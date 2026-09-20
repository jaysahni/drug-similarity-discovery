# 07 — Biologics: making the co-folding path accept a peptide or protein binder

Scope of this document: the **co-folding submission path** —
`scripts/repurpose.py` `shortlist` and `submit`. Site definition for epitopes
(`scripts/epitope_gate.py`), polymer contact extraction (`scripts/interfaces.py`),
and the corpus itself (`data/biologic_drugs.csv`, `scripts/build_biologic_library.py`)
are separate pieces of work and are referred to here only where this path depends
on them.

Every number below is labelled **measured**, **computed here**, or **EXTRAPOLATED**.
**Nothing in this document was produced by a Rowan run. Zero credits were spent.**
Verified 2026-09-19 against `rowan-python 3.2.0`, `stjames 0.0.275`, Python 3.14.7.

## Why the target-side formulation is stronger here than it was for small molecules

For small molecules there is a cheap ligand-side baseline — ECFP4 chemical
similarity — and this repo has measured that it nearly matches the structural
method. For biologics that baseline does not exist: two peptides with unrelated
sequences can engage the same epitope, and two near-identical sequences can behave
differently. Matching in **target-side residues** is not a preference here, it is
the only formulation available. That is also why `interfaces.py`'s
`consensus_signature()`, `weighted_jaccard()` and `coverage()` needed no change —
they already work in residue space and never look at what the binder is made of.

## What changed

A shortlist row now carries a `modality`, and that field decides **which input slot
the binder occupies** — nothing else. Same model, same MSA setting, same scoring
downstream, and the same prediction head left switched off
(`ligand_binding_affinity_index=None`, never read — `PROJECT_GOAL.md` §4.4).

| modality | `initial_protein_sequences` | `initial_smiles_list` | binder a constraint names |
|---|---|---|---|
| `small_molecule` (or **absent**) | `[target]` | `[drug]` | `ligand` input `0` |
| `peptide` / `protein` / `antibody` / `biologic` | `[target, binder_chain, …]` | `[]` | `protein` input `1` |

An **absent** `modality` means small molecule, so every shortlist written before
this existed keeps its meaning. Checked, not assumed: the payload built for a
small-molecule row is **byte-identical** to the one the pre-change code built
(dict comparison of `workflow_data` against a `ProteinCofoldingWorkflow`
constructed the old way), and the small-molecule selection code in `cmd_shortlist`
is unchanged apart from two lines moving above the `--biologics` dispatch (`git
diff` of the function body). `scripts/test_pipeline.py` passes 10/10,
`scripts/test_interfaces.py` 36/36, and `scripts/lint_language.py` reports 0
violations on its default targets.

## What the installed package actually supports

Read from the installed source, not from documentation:

- `stjames/workflows/workflow.py:45` — `initial_protein_sequences: list[ProteinSequence] | list[str]`.
  **Several protein sequences are accepted**, so a peptide binder is simply a second
  entry. It is a **union of two homogeneous lists**: one `ProteinSequence` forces
  every entry to be a `ProteinSequence`, which is why a `cyclic` binder makes this
  path wrap the target too.
- `stjames/workflows/protein_cofolding.py:27,45` — `Token(input_type, input_index,
  token_index)` and `PocketConstraint(input_type, input_index, contacts, …)`.
- **`input_index` indexes the list for that `input_type`, not a global input
  counter.** Evidence: `stjames/workflows/batch_protein_cofolding.py:84,91`
  range-checks a ligand-side pocket constraint against `ligand_slot_count`, a count
  of **ligand slots alone**. Supporting evidence from this repo's own paid runs: the
  two constrained probes in `results/cofold_constraint_evidence.json` used
  `input_type="ligand", input_index=0` alongside one protein and one ligand, and the
  constraint moved the ligand: aspirin went from 6/17 to 11/17 pocket residues
  engaged (axitinib was already 16/17 unconstrained and stayed there). A constraint
  that had addressed the protein instead could not have moved the ligand's pose that
  way. n=2 probes, so this is supporting evidence, not the primary one.
- **`token_index` is the 0-based offset into that chain**, so the site's author
  numbering has to be converted (`auth_resid − seq_start`). For the KDR construct
  (P35968, residues 834–1162, 329 aa) the rank-1 pocket's 17 residues map to token
  indices 6 … 213, all inside the construct. `pocket_tokens()` now hard-stops if any
  residue falls outside it instead of quietly building a constraint that points at
  nothing.
- **Therefore the target must keep protein input 0.** The binder is *appended*. Every
  contact token this file emits stays `("protein", 0, token_index)` whatever the
  modality, so adding a biologic cannot silently renumber the site.
- **`stjames` does not range-check constraint indices for a single co-fold.**
  `ProteinCofoldingWorkflow.validate_model_compatibility`
  (`protein_cofolding.py:190`) checks model/feature compatibility only; the one range
  check in the library lives in the *batch* workflow and does not run here. Verified
  by constructing a workflow with one protein input and a pocket constraint naming
  protein input 1 — it validates. `check_indices()` in `repurpose.py` now refuses it
  locally.
- `rowan-python` exposes no batch-cofolding submitter (`rowan/workflows/` contains
  `protein_cofolding.py` and no batch equivalent; `grep -l BatchProteinCofolding
  rowan/` returns nothing), so one job per binder is the only option from this client.

## The biologic shortlist

`shortlist --biologics` reads `data/biologic_drugs.csv` (path overridable with
`--biologics-csv`; column names resolved through aliases so the corpus can be
rebuilt without breaking this). It emits the **same role labels** as the
small-molecule path: `known_binder`, `decoy`, `hard_decoy`.

**Decoy design, and why.** MW/cLogP matching is meaningless for a polymer, so the
match is on **(length ⁄ 10, net charge K+R−D−E)** — the peptide analogue of
(MW ⁄ 100, cLogP):

- *Length* is the confound that actually threatens an interface score: a longer
  binder buries more target residues for reasons that have nothing to do with binding
  (`PROJECT_GOAL.md` §1.4a). A decoy the same size as the positive it replaces takes
  that explanation away.
- *Net charge* stands in for cLogP because peptide–epitope recognition is
  electrostatics-heavy and a 30-mer has no useful logP. It ignores histidine and the
  termini: a rank statistic for matching, not a titration curve.
- `--hard-decoys` narrows the pool to biologics annotated against **some other**
  target — the peptide analogue of "a `-tinib` that does not hit this kinase". The
  same caveat carries over: annotation is incomplete, so "not annotated against X" is
  not "does not bind X".
- **One confound has no small-molecule counterpart.** A decoy that is a *sequence
  relative* of a positive may well bind the same site. Anything at or above
  `--max-identity` (default 0.50, difflib ratio — a crude proxy used only to
  **exclude**, so a loose measure errs safe) is dropped, and the identity actually
  achieved is recorded on every decoy row.
- Decoys beyond one-per-positive are topped up **by distance to the nearest
  positive**, not at random (the small-molecule path draws randomly). Measured
  consequence of getting this wrong: file order put 664- and 1,334-residue antibodies
  against 31–44-residue peptide positives — the exact size confound the matching
  exists to remove, at 10–40× the co-folding cost each.

**Worked example, computed here** (`--target GLP1R --hard-decoys --n-positives 4
--n-decoys 8`, nothing submitted):

| role | drug | corpus modality | len | net charge | identity to nearest positive |
|---|---|---|---:|---:|---:|
| known_binder | DULAGLUTIDE | protein | 275 | −7 | — |
| known_binder | EXENATIDE | peptide | 39 | −3 | — |
| known_binder | LIXISENATIDE | peptide | 44 | +3 | — |
| known_binder | SEMAGLUTIDE | peptide | 31 | −2 | — |
| hard_decoy | BROLUCIZUMAB | antibody | 252 | −5 | 0.120 |
| hard_decoy | INSULIN DEGLUDEC | peptide | 50 | −3 | 0.247 |
| hard_decoy | TESAMORELIN | peptide | 44 | +4 | 0.250 |
| hard_decoy | INSULIN DETEMIR | peptide | 50 | −2 | 0.247 |
| hard_decoy | INSULIN ASPART | protein | 51 | −3 | 0.222 |
| hard_decoy | INSULIN GLULISINE | protein | 51 | −3 | 0.267 |
| hard_decoy | sermorelin | peptide | 29 | +3 | 0.333 |
| hard_decoy | TEDUGLUTIDE | peptide | 33 | −4 | 0.406 |

Decoy similarity to the nearest positive: max 0.41, median 0.25 over the eight
hard decoys in the table above (difflib ratio, n=8); the length/charge-matched arm
(`--n-decoys 4`, no `--hard-decoys`) gives max 0.31, median 0.24 (n=4). TEDUGLUTIDE — a GLP-2
analogue that is annotated against GLP2R, not GLP1R — is the kind of hard decoy this
was built to produce.

**Corpus coverage, computed here** from `data/biologic_drugs.csv` (447 rows,
`sha256:7592edc4071c`). These counts are only valid for that build of the corpus - an
earlier draft of this table was computed against a previous build and was wrong by
2 peptides; re-derive them with the snippet in `scripts/build_biologic_library.py`
if the sha above does not match the file on disk.

| | n |
|---|---:|
| rows | 447 |
| co-foldable sequence (corpus `usable_sequence=1`, modality ≠ `other`) | 301 |
| of those, parse as chains of the 20 standard residues | 301 / 301 |
| antibody / protein / peptide among them | 118 / 109 / 74 |
| not evaluated: no sequence in source | 77 |
| not evaluated: not a polypeptide (oligonucleotide, gene, cell) | 51 |
| not evaluated: non-standard residue unmappable | 13 |
| not evaluated: unresolved residue letter | 5 |
| distinct target symbols with ≥1 co-foldable approved biologic | 145 |
| target symbols with ≥2 approved **peptide** binders (a screenable positive set) | 22 |

Those exclusions are the **corpus's own per-row verdict**, used rather than
re-derived, so the decision about each non-standard residue stays auditable where it
was made. `binder_chains()` independently refuses anything outside the 20 standard
residues: a lipidated or Aib-containing peptide has to be declared as a stjames
`ResidueModification` with a CCD code, and submitting it without the modification
would be co-folding a different molecule from the drug.

## What it costs

**Measured** (small molecule, this receptor, from
`results/pipeline/colorectal-cancer/repurpose_state.json` and
`results/cofold_constraint_evidence.json`):

- 42 collected co-folds of the same 329-residue receptor: **4.55 credits each**
  (sd 1.02, 95% CI 4.23–4.87, range 3.06–6.96, **n=42**).
- Cost is **flat in ligand size** over 2–85 heavy atoms: Pearson r = −0.090,
  p = 0.573, n=42 (Spearman ρ = −0.037, p = 0.815). **The receptor dominates.**
- Four separately timed probes billed at **0.0500 credits/second**, with no spread
  (n=4). A credit is wall-clock time; 4.55 credits ≈ 91 s.

**EXTRAPOLATED** (peptide — nothing with a polymer binder has been run):

A 9–50 residue binder adds 9–50 tokens to a 329-token receptor (+3% to +15%).
Ligand-size independence says the receptor dominates, so the two honest brackets are:

| binder | tokens | flat (cost ∝ receptor) | quadratic (cost ∝ tokens²) |
|---|---:|---:|---:|
| 15 aa (median approved peptide) | 344 | 4.55 | 4.97 |
| 39 aa (exenatide) | 368 | 4.55 | 5.69 |
| 50 aa (longest "peptide" in the corpus) | 379 | 4.55 | 6.04 |
| 664 aa (a 2-chain antibody) | 993 | 4.55 | 41.45 |

**≈5–6 credits for a peptide co-fold, confidence LOW on the exact number and moderate
on the order of magnitude.** The unquantified term is the **second MSA**:
`use_msa_server=True` now has two chains to search and nothing here measures what
that costs, and short peptides are exactly the case where an MSA server is slow and
unhelpful. The exponent is an **assumption**, exposed as `--token-exponent` (default
2.0) and used **only to size the budget guard** — no number produced by it should be
reported as a cost. One real peptide co-fold would settle it, and `max_credits` caps
the downside at 3× the guard.

For the 12-row GLP1R shortlist above the dry run totals **115.4 credits** under the
quadratic guard versus 12 × 4.55 = **54.6** at the measured flat rate. Against the
~103 credits remaining, **a 12-row biologic screen that includes the two
antibody-sized rows does not fit the budget under the pessimistic assumption**. The
10 rows at ≤ 51 aa come to **76.4 guarded / 45.5 flat** — which does fit, and is the
shortlist worth running first.

## Implemented vs planned

| | status |
|---|---|
| peptide / protein binder as a second protein chain | **implemented**, dry-run validated, never submitted |
| multi-chain binder (insulin A+B, a Fab) — `/`, `;` or `+` separated | **implemented**; insulin degludec builds 3 protein inputs (329 + 29 + 21) |
| cyclic binder (`ProteinSequence(cyclic=True)`) | **implemented**; forces the target into a `ProteinSequence` too, flag comes from the corpus's `is_cyclic` |
| pocket constraint on a polymer binder (`protein` input 1) | **implemented**, unconstrained still the default |
| pocket constraint on a **multi-chain** binder | **refused with a reason** — one `PocketConstraint` names one binder input, and which chain sits on the epitope is a decision nothing here has measured |
| local index range checks stjames does not do | **implemented** (`check_indices`) |
| `--dry-run` payload printing + local pydantic validation | **implemented** |
| target-vs-shortlist mismatch guard | **implemented** — refuses to co-fold a GLP1R shortlist against the KDR construct |
| shortlist-replacement guard | **implemented** — a plain `shortlist` over a state with submitted jobs would orphan them; needs `--append`, `--out` or `--replace-shortlist` |
| non-standard residues via `ResidueModification` | **planned** — refused with a reason today; 16 corpus rows blocked on it |
| disulfides via `BondConstraint` | **planned** — the corpus records `n_disulfide`; nothing is emitted, so a cyclic peptide is submitted with its backbone only |
| **scoring a polymer pose** | **NOT in this path.** `cmd_score` still extracts contacts through `drug_like_ligands()` / `ligand_contacts()` (`repurpose.py:991,997`), which see HET groups only. `interfaces.py` has gained `chain_contacts()` and polymer chain detection; wiring `score` to it is the next step and is owned elsewhere |
| batch co-folding (one workflow, many binders) | **not available** in `rowan-python 3.2.0` |

## Not evaluated / why

| item | why |
|---|---|
| Any real peptide co-fold | Would spend credits. Budget is ~103 and the goal was to spend nothing. |
| The actual credit cost of a polymer binder | Same. The figure above is an extrapolation from small-molecule runs, explicitly bracketed. |
| Whether a second chain's MSA adds materially to the cost | Nothing measured it; it is the largest unknown in the estimate. |
| How Boltz-2 names the binder chain in the returned pose | No pose has been returned. The scorer will need to identify the binder chain, and this is the first thing to check on the first real run. |
| Whether a pocket constraint helps or hurts a **peptide** | The constrained-vs-free measurement is small-molecule only (n=2 drugs, one target). A peptide's footprint is the size of the site itself, so that result does not transfer. Default stays unconstrained. |
| A biologic screen on either existing pipeline target | **Verified programmatically:** CDK2 has **0** approved biologics annotated against it in the corpus, and KDR has exactly **1** — ramucirumab, a 1,320-residue antibody. One positive is not a validation set, and the KDR pipeline's construct is the intracellular kinase domain (P35968, 834–1162), which an antibody raised against the receptor's extracellular region cannot engage. A biologic screen needs a new target whose site is the binder's actual site; the 22 symbols with ≥2 approved peptide binders are the candidate list. |
| Whether P2Rank's site is even the right site for a biologic | `scripts/epitope_gate.py` (separate work) is the measurement; `PROJECT_GOAL.md` §1.4b predicts it does markedly worse on a flat epitope, and `sm_addressable` is still not computed. |

## Reproduce, offline, for nothing

```bash
# a biologic shortlist, written to a file so the pipeline state is untouched
./env/bin/python scripts/repurpose.py shortlist --biologics --target GLP1R \
    --hard-decoys --n-positives 4 --n-decoys 8 --out /tmp/bio_glp1r.json

# every payload built and validated against the stjames models; nothing sent
./env/bin/python scripts/repurpose.py submit --dry-run --shortlist /tmp/bio_glp1r.json \
    --allow-target-mismatch --max-credits 400
```

`--allow-target-mismatch` is needed **only** because this repo has no GLP1R pipeline:
the run pairs GLP-1 binders with the KDR construct to exercise the payload builder,
prints `!! TARGET MISMATCH` on every row, and is not a biological claim. Without the
flag the mismatch is a hard stop.

Console excerpt (constrained variant, to show where the indices land):

```
  DRY RUN  known_binder EXENATIDE                              modality=peptide
      protein input 0:  329 aa  TARGET KDR
      protein input 1:   39 aa  BINDER
      constraints: pocket on protein input 1, 17 target contact tokens, <= 6.0 A
      est. credits 7.5 (EXTRAPOLATED: 329 target + 39 binder tokens, x1.25 on an
                        assumed tokens^2.0 scaling); max_credits=22;
                        payload validated by stjames, nothing sent
      {
       "workflow_type": "protein_cofolding",
       "workflow_data": {
        "initial_protein_sequences": ["LKLGKPLGRGAFGQVIE…VEHLGN",
                                      "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS"],
        "initial_smiles_list": [],
        "use_msa_server": true,
        "pocket_constraints": [
         {"input_type": "protein", "input_index": 1,
          "contacts": [{"input_type": "protein", "input_index": 0, "token_index": 6},
                       {"input_type": "protein", "input_index": 0, "token_index": 14}, …],
          "max_distance": 6.0, "force": false}],
        "model": "boltz_2"
       },
       "name": "GLP1R cofold - EXENATIDE", "max_credits": 22, "is_draft": false
      }
```

(The sequences are elided here with `…`; the command prints them in full — that is
the exact body `rowan.submit_protein_cofolding_workflow` would POST, assembled the
same way `rowan/workflows/protein_cofolding.py` assembles it.)

**Proof that a dry run cannot spend anything**, run in-process: after
`repurpose.main()` builds all 12 payloads, `'rowan' in sys.modules` is `False` and
`ROWAN_API_KEY` was never read into the environment. `rowan` is the only HTTP client
on this path and it is imported only in the non-dry-run branch.
