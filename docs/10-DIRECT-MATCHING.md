# 10 — Direct matching: generate a candidate, find the FDA drug it resembles

The pipeline in `docs/06` answers a different question from the one intended, and
this note records the switch, what it cost, and what it found.

## The two questions

**Footprint matching** (PROJECT_GOAL.md v0.2, built first): BoltzGen designs a
binder, the binder is *discarded*, and what is kept is the set of target residues
it touched. Approved drugs are co-folded and ranked by how much of that residue
set they also touch. Similarity is **footprint ↔ footprint**.

**Direct matching** (this note): generate a candidate, and find the approved drug
most similar **to the candidate itself**. Similarity is **molecule ↔ molecule**.

Direct matching only works when the generated thing and the retrieved thing are
the same *kind* of object. That is why it needed two arms and a new target.

## Why the switch

Three measured results, not a preference:

| Finding | Number | Source |
|---|---|---|
| BoltzGen consensus loses to a 0.47 s pocket finder | Jaccard 0.338 vs 0.494, Δ −0.155, Holm p=0.0016, n=31 | `results/m2_gate_CDK2.json` |
| Hotspot coverage does not separate CDK2 binders from decoys | 0.507 vs 0.427, Mann-Whitney p=0.229, n=5/tier | `results/demo/cdk2/05_ranked.json` |
| The top-ranked drug is a decoy | atorvastatin 1st, palbociclib 9th of 14 | same |

Boltz-2 pose confidence *did* separate the tiers (ipTM 0.962 vs 0.878, p=0.0159),
but that is pose confidence rather than hotspot coverage, PROJECT_GOAL.md §4.4
deliberately keeps confidence out of the ranking path, and it is the quantity
most exposed to leakage. Recorded as an observation, not promoted to the ranking.

## Target: thrombin (F2, P00734), not CDK2

CDK2 cannot host this pipeline. No approved peptide binds it, and its ATP slot is
too enclosed for a designed peptide to enter. Thrombin is the one target where
**both** arms have an approved positive control, verified against the corpora:

| Arm | Positive control | Evidence |
|---|---|---|
| P (peptide) | bivalirudin (20 aa), lepirudin (65 aa) | `data/approved_biologics.csv`, `target_gene_symbol = F2` |
| S (small molecule) | argatroban, ximelagatran, dabigatran etexilate | `data/approved_drugs.csv`, `moa_targets = F2`; 18 drugs annotated to F2 |

Structure: **1PPB chain B**. 1PPB was chosen over the higher-resolution
4UD9/4UE7/5AFY because those carry a hirudin fragment as a **polymer** chain,
which `remove_heterogens` does not strip — designing against one would repeat the
occupied-pocket mistake with a peptide instead of a ligand.

## Three things that silently rename what you asked for

Each of these produces output that looks entirely ordinary.

**1. Insertion codes.** Thrombin uses chymotrypsin numbering. 1PPB chain H holds
**259 residues under only 231 distinct residue numbers** — 28 carry an insertion
code, including the whole 60-loop (60A–60I) that lines the active site. Keying
residues on `int(line[22:26])` merged them and lost 28, landing exactly on the
residues that matter. `contacts.residue_label` now keys on `(resseq, icode)`;
CDK2 has no insertion codes so its labels stay ints and its artifacts stay
comparable. Guarded by six tests: PPACK's contacts must include Ser195, Asp189
and the 60-loop cage.

**2. `Protein.prepare` renames chains.** 1PPB's H (catalytic) and L (light) come
back as **B** and **A**. The chain letter in `TARGETS` is a hint about the
original PDB, not a fact about the prepared one. `pipeline.resolve_chain` now
finds the target chain by aligning each chain to the canonical sequence.

**3. It also renumbers.** Residues come back numbered from 1, keeping insertion
letters (`22A`, `47A`), so no constant offset recovers author numbering.
`contacts.align_to_reference` does Needleman-Wunsch against the canonical
sequence and refuses below 95% identity.

## Pocket detection cannot find a protease active site

**MEASURED**: Rowan returns 6 pockets for 1PPB chain B and **not one contains
Ser195 or Asp102**. The best overlap with the UniProt-annotated catalytic triad
is 1 of 3, and two pockets tie at it.

This is a fact about the target, not a bug in the detector. A protease active
site is a shallow groove across S1/S2/S3 subsites, not an enclosed cavity, so a
cavity finder fragments it. CDK2's ATP site, being a real cavity, overlapped 11
of 20.

Stage 2 therefore falls back to a **geometric site** when pocket overlap is below
threshold, and records that it did. The anchor is the catalytic nucleophile
Ser195 alone — His57 and Asp102 sit *behind* the serine, so a sphere around all
three reaches into the core and away from where a ligand binds. Radius 6 Å gives
20 residues, sized for the 8–16-mer we then ask for.

Both anchor and radius are set from enzymology and from the intended binder's
size, not tuned against the co-crystal ligand. **Validation, computed
afterwards**: the site recovers 10 of PPACK's 23 contact residues, Jaccard 0.303.

## The null is the whole point

With thousands of approved entries, **every query has a nearest neighbour**. A
bare "0.42 Tanimoto to X" is not a result. `demo/match_direct.py` supplies the
missing piece — `visualization/fixture.py:35` had been shipping the empty slot
for it (`"No scientific null was computed"`).

Precedent, and the significance test it hands us:

- **Davis FP & Sali A**, *PLOS Comput Biol* 6(2):e1000668, 2010
  ([doi](https://doi.org/10.1371/journal.pcbi.1000668)) — binding-site overlap as
  the fraction of interface residues aligning to ligand-binding residues, tested
  by Fisher's exact test against random independent placement.
  `demo/embed.py:overlap_pvalue` implements that null exactly.
- **Rácz, Bajusz & Héberger**, *J Cheminform* 10:48, 2018
  ([doi](https://doi.org/10.1186/s13321-018-0302-y)) — 44 similarity measures
  benchmarked on interaction fingerprints. Six beat Tanimoto; **cosine is not
  among the recommended ones**, so `tanimoto_core` is reported beside it.

### The calibration test earned its place immediately

Held-out queries drawn from the decoy distribution must score uniform
percentiles. A KS test said otherwise for Tanimoto: **D=0.187, p<0.001, mean
percentile 0.457** instead of 0.5. Cause: Tanimoto over sparse fingerprints is
heavily discrete, and `percentile_of` used a strict `<`, so every exact tie was
scored against the candidate. `scripts/metrics.py:10-20` warns about the same tie
problem for 2048-bit ranking.

Now mid-rank — `P(null<v) + 0.5·P(null==v)` — and both metrics calibrate: cosine
D=0.062 p=0.416, Tanimoto D=0.090 p=0.074. Without the fix every candidate would
have looked less remarkable than it is, consistently and invisibly.

**A null must also be independent of what it scores.** The first cut of arm S
built the null from the same molecules it was scoring, which only recovers each
molecule's rank within its own set. Arm S now splits by pocket evidence: docked
molecules are queries, undocked ones the null. When docking has not run the
output says `DEGENERATE` rather than reporting ranks as a test.

## Arm P — peptide

`demo/peptide_arm.py`. BoltzGen `peptide-anything` at 8..16 residues → ESM-C
cosine against approved biologics ≤100 aa (37 with usable sequences, 2 annotated
to F2, base rate 5.4%).

**The protocol change fixed the designs.**

| | CDK2, `protein-anything` 60..90 | thrombin, `peptide-anything` 8..16 |
|---|---|---|
| ipTM | 0.134–0.197 (mean 0.162) | **0.525–0.640 (mean 0.598)** |
| filter funnel (`budget`) | off, 1:1 | **on, 5:1 — 20 designs → 4** |

That confirms the geometric mismatch in `docs/07` §3.4: a 60–90-mer scores badly
on an enclosed ATP pocket because it cannot enter it. 0.60 is also roughly where
published miniprotein campaigns top out, so it is a normal number, not a good one.

**The match nonetheless fails.**

- Every design's nearest neighbour is **abarelix**, a GnRH antagonist
- Null percentiles **0.947, 0.463, 0.423, 0.287** — three of four sit inside the
  null of composition-preserving shuffles
- **Enrichment factor 0.0** for all four; bivalirudin and lepirudin rank 10–17 of 37
- ESM-C and BLOSUM62 identity disagree on the nearest neighbour for all four

**Why that is credible rather than a shrug — the positive control passes.**
bivalirudin → lepirudin ranks **1 of 36** (cosine 0.9179); lepirudin →
bivalirudin ranks 5 of 36. ESM-C cosine *does* retrieve co-target peptides in
this corpus. The metric works; the designs do not pass it. Without the control,
"our designs failed" and "our metric is broken" are indistinguishable.

**And the hub check explains the abarelix result.** 77% of randomly shuffled
sequences also land on abarelix — it is an embedding hub, so a design landing
there carries no information about the design.

A correction worth recording: ESM-C cosine is **not** compressed into a narrow
band. Per-pair cosines span 0.211–0.985 (sd 0.242) for one design. It is the
*null* that sits near 0.97, because the null is a distribution of **maxima over
37 entries**. Read the percentile, not the raw cosine.

## Arm S — small molecule

`demo/generate.py`. `entropy/gpt2_zinc_87m` (MIT, GPT-2 87M on ZINC) → RDKit
drug-likeness filter → Rowan batch docking into the thrombin site → ECFP4-2048
Tanimoto against 4,099 approved drugs (18 annotated to F2, base rate 0.44%).

591 of 600 generated molecules survive filtering (1 invalid, 8 below QED, 0
already approved). Fingerprints come from `scripts/representations.py`'s registry,
so the arm uses the exact representation the benchmark found best of 16.

### Result: nothing found, and both controls say why

**Nothing resembles a thrombin drug.** Best nearest-neighbour Tanimoto is 0.352
(idrocilamide, no F2 annotation) at null percentile 0.908. The whole top 10 is
`novel_pairing`. Null over the 567 undocked generated molecules: mean 0.283,
p95 0.379.

**Docking adds nothing measurable.** Among the 23 generated molecules that
docked, Vina score and similarity-to-known-F2-binders are uncorrelated:
**Spearman ρ = −0.084, p = 0.703, n = 23**.

**But Vina works.** The three known thrombin drugs, docked in the *same batch*,
beat the generated median: argatroban −6.217, ximelagatran −5.899, dabigatran
etexilate −6.414.

**And ECFP4 works — which took a correction to see.** Reading the positive
control as a top-5 list nearly produced a false negative: argatroban's five
nearest approved drugs are all peptidomimetics (angiotensin II, icatibant,
lisinopril…) and not one is a thrombin drug. The *rank* of the co-mechanism
drugs says the opposite. Out of 4,098, against ~2,050 expected by chance:

| query | ranks of the other F2-mechanism drugs |
|---|---|
| argatroban | bivalirudin **10**, ximelagatran **48**, dabigatran etexilate 1138 |
| ximelagatran | bivalirudin **29**, dabigatran etexilate **30**, argatroban **51** |
| bivalirudin | argatroban **66**, ximelagatran **82**, dabigatran etexilate 632 |
| dabigatran etexilate | ximelagatran **17**, argatroban 988, bivalirudin 1033 |

Dabigatran etexilate is the outlier in both directions, which is what a prodrug
ester should be. The control now reports ranks alongside the top-5, because the
top-5 alone is misleading.

This **replicates README result #3 on a second target**: chemical similarity
finds the drugs that *look* like the query — argatroban retrieves peptidomimetics
generally — and the co-target drugs it does find, it finds at rank 10–82 rather
than rank 1.

**What this is: virtual screening with a generative front end.** The generator is
not conditioned on the binding site; the pocket enters only as a docking filter.
It must not be described as structure-based de novo design.

That was forced. Every pocket-conditioned generator — DiffSBDD, Pocket2Mol,
TargetDiff, PocketFlow, DecompDiff, LiGAN — depends on `torch-scatter` /
`torch-sparse` / `torch-cluster`, which have **no build for torch 2.14 and have
never shipped a native macOS arm64 wheel** (verified against the PyG wheel index;
upstream issue open and unanswered). Rowan hosts no de novo generator either: all
41 workflows take molecules as input, and `analogue_docking` requires the caller
to supply the analogues.

## Not evaluated, and why

| Approach | Why not |
|---|---|
| ESM-C on the design vs fingerprints on the drugs | No protein language model embeds small molecules. ESM-C is protein-sequence-only in every variant (300M/960-d, 600M/1152-d, 6B/2560-d). Undefined, not inaccurate. |
| ESM-2 (the toolkit's model) on a design | ESM-2's training set **explicitly excluded de novo designs** — 1,027 sequences tagged "artificial sequence" plus 58,462 similar to 81 known designs. Out of distribution by construction. |
| ConPLex | A genuinely shared 1024-d space, but trained on a *binding* objective with DUD-E decoys chosen to be chemically similar yet non-binding — the objective deliberately decorrelates chemical similarity from embedding distance. It also takes the protein as the *target*, inverting our query. Feb-2024 pre-release, Python ≤3.11. |
| DrugCLIP | Needs a 3D pocket to encode; a designed peptide has none. Google-Drive checkpoint, self-described "raw" code. |
| MolTrans, HyperAttentionDTI | No separable per-entity embedding — they score a drug–protein *pair* through a joint network, so there is nothing to retrieve against. |
| Antibody arm on thrombin | Thrombin is extracellular so biologics can reach it, but the corpus has only 3 F2-annotated entries and 2 are short enough to compare against an 8–16-mer. Reported as the peptide arm instead. |

## Environment

`novakit` requires Python ≥3.12,<3.13 and the project `env/` is 3.14.7, so all of
this runs in the gitignored `env-kit/` (3.12.13, `uv`). Added there for these
arms: `biopython` (BLOSUM62 local alignment), `torch`/`transformers` (ESM-C and
the SMILES generator), `scipy`.

## Standing caveats

- Similarity is not affinity. Nothing here predicts binding.
- Every output is a computationally-ranked hypothesis requiring experimental
  validation.
- n is small throughout: 4 designs in arm P, 37 peptides in its corpus, one
  target, one structure. No cross-target claim is made.
