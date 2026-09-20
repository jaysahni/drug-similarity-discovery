# 08 — Making the pipeline accurate on all modalities

**Scope of this document.** An executable plan for extending the interface-signature
pipeline from small molecules to peptides and protein/antibody biologics: what
defines the site, what folds the complex, what scores it, what gates it, and what
validates it — per modality. Plus a ranked experiment queue, an explicit
scoped-out list, and the places where the field has no answer and our own run is
the contribution.

**Evidence labelling used throughout.** **[M]** measured in this repo, with the
file that holds the number. **[L-a]** a number reported on a named benchmark in a
primary source. **[L-b]** an abstract-level or secondary-source claim. **[?]**
unverified — do not build on it without checking first.

Everything under **[M]** is reproducible from files on disk today. Nothing in this
document is an expected result written ahead of a run, except where a section is
explicitly labelled *prediction* — and each of those states its own falsifier.

---

## 0. The one-paragraph answer

Site definition has to be replaced for biologics, and we have measured that
ourselves rather than assumed it: P2Rank top-1 recovers a biologic interface at
Jaccard 0.1042 (n=200) against 0.6355 on the KDR small-molecule site (n=38), and
on antibody epitopes specifically it is **indistinguishable from a size-matched
random surface patch** (0.0255 vs 0.0270, n=58, Holm p=1) **[M**
`results/epitope_gate.json`**]**. But the modality label turns out to be the wrong
axis. Peptide sites sit at 0.1836 and PPI sites at 0.0903 — a concavity gradient,
not a modality cliff — and for peptides P2Rank usually *offers* a good pocket and
then mis-ranks it (best-at-any-rank oracle 0.3112, precision 0.619). Separately,
ranking: on the KDR board, unconstrained Boltz-2 **ipTM separates known binders
from unrelated drugs at AUC 0.995 and from other kinase inhibitors at only 0.696
(p=0.127, n=41)** **[M**, computed for this plan from
`results/pipeline/colorectal-cancer/repurpose_state.json`**]** — the same hard-null
ceiling our interface-overlap score hits at 0.717 (p=0.092). So the honest scope
is: **rank on small molecules; rank on peptides only after the hard-null test in
E7 passes; do not rank on antibodies at all** — look their epitopes up instead of
predicting them.

---

## 1. Constraints verified today, before anything is built on them

| Fact | Status | How it was checked |
|---|---|---|
| Rowan bills **exactly 3.0 credits/GPU-minute** | **[M]** | `results/cofold_constraint_evidence.json`: 102.0 s → 5.10 cr; 222.1 s → 11.11 cr; 261.8 s → 13.09 cr. Three runs, three exact matches. Costs below are derived from wall-time estimates through this constant. |
| A small-molecule co-fold costs **~5 cr / ~1.7–4 min** | **[M]** | same file; 41-drug board averaged ~4.4 cr/drug. |
| Rowan's returned score object is **`{ptm, iptm, avg_lddt, confidence_score}`** — **no PAE matrix** | **[M]**, for small-molecule co-folds only | every `jobs[*].scores` in `repurpose_state.json` and every probe in `cofold_constraint_evidence.json` has exactly these four keys. **Consequence: ipSAE, Pinc, LIS, MiniPAE, pDockQ2 and AF2-style `pae_interaction` are all unavailable to us as written.** Whether a *protein–protein* co-fold returns more is untested — that is experiment **E6**, and it gates the entire interface-restricted-confidence branch. |
| `torch 2.14.0` is installed and a **cp314 arm64 `gemmi` wheel exists** (0.7.5) | **[M]** | `env/bin/python -c "import torch"`; `pip download gemmi --no-deps` fetched `gemmi-0.7.5-cp314-cp314-macosx_11_0_arm64.whl`. **PeSTo's two hard dependencies are satisfiable locally.** This closes the single biggest open question about the E1 path. |
| `cc`, `gcc`, `make`, `java` present; **`cmake` absent** | **[M]** | `which`. fpocket ships a plain Makefile, so it should build; PocketMiner and anything CMake-based would need cmake first. |
| P2Rank 2.5 lives under the **session scratchpad**, not in the repo | **[M]** | `scripts/run_p2rank.py:45`. This is fragile across sessions — anything in the queue below that re-runs P2Rank must re-materialise it. |
| Modal: **CPU works, every GPU tier is payment-blocked** | **[M]** | `docs/03-SCOPE-AND-CONSTRAINTS.md:20-21`, four tiers probed. **Any Modal GPU method in §5 needs a payment method added, not just a key.** |
| Rowan credits | ~103 at the time `epitope_gate.json` was written; **2000+ now** per the session brief | the not-evaluated entries in `epitope_gate.json` that say "no credits remain" are stale and should be re-read as "not yet run". |
| PeSTo weights are **CC BY-NC-SA 4.0** | **[L-b]** | conflicts with this repo's MIT posture **if weights are vendored**. Plan: script against them, report numbers, never commit weights, and state the dependency in the Regeneron submission. |

---

## 2. The metric change that comes before any new modality

Do not report cross-modality Jaccard. This is arithmetic, not taste: a 12-residue
drug footprint *perfectly contained* in a 40-residue epitope scores J = 0.30, which
is indistinguishable from genuine partial overlap. Our small-molecule numbers
(0.635) and any biologics number are therefore on different scales and must never
be pooled or compared.

We already have the right function. `scripts/interfaces.py:1314 coverage()`
returns `core_coverage`, and its own docstring calls it *the primary score*. Make
it the headline everywhere and demote `weighted_jaccard` to a same-modality
diagnostic. Precedent for the asymmetric choice: Davis & Sali (PLoS Comput Biol
2010) quantified small-molecule-vs-protein site overlap as the fraction of
interface residues aligned to ligand-site residues **[L-b]**.

Three additions, in cost order:

1. **Overlap coefficient** `|A∩B| / min(|A|,|B|)` alongside coverage. ~5 lines.
2. **Min-BSA-weighted overlap**: per-residue `min(BSA_a, BSA_b)` normalised by the
   smaller interface, the AbLang-PDB convention **[L-b]**. We already compute
   `buried_area_by_residue` and `buried_area_by_residue_chains`
   (`interfaces.py:486,928`), so this is ~20 lines and it stops a residue
   contributing 5 Å² counting the same as one contributing 60 Å².
3. **DockQ v2 `fnat`** as an external cross-check. pip-installable, CPU-only,
   CAPRI-standard, covers multimers and small molecules **[L-a]**. Reporting it
   next to our own coverage costs nothing and pre-empts "you invented a metric".

And: report **three tables, one per modality, each with its own n and its own
null**. The cross-modality claim we are allowed to make is "the same target-side
residues are engaged", never "these scores are comparable".

---

## 3. The plan, per modality

### 3.1 Small molecule — ship it, with the hard null named as hard

| Stage | Decision |
|---|---|
| **Site definition** | **Keep P2Rank 2.5.** Measured J 0.635 / recall 0.744 on 38 KDR co-crystals, 0.494 on 31 CDK2 **[M** `results/m2_gate_KDR.json`, `m2_gate_CDK2.json`**]**. Revisit only if E1 shows PeSTo's ligand head beats it. |
| **Co-folding** | Boltz-2 on Rowan, **unconstrained**. The pocket constraint is measured to manufacture confidence: aspirin goes 6/17 → 11/17 pocket residues and ipTM 0.847 → 0.965 when constrained, while costing 2× **[M** `cofold_constraint_evidence.json`**]**. |
| **Scoring** | `precision_in_core` / `core_coverage` against the P2Rank signature. Unchanged. |
| **Confidence gate** | ipTM as a **structure-sanity filter only**. See E0: it is a superb easy-null discriminator (0.995) and a poor hard-null one (0.696). Candidate upgrade, **not adopted and not runnable without a written ruling** — see E8. <!--LANG-EXEMPT-->Boltz-2 exposes a binder-vs-decoy head (`affinity_probability_binary`, MF-PCBA AUROC 0.812 vs ipTM 0.566 **[L-a]**). `PROJECT_GOAL.md` §4.4 removes affinity from the ranking path, G8 bans the vocabulary, and `scripts/repurpose.py:291` hard-codes `ligand_binding_affinity_index=None`; nothing here changes that. A *gate* is not a *ranker*, but the distinction has to be decided by Evan in writing, not slipped in.<!--/LANG-EXEMPT--> |
| **Validation** | KDR 41-drug board (10 known binders / 12 hard decoys / 19 decoys) + CDK2 second target. Both nulls always reported, and the easy one labelled easy. |

### 3.2 Peptide — keep in scope, gate it, and do not promise a leaderboard yet

| Stage | Decision |
|---|---|
| **Site definition** | **P2Rank + a re-ranker**, not a replacement — provisional, pending E2. The measurement says the detector is not the problem: top-1 J 0.1836 → best-at-any-rank 0.3112 with precision 0.619, and a pocket touches the true interface in 78.6% of cases (n=70) **[M** `epitope_gate.json`**]**. Fallback if the re-ranker fails: PeSTo's protein-protein head (E1). |
| **Co-folding** | Boltz-2 on Rowan, same path as small molecules. `docs/07-BIOLOGICS.md` already shows the binder occupies `initial_protein_sequences[1]`. **[?]** Boltz-2 has **no published protein–peptide DockQ benchmark at all** — every peptide accuracy number available is AF3's or Boltz-1's. This is a gap in the field, not just in our search. |
| **Scoring** | `core_coverage` + overlap coefficient. **Never** the protein–protein DockQ bands: report under CAPRI *peptide* criteria (fnat 0.2/0.5/0.8 with matched L-RMSD / I-RMSD) **[L-a]**. |
| **Confidence gate** | **ipTM ≥ 0.8**, recomputed from PepPCBench's 52,200 released AF3 predictions: the gate keeps 67% of complexes and 90.9% of those recover ≥ half the native contacts, vs 73.6% ungated **[L-a**, recomputed from primary released data; tables at `…/scratchpad/peppcbench/`**]**. Everything below the gate goes in the *not evaluated / why* table with its ipTM — never into the ranking with a fabricated engagement score. Justification for gating hard: median fnat is **0.000** in the incorrect-DockQ band and 0.259 in the acceptable band (Spearman fnat–DockQ 0.919) — there is no regime where the pose is wrong but the contacts are still informative. |
| **Validation** | **GLP-1R (P43220)**, 45 peptide-bound entries **[L-a**, RCSB query run by a researcher 2026-09-19**]**: exenatide (7LLL), semaglutide (7KI0) and tirzepatide (7FIM/7RGP) as solved positives; **liraglutide, lixisenatide, dulaglutide, albiglutide as a genuine zero-structure holdout**; and 12+ small-molecule agonist complexes at a *different* site on the same receptor, which makes it a cross-modality site-discrimination test in one target. Second target PTH1R (51 entries; teriparatide, abaloparatide). **KDR is confirmed vacuous for peptides** — 3 entries with any 5–50-residue chain, all NMR TM fragments **[L-a]**. |
| **Negatives** | Three tiers, mirroring the small-molecule design. Tier 0 unrelated approved peptides (eptifibatide, oxytocin, setmelanotide). Tier 1 length-matched, unrelated receptor (tesamorelin, corticorelin, sermorelin). **Tier 2, the real test: approved class-B1 peptides for sibling receptors** (secretin, teriparatide, abaloparatide, native GIP). Match the length distribution explicitly and report the match, or length alone wins the AUC. Treat glucagon and oxyntomodulin as **graded, not binary** negatives — both have real GLP-1R activity. |

### 3.3 Antibody / nanobody — lookup, not prediction; and no ranking

**Scoped out of the ranking path. Kept as a site-definition and coverage arm
only.** Three independent lines of evidence, one of them ours:

1. **Ours [M]:** P2Rank on antibody epitopes is chance — J 0.0255 vs a random
   surface floor of 0.0270 (n=58), and the best-at-any-rank oracle only reaches
   0.062, so no re-ranker rescues it. 32.8% of antibody receptors returned **no
   pocket at all**.
2. **Ranking [L-a]:** 519 therapeutic mAbs from Thera-SAbDab, cognate vs
   non-cognate antigen, Boltz-2 — **every confidence metric scored ROC-AUC between
   0.50 and 0.60**, and feature aggregation did not help.
3. **Poses [L-a]:** on post-cutoff SAbDab complexes, Boltz-1 reaches 4.1%
   high-accuracy / 20.4% acceptable (n=49 antibodies). A leaderboard built on
   co-folded antibodies would be mostly docking failures wearing ranks.

| Stage | Decision |
|---|---|
| **Site definition** | **Experimental lookup.** Thera-SAbDab → SAbDab → the deposited complex → 4.5 Å heavy-atom contacts through `interfaces.py:947 chain_contacts`. Zero prediction error, zero GPU, zero credits. `docs/03` already makes `source` a discriminator on a signature, so this is a swap, not a new concept. |
| **Co-folding** | **None in the ranking path.** |
| **Scoring** | Epitope-overlap coverage against the *annotated* epitope, reported as a coverage table, not a leaderboard. |
| **Confidence gate** | N/A — nothing is being gated, because nothing is being ranked. |
| **Validation** | The 58 antibody complexes already in `epitope_gate.json`; optionally the Pierce-lab 67-case docking benchmark (**CC BY 4.0, confirmed**) for bound-vs-unbound. |
| **Best available discrimination test** | HER2: trastuzumab (domain IV) vs pertuzumab (domain II) — two approved antibodies, **non-overlapping epitopes on one antigen**, co-complexes solved (6OGE, 8PWH, 8Q6J) **[L-a]**. If our coverage metric cannot tell those two epitopes apart, the biologics arm has no signal at all. |

**Discrepancy to resolve before citing AsEP.** Two researchers reported different
best-method numbers on the same benchmark: MCC 0.305 (ratio) / 0.152 (group) in
one report, 0.210 / 0.077 in the other. **[?]** Neither is safe to quote until the
AsEP table is read directly. The qualitative conclusion — antibody-specific epitope
prediction is in the MCC 0.1–0.3 range, i.e. not rankable — survives either way.

---

## 4. Ranked experiment queue

Ordered by (decisiveness ÷ cost). The first five cost **zero credits**.

### E0 — Does co-folding confidence discriminate better than our own score? · **0 credits · DONE**

*Already computed while writing this plan*, from `repurpose_state.json` +
`repurpose_colorectal-cancer.json`, on exactly the 41 interface-scored drugs
(10 known binders / 12 hard decoys / 19 decoys), AUC with a 20,000-draw
two-sided permutation test:

| score | vs unrelated decoys | vs other kinase inhibitors |
|---|---|---|
| **ipTM** | **0.995** (p < 0.0001) | **0.696** (p = 0.127) |
| confidence_score | 0.845 (p = 0.002) | 0.608 (p = 0.403) |
| avg_lddt | 0.768 (p = 0.019) | 0.529 (p = 0.833) |
| pTM | 0.589 (p = 0.453) | 0.608 (p = 0.409) |
| *interface overlap (`precision_in_core`)* | *0.913 (p = 0.0003)* | *0.717 (p = 0.092)* |

**This is uncomfortable and it goes in the README.** Against the easy null, a
confidence field we already get for free **beats** our structural score. Against
the hard null, both fail, at the same place, with overlapping numbers. The
defensible claims that survive are: (i) interface overlap is *interpretable* —
it names the residues, ipTM does not; (ii) the hard-null ceiling is shared, so it
is a property of the co-folding backend, not of our scorer; (iii) our earlier
"confidence does not discriminate" framing was drawn from the **constrained**
aspirin probe (ipTM 0.965 vs axitinib 0.992) and does **not** generalise to
unconstrained runs — the constraint is what destroyed the signal.
**Action:** commit this as `scripts/confidence_gate.py` so it is reproducible, and
correct the README claim.

### E1 — PeSTo's ligand head vs P2Rank, on ground truth we already own · **0 credits · CPU · ~1 day**

*Question.* Does one structure-only model (PeSTo, Nat Commun 2023 — separate
protein-protein and small-molecule-ligand heads from the same forward pass,
reported ROC AUC 0.93 PPI on n=417 and 0.86 ligand **[L-a]**) match P2Rank on
small molecules? If yes, the whole multi-modality claim rests on one validated
model and the biologics site stage becomes the *other head of a model we already
trust*.

*Method.* Run PeSTo's ligand head on the same 38 KDR and 31 CDK2 apo receptors,
score against the identical 4.5 Å contact truth, paired per structure, same Holm
correction used for the BoltzGen comparison. **Define the residues→pocket rule
before running** (PeSTo emits per-residue scores, not pockets): fix it as
"score ≥ 0.5, largest connected surface component", and pre-register that choice
so it cannot be tuned to the answer. Then run the PPI head over the 200
`epitope_gate` receptors against the measured floor (0.1042 all / 0.1836 peptide /
0.0903 PPI / 0.0255 antibody).

*Pass.* PeSTo-ligand Jaccard ≥ P2Rank's 0.635 within the paired CI on KDR **and**
≥ 0.494 on CDK2 → adopt PeSTo as the single site-definition stage.
*Fail.* Keep P2Rank for pockets, and use PeSTo only where it beats the epitope
floor. **Either outcome is reportable**, and it is a genuinely new number: **[?]**
no published head-to-head of PeSTo's ligand head against P2Rank or fpocket exists.

*Feasibility.* Confirmed today: torch 2.14.0 installed, gemmi cp314 arm64 wheel
available, paper reports ~5.3 s/structure on CPU. **Do not vendor the weights**
(CC BY-NC-SA 4.0 vs our MIT).

*Unblocks.* E2, and the entire biologics site-definition decision.

### E2 — Re-rank P2Rank's peptide pockets · **0 credits · CPU · ~half a day**

*Question.* For peptides, P2Rank offers a usable pocket and picks the wrong one
(top-1 0.1836 → oracle 0.3112, precision 0.619, a pocket touches the interface
78.6% of the time). Can a re-ranking rule close the gap?

*Method.* Score each offered pocket by candidate rules — PeSTo PPI-head score mass
inside the pocket (if E1 runs), fpocket druggability, pocket surface area, burial —
on a 50/50 split of the 70 peptide complexes, choosing the rule on the first half
and reporting on the second.

*Pass.* ≥ half the oracle gap closed on the held-out half, paired Wilcoxon.
*Fail.* Peptide site definition moves to PeSTo's PPI head or to known-site lookup.

*Why it ranks this high.* It is the single largest measured headroom in the whole
pipeline (+0.128 Jaccard available for free), and nobody has to be asked for a key.

### E3 — Prove the Jaccard artifact on real structures · **0 credits · CPU · ~2 hours**

*Question.* Does cross-modality Jaccard fall as the size ratio grows, on pairs
crystallographers describe as *the same site*?

*Method.* 10–30 same-target pairs with both a small-molecule and a
peptide/protein co-crystal. Verified starters: **CXCR4 3ODU (IT1t) vs 3OE0 (CVX15,
16-mer cyclic)** — the primary paper states the sites significantly overlap and
the peptide fills most of the pocket volume, i.e. containment, exactly the case
Jaccard breaks on; **MDM2 1YCR (p53 peptide)** vs a nutlin complex; **BCL-2 6O0K
(venetoclax)** vs a BH3-peptide complex. **[?]** The nutlin and BH3 partner PDB IDs
still need looking up. 2P2Idb supplies more pairs mechanically. Use only existing
code: `ligand_contacts`, `chain_contacts`, `buried_area_by_residue*`.

*Pass.* Jaccard correlates strongly and negatively with |B|/|A| while coverage and
min-BSA overlap do not → the metric change is justified by our own data, not by an
argument. *Fail.* A genuine surprise, and cheaper to learn now than after E7.

### E4 — Is a peptide-**drug** site more pocket-like than a generic peptide site? · **0 credits · CPU + downloads · ~half a day**

*Question.* The current peptide class (n=70, J 0.1836) is dominated by flat
coactivator-groove motifs. Class-B1 GPCR peptide drugs insert their N-terminus
into the seven-TM core — an enclosed concave volume. Is the drug-relevant peptide
site materially easier?

*Method.* `scripts/epitope_gate.py` already computes the right metric; only its
selection filters block this. Today it requires X-ray, ≤ 3.0 Å, and
**exactly 2 polymer entities** — which excludes every GPCR–peptide-drug complex,
since they are cryo-EM with Gs α/β/γ plus a nanobody (5–6 entities). Relax to
allow cryo-EM and > 2 entities while restricting the binder to the annotated
peptide entity, and point it at ~120 fixed PDB IDs (45 GLP-1R, 51 PTH1R, 13 SSTR2,
10 MC4R). At the measured 0.043 s/structure this is under 10 s of P2Rank.

*Pass.* Drug-site Jaccard significantly above the generic peptide 0.1836
(Mann-Whitney, with n). *Fail.* Peptide site definition needs replacing too, and
E7 should not be funded until it is.

*Why it matters beyond the number.* If it passes, the honest conclusion is
narrower and more defensible than "pocket finders fail on peptides": **pocket
finders fail on flat interfaces, peptide sites span the range, and `site_kind`
should be assigned by measured burial of the true interface, not by what the
binder is made of.** That reframing is worth more than another AUC.

### E5 — Fix the biologics library before any arm is built on it · **0 credits · CPU · ~2 hours**

*Current measured state* of `data/biologic_drugs.csv` (447 rows, re-counted
today — note these differ from an earlier researcher's counts, so the file has
moved): antibody 151, protein 145, peptide 96, other 55. Usable sequences:
peptide 72/96, antibody 118/151, protein 109/145, other 0/55.

Concrete defects still present: **liraglutide, octreptide and glucagon are
recorded as `modality=protein` with a blank `n_residues`**; **abaloparatide and
tirzepatide are `usable_sequence=0`** (excluded on non-standard residues — Aib);
**secretin is absent** although it is the natural Tier-2 decoy; and semaglutide is
recorded at 31 residues with `has_nonstandard_residues=0`, which **[?]** may mean
its Aib8 was silently substituted away. Several peptides reach usability only by
substitution (`residue_decision=substitute`: exenatide, setmelanotide,
pramlintide) — that is a stated approximation, not a free pass.

*Action.* Fix the modality column, resolve the semaglutide question by tracing it
to source, and write the natively-standard-residue ceiling into
`docs/03-SCOPE-AND-CONSTRAINTS.md` as a stated constraint up front. Declare
**lipidation and D-amino acids out of scope in writing**: AF3-family models do not
support covalent bonds into polymer entities, and D-residue chirality violation is
reported around 50%, i.e. a coin flip per residue **[L-b]**. If semaglutide is
modelled at all, model the unlipidated backbone and say so on the row — never
silently.

### E6 — The two-minute probe that gates every biologics gate · **~10–30 credits · ~15 min**

*Question, three parts.* (1) Do protein-chain co-folds return a **PAE matrix**?
(2) What does a peptide co-fold actually cost? (3) Does Rowan's Boltz-2 endpoint
accept three chains (VH + VL + antigen), or require an scFv linker?

*Method.* Submit **one** GLP-1R + exenatide co-fold and **one** Fv + antigen
co-fold. Record credits, wall seconds, the full returned score object, and whether
any artifact contains PAE.

*Pass/fail.* If PAE is returned, ipSAE / Pinc / LIS / MiniPAE all become available
and the biologics confidence gate is a solved problem. **If it is not — which is
what the four small-molecule runs on disk suggest — then pDockQ (interface pLDDT +
contact count, both derivable from the pose) is the only interface-aware gate
available to us**, and that constraint goes into `docs/03` in writing rather than
being discovered at 3am.

*Unblocks.* E7's budget, and every cost estimate below. **Run this before
committing to any multi-job biologics spend** — our 5.4 cr / 4 min figure is a
small-molecule number and will not transfer.

### E7 — The GLP-1R peptide specificity test · **~150–400 credits · half a day**

*Question.* The thing nobody has answered for us: can interface overlap separate a
true peptide binder from a **class-matched** decoy?

*Method.* GLP-1R. 7 approved agonists (3 with solved complexes to define the site,
4 as the zero-structure holdout) + ~12 decoys across the three tiers of §3.2. One
seed each — the sampling budget belongs on decoys, not seeds (1 seed gives 70.4%
DockQ ≥ 0.49 vs 73.6% for 20 seeds **[L-a]**, so 20× the credits buys 3 points).
Compute interface overlap against the crystal site, record ipTM, apply the ≥ 0.8
gate, and report the abstention rate as a first-class number.

*Pass criterion, stated before the run.* Overlap must separate positives from
**Tier-2 class-B siblings**, not merely from shuffled or unrelated peptides.
Report AUC with n and a permutation p either way.
*Fail.* The peptide arm has the same ceiling as the small-molecule arm, it goes in
the results table next to 0.717/0.696, and we say so rather than shipping it as a
discovery engine.

*Cost.* 19 jobs. At 3 cr/min a 40-mer peptide co-fold is plausibly 3–5× a small
molecule → budget 300–500 and let E6's measurement replace this estimate.

<!--LANG-EXEMPT-->
### E8 — The affinity-head arm on KDR · **~220 credits · needs a §4.4 ruling first**

*Question.* Does Boltz-2's `affinity_probability_binary` — the head actually
trained for binder-vs-decoy (MF-PCBA AUROC 0.812, AP 0.0248, EF@0.5% 18.4, vs
ipTM 0.566 **[L-a]**) — separate the **hard** kinase-inhibitor null where
interface overlap (0.717) and ipTM (0.696) both fail?

*Blocker, not a technicality.* `PROJECT_GOAL.md` §4.4 removes affinity from the
ranking path and `repurpose.py` enforces it. Using the head as a **gate** is
arguably compatible; using it as a **ranker** is not. **Evan decides, in writing,
before this runs.**

*Method.* Re-run the same 41 drugs against KDR with affinity enabled,
unconstrained. One table, n=41, four columns — interface overlap, ipTM,
ligand_ipTM, `affinity_probability_binary` — each with AUC against both nulls and
a significance test. Aspirin gets its own labelled row.

*Pass.* Affinity head separates the hard null → it becomes the gate and interface
overlap becomes the interpretable ranker downstream.
*Fail.* We have measured that **no available score discriminates chemically
related actives on this target** — which matches the published dopamine-D4 result
(AF3 ligand-pLDDT AUC 0.46 **[L-a]**) and is a publishable finding, not a gap.
<!--/LANG-EXEMPT-->

### E9 — Only if E1 fails: an actual epitope predictor · **Modal, GPU, blocked**

Options in order: ScanNet `--mode epitope` (Apache-2.0, weights included, AUCPR
0.178 on B-cell epitopes **[L-a]** — note how low that is in absolute terms);
DiscoTope-3.0 (AUC-PR 0.232 on a 24-antigen external set **[L-a]**, **[?]** licence
likely academic-only); dMaSIF. **ScanNet cannot run on this Mac** — it pins
`tensorflow==1.14.0` / `keras==2.2.5` / `numpy==1.19.5`, which has no arm64 or
py3.14 path at all, so it needs a linux/amd64 py3.7 container. **This is the one
place Modal is genuinely justified** — and it needs a payment method added, since
all four GPU tiers are currently refused **[M]**. Given AsEP's ceiling, I would not
spend the hackathon on it.

---

## 5. What needs Modal rather than Rowan (i.e. what a key would be for)

| Method | Why Rowan cannot | Worth asking for? |
|---|---|---|
| **ScanNet** (epitope mode) | TF 1.14 / py3.7 / linux-amd64 container | Only if E1 fails. Note the account also needs a **payment method**, not just a key — every GPU tier is currently refused **[M]**. |
| AlphaFold2-Multimer / ColabFold | Rowan hosts OpenFold-3, not AF2-Multimer **[?]** — confirm | Interesting: AF2 is the surprise winner on leakage-controlled peptide sets (94% correct mode on 67 post-cutoff GPCR–peptide complexes, beating AF3 **[L-a]**). Second priority. |
| CyclicBoltz1 / HighFold3 / HighFold-MeD2 | research code, GPU | **No.** Cyclic/D-residue chemistry is a second project; declare it out of scope instead. |
| dMaSIF, DeeplyTough, PocketMiner | GPU-preferred, some need cmake | Low priority. PocketMiner matters only if we publish "flat, biologics only" and need to rule out a cryptic pocket. |
| MM-GBSA / funnel metadynamics / ABFE | CPU-hours to > 20 GPU-hours per complex | **No.** Goes in *not evaluated / why*, cited as the reason confidence cannot be trusted. |

**Recommended ask:** nothing yet. E0–E5 are free and E6–E8 are Rowan. Revisit
after E1.

---

## 6. What the literature does not settle — where our run is the contribution

These are the places where a search found no answer, so our number would be new
rather than a reproduction. They are the submission's technical-complexity story.

1. **No published head-to-head of a pocket finder on PPI epitopes vs
   small-molecule sites.** We have it: n=200, four arms, a random-surface floor,
   Holm-corrected, 8.6 s of compute, 0 credits (`results/epitope_gate.json`). The
   concavity gradient (peptide 0.1836 > PPI 0.0903 > antibody 0.0255 ≈ random) is,
   as far as six searches could tell, unpublished.
2. **No comparison of PeSTo's ligand head against P2Rank or fpocket.** E1 produces
   it.
3. **No benchmark of interface-overlap ranking against within-family hard
   negatives.** Our 0.717 (overlap) and 0.696 (ipTM) on other kinase inhibitors
   have nothing to be compared against — which also means we cannot say whether
   they are good or bad for the method class. Genuinely unsettled; say so.
4. **No protein–peptide DockQ benchmark for Boltz-2 anywhere.** Its own paper uses
   DockQ only for antibodies. E7 gives a target-specific datapoint.
5. **Whether Rowan exports PAE for a protein–protein co-fold.** E6, 15 minutes.
6. **Whether a generic interface predictor transfers to antibody epitopes.**
   MaSIF-site scores ROC AUC ~0.84 on generic interfaces and MCC 0.037 on AsEP
   **[L-a]** — so PeSTo's 0.93 should **not** be assumed to transfer, and we should
   not claim it does without measuring it.

---

## 7. Honest ledger — our own results that cut against the project

Kept here so they cannot be quietly dropped, per the working agreement.

- **Our design-based signature lost to a pocket finder.** BoltzGen design consensus
  vs P2Rank: −0.324 Jaccard on KDR and −0.155 on CDK2, Holm p=0.0016 **[M**
  `m2_gate_KDR.json`, `m2_gate_CDK2.json`**]**. The v0.2 thesis in `PROJECT_GOAL.md`
  §1.1 — that the designed binder reveals the hotspot map — is **not supported by
  our own measurement**. P2Rank geometry also ties the known-ligand signature
  (p=0.35 on Jaccard), so the cheapest arm is the best arm.
- **A pocket constraint manufactures confidence.** Constrained aspirin: 11/17
  pocket residues, ipTM 0.965, 11.1 cr. Free aspirin: 6/17, ipTM 0.847, 5.7 cr
  **[M]**. We removed the constraint because of this, against `PROJECT_GOAL.md` F3.
- **Confidence outperformed our score on the easy null** (E0). ipTM 0.995 vs our
  0.913. The README's current framing of confidence as uninformative was drawn from
  the constrained probe and does not survive the unconstrained measurement.
- **Both scores fail on the hard null**, at 0.696 and 0.717, neither significant at
  n=41. The method's real limit is within-family discrimination, and no metric
  change in §2 is expected to rescue it.
- **P2Rank returns no pocket at all for 19.5% of biologic receptors** (32.8% of
  antibody ones) **[M]**. Those are not low scores, they are absent predictions.
- **What this method structurally cannot find.** PD-L1/BMS-202 (5J89): the small
  molecule binds a cryptic cleft at the PD-L1 **homodimer** interface and blocks
  PD-1 allosterically — it does not sit on the PD-1 epitope **[L-b]**. An
  interface-overlap ranking scores it ~0 and misses a real, clinically pursued
  mechanism. This belongs in the README's limitations section, stated up front.
- **Every receptor in the epitope benchmark is the bound conformation with the
  binder deleted** — the generous case for a pocket finder. A prospective run
  starts from apo or predicted structure, and PeSTo's own PPDB5 numbers drop from
  0.85 bound to 0.78 unbound **[L-a]**.

---

## 8. What we are allowed to claim

| Modality | Rank on it? | Basis |
|---|---|---|
| **Small molecule** | **Yes**, with both nulls reported and the easy one labelled easy | J 0.635 site recovery (n=38); AUC 0.913 easy / 0.717 hard (n=41) **[M]** |
| **Peptide** | **Not yet** — site definition and specificity are both unproven. Ship the *coverage table*, not a leaderboard, until E4 and E7 pass | P2Rank top-1 J 0.1836 with 0.3112 oracle headroom (n=70) **[M]**; no Boltz-2 peptide benchmark exists **[?]** |
| **Antibody / nanobody** | **No.** Site definition by **lookup**; no ranking, no co-folding in the ranking path | P2Rank ≈ random (n=58) **[M]**; cognate-vs-non-cognate AUC 0.50–0.60 (n=519) and Boltz-1 4.1% high-accuracy poses (n=49) **[L-a]** |

*Not evaluated / why* (to be carried into the README): antibody epitope
**prediction** — best published binder-vs-non-binder ROC-AUC on therapeutic
antibodies is 0.50–0.60 (n=519, Thera-SAbDab), and confidence-based selection is
shown not to recover the good predictions. That is a finding with a citation and
an n, not a gap.
