# 11 — Findings: four negatives, one positive, and the controls that tell them apart

Final report for the lightweight demo. Run of 2026-09-20; targets CDK2 (6Q4G) and
thrombin (1PPB); 319 of 2,586 Rowan credits spent.

Rendered version: <https://claude.ai/artifact/AtzQNwB4fqYkZSV8otfcpq>

Every number here was computed by a script in this repository and is recorded
under `results/demo/`. Nothing is estimated.

---

## The ledger

### A — Hotspot coverage does not separate known binders from decoys
**NEGATIVE · CDK2 · n = 14 drugs** · `results/demo/cdk2/05_ranked.json`

Fourteen approved drugs co-folded against CDK2 with Boltz-2, ranked by how much
of the designed binder's hotspot set each engages.

| tier | n | mean coverage |
|---|---|---|
| positive (CDK2-annotated) | 5 | 0.507 |
| candidate (other kinase drugs) | 4 | 0.417 |
| decoy (non-kinase) | 5 | 0.427 |

Mann-Whitney positive > decoy: **U = 16.5, p = 0.229**. The top-ranked drug is
**atorvastatin**, a statin with no CDK2 relationship. The five CDK2-annotated
drugs land at ranks 2, 3, 4, 7, 9 of 14; palbociclib, an approved CDK4/6
inhibitor, ranks 9th.

Boltz-2 *pose confidence* does separate the tiers (ipTM 0.962 vs 0.878,
p = 0.0159), but that is confidence rather than hotspot coverage,
PROJECT_GOAL.md §4.4 deliberately keeps confidence out of the ranking path, and
it is the quantity most exposed to training-set leakage. Recorded as an
observation; not promoted to the ranking.

### B — Right-sizing the binder nearly triples interface confidence
**METHOD · CDK2 → thrombin** · `results/demo/*/03_designs.json`

| | protocol | binder | ipTM |
|---|---|---|---|
| CDK2, ATP pocket | `protein-anything` | 60–90 aa | 0.183–0.288, mean **0.213** |
| thrombin, protease groove | `peptide-anything` | 8–16 aa | 0.489–0.642, mean **0.572** |

A 60–90-mer cannot enter an ATP slot; the heavy pipeline's binder centroids sit a
median 26.4 Å from the pocket. An 8–16-mer fits a protease groove.

Both earlier runs had `budget = num_designs`, which switches BoltzGen's own
filter funnel off and admits every design, worst included, to the consensus. It
now runs at the documented 5:1 ratio.

Separately: PROJECT_GOAL.md §4.3's ipTM ≥ 0.85 gate has no BoltzGen provenance.
The paper publishes no ipTM distribution at all (verified, PMC12697729), and
published miniprotein work filters at > 0.5. **The plan document is wrong on this
point and should be corrected** — see `docs/07`.

### C — No designed peptide resembles an approved thrombin peptide
**NEGATIVE · thrombin · n = 10 designs** · `results/demo/thrombin/06_peptide_match.json`

Every design's nearest neighbour is **abarelix** (a GnRH antagonist) or
**afamelanotide** (a melanocortin agonist). Null percentiles span 0.223–0.983,
seven of ten inside the null of composition-preserving shuffles. ESM-C and
BLOSUM62 identity disagree on the nearest neighbour for **all ten**. Corpus: 37
approved biologics ≤ 100 aa, 2 annotated to F2, base rate 5.4%.

**Why this is readable rather than a shrug:**

- **Positive control passes.** bivalirudin → lepirudin at **rank 1 of 36**
  (cosine 0.9179); lepirudin → bivalirudin at rank 5. ESM-C cosine does retrieve
  co-target peptides in this corpus. The metric works; the designs fail it.
- **Hub check explains abarelix.** **53%** of randomly shuffled sequences land on
  it too — it attracts anything that looks like noise, so a design landing there
  carries no information about the design.

A correction worth recording: ESM-C cosine is **not** compressed into a narrow
band. Per-pair cosines span 0.211–0.985 (sd 0.242). It is the *null* that sits
near 0.97, because the null is a distribution of maxima over 37 entries. Read the
percentile, not the raw cosine.

### D — Docking does enrich for known-binder chemistry, but only with power
**POSITIVE · thrombin · n = 199 docked** · `results/demo/thrombin/09_smallmol_match.json`

Among generated molecules docked into the thrombin site, a better Vina score goes
with greater ECFP4 similarity to the drugs already known to bind thrombin. Vina
is negative-is-better, so the negative correlation points the right way.

| | ρ | p |
|---|---|---|
| n = 23 (first pass) | −0.084 | 0.703 |
| **n = 199** | **−0.146** | **0.040** |

Quartile test agrees: best 49 dockers mean 0.1835 similarity to known binders vs
worst 49 at 0.1662, Mann-Whitney **p = 0.021**.

**This result did not exist at n = 23.** Docking 200 molecules instead of 24 cost
19.23 credits and turned a reported null into a significant effect. The effect
size is small and stays small; what changed was the ability to see it. The first
pass was reported as "docking adds nothing" — that was an underpowered null
mislabelled as a negative.

### E — No generated molecule resembles a thrombin drug
**NEGATIVE · thrombin · 591 generated vs 2,153 approved** · `09_smallmol_match.json`

591 drug-like molecules sampled from `entropy/gpt2_zinc_87m` (MIT), matched by
ECFP4-2048 Tanimoto. Best match is **dithranol**, a psoriasis anthralin, at
Tanimoto 0.579 — a real chemical resemblance to an approved drug, but not to a
thrombin drug. Nothing in the top ten carries an F2 annotation. Null over the 391
undocked molecules: mean 0.270, p95 0.350. F2-annotated drugs in the corpus:
12 of 2,153, base rate 0.56%.

**What this is:** virtual screening with a generative front end. The generator is
not conditioned on the binding site; the pocket enters only as a docking filter.
Every pocket-conditioned generator (DiffSBDD, Pocket2Mol, TargetDiff, PocketFlow,
DecompDiff, LiGAN) depends on `torch-scatter`/`torch-sparse`/`torch-cluster`,
which have no build for torch 2.14 and have never shipped a native macOS arm64
wheel. Rowan hosts no de novo generator either.

---

## Second round: improving the similarity itself

Five further experiments, run after the first report, aimed at the matching step
rather than the pipeline around it.

### F — Rank fusion wins global ranking and loses the top of the list
**MIXED · n = 2,114 drugs** · `results/demo/fusion/fusion_targets.json`

Reciprocal Rank Fusion over ECFP4 + USRCAT + pharmacophore, on the E1 retrieval
task, using `scripts/metrics.py` for CIs, paired Wilcoxon and Holm correction:

| representation | AUROC | p@1 | nDCG@10 | ef@5pct |
|---|---|---|---|---|
| morgan (ECFP4) | 0.6692 | **0.6597** | **0.5122** | 5.2437 |
| usrcat | 0.6227 | 0.3504 | 0.2580 | 3.2134 |
| gobbi_pharm2d | 0.6627 | 0.6133 | 0.4648 | 4.8434 |
| **RRF fusion** | **0.6839** | 0.6240 | 0.4979 | **5.3441** |

All 12 comparisons Holm-significant. Fusion beats ECFP4 on AUROC (+0.0147,
p = 0.0024) and ef@5pct (+0.100, p = 0.0172); loses on p@1 (−0.0358) and nDCG@10
(−0.0143). README result #2 a third time — report one metric and you hide the
disagreement. USRCAT alone is catastrophic at the top (p@1 0.350 vs 0.660), which
plausibly drags the fusion there.

RRF rather than score averaging because Tanimoto and cosine are not on a common
scale, and normalising them invents a calibration nobody measured.

### G — Multi-query retrieval is worse, not better
**NEGATIVE · n = 300 drugs**

Querying by the whole co-target group instead of one molecule: p@1 −0.158
(Holm p = 0.0032), nDCG@10 −0.074 (p = 0.0056), ef@5pct −1.05 (p = 0.026), AUROC
flat. A max over the group means an item scores high if it resembles *any*
member, so one promiscuous binder drags unrelated molecules up and flattens
exactly the top of the list.

Limitation: the group is "shares any target", and the corpus maximum is 275
annotated targets per drug, so those groups are large and heterogeneous. A
per-named-target version would be a fairer test and is not what this measured.

A leak was caught mid-implementation: the first version scored each relevant item
against a group containing itself, which would have produced near-perfect
retrieval measuring nothing.

### H — The two arms are uncorrelated, and fusing them beats either alone
**POSITIVE · CDK2 · n = 14 drugs, 5 positives** · `results/demo/cdk2/10_arm_agreement.json`

| arm | mean positive | mean decoy | p | top-ranked | positive ranks |
|---|---|---|---|---|---|
| footprint | 0.507 | 0.427 | 0.2289 | atorvastatin *(decoy)* | 2, 3, 4, 7, 9 |
| chemical | 0.408 | 0.175 | **0.0297** | ribociclib | 1, 2, 3, 5, 11 |
| **fused** | — | — | **0.0159** | trilaciclib | 1, 2, 3, 5, **7** |

Kendall τ between arms = **−0.083, p = 0.695** — they order the same drugs
essentially independently, which is the precondition for combining them being
worth anything. The footprint arm is not significant alone yet still contributes.

Caveats: n = 14 with 5 positives, so read the direction not the p-value; and the
chemical arm scores against curated annotations, so it is closer to a lookup than
the footprint arm and should be expected to do better for that reason alone.

### I — Scaffold seeding: a decoy control destroys the headline
**NEGATIVE · thrombin** · `results/demo/seeded/11_seeded_vs_unseeded.json`

| set | n | best NN | median | NN is F2-annotated |
|---|---|---|---|---|
| seeded (from the 4 F2 drugs) | 378 | 0.6026 | 0.3066 | 151 (40.0%) |
| seeded, seeds deleted | 378 | 0.4182 | 0.2872 | 23 (6.1%) |
| LibInvent (REINVENT4) | 81 | 0.5789 | 0.3725 | 1 (1.2%) |
| **decoy-seeded**, F2-free drugs | 726 | **0.8200** | **0.3623** | 1 (0.14%) |
| unseeded | 591 | 0.5789 | 0.2639 | 5 (0.85%) |

Seeding beats unseeded on Tanimoto at p = 1.25e-45 — and **BRICS seeded from
random, MW-matched, thrombin-free approved drugs does better still**
(rank-biserial −0.452 against the thrombin arm). Recombining approved-drug
fragments produces approved-drug-like molecules regardless of the seeds, so
**the Tanimoto lift is 0% target-specific**. Without the decoy panel, "0.60 to
argatroban beats 0.58 to dithranol" would have read as a result.

The annotation enrichment is target-specific but ~89% circular: 135 of 151
annotated nearest neighbours are the seeds themselves (argatroban ×115,
ximelagatran ×20), and all ten of the top ten are argatroban. Leave-seeds-out
leaves 23, every one captopril — an off-target annotation. Structural caveat: all
four F2-*mechanism* drugs are the seeds, so leave-seeds-out makes a
mechanism-level hit impossible by construction.

**One non-circular hit, n = 1.** LibInvent proposed **nafamostat**, a
guanidinobenzoate serine-protease inhibitor — chemically right for an Asp189 S1
pocket, not a seed, and verified F2-annotated in the corpus (off-target). 1/81 vs
5/591, **p = 0.539**. Not significant; reported as what it is.

**Where seeding unambiguously works:** the arginine mimetic, the pharmacophore
this target needs. Unseeded proposes it 1 time in 591 (0.17%); seeded 38.6%.
Fisher OR = 371, p = 1.66e-67.

### J — Annotation enrichment helps coverage, not thrombin
**MIXED · n = 28,337 annotation rows** · `docs/12-TARGET-ANNOTATIONS.md`

Approved-library coverage 77.5% → 80.4%; 63 of 484 unannotated approved drugs
rescued. Agreement with the curated column, n = 1,167: 96.1% share ≥1 target.

But **F2 coverage does not improve**: 12 → 13, and the 13th is liothyronine, a
false positive from a `THR` abbreviation collision assigning thyroid hormone
receptor data to Prothrombin. Inspection turns up further suspects at threshold —
captopril (an ACE inhibitor) at Ki 7.32, betrixaban (a selective FXa drug) at
IC50 7.75. **The enriched set is noisier than the curated one for the one target
this project uses for controls, so it is not wired into them.**

The cause is the corpus, not the annotations: dabigatran's active form,
melagatran, lepirudin, desirudin, hirudin and heparin are all absent from
`approved_drugs.csv`. No annotation source fixes that; only adding drugs does.

Three findings that generalise: an unfiltered ChEMBL pull is actively harmful
(aspirin→TSHR; restricting to Ki/Kd/IC50/EC50 and assay type B/F moved agreement
83.1% → 96.1% with recall *up*); human-only was the wrong filter (27% of a random
120 unannotated drugs get a hit once any organism is allowed — the tail is
anti-infectives); and 212 approved drugs have no activity data anywhere, many of
them imaging agents, sunscreens and excipients that should be flagged **out of
the retrieval denominator** rather than counted as misses.

---

## Controls

A failed search and a broken instrument look identical from outside. Every arm
runs a positive control on the same data with the same metric.

| control | query → target | result | verdict |
|---|---|---|---|
| ESM-C, peptides | bivalirudin → lepirudin | rank 1 of 36 | passes |
| ESM-C, peptides | lepirudin → bivalirudin | rank 5 of 36 | passes |
| ECFP4, small molecules | argatroban → bivalirudin | rank 9 of 2,152 | passes |
| ECFP4, small molecules | ximelagatran → dabigatran | rank 21 of 2,152 | passes |
| Vina docking | 3 known drugs vs 199 generated | beat the median | passes |
| Hub check | random shuffles → abarelix | 53% of 60 | hub found |

**One control nearly reported itself as broken.** Read as a top-5 list, ECFP4
looked useless: argatroban's five nearest approved drugs are all peptidomimetics
— angiotensin II, icatibant, lisinopril — and not one is a thrombin drug. The
*rank* of the co-mechanism drugs says the opposite. Chemical similarity groups by
chemotype; the co-target drugs it finds, it finds at rank 10–82 rather than rank
1. This replicates README result #3 on a second target.

---

## Six faults that produced ordinary-looking wrong answers

None threw an error. Each is now guarded by a test.

| fault | consequence | evidence |
|---|---|---|
| Co-crystal ligand left in the pocket | designs bound elsewhere on the protein entirely | Jaccard 0.000 vs known ligand |
| PDB parser ignored insertion codes | 28 of thrombin's 259 residues vanish, including the whole 60-loop lining the active site | 259 → 231 |
| `Protein.prepare` renames chains | catalytic chain H becomes B; every later stage points at nothing | H→B, L→A |
| Design structures cached without a uuid | a re-run scored against the *previous* run's complexes | 63% identity, refused |
| Residues mapped by a constant offset | correct by accident on one run, wrong on a scattered subset | 123 of 298 |
| Percentile used a strict `<` | every exact Tanimoto tie scored against the candidate — all results biased low, invisibly | KS p < 0.001 |

**The last was caught by its own test.** Held-out queries from the null must
score uniform percentiles; KS said D = 0.187, mean 0.457 instead of 0.500.
Tanimoto over sparse fingerprints is heavily discrete, so ties are common and the
strict comparison charged all of them to the candidate. Corrected to mid-rank;
both metrics now calibrate at p = 0.416 and p = 0.074.

They share one root cause: **trusting an identifier that an upstream tool is free
to change** — a chain letter, a residue number, a cache key. The fix each time
was to re-derive identity from content.

---

## Incidental: pocket detection cannot find a protease active site

Rowan returns six pockets for thrombin chain B and **not one contains Ser195 or
Asp102**. Best overlap with the curated catalytic triad is 1 of 3, with two
pockets tied.

This is a fact about the target, not a defect in the detector: a protease active
site is a shallow groove across S1/S2/S3 subsites, not an enclosed cavity, so a
cavity finder fragments it. CDK2's ATP site, a real cavity, overlapped 11 of 20.

The pipeline falls back to a geometric site anchored on the catalytic nucleophile
(His57 and Asp102 sit *behind* Ser195, so a sphere around all three reaches into
the core) and records that it did. Validation computed afterwards and never used
to choose: the site recovers 10 of PPACK's 23 contact residues, Jaccard 0.303.

---

## What this does not claim

- **Similarity is not affinity.** Nothing here predicts binding. Every output is
  a computationally ranked hypothesis requiring experimental validation.
- **n is small throughout.** Ten designs in arm P, five drugs per tier on CDK2,
  one target per arm, one structure per target. No cross-target claim.
- **Thrombin is a friendly case**, chosen *because* both arms have approved
  positive controls. A target where the answer is less knowable would test the
  method harder.
- **The generated arm was already known to lose.** `results/m2_gate_CDK2.json`
  measured it against a pocket finder over 31 held-out CDK2 co-crystals: Jaccard
  0.338 vs 0.494, Holm p = 0.0016. This run replicates that rather than
  overturning it — after the stripping fix, the design consensus reaches 0.083
  against known ligand contacts where pocket geometry reaches 0.455.
- **The biologics corpus is incomplete.** 45 approved antibodies carry a ChEMBL
  record with no sequence, so TNF coverage is 1 of 5. Listed, not hidden.

---

## Retrospective

### What went wrong in how this was built

1. **The wrong pipeline was built first.** Footprint matching (generated binder as
   a discarded *probe*) was built before establishing that the intent was direct
   matching (generated binder as the *query*). One question at the outset — probe
   or query? — would have saved roughly half the working time.
2. **A specified step was skipped.** PROJECT_GOAL §3 says strip ligands. It was
   not done, so BoltzGen designed around an occupied ATP pocket: 52 credits and a
   full run, with designs that bound somewhere else entirely.
3. **A number was reported wrong.** Our ipTM 0.16 was compared against the heavy
   run's *top few* designs rather than its median (0.1685). The two runs are
   statistically indistinguishable (Mann-Whitney p = 0.14).
4. **An underpowered null was reported as a negative.** Finding D at n = 23 was
   written up as "docking adds nothing." At n = 199 it is significant. The result
   existed the whole time.
5. **A self-referential null.** Arm S's first cut scored molecules against a null
   built from those same molecules, which only recovers each molecule's rank
   within its own set.
6. **A broad `git add`** swept in another agent's mid-edit file, committing a
   stale intermediate version.
7. **The site definition omitted the pocket's defining residue.** The thrombin
   site was built as a 6 A shell around catalytic Ser195, and **Asp189 is not in
   it** — UniProt 562, structure label 199, the base of the S1 pocket and the
   salt-bridge partner for essentially every thrombin inhibitor. The anchor
   choice was justified in writing; what was never checked is whether the
   resulting shell contained S1. The peptide designs were aimed at a site
   missing it.
8. **"48% of the library is unannotated" was over the wrong denominator.** 1,501
   of the 1,985 unannotated rows are `approved != 1`. The approved library is
   **77.5%** annotated, so the gap called "the binding constraint" is about half
   the size reported.
9. **The corpus was mis-described as "4,099 approved drugs" throughout.** It is
   4,099 DrugCentral *structures*, of which only **2,153** carry `approved == 1`;
   the rest have an empty `approval_agencies`. The matching code did not filter,
   so two reported top hits (idrocilamide, pamaquine) were not approved drugs at
   all, and the target base rate was understated by a diluted denominator.
   Corrected: the effect in finding D survives but weakens, from ρ = −0.182,
   p = 0.0099 to ρ = −0.146, p = 0.040.

### Ideas, ranked by value per hour

1. **Run a power check before reporting any null.** Compute the n needed to detect
   a small effect; either reach it, or label the result *underpowered* rather than
   *negative*. Docking costs ~0.1 credit per molecule — power was affordable and
   was not bought.
2. **Get a genuinely pocket-conditioned generator.** DiffSBDD is reachable: MIT,
   Zenodo checkpoints, and its only blocker (`torch-scatter`) needs a ~10-line
   shim. Or REINVENT4 with Rowan docking as the RL reward.
3. **Make content-derived identity a stated invariant with a test.** All six
   faults above share one root cause; they were fixed six separate times.
4. **Build the golden fixtures** PROJECT_GOAL §5 rule 4 already mandates. Their
   absence is why a statistics block vanished the first time a stage re-ran.
5. **Compare the two arms.** Two independent rankings of approved drugs for the
   same target now exist and have never been compared. Where they agree is more
   interesting than either alone, and the data is already on disk.
6. **Replicate the direct arms on a second target.** The M2 gate was not trusted
   until it replicated; the new arms have n = 1 target.
7. **Correct PROJECT_GOAL §4.3.** The document specifies a filter this repo's own
   investigation has shown to be unreachable by construction.
