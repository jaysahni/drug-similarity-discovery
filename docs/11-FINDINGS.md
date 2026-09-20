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
7. **The corpus was mis-described as "4,099 approved drugs" throughout.** It is
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
