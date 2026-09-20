# 05 — E2 design critique: what could make the held-out-ligand result wrong

**Status:** adversarial design review, written *before* the experiment reports a number.
**Scope:** the E2 experiment — build a consensus interface signature from all-but-one
ligand of a target, predict the held-out ligand's contact residues, compare against
P2Rank's pocket on the same event.
**Author's role:** critic only. Nothing here implements E2.

Every claim below is labelled **CONFIRMED** (measured on this repo's data, number given)
or **HYPOTHETICAL** (argued, not yet measured). Measurement provenance is in
[§10](#10-how-the-numbers-here-were-produced).

---

## 0. The one-paragraph verdict

I ran a 6-target, 173-trial pilot of the exact experiment being implemented. **As
specified — consensus core at `core_threshold=0.6`, the PROJECT_GOAL 4.3 default — the
result is null**: consensus vs P2Rank is Δ weighted-Jaccard **+0.007, 95% CI
[−0.049, +0.059], p = 0.81, n = 161**. Change one representational choice — truncate the
consensus to the held-out contact set's size instead of thresholding at 0.6 — and the
same data gives **Δ +0.138, CI [+0.100, +0.173], p = 0.0002**. The headline answer to the
decisive question is therefore currently **determined by a free parameter nobody has
registered**, not by information content. That is finding C1 and it outranks everything
else in this document. Worse, a null that costs nothing — *pick one other ligand of this
target at random* — scores median J = 0.469, and the consensus core does **not** beat it
(Δ −0.026, p = 0.10). Fix C1–C3 before any number leaves this repo.

---

## 1. Ranked findings

| # | Finding | Status | Severity | Kills the result? |
|---|---|---|---|---|
| **C1** | Answer flips sign-of-conclusion on the consensus→set rule | **CONFIRMED** | **Fatal** | Yes |
| **C2** | The n=1 null is not beaten; "aggregation" premise unsupported | **CONFIRMED** | **Fatal** | Yes |
| **C3** | This is I6.2 (the ceiling), not the I6.1 M2 gate | **CONFIRMED** | **Fatal to interpretation** | Yes, as stated |
| **C4** | Trial-level pooling violates independence; CIs ~3.4× too narrow | **CONFIRMED** | High | Inflates significance |
| **C5** | Size asymmetry caps P2Rank's Jaccard at 0.80 before any biology | **CONFIRMED** | High | Biases the comparison |
| **C6** | Lipids/detergents/mercurials pollute the held-out trial set | **CONFIRMED** | High | Biases toward consensus |
| **C7** | 9.8% of trials leave an *identical* ligand in the consensus | **CONFIRMED** | High | Direct leakage |
| **C8** | Insertion-code collapse corrupts residue ids in serine proteases | **CONFIRMED** | High | Manufactures a top score |
| **C9** | Author ids pool across chains in 15/60 targets | **CONFIRMED** | Medium | Silent id collision |
| **C10** | 4.5 Å cutoff churns ~18% of the ground truth | **CONFIRMED** | Medium | Ground-truth noise |
| **C11** | Near-duplicate chemistry — *weaker* than the brief assumed | **CONFIRMED** | Low-Medium | Hedge, don't panic |
| **C12** | Resolution confound, 3× spread in median target resolution | **CONFIRMED** | Medium | Shared confound |
| **C13** | Greedy diversity selection makes the set unrepresentative | **CONFIRMED** | Medium | External validity |
| **C14** | 32 structures on disk are not in the manifest | **CONFIRMED** | Low | Join hazard |
| **C15** | Non-primary ligand copies silently dropped | **CONFIRMED** | Low | Small bias |

---

## C1 — The conclusion is a free parameter. **CONFIRMED. Fatal.**

The consensus is a *frequency map over residues*. Turning it into a predicted set requires
a rule, and the two obvious rules give opposite answers to the M2-style question on
identical data.

Pilot: 6 targets (ESR1, BRD4, CDK2, F2, CA2, MAPK14), 173 leave-one-ligand-out trials,
161 paired with P2Rank. Metric = Jaccard against the held-out ligand's observed 4.5 Å
contact residues. Tests via `scripts/metrics.py`.

| Predictor | median J | mean J | Δ vs P2Rank (paired) | 95% CI | p |
|---|---|---|---|---|---|
| **Null: one random other ligand** | 0.469 | 0.494 | — | — | — |
| **P2Rank top pocket** | 0.533 | 0.461 | — | — | — |
| **Consensus core, freq ≥ 0.6** | 0.545 | 0.468 | **+0.007** | [−0.049, +0.059] | **0.81** |
| **Consensus top-\|O\| (size-matched)** | 0.692 | 0.599 | **+0.138** | [+0.100, +0.173] | **0.0002** |

Both rows are "the consensus interface signature". One says the generated-signature idea
adds nothing over a 1.7-second pocket finder; the other says it adds a lot. The
difference is entirely that `core_threshold=0.6` emits a **median of 10 residues** against
a **median 15-residue** truth — it under-predicts by construction and pays recall
(0.625) for precision (0.889).

**Severity: fatal.** Section 7 says the M2 answer decides whether the project ships the
design stack or the cheap pocket pipeline. That decision currently rests on an
unregistered default.

**Mitigation (required):**
1. **Pre-register the primary rule before running.** I recommend the size-matched
   top-|O| truncation as primary — it removes the size confound of C5 symmetrically — and
   `core@0.6` as the named secondary, because 4.3 defines it.
2. **Report both, always, in the same table.** Dropping the null one would violate
   "keep the negative results".
3. Report a **threshold sweep** (core ∈ {0.3, 0.4, 0.5, 0.6, 0.7, 0.8}) as a sensitivity
   curve so the reader sees the conclusion's dependence rather than inheriting a choice.

---

## C2 — The null that matters is not "random surface residues". **CONFIRMED. Fatal.**

The brief asks whether random-surface is too weak a floor. It is far too weak, but the
important floor is cheaper still and more damning:

> **Null-1: predict the held-out ligand's contacts using the contacts of a single other
> ligand of the same target, chosen at random.**

Measured (40 resamples/trial, seed 0): **median J = 0.469, mean 0.494**, n = 161.

Paired against it:

| Arm | Δ vs Null-1 | 95% CI | p |
|---|---|---|---|
| Consensus core@0.6 | **−0.026** | [−0.057, +0.005] | 0.10 |
| Consensus top-\|O\| | +0.105 | [+0.074, +0.136] | 0.0002 |
| P2Rank top pocket | −0.033 (derived) | — | — |

The as-specified consensus is **nominally worse than a single arbitrary crystal
structure**, and cannot be distinguished from it. PROJECT_GOAL 1.2's central cost argument
— that aggregating over *thousands* of designs buys a converged contact map — receives
**no support** from this arm. The reason is mechanical: distinct ligands of one target
already share most of their contacts (**median pairwise contact-set Jaccard = 0.542**,
n = 3,875 within-target ligand pairs across the 6 pilot targets). Most of what any arm
"predicts" is simply *the pocket is the pocket*, which is exactly the information P2Rank
returns for free.

**Severity: fatal to the headline.** A result reported only against a surface-random floor
would look impressive and mean nothing.

**Mitigation (required):** report a null ladder, every rung with its own n and paired
test, ordered from weakest to strongest:

1. Random surface residues, size-matched to |O| — the trivial floor (expected J ≪ 0.1).
2. The k most-buried residues of the chain, size-matched — the "any concave spot" floor.
3. **Null-1 above (one random other ligand)** — the floor that actually binds.
4. **Null-1b: the single other ligand most chemically dissimilar** to the held-out one —
   the strictly hardest version, and the one that speaks to the cross-chemistry claim.

An arm that does not clear rung 3 with a CI excluding zero has not demonstrated anything,
and the README must say so in those words.

---

## C3 — This experiment cannot answer the M2 gate. **CONFIRMED. Fatal to interpretation.**

PROJECT_GOAL 8.3 reads:

> `source="boltzgen_consensus"` vs `source="p2rank_geometry"` … **`source="known_ligand"`
> gives the upper bound.**

and WS-I lists them as *different ablations*: I6.1 is BoltzGen-vs-P2Rank; I6.2 is
signature-vs-`known_ligand`, "the ceiling". The task brief proposes running
`known_ligand` vs P2Rank and treating the outcome as the **I6.1** answer. It is not.
It is **I6.2 with the arm of interest deleted**, and the two have different logic:

| | Question | What a win means | What a loss means |
|---|---|---|---|
| **I6.1** (the gate) | Does BoltzGen beat P2Rank? | Build the design stack | Ship the pocket pipeline (§7) |
| **E2 as designed** | Does a crystal consensus beat P2Rank? | The *ceiling* is above P2Rank — BoltzGen has headroom to compete for. Says nothing about whether BoltzGen reaches it. | **The ceiling itself is at or below P2Rank** — so no generator can win, and §7's fallback is forced *a fortiori*. |

Note the asymmetry, because it is the single most useful thing this experiment can do:
**a negative E2 is far more informative than a positive one.** A null result here
(which C1 shows is the current as-specified outcome) genuinely forecloses the BoltzGen
arm without ever running BoltzGen, because you cannot beat P2Rank through a ceiling that
does not itself beat P2Rank. A positive E2 licenses almost nothing about BoltzGen.

**Mitigation (required):** state verbatim in the README and in the results artifact:

> This is ablation I6.2 (`known_ligand` ceiling) run against `p2rank_geometry`, not
> ablation I6.1. The `boltzgen_consensus` arm was not run — Modal authenticates but every
> GPU tier (T4, L4, A10G, A100) is refused pending a payment method. **The M2 gate is not
> passed and not failed by this result.** What it bounds is the ceiling the absent arm
> would have had to reach.

Add BoltzGen to the "not evaluated / why" table (§9), not to the results table.

---

## C4 — Per-trial pooling overstates significance. **CONFIRMED. High.**

Trials within a target share a protein, a pocket, a numbering scheme and a crystal form.
They are not independent draws. Measured effect on the interval, same data, same
estimator:

| Comparison | Trial-level CI (n=161) | Target-level CI (n=6) | Width ratio |
|---|---|---|---|
| core vs P2Rank | [−0.049, +0.059] | [−0.135, +0.227] | **3.4×** |
| top-\|O\| vs P2Rank | [+0.100, +0.173] | [+0.072, +0.236] | 2.2× |

The trial-level interval is roughly a third of its honest width. Note the effect survives
clustering for top-|O| (CI still excludes 0) but the *precision* claim does not.

On the imbalance the brief expected: **it is not there, and this is worth recording as a
corrected assumption.** The greedy ≤25-entry cap equalised the targets — distinct ligands
per target run **min 25, median 26, max 31** across all 60 targets, top-5 share of all
1,566 ligands = 9.4%. Zero targets have fewer than 5 ligands. **The damaging imbalance is
in per-target *difficulty*, not per-target *count*:** median J for the core arm ranges
from **0.000 (CA2)** to **0.895 (F2)** — and CA2 alone supplies 41 of 173 pilot trials
(23.7%). Per-trial means are therefore dominated by whichever hard target happens to
carry the most structures.

**Mitigation (required):**
- Primary estimand = **mean of per-target means**, with `paired_bootstrap` resampling
  **targets**, not trials. Report n_targets alongside n_trials in every headline.
- Publish the **per-target table in full** (all 60 rows, with n). The sign of the
  comparison flips across targets — in the pilot P2Rank wins CA2 (0.500 vs 0.000) and
  MAPK14 (0.621 vs 0.417) while the consensus wins BRD4 (0.867 vs 0.556) and F2 (0.895 vs
  0.604). A single pooled number hides a genuinely split result, and hiding it would
  breach the "keep the negative results" rule.
- `wilcoxon` at the target level too; `holm_bonferroni` across the arm family.

---

## C5 — P2Rank is strawmanned by object size. **CONFIRMED. High.**

P2Rank returns a *pocket*; the ground truth is a *contact set*. They are different-sized
objects, so Jaccard penalises P2Rank for reasons unrelated to information.

Measured over all 1,531 structures (final `results/p2rank_pockets.json`, n_parsed = 1531,
n_failed = 0) and the 161 paired pilot trials:

- P2Rank top pocket: **median 25 residues** (mean 28.0, p10 16, p90 43, max 106),
  median **5 pockets offered** per structure.
- Observed contact sets: **median 15 residues** (mean 14.8, min 3, max 28, n = 173).
- Paired ratio |P|/|O|: **median 1.23**; P2Rank's pocket is **larger in 77%** of trials.
- **Maximum Jaccard achievable from size alone**, assuming perfect nesting:
  **median 0.800, mean 0.731.** P2Rank's observed 0.533 is 67% of its own structural
  ceiling.
- Decomposed: P2Rank **recall 0.765**, **precision 0.647**. It is finding the right
  residues and being charged for volunteering extra ones.

Two further fairness points, both measured:

- P2Rank is scored on its **rank-1** pocket while the consensus gets the whole target's
  ligand history. Giving P2Rank its best-of-all-pockets oracle moves it only 0.533 → 0.550
  (best pocket *was* rank 1 in 94% of trials), so rank-1 is a fair choice — record this,
  as it pre-empts the obvious objection.
- The top pocket spans **>1 chain in 201 of 1,531 structures** (180 span 2, 11 span 3,
  2 span 4), so on multimers P2Rank is sometimes describing an assembly-wide site while
  the ligand sits in one chain.

**Mitigation (required):** the primary comparison must be **size-matched in both
directions** — truncate the consensus to |O| (C1) *and* truncate P2Rank's pocket to its
|O| highest-scoring residues, using the per-residue scores in each run's
`*_residues.csv`. Report precision, recall and Jaccard **separately** — a single Jaccard
conflates a size artefact with an information difference. Also report the size-free
comparison: **rank-based AUC / average precision** over each method's residue ranking
against the binary contact label, which is invariant to where you cut.

---

## C6 — A fifth of the held-out "drug" trials are not drugs. **CONFIRMED. High.**

`drug_like_ligands` keeps anything with ≥10 heavy atoms that is not on `EXCLUDED_HET`.
Applying the real exclusion list to the manifest's 1,980 (entry, ligand) instances:
**1,649 survive, 331 are excluded.** Of the survivors, these are not drug-like candidates:

| Comp id | Instances | What it is |
|---|---|---|
| CLR | 25 | cholesterol |
| OLA | 25 | oleic acid |
| OLC | 20 | monoolein (lipidic cubic phase) |
| B7G | 19 | heptyl glucoside detergent |
| TEP | 13 | theophylline |
| MBO | 8 | **organomercurial** (phasing heavy atom) |
| MMC | 6 | **mercury compound** |
| OLB / C15 / PLM / A6L | 10 | lipids |

**119 of 1,649 instances (7.2%) are lipid or detergent**, spread over 6 targets, and the
concentration is extreme where it occurs:

| Target | Surviving instances | Lipid/detergent | Share |
|---|---|---|---|
| **ADORA2A** | 113 | **85** | **75.2%** |
| PPARD | 45 | 19 | 42.2% |
| CA2 | 41 | 8 (+14 mercurials) | 19.5% (53.7% w/ Hg) |

This is not cosmetic, and it biases **toward the consensus arm**:

- Membrane lipids bind the **lipid-facing outer surface** of a GPCR, not the orthosteric
  pocket. P2Rank returns one concave pocket and will miss them entirely. The consensus,
  built from 24 other cholesterol/monoolein copies, will nail them. **ADORA2A alone would
  contribute 113 of ~1,649 trials (6.9%), 75% of them rigged in the consensus's favour.**
- CA2's mercurials bind surface cysteines. This is visible in the pilot as the reason
  CA2's core arm scores **J = 0.000** (no residue reaches frequency 0.6 because the
  "ligands" occupy unrelated sites) and its within-target pairwise contact Jaccard is
  **0.048** against 0.542 overall.

**Mitigation (required):**
1. Extend `EXCLUDED_HET` with a lipid/detergent/heavy-atom class before building trials:
   `CLR OLA OLC OLB PLM A6L MYR STE C15 B7G LDA D12 DAO LMT DMU CPS CHD MBO MMC HGB` at
   minimum. Do this in `interfaces.py` so both arms see the same set, and record the
   addition in a doc — it changes a published filter.
2. **Report the count of dropped instances and the per-target breakdown**, don't just
   drop them silently ("verify absence programmatically").
3. Cross-check survivors against DrugCentral/ChEMBL approved-drug identity where possible,
   and put anything unresolvable in the "not evaluated" table.
4. Run the headline **both with and without ADORA2A**. If the sign changes, the finding is
   about cholesterol, not about interface signatures.

---

## C7 — One trial in ten leaves an identical ligand in the consensus. **CONFIRMED. High.**

If trials are enumerated at the (entry, ligand) level — the natural reading of "all but
one ligand" — then holding out one copy of a ligand that was solved several times against
the same target leaves its twins in the training half.

Measured over the 1,649 surviving instances:

- **162 trials (9.8%)** have the held-out comp id **still present** in the consensus.
- Distribution of identical copies remaining: 0 → 1,487 trials; 1 → 32; 2 → 15;
  **≥5 → 115 trials**.

Those 115 are not predictions, they are lookups: 24 other cholesterols determine the
answer exactly. And this is the *floor*, because it counts only exact comp-id identity —
it misses the same molecule deposited under two component ids, and the same ligand in two
crystal forms.

**Mitigation (required):**
1. **Deduplicate to the (target, comp_id) level**: one trial per distinct ligand, its
   contacts averaged over copies, and *all* copies removed from the consensus when it is
   held out. This is the single highest-value, lowest-cost fix in this document.
2. Additionally hold out by **chemical cluster**, not by molecule — see C11.
3. Report the exact-duplicate rate (9.8%) as a measured property of the set.

---

## C8 — Insertion-code collapse corrupts ids where the consensus scores best. **CONFIRMED. High.**

`interfaces.py` documents that insertion codes are dropped from residue ids and warned
about. In the serine proteases, which use chymotrypsin numbering, this is not a corner
case. Sampling one entry per target across all 60:

| Target | Entry | Polymer residues with an insertion code | Within-chain id collisions |
|---|---|---|---|
| **F2** (thrombin) | 5AFY | **37 / 291 = 12.7%** | 37 |
| **PLAU** | 5YC6 | 19 / 246 = 7.7% | 19 |
| **F10** | 9I24 | 7 / 296 = 2.4% | 7 |

3 of 60 targets affected. In the pilot, **20 of 173 instances** raised
`residue-id collision across insertion codes for [60]` — the thrombin 60-loop
(60A…60I) collapsing to a single id `60`.

Why this is worse than a rounding error: **F2 is the target where the consensus posts its
highest score (J = 0.895)**. Collapsing 9 distinct 60-loop residues into one id that every
ligand touches inflates agreement for *both* arms, but it inflates the frequency-weighted
consensus more, because a collapsed id reaches frequency 1.0 trivially. One of the two
targets carrying the positive result is the one with corrupted ids.

**Mitigation (required):**
- Make the residue key `(chain_id, auth_seq_id, insertion_code)` throughout, or
- exclude F2/PLAU/F10 from the primary analysis and list them in "not evaluated / why",
  and
- either way, **re-run the per-target table and check whether the headline survives their
  removal.** Report it.

---

## C9 — Author ids pool across chains. **CONFIRMED. Medium.**

Contact sets are sets of author sequence numbers with no chain component, so in a
homodimer residue 350 of chain A and residue 350 of chain B are the same element.

- **15 of 60 targets** have >1 polymer chain in their sampled entry (NQO2 has 4; BACE1,
  PTPN1, TTR, AURKA, F10, PDE4, IRAK4, DPP4, FGFR1, PPARD, PIK3CA, MAP2K1, CA1 have 2).
- **22 of 173 pilot instances** warned `contacts span polymer chains ['A','B']`.

Effect is two-sided and unquantified: pooling can inflate agreement (a ligand in chain B
appears to match a chain-A consensus) or deflate it (a genuine cross-chain interface
looks like noise). P2Rank's residues *do* carry a chain, so the two arms are keyed
differently — which is itself an apples-to-oranges hazard.

**Mitigation:** key both arms on `(chain, resnum)`; where a target is a symmetric
homomultimer, canonicalise to one protomer explicitly and say so.

---

## C10 — The ground truth moves with the cutoff. **CONFIRMED. Medium.**

Contacts recomputed at 4.0 / 4.5 / 5.0 Å over 173 instances:

| Target | mean \|O\| @4.0 | @4.5 | @5.0 | J(4.0, 4.5) | J(4.5, 5.0) |
|---|---|---|---|---|---|
| ESR1 | 9.5 | 11.5 | 12.8 | 0.835 | 0.903 |
| BRD4 | 10.6 | 13.6 | 14.2 | 0.788 | 0.958 |
| CDK2 | 12.6 | 15.4 | 16.6 | 0.817 | 0.930 |
| F2 | 15.9 | 18.1 | 20.0 | 0.882 | 0.910 |
| CA2 | 10.0 | 12.1 | 13.9 | 0.843 | 0.851 |
| MAPK14 | 15.9 | 19.4 | 21.5 | 0.813 | 0.902 |

A 0.5 Å change churns **~18% of the ground-truth set** (J ≈ 0.79–0.88). That is the same
order as the entire top-|O| effect (+0.138) and an order larger than the core effect
(+0.007). The cutoff is therefore not a detail; **at the core threshold the effect is
smaller than the measurement's own definitional noise.**

**Mitigation (required):** run the headline at 4.0 / 4.5 / 5.0 and report the Δ at each.
A conclusion that does not hold across all three is not a conclusion. Also report, per
CLAUDE.md, whether altloc and occupancy handling is uniform — `ligand_contacts` should be
checked for whether it takes altloc A only, and the answer stated.

---

## C11 — Near-duplicate chemistry: the brief's fear is *not* borne out. **CONFIRMED. Low-Medium.**

The brief expects a target's ligands to be "mostly one chemical series". Measured — ECFP4
(Morgan r=2, 2048 bits), 12 targets, 310 ligands fingerprinted from the CCD (4 failed:
`25X 34B 39B HGB`), 3,875 within-target pairs:

- **Median pairwise Tanimoto = 0.157.**
- **Fraction of pairs ≥ 0.9: 0.0018 (0.18%).** ≥0.7: 1.5%. ≥0.5: 5.2%.
- Per-ligand nearest neighbour within its own target: **3.5% have a NN ≥ 0.9**, 22.6%
  ≥ 0.7, 45.8% ≥ 0.5. Median NN Tanimoto 0.30–0.73 by target.

| Target | n | median pairwise T | median max-T to any sibling | n with NN ≥ 0.7 |
|---|---|---|---|---|
| F2 | 24 | 0.236 | 0.732 | 14 |
| ESR1 | 26 | 0.165 | 0.657 | 11 |
| DHFR | 23 | 0.231 | 0.645 | 6 |
| PPARG | 25 | 0.162 | 0.562 | 5 |
| CDK2 | 25 | 0.138 | 0.354 | 4 |
| MAPK14 | 29 | 0.146 | 0.304 | 3 |

The chemical-series worry is **real but bounded**, and concentrated: F2 has 14 of 24
ligands with a ≥0.7 sibling. The reason the rate is low is C13 — the greedy selection
explicitly maximised distinct ligands.

**But chemical dissimilarity does not buy site dissimilarity, and that is the point that
matters.** Median within-target *contact-set* Jaccard is **0.542** while median
within-target *chemical* Tanimoto is **0.157**. Ligands that share almost no chemistry
still share half their contact residues, because they occupy the same pocket. So the
claim "predicts a held-out ligand's contacts across chemotypes" is largely **"recognises
the pocket"** — which is precisely P2Rank's job description, and precisely why C1/C2 come
out near-null.

**Mitigation:** (a) hold out by **Butina cluster at Tanimoto 0.5**, not by molecule, and
report n_clusters as the real n; (b) report the headline **stratified by the held-out
ligand's max Tanimoto to the training half** (<0.3 / 0.3–0.5 / ≥0.5) — the <0.3 stratum
is the only one that supports a cross-chemistry claim; (c) never describe the result as
cross-chemotype generalisation without that stratum's number and its n.

---

## C12 — Resolution confound. **CONFIRMED. Medium.**

- Across all 1,500 manifest entries: median resolution 1.72 Å; only 2.1% worse than 2.5 Å,
  0.3% worse than 2.8 Å. The ≤3.0 Å filter did its job — good.
- But **median resolution per target spans 0.93 Å (AKR1B1) to 2.74 Å (KDM1A)**, a ~3×
  range, and it is confounded with target identity.
- Within the pilot (a high-resolution slice, 0.90–1.80 Å), **Pearson r(resolution,
  |O|) = +0.551**, n = 173. The size of the ground truth itself tracks resolution.

There is no *between-arm* confound in the usual sense — both arms read the same
structures — but there is a **between-target** one: the per-target effect will correlate
with crystallographic quality, and since C4 makes the estimand a per-target mean, that
correlation propagates into the headline.

**Mitigation:** report the per-target effect against median resolution (scatter + Spearman
with n); if |rho| is material, report the headline additionally on the ≤2.0 Å subset.
State occupancy/altloc handling explicitly rather than assuming it is neutral.

---

## C13 — The benchmark set is deliberately unrepresentative. **CONFIRMED. Medium.**

`structures_manifest.json` selected entries greedily to **maximise distinct ligands**,
capped at ≤25 entries per target. Consequences, measured:

- ESR1: **434 distinct ligands available, 26 selected.** BRD4 410→25, BACE1 378→26,
  CDK2 378→25. Median available across the 60 targets is 79.5; median selected 26.
- This is why C11's duplicate rate is low and why C4's count imbalance vanished.

Both effects cut **against** the consensus arm (it is handed the most diverse, hardest
possible training set), so the design is *conservative* — a positive result is not an
artefact of this choice. But two honest caveats follow: (a) the measured 3.5% near-
duplicate rate is a **lower bound** on what a naive PDB-wide run would show, so it cannot
be quoted as a property of the PDB; (b) a *negative* result may partly reflect the
deliberate hardness, so a null must be reported with the selection rule stated beside it.

**Mitigation:** state the selection rule in the results artifact. Optionally run one
target (ESR1) with all 434 ligands as a sensitivity check on how much diversity
maximisation costs.

---

## C14 & C15 — Join and copy-selection hazards. **CONFIRMED. Low.**

- **C14.** `data/raw/structures/` holds **1,531** `.cif`; the manifest names **1,499**.
  The 32 extras (VEGFR2 self-test entries and similar) are *not* part of the benchmark.
  `results/p2rank_pockets.json` covers all 1,531. **Join trials on the manifest, not on
  the structures directory**, or 32 off-target structures enter the benchmark silently.
  Verified: manifest entries missing a cif = 0; cif without P2Rank = 0; P2Rank
  n_failed = 0, n_structures_with_no_pocket = **8** → those 8 need a "not evaluated" row,
  not a zero.
  *Timing note:* while this review was being written, `results/p2rank_pockets.json`
  contained only **8 of 1,531** structures (the run was mid-flight). It is now complete
  (n_parsed = 1531, 720.9 s). Any intermediate read would have silently produced an
  8-structure benchmark. **Assert `n_structures_parsed == 1531` before use.**
- **C15.** `ligand_contacts` uses one copy and warns, e.g. `6 copies of HHT; using chain A
  resseq 301, ignoring [A/302, A/303, …]`. 8 such warnings in 173 instances (~4.6%). Where
  the biologically relevant copy is not the first, the ground truth is the wrong site.
  **Mitigation:** select the copy with the largest buried area, or average over copies;
  either way report how many trials were affected.

---

## 2. What the result will and will not license

**It can say:**
- Whether a consensus built from real co-crystals predicts a held-out ligand's contact
  residues better than P2Rank's pocket, **on this 60-target set, under a stated
  consensus→set rule, at a stated cutoff, against a stated null.**
- An **upper bound on the value of any interface-signature builder** on this data —
  including BoltzGen's. This is the genuinely valuable output (C3).
- Whether the §7 fallback is forced. If the crystal ceiling does **not** clear P2Rank and
  the n=1 null, no generator can, and the cheap pipeline is the honest one — a real,
  reportable conclusion reached without a GPU.

**It cannot say:**
- That the **M2 gate passed**. The `boltzgen_consensus` arm does not exist here. A
  positive E2 shows headroom, not that BoltzGen reaches it.
- Anything about **convergence** (PROJECT_GOAL 1.2 / task C5). With 25–31 ligands per
  target you cannot speak to the behaviour of 10,000 designs. In fact C2 shows the
  consensus does not beat n=1 on this data, which is evidence *against* the aggregation
  premise, not for it.
- Anything about **cross-chemotype transfer**, unless the low-Tanimoto stratum of C11 is
  reported with its own n.
- Anything about **affinity, repurposing or ranking**. This is I1/I6 geometry only;
  I2/I3 (recall@k, temporal holdout) are untouched. Per PROJECT_GOAL G8, no
  affinity/potency language may appear near these numbers.
- That the signature is **small-molecule addressable** (`sm_addressable`, §1.4b) — not
  evaluated.

---

## 3. Not evaluated / why

| Item | Why |
|---|---|
| `source="boltzgen_consensus"` (the M2 gate arm) | Modal authenticates; every GPU tier (T4, L4, A10G, A100) refused pending a payment method. Verified, commit 2a13b07. |
| Convergence vs design count (C5 / §1.2) | Requires the generated ensemble. 25–31 crystal ligands cannot stand in. |
| I2 recall@k, EF, BEDROC; I3 temporal holdout | Out of scope for E2; no co-folding run. |
| 8 structures with no P2Rank pocket | P2Rank returned zero pockets; must be rows, not zeros. |
| 4 ligands unfingerprintable (`25X 34B 39B HGB`) | CCD entry unavailable or RDKit refuses to sanitise; excluded from C11 only. |
| Altloc / occupancy sensitivity | Not measured here; `ligand_contacts` behaviour should be stated before the run. |
| 32 cif files not in the manifest | Not part of the benchmark set (C14). |

---

## 4. Minimum bar before any number is published

Ordered. 1–4 are blocking.

1. **Pre-register the consensus→set rule** and report both rules side by side (C1).
2. **Report the n=1 null ladder**; no arm claims anything it cannot beat at rung 3 (C2).
3. **Relabel as I6.2-with-a-missing-arm**, with the verbatim disclaimer in §C3 (C3).
4. **Cluster the statistics by target**; headline carries n_targets and n_trials (C4).
5. Size-match both arms; report precision/recall separately and a rank-based AUC (C5).
6. Filter lipids/detergents/mercurials; report what was dropped; run ±ADORA2A (C6).
7. Deduplicate trials to (target, comp_id) and hold out whole chemical clusters (C7, C11).
8. Key residues by `(chain, seq, icode)`, or exclude F2/PLAU/F10 and re-check (C8, C9).
9. Sweep the cutoff at 4.0/4.5/5.0; a conclusion must survive all three (C10).
10. Publish the **full 60-row per-target table**, negative rows included (CLAUDE.md).
11. Assert `n_structures_parsed == 1531` and join on the manifest (C14).

---

## 10. How the numbers here were produced

Pilot run by this review, not by the E2 implementation, using the repo's own verified
modules (`scripts/interfaces.py`, `scripts/metrics.py`) via `./env/bin/python`.

- **Targets:** ESR1, BRD4, CDK2, F2, CA2, MAPK14 (6 of 60) — chosen for high ligand count
  and resolution spread, not for outcome.
- **Trials:** 173 leave-one-ligand-instance-out; 161 paired with P2Rank.
- **Contacts:** `ligand_contacts(..., cutoff ∈ {4.0, 4.5, 5.0}, with_buried_area=False)`.
- **Consensus:** frequency map over the held-in contact sets; `core` = freq ≥ 0.6
  (PROJECT_GOAL 4.3 default); `top-|O|` = the |O| highest-frequency residues.
- **P2Rank:** rank-1 pocket parsed from `data/raw/p2rank/chunk*/out/**/*_predictions.csv`
  (1,531 structures), cross-checked against the completed
  `results/p2rank_pockets.json`.
- **Chemistry:** `ccd_mol` → RDKit Morgan r=2/2048; 12 targets, 310 ligands, 3,875 pairs.
- **Statistics:** `metrics.paired_bootstrap`, `metrics.wilcoxon`. Null-1 = 40 resamples
  per trial, `random.seed(0)`.

**Caveat on this critique's own numbers:** the pilot is 6 targets, and C4 is the reason to
distrust any 6-target estimate as a point value. The pilot establishes that the C1 flip,
the C2 null failure and the C4 clustering inflation **exist and are large**; it does not
establish their magnitude on the full 60-target set. That is the implementation's job —
which is exactly why these mitigations belong in the design rather than in a post-hoc
defence of a number already published.
