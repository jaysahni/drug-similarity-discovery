# 07 — BoltzGen behaviour: is ipTM 0.16 expected, our fault, or fixable?

Investigation prompted by the lightweight demo's first CDK2 run (workflow
`e715cb2e-c3d8-45c5-a9a1-4bcf248b4cf7`, 12 designs, 51.95 credits) coming back at
ipTM 0.134–0.197 with a 3-residue consensus core that recovered **none** of the
known ligand contacts.

Every claim below is labelled **MEASURED** (computed here, on this repo's data,
number and n given), **PUBLISHED** (cited URL), or **INFERRED** (argued, not
measured). No Rowan workflow was submitted for this document; all Rowan data is
read from artifacts already in the repo.

---

## 0. Verdict in one paragraph

**Three separate things were wrong, and only one of them is ipTM.** The holo
target cost us *targeting* — 4/12 designs on the requested pocket versus 24/24 in
the heavy pipeline's stripped run (Fisher p = 1.6 × 10⁻⁵) — while leaving ipTM
statistically unchanged (Mann-Whitney p = 0.14). The plan's ipTM > 0.85 filter
is **community folklore imported from AlphaFold2/BindCraft tooling** and has no
BoltzGen provenance; the paper publishes no ipTM distribution at all. And
`budget == num_designs` in both our run and the heavy run **switched BoltzGen's
own filtering funnel off**, so every design entered the consensus. Design *count*
is measurably not the problem. But none of these fixes changes the conclusion
already in `results/m2_gate_CDK2.json`: a correctly stripped, 24-design consensus
still loses to pocket geometry on 31 held-out structures.

---

## 1. What ipTM should we expect?

### 1.1 The 0.85 threshold has no BoltzGen provenance

`PROJECT_GOAL.md:492` says, of the filter thresholds:

> Sensible defaults from the community: ipTM > 0.85, pTM > 0.8, complex RMSD < 2.5 Å

**PUBLISHED.** That exact triple traces to a review of scoring across binder
design tools "including BindCraft and others" ([search result summary, Predicting
Experimental Success in De Novo Binder
Design](https://www.biorxiv.org/content/10.1101/2025.08.14.670059.full.pdf)).
BindCraft-family tools compute ipTM with **AlphaFold2-Multimer on a hallucinated
complex**. BoltzGen's `iptm` comes from **Boltz-2 refolding a diffusion-generated
design**. These are different models scoring different objects; there is no
published mapping between their scales.

**MEASURED.** The number is also not reachable in either of this repo's two CDK2
runs. Heavy run max = 0.257 (n = 24). Demo run max = 0.197 (n = 12). Neither run
produced a single design within 0.59 of the threshold.

### 1.2 The BoltzGen paper publishes no ipTM distribution

**PUBLISHED.** Fetched [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729/):
the paper reports "confidence metrics (pTMs, pAEs)" as inputs to filtering but
gives **no ipTM values or thresholds**. Its only stated numeric filter is
`RMSD < 2.5 Å` against Boltz-2 refolded models, for the small-molecule campaign.

So **there is no published BoltzGen baseline to compare 0.16 against.** Anyone
asserting 0.16 is "bad for BoltzGen" is extrapolating, this document included.

### 1.3 What the wider literature uses

**PUBLISHED.** De novo miniprotein work commonly filters at **ipTM > 0.5**, not
0.85 ([ProtDBench / sweet-receptor binder
design](https://arxiv.org/pdf/2601.14574)). In that TAS1R2 campaign the *most
native-like* designs reached only **ipTM ≈ 0.55–0.60**. A meta-analysis of 3,766
experimentally characterised binders puts the F1-optimal **iPSAE** cutoff at 0.61
with "elite" at ≥ 0.8 — a different metric, quoted here only to show that 0.85 is
at the top of the range even for the metric it was coined for.

### 1.4 Reading

**INFERRED.** `ipTM > 0.85` is the wrong gate for this stack and should be
restated in `PROJECT_GOAL.md` as what it is — an AlphaFold-era heuristic — rather
than as a BoltzGen specification. **But our designs are low-confidence by any
published threshold**, including the permissive 0.5. Both things are true.

The more useful finding is §2.3: in our data ipTM did not predict whether a
design landed on the requested pocket at all.

---

## 2. Apo vs holo: it matters, and here is the number

### 2.1 The heavy pipeline stripped. We did not.

**MEASURED.** `results/pipeline/cdk2-second-target/target/pocket.json` records:

```json
"structure_path": "…/structure/6Q4G_A_stripped.pdb",
"apo_or_holo": "stripped_holo",
"detector": "P2Rank 2.5 (default model), run on the ligand-stripped structure"
```

Both `6Q4G_A_holo.pdb` and `6Q4G_A_stripped.pdb` exist in that directory; the
pipeline used the stripped one. The demo's run 1 called
`rowan.create_protein_from_pdb_id("6Q4G")`, which **retains non-polymer
entities** — verified by counting the downloaded structure: **HJK 40 atoms (one
ligand at two altlocs) and 214 waters**.

### 2.2 The effect is on targeting, and it is large

**MEASURED**, CDK2, same target, same `include_proximity` radius 12:

| | heavy (stripped), n = 24 | demo run 1 (holo), n = 12 |
|---|---|---|
| designs overlapping the requested pocket | **24 / 24** | **4 / 12** |
| mean pairwise Jaccard between designs | 0.397 | 0.169 |
| consensus core residues | 19 | 3 |
| Jaccard, core vs known-ligand contacts | 0.338 | **0.000** |
| median contacts per design | 28.5 | 17.0 |

Fisher exact on the on-target counts: **p = 1.64 × 10⁻⁵**.

### 2.3 But ipTM is blind to it

**MEASURED.** ipTM distributions of the two runs are not distinguishable:

| run | n | min | median | mean | max |
|---|---|---|---|---|---|
| heavy, stripped | 24 | 0.114 | **0.1685** | 0.1819 | 0.257 |
| demo, holo | 12 | 0.134 | **0.1615** | 0.1620 | 0.197 |

Mann-Whitney U = 188.5, **p = 0.1397** — no significant difference.

So the holo ligand degraded *where the designs bound* without degrading *how
confident the model was about them*. **ipTM is not a targeting filter.** Any
gating strategy built on ipTM alone — including PROJECT_GOAL.md §4.3's — would
have passed the holo run and the stripped run equally.

The heavy run's own data says the same thing internally: **all 24 of its designs
overlapped the requested pocket** (median 16 of 21 pocket residues engaged), so
within that run there is no on-target/off-target contrast for ipTM to explain.

### 2.4 No published source requires stripping

**PUBLISHED.** The paper says only "Unless mentioned otherwise, we provide the
structure of the targets as input to BoltzGen" — no apo/holo statement. [Rowan's
BoltzGen guide](https://www.rowansci.com/blog/how-to-run-boltzgen) notes its
protein editor "makes it easy to sanitize PDB files and remove existing chains,
binders, or small-molecule ligands" but does not mandate it. A targeted search of
the BoltzGen repo docs surfaced nothing on HETATM handling.

**So §2.2 is the evidence.** This repo's measurement is, as far as this search
found, the only quantitative statement on the question.

---

## 3. Parameters

### 3.1 `budget` is a filter funnel, and we switched it off

**MEASURED**, from `env-kit/lib/python3.12/site-packages/stjames/workflows/protein_binder_design.py:284-286`:

```python
protocol: BoltzGenProtocol = BoltzGenProtocol.PROTEIN_ANYTHING
num_designs: int = 100
budget: int = 20
```

The class docstring documents `num_designs`, `protocol` and `binding_residue` but
**not `budget`**. The defaults imply a **5:1 generate-to-keep ratio**.

**PUBLISHED.** [Rowan's guide](https://www.rowansci.com/blog/how-to-run-boltzgen)
describes a test run using "10 designs" and "a final budget of two designs" —
also 5:1. So `budget` is the count retained after BoltzGen's internal filtering.

**MEASURED.** Demo run 1 used `num_designs=12, budget=12`. The heavy CDK2 run
used `num_designs=24, budget=24` (`results/m2_gate_CDK2.json`, and the recovered
`binder_design_settings` shows `{"budget": 24, "num_designs": 24}`). **Both ran a
1:1 ratio — no filtering at all.** Every design, including the worst, entered the
consensus. `num_filters_passed` is **1 for all 12** demo designs, a constant that
discriminates nothing at this setting.

This is the concrete mechanism behind README.md:295's "the filtered one is
genuinely untested".

### 3.2 `quality_score` is a rank, not a quality

**MEASURED**, n = 12. Reported values are exactly `(n-1-i)/(n-1)` truncated to
3 dp — `[1.0, 0.909, 0.818, …, 0.091, 0.0]`, max deviation 4.5 × 10⁻⁴ — and
designs are returned pre-sorted by it.

It carries **no absolute information**, and in this run it is *anti*-correlated
with interface confidence: Spearman ρ(quality_score, iptm) = **−0.445**, p = 0.147
(n = 12, not significant). Concretely:

- design with the highest ipTM (0.197) → quality_score **0.273**
- design with the highest quality_score (1.0) → ipTM **0.141**

**Consequence for `demo/embed.py`:** it selects the "top design" by
`quality_score`, which at `budget == num_designs` means "first in an arbitrary
within-batch ordering that does not track interface confidence." The parent
should either pick by `design_to_target_iptm`, or use the consensus vector, which
`embed.py` already computes and which is the defensible object anyway.

### 3.3 `include_proximity` is what aims the design

**MEASURED** (previous session, KDR, README.md:264-273) — already correct in the
demo, recorded here so it is not lost:

| | `binding_types` only | + 12 Å `include_proximity` |
|---|---|---|
| designs hitting requested pocket | 0, 2, 0, 0 of 17 | 6, 12, 12, 15 of 17 |
| mutual contact-set Jaccard | 0.087 | 0.441 |
| consensus core residues | **0** | 29 |

### 3.4 Protocol and binder length

**PUBLISHED.** [Rowan's guide](https://www.rowansci.com/blog/how-to-run-boltzgen)
shows `140..180` residues for protein protocols and `8..16` for peptides; the
four protocols are `protein-anything`, `peptide-anything`,
`protein-small_molecule`, `nanobody-anything`.

**INFERRED, not measured.** We ask `protein-anything` for a **60–90 residue**
binder against a **15-residue enclosed ATP pocket**. A 60–90mer buries far more
surface than that pocket presents, so most of its interface must land on flat
surface around the cleft — which is consistent with the heavy run's median binder
centroid sitting **26.4 Å** from the pocket (range 14.6–39.2 Å) even when all 24
designs technically overlapped it. `peptide-anything` at `8..16` is the better
geometric match for an ATP site.

`protein-small_molecule` is **not** the answer despite the paper's rucaparib
campaign using it: it designs a protein that binds *the small molecule*, the
reverse of a target-side signature. `docs/04-BOLTZGEN-MODAL.md` records that
`run_design.py` refuses it without `--force`, for this reason.

**No measurement supports switching protocol.** It is a hypothesis costing one
run to test.

---

## 4. Is n = 12 the problem? No — measured

**MEASURED.** `results/m2_gate_CDK2.json`'s `convergence_c5` block, computed by
sub-sampling the heavy run's own 24 designs:

| n designs | mean Jaccard to full consensus | sd | mean core size |
|---|---|---|---|
| 2 | 0.503 | 0.121 | 17.9 |
| 4 | 0.616 | 0.107 | 19.3 |
| 8 | 0.716 | 0.095 | 21.6 |
| **12** | **0.781** | 0.098 | 17.5 |
| 16 | 0.853 | 0.055 | 18.9 |
| 20 | 0.889 | 0.061 | 20.4 |
| 24 | 1.000 | 0.000 | 19.0 |

At n = 12 the contact map is already **78% of the way** to the n = 24 consensus,
and core size is flat (17–22) across the whole ladder. README.md:243 states the
same conclusion: "The contact map had largely converged before the budget ran
out."

The paper's 10,000–60,000 designs per target exist to find **wet-lab-ready
candidates** (10,000 → 100 → 6 for rucaparib; 60,000 → ≤15 nanobodies per target
— [PMC12697729](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729/)). That is a
1-in-1,000 extremes problem. A **converged contact map** is a central-tendency
problem and converges three orders of magnitude sooner. PROJECT_GOAL.md §1.2
anticipated exactly this and C5 measured it.

**Caveat, and it matters:** that ladder was measured on the *stripped* run where
24/24 designs were on target. Our holo designs disagreed with each other (pairwise
Jaccard 0.169), so their consensus was converging on noise. Convergence of a
mis-aimed ensemble is not evidence about a correctly-aimed one.

---

## 5. Bottom line

### 5.1 Ranked by measured effect size

1. **Strip the target.** Effect measured at Fisher p = 1.6 × 10⁻⁵ on on-target
   rate. Costs nothing. *The parent has already fixed this* via
   `Protein.prepare(remove_heterogens=True, keep_waters=False)`.
2. **Restore the filter funnel** — `num_designs=100, budget=20`. This is the arm
   README.md:295 lists as "genuinely untested", and the only one of these changes
   that could plausibly move the M2 gate result.
3. **Stop selecting designs by `quality_score`** (§3.2). Free.
4. Protocol/length change to `peptide-anything`, `8..16`. Hypothesis only.

### 5.2 Cost of the funnel run

**MEASURED**, two points: heavy 24 designs = 55.14 credits (2.30/design); demo 12
designs = 51.95 credits (4.33/design). Per-design cost is *not* constant, and a
two-point linear fit gives ≈ 48.8 fixed + 0.27/design — implying 100 designs might
cost ≈ 75 credits rather than 230–430.

**That fit is n = 2, across runs that differed in target preparation.** Treat
75–430 credits as the range and budget for the top of it. Either way it fits
inside ~2,400 credits.

### 5.3 The honest answer

Fixing all four items above gets us back to the heavy pipeline's stripped
baseline. **That baseline already lost.** `results/m2_gate_CDK2.json`, n = 31
held-out CDK2 co-crystals:

| arm | Jaccard vs known ligand contacts |
|---|---|
| known ligand (ceiling) | 0.679 |
| **p2rank geometry** | **0.494** |
| boltzgen consensus | 0.338 |
| random | 0.033 |

Δ boltzgen − p2rank = **−0.155**, Holm-corrected **p = 0.0016**. The same gate ran
on KDR (n = 38) with Δ = **−0.324**, Holm p = 0.0016.

So: **the pocket-geometry arm is the one to ship.** The highest-value use of the
remaining credits is not rescuing BoltzGen — it is the filtered-consensus run at
`num_designs=100, budget=20`, because that single experiment closes the one
caveat this project has repeatedly flagged as open and cannot currently answer.
If it also loses, the negative result is finally unqualified.

---

## 6. Not investigated / why

| Item | Why |
|---|---|
| Whether BoltzGen's own `parse_pdb` includes HETATM in the design context | Rowan's hosted path is opaque; the local `scripts/modal_boltzgen.py:185` skips HETATM when building the residue map, but that code has never run. §2.2 measures the *effect* without establishing the mechanism |
| Whether `peptide-anything` improves ATP-site targeting | costs one run; no measurement exists either way |
| ipTM calibration between Boltz-2 and AlphaFold2-Multimer | no published mapping found |
| Re-running the heavy KDR comparison under a filter funnel | out of scope; KDR's Δ is twice CDK2's, so CDK2 is the friendlier test |
| `alpha` (Rowan default 0.01 peptides / 0.001 otherwise) | surfaced in Rowan's guide, not exposed by `submit_protein_binder_design_workflow`'s signature |

## 7. How the numbers here were produced

All computed with `env-kit/bin/python` (3.12.13) against committed artifacts;
`scipy` was added to `env-kit/` for the two hypothesis tests.

| Number | Source |
|---|---|
| heavy ipTM distribution, n = 24 | `results/boltzgen_signature_cdk2.json` → `per_design[].iptm` |
| demo ipTM distribution, n = 12 | `results/demo/cdk2/holo-run/03_designs.json` → `designs[].iptm` |
| Mann-Whitney p = 0.1397 | `scipy.stats.mannwhitneyu`, two-sided, the two arrays above |
| on-target 24/24 vs 4/12, Fisher p = 1.64e-05 | heavy: `per_design[].overlap_with_requested_pocket`; demo: intersect `holo-run/04_signature.json` `per_design[].residues_author` with `holo-run/02_site.json` `chosen.residue_ids_author` |
| `quality_score` is a rank | `holo-run/03_designs.json`, compared against `(n-1-i)/(n-1)` |
| Spearman ρ = −0.445 | `scipy.stats.spearmanr(quality_score, iptm)`, n = 12 |
| convergence ladder | `results/m2_gate_CDK2.json` → `convergence_c5` |
| m2 gate arms | `results/m2_gate_CDK2.json` → `arms`, `comparisons` |
| stripped-structure provenance | `results/pipeline/cdk2-second-target/target/pocket.json` → `apo_or_holo`, `structure_path` |
| HJK 40 atoms + 214 waters | HETATM count of the downloaded run-1 target structure |

**Sources:**

- [BoltzGen: Toward Universal Binder Design (PMC12697729)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697729/)
- [BoltzGen preprint (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.11.20.689494v1)
- [How to Design Protein Binders with BoltzGen — Rowan](https://www.rowansci.com/blog/how-to-run-boltzgen)
- [Interpreting BoltzGen Metrics and Filtering — Neurosnap](https://neurosnap.ai/blog/post/interpreting-boltzgen-metrics-and-filtering-in-protein-design/691cccde8b9522d6ffeff2d9)
- [Predicting Experimental Success in De Novo Binder Design: Meta-Analysis of 3,766 Binders](https://www.biorxiv.org/content/10.1101/2025.08.14.670059.full.pdf)
- [De novo design of protein binders targeting the human sweet taste receptor](https://arxiv.org/pdf/2601.14574)
- [Introducing BoltzGen — MIT Jameel Clinic](https://jclinic.mit.edu/boltzgen/)
