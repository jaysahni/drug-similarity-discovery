# drug-similarity-discovery

Benchmarking representations for **drug similarity matching** — given a drug, find the
most similar drugs — and then asking whether the same question is better posed on the
*target* side: not "what does this molecule look like" but "which residues does it
engage".

Built at HackMIT 2026. MIT licensed. Every number below was computed by a script in
this repo and can be reproduced by running it.

---

## The three results

**1. Nothing beats ECFP4.** Across 16 representations on a 2,069-query retrieval task,
a 2048-bit Morgan fingerprint is the best thing tested. Three others tie it; the
learned chemical-language embedding is significantly *worse*.

**2. But the metrics disagree, and that matters.** ECFP4 wins the top of the ranking
(p@1) and loses on global ranking quality (AUROC) to four other representations.
Reporting one number would have hidden that.

**3. And where it fails is the interesting part.** Ranking approved drugs by similarity
to a potent VEGFR2 ligand enriches known VEGFR2 binders **10.2×** — and still buries
sunitinib at rank 464 and pazopanib at 764. Chemical similarity finds the drugs that
*look* like the query. The ones that engage the same residues without resembling it are
exactly what it cannot see.

**4. Target-side beats ligand-side — but for a smaller reason than it first appears.**
A consensus of which residues a target's known ligands engage predicts a *held-out*
ligand's contacts far better than a pocket finder does (+0.27 precision, n=1,531). But
a single randomly chosen other ligand already gets most of the way there. Aggregation
buys **+0.11**, not +0.27. The rest is just knowing where the pocket is.

---

## E1 — which representation actually retrieves target-mates?

Given a query drug, rank all others. A hit shares a human protein target.
`n = 2,069` queries over 2,114 DrugCentral drugs (median 87 target-sharing
neighbours; random-pair relevance 0.0650). Floor = random vectors. Ceiling = perfect
ranking. 95% bootstrap CIs; comparisons are paired bootstrap + Wilcoxon with
Holm-Bonferroni over the 32-comparison family.

| representation | p@1 | p@10 | AUROC | BEDROC(α=20) | vs ECFP4 (Holm) |
|---|---|---|---|---|---|
| **morgan (ECFP4)** | **0.6597** [0.6394, 0.6796] | 0.4503 | 0.6692 | 0.3528 | — baseline |
| atom_pair | 0.6588 | 0.4362 | 0.6845 | 0.3561 | −0.0010, p=1.00 **n.s.** |
| fcfp4 | 0.6573 | 0.4505 | **0.6907** | **0.3637** | −0.0024, p=1.00 **n.s.** |
| topological_torsion | 0.6520 | 0.4319 | 0.6695 | 0.3352 | −0.0077, p=0.78 **n.s.** |
| morgan_count | 0.6341 | 0.3999 | 0.6763 | 0.3278 | −0.0256, p=6.4e-03 |
| layered_fp | 0.6303 | 0.4134 | 0.6560 | 0.3190 | −0.0295 |
| avalon | 0.6249 | 0.4129 | 0.6514 | 0.3279 | −0.0348 |
| gobbi_pharm2d | 0.6133 | 0.4099 | 0.6627 | 0.3238 | −0.0464 |
| rdkit_fp | 0.6100 | 0.3901 | 0.6349 | 0.2926 | −0.0498 |
| rdkit_descriptors | 0.6037 | 0.3875 | 0.6865 | 0.3465 | −0.0561 |
| maccs | 0.5988 | 0.3962 | 0.6726 | 0.3382 | −0.0609 |
| pattern_fp | 0.5660 | 0.3713 | 0.6657 | 0.3072 | −0.0938 |
| **chemberta** | 0.5607 | 0.3690 | 0.6849 | 0.3287 | **−0.0991, p=6.4e-03** |
| murcko_scaffold | 0.3586 | 0.1485 | 0.5610 | 0.1394 | −0.3011 |
| usrcat (3D shape) | 0.3422 | 0.2262 | 0.6213 | 0.2236 | −0.3175 |
| mw_logp *(weak control)* | 0.1547 | 0.1046 | 0.6152 | 0.1687 | −0.5051 |
| random *(floor)* | 0.0691 | 0.0665 | 0.4990 | 0.0969 | −0.5906 |
| *ceiling* | *1.0000* | *0.9202* | *1.0000* | *1.0000* | — |

Read this way: ECFP4 gets **66% of the way to ceiling on p@1** from a floor of 7%.

**What's interesting is the disagreement.** On AUROC, `fcfp4` (0.6907),
`rdkit_descriptors` (0.6865), `chemberta` (0.6849) and `atom_pair` (0.6845) all beat
`morgan` (0.6692). ECFP4 is the best at putting *one* correct drug first and is not the
best at ordering the whole list. A benchmark that reported only p@1 would have called
this a clean win.

**ChemBERTa losing is a real result, not a bug.** It was mean-pooled over non-pad
tokens, ran on the full corpus, and its determinism was verified by recompute. A
77M-parameter transformer trained on SMILES is worse than a hash of substructures at
this task — while being *better* on AUROC, which is the same disagreement again.

`mw_logp` is a deliberately weak control (molecular weight + cLogP only). It lands near
the bottom, which is the evidence that the benchmark can distinguish anything at all.

**It replicates on a stricter ground truth.** Re-run against shared *mechanism-of-action*
target (n=1,319 queries, 1,408 drugs, base rate 0.0171): morgan p@1 **0.6960**
[0.6725, 0.7202] is again first, with topological_torsion (−0.011), fcfp4 (−0.011) and
atom_pair (−0.017) again tying it (Holm p = 0.42, 0.42, 0.21). The AUROC inversion
replicates too — fcfp4 reaches 0.8299 against morgan's 0.8139.

```bash
./env/bin/python scripts/benchmark.py            # full panel, ~10 min first run
./env/bin/python scripts/benchmark.py --quick    # skips conformers + neural
./env/bin/python scripts/benchmark.py --truth moa_targets
```

---

## E3 — where chemical similarity fails, and why that motivates E2

The pipeline end to end: disease → target → best literature ligand → similar approved
drugs. For colorectal cancer the target is KDR/VEGFR2 and the ligand resolved from
ChEMBL is **axitinib** (IC50 0.02 nM, pChEMBL 10.70, assay CHEMBL5621867 confidence 9,
PMID 37285684). Rank the 4,099-structure approved corpus by ECFP4 similarity to it —
ECFP4 because E1 measured it to be the best, not because it is conventional — and ask
where the **32 approved drugs DrugCentral annotates against KDR** land.

| | |
|---|---|
| known binders in top 25 | 2 of 25, base rate 0.0078 → **enrichment 10.2×** |
| axitinib | rank 1 — *it is the query; a sanity check, not a discovery* |
| sorafenib / regorafenib / apatinib | 21 / 39 / 51 |
| **sunitinib** | **464** |
| **nintedanib / pazopanib** | **755 / 764** |
| **vandetanib** | **1,760** |
| **gefitinib** | **2,363** |

So similarity-to-axitinib is genuinely enriched — 10× is not nothing — and it is also
useless for most of the drugs that actually bind this target. Every drug in the lower
half of that table is a real VEGFR2 binder that a chemical-similarity shortlist would
never surface.

That is the whole motivation for asking the question on the target side instead. Two
drugs that share a binding site need not share any chemistry, and E2 measures exactly
how much is recoverable from the site rather than the molecule.

*No affinity, potency or binding strength is predicted or implied anywhere in this
ranking; it is chemical similarity only.*

---

## E2 — the same question, asked on the target side

`PROJECT_GOAL.md` argues a binding site is better described by *which target residues a
ligand engages* than by the ligand's own chemistry, because that description is
comparable across molecules sharing no chemistry at all.

**The experiment.** Take targets with many distinct co-crystallised ligands. Hold one
ligand out. Build a consensus "interface signature" from the others. Predict the
held-out ligand's actual contact residues. Compare against P2Rank, a pocket finder that
takes 0.47 s per structure.

**Scale.** 233 targets qualified from DrugCentral's drug-bound human proteins; the top
60 were taken forward — 1,531 structures, 1,566 distinct ligands, 0 fetch failures,
10,735 P2Rank pockets. `n = 1,531` trials, one per held-out ligand.

| arm | precision | recall | F1 | mean size |
|---|---|---|---|---|
| consensus (core ≥ 0.6) | **0.7322** [0.7137, 0.7506] | 0.5930 | 0.6394 | 13.1 |
| p2rank_top1 | 0.4616 [0.4450, 0.4783] | 0.6121 | 0.5047 | 23.4 |
| p2rank_best_of_3 *(oracle, generous)* | 0.6033 | 0.7771 | 0.6520 | 25.5 |
| union_all *(ceiling)* | 0.3584 | 0.9470 | 0.5002 | 57.2 |
| random *(floor)* | 0.0737 | 0.0608 | 0.0645 | 13.1 |

Δ consensus − p2rank_top1 = **+0.2706** [+0.2498, +0.2913], Wilcoxon p=7.9e-113
(n_eff 1372), Holm-adjusted 8e-04. Per-target-weighted (every target counts once):
**+0.2736** [+0.2034, +0.3422], n=60 — so it is not carried by the ligand-rich targets.

### The number that matters more than the headline

A null that costs nothing — predict the held-out ligand's contacts from **one randomly
chosen other ligand** — already scores precision **0.6193**, beating P2Rank by +0.158.
The consensus beats *that* by **+0.1130** [+0.1036, +0.1225].

So of the +0.27 over a pocket finder, roughly **+0.16 is "knowing that ligands bind
here at all"** and only **+0.11 is what aggregating many ligands buys.** Most of what
any arm predicts is that the pocket is the pocket, which is exactly what P2Rank returns
for free.

### What this is *not*

It is **not** the M2 gate. `PROJECT_GOAL.md` §8.3 defines ablation I6.1 as
*BoltzGen-generated* signature vs pocket geometry, and explicitly calls the
`known_ligand` builder used here the **upper bound** — ablation I6.2. Running I6.2 and
reading it as I6.1 would be an overclaim, so the artifact records
`m2_gate_status: "neither passed nor failed"`.

The asymmetry is the useful part: a *negative* result here would have foreclosed the
BoltzGen arm a fortiori — no generator beats P2Rank through a ceiling that doesn't. A
positive result, which is what we got, licenses nothing about BoltzGen and only shows
the headroom it would have to reach.

### Robustness, and what fails

An adversarial critique (`docs/05-E2-DESIGN-CRITIQUE.md`, 15 findings) was run
*against* this result before it was believed. It survived every blocking check:

| check | result |
|---|---|
| alternative set rule (size-matched top-\|O\|) | +0.2618 precision — agrees with core@0.6 |
| contact cutoff 4.0 / 4.5 / 5.0 Å | +0.300 / +0.271 / +0.248 |
| resolution confound | Spearman ρ=+0.138, p=0.29, n=60 — none |
| drop all lipid/detergent held-out trials | +0.2725 |
| drop insertion-code targets (F2, PLAU, F10) | +0.2783 |
| drop trials where either arm predicts nothing | +0.2617, n=1,214 |

Kept negative results:

- **The consensus loses on recall** (−0.019), and the two tests disagree — bootstrap
  p=0.11 (n.s.), Wilcoxon p=4.3e-05 (significant). Both are reported.
- **It ties the generous oracle baseline on F1** at the plan's default threshold
  (−0.0125, n.s.); 25 of 60 targets have consensus F1 below it.
- **Empty-core failure in 151/1,531 trials (9.9%)**, concentrated in exactly 6 targets
  — and those same 6 are the *only* targets where the consensus loses. Diagnosed: their
  PDB ligands occupy more than one site (ADORA2A's set includes cholesterol and oleic
  acid on the membrane-facing surface). The `known_ligand` builder assumes one site per
  target and that assumption is what breaks.
- **The threshold costs recall.** `union_all` recovers 94.7% of held-out contacts, so
  the information is there; the 0.6 core keeps 0.593 of it. The frequency threshold buys
  precision by discarding about a third of the recoverable signal.

---

## BoltzGen — the generated arm

Reachable via Rowan (no local GPU). `protein-anything` against KDR/VEGFR2 (3VHE chain A,
ligand-stripped), A100-80GB. Two findings, both from running it:

**1. BoltzGen does not use author numbering.** It indexes residues `1..N` over residues
*present in the file*. Our 303-residue construct is authored 811–1169, so every P2Rank
pocket residue (840…1047) was out of range and the first run died with
`BoltzGen exited with code 1`. Cost: ~1 credit.

**2. `binding_types` alone does not aim the design; `include_proximity` does.**

| | `binding_types` only | + 12 Å `include_proximity` |
|---|---|---|
| designs hitting the requested pocket | 0, 2, 0, 0 of 17 | **6, 12, 12, 15 of 17** |
| binder centroid distance to pocket | 28–52 Å | 26–36 Å |
| mutual contact-set Jaccard | 0.087 | **0.441** |
| consensus core residues | **0** | **29** |
| core residues on the requested pocket | 0 | **11 of 17** |

Naming the hotspot residues as a binding constraint left the designs scattered across
the protein surface with an *empty* consensus. Restricting the presented target surface
is what worked. Measured cost: **~5–7 credits per design**, 396–564 s per 4-design batch.

---

## Not evaluated, and why

| Item | Why |
|---|---|
| Ablation I6.1 proper (BoltzGen consensus vs pocket) at scale | ~5–7 credits/design and a published ~3-in-100 filter pass rate put a *filtered* consensus beyond the free tier. Targeting is now solved, so this is a budget limit, not a method limit |
| Boltz-2 co-folding of the approved library | needs GPU-hours; `boltz` also cannot install on Python 3.14 |
| Predicted affinity (pIC50, Kd, IC50) | no affinity model is run anywhere, so no affinity number is reported anywhere |
| Approved peptide + biologic tiers | only scoreable with a co-folding backend |
| Docking as a co-folding substitute | checked, not assumed: `vina` 1.2.7 publishes no cp314 wheel and its sdist fails to build on Python 3.14 |
| Multi-conformer USRCAT; pharmacophore >150 heavy atoms | cost; recorded in `results/representations_selftest.json` |
| Modal GPU | account authenticates and runs CPU, but all four GPU tiers are refused pending a payment method |

---

## Layout

| Path | Contents |
|---|---|
| `scripts/benchmark.py` | E1 — the retrieval benchmark |
| `scripts/representations.py` | the 17-representation zoo |
| `scripts/metrics.py` | ranking metrics, bootstrap CIs, paired tests, Holm |
| `scripts/interfaces.py` | target-side contact extraction + InterfaceSignature |
| `scripts/hotspot_recovery.py` | E2 — leave-one-ligand-out hotspot recovery |
| `scripts/run_p2rank.py`, `build_interface_set.py`, `fetch_structures.py` | the E2 baseline arm and its data |
| `scripts/match_candidates.py` | E3 — similarity shortlist + the rediscovery check |
| `scripts/boltzgen_signature.py` | BoltzGen run → InterfaceSignature + targeting diagnostics |
| `docs/03-SCOPE-AND-CONSTRAINTS.md` | what this machine can and cannot do, stated before building on it |
| `docs/05-E2-DESIGN-CRITIQUE.md` | the adversarial critique E2 had to survive |
| `results/` | computed outputs only — never hand-edited |

## Reproduce

```bash
python3.14 -m venv env && ./env/bin/pip install -r requirements.txt
./env/bin/python scripts/benchmark.py                 # E1
./env/bin/python scripts/build_interface_set.py       # E2 data  (~50 min, network)
./env/bin/python scripts/fetch_structures.py          # ~1 GB of mmCIF
./env/bin/python scripts/run_p2rank.py                # baseline arm (~12 min)
./env/bin/python scripts/hotspot_recovery.py          # E2
./env/bin/python scripts/match_candidates.py          # E3
./env/bin/python scripts/test_metrics.py              # 23 tests
```

P2Rank needs a JVM; `docs/03-SCOPE-AND-CONSTRAINTS.md` records how one was obtained
without Homebrew. BoltzGen needs `ROWAN_API_KEY` in a gitignored `.env`.
