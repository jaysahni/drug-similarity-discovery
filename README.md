# drug-similarity-discovery

Benchmarking representations for **drug similarity matching** — given a drug, find the
most similar drugs — and then asking whether the same question is better posed on the
*target* side: not "what does this molecule look like" but "which residues does it
engage".

Built at HackMIT 2026. MIT licensed. Every number below was computed by a script in
this repo and can be reproduced by running it.

---

## The pipeline

```
binding site  ->  BoltzGen designs binders  ->  consensus INTERFACE SIGNATURE
              ->  co-fold each approved drug WITH the target, unconstrained (Boltz-2)
              ->  score the resulting pose on engagement with that site
              ->  rank by how much of the signature each drug engages
```

The match is computed in **target-side coordinates** — which residues are engaged —
not in chemical space, because a designed miniprotein and a small molecule have no
chemistry in common to compare.
<!--LANG-EXEMPT-->
**No affinity, Kd, IC50 or potency is predicted or read anywhere in this pipeline.**
Boltz-2 can emit an affinity score; it is deliberately never read. Ranking is
interface overlap only.
<!--/LANG-EXEMPT-->

**Ask it something** — instant, offline, no credits:

```bash
./env/bin/python scripts/ask.py --list             # what can be asked about
./env/bin/python scripts/ask.py "colorectal cancer"
./env/bin/python scripts/ask.py --target KDR --json
./env/bin/python scripts/ask.py --explain sunitinib
```

Run the pipeline live for a new target (needs network + Rowan credits):

```bash
./env/bin/python scripts/autorepurpose.py run --target KDR --uniprot P35968
```

Or drive the stages by hand:

```bash
./env/bin/python scripts/repurpose.py shortlist --target KDR
./env/bin/python scripts/repurpose.py submit --max-credits 190
./env/bin/python scripts/repurpose.py collect
./env/bin/python scripts/repurpose.py score --designs v4
```

### It works — and the free pocket finder wins

KDR/VEGFR2. **42 approved drugs co-folded, 41 scored** (cupric oxide produced no drug-like ligand in its pose and is kept in the board as unscored, not dropped) with Boltz-2 and ranked by interface
overlap: 10 known binders, 19 decoys (matched on MW/cLogP where a near neighbour existed; the rest drawn at random, so 6 of the 19 sit far from any positive), and **12 hard decoys** —
approved kinase inhibitors that DrugCentral does not annotate against KDR, each of
which binds *some* ATP pocket by construction. 191 Rowan credits.

| signature | enrichment @25% | median known-binder rank |
|---|---|---|
| **p2rank_geometry** *(free; see timing note)* | **2.87×** | **8.0** |
| boltzgen_consensus *(85 credits)* | 2.05× | 12.5 |
| known_ligand | 0.82× | 16.5 |

*(4.10× is the arithmetic maximum at this base rate.)* 7 of the top 10 are known
binders; known-binder ranks 3, 4, 5, 6, 7, 9, 10, 18, 20, 22 of 41.

| rank | drug | role | score |
|---|---|---|---|
| 1 | *proxibarbal* | *decoy* | 0.923 |
| 2 | *lapatinib* | *hard decoy* | 0.882 |
| 3 | quercetin | known binder | 0.875 |
| 4 | **pazopanib** | known binder | 0.850 |
| 5 | **sunitinib** | known binder | 0.812 |
| 6 | neratinib | known binder | 0.810 |
| 7 | dasatinib | known binder | 0.800 |
| 8 | *pemigatinib* | *hard decoy* | 0.800 |
| 9 | **nintedanib** | known binder | 0.762 |
| 10 | **vandetanib** | known binder | 0.762 |

### Does it beat chemical similarity? Yes — modestly, and only like-for-like

An earlier version of this README claimed these drugs rank "4, 5, 9, 10 by interface
overlap versus 464, 755, 764, 1760 by chemical similarity". **That comparison was
wrong twice over** and it is corrected here rather than quietly dropped.

**It compared different populations.** 4 is a rank out of the 41-drug co-folded board;
764 is a rank out of the whole 4,099-drug approved corpus. The ratio 4099/41 ≈ 100 is
almost exactly the apparent gap — the comparison was measuring the population sizes.

**And those four were the four with the largest raw gap.** Over all 10 known binders,
on a percentile basis against the full corpus, interface overlap is better for **4 of
10**, paired Wilcoxon **p=0.85** — no detectable difference. Sorafenib is 53.7%
structurally against 0.5% chemically.

The honest comparison re-ranks **the same 41 drugs** by chemical similarity
(`scripts/compare_rankings.py`):

| | structural | chemical |
|---|---|---|
| median rank of a known binder (of 41) | **8.0** | 12.5 |
| AUC vs property-matched decoys | **0.913** (p=0.0003) | 0.795 (p=0.011) |
| AUC vs hard decoys | **0.717** (p=0.092) | 0.621 (p=0.356) |

So structural ranking genuinely beats the chemical baseline on the same drugs — by
roughly 1.5×, not 100×. Per drug: pazopanib 4 vs 13, sunitinib 5 vs 9, nintedanib
9 vs 12, vandetanib 10 vs 20.

### The hard decoys break it, and that is the real result

The enrichment above is carried by the *easy* decoys. Mann-Whitney U on the ranking
metric, two-sided, separates the two nulls cleanly:

| signature | known vs easy decoy | known vs **hard** decoy |
|---|---|---|
| p2rank_geometry | U=173.5, **p=0.0003**, AUC 0.913 | U=86.0, **p=0.092**, AUC 0.717 |
| boltzgen_consensus | U=178.0, **p=0.0002**, AUC 0.937 | U=63.5, **p=0.843**, AUC **0.529** |

*n = 10 known binders vs 19 easy / 12 hard.*

**Against the only null that is actually hard, neither arm is significant.** The
pocket signature is marginal (AUC 0.717, p=0.092); the BoltzGen design signature is at
**AUC 0.529 — a coin flip.** Dropping the easy decoys entirely, enrichment falls to
**1.83×** and **1.10×** against a 2.20× maximum.

So the honest verdict: this pipeline separates known binders from *chemically
unrelated* drugs very well, and has **not been shown** to separate them from other
kinase inhibitors. The 2.87× headline is real but it is measuring the easier task.

**The caveat cuts the other way and belongs right here.** lapatinib ranks 2,
pemigatinib 8, bosutinib 11 — all "hard decoys", all promiscuous kinase inhibitors,
and DrugCentral's annotation is incomplete. A well-scoring hard decoy may be
**mislabelled rather than a false positive**, which makes AUC 0.529 a *lower bound*
on true discrimination rather than a measurement of it. Both statements have to sit
together; either alone misleads.

### What had to be measured rather than assumed

**Co-folding is unconstrained**, departing from §F3. F3 wants a pocket constraint so
promiscuous drugs don't dock somewhere irrelevant. Measured, it does the opposite
damage — it forces *everything* into the site:

| | constrained | unconstrained |
|---|---|---|
| axitinib *(binder)* | 16/17 pocket, 13.1 cr | 16/17 pocket, **5.1 cr** |
| aspirin *(control)* | 11/17 pocket | **6/17** pocket, 5.7 cr |

**Size normalisation** (§1.4a/E4). Raw core coverage ranks a 39-residue lipopeptide
decoy **first**, purely for being large (r(score, n_engaged) = +0.641). Five scores
were compared; `precision_in_core` drops that to **+0.135** and the decoy to rank 14.
Enrichment and median rank are *identical* under all five, so the separation is not an
artefact of the choice — only the oversized decoy moves.

### What does not work

**The design step loses, in the product itself.** Against easy decoys the BoltzGen
signature and the free pocket finder track each other closely (Spearman **0.776** on
the ranking metric, 0.855 on core coverage). The harder null separates them — and
separates them *against* BoltzGen: **2.87× vs 2.05×**, median rank 8 vs 12.5, and
AUC 0.529 vs 0.717 on the known-vs-hard-decoy test. An 85-credit design run is beaten
by a pocket prediction costing well under a second per structure in batch. That is
ablation I6.1 appearing a third time.

**A single known ligand is worse than chance** (0.82×). One holo ligand's contacts
describe that ligand, not the site — so aggregating across many binding events does do
real work. It just doesn't need a GPU to do it.

**The designs are not high-confidence.** All 24 sit at ipTM 0.12–0.42 against the
plan's 0.85 filter; selecting the best 6 barely changes the signature. A filtered
ensemble needs ~333 designs ≈ 1,188 credits against a 500-credit tier.

**The decoys are not yet hard.** MW/cLogP matching returns steroids and
perfluorocarbons — matched on bulk properties, but not plausible ATP-site binders.
`--hard-decoys` selects approved kinase inhibitors that miss KDR (34 available); that
is the test this result still needs.

---

## Why the design step is the weak link

**The project's own decisive experiment fails.**
`PROJECT_GOAL.md` is built around one gate: does a BoltzGen design consensus beat a
pocket finder at saying where a ligand binds? It does not. A sub-second-per-structure pocket
prediction reaches Jaccard **0.6355** against held-out ligand contacts; a 24-design,
85-credit, 28-minute-of-A100 BoltzGen consensus reaches **0.3115** (Δ −0.324, Holm
p=0.0016, n=38 paired). The plan says to be willing to accept that outcome; it is
accepted below, with its caveats quantified rather than asserted.

**1. Nothing beats ECFP4.** Across 16 representations on a 2,069-query retrieval task,
a 2048-bit Morgan fingerprint is the best thing tested. Three others tie it; the
learned chemical-language embedding is significantly *worse*.

**2. But the metrics disagree, and that matters.** ECFP4 wins the top of the ranking
(p@1) and loses on global ranking quality (AUROC) to four other representations.
Reporting one number would have hidden that.

**3. And where it fails is the interesting part.** Ranking approved drugs by similarity
to the target's best literature ligand enriches known VEGFR2 binders **10.2×** — and still buries
sunitinib at rank 464 and pazopanib at 764. Chemical similarity finds the drugs that
*look* like the query. The ones that engage the same residues without resembling it are
exactly what it cannot see.

**4. Target-side beats ligand-side — but for a smaller reason than it first appears.**
A consensus of which residues a target's known ligands engage predicts a *held-out*
ligand's contacts far better than a pocket finder does (+0.27 precision, n=1,531). But
a single randomly chosen other ligand already gets most of the way there. Aggregation
buys **+0.11**, not +0.27. The rest is just knowing where the pocket is.

---

## The end-to-end result

Start from a disease, end at a structurally-validated repurposing candidate.
Plain-language write-up: **[`SUCCESSES.md`](SUCCESSES.md)**. Submission summary:
**[`PLUME.md`](PLUME.md)**.

**9. The pipeline recovers a documented repurposing case blind.** Searching 2,382
approved drugs for VEGFR2 by chemical structure alone — no annotation visible to
the ranking — returns the known VEGFR2 cancer drugs at ranks 3, 6, 17 and 22, and
**mebendazole, a 1974 anthelmintic, at rank 24** (top 1.01%, null percentile
0.992). Enrichment of VEGFR2-linked drugs in the top 50: **7.44×**.

Co-folding puts mebendazole in the same pocket as approved VEGFR2 drug axitinib,
sharing **19 of 19 contacts** (Jaccard 0.950). Four unrelated approved drugs
co-folded as controls reach 0.409–0.636. Note that ipTM does *not* discriminate
here — paracetamol scores 0.970 against mebendazole's 0.990 — so the contact
overlap is the evidence, not the confidence score.

**10. EGFR is the cleanest control.** Query with erlotinib and the other approved
EGFR drugs return at ranks **1, 2, 3, 4, 5, 7, 12** of 2,382, enrichment
**10.79×**. It also marks the boundary: chlorpromazine, which touches EGFR on
unrelated chemistry, sits at rank 1,013 — result #3 again, on a third target.

**11. And what an undrugged target looks like.** PADI4 drives rheumatoid arthritis
and has **zero** approved drugs, so there is no query molecule and no positive
control. Seeded with the arginine-mimetic chemotype PAD4 inhibitor programmes
use, the search returns **pentamidine** and **hydroxystilbamidine** — both
bis-amidines, chemically the right class — alongside benzoic acid and phenol,
which resemble a small query only by being small. Nothing in the run can separate
the two. `results/demo/rediscover/PADI4_prospective.json`

---

## Final report

**[`docs/11-FINDINGS.md`](docs/11-FINDINGS.md)** — the five demo findings with their
n and significance tests, the positive controls that make the negatives readable,
the six silent-corruption faults caught along the way, and a retrospective on what
was built wrong and what to do next.

---

## The lightweight demo — the whole pipeline on a toolkit, and what it found

`demo/` runs the PROJECT_GOAL.md pipeline end to end in ~1,500 lines, with the
scientific work done by [`cheminformatics-kit`](https://github.com/jaysahni/cheminformatics-kit)
(`novakit`) rather than by this repo's 12,305-line `scripts/`. Three Rowan-hosted
stages, everything else local on a Mac. See `docs/06-LIGHTWEIGHT-DEMO.md` and
`docs/10-DIRECT-MATCHING.md`.

It was built to answer a specific question — *given a designed binder, which
approved drug resembles it?* — and it answers it in the negative, three times,
each with a control that says the negative is real.

**5. Hotspot coverage does not separate known binders from decoys.** 14 approved
drugs co-folded against CDK2 with Boltz-2. Positives 0.507 vs decoys 0.427,
Mann-Whitney **p = 0.229**, n = 5 per tier. The **top-ranked drug is atorvastatin**
— a statin with no CDK2 relationship. Palbociclib, an approved CDK4/6 inhibitor,
ranks 9th of 14. Boltz-2 *pose confidence* does separate the tiers (ipTM 0.962 vs
0.878, p = 0.0159), but that is not hotspot coverage, PROJECT_GOAL.md §4.4 keeps
confidence out of the ranking path, and it is the quantity most exposed to
training-set leakage. Recorded, not promoted.

**6. Right-sizing the binder fixes the designs; the match still fails.** Switching
BoltzGen from `protein-anything` at 60–90 residues to `peptide-anything` at 8–16
raised interface confidence from ipTM 0.162 to **0.598** — a 60–90-mer cannot
enter an ATP slot, an 8–16-mer fits a protease groove. But none of the resulting
thrombin peptides resembles an approved thrombin peptide: every one's nearest
neighbour is **abarelix**, a GnRH antagonist, at null percentiles 0.29–0.95, with
enrichment factor 0.0.

That negative is believable because **the positive control passes**: the two
approved thrombin peptides retrieve each other, bivalirudin → lepirudin at
**rank 1 of 36**. The metric works; the designs do not pass it. And a hub check
explains abarelix — **77% of randomly shuffled sequences land on it too**.

**7. The null was wrong, and a calibration test caught it.** Every headline number
here is a percentile against a null. Held-out queries must score uniform
percentiles; a KS test said otherwise for Tanimoto — **D = 0.187, p < 0.001**, mean
0.457 instead of 0.5 — because Tanimoto over sparse fingerprints is heavily
discrete and the percentile used a strict `<`, scoring every exact tie against the
candidate. Fixed to mid-rank; both metrics now calibrate (p = 0.416, p = 0.074).
Without it every candidate would have looked less remarkable than it is,
consistently and invisibly.

**8. The small-molecule arm finds nothing either — and replicates result #3.**
591 molecules generated with a ZINC-trained GPT-2, 200 docked into the thrombin
site, all matched against the 2,153-drug approved library by ECFP4. Best
nearest-neighbour Tanimoto **0.579** (dithranol, a psoriasis anthralin), nothing
with an F2 annotation. Docking score *does* track similarity to known thrombin
binders — **Spearman ρ = −0.146, p = 0.040, n = 199** — the project's only
positive result, and one that was invisible at n = 23 (p = 0.703).

Both controls pass, so the negative is about the molecules rather than the method.
Vina ranks the three real thrombin drugs above the generated median, docked in the
same batch. And ECFP4 does group thrombin drugs — argatroban retrieves bivalirudin
at **rank 9 of 2,152** and ximelagatran at 34, against ~1,076 expected by chance —
but its *top five* are all peptidomimetics with no thrombin annotation. Reading that
top-five as the control nearly produced a false negative. Chemical similarity finds
what looks like the query; the co-target drugs it finds, it finds at rank 10–82.

Three bugs of the same shape were found and fixed along the way, each of which
produced ordinary-looking output: a PDB parser that ignored insertion codes and so
**lost 28 of thrombin's 259 residues**, including the entire 60-loop lining the
active site; structure preparation that silently **renames chains** (H→B, L→A);
and a ligand left in the pocket during design, which sent every design to the
wrong site (Jaccard 0.0 against the known ligand contacts). All three are now
guarded by tests.

One finding is about the tools rather than the drugs: **pocket detection cannot
find a serine protease active site.** Six pockets on thrombin and not one contains
Ser195 or Asp102 — a protease active site is a shallow groove across subsites, not
an enclosed cavity, so a cavity finder fragments it.

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
ChEMBL is **axitinib**.
<!--LANG-EXEMPT-->
(Selection input, not a pipeline output: the published assay value is IC50 0.02 nM,
pChEMBL 10.70, assay CHEMBL5621867 confidence 9, PMID 37285684. It is a literature
measurement used to pick the reference ligand — nothing here predicts it.)
<!--/LANG-EXEMPT--> Rank the 4,099-structure approved corpus by ECFP4 similarity to it —
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

<!--LANG-EXEMPT-->
*No affinity, potency or binding strength is predicted or implied anywhere in this
ranking; it is chemical similarity only.*
<!--/LANG-EXEMPT-->

---

## E2 — the same question, asked on the target side

`PROJECT_GOAL.md` argues a binding site is better described by *which target residues a
ligand engages* than by the ligand's own chemistry, because that description is
comparable across molecules sharing no chemistry at all.

**The experiment.** Take targets with many distinct co-crystallised ligands. Hold one
ligand out. Build a consensus "interface signature" from the others. Predict the
held-out ligand's actual contact residues. Compare against P2Rank, a pocket finder that
takes 0.47 s/structure amortised over a 1,531-structure batch at 12 threads, 250 per JVM; a single structure in isolation measured 2.1-2.5 s wall (P2Rank self-reports 1.87 s), so the batch figure is not a single-run cost.

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

## The M2 gate — ablation I6.1, in full

`PROJECT_GOAL.md` §8.3 calls this "the decisive ablation" and §7 makes it the gate the
project turns on: *if a generated consensus gives you what P2Rank gives you for free,
ship the cheap pipeline.*

Ground truth: the residues the known ligand actually contacts (4.5 Å heavy atom) in each
of **38 KDR co-crystals**. All arms scored on the same 38, so comparisons are paired.

| arm | precision | recall | F1 | Jaccard | \|pred\| |
|---|---|---|---|---|---|
| **p2rank_geometry** | 0.8267 | 0.7441 | 0.7718 | **0.6355** | 19.5 |
| known_ligand *(ceiling, I6.2)* | 0.8366 | 0.7321 | 0.7771 | 0.6630 | 19.0 |
| **boltzgen_consensus** | 0.4565 | 0.4890 | 0.4698 | **0.3115** | 23.0 |
| random *(floor)* | 0.0745 | 0.0791 | 0.0763 | 0.0405 | 23.0 |

Δ boltzgen − p2rank, Jaccard **−0.3240** [−0.3526, −0.2941], Holm p=0.0016. It also
loses on precision (−0.370), recall (−0.255) and F1 (−0.302).

**P2Rank lands within 0.03 Jaccard of the known-ligand ceiling.** There is almost no
headroom above pocket geometry on this target for anything to occupy.

### Why this isn't a cheap dismissal

**Design count is not the excuse.** Task C5's convergence study, by sub-sampling the
same 24 designs: Jaccard to the final consensus is 0.850 at N=12, 0.891 at N=16, 0.897
at N=20. The contact map had largely converged before the budget ran out — which is what
C5 predicted, and it means more designs would not obviously rescue it.

**Design quality is the real caveat, and it is measured.** **Zero of 24 designs pass the
ipTM > 0.85 filter the plan suggests. Zero pass even 0.5** (min 0.119, median 0.210, max
0.420). This is an *unfiltered* consensus, and the plan is explicit that the filters do
enormous work. Reaching 10 filtered survivors at the published ~3-in-100 rate needs ~333
designs at 3.56 credits each ≈ **1,188 credits against a 500-credit tier** — 2.4× beyond
what was available.

So the bounded claim: **BoltzGen, as we could afford to run it, is decisively worse than
pocket geometry at locating a ligand's contacts** on one target with 38 held-out ligands.
Whether a properly filtered ensemble closes a 0.32 Jaccard gap is *not answered here*.

```bash
./env/bin/python scripts/m2_gate.py
```

---

## BoltzGen — getting the generated arm to run at all

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
| A **filtered** BoltzGen consensus (ipTM > 0.85) | 0 of 24 designs reached even ipTM 0.5; ~333 designs ≈ 1,188 credits would be needed for 10 survivors, against a 500-credit tier. The unfiltered arm was run and lost; the filtered one is genuinely untested |
| I6.1 across multiple targets | run on KDR only (n=38 held-out ligands, 1 target). `PROJECT_GOAL.md` anticipates several |
| Boltz-2 co-folding of the approved library | needs GPU-hours; `boltz` also cannot install on Python 3.14 |
| Predicted binding strength | <!--LANG-EXEMPT-->no affinity model is run anywhere, so no affinity, Kd or IC50 number is reported anywhere<!--/LANG-EXEMPT--> |
| Approved peptide + biologic tiers | only scoreable with a co-folding backend |
| Docking as a co-folding substitute | checked, not assumed: `vina` 1.2.7 publishes no cp314 wheel and its sdist fails to build on Python 3.14 |
| Multi-conformer USRCAT; pharmacophore >150 heavy atoms | cost; recorded in `results/representations_selftest.json` |
| Modal GPU | account authenticates and runs CPU, but all four GPU tiers are refused pending a payment method |

---

## Does our score beat a free baseline? On one null, no

Every co-folded pose comes back with Boltz-2 confidence fields we did not ask for.
Tested as binder discriminators on the same 41 drugs (`scripts/confidence_gate.py`,
costs nothing — the numbers were already cached):

| signal | vs unrelated drugs | vs other kinase inhibitors |
|---|---|---|
| **ipTM** | **0.995** (p<1e-4) | 0.696 (p=0.129) |
| interface overlap *(ours)* | 0.913 (p=0.0003) | 0.717 (p=0.092) |
| confidence_score | 0.845 (p=0.003) | 0.608 (p=0.41) |
| avg_lddt | 0.768 (p=0.02) | 0.529 (p=0.84) |
| pTM | 0.589 (p=0.45) | 0.608 (p=0.41) |

*n = 10 known binders vs 19 unrelated / 12 hard decoys.*

**A free field beats our structural score on the easy null.** That is reported
because it is true, not because it helps. Two things survive it:

- **On the hard null nothing works** — not ours, not any confidence field, none
  significant. That ceiling is shared with the co-folding backend rather than caused
  by our scorer, which is a more useful diagnosis than "our score is weak".
- **Interface overlap is interpretable and ipTM is not.** Ours names the residues
  engaged and missed, so a hit can be inspected and argued with. ipTM is one number
  with no mechanism attached.

An earlier version of this claim said confidence *does not* discriminate at all. That
was drawn from a single **constrained** probe — aspirin forced into the pocket reaches
ipTM 0.965 — and it does not generalise to the unconstrained runs the board is built
from. The constrained observation is still why constraints stay off; it was just never
evidence about confidence in general.

## The scoring definition

`precision_in_core` — **of the residues this drug engages, the fraction that lie in the
site core.**

    score = |engaged ∩ core| / |engaged|

- **engaged**: target residues with any heavy atom within **4.5 Å** of any ligand heavy
  atom. Hydrogens excluded; waters, ions and crystallisation additives excluded by a
  named list in `scripts/interfaces.py`; contacts restricted to the ligand's own chain.
- **core**: the residues of the chosen signature (P2Rank rank-1 pocket by default).
- **numbering**: author numbering (`auth_seq_id`) throughout. Co-folded poses are
  numbered 1..N over the construct and converted back with the kinase-domain offset
  (834 for KDR), checked against known identities — Cys919, Asp1046, Phe1047.
- **missing data**: a drug whose pose has no drug-like ligand keeps its row with a
  `status` and **no score and no rank**. It is counted out of `n_scored`, never given a
  plausible-looking number.

**Why this normalisation.** Raw core coverage ranks bigger ligands higher for reasons
unrelated to binding: it correlates with the number of residues engaged at r=+0.64 and
put a 39-residue lipopeptide decoy first. Five variants were compared;
`precision_in_core` drops that to r=+0.14. All five are stored per drug so the choice
can be re-checked, and enrichment is identical under all of them.

## Output schema

`results/repurpose_<pipeline>.json`

| field | meaning |
|---|---|
| `ranked_by`, `ranking_metric` | which signature and score ordered the board |
| `cofolding.pocket_constrained` | **false** for the shipped run — see the co-folding note above |
| `versions` | python, rowan-python, stjames, rdkit, biopython, numpy, scipy, P2Rank |
| `signatures{}` | each signature's core residues and provenance |
| `validation_by_signature{}` | enrichment and known-binder ranks per signature |
| `results[]` | per drug: `name`, `role`, `status`, `rank`, `n_engaged`, `engaged_residues`, and `by_signature{}` with all five scores plus `engaged_core` / `missed_core` |
| `n_scored` | excludes failed candidates, which remain in `results[]` |

Companion artifacts: `role_separation.json` (significance tests),
`ranking_comparison_<pipeline>.json` (structural vs chemical, same drugs),
`report_<pipeline>.html` (self-contained report).

## Implemented vs planned

| stage | status |
|---|---|
| disease → ranked targets, ambiguity exposed | **implemented** (`autorepurpose.py targets`) |
| target → structure → pocket (P2Rank) | **implemented** |
| co-folding of approved drugs (Boltz-2 via Rowan) | **implemented**, needs credits |
| engagement scoring + ranking | **implemented** |
| query interface + HTML report | **implemented** (report renders KDR only) |
| BoltzGen design signature | **optional**, and measured to lose to the pocket finder |
| second target end-to-end | **planned** — CDK2 is prepared but never screened |
| report for an arbitrary target | **planned** — `report.py` inputs are hardcoded to KDR |
| filtered design ensemble (ipTM > 0.85) | **planned** — 0 of 24 designs reached 0.5 |

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
| `scripts/repurpose.py` | **THE PIPELINE — site → designs → signature → co-folded drugs → ranking** |
| `scripts/m2_gate.py` | the M2 gate — ablation I6.1, plus the C5 convergence study |
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
./env/bin/python scripts/m2_gate.py                   # the M2 gate (needs the design run)
./env/bin/python scripts/test_metrics.py              # 23 tests
```

P2Rank needs a JVM; `docs/03-SCOPE-AND-CONSTRAINTS.md` records how one was obtained
without Homebrew. BoltzGen needs `ROWAN_API_KEY` in a gitignored `.env`.
