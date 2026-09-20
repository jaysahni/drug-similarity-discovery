# CDK2 — lightweight AutoRepurpose demo

Target **CDK2** (P24941, 298 aa), design structure **6Q4G**. Computed with the `novakit` toolkit; every number below comes from a file in this directory.

## What ran

| Stage | Toolkit call | Result | Credits |
|---|---|---|---|
| research | `database.fetch_uniprot_entry`, `literature.search_pubmed` | 20 annotated site residues, 10 PubMed hits | 0 |
| site | `rowan.detect_pockets` | 5 pockets, chose rank 3 | 0.05 |
| design | `rowan.design_protein_binder` (BoltzGen) | 12 designs | 51.95 |
| signature | `demo/contacts.py` | 3 core residues from 12 designs | 0 |
| match | `rowan.cofold` | 14 approved drugs co-folded | 51.42 |

## The ATP site is not the top-scoring pocket

| Pocket | Score | Volume Å³ | Residues | Overlap with annotated site |
|---|---|---|---|---|
| 0 | 5.85 | 1336 | 12 | **0** |
| 1 | 4.64 | 932 | 9 | **0** |
| 2 | 3.41 | 1007 | 17 | **1** |
| 3 ←chosen | 3.00 | 918 | 15 | **11** |
| 4 | 1.96 | 607 | 12 | **0** |

Rowan's pocket score ranks the ATP site below pockets with no annotated-site overlap at all. Taking the top-scoring pocket would have designed against the wrong site.

Chosen site, author numbering: `[13, 14, 15, 16, 17, 18, 33, 34, 35, 129, 131, 132, 134, 144, 145]`. Selection rule: the detected pocket with the most residues in common with the UniProt-annotated ATP/Mg site. NOT the top-scoring pocket -- see finding below.

## Designs

12 BoltzGen designs, ipTM 0.134–0.197 (mean 0.162). **0 of 12** reached the ipTM 0.85 gate PROJECT_GOAL.md §4.3 suggests filtering at.

Consensus over 12 designs at frequency ≥ 0.6: **3 core residues**. Mean pairwise Jaccard between designs: **0.1687** — how much the designs agree with each other at all.

Core residues: `[25, 26, 28]`

### Scored against the known ligand (validation only)

| Signature | Jaccard vs known-ligand contacts |
|---|---|
| BoltzGen consensus | 0.0 |
| chosen pocket geometry | 0.4545 |
| UniProt annotated site | 0.48 |

Known-ligand contacts (n=17) are read **after** the signature is built and never used to construct it (task B12).

## Approved drugs, ranked

14 drugs of 4099 in `data/approved_drugs.csv`, by tier: {'positive': 5, 'candidate': 4, 'decoy': 5}. Co-folded with Boltz-2 via `rowan.cofold`; score is the fraction of core hotspot residues the docked pose engages within 4.5 Å.

| # | Drug | Tier | CDK2 annotated | Coverage (pocket) | Coverage (BoltzGen) | ipTM | Engaged residues |
|---|---|---|---|---|---|---|---|
| 1 | atorvastatin | decoy | no | **0.67** | 0.00 | 0.882 | 26 |
| 2 | lapatinib | positive | yes | **0.60** | 0.00 | 0.888 | 24 |
| 3 | trilaciclib | positive | yes | **0.53** | 0.00 | 0.985 | 21 |
| 4 | ceritinib | positive | yes | **0.53** | 0.00 | 0.98 | 22 |
| 5 | baricitinib | candidate | no | **0.53** | 0.00 | 0.964 | 17 |
| 6 | furosemide | decoy | no | **0.53** | 0.00 | 0.849 | 19 |
| 7 | ribociclib | positive | yes | **0.47** | 0.00 | 0.978 | 19 |
| 8 | erlotinib | candidate | no | **0.47** | 0.00 | 0.939 | 21 |
| 9 | palbociclib | positive | yes | **0.40** | 0.00 | 0.979 | 18 |
| 10 | loratadine | decoy | no | **0.40** | 0.00 | 0.833 | 17 |
| 11 | dasatinib | candidate | no | **0.33** | 0.00 | 0.945 | 19 |
| 12 | tofacitinib | candidate | no | **0.33** | 0.00 | 0.951 | 17 |
| 13 | warfarin | decoy | no | **0.33** | 0.00 | 0.889 | 15 |
| 14 | paracetamol | decoy | no | **0.20** | 0.00 | 0.937 | 8 |

### Does the score separate the tiers?

| Tier | n | Mean coverage (pocket) | Mean coverage (BoltzGen) | Mean ipTM |
|---|---|---|---|---|
| positive | 5 | 0.507 | 0.000 | 0.962 |
| candidate | 4 | 0.417 | 0.000 | 0.950 |
| decoy | 5 | 0.427 | 0.000 | 0.878 |

## Caveats

- UNFILTERED consensus: no design reached the ipTM 0.85 gate.
- consensus over 12 designs from ONE structure (6Q4G).
- results/m2_gate_CDK2.json measured this arm against pocket geometry on 31 held-out CDK2 co-crystals and BoltzGen LOST (Jaccard 0.338 vs 0.494, Holm p=0.0016). Stage 5 therefore scores drugs against both signatures, not just this one.
- **Leakage.** Boltz-2 co-folding an approved drug against a target it is already annotated to bind is very likely reproducing a complex in its training set. The positive tier's scores are therefore an upper bound, not a blind prediction.
- **n = 14 drugs, one target, one structure.** No significance test is reported because none is meaningful at this n; the tier means above are descriptive.
- Output is a computationally-ranked, evidence-linked repurposing hypothesis. It requires experimental validation.

