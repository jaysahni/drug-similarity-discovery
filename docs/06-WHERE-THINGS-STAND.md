# 06 — Where things stand

Written as a handoff. If you are picking this up cold, read this file and the
README, in that order, and you will know what exists and what it is worth.

## The one-paragraph version

The repo does two separate things. It **benchmarks drug-similarity
representations** (the original brief in `CLAUDE.md`), and it **implements
`PROJECT_GOAL.md`'s v0.2 repurposing pipeline** — binding site → designed
binders → interface signature → co-folded approved drugs → ranked matches. Both
produce real, computed numbers. The pipeline demonstrably ranks known binders
above chemically unrelated drugs and finds four that chemical similarity buries.
It does **not** separate them from other kinase inhibitors, and the BoltzGen
design step is measurably worse than a free pocket finder.

## What you can run right now, offline, with no credits

```bash
./env/bin/python scripts/ask.py --list            # what can be asked about
./env/bin/python scripts/ask.py --target KDR      # ranked possible matches
./env/bin/python scripts/report.py                # regenerate the HTML report
./env/bin/python scripts/benchmark.py --quick     # the representation benchmark
./env/bin/python scripts/test_metrics.py          # 23 tests
./env/bin/python scripts/lint_language.py         # the G8 language ban
```

`scripts/autorepurpose.py run --target X --uniprot Y` runs the live pipeline.
It needs network and Rowan credits and takes ~5.4 credits per drug.

## What was measured, with the number that matters

| Question | Answer |
|---|---|
| Best representation for target-sharing retrieval? | **ECFP4.** Nothing beat it, n=2,069. ChemBERTa is significantly worse on p@1 and better on AUROC — the metric picks the winner |
| Does the pipeline rank known binders up? | **Yes** against chemically unrelated drugs: AUC 0.913, p=0.0003 |
| Does it beat chemical similarity? | **Yes, and this is the point.** pazopanib/sunitinib/nintedanib/vandetanib rank 4/5/9/10 by interface overlap and 764/464/755/1760 by chemistry |
| Does it separate *hard* decoys (other kinase inhibitors)? | **No.** Best AUC 0.717, p=0.092. The design arm is 0.529 — a coin flip |
| Does the BoltzGen design step earn its cost? | **No.** It loses ablation I6.1 on two targets (KDR −0.324, CDK2 −0.155 Jaccard, Holm p=0.0016) and loses inside the product (2.05× vs 2.87×) |

## The three things most likely to mislead you

1. **The 2.87× enrichment is carried by the easy decoys.** Drop them and it is
   1.83×. `results/role_separation.json` has the significance tests; quote them
   with the enrichment, never the enrichment alone.
2. **Hard decoys near the top may be labelling gaps, not errors.** lapatinib (2),
   pemigatinib (8) and bosutinib (11) are promiscuous kinase inhibitors and
   DrugCentral's annotation is incomplete. This run cannot tell a false positive
   from a missing annotation, so AUC 0.529 is a *lower bound* on discrimination.
3. **Known binders in the board are positive controls, not discoveries.** They
   were put there to see whether the ranking works. A demo that presents
   sunitinib as a novel finding is misreading its own output.

## Not evaluated, and why

| Item | Why |
|---|---|
| A **filtered** BoltzGen ensemble (ipTM > 0.85) | 0 of 24 designs reached even 0.5. ~333 designs ≈ 1,188 credits against a 500-credit tier |
| CDK2 drug board | the target is prepared (pocket, structures, designs) but no co-folding was run |
| Any second disease | only KDR/VEGFR2 has a full screen |
| Predicted binding strength | <!--LANG-EXEMPT-->no affinity model runs anywhere; PROJECT_GOAL.md 4.4 removes it from the ranking path and G8 bans the vocabulary<!--/LANG-EXEMPT--> |
| Visual rendering of the HTML reports | no browser in the dev sandbox; verified structurally only. **Open `results/report_colorectal-cancer.html` before showing it to anyone** |

## Obvious next moves, in order of value per hour

1. **Screen a second target end to end** (CDK2 is already prepared — it needs
   only `repurpose.py shortlist/submit/collect/score`). One target is the
   weakest thing about every claim here.
2. **Fix the hard-null result or state it louder.** Either find a scoring
   variant that separates kinase inhibitors, or lead with the limitation.
3. **Resolve the labelling ambiguity** — check ChEMBL for measured KDR activity
   of lapatinib/bosutinib/pemigatinib. If they do bind, they are not decoys and
   the AUC is being unfairly penalised.
4. **Rotate the Rowan key.** It was pasted into a session transcript.

## Environment gotchas that cost time

All in `docs/03-SCOPE-AND-CONSTRAINTS.md` and the session memory, but the two
that bite hardest: **BoltzGen indexes residues 1..N over residues present in the
file, not author numbering** (a mismatch fails the run outright), and
**`binding_types` alone does not aim a design — `include_proximity` does.**
