# AutoRepurpose — HackMIT 2026 submission

**Finding approved drugs that fit a disease target they were never designed for.**

Tracks: Healthcare · Regeneron · Voloridge
Repo: <https://github.com/jaysahni/drug-similarity-discovery> (MIT)

---

## What it does

Type a disease name. The system finds the protein driving it, searches every
approved drug for one whose shape fits that protein's binding site, and predicts
the 3D structure of the two locked together.

```
disease  →  target  →  sequence  →  search 2,382 approved drugs  →  co-fold  →  candidate
```

Repurposing matters because most of a drug's decade-long, billion-dollar path to
market is spent proving it is safe in humans. A drug already approved for one
illness has done that work. Finding which approved drug fits a new target is a
search problem — and this automates it end to end.

## The headline result

Given **VEGFR2**, a cancer target, and 2,382 approved drugs to search, with every
label hidden:

| what came back | rank |
|---|---|
| known VEGFR2 cancer drugs *(the control)* | **3, 6, 17, 22** |
| **mebendazole** — a 1974 deworming pill | **24 of 2,382** |

Co-folding mebendazole against the target puts it in the same pocket as the
approved cancer drug, touching **19 of the same 19 amino acids** (overlap 0.950).
Four unrelated approved drugs run as controls reach 0.41–0.64.

Mebendazole-for-cancer is a real line of research, currently in human trials. The
pipeline had no access to that — it was given a protein and a library of
structures and surfaced the right molecule in the top 1% on shape alone.

## Two more examples

**EGFR — the method at its clearest.** Query with one approved EGFR drug and the
others return at ranks **1, 2, 3, 4, 5, 7 and 12** of 2,382; EGFR-linked drugs
are 10.8× enriched in the top 50.

**PADI4 — an undrugged rheumatoid arthritis target.** With no approved drug to
query, the search uses the chemistry the enzyme recognises and returns
**pentamidine** and **hydroxystilbamidine** — both bis-amidines, the exact class
published PAD4 inhibitor programmes are built on. A chemically coherent
hypothesis, and no way to grade it without a lab.

## How it is built

| stage | what runs |
|---|---|
| disease → target | Open Targets, with structure selection from the PDB |
| sequence | UniProt, via the `novakit` toolkit |
| search | ECFP4 fingerprints, best of 16 representations benchmarked here |
| co-folding | Boltz-2 on Rowan |
| contacts | 4.5 Å heavy-atom, validated against published crystallography |

Runs on a laptop plus API credits. No wet lab, no cluster. The whole project
used **under 450 compute credits**.

## What makes it trustworthy

**Every result carries a control.** A ranking is worthless unless you know what a
wrong answer looks like. Known drugs must return first; unrelated drugs must
fail. Both are measured on the same data in the same run.

**Negative results are reported, not hidden.** Several approaches were measured
and found not to work, and they are in the repository with their numbers —
including one where a control we built showed our own apparent success was an
artifact.

**The pipeline refuses work it cannot do.** It rejects targets with no approved
drug to validate against, protein structures with the binding site already
occupied, and targets where inhibition would be the wrong direction — a check
added after a run selected a tumour suppressor, where blocking the protein would
make the disease worse.

## Limits

Everything here is a computational hypothesis and a starting point for
laboratory work. Mebendazole's activity is documented, which is what made it a
valid test of the method rather than a discovery. Confidence scores measure how
sure a model is about a predicted pose, not how tightly a molecule binds.

## Where to look

| | |
|---|---|
| plain-language write-up | `SUCCESSES.md` |
| all findings, with n and significance tests | `docs/11-FINDINGS.md` |
| how to run it | `demo/README.md` |
| methodology | `docs/10-DIRECT-MATCHING.md`, `docs/15-AUTORESEARCH.md` |

Every number in this submission was computed by a script in the repository and
is recorded under `results/`.
