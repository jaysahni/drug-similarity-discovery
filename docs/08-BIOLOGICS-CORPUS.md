# 08 — The approved-biologics corpus

`data/approved_biologics.csv`, built by `demo/build_biologics.py`.

## Why it exists

`data/approved_drugs.csv` has 4,099 approved entries and **zero biologics** —
verified programmatically, not assumed: adalimumab, pembrolizumab, etanercept,
somatropin, trastuzumab, bevacizumab and insulin glargine are all absent. It is
SMILES-derived, so anything too large to write as SMILES simply is not in it:

```
4,099 approved entries in data/approved_drugs.csv
  4,065  small molecules
     34  peptide-scale (semaglutide, exenatide, teriparatide, enfuvirtide…)
      0  antibodies or large biologics
```

A SMILES string cannot be embedded by a protein language model, so the biologic
half of the approved library had no representation at all. `demo/esmc.py` needs
sequences; this file supplies them.

## Contents, n = 262

| modality | n |
|---|---|
| mab | 109 |
| hormone_cytokine | 47 |
| other | 34 |
| enzyme | 28 |
| peptide | 26 |
| fusion_protein | 9 |
| antibody_fragment | 9 |

222 of 262 carry at least one UniProt target accession; 251 of 262 a
`first_approval` year. Total 129,913 residues.

**Chain format**: `chains` holds `ID:SEQUENCE|ID:SEQUENCE` with `H*`/`L*`/`A*`
identifiers, plus separate `heavy_chain_sequence` / `light_chain_sequence`
columns. `demo/esmc.py:embed_entity` combines multi-chain entities by
equal-weighted mean of per-chain vectors (see `docs/09`).

## Verification

All computed, none estimated.

| Check | Result |
|---|---|
| A — INN N-terminal motifs, independent of ChEMBL, n=11 | 10 pass, 0 fail, 1 not in corpus (adalimumab) |
| B — chain is a substring of its UniProt parent, n=5 | 5 pass/partial, 0 fail |
| C — whole-corpus consistency, n=262 | **0 violations** |

Check C covers: heavy chain longer than light 106/106; kappa/lambda C-terminus
98/106; IgG CH3 C-terminus 101/106; all residues in the standard alphabet
262/262; names unique.

Two partials in check B are expected and remain unproven rather than asserted:
etanercept matches the first 235 aa of TNFRSF1B then diverges (it is an Fc
fusion), and aldesleukin matches 123 of 132 of IL-2 (des-Ala1, C125S).

## Two bugs found while building it

1. **ChEMBL labels 160 of 400 protein components with the literal string
   `Sequence`**, which filed trastuzumab, bevacizumab and ranibizumab heavy
   chains as neutral `C1` instead of `H1`. `assign_chain_ids()` now infers H/L
   only under unambiguous conditions. Verification A went 4/8 → 10/11 as a result.
2. **`approval_agencies` by DrugCentral name-match filled only 16 of 262** —
   because `data/approved_drugs.csv` contains no biologics to match against, which
   is the premise of this file. Replaced with ChEMBL provenance evidence
   (DailyMed/FDA xref → `FDA`, EMA xref → `EMA`) plus an
   `approval_agencies_source` column.

## Licence — read this before redistributing

**The committed CSV is CC BY-SA 3.0, not MIT.** All sequences come from ChEMBL,
which is CC BY-SA 3.0: redistributable with attribution and share-alike. The
repo's MIT licence covers the *code*; this data file does not inherit it.

Other sources and their status:

| Source | Licence | Used? |
|---|---|---|
| ChEMBL | CC BY-SA 3.0 | yes — all committed sequences |
| RCSB PDB | CC0 | stage written, **not executed** |
| UniProt | CC BY 4.0 | gene symbols and verification |
| DrugBank | CC BY-NC 4.0 | **not used, not committed** |
| Thera-SAbDab | no published licence | not committed; opt-in `--therasabdab audit\|merge`, where `merge` refuses to write to the committed path |

## Gaps, stated plainly

- **45 approved antibodies have a ChEMBL biotherapeutic record carrying zero
  sequences**, so they are absent. These include adalimumab, infliximab,
  golimumab, certolizumab pegol, ipilimumab, denosumab, pertuzumab, ustekinumab,
  tocilizumab, eculizumab, dupilumab, natalizumab, palivizumab, panitumumab,
  ocrelizumab, ofatumumab and sotrovimab. **TNF coverage is therefore 1 of 5**
  (etanercept only).
- **251 of 513 ChEMBL candidates are excluded**; the largest bucket (146) has no
  biotherapeutic record at all, including alteplase, asparaginase, botulinum
  toxin A, albiglutide and calcitonin. ChEMBL is simply incomplete for biologics
  and no fully-redistributable source found here closes that.
- An RCSB PDB stage was written to close the antibody gap — a probe showed RCSB
  returns cleanly labelled entities ("Adalimumab Fab Heavy chain", 3WD5/4NYL) —
  but **it has never been run or validated**, so no claim is made for it.
- A Thera-SAbDab audit found heavy chains for **31 of the 45** ChEMBL antibody
  gaps, so that route would close most of it if the licence were resolved.

## Target coverage

| Target | Coverage |
|---|---|
| PDCD1 (PD-1) | 6/6 |
| CD274 (PD-L1) | 3/3 |
| VEGFA | 3/3 |
| KDR (VEGFR2) | 1/1 |
| **F2 (thrombin)** | 3 — bivalirudin (20 aa), lepirudin (65 aa), antithrombin alfa (432 aa) |
| TNF | **1/5** — etanercept only |

F2 is the target `docs/10` uses, and two of its three entries are short enough to
sit in the ≤100 aa matching corpus, giving arm P two real positive controls.
