# 01 - Data

Source: **DrugCentral 2021_09_01 release** (unmtid-dbs.net/download/DrugCentral).
Raw TSVs are gitignored (size); `data/drugs.csv` is committed.

| File | Rows | Contents |
|---|---|---|
| `data/structures.smiles.tsv` | ~4.6k | struct_id, INN name, SMILES (raw, gitignored) |
| `data/drug.target.interaction.tsv` | ~30k | struct_id → gene, organism, MOA flag (raw, gitignored) |
| `data/drugs.csv` | 2,114 | struct_id, name, smiles, targets, moa_targets — built by `scripts/prepare_data.py`; **only drugs with ≥1 human target** |

Caveat for the pipeline: `drugs.csv` is filtered to drugs with known human
targets. The embedding-matching corpus (task T1) needs *all* approved
structures with SMILES, not just target-annotated ones.

`scripts/fetch_chembl.py` is an unused fallback — ChEMBL REST was down when
written; retry before assuming it works.

## `data/approved_drugs.csv` — the unfiltered corpus (task T1)

`drugs.csv` above keeps only drugs with >=1 annotated human target, which is
right for the retrieval benchmark and wrong for repurposing: a candidate may
have no target record at all. `scripts/build_drug_corpus.py` builds the
unfiltered corpus.

| | |
|---|---|
| rows | 4,099 (every data row of `structures.smiles.tsv`; none lost) |
| rows with a parseable SMILES | 4,099 / 4,099 = **100%** |
| `approved = 1` | 2,153 |
| `chembl_max_phase4 = 1` | 2,082 |
| flagged `charged_parent` | 140 |
| flagged `inorganic` | 32, plus 12 `inorganic;element` |

Columns: `struct_id, name, smiles, inchikey, approved, approval_agencies,
first_approval, chembl_max_phase4, n_targets, targets, moa_targets, flags`.

### Caveats that downstream code must know about

These were found by an independent verification pass and confirmed here by
recount, not taken on trust.

- **`first_approval` is blank for 290 of the 2,153 approved rows (13.5%).**
  `approved` means presence in DrugCentral's approval table; the date column is
  optional. Anything keying on approval date — a temporal holdout, for instance
  — silently loses an eighth of the corpus unless it handles the blank.
- **The ChEMBL cross-check matches on the 14-character InChIKey skeleton**, so
  stereoisomers and charge-layer variants collapse together. Matching on the
  full 27-character key instead gives 2,002 rather than 2,082 hits: **80 of the
  2,082 depend on the loose match.** The looseness is deliberate (a salt or a
  single enantiomer should still match its parent) but it is not free, and the
  DrugCentral/ChEMBL agreement table is a skeleton-level comparison.
- **Inorganics and elements are flagged, not deleted** (44 rows). A
  fingerprint-based representation will produce a near-empty vector for these;
  they are kept so the count is visible rather than quietly absent.
- `approval_agencies` stores full agency strings, e.g.
  `Korean Food and Drug Administration (KFDA)`, not the short label. A filter
  written against `KFDA` will miss it.
