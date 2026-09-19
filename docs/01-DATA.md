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
