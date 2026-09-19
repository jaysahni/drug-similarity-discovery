"""Build the benchmark dataset from the DrugCentral release files.

Inputs (downloaded into data/ from unmtid-dbs.net/download/DrugCentral/2021_09_01):
    structures.smiles.tsv          SMILES + INN name per DrugCentral struct id
    drug.target.interaction.tsv    drug -> protein target, with a mechanism flag

Output: data/drugs.csv with columns
    struct_id, name, smiles, targets, moa_targets

`targets` is every human gene the drug is recorded as binding; `moa_targets` is
the subset flagged as the drug's mechanism of action. The benchmark's ground
truth is built from these, so the two columns give a loose and a strict
definition of "similar" to test against.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
STRUCTURES = DATA / "structures.smiles.tsv"
INTERACTIONS = DATA / "drug.target.interaction.tsv"
OUT = DATA / "drugs.csv"


def load_structures():
    with STRUCTURES.open(newline="") as fh:
        return {
            r["ID"]: {"struct_id": r["ID"], "name": r["INN"], "smiles": r["SMILES"]}
            for r in csv.DictReader(fh, delimiter="\t")
            if r.get("SMILES") and r.get("ID")
        }


def load_targets():
    """struct_id -> (all human gene symbols, mechanism-of-action gene symbols)."""
    allt = defaultdict(set)
    moat = defaultdict(set)
    with INTERACTIONS.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid, gene = r.get("STRUCT_ID"), (r.get("GENE") or "").strip()
            if not sid or not gene:
                continue
            if r.get("ORGANISM") != "Homo sapiens":
                continue
            # A single row can list several genes for a protein complex.
            genes = {g.strip() for g in gene.split("|") if g.strip()}
            allt[sid] |= genes
            if (r.get("MOA") or "").strip() == "1":
                moat[sid] |= genes
    return allt, moat


def main():
    structures = load_structures()
    allt, moat = load_targets()

    rows = []
    for sid, rec in structures.items():
        targets = sorted(allt.get(sid, ()))
        if not targets:
            continue  # no ground truth available for this drug
        rec["targets"] = ";".join(targets)
        rec["moa_targets"] = ";".join(sorted(moat.get(sid, ())))
        rows.append(rec)

    rows.sort(key=lambda r: int(r["struct_id"]))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["struct_id", "name", "smiles", "targets", "moa_targets"]
        )
        w.writeheader()
        w.writerows(rows)

    n_moa = sum(1 for r in rows if r["moa_targets"])
    genes = {g for r in rows for g in r["targets"].split(";")}
    print(f"structures with SMILES:        {len(structures)}")
    print(f"  ...with >=1 human target:    {len(rows)}")
    print(f"  ...with >=1 MOA target:      {n_moa}")
    print(f"distinct human genes targeted: {len(genes)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
