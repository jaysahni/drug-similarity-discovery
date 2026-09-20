"""The approved-drug shortlist the demo co-folds against the target.

Co-folding costs real credits (~13 per drug against a ~300-residue kinase), so
the full 4,099-entry approved library in data/approved_drugs.csv is out of reach
for a demo. This module makes that cut explicit, and makes it *interpretable*:
the shortlist is three tiers, not a top-k from some upstream score.

  positive   annotated CDK2 binders. If the pipeline works, these rank high.
             If they don't, the run has told us something and we report it.
  candidate  approved kinase inhibitors with no CDK2 annotation. These are the
             actual repurposing question.
  decoy      non-kinase drugs. These calibrate the score: without them a
             coverage number has no scale.

Selection is by name, fixed in source, so a re-run scores the same molecules.
A shortlist chosen by an upstream model would make the leaderboard a measure of
that model instead.
"""

from __future__ import annotations

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LIBRARY = REPO / "data" / "approved_drugs.csv"

TIERS: dict[str, list[str]] = {
    "positive": [
        "palbociclib",
        "ribociclib",
        "trilaciclib",
        "ceritinib",
        "lapatinib",
    ],
    "candidate": [
        "dasatinib",
        "erlotinib",
        "baricitinib",
        "tofacitinib",
    ],
    "decoy": [
        "paracetamol",
        "atorvastatin",
        "warfarin",
        "furosemide",
        "loratadine",
    ],
}


def load(target_symbol: str = "CDK2") -> list[dict]:
    """Return the shortlist, annotated with tier and known-target status.

    Raises if a named drug is absent from the library rather than quietly
    shipping a shorter list than the source says.
    """
    rows = {r["name"].lower(): r for r in csv.DictReader(LIBRARY.open())}

    shortlist = []
    missing = []
    for tier, names in TIERS.items():
        for name in names:
            row = rows.get(name)
            if row is None:
                missing.append(name)
                continue
            targets = [t for t in row["targets"].split(";") if t]
            moa = [t for t in row["moa_targets"].split(";") if t]
            shortlist.append(
                {
                    "name": name,
                    "tier": tier,
                    "smiles": row["smiles"],
                    "struct_id": row["struct_id"],
                    "n_annotated_targets": len(targets),
                    "hits_target": target_symbol in targets,
                    "target_is_moa": target_symbol in moa,
                    "first_approval": row.get("first_approval") or None,
                }
            )

    if missing:
        raise SystemExit(
            f"shortlist names absent from {LIBRARY.name}: {missing}. "
            "Fix the name or drop it from TIERS -- do not run a partial library silently."
        )
    return shortlist


def summarise(shortlist: list[dict]) -> dict:
    counts: dict[str, int] = {}
    for drug in shortlist:
        counts[drug["tier"]] = counts.get(drug["tier"], 0) + 1
    return {
        "n": len(shortlist),
        "by_tier": counts,
        "library_file": str(LIBRARY.relative_to(REPO)),
        "library_size": sum(1 for _ in LIBRARY.open()) - 1,
    }
