"""Rank approved drugs by similarity to the target's best literature ligand, then
check whether that ranking rediscovers the target's KNOWN approved binders.

This closes the pipeline loop (tasks T7 + T10 of docs/02-PIPELINE-PLAN.md) using the
representation the benchmark actually chose, rather than one picked by assumption:
scripts/benchmark.py found nothing beats ECFP4 on this corpus, so ECFP4 is the default
here and the choice is traceable to a measured result.

The rediscovery check is the honest signal. The "ideal binder" for KDR is axitinib, an
approved VEGFR2 inhibitor, and several other approved drugs are annotated against KDR
in DrugCentral. If similarity-to-axitinib is worth anything as a repurposing shortlist,
those drugs should surface near the top; the enrichment over their base rate in the
corpus says by how much. A miss is reported with its rank, not hidden.

NO AFFINITY IS PREDICTED OR REPORTED ANYWHERE. This ranks by chemical similarity only.
PROJECT_GOAL.md G8 bans affinity language around numbers like these and that ban is kept.

Usage:
    ./env/bin/python scripts/match_candidates.py
    ./env/bin/python scripts/match_candidates.py --rep fcfp4 --top 25
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from representations import REPRESENTATIONS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "data" / "approved_drugs.csv"
PIPE = ROOT / "results" / "pipeline" / "colorectal-cancer"


def tanimoto(query_vec, mat):
    inter = mat @ query_vec
    card = mat.sum(1)
    union = card + query_vec.sum() - inter
    union[union == 0] = 1.0
    return inter / union


def cosine(query_vec, mat):
    qn = np.linalg.norm(query_vec) or 1.0
    mn = np.linalg.norm(mat, axis=1)
    mn[mn == 0] = 1.0
    return (mat @ query_vec) / (mn * qn)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep", default="morgan", help="representation (default: the benchmark winner)")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--target", default="KDR")
    args = ap.parse_args()
    if args.rep not in REPRESENTATIONS:
        sys.exit(f"unknown representation {args.rep!r}; known: {list(REPRESENTATIONS)}")

    rows = [r for r in csv.DictReader(CORPUS.open()) if r["smiles"]]
    ideal = [r for r in csv.DictReader((PIPE / "ideal_binder.csv").open())
             if r["selected_as_ideal_binder"] == "True"]
    if not ideal:
        sys.exit("no row flagged selected_as_ideal_binder in ideal_binder.csv")
    ideal = ideal[0]
    print(f"ideal binder: {ideal['pref_name'] or ideal['chembl_id']} "
          f"({ideal['standard_type']} {ideal['standard_value_nm']} nM, "
          f"assay {ideal['assay_chembl_id']})")
    print(f"corpus: {len(rows)} approved-drug structures | representation: {args.rep}\n")

    fn = REPRESENTATIONS[args.rep]
    mat = fn([r["smiles"] for r in rows] + [ideal["smiles"]])
    corpus_mat, query = mat[:-1], mat[-1]
    if not query.any():
        sys.exit("the ideal binder did not encode under this representation")
    sim = (tanimoto(query, corpus_mat) if getattr(fn, "binary", False)
           else cosine(query, corpus_mat))

    # ground truth for the rediscovery check: approved drugs DrugCentral annotates
    # against this target. Verified present rather than assumed.
    known = {r["struct_id"] for r in rows
             if args.target in (r.get("targets") or "").split(";")}
    print(f"{len(known)} approved drugs in the corpus are annotated against {args.target}")
    if not known:
        sys.exit(f"no corpus drug is annotated against {args.target}; nothing to rediscover")

    order = np.argsort(-sim, kind="stable")
    ranks = {}
    for pos, i in enumerate(order, 1):
        if rows[i]["struct_id"] in known:
            ranks[rows[i]["struct_id"]] = (pos, rows[i]["name"], float(sim[i]))

    top = []
    print(f"\ntop {args.top} by similarity to the ideal binder:")
    for pos, i in enumerate(order[:args.top], 1):
        r = rows[i]
        hit = r["struct_id"] in known
        top.append({"rank": pos, "struct_id": r["struct_id"], "name": r["name"],
                    "similarity": round(float(sim[i]), 4),
                    "known_target_binder": hit,
                    "n_annotated_targets": int(r.get("n_targets") or 0)})
        print(f"  {pos:3d}. {r['name'][:38]:<38} {sim[i]:.4f}"
              + ("   <- known " + args.target + " binder" if hit else ""))

    n_hits_top = sum(1 for t in top if t["known_target_binder"])
    base = len(known) / len(rows)
    ef = (n_hits_top / args.top) / base if base else float("nan")
    found = sorted(ranks.values())
    print(f"\nREDISCOVERY ({args.target}):")
    print(f"  {n_hits_top}/{args.top} of the top {args.top} are known binders "
          f"(base rate {base:.4f}) -> enrichment {ef:.1f}x")
    print(f"  ranks of all {len(known)} known binders:")
    for pos, name, s in found:
        print(f"    {pos:5d}  {name[:40]:<40} sim={s:.4f}")
    missed = [k for k in known if k not in ranks]
    if missed:
        print(f"  NOT RANKED (no similarity computed): {missed}")

    out = {
        "representation": args.rep,
        "representation_choice_rationale": (
            "benchmark.py measured 16 representations on target-sharing retrieval; "
            "nothing beat ECFP4 on p@1 (n=2069), so it is the default here"),
        "ideal_binder": {k: ideal[k] for k in
                         ("chembl_id", "pref_name", "smiles", "standard_type",
                          "standard_value_nm", "pchembl", "assay_chembl_id", "pmid")},
        "target": args.target,
        "n_corpus": len(rows),
        "n_known_binders_in_corpus": len(known),
        "base_rate": base,
        "top_k": args.top,
        "n_known_in_top_k": n_hits_top,
        "enrichment_factor": ef,
        "known_binder_ranks": [{"rank": p, "name": n, "similarity": s} for p, n, s in found],
        "unranked_known_binders": missed,
        "top": top,
        "caveat": ("ranked by chemical similarity to the best literature ligand only. "
                   "No affinity, potency or binding strength is predicted or implied. "
                   "The ideal binder is itself an approved drug, so its own appearance "
                   "near rank 1 is a sanity check, not a discovery."),
    }
    dest = PIPE / "candidates.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
