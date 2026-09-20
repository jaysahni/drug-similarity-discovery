"""Compare the structural ranking against a chemical baseline ON THE SAME DRUGS.

This exists because the obvious comparison is wrong in a way that flatters the
structural method. The repurposing board ranks ~41 co-folded drugs; the chemical
shortlist in candidates.json ranks the whole ~4,099-drug approved corpus. Quoting
"pazopanib 4th structurally versus 764th chemically" compares a rank out of 41
against a rank out of 4,099 - different populations and different tasks. In
percentile terms that headline partly inverts: sunitinib is 12.2% structurally and
11.3% chemically, nintedanib 22.0% versus 18.4%.

The honest comparison re-ranks the SAME board by chemical similarity to the same
reference ligand and compares like with like. The structural method still wins,
by a real but smaller margin, and that is what this script computes.

Both rankings are kept strictly separate from the evaluation labels: roles
(known_binder / decoy / hard_decoy) come from the board, the chemical score comes
only from SMILES and the reference ligand, and neither ranking sees the other.

Usage:
    ./env/bin/python scripts/compare_rankings.py
    ./env/bin/python scripts/compare_rankings.py --pipeline colorectal-cancer
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parent.parent
POSITIVE_ROLE = "known_binder"


def load(path: Path):
    if not path.exists():
        sys.exit(f"missing input: {path}\nrun the pipeline first "
                 f"(scripts/autorepurpose.py) or see docs/06-WHERE-THINGS-STAND.md")
    return json.loads(path.read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pipeline", default="colorectal-cancer")
    ap.add_argument("--representation", default="morgan",
                    help="chemical baseline; benchmark.py measured ECFP4 (morgan) best")
    ap.add_argument("--signature", default=None,
                    help="structural signature to compare (default: whatever the board ranked by)")
    args = ap.parse_args()

    pipe = ROOT / "results" / "pipeline" / args.pipeline
    board = load(ROOT / "results" / f"repurpose_{args.pipeline}.json")
    state = load(pipe / "repurpose_state.json")
    ideal_path = pipe / "ideal_binder.csv"
    if not ideal_path.exists():
        sys.exit(f"missing {ideal_path}; the chemical baseline needs a reference ligand")
    ideal = next((r for r in csv.DictReader(ideal_path.open())
                  if r["selected_as_ideal_binder"] == "True"), None)
    if ideal is None:
        sys.exit(f"no row flagged selected_as_ideal_binder in {ideal_path}")

    sig = args.signature or board["ranked_by"]
    metric = board["ranking_metric"]
    scored = [r for r in board["results"] if r.get("status") == "scored"]
    if not scored:
        sys.exit("no scored rows in the board")

    smiles = {d["struct_id"]: d["smiles"] for d in state["shortlist"]}
    missing = [r["name"] for r in scored if r["struct_id"] not in smiles]
    ids = [r["struct_id"] for r in scored if r["struct_id"] in smiles]

    from representations import REPRESENTATIONS
    fn = REPRESENTATIONS[args.representation]
    mat = fn([smiles[i] for i in ids] + [ideal["smiles"]])
    corpus, query = mat[:-1], mat[-1]
    if getattr(fn, "binary", False):
        inter = corpus @ query
        union = corpus.sum(1) + query.sum() - inter
        union[union == 0] = 1.0
        chem_score = inter / union
    else:
        qn = np.linalg.norm(query) or 1.0
        cn = np.linalg.norm(corpus, axis=1)
        cn[cn == 0] = 1.0
        chem_score = (corpus @ query) / (cn * qn)

    chem_rank = {ids[i]: p for p, i in enumerate(np.argsort(-chem_score), 1)}
    rows = []
    for r in scored:
        sid = r["struct_id"]
        rows.append({
            "struct_id": sid, "name": r["name"], "role": r["role"],
            "structural_rank": r["rank"],
            "structural_score": r["by_signature"][sig][metric],
            "chemical_rank": chem_rank.get(sid),
            "chemical_score": (float(chem_score[ids.index(sid)])
                               if sid in chem_rank else None),
            "rank_delta": (chem_rank[sid] - r["rank"]) if sid in chem_rank else None,
        })

    from scipy.stats import mannwhitneyu
    n = len(ids)
    roles = sorted({r["role"] for r in rows} - {POSITIVE_ROLE})
    tests = {}
    for label, key in (("structural", "structural_score"), ("chemical", "chemical_score")):
        for neg in roles:
            a = np.array([r[key] for r in rows if r["role"] == POSITIVE_ROLE
                          and r[key] is not None])
            b = np.array([r[key] for r in rows if r["role"] == neg
                          and r[key] is not None])
            if not a.size or not b.size:
                continue
            u = mannwhitneyu(a, b, alternative="two-sided")
            tests[f"{label} vs {neg}"] = {
                "auc": float(u.statistic / (a.size * b.size)),
                "p_two_sided": float(u.pvalue),
                "n_positive": int(a.size), "n_negative": int(b.size),
            }

    def med(key):
        v = [r[key] for r in rows if r["role"] == POSITIVE_ROLE and r[key] is not None]
        return float(np.median(v)) if v else None

    out = {
        "what": "structural vs chemical ranking of the SAME drugs",
        "why": ("the board ranks a shortlist; candidates.json ranks the whole approved "
                "corpus. Comparing a rank out of the shortlist against a rank out of the "
                "corpus is not like-for-like, so the chemical baseline is recomputed here "
                "over exactly the drugs the board contains."),
        "population": {"n_drugs": n, "n_excluded_no_smiles": len(missing),
                       "excluded": missing,
                       "roles": {r: sum(1 for x in rows if x["role"] == r)
                                 for r in sorted({x["role"] for x in rows})}},
        "structural": {"signature": sig, "metric": metric},
        "chemical": {"representation": args.representation,
                     "reference_ligand": ideal.get("pref_name") or ideal.get("chembl_id"),
                     "reference_chembl_id": ideal.get("chembl_id"),
                     "note": "similarity to the reference ligand only; no structure used"},
        "median_rank_of_a_known_binder": {"structural": med("structural_rank"),
                                          "chemical": med("chemical_rank"),
                                          "out_of": n},
        "separation": tests,
        "per_drug": sorted(rows, key=lambda r: r["structural_rank"]),
        "caveat": ("site engagement is a structural prioritisation signal. It is not "
                   "evidence of binding, efficacy or clinical suitability."),
    }
    dest = ROOT / "results" / f"ranking_comparison_{args.pipeline}.json"
    dest.write_text(json.dumps(out, indent=1))

    print(f"same-population comparison, n={n} drugs "
          f"({out['population']['roles']})")
    if missing:
        print(f"  excluded for lack of SMILES: {missing}")
    print(f"\nmedian rank of a known binder (of {n}): "
          f"structural {med('structural_rank'):.1f}   chemical {med('chemical_rank'):.1f}")
    print(f"\n{'signal':<26}{'vs':<13}{'AUC':>7}{'p':>10}   n")
    for k, v in tests.items():
        label, neg = k.split(" vs ")
        print(f"{label:<26}{neg:<13}{v['auc']:>7.3f}{v['p_two_sided']:>10.4f}   "
              f"{v['n_positive']}v{v['n_negative']}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
