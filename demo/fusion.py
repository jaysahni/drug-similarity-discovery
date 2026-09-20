"""Rank fusion and multi-query retrieval, scored on the repo's own benchmark task.

Two improvements to similarity, measured against ECFP4 alone on the E1 retrieval
task rather than asserted:

  FUSION        combine several representations by Reciprocal Rank Fusion
  MULTI-QUERY   score a drug by its best similarity to ANY known binder of the
                target, instead of to one arbitrary query molecule

Both exist because of one measured failure. ECFP4 is the best single
representation in this repo's 16-way benchmark, but it ranks by *chemotype*:
argatroban's five nearest approved drugs are all peptidomimetics -- angiotensin
II, icatibant, lisinopril -- and not one is a thrombin drug, even though
argatroban is one. Shape and pharmacophore representations fail differently, so
fusing them should recover co-target drugs that 2D topology alone puts out of
reach.

Reciprocal Rank Fusion rather than score averaging: Tanimoto over ECFP4 and
cosine over z-scored USRCAT are not on a common scale, and normalising them onto
one invents a calibration nobody measured. RRF only uses each list's ORDER, so it
needs no such assumption -- score = sum over representations of 1/(k + rank),
with k=60 as in Cormack et al. 2009.

Nothing here is new statistics: per-query metrics, bootstrap CIs, the paired
Wilcoxon and the Holm correction all come from scripts/metrics.py, so these
numbers sit on the same footing as the benchmark's existing table.

Run:  ./env-kit/bin/python -m demo.fusion --reps morgan usrcat gobbi_pharm2d
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

OUT = REPO / "results" / "demo" / "fusion"
RRF_K = 60  # Cormack, Clarke & Buettcher 2009
REPORTED = ("auroc", "p@1", "ndcg@10", "ef@5pct")


def rank_matrix(sim: np.ndarray) -> np.ndarray:
    """Per-row ranks, 1 = most similar. Ties broken by ascending index.

    The tie-break matters: 2048-bit Tanimoto produces many exact ties, and
    scripts/metrics.py:10-20 documents the same convention for the benchmark.
    """
    order = np.argsort(-sim, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(sim.shape[0])[:, None]
    ranks[rows, order] = np.arange(1, sim.shape[1] + 1)[None, :]
    return ranks


def reciprocal_rank_fusion(sims: list[np.ndarray], k: int = RRF_K) -> np.ndarray:
    """Fuse similarity matrices by rank, not by score.

    Returns a matrix on the RRF scale (higher is better), which is all the
    downstream metrics need -- they only ever rank it.
    """
    if not sims:
        raise ValueError("no similarity matrices to fuse")
    fused = np.zeros_like(sims[0], dtype=float)
    for sim in sims:
        fused += 1.0 / (k + rank_matrix(sim))
    return fused


def multi_query(sim: np.ndarray, rel: np.ndarray) -> np.ndarray:
    """Leave-one-out multi-query: score each item by its best similarity to any
    OTHER known binder of the query's target.

    Single-query retrieval asks "what looks like this one molecule", which lets
    one arbitrary query fix the chemotype. Querying by a group instead should be
    more robust -- and for repurposing it is also the realistic setting, because
    a target usually has several known binders and there is no reason to pick one.

    Leave-one-out is what keeps it honest. The query group for row i is the set
    of i's co-target drugs EXCLUDING i itself, and row i's own column is masked
    out of the result. Without that, an item would be scored by its similarity to
    itself and every number would be meaningless.

    **This is not a fair fight, and the comparison must say so.** Multi-query
    sees the target labels at query time; single-query does not. It should be
    expected to win.

    MEASURED: it loses, and badly. On 300 drugs it is worse than single-query
    ECFP4 at p@1 (-0.158, Holm p = 0.0032), nDCG@10 (-0.074, p = 0.0056) and
    ef@5pct (-1.05, p = 0.026), with AUROC flat. Taking a MAX over the group
    means an item scores high if it resembles ANY group member, so one
    promiscuous or atypical binder drags unrelated molecules up the ranking and
    flattens the top of the list -- which is exactly where p@1 and nDCG@10 look.

    A caveat on that result, because it limits what it licenses: the group here
    is "drugs sharing ANY target with the query", which for a promiscuous drug
    (the corpus maximum is 275 annotated targets) is a large, chemically
    heterogeneous set. A version that queried by the binders of ONE named target
    would be a different and fairer test, and is not what this measures.
    """
    n = sim.shape[0]
    out = np.full_like(sim, -np.inf, dtype=float)
    for i in range(n):
        group = np.flatnonzero(rel[i])          # co-target drugs of i, diagonal already False
        if group.size < 2:
            continue                            # nothing left after leaving one out
        block = sim[group].copy()               # (g, n): each group member vs everything
        # Remove each group member's similarity to ITSELF. Without this, every
        # relevant item scores 1.0 (it is in its own query group) and retrieval
        # comes out perfect -- a spectacular result that measures nothing.
        block[np.arange(group.size), group] = -np.inf
        scores = block.max(axis=0)
        scores[i] = -np.inf                     # and never score the query itself
        out[i] = scores
    return out


def evaluate(sim: np.ndarray, rel: np.ndarray) -> list[dict]:
    """Per-query retrieval metrics, via the repo's own implementation."""
    import metrics

    return metrics.per_query_metrics(sim, rel)


def baseline_name(args) -> str:
    return args.reps[0]


def main() -> None:
    import metrics
    import representations
    from benchmark import load_drugs, relevance_matrix, similarity

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", nargs="+", default=["morgan", "usrcat", "gobbi_pharm2d"],
                        help="representations to fuse; the first is the baseline")
    parser.add_argument("--truth", default="targets", choices=["targets", "moa_targets"])
    parser.add_argument("--limit", type=int, default=0, help="subsample the corpus for a fast check")
    parser.add_argument("--multi-query", action="store_true",
                        help="also evaluate leave-one-out group queries (sees target labels; see docstring)")
    args = parser.parse_args()

    rows = load_drugs(args.truth)
    if args.limit:
        rows = rows[: args.limit]
    smiles = [r["smiles"] for r in rows]
    rel = relevance_matrix(rows)
    print(f"corpus {len(rows)} drugs, truth={args.truth}, "
          f"relevant pairs {int(rel.sum())} ({rel.mean():.4%} of pairs)")

    sims, built = {}, {}
    for name in args.reps:
        fn = representations.REPRESENTATIONS.get(name)
        if fn is None:
            raise SystemExit(f"unknown representation {name!r}")
        t0 = time.time()
        mat = fn(smiles)
        sims[name] = similarity(mat, fn.binary)
        built[name] = round(time.time() - t0, 1)
        n_failed = len(representations.FAILURES.get(name, (0, []))[1]) if name in representations.FAILURES else 0
        print(f"  {name:16s} dim={mat.shape[1]:5d} binary={fn.binary} "
              f"{built[name]:6.1f}s failed={n_failed}")

    sims["RRF_fusion"] = reciprocal_rank_fusion([sims[n] for n in args.reps])
    if args.multi_query:
        sims[f"{baseline_name(args)}_multiquery"] = multi_query(sims[args.reps[0]], rel)
        sims["RRF_fusion_multiquery"] = multi_query(sims["RRF_fusion"], rel)

    per_query = {name: evaluate(sim, rel) for name, sim in sims.items()}
    # per_query_metrics returns dict[metric -> ndarray over queries].
    table = {}
    for name, pq in per_query.items():
        entry = {"n_queries": int(len(pq["auroc"]))}
        for metric in REPORTED:
            mean, lo, hi = metrics.bootstrap_ci(pq[metric])
            entry[metric] = {"mean": round(float(mean), 5),
                             "ci_lo": round(float(lo), 5), "ci_hi": round(float(hi), 5)}
        table[name] = entry

    baseline = args.reps[0]
    comparisons, pvals = {}, {}
    for name in sims:
        if name == baseline:
            continue
        for metric in REPORTED:
            a = np.asarray(per_query[name][metric], dtype=float)
            b = np.asarray(per_query[baseline][metric], dtype=float)
            keep = np.isfinite(a) & np.isfinite(b)
            a, b = a[keep], b[keep]
            boot = metrics.paired_bootstrap(a, b)
            w = metrics.wilcoxon(a, b)
            key = f"{name} - {baseline} [{metric}]"
            comparisons[key] = {
                "delta": round(float(a.mean() - b.mean()), 5),
                "ci": [round(float(x), 5) for x in (boot["ci_lo"], boot["ci_hi"])],
                "p_value": float(f'{boot["p_value"]:.4g}'),
                "wilcoxon_p": float(f'{w["p_value"]:.4g}'),
                "n": int(len(a)),
            }
            pvals[key] = boot["p_value"]

    # holm_bonferroni takes a dict keyed by comparison name and returns one.
    for key, p_holm in metrics.holm_bonferroni(pvals).items():
        comparisons[key]["p_holm"] = float(f"{p_holm:.4g}")

    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "task": "E1 retrieval -- rank every other drug; a hit shares a target",
        "corpus": {"n": len(rows), "truth": args.truth,
                   "relevant_pair_rate": round(float(rel.mean()), 5)},
        "representations_fused": args.reps,
        "fusion": {"method": "Reciprocal Rank Fusion", "k": RRF_K,
                   "why": "Tanimoto and cosine are not on a common scale; RRF uses order only"},
        "build_seconds": built,
        "results": table,
        "comparisons_vs_baseline": comparisons,
        "baseline": baseline,
        "multi_query": {
            "evaluated": bool(args.multi_query),
            "construction": "leave-one-out: query group is the co-target drugs of the "
                            "query, each member's self-similarity removed",
            "not_a_fair_fight": "multi-query sees the target labels at query time and "
                                "single-query does not, so it is EXPECTED to win. The "
                                "number to read is how large the gap is.",
        },
        "caveats": [
            "Ranking quality on a target-sharing task is not the same as finding a "
            "repurposing candidate; a fused ranker that wins here still has to be "
            "checked on the arm-S and arm-P tasks separately.",
            "Representations that fail to parse a SMILES contribute an all-zero row "
            "rather than being dropped, so their similarities to that drug are 0.",
        ],
    }
    path = OUT / f"fusion_{args.truth}.json"
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"\n-> {path.relative_to(REPO)}")

    head = "".join(f"{m:>22s}" for m in REPORTED)
    print(f'\n{"representation":18s}{head}')
    for name, e in table.items():
        mark = "  <- fused" if name == "RRF_fusion" else ""
        cells = "".join(f'{e[m]["mean"]:>9.4f} [{e[m]["ci_lo"]:.3f},{e[m]["ci_hi"]:.3f}]' for m in REPORTED)
        print(f"{name:18s}{cells}{mark}")

    print(f"\nvs {baseline}, Holm-corrected over {len(pvals)} comparisons:")
    for key, c in comparisons.items():
        flag = "*" if c["p_holm"] < 0.05 else " "
        print(f'  {flag} {key:44s} delta {c["delta"]:+.5f}  CI [{c["ci"][0]:+.4f},{c["ci"][1]:+.4f}]  p_holm {c["p_holm"]:.4g}')


if __name__ == "__main__":
    main()
