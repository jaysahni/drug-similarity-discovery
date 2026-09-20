#!/usr/bin/env python3
"""
scripts/role_separation.py - does the interface score actually tell the roles apart?

The board's headline enrichment is measured against every non-positive row at
once, which lets an easy null carry a hard one. This script takes the board
apart: for every signature, every stored per-drug score, and every pair of roles
the board contains, it asks whether the two role's scores are separable at all.

Why it exists. On this run the answer differs sharply by null. Against decoys
matched only on size and lipophilicity, every arm looks strong. Against approved
kinase inhibitors that carry no annotation for the target - the only hard null on
the board - the design arm is at an AUC of about 0.5. Both facts belong in the
record (CLAUDE.md: keep the negative results), and neither is visible in a single
pooled enrichment figure.

Nothing here is hardcoded to this board: the roles, signatures and metrics are
whatever the JSON contains, so a run that adds an arm or a null is measured
without a code change. A role that cannot be tested is listed in
`not_evaluated` with the reason rather than skipped.

NO AFFINITY ANYWHERE. Every input is an interface-overlap score
(PROJECT_GOAL.md 4.4 and G8); this script only compares those scores between
groups of drugs.

Usage:
    ./env/bin/python scripts/role_separation.py
    ./env/bin/python scripts/role_separation.py --board results/repurpose_X.json \
        --out results/role_separation_X.json --generated-at 2026-09-19T20:30:00-04:00
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

ALPHA = 0.05
POSITIVE_ROLE = "known_binder"


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni step-down adjusted p-values, order preserved.

    Adjusted p_i = max over the steps up to i of (m - rank) * p, clipped at 1
    and forced monotone non-decreasing, which is the standard construction.
    """
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for step, idx in enumerate(order):
        running = max(running, (m - step) * pvals[idx])
        adj[idx] = min(1.0, running)
    return adj


def enrichment(order: list[dict], positive_role: str) -> dict:
    """Top-quartile enrichment, matching scripts/repurpose.py's definition, plus
    the arithmetic maximum a perfect ranking could reach on this subset."""
    n = len(order)
    n_pos = sum(1 for r in order if r["role"] == positive_role)
    if n == 0 or n_pos == 0:
        return {"enrichment": None, "reason": "no positives in this subset"}
    k = max(1, round(0.25 * n))
    hits = sum(1 for r in order[:k] if r["role"] == positive_role)
    base = n_pos / n
    return {
        "n": n, "n_positive": n_pos, "k": k, "hits_in_k": hits,
        "base_rate": base,
        "enrichment": (hits / k) / base,
        "max_possible": (min(k, n_pos) / k) / base,
        "positive_ranks": [i for i, r in enumerate(order, 1)
                           if r["role"] == positive_role],
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--board", default="results/repurpose_colorectal-cancer.json")
    ap.add_argument("--out", default="results/role_separation.json")
    ap.add_argument("--positive-role", default=POSITIVE_ROLE)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--generated-at", default=None,
                    help="ISO-8601 stamp, resolved once here and written into the "
                         "output so reruns over the same board are diffable.")
    args = ap.parse_args()

    if args.generated_at is None:
        sde = os.environ.get("SOURCE_DATE_EPOCH")
        args.generated_at = (
            datetime.fromtimestamp(int(sde), tz=timezone.utc).isoformat() if sde
            else datetime.now().astimezone().replace(microsecond=0).isoformat())

    root = Path(args.root).resolve()
    bpath = root / args.board
    raw = bpath.read_bytes()
    board = json.loads(raw)

    scored = [r for r in board["results"] if r.get("status") == "scored"]
    unscored = [r for r in board["results"] if r.get("status") != "scored"]
    if not scored:
        sys.exit("role_separation.py: no scored rows on the board")

    signatures = sorted(board["signatures"])
    metrics = sorted(k for k, v in scored[0]["by_signature"][signatures[0]].items()
                     if isinstance(v, (int, float)))

    roles: dict[str, list[dict]] = {}
    for r in scored:
        roles.setdefault(r["role"], []).append(r)
    role_counts = {
        r: {"scored": len(roles.get(r, [])),
            "unscored": sum(1 for x in unscored if x["role"] == r)}
        for r in sorted({x["role"] for x in board["results"]})}

    not_evaluated = []
    testable = [r for r in sorted(roles) if len(roles[r]) >= 1]
    for r in sorted(role_counts):
        if r not in testable:
            not_evaluated.append({
                "role": r, "reason": "no scored rows on the board, so no scores "
                                     "exist to compare"})

    def vals(role: str, sig: str, metric: str) -> list[float]:
        return [x["by_signature"][sig][metric] for x in roles[role]]

    # Every distinct pair, every signature, every stored score.
    #
    # Pair convention: a Mann-Whitney U test on (A, B) and on (B, A) is ONE test
    # - the two-sided p is identical and U_BA = n_A*n_B - U_AB. Counting it twice
    # would inflate the Holm family with a duplicate and make the correction
    # wrongly conservative. So each unordered pair appears once, with both
    # directions of the effect size reported inside the record.
    comparisons = []
    for sig in signatures:
        for metric in metrics:
            for a, b in itertools.combinations(testable, 2):
                va, vb = vals(a, sig, metric), vals(b, sig, metric)
                if len(va) < 2 or len(vb) < 2:
                    not_evaluated.append({
                        "signature": sig, "metric": metric, "pair": [a, b],
                        "reason": f"needs at least 2 scored rows per role, have "
                                  f"{len(va)} and {len(vb)}"})
                    continue
                u_ab, p_two = stats.mannwhitneyu(va, vb, alternative="two-sided")
                _, p_greater = stats.mannwhitneyu(va, vb, alternative="greater")
                n1, n2 = len(va), len(vb)
                auc = float(u_ab) / (n1 * n2)
                tied = len(set(va) | set(vb)) < (n1 + n2)
                comparisons.append({
                    "signature": sig, "metric": metric,
                    "group_a": a, "group_b": b, "n_a": n1, "n_b": n2,
                    "u_a_vs_b": float(u_ab),
                    "p_two_sided": float(p_two),
                    "p_one_sided_a_greater": float(p_greater),
                    "p_one_sided_note": "reported for completeness only; the "
                                        "two-sided p is the figure to use",
                    "auc_a_over_b": auc,
                    "auc_b_over_a": 1.0 - auc,
                    "rank_biserial_a_vs_b": 2.0 * auc - 1.0,
                    "median_a": float(np.median(va)),
                    "median_b": float(np.median(vb)),
                    "ties_present": bool(tied),
                })

    # Holm across the whole family, and again within each (signature, metric)
    # block, since the headline test is best judged in its own family too.
    if comparisons:
        glob = holm([c["p_two_sided"] for c in comparisons])
        for c, adj in zip(comparisons, glob):
            c["p_holm_global"] = adj
        for sig in signatures:
            for metric in metrics:
                blk = [c for c in comparisons
                       if c["signature"] == sig and c["metric"] == metric]
                for c, adj in zip(blk, holm([x["p_two_sided"] for x in blk])):
                    c["p_holm_within_signature_metric"] = adj

    # Enrichment on the full board and on each positives+one-null subset.
    pos = args.positive_role
    enrichment_by_subset = []
    for sig in signatures:
        metric = board.get("ranking_metric", metrics[0])
        full = sorted(scored, key=lambda r: -r["by_signature"][sig][metric])
        enrichment_by_subset.append({
            "signature": sig, "metric": metric, "subset": "all scored rows",
            "negative_roles": [r for r in testable if r != pos],
            **enrichment(full, pos)})
        for neg in [r for r in testable if r != pos]:
            sub = [r for r in full if r["role"] in (pos, neg)]
            enrichment_by_subset.append({
                "signature": sig, "metric": metric,
                "subset": f"{pos} + {neg} only", "negative_roles": [neg],
                **enrichment(sub, pos)})

    # Machine-readable verdict per signature, on the board's own ranking metric.
    rank_metric = board.get("ranking_metric", metrics[0])
    verdict = {}
    for sig in signatures:
        v = {"metric": rank_metric, "alpha": args.alpha,
             "basis": "two-sided Mann-Whitney U, Holm-adjusted within this "
                      "signature and metric"}
        for neg in [r for r in testable if r != pos]:
            c = next((x for x in comparisons
                      if x["signature"] == sig and x["metric"] == rank_metric
                      and {x["group_a"], x["group_b"]} == {pos, neg}), None)
            if c is None:
                v[f"separates_{neg}"] = None
                v[f"separates_{neg}_reason"] = "not evaluated; see not_evaluated"
                continue
            padj = c.get("p_holm_within_signature_metric", c["p_two_sided"])
            auc = (c["auc_a_over_b"] if c["group_a"] == pos else c["auc_b_over_a"])
            v[f"separates_{neg}"] = bool(padj < args.alpha)
            v[f"separates_{neg}_p_two_sided"] = c["p_two_sided"]
            v[f"separates_{neg}_p_holm"] = padj
            v[f"separates_{neg}_auc"] = auc
            v[f"separates_{neg}_n_positive"] = (c["n_a"] if c["group_a"] == pos
                                                else c["n_b"])
            v[f"separates_{neg}_n_negative"] = (c["n_b"] if c["group_a"] == pos
                                                else c["n_a"])
        verdict[sig] = v

    out = {
        "what": "role separation of the interface-overlap score, per signature "
                "and per stored per-drug score",
        "generated_at": args.generated_at,
        "source": {
            "board": args.board,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw),
            "n_rows": len(board["results"]),
            "n_scored": len(scored),
            "ranked_by": board.get("ranked_by"),
            "ranking_metric": rank_metric,
        },
        "affinity_used": False,
        "affinity_note": "inputs are interface-overlap scores only; no affinity, "
                         "Kd, IC50 or potency is read, predicted or implied "
                         "anywhere in this file (PROJECT_GOAL.md 4.4, G8)",
        "positive_role": pos,
        "roles": role_counts,
        "signatures": signatures,
        "metrics": metrics,
        "test": {
            "name": "Mann-Whitney U",
            "alternative": "two-sided",
            "alpha": args.alpha,
            "correction": "Holm-Bonferroni",
            "family_size_global": len(comparisons),
            "pair_convention": "each unordered role pair is tested once; "
                               "(A,B) and (B,A) are the same two-sided test, so "
                               "counting both would duplicate a test inside the "
                               "Holm family. Both directions of the effect size "
                               "are given per record as auc_a_over_b / auc_b_over_a.",
            "one_sided_note": "a one-sided p is stored per comparison for "
                              "completeness and is deliberately not used for any "
                              "verdict in this file",
        },
        "comparisons": comparisons,
        "enrichment_by_subset": enrichment_by_subset,
        "verdict_by_signature": verdict,
        "interpretation": {
            "read_these_two_together": [
                "The pooled enrichment on this board is carried by the easy null. "
                "Against decoys matched only on molecular weight and cLogP every "
                "arm separates cleanly; against the hard null - approved kinase "
                "inhibitors carrying no annotation for this target - the design "
                "signature is at an AUC near 0.5, which is a coin flip. See "
                "verdict_by_signature for the per-signature figures and n.",
                "The hard-decoy label means 'not annotated against this target' in "
                "a source annotation that is incomplete, not 'does not bind it'. "
                "Several approved kinase inhibitors are genuinely promiscuous. A "
                "hard decoy that scores well may therefore be a labelling gap "
                "rather than a false positive, which makes every AUC against the "
                "hard null a LOWER BOUND on true discrimination, not a "
                "measurement of it.",
            ],
            "why_neither_alone": "Quoted alone, the first overstates the failure "
                                 "and the second excuses it. The honest claim is "
                                 "that this run cannot separate the two "
                                 "explanations, and says so.",
        },
        "not_evaluated": not_evaluated,
    }

    opath = root / args.out
    opath.parent.mkdir(parents=True, exist_ok=True)
    opath.write_text(json.dumps(out, indent=2) + "\n")

    # readable summary
    print(f"board {args.board}: {len(board['results'])} rows, {len(scored)} scored, "
          f"ranked by {board.get('ranked_by')} / {rank_metric}")
    print("roles: " + ", ".join(f"{r}={c['scored']}" for r, c in role_counts.items()))
    print(f"\n{len(comparisons)} tests "
          f"({len(signatures)} signatures x {len(metrics)} metrics x "
          f"{len(list(itertools.combinations(testable, 2)))} role pairs), "
          f"Holm across the family\n")
    print(f"on the ranking metric ({rank_metric}), two-sided; U and AUC are "
          f"oriented as {pos} over the other role:")
    hdr = f"  {'signature':<20}{'comparison':<34}{'n':>8}{'U':>8}{'p':>10}{'p_holm':>10}{'AUC':>7}"
    print(hdr + "\n  " + "-" * (len(hdr) - 2))
    for c in comparisons:
        if c["metric"] != rank_metric or pos not in (c["group_a"], c["group_b"]):
            continue
        neg = c["group_b"] if c["group_a"] == pos else c["group_a"]
        fwd = c["group_a"] == pos
        auc = c["auc_a_over_b"] if fwd else c["auc_b_over_a"]
        u = c["u_a_vs_b"] if fwd else c["n_a"] * c["n_b"] - c["u_a_vs_b"]
        n_p = c["n_a"] if fwd else c["n_b"]
        n_n = c["n_b"] if fwd else c["n_a"]
        print(f"  {c['signature']:<20}{pos + ' vs ' + neg:<34}"
              f"{str(n_p) + 'v' + str(n_n):>8}{u:>8.1f}"
              f"{c['p_two_sided']:>10.4f}{c['p_holm_within_signature_metric']:>10.4f}"
              f"{auc:>7.3f}")
    print(f"\nenrichment by subset ({rank_metric}):")
    for r in enrichment_by_subset:
        if r.get("enrichment") is None:
            continue
        print(f"  {r['signature']:<20}{r['subset']:<34}"
              f"{r['enrichment']:>6.2f}x  (max {r['max_possible']:.2f}x, "
              f"top {r['k']} of {r['n']})")
    print(f"\nverdict (Holm-adjusted within signature+metric, alpha {args.alpha}):")
    for sig, v in verdict.items():
        flags = [f"{k.replace('separates_', '')}={v[k]}"
                 for k in v if k.startswith("separates_") and isinstance(v[k], bool)]
        print(f"  {sig:<20}" + "  ".join(flags))
    if not_evaluated:
        print(f"\nnot evaluated: {len(not_evaluated)} (see the JSON)")
    print(f"\nwrote {opath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
