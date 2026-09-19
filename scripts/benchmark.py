"""Retrieval benchmark for drug-similarity representations.

Task: given a query drug, rank every other drug by similarity. A retrieved drug
counts as relevant if it shares a protein target with the query.

A representation is any function  list[smiles] -> ndarray (n, d)  registered in
scripts/representations.py. Similarity is cosine on that matrix, except for
fingerprint representations which declare `binary=True` and use Tanimoto.

Every number carries n and a 95% bootstrap CI, and every comparison carries a
paired test (bootstrap + Wilcoxon) with a Holm-Bonferroni correction over the
family -- a panel this wide will otherwise manufacture a winner. Two references
frame the table: a random-ranking FLOOR and the perfect-ranking CEILING, so a
score reads as a fraction of the achievable range rather than as a bare number.

Usage:
    ./env/bin/python scripts/benchmark.py                    # every representation
    ./env/bin/python scripts/benchmark.py morgan maccs       # a subset
    ./env/bin/python scripts/benchmark.py --truth moa_targets
    ./env/bin/python scripts/benchmark.py --quick            # skip the slow ones
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import metrics as M                      # noqa: E402
from representations import REPRESENTATIONS, FAILURES  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "drugs.csv"
RESULTS = ROOT / "results"
KS = (1, 5, 10)
SEED = 0

# Representations whose cost is dominated by conformer generation or a neural
# forward pass. --quick drops them; they are cached, so the second run is fast.
SLOW = {"usrcat", "chemberta", "gobbi_pharm2d"}

# The comparison every other representation is measured against. ECFP4 is the
# standard chemical-similarity baseline, so "beats morgan" is the claim that matters.
BASELINE = "morgan"
FLOOR = "random"


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def load_drugs(truth):
    """Drugs with a structure and at least one target under the chosen ground truth."""
    rows = []
    with DATA.open() as fh:
        for r in csv.DictReader(fh):
            labels = [t for t in r[truth].split(";") if t]
            if r["smiles"] and labels:
                r["labels"] = labels
                rows.append(r)
    return rows


def relevance_matrix(rows):
    """Boolean (n, n): True where i and j share a target. Diagonal False."""
    labels = sorted({t for r in rows for t in r["labels"]})
    idx = {t: i for i, t in enumerate(labels)}
    memb = np.zeros((len(rows), len(labels)), dtype=bool)
    for i, r in enumerate(rows):
        for t in r["labels"]:
            memb[i, idx[t]] = True
    rel = memb.astype(np.int16) @ memb.astype(np.int16).T > 0
    np.fill_diagonal(rel, False)
    return rel


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def similarity(mat, binary):
    if binary:
        # Tanimoto for bit vectors: |a&b| / (|a|+|b|-|a&b|)
        inter = mat @ mat.T
        card = mat.sum(1)
        union = card[:, None] + card[None, :] - inter
        union[union == 0] = 1.0
        return inter / union
    norm = np.linalg.norm(mat, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    unit = mat / norm
    return unit @ unit.T


def evaluate(sim, rel):
    """Mean p@k, r@k, MRR and AUROC. Kept for callers that predate metrics.py."""
    pq = M.per_query_metrics(sim, rel, ks=KS)
    s = M.summarise(pq)
    return {k: s[k] for k in
            ["n_queries"] + [f"p@{k}" for k in KS] + [f"r@{k}" for k in KS] + ["mrr", "auroc"]}


def ceiling(rel):
    """Perfect ranking: all relevant drugs first. The maximum any method could score."""
    n = rel.shape[0]
    sim = rel.astype(np.float64) + np.random.default_rng(SEED).uniform(0, 1e-9, (n, n))
    return M.summarise(M.per_query_metrics(sim, rel, ks=KS))


# --------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    truth = "targets"
    quick = "--quick" in argv
    if quick:
        argv.remove("--quick")
    if "--truth" in argv:
        i = argv.index("--truth")
        truth = argv[i + 1]
        del argv[i:i + 2]
    if truth not in ("targets", "moa_targets"):
        sys.exit(f"--truth must be 'targets' or 'moa_targets', got {truth!r}")

    wanted = argv or [n for n in REPRESENTATIONS if not (quick and n in SLOW)]
    unknown = [w for w in wanted if w not in REPRESENTATIONS]
    if unknown:
        sys.exit(f"unknown representation(s): {unknown}\nknown: {list(REPRESENTATIONS)}")

    rows = load_drugs(truth)
    rel = relevance_matrix(rows)
    n_rel = rel.sum(1)
    print(f"ground truth: shared {truth}")
    print(f"{len(rows)} drugs with a structure and >=1 target")
    print(f"{int((n_rel > 0).sum())} have >=1 target-sharing neighbour "
          f"(median {int(np.median(n_rel[n_rel > 0]))} neighbours)")
    print(f"random-pair relevance rate: {rel.mean():.4f}\n")

    smiles = [r["smiles"] for r in rows]
    per_query, summary = {}, {}
    for name in wanted:
        fn = REPRESENTATIONS[name]
        t0 = time.time()
        print(f"computing {name}...", end=" ", flush=True)
        mat = fn(smiles)
        pq = M.per_query_metrics(similarity(mat, getattr(fn, "binary", False)), rel, ks=KS)
        per_query[name] = pq
        s = M.summarise(pq)
        n_failed = FAILURES.get(name, (0, []))[0]
        # a row of zeros retrieves nothing and scores 0 on every metric; that is an
        # encoding failure, not a miss, so it is counted separately either way
        n_empty = int((np.abs(mat).sum(1) == 0).sum())
        s.update({"dim": int(mat.shape[1]), "binary": bool(getattr(fn, "binary", False)),
                  "blurb": getattr(fn, "blurb", ""), "seconds": round(time.time() - t0, 1),
                  "n_failed_to_encode": n_failed, "n_all_zero_rows": n_empty})
        for k in ("p@1", "p@10", "auroc", "ap", "bedroc_a20"):
            mean, lo, hi = M.bootstrap_ci(pq[k], n_boot=2000, seed=SEED)
            s[f"{k}_ci"] = [lo, hi]
        summary[name] = s
        print(f"{s['seconds']}s  p@1={s['p@1']:.4f}  (failed {n_failed}, all-zero {n_empty})")

    # --- comparisons, corrected over the whole family --------------------
    comparisons, raw_p = {}, {}
    for ref in (FLOOR, BASELINE):
        if ref not in per_query:
            continue
        for name in wanted:
            if name == ref:
                continue
            key = f"{name} - {ref}"
            a, b = per_query[name]["p@1"], per_query[ref]["p@1"]
            st = M.paired_bootstrap(a, b, n_boot=10000, seed=SEED)
            st["wilcoxon"] = M.wilcoxon(a, b)
            comparisons[key] = st
            raw_p[key] = st["p_value"]
    if raw_p:
        for key, adj in M.holm_bonferroni(raw_p).items():
            comparisons[key]["p_holm"] = adj

    results = {
        "_meta": {"truth": truth, "n_drugs": len(rows),
                  "n_queries": int((n_rel > 0).sum()),
                  "random_pair_relevance": float(rel.mean()),
                  "primary_metric": "p@1", "baseline": BASELINE, "floor": FLOOR,
                  "ks": list(KS), "seed": SEED,
                  "family_size_for_holm": len(raw_p)},
        "_ceiling": ceiling(rel),
        "representations": summary,
        "comparisons": comparisons,
    }
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"benchmark_{truth}.json"
    out.write_text(json.dumps(results, indent=2))

    # --- table ------------------------------------------------------------
    cols = ["p@1", "p@5", "p@10", "r@10", "mrr", "auroc", "ap", "bedroc_a20"]
    order = sorted((n for n in wanted if n != FLOOR),
                   key=lambda n: -summary[n]["p@1"])
    order = ([FLOOR] if FLOOR in wanted else []) + order
    width = max(len(n) for n in wanted) + 2
    print(f"\n{'method'.ljust(width)}" + "".join(c.rjust(9) for c in cols) + "   p@1 95% CI")
    print("-" * (width + 9 * len(cols) + 20))
    for name in order:
        r = summary[name]
        lo, hi = r["p@1_ci"]
        print(name.ljust(width) + "".join(f"{r[c]:9.4f}" for c in cols)
              + f"   [{lo:.4f}, {hi:.4f}]")
    c = results["_ceiling"]
    print("_ceiling".ljust(width) + "".join(f"{c[col]:9.4f}" for col in cols))

    if BASELINE in per_query:
        print(f"\nvs {BASELINE} on p@1 (Holm-adjusted over {len(raw_p)} comparisons):")
        for name in order:
            key = f"{name} - {BASELINE}"
            if key not in comparisons:
                continue
            st = comparisons[key]
            flag = "" if st["p_holm"] < 0.05 else "  (n.s.)"
            print(f"  {name.ljust(width)} delta {st['delta']:+.4f} "
                  f"[{st['ci_lo']:+.4f}, {st['ci_hi']:+.4f}]  "
                  f"Holm p={st['p_holm']:.2e}{flag}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
