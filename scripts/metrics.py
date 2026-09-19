"""Per-query retrieval metrics, and the statistics needed to compare them.

benchmark.py's evaluate() returns bare means. That is not enough to say one
representation beats another: the repo's rule is that every headline claim
carries n and, where it is a comparison, a significance test. This module
returns the SAME metrics as PER-QUERY VECTORS, so a mean can be given a
confidence interval and two representations can be compared on the queries
they share.

Ranking conventions, kept identical to benchmark.py so numbers are comparable:
  - the query itself is removed from its own ranking (diagonal of `sim` set to
    -inf, then index i dropped from the order);
  - ties are broken by ascending index via a stable sort. Tanimoto on 2048-bit
    fingerprints produces a great many exact ties, so this matters: sklearn's
    roc_auc_score averages over tied ranks instead, and the two agree exactly
    only when the scores are tie-free. benchmark.py's ceiling() adds ~1e-9 of
    noise for the same reason.
  - only queries with >=1 relevant item are scored, and every returned vector
    has one entry per such query, in ascending query index (see query_index()).

Usage:
    ./env/bin/python scripts/metrics.py      # demo: runs on the repo's own data
    from metrics import per_query_metrics, summarise, paired_bootstrap
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

KS = (1, 5, 10)
BEDROC_ALPHA = 20.0
EF_FRACTIONS = {"ef@1pct": 0.01, "ef@5pct": 0.05}


# --------------------------------------------------------------------------
# single-ranking metrics.  each takes `hits`: a boolean vector in rank order,
# self already removed, True where the retrieved item is relevant.
# --------------------------------------------------------------------------
def auroc(hits):
    """P(relevant ranked above irrelevant), via the Mann-Whitney rank sum.

    NaN when there is nothing to rank against (every item relevant).
    """
    n_pos = int(hits.sum())
    n_neg = hits.size - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    rank_of_pos = np.flatnonzero(hits) + 1
    u = rank_of_pos.sum() - n_pos * (n_pos + 1) / 2.0
    return float(1.0 - u / (n_pos * n_neg))


def average_precision(hits):
    """Mean of precision@rank taken at each relevant item.

    Agrees with sklearn's average_precision_score on tie-free rankings. Tied
    scores are broken by the caller's sort order here, and Tanimoto on 2048-bit
    fingerprints produces a great many exact ties, so per-query values can
    differ from sklearn's tie handling (measured on the real morgan run:
    mean ap 0.2502 here vs 0.2494 from sklearn, per-query deviation up to 0.50).
    """
    n_pos = int(hits.sum())
    if n_pos == 0:
        return float("nan")
    positions = np.flatnonzero(hits) + 1
    return float(np.mean(np.arange(1, n_pos + 1) / positions))


def ndcg(hits, k):
    """nDCG@k with binary gains; ideal ranking puts min(k, n_relevant) hits first."""
    n_pos = int(hits.sum())
    if n_pos == 0:
        return float("nan")
    disc = 1.0 / np.log2(np.arange(2, hits.size + 2))
    dcg = float((hits[:k] * disc[:k]).sum())
    ideal_n = min(k, n_pos, hits.size)
    idcg = float(disc[:ideal_n].sum())
    return dcg / idcg


def enrichment_factor(hits, fraction):
    """Hit rate in the top `fraction` of the ranking over the overall hit rate.

    A random ranker scores 1.0 in expectation. The cut is rounded to the
    nearest item and never smaller than 1, so it stays defined on short lists
    - but that also means every fraction collapses to the same top-1 cut once
    the list is short enough (n <= 29 makes the 1% and 5% cuts both 1 item, so
    ef@1pct and ef@5pct become the same number under two names). Check n
    before reading them as two metrics.
    """
    n_total = hits.size
    n_pos = int(hits.sum())
    if n_pos == 0 or n_total == 0:
        return float("nan")
    n_top = max(1, int(round(fraction * n_total)))
    return float((hits[:n_top].sum() / n_top) / (n_pos / n_total))


def bedroc(hits, alpha=BEDROC_ALPHA):
    """BEDROC: early-recognition-weighted ranking score in [0, 1].

    Truchon & Bayly, "Evaluating Virtual Screening Methods: Good and Bad
    Metrics for the 'Early Recognition' Problem", J. Chem. Inf. Model. 2007,
    47, 488-508 — RIE (eq. 15) rescaled to [0, 1] by eq. 36. alpha=20 puts
    ~80% of the weight in the top 8% of the ranking.

    Note its random-ranking expectation is NOT 0.5: it depends on alpha and on
    the fraction of relevant items, so BEDROC is only comparable across methods
    scored on the same query set.

    Degenerate cases: no relevant items -> 0.0; all items relevant -> NaN.
    The all-relevant case makes eq. 36 a 0/0, and there are no decoys to
    recognise early, so it carries no information. It returns NaN rather than
    1.0 so that bedroc and auroc are averaged over exactly the same set of
    queries - otherwise a corpus with many all-relevant neighbourhoods pulls
    the bedroc mean toward 1 with no record of the different denominator.
    """
    n_total = hits.size
    n_pos = int(hits.sum())
    if n_total == 0:
        return float("nan")
    if n_pos == 0:
        return 0.0
    if n_pos == n_total:
        return float("nan")
    ra = n_pos / n_total
    ranks = np.flatnonzero(hits) + 1.0
    observed = float(np.exp(-alpha * ranks / n_total).sum())
    expected = n_pos * (1.0 - math.exp(-alpha)) / (n_total * math.expm1(alpha / n_total))
    rie = observed / expected
    scale = ra * math.sinh(alpha / 2.0) / (
        math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)
    )
    return float(rie * scale - 1.0 / math.expm1(alpha * (1.0 - ra)))


# --------------------------------------------------------------------------
# per-query evaluation
# --------------------------------------------------------------------------
def query_index(rel):
    """Indices of the queries that per_query_metrics scores, in vector order.

    Validated the same way per_query_metrics validates, because this function
    is the authority on how per-query vectors line up with drugs: a caller that
    got [0,1,2] back from a non-square rel here would silently misalign its
    labels instead of seeing the error.
    """
    rel = np.asarray(rel).astype(bool)
    if rel.ndim != 2 or rel.shape[0] != rel.shape[1]:
        raise ValueError(f"rel must be square (n, n), got {rel.shape}")
    np.fill_diagonal(rel, False)
    return np.flatnonzero(rel.sum(1) > 0)


def per_query_metrics(sim, rel, ks=KS):
    """Every metric as a length-n_queries vector, one entry per scored query.

    sim: (n, n) similarity, higher is better. rel: (n, n) boolean ground truth.
    Returns dict[str, ndarray] keyed p@k, r@k, mrr, auroc, ap, ndcg@k,
    ef@1pct, ef@5pct, bedroc_a20.

    Entries are NaN where a metric is undefined for that query (auroc and
    bedroc when every other item is relevant); summarise() reports the
    affected counts rather than silently averaging over a different
    denominator.

    Denominator note: p@k divides by the number of items actually ranked,
    min(k, n-1), not by k. On this corpus (n=2114, k<=10) the two agree. On a
    corpus smaller than k they do not, and p@k is then inflated relative to
    the divide-by-k convention - so per_query_metrics warns once when any k
    exceeds the ranked-list length rather than reporting it silently.
    """
    sim = np.array(sim, dtype=np.float64)   # copy: the caller's matrix is not touched
    rel = np.asarray(rel).astype(bool)
    if sim.ndim != 2 or sim.shape[0] != sim.shape[1] or sim.shape != rel.shape:
        raise ValueError(f"need square, matching sim/rel; got {sim.shape} and {rel.shape}")

    # self-exclusion, exactly as benchmark.py does it. relevance_matrix()
    # already clears rel's diagonal; clear it again so n_rel cannot count self.
    np.fill_diagonal(sim, -np.inf)
    np.fill_diagonal(rel, False)
    n_rel = rel.sum(1)
    queries = query_index(rel)   # single source of truth for the vector order

    keys = (
        [f"p@{k}" for k in ks]
        + [f"r@{k}" for k in ks]
        + [f"ndcg@{k}" for k in ks]
        + ["mrr", "auroc", "ap", "bedroc_a20"]
        + list(EF_FRACTIONS)
    )
    out = {key: np.full(queries.size, np.nan) for key in keys}

    n_ranked = sim.shape[0] - 1          # self is removed from every ranking
    if n_ranked and max(ks) > n_ranked:
        print(f"warning: k={max(ks)} exceeds the {n_ranked} items ranked per "
              f"query; p@k and ndcg@k divide by {n_ranked}, not by k",
              file=sys.stderr)

    for row, i in enumerate(queries):
        order = np.argsort(-sim[i], kind="stable")
        order = order[order != i]
        hits = rel[i][order]

        for k in ks:
            out[f"p@{k}"][row] = hits[:k].mean()
            out[f"r@{k}"][row] = hits[:k].sum() / n_rel[i]
            out[f"ndcg@{k}"][row] = ndcg(hits, k)

        first = np.flatnonzero(hits)
        out["mrr"][row] = 1.0 / (first[0] + 1) if first.size else 0.0
        out["auroc"][row] = auroc(hits)
        out["ap"][row] = average_precision(hits)
        out["bedroc_a20"][row] = bedroc(hits, BEDROC_ALPHA)
        for name, frac in EF_FRACTIONS.items():
            out[name][row] = enrichment_factor(hits, frac)

    return out


def summarise(per_query):
    """Means over queries, plus n_queries. NaNs are skipped, and any metric
    whose denominator is therefore smaller gets its own n_<metric> key."""
    n_queries = int(len(next(iter(per_query.values())))) if per_query else 0
    out = {"n_queries": n_queries}
    for key, values in per_query.items():
        values = np.asarray(values, dtype=np.float64)
        finite = np.isfinite(values)
        out[key] = float(np.mean(values[finite])) if finite.any() else float("nan")
        if int(finite.sum()) != n_queries:
            out[f"n_{key}"] = int(finite.sum())
    return out


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
# Rows of the resample matrix built at once. Peak memory is ~chunk * n_queries
# floats instead of n_boot * n_queries (measured 331 MB -> 33 MB at n=2069, and bounded
# as the query set grows). Verified: chunking one Generator's draws reproduces
# the single-call stream exactly, so no reported number changes.
_BOOT_CHUNK = 1000


def _resample_means(rng, arrays, n_boot, n):
    """Bootstrap means of each array under a SHARED resampling of the n queries.

    Sharing the row indices across arrays is what makes paired_bootstrap paired.
    """
    out = [np.empty(n_boot) for _ in arrays]
    for start in range(0, n_boot, _BOOT_CHUNK):
        rows = min(_BOOT_CHUNK, n_boot - start)
        idx = rng.integers(0, n, size=(rows, n))
        for dest, values in zip(out, arrays):
            dest[start:start + rows] = values[idx].mean(axis=1)
    return out


def bootstrap_ci(values, n_boot=10000, seed=0, alpha=0.05):
    """Percentile bootstrap over queries. Returns (mean, lo, hi).

    Resamples queries with replacement, so the interval reflects the query
    sample actually benchmarked. NaN entries are dropped first.
    """
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    (means,) = _resample_means(rng, [values], n_boot, values.size)
    lo, hi = np.percentile(means, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return (float(values.mean()), float(lo), float(hi))


def paired_bootstrap(a, b, n_boot=10000, seed=0, alpha=0.05):
    """PAIRED percentile bootstrap on a - b. This is the test to use here.

    a and b are per-query vectors for two representations scored on THE SAME
    queries, in the same order. Pairing matters: per-query scores vary far more
    across queries (a drug with 40 target-sharing neighbours is easy for every
    method, one with 1 is hard for all of them) than between representations,
    so an unpaired test spends nearly all its power on that shared query
    difficulty and will miss a real difference between methods. Resampling
    QUERIES and taking the difference within each resampled query cancels it.

    Returns delta (mean a - mean b), the CI of delta, a two-sided p, and n.
    The p is the bootstrap achieved significance level: the fraction of
    resampled deltas on the far side of zero, doubled, with add-one smoothing
    so it is never exactly 0 (the floor is 2/(n_boot+1)).

    That floor is a resolution limit, not a measurement. At the default
    n_boot=10000 it is 2.0e-4, and a Holm family of 594 comparisons (66 pairs
    x 9 metrics) multiplies the smallest rank by 594, giving adj p 0.119 - so
    every comparison would read as non-significant purely because the
    bootstrap could not resolve a smaller p. When a whole family sits on the
    floor, raise n_boot or use the Wilcoxon p, which has no floor.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"paired vectors must match: {a.shape} vs {b.shape}")
    keep = np.isfinite(a) & np.isfinite(b)  # pairwise-complete
    a, b = a[keep], b[keep]
    if a.size == 0:
        raise ValueError("no queries with a finite score in both vectors")

    rng = np.random.default_rng(seed)
    mean_a, mean_b = _resample_means(rng, [a, b], n_boot, a.size)
    deltas = mean_a - mean_b
    lo, hi = np.percentile(deltas, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    n_le = int((deltas <= 0).sum())
    n_ge = int((deltas >= 0).sum())
    p = 2.0 * min(n_le + 1, n_ge + 1) / (n_boot + 1)
    return {
        "delta": float(a.mean() - b.mean()),
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "p_value": float(min(1.0, p)),
        "n": int(a.size),
    }


def wilcoxon(a, b):
    """Wilcoxon signed-rank on the paired differences - a distribution-free
    cross-check on paired_bootstrap.

    Returns statistic, p_value, n (pairs finite in both vectors) and
    n_effective (pairs that actually enter the test). scipy's default
    zero_method="wilcox" discards zero differences, so for two near-identical
    representations n can be thousands while the test ran on a handful. The
    repo's rule is that every claim carries its n, so report the n that the
    test used, not the n that was offered to it.

    All-zero differences are not a valid input to scipy.stats.wilcoxon, and are
    reported here as statistic 0, p 1: identical vectors are no evidence of a
    difference.
    """
    from scipy import stats

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"paired vectors must match: {a.shape} vs {b.shape}")
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    n_effective = int(np.count_nonzero(a - b))
    if a.size == 0 or n_effective == 0:
        return {"statistic": 0.0, "p_value": 1.0, "n": int(a.size), "n_effective": 0}
    res = stats.wilcoxon(a, b)
    return {"statistic": float(res.statistic), "p_value": float(res.pvalue),
            "n": int(a.size), "n_effective": n_effective}


def holm_bonferroni(pvalues):
    """Holm-Bonferroni step-down adjustment over a family of comparisons.

    A panel of ~12 representations is ~66 pairwise tests, so uncorrected p
    values will manufacture winners. Holm controls the family-wise error rate
    and is uniformly more powerful than plain Bonferroni. Adjusted p values are
    capped at 1 and forced non-decreasing in rank, so they can be read against
    the same threshold as the raw ones.

    Non-finite or out-of-range p values are rejected rather than adjusted. They
    used to pass straight through: max(running, nan) returns running, so a NaN
    came out as adjusted p 0.0 - maximally significant - and shifted every
    later entry's rank multiplier, making the whole family order-dependent. A
    caller that catches paired_bootstrap's ValueError and records NaN would
    then have manufactured a winner at p=0.000.
    """
    for name, p in pvalues.items():
        p = float(p)
        if not math.isfinite(p) or not 0.0 <= p <= 1.0:
            raise ValueError(
                f"p value for {name!r} is {p!r}; holm_bonferroni needs a finite "
                "p in [0, 1]. Drop the comparison or fix it upstream - do not "
                "pass NaN through a correction.")
    items = sorted(pvalues.items(), key=lambda kv: kv[1])
    m = len(items)
    out, running = {}, 0.0
    for rank, (name, p) in enumerate(items):
        running = max(running, (m - rank) * float(p))
        out[name] = min(1.0, running)
    return {name: out[name] for name in pvalues}


# --------------------------------------------------------------------------
# demo: the module applied to the repo's own benchmark, end to end
# --------------------------------------------------------------------------
def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import benchmark as bm

    rows = bm.load_drugs("targets")
    rel = bm.relevance_matrix(rows)
    smiles = [r["smiles"] for r in rows]
    print(f"{len(rows)} drugs, {int(query_index(rel).size)} scorable queries\n")

    out = {"n_drugs": len(rows), "n_queries": int(query_index(rel).size),
           "representations": {}, "comparisons": {}}
    per_query = {}
    for name in ("morgan", "random"):
        fn = bm.REPRESENTATIONS[name]
        sim = bm.similarity(fn(smiles), getattr(fn, "binary", False))
        per_query[name] = per_query_metrics(sim, rel)

    cols = ["p@1", "p@10", "r@10", "mrr", "auroc", "ap", "ndcg@10", "ef@1pct", "bedroc_a20"]
    for name, pq in per_query.items():
        s = summarise(pq)
        print(f"{name}  (n={s['n_queries']})")
        entry = {"n_queries": s["n_queries"]}
        for c in cols:
            mean, lo, hi = bootstrap_ci(pq[c])
            entry[c] = {"mean": mean, "ci_lo": lo, "ci_hi": hi}
            print(f"  {c:<11} {mean:7.4f}  95% CI [{lo:.4f}, {hi:.4f}]")
        out["representations"][name] = entry
        print()

    raw = {}
    for c in cols:
        stat = paired_bootstrap(per_query["morgan"][c], per_query["random"][c])
        wil = wilcoxon(per_query["morgan"][c], per_query["random"][c])
        raw[c] = stat["p_value"]
        out["comparisons"][c] = {"pair": "morgan - random", **stat, "wilcoxon": wil}
        print(f"morgan - random  {c:<11} delta {stat['delta']:+7.4f} "
              f"[{stat['ci_lo']:+.4f}, {stat['ci_hi']:+.4f}] "
              f"p={stat['p_value']:.2e}  wilcoxon p={wil['p_value']:.2e}  "
              f"n={stat['n']} (wilcoxon n={wil['n_effective']})")
    adj = holm_bonferroni(raw)
    print("\nHolm-Bonferroni adjusted p over the " f"{len(raw)} metrics above:")
    for c in cols:
        out["comparisons"][c]["p_holm"] = adj[c]
        print(f"  {c:<11} raw {raw[c]:.2e} -> adj {adj[c]:.2e}")

    dest = Path(__file__).resolve().parent.parent / "results" / "metrics_demo.json"
    dest.parent.mkdir(exist_ok=True)
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
