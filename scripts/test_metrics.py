"""Correctness tests for scripts/metrics.py, against analytically known answers.

Every case here has an answer that can be derived on paper (perfect ranking,
worst ranking, the expectation of a random ranker) or from an independent
implementation (sklearn's roc_auc_score / average_precision_score / ndcg_score,
and benchmark.py's own evaluate()). Tolerances on the random-ranker tests are
derived from the standard error of the trials actually run, not guessed.

Plain asserts, no pytest (not installed in env/).

Usage:
    ./env/bin/python scripts/test_metrics.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import metrics as M

PASSED = []


def check(name, fn):
    fn()
    PASSED.append(name)
    print(f"  ok   {name}")


# --------------------------------------------------------------------------
# fixtures: 60 drugs in 5 disjoint groups of 12, so every query has exactly
# 11 relevant and 48 irrelevant items among the 59 it ranks.
# --------------------------------------------------------------------------
N_DRUGS, GROUP = 60, 12
N_REL, N_RANKED = GROUP - 1, N_DRUGS - 1


def grouped_rel():
    g = np.arange(N_DRUGS) // GROUP
    rel = g[:, None] == g[None, :]
    np.fill_diagonal(rel, False)
    return rel


def perfect_sim(rel, seed=0):
    """Relevant items strictly first, ties broken by ~1e-9 noise as ceiling() does."""
    rng = np.random.default_rng(seed)
    return rel.astype(np.float64) + rng.uniform(0, 1e-9, rel.shape)


def worst_sim(rel, seed=0):
    rng = np.random.default_rng(seed)
    return -rel.astype(np.float64) + rng.uniform(0, 1e-9, rel.shape)


# --------------------------------------------------------------------------
# 1. perfect ranking
# --------------------------------------------------------------------------
def test_perfect_ranking():
    rel = grouped_rel()
    pq = M.per_query_metrics(perfect_sim(rel), rel)
    s = M.summarise(pq)
    assert s["n_queries"] == N_DRUGS, s["n_queries"]
    # every query has 11 relevant, so the top 10 are all relevant
    for k in (1, 5, 10):
        assert np.allclose(pq[f"p@{k}"], 1.0), (k, pq[f"p@{k}"].min())
        assert np.allclose(pq[f"ndcg@{k}"], 1.0), (k, pq[f"ndcg@{k}"].min())
        assert np.allclose(pq[f"r@{k}"], k / N_REL), (k, pq[f"r@{k}"][0])
    assert np.allclose(pq["mrr"], 1.0)
    assert np.allclose(pq["auroc"], 1.0)
    assert np.allclose(pq["ap"], 1.0)
    assert np.allclose(pq["bedroc_a20"], 1.0), pq["bedroc_a20"].min()
    # EF is capped by the best achievable rate: top-1 is a hit for every query,
    # so ef@1pct == 1 / (11/59) == N_RANKED / N_REL
    assert np.allclose(pq["ef@1pct"], N_RANKED / N_REL)


# --------------------------------------------------------------------------
# 2. worst ranking
# --------------------------------------------------------------------------
def test_worst_ranking():
    rel = grouped_rel()
    pq = M.per_query_metrics(worst_sim(rel), rel)
    assert np.allclose(pq["auroc"], 0.0), pq["auroc"].max()
    assert np.allclose(pq["bedroc_a20"], 0.0, atol=1e-9), pq["bedroc_a20"].max()
    for k in (1, 5, 10):
        assert np.allclose(pq[f"p@{k}"], 0.0)
        assert np.allclose(pq[f"ndcg@{k}"], 0.0)
    # first hit sits at rank 49 (48 irrelevant items ahead of it)
    n_neg = N_RANKED - N_REL
    assert np.allclose(pq["mrr"], 1.0 / (n_neg + 1)), pq["mrr"][0]
    assert np.allclose(pq["ef@1pct"], 0.0)
    positions = np.arange(n_neg + 1, N_RANKED + 1)
    assert np.allclose(pq["ap"], np.mean(np.arange(1, N_REL + 1) / positions))


# --------------------------------------------------------------------------
# 3. random ranking: AUROC -> 0.5 and EF -> 1.0 in expectation
# --------------------------------------------------------------------------
def _random_trials(n_trials=200):
    rel = grouped_rel()
    rng = np.random.default_rng(1234)
    means = {"auroc": [], "ef@1pct": [], "ef@5pct": [], "bedroc_a20": []}
    for _ in range(n_trials):
        sim = rng.standard_normal((N_DRUGS, N_DRUGS))
        pq = M.per_query_metrics(sim, rel)
        for key in means:
            means[key].append(float(np.mean(pq[key])))
    return {k: np.asarray(v) for k, v in means.items()}, n_trials


RANDOM_TRIALS, N_TRIALS = _random_trials()


def _within_4se(values, expected, label):
    """4 standard errors of the trial means: a ~1-in-16000 two-sided false alarm."""
    se = values.std(ddof=1) / np.sqrt(values.size)
    err = abs(values.mean() - expected)
    assert err < 4 * se, f"{label}: mean {values.mean():.5f} vs {expected}, err {err:.5f} >= 4se {4*se:.5f}"
    print(f"       {label}: mean {values.mean():.5f}  expected {expected}  "
          f"|err| {err:.5f}  4se {4 * se:.5f}  ({values.size} trials)")


def test_random_auroc_is_half():
    _within_4se(RANDOM_TRIALS["auroc"], 0.5, "random auroc")


def test_random_enrichment_is_one():
    _within_4se(RANDOM_TRIALS["ef@1pct"], 1.0, "random ef@1pct")
    _within_4se(RANDOM_TRIALS["ef@5pct"], 1.0, "random ef@5pct")


# --------------------------------------------------------------------------
# 4. cross-checks against sklearn (tie-free scores, so ranking is unambiguous)
# --------------------------------------------------------------------------
def _random_case(n=80, seed=7):
    rng = np.random.default_rng(seed)
    sim = rng.random((n, n))
    rel = rng.random((n, n)) < 0.12
    np.fill_diagonal(rel, False)
    for i in range(n):  # guarantee >=1 relevant and >=1 irrelevant per query
        rel[i, (i + 1) % n] = True
        rel[i, (i + 2) % n] = False
    return sim, rel


def test_auroc_matches_sklearn():
    from sklearn.metrics import roc_auc_score

    sim, rel = _random_case()
    pq = M.per_query_metrics(sim, rel)
    worst = 0.0
    for row, i in enumerate(M.query_index(rel)):
        mask = np.arange(sim.shape[0]) != i
        ref = roc_auc_score(rel[i][mask], sim[i][mask])
        worst = max(worst, abs(pq["auroc"][row] - ref))
    assert worst < 1e-9, worst
    print(f"       max |auroc - sklearn| = {worst:.3e} over {len(M.query_index(rel))} queries")


def test_ap_matches_sklearn():
    from sklearn.metrics import average_precision_score

    sim, rel = _random_case()
    pq = M.per_query_metrics(sim, rel)
    worst = 0.0
    for row, i in enumerate(M.query_index(rel)):
        mask = np.arange(sim.shape[0]) != i
        ref = average_precision_score(rel[i][mask], sim[i][mask])
        worst = max(worst, abs(pq["ap"][row] - ref))
    assert worst < 1e-9, worst
    print(f"       max |ap - sklearn| = {worst:.3e}")


def test_ndcg_matches_sklearn():
    from sklearn.metrics import ndcg_score

    sim, rel = _random_case()
    pq = M.per_query_metrics(sim, rel)
    worst = 0.0
    for row, i in enumerate(M.query_index(rel)):
        mask = np.arange(sim.shape[0]) != i
        for k in (1, 5, 10):
            ref = ndcg_score([rel[i][mask].astype(float)], [sim[i][mask]], k=k)
            worst = max(worst, abs(pq[f"ndcg@{k}"][row] - ref))
    assert worst < 1e-9, worst
    print(f"       max |ndcg@k - sklearn| = {worst:.3e}")


def test_matches_benchmark_evaluate():
    """The means must reproduce benchmark.py's evaluate() exactly."""
    import benchmark as bm

    sim, rel = _random_case()
    ref = bm.evaluate(sim, rel)
    got = M.summarise(M.per_query_metrics(sim, rel))
    assert got["n_queries"] == ref["n_queries"], (got["n_queries"], ref["n_queries"])
    for key in ["p@1", "p@5", "p@10", "r@1", "r@5", "r@10", "mrr", "auroc"]:
        assert abs(got[key] - ref[key]) < 1e-12, (key, got[key], ref[key])
    print(f"       reproduces benchmark.evaluate() on 8 metrics, n={ref['n_queries']}")


# --------------------------------------------------------------------------
# 5. BEDROC edge cases
# --------------------------------------------------------------------------
def test_bedroc_edges():
    # all-relevant is NaN, not 1.0: there are no decoys to recognise early, so
    # the query carries no information. Matching auroc here keeps the two means
    # over the same set of queries instead of quietly different denominators.
    assert np.isnan(M.bedroc(np.ones(50, dtype=bool)))
    assert M.bedroc(np.zeros(50, dtype=bool)) == 0.0       # nothing to find
    assert np.isnan(M.bedroc(np.zeros(0, dtype=bool)))     # empty ranking
    best = np.zeros(200, dtype=bool); best[:20] = True
    worst = np.zeros(200, dtype=bool); worst[-20:] = True
    assert abs(M.bedroc(best) - 1.0) < 1e-9, M.bedroc(best)
    assert abs(M.bedroc(worst) - 0.0) < 1e-9, M.bedroc(worst)
    # early recognition: the same 20 hits, moved earlier, must score higher
    mid = np.zeros(200, dtype=bool); mid[90:110] = True
    assert M.bedroc(worst) < M.bedroc(mid) < M.bedroc(best)
    print(f"       bedroc best={M.bedroc(best):.6f} mid={M.bedroc(mid):.6f} "
          f"worst={M.bedroc(worst):.6f}")


def test_undefined_metrics_are_nan_not_dropped():
    """A query where every other item is relevant has no AUROC; it must be NaN
    and summarise() must say the denominator shrank."""
    rel = np.ones((6, 6), dtype=bool)
    sim = np.random.default_rng(0).random((6, 6))
    pq = M.per_query_metrics(sim, rel)
    assert np.all(np.isnan(pq["auroc"])), pq["auroc"]
    assert np.all(np.isnan(pq["bedroc_a20"])), pq["bedroc_a20"]
    s = M.summarise(pq)
    assert s["n_queries"] == 6 and s["n_auroc"] == 0 and np.isnan(s["auroc"])
    assert s["n_bedroc_a20"] == 0 and np.isnan(s["bedroc_a20"])


# --------------------------------------------------------------------------
# 6. statistics
# --------------------------------------------------------------------------
def test_paired_bootstrap_identical():
    rng = np.random.default_rng(3)
    a = rng.random(120)
    out = M.paired_bootstrap(a, a.copy(), n_boot=2000, seed=0)
    assert out["delta"] == 0.0, out
    assert out["p_value"] == 1.0, out
    assert out["ci_lo"] == 0.0 and out["ci_hi"] == 0.0, out
    assert out["n"] == 120
    w = M.wilcoxon(a, a.copy())
    assert w["p_value"] == 1.0, w


def test_paired_bootstrap_clear_difference():
    rng = np.random.default_rng(5)
    base = rng.random(200)
    a = base + 0.20 + rng.normal(0, 0.01, 200)   # a beats b on essentially every query
    b = base
    out = M.paired_bootstrap(a, b, n_boot=10000, seed=0)
    assert out["p_value"] < 0.01, out
    assert abs(out["delta"] - (a.mean() - b.mean())) < 1e-12
    assert out["ci_lo"] > 0, out
    w = M.wilcoxon(a, b)
    assert w["p_value"] < 0.01, w
    print(f"       delta {out['delta']:.4f} CI [{out['ci_lo']:.4f}, {out['ci_hi']:.4f}] "
          f"p={out['p_value']:.2e}  wilcoxon p={w['p_value']:.2e}")


def test_paired_bootstrap_beats_unpaired_on_correlated_data():
    """Why pairing matters: a small consistent edge buried in huge query-to-query
    variance is invisible unless the shared difficulty is cancelled out."""
    rng = np.random.default_rng(11)
    difficulty = rng.normal(0, 1.0, 300)      # dominates both vectors
    a = difficulty + 0.05 + rng.normal(0, 0.02, 300)
    b = difficulty + rng.normal(0, 0.02, 300)
    paired = M.paired_bootstrap(a, b, n_boot=10000, seed=0)
    # unpaired: resample each vector independently
    r = np.random.default_rng(0)
    deltas = a[r.integers(0, 300, (10000, 300))].mean(1) - b[r.integers(0, 300, (10000, 300))].mean(1)
    p_unpaired = min(1.0, 2.0 * min((deltas <= 0).sum() + 1, (deltas >= 0).sum() + 1) / 10001)
    assert paired["p_value"] < 0.01 < p_unpaired, (paired["p_value"], p_unpaired)
    print(f"       paired p={paired['p_value']:.2e}  unpaired p={p_unpaired:.3f}")


def test_bootstrap_ci_covers_true_mean():
    """Nominal 95% coverage, measured over independent samples from a known mean."""
    true_mean, n_rep, n_sample = 0.3, 200, 200
    rng = np.random.default_rng(42)
    covered = 0
    for rep in range(n_rep):
        sample = rng.normal(true_mean, 0.1, n_sample)
        _, lo, hi = M.bootstrap_ci(sample, n_boot=1500, seed=rep)
        covered += lo <= true_mean <= hi
    rate = covered / n_rep
    se = np.sqrt(0.95 * 0.05 / n_rep)         # 4se ~= 0.062 around the nominal 0.95
    assert abs(rate - 0.95) < 4 * se, (rate, 4 * se)
    mean, lo, hi = M.bootstrap_ci(rng.normal(true_mean, 0.1, 500), n_boot=10000)
    assert lo < true_mean < hi and lo < mean < hi
    print(f"       coverage {rate:.3f} over {n_rep} samples (nominal 0.950, 4se {4*se:.3f})")


def test_bootstrap_ci_ignores_nan():
    v = np.array([1.0, 2.0, 3.0, np.nan])
    mean, lo, hi = M.bootstrap_ci(v, n_boot=1000)
    assert abs(mean - 2.0) < 1e-12 and lo <= 2.0 <= hi
    assert all(np.isnan(x) for x in M.bootstrap_ci(np.array([np.nan, np.nan])))


def test_holm_bonferroni():
    # worked by hand: m=4, sorted .001 .008 .039 .041
    #   .001*4=.004 | .008*3=.024 | .039*2=.078 | .041*1=.041 -> forced up to .078
    adj = M.holm_bonferroni({"a": 0.001, "b": 0.008, "c": 0.039, "d": 0.041})
    assert abs(adj["a"] - 0.004) < 1e-12, adj
    assert abs(adj["b"] - 0.024) < 1e-12, adj
    assert abs(adj["c"] - 0.078) < 1e-12, adj
    assert abs(adj["d"] - 0.078) < 1e-12, adj        # monotonicity, not 0.041
    assert list(adj) == ["a", "b", "c", "d"]          # input order preserved
    assert M.holm_bonferroni({"x": 0.5, "y": 0.9}) == {"x": 1.0, "y": 1.0}  # capped at 1
    assert M.holm_bonferroni({}) == {}


def test_summarise_shapes():
    rel = grouped_rel()
    pq = M.per_query_metrics(perfect_sim(rel), rel, ks=(1, 3))
    assert set(pq) == {"p@1", "p@3", "r@1", "r@3", "ndcg@1", "ndcg@3", "mrr",
                       "auroc", "ap", "bedroc_a20", "ef@1pct", "ef@5pct"}, set(pq)
    assert all(v.shape == (N_DRUGS,) for v in pq.values())
    assert M.summarise(pq)["n_queries"] == N_DRUGS
    assert np.array_equal(M.query_index(rel), np.arange(N_DRUGS))


def test_queries_without_relevant_items_are_excluded():
    rel = grouped_rel().copy()
    rel[0, :] = False          # drug 0 now has no relevant neighbour
    rel[:, 0] = False
    sim = np.random.default_rng(0).random((N_DRUGS, N_DRUGS))
    pq = M.per_query_metrics(sim, rel)
    assert M.summarise(pq)["n_queries"] == N_DRUGS - 1
    assert np.array_equal(M.query_index(rel), np.arange(1, N_DRUGS))


def test_self_is_excluded():
    """A representation that ranks the query itself first gains nothing from it."""
    rel = grouped_rel()
    sim = perfect_sim(rel)
    np.fill_diagonal(sim, 1e6)      # self scores far above everything
    pq = M.per_query_metrics(sim, rel)
    assert np.allclose(pq["p@1"], 1.0) and np.allclose(pq["auroc"], 1.0)
    # and a True diagonal in rel must not inflate n_rel
    rel_self = rel.copy()
    np.fill_diagonal(rel_self, True)
    assert np.allclose(M.per_query_metrics(sim, rel_self)["r@10"], 10 / N_REL)


def test_bad_shapes_raise():
    for fn, args in [
        (M.per_query_metrics, (np.zeros((3, 4)), np.zeros((3, 4), bool))),
        (M.paired_bootstrap, (np.zeros(3), np.zeros(4))),
        (M.wilcoxon, (np.zeros(3), np.zeros(4))),
    ]:
        try:
            fn(*args)
        except ValueError:
            continue
        raise AssertionError(f"{fn.__name__} accepted mismatched shapes")


def test_holm_rejects_non_finite_p():
    """Regression: a NaN p used to come out of Holm as adjusted p 0.0.

    max(running, nan) returns running, so the NaN entry was assigned the
    running maximum (0.0 at rank 0) and every later entry's rank multiplier
    shifted - the family became order-dependent and manufactured a winner at
    p=0.000. Non-finite and out-of-range p must be rejected instead.
    """
    for bad in ({"a": float("nan"), "b": 0.01}, {"a": -0.1, "b": 0.01},
                {"a": 1.5, "b": 0.01}, {"a": float("inf"), "b": 0.01}):
        try:
            M.holm_bonferroni(bad)
        except ValueError:
            continue
        raise AssertionError(f"holm_bonferroni accepted {bad}")
    # and the order-dependence it caused is gone for valid families
    fam = {"a": 0.04, "b": 0.01, "c": 0.5}
    assert M.holm_bonferroni(fam) == M.holm_bonferroni(dict(reversed(list(fam.items()))))


def test_wilcoxon_reports_effective_n():
    """scipy drops zero differences, so n and the test's real n can diverge."""
    a = np.zeros(10)
    b = np.zeros(10); b[0] = 1.0
    out = M.wilcoxon(a, b)
    assert out["n"] == 10, out
    assert out["n_effective"] == 1, out
    same = M.wilcoxon(a, a.copy())
    assert same["n_effective"] == 0 and same["p_value"] == 1.0, same


def test_query_index_validates_shape():
    """query_index is the authority on vector ordering, so a non-square rel
    must raise rather than return indices that per_query_metrics would reject."""
    try:
        M.query_index(np.ones((3, 4), dtype=bool))
    except ValueError:
        return
    raise AssertionError("query_index accepted a non-square rel")


TESTS = [
    ("perfect ranking -> p@k=1, mrr=1, auroc=1, ap=1, ndcg=1, bedroc=1", test_perfect_ranking),
    ("worst ranking -> auroc=0, bedroc=0, mrr=1/49", test_worst_ranking),
    ("random ranking -> auroc ~ 0.5", test_random_auroc_is_half),
    ("random ranking -> enrichment factor ~ 1.0", test_random_enrichment_is_one),
    ("auroc == sklearn.roc_auc_score", test_auroc_matches_sklearn),
    ("ap == sklearn.average_precision_score", test_ap_matches_sklearn),
    ("ndcg == sklearn.ndcg_score", test_ndcg_matches_sklearn),
    ("means == benchmark.evaluate()", test_matches_benchmark_evaluate),
    ("bedroc edge cases (all-relevant, none-relevant, empty)", test_bedroc_edges),
    ("undefined metrics are NaN and counted, not dropped", test_undefined_metrics_are_nan_not_dropped),
    ("paired_bootstrap on identical input -> delta 0, p 1", test_paired_bootstrap_identical),
    ("paired_bootstrap on a clear difference -> p < 0.01", test_paired_bootstrap_clear_difference),
    ("pairing recovers a signal an unpaired test misses", test_paired_bootstrap_beats_unpaired_on_correlated_data),
    ("bootstrap_ci covers the true mean at ~95%", test_bootstrap_ci_covers_true_mean),
    ("bootstrap_ci drops NaN", test_bootstrap_ci_ignores_nan),
    ("holm_bonferroni matches a hand-worked family", test_holm_bonferroni),
    ("summarise / query_index shapes and keys", test_summarise_shapes),
    ("queries with no relevant item are excluded", test_queries_without_relevant_items_are_excluded),
    ("self-match is excluded from the ranking and from n_rel", test_self_is_excluded),
    ("mismatched shapes raise ValueError", test_bad_shapes_raise),
    ("holm_bonferroni rejects non-finite p (regression)", test_holm_rejects_non_finite_p),
    ("wilcoxon reports the n the test actually used", test_wilcoxon_reports_effective_n),
    ("query_index validates its shape", test_query_index_validates_shape),
]


def main():
    print(f"scripts/test_metrics.py - {len(TESTS)} tests\n")
    failed = []
    for name, fn in TESTS:
        try:
            check(name, fn)
        except Exception as exc:  # noqa: BLE001 - a crashing test is a failing
            # test, not a reason to abandon the remaining ones. Catching only
            # AssertionError here meant an ImportError in test 2 silently took
            # tests 3-20 with it and printed no summary at all.
            failed.append((name, exc))
            print(f"  FAIL {name}\n       {exc}")
    print(f"\n{len(PASSED)} passed, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
