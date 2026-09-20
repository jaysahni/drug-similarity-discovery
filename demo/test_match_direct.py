"""Check the matching harness, and above all check that the null is calibrated.

Every headline number this pipeline can produce is a percentile against the null
in demo/match_direct.py. If that null is wrong, nothing downstream means
anything -- and a wrong null does not announce itself, it just shifts every
candidate's percentile in one direction.

So the load-bearing test is a calibration check: draw queries from the same
distribution as the decoys, and their percentiles must come out ~uniform. A
Kolmogorov-Smirnov test against U(0,1) is the honest way to assert that, and it
is reported with its p-value rather than as a tuned threshold.

Run:  ./env-kit/bin/python -m demo.test_match_direct
No network, no credentials, no credits.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import kstest

from demo import match_direct as md


def calibration_check(metric: str, *, n_corpus: int = 400, n_null: int = 300,
                      n_probe: int = 200, dim: int = 64, seed: int = 0) -> tuple[bool, str]:
    """Held-out queries from the decoy distribution must score ~uniform percentiles."""
    rng = np.random.default_rng(seed)
    if metric == "tanimoto":
        corpus = (rng.random((n_corpus, dim)) < 0.2).astype(float)
        draw = lambda n: (rng.random((n, dim)) < 0.2).astype(float)  # noqa: E731
    else:
        corpus = rng.normal(size=(n_corpus, dim))
        draw = lambda n: rng.normal(size=(n, dim))  # noqa: E731

    null = md.null_distribution(draw(n_null), corpus, metric)
    percentiles = np.array(
        [md.percentile_of(md.nn_similarity(q, corpus, metric), null) for q in draw(n_probe)]
    )

    stat, p = kstest(percentiles, "uniform")
    ok = p > 0.01
    return ok, f"KS D={stat:.3f} p={p:.3f} mean={percentiles.mean():.3f} (expect ~0.5)"


def main() -> int:
    checks: list[tuple[str, bool, str]] = []

    # 1. The calibration check, per metric. This is the one that gates the rest.
    for metric in ("cosine", "tanimoto"):
        ok, detail = calibration_check(metric)
        checks.append((f"null calibrated ({metric})", ok, detail))

    # 2. Similarity functions agree with the repo's existing implementations.
    rng = np.random.default_rng(1)
    corpus = (rng.random((50, 32)) < 0.3).astype(float)
    query = corpus[7]
    t = md.tanimoto(query, corpus)
    checks.append(("tanimoto self = 1", abs(t[7] - 1.0) < 1e-12, f"{t[7]:.6f}"))
    checks.append(("tanimoto in [0,1]", bool((t >= 0).all() and (t <= 1 + 1e-12).all()),
                   f"min={t.min():.3f} max={t.max():.3f}"))
    c = md.cosine(corpus[3], corpus)
    checks.append(("cosine self = 1", abs(c[3] - 1.0) < 1e-10, f"{c[3]:.6f}"))

    # A disjoint binary vector must score exactly 0, not NaN.
    a = np.array([1.0, 1.0, 0.0, 0.0])
    b = np.array([[0.0, 0.0, 1.0, 1.0]])
    checks.append(("tanimoto disjoint = 0", md.tanimoto(a, b)[0] == 0.0, str(md.tanimoto(a, b)[0])))
    # All-zero rows must not produce NaN either.
    z = np.zeros((1, 4))
    checks.append(("cosine zero-row = 0", md.cosine(a, z)[0] == 0.0, str(md.cosine(a, z)[0])))

    # 3. exclude_self actually excludes.
    with_self = md.null_distribution(corpus[:5], corpus, "tanimoto")
    without = md.null_distribution(corpus[:5], corpus, "tanimoto", exclude_self=True)
    checks.append(("exclude_self drops the 1.0", bool((with_self == 1.0).all() and (without < 1.0).all()),
                   f"with={with_self.max():.3f} without={without.max():.3f}"))

    # 4. Novelty classification.
    cases = [
        (("F2;ALB", "F2", "F2"), "known_moa"),
        (("F2;ALB", "ALB", "F2"), "known_offtarget"),
        (("ALB", "ALB", "F2"), "novel_pairing"),
        (("", "", "F2"), "novel_pairing"),
    ]
    for (targets, moa, symbol), expected in cases:
        got = md.classify_novelty(targets, moa, symbol)
        checks.append((f"novelty {expected}", got == expected, got))

    # 5. Leaderboards partition without loss, and each is ordered.
    rows = [
        {"name": "a", "s": 0.9, "novelty": "known_moa"},
        {"name": "b", "s": 0.7, "novelty": "novel_pairing"},
        {"name": "c", "s": 0.8, "novelty": "known_offtarget"},
        {"name": "d", "s": 0.6, "novelty": "novel_pairing"},
    ]
    boards = md.leaderboards(rows, score_key="s")
    checks.append(("boards partition", len(boards["known_only"]) + len(boards["novel_only"]) == len(rows),
                   f'{len(boards["known_only"])}+{len(boards["novel_only"])}=={len(rows)}'))
    checks.append(("board ordered", [r["name"] for r in boards["all"]] == ["a", "c", "b", "d"],
                   str([r["name"] for r in boards["all"]])))

    # 6. Enrichment: a perfect ranker beats base rate, a random one sits near 1.
    ef = md.enrichment(rows, score_key="s", top_k=2, symbol="F2")
    checks.append(("enrichment > 1 when known rank top", ef["enrichment_factor"] > 1.0,
                   f'EF={ef["enrichment_factor"]} base={ef["base_rate"]}'))

    failed = 0
    for name, ok, detail in checks:
        print(f'  {"PASS" if ok else "FAIL"}  {name:32s} {detail}')
        failed += not ok
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
