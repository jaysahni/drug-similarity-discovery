"""Direct matching: rank approved drugs by similarity to a generated candidate.

This is the harness both arms share. Arm P compares a designed peptide to
approved peptides by ESM-C cosine; arm S compares a generated small molecule to
approved small molecules by ECFP4 Tanimoto. The ranking is the easy part. The
three things below are what make a ranking mean anything, and none of them
existed in this repo before.

**The null is the point.** The approved library has thousands of entries, so
*every* query has a nearest neighbour and a bare "0.42 Tanimoto to X" says
nothing. What matters is whether 0.42 is unusual: draw decoy queries, compute
each one's best match against the same library, and report the real query's
nearest-neighbour similarity as a percentile of that distribution. A candidate
sitting at the 50th percentile has found nothing.
`visualization/fixture.py:35` has been shipping an empty slot for exactly this
("No scientific null was computed"); this fills it.

**Novelty, because rediscovery is the failure mode.** PROJECT_GOAL.md rates
"circular rediscovery presented as novel" as high-likelihood and fatal to
credibility (§9), and §6.G3 requires three leaderboards rather than one. Finding
that a thrombin-conditioned peptide resembles bivalirudin is a database lookup
wearing a lab coat unless it is labelled as one.

**Nothing here is an affinity.** These are similarity scores between
representations. No Kd, no potency, no ranking claim beyond "resembles".
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# similarity
# ---------------------------------------------------------------------------


def cosine(query: np.ndarray, corpus: np.ndarray) -> np.ndarray:
    """One-vs-many cosine. Same convention as scripts/match_candidates.py:49."""
    denominator = np.linalg.norm(corpus, axis=1) * np.linalg.norm(query)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(denominator > 0, (corpus @ query) / denominator, 0.0)
    return out


def tanimoto(query: np.ndarray, corpus: np.ndarray) -> np.ndarray:
    """One-vs-many Tanimoto over binary vectors.

    Lifted from scripts/match_candidates.py:41 so the two agree by construction.
    """
    intersection = corpus @ query
    union = corpus.sum(1) + query.sum() - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(union > 0, intersection / union, 0.0)


METRICS = {"cosine": cosine, "tanimoto": tanimoto}


# ---------------------------------------------------------------------------
# the null
# ---------------------------------------------------------------------------


def nn_similarity(query: np.ndarray, corpus: np.ndarray, metric: str) -> float:
    """Best match this query finds anywhere in the corpus."""
    return float(METRICS[metric](query, corpus).max())


def null_distribution(
    decoys: np.ndarray, corpus: np.ndarray, metric: str, *, exclude_self: bool = False
) -> np.ndarray:
    """Nearest-neighbour similarity for each decoy query: the chance baseline.

    `exclude_self=True` when the decoys are themselves rows of `corpus` -- without
    it every decoy matches itself at 1.0 and the null becomes a wall of ones,
    which would make any real candidate look terrible rather than merely
    unremarkable.
    """
    out = np.empty(len(decoys), dtype=float)
    for i, decoy in enumerate(decoys):
        sims = METRICS[metric](decoy, corpus)
        if exclude_self:
            sims[np.argmax(sims)] = -np.inf  # drop the self-match
        out[i] = float(sims.max())
    return out


def percentile_of(value: float, null: np.ndarray) -> float:
    """Fraction of the null the observed value beats. 0.5 means unremarkable.

    Mid-rank, i.e. P(null < v) + 0.5 * P(null == v), not a strict `<`.

    This is not pedantry. Tanimoto over sparse fingerprints produces a heavily
    discrete distribution with many exact ties -- scripts/metrics.py:10-20 warns
    about the same thing for 2048-bit ranking. With a strict `<` those ties are
    all scored against the candidate, and the calibration check in
    demo/test_match_direct.py caught it: KS D=0.187, p<0.001 against uniform,
    mean percentile 0.457 instead of 0.5. Every candidate would have looked less
    remarkable than it is, consistently and invisibly.
    """
    if null.size == 0:
        return float("nan")
    below = float((null < value).mean())
    tied = float((null == value).mean())
    return below + 0.5 * tied


def calibration(null: np.ndarray) -> dict:
    """Summary of the null, so a reader can see its shape rather than trust it."""
    return {
        "n": int(null.size),
        "mean": round(float(null.mean()), 4),
        "sd": round(float(null.std(ddof=1)), 4) if null.size > 1 else None,
        "min": round(float(null.min()), 4),
        "p50": round(float(np.percentile(null, 50)), 4),
        "p95": round(float(np.percentile(null, 95)), 4),
        "max": round(float(null.max()), 4),
    }


# ---------------------------------------------------------------------------
# novelty
# ---------------------------------------------------------------------------


def classify_novelty(targets: str, moa_targets: str, symbol: str) -> str:
    """known_moa / known_offtarget / novel_pairing, per PROJECT_GOAL.md 4.6.

    The repo specifies this field (PROJECT_GOAL.md:380) and renders it
    (visualization/report.py:37), but until now the only thing producing it was a
    test fixture assigning values by list index. Both annotation columns already
    exist in the corpora, so this is set membership, not inference.
    """
    known = {t for t in (targets or "").split(";") if t}
    moa = {t for t in (moa_targets or "").split(";") if t}
    if symbol in moa:
        return "known_moa"
    if symbol in known:
        return "known_offtarget"
    return "novel_pairing"


def leaderboards(rows: list[dict], *, score_key: str) -> dict:
    """all / novel-only / known-only, each ranked independently (§6.G3).

    The known-only board is the sanity check and renders every time: if the
    drugs already known to hit this target do NOT come out near the top of it,
    the method is not working and no novel hit from the same run is believable.
    """
    ordered = sorted(rows, key=lambda r: -r[score_key])
    boards = {
        "all": ordered,
        "known_only": [r for r in ordered if r["novelty"] != "novel_pairing"],
        "novel_only": [r for r in ordered if r["novelty"] == "novel_pairing"],
    }
    return {
        name: [dict(r, rank_in_board=i + 1) for i, r in enumerate(board)]
        for name, board in boards.items()
    }


def enrichment(rows: list[dict], *, score_key: str, top_k: int, symbol: str) -> dict:
    """Enrichment of target-annotated entries in the top k, over their base rate.

    Same idiom as scripts/match_candidates.py:110-112. A random ranker gives 1.0.
    """
    ordered = sorted(rows, key=lambda r: -r[score_key])
    known = [r for r in rows if r["novelty"] != "novel_pairing"]
    if not rows or not known or top_k <= 0:
        return {"enrichment_factor": None, "reason": "no known binders in corpus"}
    hits = sum(1 for r in ordered[:top_k] if r["novelty"] != "novel_pairing")
    base = len(known) / len(rows)
    return {
        "symbol": symbol,
        "top_k": top_k,
        "n_known_in_top_k": hits,
        "n_known_total": len(known),
        "n_corpus": len(rows),
        "base_rate": round(base, 4),
        "enrichment_factor": round((hits / top_k) / base, 3) if base else None,
        "known_ranks": [i + 1 for i, r in enumerate(ordered) if r["novelty"] != "novel_pairing"],
    }
