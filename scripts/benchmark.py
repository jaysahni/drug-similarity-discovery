"""Retrieval benchmark for drug-similarity representations.

Task: given a query drug, rank every other drug by similarity. A retrieved drug
counts as relevant if it shares a protein target with the query.

A representation is any function  list[smiles] -> ndarray (n, d)  registered in
REPRESENTATIONS. Similarity is cosine on that matrix, except for fingerprint
representations which declare `binary=True` and use Tanimoto.

Reports precision@k, recall@k, MRR and AUROC alongside a random-ranking floor
and the perfect-ranking ceiling, so a score can be read as a fraction of the
achievable range rather than as a bare number.

Usage:
    ./env/bin/python scripts/benchmark.py
    ./env/bin/python scripts/benchmark.py --truth moa_targets
    ./env/bin/python scripts/benchmark.py morgan random
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "drugs.csv"
RESULTS = ROOT / "results"
KS = (1, 5, 10)
SEED = 0


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
# representations
# --------------------------------------------------------------------------
REPRESENTATIONS = {}


def representation(name, binary=False):
    def deco(fn):
        fn.binary = binary
        REPRESENTATIONS[name] = fn
        return fn

    return deco


@representation("morgan", binary=True)
def morgan(smiles):
    """Morgan (ECFP4) 2048-bit fingerprints - the baseline every model must beat."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        out[i] = np.asarray(gen.GetFingerprint(m), dtype=np.float32)
    return out


@representation("rdkit_descriptors")
def rdkit_descriptors(smiles):
    """Physicochemical descriptors, z-scored. Tests whether bulk properties suffice."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")
    names = [n for n, _ in Descriptors.descList]
    calc = dict(Descriptors.descList)
    out = np.zeros((len(smiles), len(names)), dtype=np.float64)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        for j, n in enumerate(names):
            try:
                v = calc[n](m)
            except Exception:  # noqa: BLE001 - individual descriptors can fail
                v = 0.0
            out[i, j] = v if np.isfinite(v) else 0.0
    out = np.nan_to_num(out, posinf=0.0, neginf=0.0)
    sd = out.std(0)
    sd[sd == 0] = 1.0
    return ((out - out.mean(0)) / sd).astype(np.float32)


@representation("random")
def random_rep(smiles):
    """Random vectors - the floor. Any representation must beat this."""
    rng = np.random.default_rng(SEED)
    return rng.standard_normal((len(smiles), 64)).astype(np.float32)


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
    """Mean precision@k, recall@k, MRR and AUROC over queries with >=1 relevant hit."""
    n = sim.shape[0]
    sim = sim.copy()
    np.fill_diagonal(sim, -np.inf)
    n_rel = rel.sum(1)
    queries = np.flatnonzero(n_rel > 0)

    prec = {k: [] for k in KS}
    recall = {k: [] for k in KS}
    rr, auroc = [], []

    for i in queries:
        order = np.argsort(-sim[i], kind="stable")
        order = order[order != i]
        hits = rel[i][order]

        for k in KS:
            prec[k].append(hits[:k].mean())
            recall[k].append(hits[:k].sum() / n_rel[i])

        first = np.flatnonzero(hits)
        rr.append(1.0 / (first[0] + 1) if first.size else 0.0)

        # AUROC over the ranking: P(relevant ranked above irrelevant)
        ranks = np.empty(hits.size)
        ranks[order.argsort()] = np.arange(hits.size)
        pos = n_rel[i]
        neg = hits.size - pos
        if neg > 0:
            rank_of_pos = np.flatnonzero(hits) + 1
            auroc.append((rank_of_pos.sum() - pos * (pos + 1) / 2) / (pos * neg))
            auroc[-1] = 1.0 - auroc[-1]

    out = {"n_queries": int(queries.size)}
    for k in KS:
        out[f"p@{k}"] = float(np.mean(prec[k]))
        out[f"r@{k}"] = float(np.mean(recall[k]))
    out["mrr"] = float(np.mean(rr))
    out["auroc"] = float(np.mean(auroc))
    return out


def ceiling(rel):
    """Perfect ranking: all relevant drugs first. The maximum any method could score."""
    n = rel.shape[0]
    sim = rel.astype(np.float64) + np.random.default_rng(SEED).uniform(0, 1e-9, (n, n))
    return evaluate(sim, rel)


def main():
    argv = sys.argv[1:]
    truth = "targets"
    if "--truth" in argv:
        i = argv.index("--truth")
        truth = argv[i + 1]
        del argv[i : i + 2]
    if truth not in ("targets", "moa_targets"):
        sys.exit(f"--truth must be 'targets' or 'moa_targets', got {truth!r}")

    wanted = argv or list(REPRESENTATIONS)
    unknown = [w for w in wanted if w not in REPRESENTATIONS]
    if unknown:
        sys.exit(f"unknown representation(s): {unknown}\nknown: {list(REPRESENTATIONS)}")

    rows = load_drugs(truth)
    print(f"ground truth: shared {truth}")
    print(f"{len(rows)} drugs with a structure and >=1 target")
    rel = relevance_matrix(rows)
    n_rel = rel.sum(1)
    print(f"{int((n_rel > 0).sum())} have >=1 target-sharing neighbour "
          f"(median {int(np.median(n_rel[n_rel > 0]))} neighbours)")
    print(f"random-pair relevance rate: {rel.mean():.4f}\n")

    smiles = [r["smiles"] for r in rows]
    results = {"_ceiling": ceiling(rel), "_meta": {"truth": truth, "n_drugs": len(rows)}}

    for name in wanted:
        fn = REPRESENTATIONS[name]
        print(f"computing {name}...", flush=True)
        mat = fn(smiles)
        results[name] = evaluate(similarity(mat, getattr(fn, "binary", False)), rel)
        results[name]["dim"] = int(mat.shape[1])

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"benchmark_{truth}.json"
    out.write_text(json.dumps(results, indent=2))

    cols = ["p@1", "p@5", "p@10", "r@10", "mrr", "auroc"]
    order = [n for n in wanted if n != "random"]
    order = ["random"] + order if "random" in wanted else order
    width = max(len(n) for n in results) + 2
    print(f"\n{'method'.ljust(width)}" + "".join(c.rjust(9) for c in cols))
    print("-" * (width + 9 * len(cols)))
    for name in order + ["_ceiling"]:
        r = results[name]
        print(name.ljust(width) + "".join(f"{r[c]:9.4f}" for c in cols))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
