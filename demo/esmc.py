"""ESM-C protein embeddings, for the biologic half of the repurposing library.

The small-molecule side of this demo is already covered: ECFP4 over SMILES, via
rdkit. Approved *biologics* -- antibodies, enzymes, hormones -- have no SMILES,
so they cannot enter that library at all. This module gives them a vector.

Why ESM-C rather than the ESM-2 already in the toolkit: ESM-C 300M reaches
ESM-2 650M's quality at half the parameters, it is a 2024/2026 model rather than
a 2022 one, and -- the part that decides it for a public repo -- the weights are
MIT (Chan Zuckerberg Biohub, 2026). See docs/09-ESMC.md for the licence check and
for the measured load and per-sequence times that chose the 300M size.

Two things here are deliberately loud rather than convenient:

  * Sequences longer than the trained 2048-token context are **refused**, not
    truncated. MEASURED: ESM-C does not raise on a 4000-residue input -- its
    rotary embeddings extrapolate and it returns finite numbers that look
    perfectly ordinary. A silent extrapolation is worse than a crash, so the
    crash is ours.
  * Anything outside ESM-C's alphabet is refused. The tokenizer maps lowercase
    'a' to <unk> without complaint, which would quietly embed a lowercased FASTA
    as a string of unknown residues.

Cosine on the mean-pooled vector is the conventional similarity for protein
language model embeddings; `similarity` does that and nothing cleverer.

Run:  ./env-kit/bin/python -m demo.esmc
"""

from __future__ import annotations

import hashlib
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "results" / "demo" / "_cache" / "esmc"

MODEL = "biohub/ESMC-300M"

# ESM-C was trained at a 2048-token context (512 for the first 1M steps, 2048
# for the last 500k -- model card). The tokenizer wraps a sequence in <cls> and
# <eos>, so 2046 residues is the longest input that stays inside it.
CONTEXT_TOKENS = 2048
MAX_RESIDUES = CONTEXT_TOKENS - 2

# The 20 standard residues plus the ambiguity codes ESM-C tokenises rather than
# maps to <unk>. '|' is the model's own chain-break token and is inserted by
# this module, never accepted from a caller inside a chain.
ALPHABET = set("LAGVSERTIDPKQNFYMHWCXBUZO")

_LOADED: dict[str, tuple] = {}


def device() -> str:
    """Best available torch device.

    MEASURED on this machine (docs/09-ESMC.md): MPS is 2.1x faster than CPU for
    a batch of 8 and agrees with it to 1.3e-5 on the hidden states, so MPS is
    preferred when present rather than assumed to help.
    """
    import torch

    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def load(model: str = MODEL, dev: str | None = None) -> tuple:
    """Return (tokenizer, model, device), loading once per process."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    dev = dev or device()
    key = f"{model}@{dev}"
    if key not in _LOADED:
        tokenizer = AutoTokenizer.from_pretrained(model)
        net = AutoModel.from_pretrained(model, dtype=torch.float32).to(dev).eval()
        _LOADED[key] = (tokenizer, net, dev)
    return _LOADED[key]


def check(sequence: str) -> str:
    """Validate one chain. Raises rather than embedding something meaningless."""
    if not isinstance(sequence, str):
        raise TypeError(f"expected a sequence string, got {type(sequence).__name__}")
    seq = "".join(sequence.split())
    if not seq:
        raise ValueError("empty sequence: refusing to return a vector for nothing")
    bad = sorted(set(seq) - ALPHABET)
    if bad:
        raise ValueError(
            f"characters outside ESM-C's alphabet: {bad}. The tokenizer maps these "
            "to <unk> silently -- lowercase FASTA is the usual cause. "
            f"Sequence starts {seq[:30]!r}"
        )
    if len(seq) > MAX_RESIDUES:
        raise ValueError(
            f"sequence is {len(seq)} residues, over ESM-C's {MAX_RESIDUES}-residue "
            f"trained context ({CONTEXT_TOKENS} tokens including <cls>/<eos>). "
            "ESM-C does NOT raise on this -- it extrapolates its rotary embeddings "
            "and returns finite, ordinary-looking numbers. Split the chain or embed "
            "a domain; do not truncate silently."
        )
    return seq


def _key(model: str, sequence: str) -> str:
    return hashlib.sha256(f"{model}\n{sequence}".encode()).hexdigest()[:24]


def embed(
    sequences: list[str],
    *,
    model: str = MODEL,
    dev: str | None = None,
    batch_size: int = 8,
) -> np.ndarray:
    """Mean-pooled ESM-C embeddings, shape (n, d). d = 960 for the 300M model.

    Pooling is the mean over *residue* positions only: <cls>, <eos> and padding
    are masked out. MEASURED: padded batches agree with one-at-a-time inference
    to 0.0 exactly, so batching costs nothing in fidelity.

    Cached per sequence under results/demo/_cache/esmc/, keyed by a hash of the
    model name and the sequence, so a re-run does no inference at all.
    """
    import torch

    if not isinstance(sequences, (list, tuple)):
        raise TypeError("embed() takes a list of sequences, not a bare string")
    if not sequences:
        raise ValueError("embed() got no sequences")

    clean = [check(s) for s in sequences]
    CACHE.mkdir(parents=True, exist_ok=True)

    out: list[np.ndarray | None] = [None] * len(clean)
    todo: list[int] = []
    for i, seq in enumerate(clean):
        path = CACHE / f"{_key(model, seq)}.npy"
        if path.exists():
            out[i] = np.load(path)
        else:
            todo.append(i)

    if todo:
        tokenizer, net, dev = load(model, dev)
        # Group by length so a short chain is not padded out to a long one.
        todo.sort(key=lambda i: len(clean[i]))
        for start in range(0, len(todo), batch_size):
            chunk = todo[start : start + batch_size]
            enc = tokenizer([clean[i] for i in chunk], return_tensors="pt", padding=True)
            enc = {k: v.to(dev) for k, v in enc.items()}
            with torch.inference_mode():
                hidden = net(**enc).last_hidden_state

            mask = enc["attention_mask"].clone()
            mask[:, 0] = 0  # <cls>
            lengths = enc["attention_mask"].sum(dim=1)
            mask[torch.arange(mask.shape[0], device=dev), lengths - 1] = 0  # <eos>
            weights = mask.unsqueeze(-1).to(hidden.dtype)
            pooled = (hidden * weights).sum(dim=1) / weights.sum(dim=1)

            for row, i in enumerate(chunk):
                vec = pooled[row].float().cpu().numpy()
                if not np.isfinite(vec).all():
                    raise RuntimeError(
                        f"non-finite embedding for a {len(clean[i])}-residue sequence"
                    )
                out[i] = vec
                np.save(CACHE / f"{_key(model, clean[i])}.npy", vec)

    return np.vstack(out)


# How to turn a multi-chain entity (an antibody's heavy + light, a haemoglobin
# tetramer's alpha + beta) into one vector. There is no settled answer, so the
# choice is a parameter and docs/09-ESMC.md carries the measured comparison.
#
#   "mean"        default. Embed each chain, average the vectors, each chain
#                 weighted equally. Keeps the result in R^d, so a two-chain
#                 antibody and a single-chain enzyme stay comparable -- which is
#                 the whole point of a library that mixes both. Equal weighting
#                 rather than length weighting so a 450-aa heavy chain does not
#                 swamp a 215-aa light chain that carries three of the six CDRs.
#   "length"      length-weighted mean; identical to pooling the concatenation's
#                 residues. Use when the entity really is one long molecule.
#   "chainbreak"  join the chains with '|', ESM-C's own chain-break token, and
#                 embed once. Model-native, and the only option where the model
#                 sees both chains at the same time -- but the total must fit in
#                 the 2046-residue context.
#   "concat"      stack the chain vectors end to end. Preserves everything, but
#                 the result is in R^(k*d): comparable only to entities with the
#                 same number of chains in the same order. Not the default for
#                 that reason.
COMBINE = ("mean", "length", "chainbreak", "concat")


def embed_entity(
    chains: list[str],
    *,
    combine: str = "mean",
    model: str = MODEL,
    dev: str | None = None,
) -> np.ndarray:
    """One vector for a multi-chain entity. See COMBINE for the four modes."""
    if combine not in COMBINE:
        raise ValueError(f"combine must be one of {COMBINE}, got {combine!r}")
    if isinstance(chains, str):
        chains = [chains]
    if not chains:
        raise ValueError("embed_entity() got no chains")

    clean = [check(c) for c in chains]

    if combine == "chainbreak":
        total = sum(len(c) for c in clean) + len(clean) - 1
        if total > MAX_RESIDUES:
            raise ValueError(
                f"{len(clean)} chains plus chain breaks come to {total} tokens, over "
                f"ESM-C's {MAX_RESIDUES}. Use combine='mean', which embeds each chain "
                "separately and never exceeds the context."
            )
        # check() rejects '|' inside a chain, so this join is unambiguous.
        joined = "|".join(clean)
        tokenizer, net, dev = load(model, dev)
        import torch

        enc = tokenizer(joined, return_tensors="pt")
        enc = {k: v.to(dev) for k, v in enc.items()}
        with torch.inference_mode():
            hidden = net(**enc).last_hidden_state
        # Pool over residues only: drop <cls>, <eos> and the chain-break tokens.
        ids = enc["input_ids"][0]
        keep = torch.ones_like(ids, dtype=torch.bool)
        keep[0] = False
        keep[-1] = False
        keep &= ids != tokenizer.convert_tokens_to_ids("|")
        return hidden[0][keep].mean(dim=0).float().cpu().numpy()

    vectors = embed(clean, model=model, dev=dev)
    if combine == "mean":
        return vectors.mean(axis=0)
    if combine == "length":
        w = np.array([len(c) for c in clean], dtype=float)
        return (vectors * w[:, None]).sum(axis=0) / w.sum()
    return vectors.reshape(-1)  # concat


def embed_entities(
    entities: list[list[str]], *, combine: str = "mean", model: str = MODEL
) -> np.ndarray:
    """embed_entity over a list. Shape (n, d) except for combine='concat'."""
    if not entities:
        raise ValueError("embed_entities() got no entities")
    rows = [embed_entity(e, combine=combine, model=model) for e in entities]
    widths = {r.shape[0] for r in rows}
    if len(widths) > 1:
        raise ValueError(
            f"entities embedded to different widths {sorted(widths)} -- combine="
            f"{combine!r} only stacks if every entity has the same chain count. "
            "Use combine='mean'."
        )
    return np.vstack(rows)


def similarity(a: np.ndarray, b: np.ndarray) -> float | np.ndarray:
    """Cosine similarity between ESM-C embeddings.

    Cosine is the conventional choice for protein language model embeddings:
    mean-pooled PLM vectors vary in norm with sequence length and composition,
    and cosine discards exactly that, comparing direction only. Euclidean
    distance on the same vectors ranks a long protein as far from a short one
    largely because it is long.

    1-D against 1-D returns a float; 2-D against 2-D returns the (n, m) matrix.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scalar = a.ndim == 1 and b.ndim == 1
    A = np.atleast_2d(a)
    B = np.atleast_2d(b)
    if A.shape[1] != B.shape[1]:
        raise ValueError(
            f"dimension mismatch: {A.shape[1]} vs {B.shape[1]}. Embeddings from "
            "different ESM-C sizes (300M=960, 600M=1152, 6B=2560) or different "
            "combine= modes do not share a space."
        )
    na = np.linalg.norm(A, axis=1, keepdims=True)
    nb = np.linalg.norm(B, axis=1, keepdims=True)
    if (na == 0).any() or (nb == 0).any():
        raise ValueError("zero-norm embedding: cosine is undefined")
    out = (A / na) @ (B / nb).T
    return float(out[0, 0]) if scalar else out


def main() -> None:
    """Load the model, embed one protein, and report what it measured."""
    seq = (
        "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVI"
        "HTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGA"
        "IKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEID"
        "QLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAA"
        "LAHPFFQDVTKPVPHLRL"
    )  # CDK2_HUMAN, UniProt P24941

    t0 = time.perf_counter()
    tokenizer, net, dev = load()
    load_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    v = embed([seq])
    embed_s = time.perf_counter() - t0

    print(f"model        {MODEL}")
    print(f"device       {dev}")
    print(f"parameters   {sum(p.numel() for p in net.parameters()) / 1e6:.0f}M")
    print(f"load         {load_s:.2f}s (weights already on disk)")
    print(f"embed        {embed_s:.3f}s for 1 sequence of {len(seq)} residues")
    print(f"shape        {v.shape}  norm {np.linalg.norm(v):.3f}")
    print(f"self-cosine  {similarity(v[0], v[0]):.6f}")
    print(f"cache        {CACHE.relative_to(REPO)}")


if __name__ == "__main__":
    main()
