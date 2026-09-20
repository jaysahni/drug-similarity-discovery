"""Arm P: design a peptide against the target, then find the approved peptide it resembles.

This is the direct-matching pipeline, peptide flavour:

    thrombin active site
      -> BoltzGen peptide-anything, 8..16 residues        [Rowan]
      -> ESM-C embedding of each design                   [local, MPS]
      -> nearest approved peptide by cosine               [local]
      -> null, novelty, three leaderboards                [demo/match_direct.py]

Unlike the footprint arm, the generated molecule is the *query* here rather than
a probe that gets discarded, so the comparison is peptide-to-peptide and a
protein language model is the right tool.

Two metrics are reported side by side, deliberately:

  esmc_cosine     mean-pooled ESM-C, the conventional PLM similarity
  identity        BLOSUM62 local alignment, normalised by the shorter sequence

ESM-C is trained on domains, not octapeptides, and its behaviour on very short
sequences is not something this repo has measured. Reporting one number would
hide that. Where the two disagree, the disagreement is the finding.

The corpus is approved biologics of at most MAX_LEN residues. Comparing a 12-mer
to a 450-residue antibody heavy chain would measure length and composition, not
binding, so the long entries are excluded rather than allowed to dominate.

Run:  DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.peptide_arm
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np

from demo import esmc, match_direct as md

REPO = Path(__file__).resolve().parent.parent
BIOLOGICS = REPO / "data" / "approved_biologics.csv"

MAX_LEN = int(os.environ.get("DEMO_PEPTIDE_MAX_LEN", "100"))
N_NULL = int(os.environ.get("DEMO_N_NULL", "300"))
SEED = 0


def load_corpus(symbol: str) -> list[dict]:
    """Approved biologics short enough to compare against a designed peptide."""
    rows = []
    for row in csv.DictReader(BIOLOGICS.open()):
        length = row.get("length_aa") or ""
        if not length.isdigit() or int(length) > MAX_LEN:
            continue
        sequence = (row.get("sequence") or "").strip().upper()
        if not sequence:
            continue
        targets = row.get("target_gene_symbol") or ""
        rows.append(
            {
                "name": row["name"],
                "modality": row.get("modality"),
                "length_aa": int(length),
                "sequence": sequence,
                "targets": targets,
                # The corpus has no MOA column, so a hit on this target counts as
                # known_moa only when the target is its sole annotation; anything
                # else that lists the target is known_offtarget. Conservative in
                # the direction that matters -- it cannot invent novelty.
                "novelty": md.classify_novelty(
                    targets, targets if targets == symbol else "", symbol
                ),
            }
        )
    if not rows:
        raise SystemExit(f"no biologics at or under {MAX_LEN} aa in {BIOLOGICS.name}")
    return rows


def identity(a: str, b: str) -> float:
    """BLOSUM62 local-alignment identity, normalised by the shorter sequence.

    A designed 12-mer against a 65-mer cannot exceed 12 identities, so
    normalising by the shorter sequence is what keeps short queries comparable.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    try:
        alignment = aligner.align(a, b)[0]
    except (IndexError, ValueError):
        return 0.0
    matches = sum(
        1
        for x, y in zip(alignment[0], alignment[1])
        if x == y and x != "-"
    )
    return matches / min(len(a), len(b))


def shuffled_decoys(sequences: list[str], n: int, seed: int = SEED) -> list[str]:
    """Composition-preserving shuffles: the null asks "is the ORDER doing work?"

    A decoy built by shuffling a design keeps its length and amino-acid
    composition exactly, so anything the real design scores above this null is
    attributable to sequence order rather than to being, say, unusually
    arginine-rich.
    """
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        base = list(sequences[len(out) % len(sequences)])
        rng.shuffle(base)
        out.append("".join(base))
    return out


def main() -> None:
    target_name = os.environ.get("DEMO_TARGET", "thrombin")
    from demo import pipeline  # imported here so DEMO_TARGET is already set

    symbol = pipeline.TARGET["symbol"]
    out_dir = pipeline.OUT
    designs = json.loads((out_dir / "03_designs.json").read_text())

    sequences = [d["sequence"] for d in designs["designs"] if d.get("sequence")]
    if not sequences:
        raise SystemExit("03_designs.json holds no binder sequences")
    print(f"[arm P] {len(sequences)} designs, {designs.get('binder_length')} residues, "
          f"protocol {designs.get('protocol')}")

    corpus = load_corpus(symbol)
    print(f"  corpus: {len(corpus)} approved biologics <= {MAX_LEN} aa; "
          f'{sum(1 for r in corpus if r["novelty"] != "novel_pairing")} annotated to {symbol}')

    corpus_vecs = esmc.embed([r["sequence"] for r in corpus])
    design_vecs = esmc.embed(sequences)

    # The null: shuffled versions of the designs themselves.
    decoys = shuffled_decoys(sequences, N_NULL)
    null = md.null_distribution(esmc.embed(decoys), corpus_vecs, "cosine")
    print(f"  null: n={null.size} mean={null.mean():.4f} p95={np.percentile(null, 95):.4f}")

    results = []
    for design, vec in zip(designs["designs"], design_vecs):
        sims = md.cosine(vec, corpus_vecs)
        rows = []
        for row, sim in zip(corpus, sims):
            rows.append(
                dict(
                    row,
                    esmc_cosine=round(float(sim), 4),
                    identity=round(identity(design["sequence"], row["sequence"]), 4),
                )
            )
        nn = max(rows, key=lambda r: r["esmc_cosine"])
        by_identity = max(rows, key=lambda r: r["identity"])
        boards = md.leaderboards(rows, score_key="esmc_cosine")
        results.append(
            {
                "design_id": design["design_id"],
                "sequence": design["sequence"],
                "length": len(design["sequence"]),
                "design_to_target_iptm": design.get("design_to_target_iptm"),
                "nearest_by_esmc": {k: nn[k] for k in ("name", "esmc_cosine", "identity", "novelty", "targets")},
                "nearest_by_identity": {k: by_identity[k] for k in ("name", "esmc_cosine", "identity", "novelty", "targets")},
                "metrics_agree": nn["name"] == by_identity["name"],
                "nn_cosine": round(float(sims.max()), 4),
                "null_percentile": round(md.percentile_of(float(sims.max()), null), 4),
                "enrichment": md.enrichment(rows, score_key="esmc_cosine", top_k=5, symbol=symbol),
                "leaderboards": {
                    name: [
                        {k: r[k] for k in ("name", "esmc_cosine", "identity", "novelty", "rank_in_board")}
                        for r in board[:10]
                    ]
                    for name, board in boards.items()
                },
            }
        )

    results.sort(key=lambda r: -r["null_percentile"])
    payload = {
        "arm": "P (peptide, direct matching)",
        "target": {"symbol": symbol, "uniprot": pipeline.TARGET["uniprot"]},
        "question": "which approved peptide does the designed binder most resemble?",
        "corpus": {
            "file": str(BIOLOGICS.relative_to(REPO)),
            "n": len(corpus),
            "max_length_aa": MAX_LEN,
            "n_annotated_to_target": sum(1 for r in corpus if r["novelty"] != "novel_pairing"),
            "annotated_names": [r["name"] for r in corpus if r["novelty"] != "novel_pairing"],
            "base_rate": round(sum(1 for r in corpus if r["novelty"] != "novel_pairing") / len(corpus), 4),
        },
        "null": {
            "construction": "composition-preserving shuffles of the designs themselves",
            **md.calibration(null),
        },
        "embedding": {"model": esmc.MODEL, "dim": int(corpus_vecs.shape[1])},
        "designs": results,
        "caveats": [
            "ESM-C is trained on domains, not on 8-16mers; its behaviour at this "
            "length is not characterised here. The identity column is reported "
            "beside it for that reason, and disagreement between them is a result.",
            "The corpus is 49-ish approved biologics under the length cap, so a "
            "nearest neighbour exists by construction. Only the null percentile "
            "says whether it is a meaningful one.",
            "Similarity is not affinity. Nothing here predicts binding.",
        ],
    }

    path = out_dir / "06_peptide_match.json"
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"  -> {path.relative_to(REPO)}")

    print(f'\n{"design":>6s} {"len":>4s} {"nearest (ESM-C)":24s} {"cos":>7s} {"ident":>6s} {"pctile":>7s} {"agree":>6s}')
    for r in results:
        n = r["nearest_by_esmc"]
        print(f'{r["design_id"]:6d} {r["length"]:4d} {n["name"][:24]:24s} '
              f'{n["esmc_cosine"]:7.4f} {n["identity"]:6.3f} {r["null_percentile"]:7.3f} '
              f'{"yes" if r["metrics_agree"] else "NO":>6s}')


if __name__ == "__main__":
    main()
