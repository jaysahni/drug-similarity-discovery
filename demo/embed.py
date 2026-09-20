"""Interface embeddings: one vector space for the design and for every drug.

The ask this answers: embed the top BoltzGen design, embed the approved drugs,
rank the drugs by similarity to the design.

The obvious way to do that -- ESM the design, fingerprint the drugs, cosine the
two -- does not work, and not for want of a better model. ESM-2 returns 1280
dimensions of learned protein-sequence features; ECFP4 returns 2048 dimensions of
hashed chemical substructures. The axes are unrelated, so the cosine is arithmetic
over quantities that are not commensurable. PROJECT_GOAL.md 1.1 states the same
thing: similarity between a de novo miniprotein and a small molecule is undefined.

So both sides are embedded on axes that DO mean the same thing: the target's
residues. A vector is one coordinate per residue of the target, holding how
strongly that residue is engaged.

  design vector   frequency, across the surviving designs, of engaging residue i
                  (or a binary contact vector, for a single design)
  drug vector     1 if the co-folded pose puts a ligand heavy atom within 4.5 A
                  of residue i, else 0

Cosine between those is well-posed: same axes, same units, same length.

Reads what the pipeline already computed; costs no credits and needs no GPU.
Run:  ./env-kit/bin/python -m demo.embed
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "demo" / "cdk2"


def load(name: str) -> dict:
    path = OUT / name
    if not path.exists():
        raise SystemExit(f"{path.relative_to(REPO)} missing -- run the pipeline first.")
    return json.loads(path.read_text())


def vector(engaged: dict[int, float], axis: list[int]) -> np.ndarray:
    """Project an engagement map onto the shared residue axis."""
    return np.array([float(engaged.get(r, 0.0)) for r in axis])


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator == 0.0:
        return 0.0
    return float(np.dot(a, b) / denominator)


def main() -> None:
    research = load("01_research.json")
    designs = load("03_designs.json")
    signature = load("04_signature.json")
    ranked = load("05_ranked.json")

    # The shared axis: every residue of the target. Fixed and explicit, so a
    # residue missing from one side is a zero rather than a silent misalignment.
    n_residues = research["uniprot"]["sequence_length"]
    axis = list(range(1, n_residues + 1))

    # -- design side ------------------------------------------------------
    # "Top design" by the quality score Rowan reports. Also embed the
    # consensus, which is the more defensible object: one design is one sample.
    scored = [d for d in designs["designs"] if d.get("quality_score") is not None]
    top = max(scored, key=lambda d: d["quality_score"]) if scored else None

    per_design = {d["design_id"]: d["residues_author"] for d in signature["per_design"]}
    top_contacts = per_design.get(top["design_id"], []) if top else []

    design_vectors = {
        "top_design": vector({r: 1.0 for r in top_contacts}, axis),
        "consensus": vector(
            {int(r): f for r, f in signature["frequency"].items()}, axis
        ),
    }

    # -- drug side --------------------------------------------------------
    rows = []
    for drug in ranked["results"]:
        drug_vector = vector({r: 1.0 for r in drug["engaged_residues"]}, axis)
        row = {
            "name": drug["name"],
            "tier": drug["tier"],
            "hits_target": drug["hits_target"],
            "iptm": drug["iptm"],
            "n_engaged_residues": drug["n_engaged_residues"],
            "core_coverage_pocket": drug["core_coverage_pocket"],
        }
        for label, design_vector in design_vectors.items():
            row[f"cosine_{label}"] = round(cosine(design_vector, drug_vector), 4)
        rows.append(row)

    for label in design_vectors:
        field = f"cosine_{label}"
        for i, row in enumerate(sorted(rows, key=lambda r: -r[field]), start=1):
            row[f"rank_{label}"] = i

    rows.sort(key=lambda r: -r["cosine_consensus"])

    tiers = {}
    for tier in ("positive", "candidate", "decoy"):
        members = [r for r in rows if r["tier"] == tier]
        if members:
            tiers[tier] = {
                "n": len(members),
                "mean_cosine_consensus": round(
                    sum(r["cosine_consensus"] for r in members) / len(members), 4
                ),
                "mean_cosine_top_design": round(
                    sum(r["cosine_top_design"] for r in members) / len(members), 4
                ),
            }

    payload = {
        "question": "which approved drug's binding footprint most resembles the BoltzGen design's?",
        "space": {
            "axes": f"the {n_residues} residues of {research['target']['symbol']} "
            f"({research['target']['uniprot']}), author numbering",
            "design_side": "engagement frequency across surviving designs (consensus), "
            "or binary contact (single design)",
            "drug_side": "binary contact, ligand heavy atom within 4.5 A, from the co-folded pose",
            "metric": "cosine",
            "why_not_esm": (
                "ESM-2 embeds sequences and ECFP4 embeds substructures; their axes "
                "are unrelated, so a cosine between them is not defined. Not a model "
                "accuracy issue -- there is no shared coordinate system. Projecting "
                "both onto target residues gives one."
            ),
        },
        "top_design": {
            "design_id": top["design_id"] if top else None,
            "quality_score": top["quality_score"] if top else None,
            "iptm": top["iptm"] if top else None,
            "sequence": top["sequence"] if top else None,
            "n_contacts": len(top_contacts),
        },
        "consensus_core": signature["core"],
        "by_tier": tiers,
        "ranked": rows,
        "caveats": [
            "A single design is one sample from a generator whose designs agree with "
            f'each other at mean pairwise Jaccard {signature["mean_pairwise_jaccard"]}. '
            "The consensus column is the more defensible of the two.",
            "Cosine here rewards engaging the same residues. It does not distinguish a "
            "hydrogen bond from a hydrophobic brush, and it does not weight by buried "
            "area, so two molecules touching the same residues differently score alike.",
            "Every drug vector comes from a Boltz-2 co-fold, so a drug already known to "
            "bind this target may be reproducing a complex from the model's training set.",
        ],
    }

    path = OUT / "06_embedding_similarity.json"
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"-> {path.relative_to(REPO)}")

    print(f'\ntop design: #{payload["top_design"]["design_id"]} '
          f'(quality {payload["top_design"]["quality_score"]}, '
          f'{payload["top_design"]["n_contacts"]} contacts)')
    print(f'\n{"drug":16s} {"tier":10s} {"cos(consensus)":>15s} {"cos(top)":>9s} {"coverage":>9s}')
    for row in rows:
        print(f'{row["name"]:16s} {row["tier"]:10s} '
              f'{row["cosine_consensus"]:15.4f} {row["cosine_top_design"]:9.4f} '
              f'{row["core_coverage_pocket"]:9.2f}')
    print("\nby tier:", json.dumps(tiers))


if __name__ == "__main__":
    main()
