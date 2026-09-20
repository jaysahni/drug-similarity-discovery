"""Do the two arms agree, and does combining them beat either alone?

This project ranks approved drugs for a target two independent ways and has
never compared them:

  FOOTPRINT   co-fold the drug, score how much of the designed binder's hotspot
              set its pose engages            (results/demo/<target>/05_ranked.json)
  CHEMICAL    ECFP4 Tanimoto to the drugs already known to hit the target
              (scripts/representations.py, the benchmark's best of 16)

They fail differently, which is the reason to look. Chemical similarity ranks by
chemotype -- argatroban retrieves peptidomimetics, not thrombin drugs. Footprint
coverage ignores chemistry entirely and ranked a statin first on CDK2. Where two
rankings with unrelated failure modes agree, the agreement is worth more than
either ranking alone; where they disagree, that is worth knowing before either is
trusted.

Combination is by Reciprocal Rank Fusion (demo/fusion.py), for the same reason it
is used there: hotspot coverage is a fraction of a residue set and Tanimoto is a
fraction of a bit set, and putting them on a common scale would invent a
calibration nobody measured. RRF uses order only.

**n is small.** The only target with both arms computed is CDK2, with 14 drugs and
5 positives. Nothing here is powered; it is reported with its n and read as a
direction, not a result.

Run:  DEMO_TARGET=cdk2 ./env-kit/bin/python -m demo.combine_arms
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, mannwhitneyu

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from demo.fusion import reciprocal_rank_fusion  # noqa: E402


def chemical_scores(shortlist: list[dict], symbol: str) -> tuple[np.ndarray, list[str]]:
    """Each shortlisted drug's best ECFP4 Tanimoto to a KNOWN binder of the target.

    The known binders are drawn from the approved library, excluding the drug
    being scored -- otherwise a positive scores 1.0 against itself and the arm is
    a lookup rather than a prediction.
    """
    import csv

    import representations
    from demo.match_direct import tanimoto

    approved = [
        r for r in csv.DictReader((REPO / "data" / "approved_drugs.csv").open())
        if r["smiles"] and r.get("approved") == "1"
    ]
    known = [r for r in approved if symbol in r["moa_targets"].split(";")]
    if len(known) < 2:
        known = [r for r in approved if symbol in r["targets"].split(";")]
    if not known:
        raise SystemExit(f"no approved drug annotated to {symbol}; cannot score chemically")

    names = [d["name"] for d in shortlist]
    mat = representations.REPRESENTATIONS["morgan"](
        [d["smiles"] for d in shortlist] + [r["smiles"] for r in known]
    )
    query_fp, known_fp = mat[: len(shortlist)], mat[len(shortlist):]
    known_names = [r["name"].lower() for r in known]

    out = np.zeros(len(shortlist))
    for i, drug in enumerate(shortlist):
        sims = tanimoto(query_fp[i], known_fp)
        for j, kn in enumerate(known_names):
            if kn == drug["name"].lower():
                sims[j] = -np.inf          # never score a drug against itself
        out[i] = float(sims.max()) if np.isfinite(sims).any() else 0.0
    return out, [r["name"] for r in known]


def main() -> None:
    from demo import pipeline

    symbol = pipeline.TARGET["symbol"]
    out_dir = pipeline.OUT
    ranked = json.loads((out_dir / "05_ranked.json").read_text())
    drugs = ranked["results"]
    print(f"[arms] {symbol}: {len(drugs)} drugs from 05_ranked.json")

    footprint = np.array([d["core_coverage_pocket"] for d in drugs])
    chemical, known_names = chemical_scores(drugs, symbol)
    print(f"  chemical reference set: {len(known_names)} approved {symbol} binders")

    # RRF needs matrices; treat each score vector as a single query's ranking.
    fused = reciprocal_rank_fusion([footprint[None, :], chemical[None, :]])[0]

    tau, tau_p = kendalltau(footprint, chemical)
    is_pos = np.array([d["tier"] == "positive" for d in drugs])
    is_dec = np.array([d["tier"] == "decoy" for d in drugs])

    arms = {"footprint": footprint, "chemical": chemical, "fused": fused}
    summary = {}
    for name, score in arms.items():
        order = np.argsort(-score)
        ranks = {drugs[i]["name"]: int(np.where(order == i)[0][0]) + 1 for i in range(len(drugs))}
        u, p = mannwhitneyu(score[is_pos], score[is_dec], alternative="greater")
        summary[name] = {
            "mean_positive": round(float(score[is_pos].mean()), 4),
            "mean_decoy": round(float(score[is_dec].mean()), 4),
            "mannwhitney_U": float(u),
            "p_positive_gt_decoy": float(f"{p:.4g}"),
            "top_ranked": drugs[int(order[0])]["name"],
            "top_ranked_tier": drugs[int(order[0])]["tier"],
            "positive_ranks": sorted(ranks[d["name"]] for d in drugs if d["tier"] == "positive"),
        }

    payload = {
        "target": symbol,
        "n_drugs": len(drugs),
        "n_positive": int(is_pos.sum()),
        "n_decoy": int(is_dec.sum()),
        "chemical_reference_set": known_names,
        "agreement": {
            "kendall_tau": round(float(tau), 4),
            "p_value": float(f"{tau_p:.4g}"),
            "reading": "how much the two arms order the shortlist alike; "
                       "0 would mean they carry independent information",
        },
        "arms": summary,
        "caveats": [
            f"n = {len(drugs)} drugs, {int(is_pos.sum())} positives. Nothing here is "
            "powered; read the direction, not the p-value.",
            "The chemical arm scores each drug against OTHER known binders, never "
            "itself, but those binders are curated annotations -- so this arm is "
            "closer to a lookup than the footprint arm is, and should be expected "
            "to do better for that reason alone.",
            "Only CDK2 has both arms computed. No cross-target claim.",
        ],
    }
    path = out_dir / "10_arm_agreement.json"
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"  -> {path.relative_to(REPO)}")

    print(f'\n  Kendall tau between arms: {tau:.4f} (p = {tau_p:.3g})')
    print(f'\n{"arm":12s} {"pos":>7s} {"decoy":>7s} {"p":>9s}  {"top-ranked":22s} positive ranks')
    for name, s in summary.items():
        print(f'{name:12s} {s["mean_positive"]:7.3f} {s["mean_decoy"]:7.3f} '
              f'{s["p_positive_gt_decoy"]:9.4f}  {s["top_ranked"][:20]:22s} {s["positive_ranks"]}')


if __name__ == "__main__":
    main()
