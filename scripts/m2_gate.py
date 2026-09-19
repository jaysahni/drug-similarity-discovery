"""Ablation I6.1 - the M2 gate. Does a BoltzGen consensus beat pocket geometry?

PROJECT_GOAL.md section 8.3 calls this "the decisive ablation" and section 7 makes it
the gate the whole project turns on: if a generated-design consensus gives you nothing
that P2Rank gives you in half a second, the honest pipeline is the cheap one.

    arms, all scored against the SAME ground truth - the residues the known ligand of
    each held-out KDR co-crystal actually contacts (4.5 A heavy atom):

      boltzgen_consensus  PROJECT_GOAL.md 4.3 source="boltzgen_consensus".
                          Built from 24 BoltzGen protein-anything designs against
                          3VHE chain A. ONE set, used unchanged for every structure -
                          it never sees any held-out ligand.
      p2rank_geometry     source="p2rank_geometry". The rank-1 pocket P2Rank predicts
                          for that structure, run blind on the ligand-stripped form.
      known_ligand        source="known_ligand". Leave-one-out consensus over the
                          OTHER 37 structures' observed contacts. The ceiling
                          (ablation I6.2), for reference.
      random              |boltzgen core| residues drawn from the observed chain. Floor.

Also runs task C5, the convergence study: sub-sample the 24 designs at
N = 2,4,6,8,12,16,20,24 and measure how fast the consensus stabilises. The plan calls
this "the main cost lever" and expects convergence far below the design counts used
for wet-lab candidate selection.

Every arm is scored on the same 38 structures, so comparisons are PAIRED.

Usage:
    ./env/bin/python scripts/m2_gate.py
    ./env/bin/python scripts/m2_gate.py --core-threshold 0.5
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "results" / "pipeline" / "colorectal-cancer"
P2RANK_DIR = ROOT / "data" / "raw" / "p2rank"
SEED = 0


def p2rank_top1(pdb_id, chain):
    """Rank-1 pocket residues (author numbering) for one structure, or None."""
    d = P2RANK_DIR / f"{pdb_id}_{chain}"
    hits = list(d.glob("*_predictions.csv"))
    if not hits:
        return None
    with hits[0].open(newline="") as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            if row.get("rank") != "1":
                continue
            out = set()
            for tok in (row.get("residue_ids") or "").split():
                _, _, num = tok.rpartition("_")
                try:
                    out.add(int(num))
                except ValueError:
                    pass
            return out
    return None


def consensus(sets, threshold):
    """Residues engaged by >= threshold of the input contact sets."""
    if not sets:
        return set()
    c = Counter(r for s in sets for r in s)
    return {r for r, k in c.items() if k / len(sets) >= threshold}


def score(pred, truth):
    if not pred:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "jaccard": 0.0, "n_pred": 0}
    tp = len(pred & truth)
    p = tp / len(pred)
    r = tp / len(truth) if truth else float("nan")
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"precision": p, "recall": r, "f1": f1,
            "jaccard": tp / len(pred | truth), "n_pred": len(pred)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-threshold", type=float, default=0.6)
    ap.add_argument("--designs", default="v4")
    args = ap.parse_args()

    sig = json.loads((ROOT / "results" / f"boltzgen_signature_{args.designs}.json").read_text())
    boltz_core = set(sig["core_residue_ids_auth"])
    rec = json.loads((PIPE / "target" / "pocket_recovery.json").read_text())
    rmap = json.loads((PIPE / "target" / "boltzgen_residue_map.json").read_text())
    idx2auth = {v: int(k) for k, v in rmap["auth_to_index"].items()}

    print(f"boltzgen_consensus: {len(boltz_core)} core residues from "
          f"{sig['n_designs']} designs (threshold {sig['core_threshold']})")
    print(f"ground truth: {rec['n_structures']} KDR co-crystals, "
          f"known-ligand contacts at 4.5 A\n")

    # ground truth per structure, reconstructed from the recorded overlap + misses
    truth, meta = {}, {}
    for s in rec["per_structure"]:
        t1 = s["top1"]
        obs = {int(x) for x in t1["overlap_residue_labels"]} | \
              {int(x) for x in t1["missed_known_contact_labels"]}
        assert len(obs) == t1["n_known_contact_residues"], s["pdb_id"]
        truth[s["pdb_id"]] = obs
        meta[s["pdb_id"]] = {"chain": s["chain_id"], "ligand": s["ligand"],
                             "resolution": s["resolution_a"],
                             "n_chain": t1["n_observed_residues_in_chain"]}

    ids = sorted(truth)
    rng = np.random.default_rng(SEED)
    arms = {k: {m: [] for m in ("precision", "recall", "f1", "jaccard")}
            for k in ("boltzgen_consensus", "p2rank_geometry", "known_ligand", "random")}
    per_structure, skipped = [], []

    for pid in ids:
        O = truth[pid]
        pk = p2rank_top1(pid, meta[pid]["chain"])
        if pk is None:
            skipped.append(pid)
            continue
        loo = consensus([truth[o] for o in ids if o != pid], args.core_threshold)
        # floor: same size as the arm under test, drawn from the observed chain
        pool = np.arange(1, meta[pid]["n_chain"] + 1)
        rnd_scores = []
        for _ in range(100):
            picked = {idx2auth.get(int(i), -1) for i in
                      rng.choice(pool, size=min(len(boltz_core), pool.size), replace=False)}
            rnd_scores.append(score(picked, O))
        row = {"pdb_id": pid, **meta[pid], "n_known_contacts": len(O),
               "boltzgen_consensus": score(boltz_core, O),
               "p2rank_geometry": score(pk, O),
               "known_ligand": score(loo, O),
               "random": {m: float(np.mean([r[m] for r in rnd_scores]))
                          for m in ("precision", "recall", "f1", "jaccard")}}
        per_structure.append(row)
        for a in arms:
            for m in arms[a]:
                arms[a][m].append(row[a][m])

    n = len(per_structure)
    print(f"{n} structures scored" + (f", {len(skipped)} skipped ({skipped})" if skipped else ""))
    print(f"\n{'arm':<20}{'precision':>11}{'recall':>9}{'F1':>9}{'Jaccard':>9}{'|pred|':>8}")
    print("-" * 66)
    summary = {}
    for a in ("boltzgen_consensus", "p2rank_geometry", "known_ligand", "random"):
        v = {m: float(np.mean(arms[a][m])) for m in arms[a]}
        lo_hi = M.bootstrap_ci(np.array(arms[a]["jaccard"]), n_boot=10000, seed=SEED)
        v["jaccard_ci"] = [lo_hi[1], lo_hi[2]]
        npred = (len(boltz_core) if a == "boltzgen_consensus" else
                 float(np.mean([r[a]["n_pred"] for r in per_structure])) if a != "random"
                 else len(boltz_core))
        v["mean_n_pred"] = npred
        summary[a] = v
        print(f"{a:<20}{v['precision']:11.4f}{v['recall']:9.4f}{v['f1']:9.4f}"
              f"{v['jaccard']:9.4f}{npred:8.1f}")

    # --- the gate ---------------------------------------------------------
    comps, raw_p = {}, {}
    pairs = [("boltzgen_consensus", "p2rank_geometry"),
             ("boltzgen_consensus", "known_ligand"),
             ("boltzgen_consensus", "random"),
             ("known_ligand", "p2rank_geometry")]
    for a, b in pairs:
        for m in ("jaccard", "precision", "recall", "f1"):
            st = M.paired_bootstrap(np.array(arms[a][m]), np.array(arms[b][m]),
                                    n_boot=20000, seed=SEED)
            st["wilcoxon"] = M.wilcoxon(np.array(arms[a][m]), np.array(arms[b][m]))
            comps[f"{a} - {b} [{m}]"] = st
            raw_p[f"{a} - {b} [{m}]"] = st["p_value"]
    for k, adj in M.holm_bonferroni(raw_p).items():
        comps[k]["p_holm"] = adj

    print(f"\n=== THE GATE: boltzgen_consensus vs p2rank_geometry (n={n} paired) ===")
    for m in ("jaccard", "precision", "recall", "f1"):
        st = comps[f"boltzgen_consensus - p2rank_geometry [{m}]"]
        verdict = "BOLTZGEN" if st["delta"] > 0 and st["p_holm"] < 0.05 else (
            "P2RANK" if st["delta"] < 0 and st["p_holm"] < 0.05 else "no difference")
        print(f"  {m:<10} delta {st['delta']:+.4f} [{st['ci_lo']:+.4f}, {st['ci_hi']:+.4f}]"
              f"  Holm p={st['p_holm']:.3g}  -> {verdict}")

    # --- C5 convergence ---------------------------------------------------
    design_sets = [set(d["contacts_auth"]) for d in sig["per_design"]]
    full = consensus(design_sets, sig["core_threshold"])
    conv = []
    for k in (2, 4, 6, 8, 12, 16, 20, len(design_sets)):
        if k > len(design_sets):
            continue
        js, sizes = [], []
        for _ in range(60):
            pick = rng.choice(len(design_sets), size=k, replace=False)
            sub = consensus([design_sets[i] for i in pick], sig["core_threshold"])
            sizes.append(len(sub))
            js.append(len(sub & full) / len(sub | full) if (sub | full) else 0.0)
        conv.append({"n_designs": k, "mean_jaccard_to_full": float(np.mean(js)),
                     "sd": float(np.std(js)), "mean_core_size": float(np.mean(sizes))})
    print(f"\n=== C5 convergence: how fast does the consensus stabilise? ===")
    for c in conv:
        print(f"  N={c['n_designs']:3d}  Jaccard to the {len(design_sets)}-design "
              f"consensus {c['mean_jaccard_to_full']:.3f} +/- {c['sd']:.3f}   "
              f"core size {c['mean_core_size']:.1f}")

    out = {
        "ablation": "I6.1 (PROJECT_GOAL.md 8.3) - the M2 gate",
        "question": "does a BoltzGen design consensus beat pocket geometry at predicting "
                    "where a real ligand binds?",
        "ground_truth": "known-ligand contact residues, 4.5 A heavy atom, per structure",
        "n_structures": n, "skipped": skipped,
        "target": {"symbol": "KDR", "uniprot": "P35968",
                   "design_structure": "3VHE chain A (ligand-stripped)"},
        "boltzgen": {"n_designs": sig["n_designs"], "workflow": sig["workflow_uuid"],
                     "core_threshold": sig["core_threshold"],
                     "n_core": len(boltz_core), "core_auth": sorted(boltz_core),
                     "credits": sig["run"]["credits_charged"],
                     "mean_pairwise_jaccard_between_designs":
                         sig["mean_pairwise_jaccard_between_designs"]},
        "arms": summary, "comparisons": comps,
        "convergence_c5": conv,
        "per_structure": per_structure,
        "caveats": [
            "the BoltzGen consensus is built on one structure (3VHE chain A) and scored "
            "against 38 co-crystals of the same protein; author numbering is shared "
            "across KDR entries, which is what makes that comparison legitimate",
            "n = 38 held-out ligands but ONE target, so this does not generalise across "
            "targets - PROJECT_GOAL.md's I6.1 anticipates several",
            "ipTM of every design is far below the 0.85 the plan suggests filtering at, "
            "so this is an UNFILTERED consensus; the plan's filter could not be applied "
            "without generating far more designs than the credit budget allowed",
        ],
    }
    dest = ROOT / "results" / "m2_gate.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
