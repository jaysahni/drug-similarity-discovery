"""Does co-folding confidence discriminate binders better than our interface score?

This exists because the README claimed it does not, and that claim was wrong.

It was drawn from a single CONSTRAINED probe: aspirin, forced into the pocket,
co-folded at ipTM 0.965 - as high as the true binder. That is a real observation
about what a pocket constraint does, and it is why constraints stay off. But it
says nothing about the UNCONSTRAINED runs the shipped board is actually built
from, and on those the picture reverses.

Every score here is already on disk - the confidences were returned by Rowan with
each pose and cached in repurpose_state.json - so this costs nothing to rerun.

Usage:
    ./env/bin/python scripts/confidence_gate.py
    ./env/bin/python scripts/confidence_gate.py --pipeline colorectal-cancer
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
POSITIVE = "known_binder"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pipeline", default="colorectal-cancer")
    ap.add_argument("--signature", default=None, help="default: whatever the board ranked by")
    args = ap.parse_args()

    board_p = ROOT / "results" / f"repurpose_{args.pipeline}.json"
    state_p = ROOT / "results" / "pipeline" / args.pipeline / "repurpose_state.json"
    for p in (board_p, state_p):
        if not p.exists():
            sys.exit(f"missing {p}")
    board = json.loads(board_p.read_text())
    state = json.loads(state_p.read_text())

    sig = args.signature or board["ranked_by"]
    metric = board["ranking_metric"]
    role, iface = {}, {}
    for r in board["results"]:
        if r.get("status") != "scored":
            continue
        role[r["struct_id"]] = r["role"]
        iface[r["struct_id"]] = r["by_signature"][sig][metric]

    conf = {sid: (j.get("scores") or {}) for sid, j in state["jobs"].items()
            if sid in role and (j.get("scores") or {}).get("iptm") is not None}
    if not conf:
        sys.exit("no cached co-folding confidences in the state file")
    constrained = sorted({bool(state["jobs"][s].get("constrained")) for s in conf})

    signals = {f"interface overlap ({metric})": {s: iface[s] for s in conf}}
    for key in ("iptm", "confidence_score", "avg_lddt", "ptm"):
        if all(key in conf[s] for s in conf):
            signals[key] = {s: conf[s][key] for s in conf}

    negatives = sorted({role[s] for s in conf} - {POSITIVE})

    def auc(score, neg):
        a = np.array([score[s] for s in conf if role[s] == POSITIVE])
        b = np.array([score[s] for s in conf if role[s] == neg])
        if not a.size or not b.size:
            return None
        u = mannwhitneyu(a, b, alternative="two-sided")
        return {"auc": float(u.statistic / (a.size * b.size)),
                "p_two_sided": float(u.pvalue),
                "n_positive": int(a.size), "n_negative": int(b.size)}

    out = {"what": "co-folding confidence vs interface overlap as a binder discriminator",
           "pipeline": args.pipeline, "signature": sig, "metric": metric,
           "n_drugs": len(conf),
           "pocket_constrained": (constrained[0] if len(constrained) == 1 else "mixed"),
           "results": {}}
    print(f"{len(conf)} drugs with both an interface score and a cached co-fold confidence")
    print(f"co-folding was {'UNCONSTRAINED' if constrained == [False] else 'constrained/mixed'}\n")
    head = f"{'signal':<30}" + "".join(f"{('vs ' + n):>22}" for n in negatives)
    print(head); print("-" * len(head))
    for name, sc in signals.items():
        row = {}
        line = f"{name:<30}"
        for neg in negatives:
            r = auc(sc, neg)
            row[neg] = r
            line += f"{r['auc']:>13.3f} p={r['p_two_sided']:<6.4f}"
        out["results"][name] = row
        print(line)

    prim = f"interface overlap ({metric})"
    easy = [n for n in negatives if n != "hard_decoy"]
    reading = []
    if "iptm" in out["results"] and easy:
        a_i = out["results"]["iptm"][easy[0]]["auc"]
        a_o = out["results"][prim][easy[0]]["auc"]
        if a_i > a_o:
            reading.append(
                f"On the easy null, ipTM ({a_i:.3f}) BEATS interface overlap ({a_o:.3f}). "
                f"A free confidence field returned with every pose is the better "
                f"discriminator there, and the README must not claim otherwise.")
    if "hard_decoy" in negatives:
        hs = {k: v["hard_decoy"] for k, v in out["results"].items()}
        if all(v["p_two_sided"] > 0.05 for v in hs.values()):
            reading.append(
                "On the hard null NOTHING reaches significance - not our score, not any "
                "confidence field. The ceiling is shared with the co-folding backend "
                "rather than caused by the scorer.")
    reading.append(
        "What survives for interface overlap is interpretability, not accuracy: it names "
        "the residues engaged and missed, so a hit can be inspected. ipTM is a single "
        "number with no mechanism attached.")
    out["reading"] = reading
    print()
    for r in reading:
        print("  * " + r)

    dest = ROOT / "results" / f"confidence_gate_{args.pipeline}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
