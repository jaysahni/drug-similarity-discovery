"""Turn a Rowan BoltzGen design run into an InterfaceSignature, and score its targeting.

This is the `source="boltzgen_consensus"` builder of PROJECT_GOAL.md 4.3 / C8 - the
arm that the P2Rank and known-ligand arms exist to be compared against (ablation
I6.1, the M2 gate of section 8.3).

It downloads each designed complex, takes the target-side residues within CUTOFF of
the designed chain, aggregates them into a per-residue engagement frequency, and
reports:
  - whether the designs actually bound WHERE THEY WERE ASKED TO (overlap with the
    requested pocket). PROJECT_GOAL.md 1.2 warns that "BoltzGen doesn't consistently
    demonstrate that designs bind where predicted", so this is measured, not assumed
  - how much the designs agree with EACH OTHER (pairwise Jaccard). A consensus over
    designs that disagree is not a signature of anything

RESIDUE NUMBERING. BoltzGen indexes residues 1..N over the residues PRESENT in the
input file, not by author numbering - verified the hard way: a first run using author
numbering (840..1047 against a 303-residue file) failed with "BoltzGen exited with
code 1". target/boltzgen_residue_map.json holds the mapping both ways.

Usage:
    ./env/bin/python scripts/boltzgen_signature.py <workflow_uuid> [--tag v3]
"""

from __future__ import annotations

import argparse
import json
import os
import warnings
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "results" / "pipeline" / "colorectal-cancer"
CUTOFF = 4.5           # same cutoff interfaces.py uses, so the arms are comparable
CORE_THRESHOLD = 0.6   # PROJECT_GOAL.md 4.3 default


def load_env():
    env = ROOT / ".env"
    if env.exists():
        for line in env.read_text().splitlines():
            if "=" in line and not line.strip().startswith("#"):
                k, _, v = line.partition("=")
                os.environ.setdefault(k.strip(), v.strip())


def download_designs(uuid, outdir):
    import rowan
    load_env()
    key = os.environ.get("ROWAN_API_KEY")
    if not key:
        raise SystemExit("no ROWAN_API_KEY - put it in .env (see .env.example)")
    rowan.api_key = key
    wf = rowan.retrieve_workflow(uuid)
    if "OK" not in str(wf.status).upper():
        raise SystemExit(f"workflow {uuid} status is {wf.status}, not completed")
    outdir.mkdir(parents=True, exist_ok=True)
    out = []
    for i, b in enumerate(wf.data["generated_binders"]):
        d = outdir / f"design_{i:03d}"
        existing = list(d.glob("*.pdb"))
        if not existing:
            rowan.retrieve_protein(b["bound_structure"]).download_pdb_file(str(d))
            existing = list(d.glob("*.pdb"))
        out.append({"path": existing[0], "scores": b["scores"],
                    "sequence": (b.get("binder_sequences") or [{}])[0].get("sequence")})
    return out, {"credits_charged": wf.credits_charged, "elapsed_s": wf.elapsed,
                 "status": str(wf.status), "settings": wf.data.get("binder_design_settings")}


def target_contacts(path, target_chain="A", design_chain="B", cutoff=CUTOFF):
    """Target-side residue ids within `cutoff` of any design heavy atom."""
    from Bio.PDB import PDBParser
    model = PDBParser(QUIET=True).get_structure("d", str(path))[0]
    tgt = [r for r in model[target_chain] if r.id[0] == " "]
    des = np.array([a.coord for r in model[design_chain] if r.id[0] == " "
                    for a in r if a.element != "H"])
    if not des.size:
        return set(), None
    hits = set()
    for r in tgt:
        ra = np.array([a.coord for a in r if a.element != "H"])
        if ra.size and np.min(np.linalg.norm(ra[:, None, :] - des[None, :, :], axis=-1)) <= cutoff:
            hits.add(r.id[1])
    return hits, des.mean(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("uuid")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--core-threshold", type=float, default=CORE_THRESHOLD)
    args = ap.parse_args()

    rmap = json.loads((PIPE / "target" / "boltzgen_residue_map.json").read_text())
    idx2auth = {v: int(k) for k, v in rmap["auth_to_index"].items()}
    pocket = set(rmap["pocket_boltzgen_idx"])
    pocket_centre = np.array(json.loads((PIPE / "target" / "pocket.json").read_text())["center"])

    designs, meta = download_designs(args.uuid, PIPE / "designs" / args.tag)
    print(f"{len(designs)} designs  ({meta['credits_charged']} credits, "
          f"{meta['elapsed_s']:.0f}s, {meta['status']})\n")

    rows, sets = [], []
    for i, d in enumerate(designs):
        hits, centroid = target_contacts(d["path"])
        sets.append(hits)
        overlap = hits & pocket
        dist = float(np.linalg.norm(centroid - pocket_centre)) if centroid is not None else None
        rows.append({
            "design": i, "n_contacts": len(hits),
            "overlap_with_requested_pocket": len(overlap),
            "pocket_size": len(pocket),
            "binder_centroid_dist_to_pocket_a": round(dist, 1) if dist else None,
            "iptm": d["scores"].get("iptm"),
            "design_ptm": d["scores"].get("design_ptm"),
            "num_filters_passed": d["scores"].get("num_filters_passed"),
            "contacts_auth": sorted(idx2auth.get(h, -1) for h in hits),
        })
        print(f"  design {i}: {len(hits):3d} contacts | {len(overlap):2d}/{len(pocket)} on the "
              f"requested pocket | {dist:5.1f} A away | ipTM {d['scores'].get('iptm')}")

    n = len(sets)
    freq = Counter(x for s in sets for x in s)
    core = sorted(r for r, k in freq.items() if k / n >= args.core_threshold)
    pairwise = [len(a & b) / len(a | b) if (a | b) else float("nan") for a, b in combinations(sets, 2)]

    out = {
        "source": "boltzgen_consensus",
        "workflow_uuid": args.uuid, "tag": args.tag,
        "n_designs": n, "cutoff_a": CUTOFF, "core_threshold": args.core_threshold,
        "run": meta,
        "requested_pocket_auth": sorted(idx2auth[i] for i in pocket),
        "per_design": rows,
        "core_residue_ids_idx": core,
        "core_residue_ids_auth": sorted(idx2auth.get(r, -1) for r in core),
        "n_core": len(core),
        "core_overlap_with_requested_pocket": len(set(core) & pocket),
        "mean_pairwise_jaccard_between_designs": (
            float(np.nanmean(pairwise)) if pairwise else None),
        "designs_agreeing_on_any_residue": sum(1 for r, k in freq.items() if k > 1),
    }
    dest = ROOT / "results" / f"boltzgen_signature_{args.tag}.json"
    dest.write_text(json.dumps(out, indent=1))

    print(f"\nconsensus at core_threshold={args.core_threshold}: {len(core)} core residues")
    print(f"  core on requested pocket: {out['core_overlap_with_requested_pocket']}")
    print(f"  mean pairwise Jaccard between designs: "
          f"{out['mean_pairwise_jaccard_between_designs']:.3f}")
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
