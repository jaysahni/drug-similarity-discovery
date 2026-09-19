"""Export the computed results as a visualization handoff bundle.

visualization/AGENTS.md defines the integration boundary between the scientific
pipeline and the viewer: a JSON manifest validated against
visualization/schema.json (1.0.0), plus immutable structure files referenced by
bundle-relative path and SHA-256. That subtree belongs to the visualization
worker; this script only PRODUCES a bundle for it and writes nothing inside it.

Everything here is read from results/. Nothing is recomputed, and nothing that
was not computed is invented: where a field has no measured value it is emitted
with status "not_evaluated" and the reason, which is what the schema's
status/reason pairs exist for.

The signature exported is the BoltzGen consensus, because that is the arm
PROJECT_GOAL.md 4.3 specifies - even though m2_gate.py shows it LOSES to pocket
geometry. Showing the arm that lost, labelled, is the point.

Usage:
    ./env/bin/python scripts/build_viz_manifest.py
    ./env/bin/python scripts/build_viz_manifest.py --out results/viz-bundle
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
PIPE = RES / "pipeline" / "colorectal-cancer"
SCHEMA_VERSION = "1.0.0"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def residue_map(pdb_path, structure_id, chain_filter=None):
    """key -> {chain, number, insertion} for every polymer residue in the file."""
    out, seen = {}, set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip() or "A"
        if chain_filter and chain != chain_filter:
            continue
        try:
            num = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        key = f"{structure_id}:{chain}:{num}{icode}"
        if key in seen:
            continue
        seen.add(key)
        out[key] = {"chain": chain, "number": num, "insertion": icode}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RES / "viz-bundle"))
    args = ap.parse_args()
    out = Path(args.out)
    (out / "structures").mkdir(parents=True, exist_ok=True)

    gate = json.loads((RES / "m2_gate.json").read_text())
    sig4 = json.loads((RES / "boltzgen_signature_v4.json").read_text())
    pocket = json.loads((PIPE / "target" / "pocket.json").read_text())
    cand = json.loads((PIPE / "candidates.json").read_text())
    bench = json.loads((RES / "benchmark_targets.json").read_text())
    hotspot = json.loads((RES / "hotspot_recovery.json").read_text()) \
        if (RES / "hotspot_recovery.json").exists() else None

    # structure file, copied in so the bundle is self-contained
    src = PIPE / "target" / "structure" / "3VHE_A_stripped.pdb"
    dest = out / "structures" / "3VHE_A_stripped.pdb"
    shutil.copyfile(src, dest)
    sid = "3VHE_A"
    rmap = residue_map(dest, sid, chain_filter="A")

    # signature residues: engagement frequency across the 24 designs
    design_sets = [set(d["contacts_auth"]) for d in sig4["per_design"]]
    freq = Counter(r for s in design_sets for r in s)
    n_des = len(design_sets)
    core = set(sig4["core_residue_ids_auth"])
    residues = []
    for r, k in sorted(freq.items()):
        key = f"{sid}:A:{r}"
        if key not in rmap:            # never emit a residue the viewer cannot resolve
            continue
        residues.append({"key": key, "label": f"KDR A:{r}",
                         "frequency": round(k / n_des, 3), "core": r in core})

    def val(label, status, reason, value, n, source):
        return {"label": label, "status": status, "reason": reason,
                "value": value, "n": n, "source": source}

    g = gate["comparisons"]["boltzgen_consensus - p2rank_geometry [jaccard]"]
    validation = [
        val("M2 gate (I6.1): BoltzGen consensus vs pocket geometry, Jaccard delta",
            "available",
            f"BoltzGen LOSES. delta {g['delta']:+.4f} [{g['ci_lo']:+.4f}, {g['ci_hi']:+.4f}], "
            f"Holm p={g['p_holm']:.3g}, paired over held-out ligands of one target",
            round(g["delta"], 4), gate["n_structures"], "results/m2_gate.json"),
        val("Hotspot recovery, BoltzGen consensus (Jaccard)", "available",
            "consensus of 24 unfiltered designs against known-ligand contacts",
            round(gate["arms"]["boltzgen_consensus"]["jaccard"], 4),
            gate["n_structures"], "results/m2_gate.json"),
        val("Hotspot recovery, P2Rank top-1 pocket (Jaccard)", "available",
            "same ground truth, pocket predicted blind on the ligand-stripped structure",
            round(gate["arms"]["p2rank_geometry"]["jaccard"], 4),
            gate["n_structures"], "results/m2_gate.json"),
        val("Known-ligand ceiling (Jaccard)", "available",
            "leave-one-out consensus over the other co-crystals; ablation I6.2",
            round(gate["arms"]["known_ligand"]["jaccard"], 4),
            gate["n_structures"], "results/m2_gate.json"),
        val("Random floor (Jaccard)", "available", "size-matched random residues, 100 draws",
            round(gate["arms"]["random"]["jaccard"], 4), gate["n_structures"],
            "results/m2_gate.json"),
        val("Best retrieval representation (p@1)", "available",
            "16 representations benchmarked; nothing beat ECFP4",
            round(bench["representations"]["morgan"]["p@1"], 4),
            bench["_meta"]["n_queries"], "results/benchmark_targets.json"),
        val("Rediscovery enrichment of known KDR binders", "available",
            "top-25 by ECFP4 similarity to the best literature ligand vs base rate; "
            "note sunitinib ranks 464 and pazopanib 764",
            round(cand["enrichment_factor"], 2), cand["top_k"],
            "results/pipeline/colorectal-cancer/candidates.json"),
        val("Interface-coverage-ranked repurposing candidates", "not_evaluated",
            "no approved drug was co-folded against this target, so none has a pose or a "
            "coverage score. The similarity shortlist in candidates.json is a LIGAND-side "
            "ordering and is deliberately not presented as a ranked interface board",
            None, 0, "none"),
        val("Predicted binding affinity", "not_evaluated",
            "no affinity model is run anywhere in this pipeline; PROJECT_GOAL.md G8 "
            "bans affinity language around these numbers", None, 0, "none"),
        val("Small-molecule addressability of the hotspot", "not_evaluated",
            "PROJECT_GOAL.md 1.4b sm_addressable was not computed", None, 0, "none"),
    ]
    if hotspot:
        h = hotspot.get("headline_consensus_vs_p2rank_top1") or {}
        if h:
            validation.insert(1, val(
                "Ablation I6.2 at scale: known-ligand consensus vs pocket (precision delta)",
                "available",
                f"delta {h['delta']:+.4f}, Wilcoxon p={h.get('wilcoxon_p'):.3g}; this is the "
                "CEILING ablation, not the M2 gate",
                round(h["delta"], 4), h["n"], "results/hotspot_recovery.json"))

    # NO CANDIDATES ARE EXPORTED, and that is a decision rather than an omission.
    # The schema requires every candidate to carry a `rank`, and it rejects a ranked
    # candidate that has no interface coverage and no passed confidence gate. We have
    # neither: co-folding was never run, so no approved drug has a pose against this
    # target, so nothing has a coverage score. What we do have is a chemical-similarity
    # ordering (results/.../candidates.json, axitinib-similarity, enrichment 10.2x).
    # Emitting that as the ranked board would present a ligand-side ordering as if it
    # were a target-side one - exactly the conflation this project exists to test.
    # The similarity shortlist stays in its own JSON and is summarised in `validation`.
    candidates = []

    evidence = []
    ev_md = PIPE / "evidence.md"
    if ev_md.exists():
        import re
        for pmid, title in re.findall(r"PMID[:\s]*(\d+)[^\n]*\n?[^\n]*?([A-Z][^\n]{20,120})",
                                      ev_md.read_text())[:8]:
            evidence.append({"id": f"PMID:{pmid}", "title": title.strip()[:160],
                             "claim": "literature evidence linking this target to the disease; "
                                      "every PMID was re-fetched and its title verified",
                             "url": f"https://europepmc.org/article/MED/{pmid}"})
    if not evidence:
        evidence = [{"id": "none", "title": "No literature evidence exported",
                     "claim": "evidence.md was not parsed into the bundle",
                     "url": None}]

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_id": "kdr-m2-gate",
        "data_origin": "computed",
        "disease": {"id": "EFO:colorectal-carcinoma", "name": "Colorectal cancer"},
        "target": {"id": "P35968", "name": "KDR / VEGFR2", "structure_id": sid,
                   "site_kind": "pocket",
                   "rationale": "selected from Open Targets association for colorectal "
                                "cancer; approved VEGFR2 inhibitors exist, which is what "
                                "makes the rediscovery check interpretable"},
        "design_reference": {
            "id": sig4["workflow_uuid"],
            "name": f"BoltzGen protein-anything, {sig4['n_designs']} designs "
                    f"({sig4['run']['credits_charged']} credits, A100-80GB)",
            "structure_id": sid,
            "protocol": "protein-anything",
            "source": "BoltzGen via Rowan protein_binder_design; targeting required "
                      "include_proximity - binding_types alone left designs 28-52 A away",
        },
        "signature": {
            "id": "kdr-boltzgen-consensus-v4",
            "source": "boltzgen_consensus",
            "status": "available",
            "reason": "consensus over 24 BoltzGen designs. UNFILTERED: zero designs reach "
                      "ipTM 0.5, against the 0.85 the plan suggests. This arm LOSES the "
                      "M2 gate to pocket geometry by 0.324 Jaccard.",
            "threshold": sig4["core_threshold"],
            "n_generated": sig4["n_designs"],
            "n_surviving": sig4["n_designs"],
            "addressability": {"status": "not_evaluated",
                               "reason": "PROJECT_GOAL.md 1.4b sm_addressable was not computed"},
            "residues": residues,
        },
        "structures": [{
            "id": sid, "path": "structures/3VHE_A_stripped.pdb",
            "sha256": sha256(dest), "format": "pdb", "origin": "experimental",
            "frame_id": "3VHE-A", "license": "PDB - free to use",
            "source": "RCSB PDB 3VHE chain A, all heteroatoms stripped before site detection",
            "residue_map": rmap,
            "ligand_residues": [],
        }],
        "candidates": candidates,
        "evidence": evidence,
        "validation": validation,
        "provenance": {
            "producer": "scripts/build_viz_manifest.py",
            "purpose": "computed results from the drug-similarity-discovery pipeline",
            "geometry": "RCSB experimental coordinates; contacts at 4.5 A heavy atom",
            "scientific_status": "M2 gate answered and NEGATIVE: the BoltzGen consensus "
                                 "loses to a 0.47 s pocket finder. All outputs are "
                                 "computational hypotheses requiring experimental "
                                 "validation. No affinity is predicted anywhere.",
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"wrote {out/'manifest.json'}")
    print(f"  signature residues: {len(residues)} ({sum(r['core'] for r in residues)} core)")
    print(f"  candidates: {len(candidates)} (none: no poses were computed) | "
          f"evidence: {len(evidence)} | "
          f"validation rows: {len(validation)} "
          f"({sum(1 for v in validation if v['status']=='not_evaluated')} not_evaluated)")
    print(f"  structure sha256: {manifest['structures'][0]['sha256'][:16]}...")


if __name__ == "__main__":
    main()
