import hashlib
import json
import math
from pathlib import Path


def make_fixture(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    mapping, lines = {}, []
    serial = 0
    for i in range(1, 13):
        angle = i * 1.75
        x, y, z = 3 * math.cos(angle), 3 * math.sin(angle), i * 1.5
        for atom, dx, dy, dz, element in [("N", -0.6, 0, -0.5, "N"), ("CA", 0, 0, 0, "C"), ("C", 0.6, 0, 0.5, "C"), ("O", 1.1, 0.5, 0.7, "O"), ("CB", 0, 1.2, 0, "C")]:
            serial += 1
            lines.append(f"ATOM  {serial:5d} {atom:^4} ALA A{i:4d}    {x+dx:8.3f}{y+dy:8.3f}{z+dz:8.3f}{1:6.2f}{20:6.2f}          {element:>2}")
        mapping[f"toy:A:{i}"] = {"chain": "A", "number": i, "insertion": ""}
    for i in range(6):
        serial += 1
        angle = i * math.pi / 3
        x, y, z = 5 + 1.4 * math.cos(angle), 1.4 * math.sin(angle), 10
        lines.append(f"HETATM{serial:5d} {('C'+str(i+1)):^4} TOY L   1    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{20:6.2f}           C")
    pdb = "\n".join(lines) + "\nEND\n"
    (destination / "toy.pdb").write_text(pdb)
    residues = [{"key": f"toy:A:{i}", "label": f"Toy A:{i}", "frequency": f, "core": f >= 0.6} for i, f in [(4, 0.4), (5, 0.8), (6, 1), (7, 0.9), (8, 0.7), (9, 0.3)]]
    core = [r["key"] for r in residues if r["core"]]
    candidates = []
    for i, (name, covered, gate, modality) in enumerate([("Fixture alpha", 3, "passed", "small_molecule"), ("Fixture beta", 2, "passed", "small_molecule"), ("Fixture gamma", 1, "failed", "small_molecule"), ("Fixture peptide", 4, "passed", "peptide")]):
        candidates.append({"id": f"fixture-{i}", "name": name, "modality": modality,
            "novelty": ["novel_pairing", "known_moa", "known_offtarget", "unknown"][i], "rank": (i + 1 if i < 2 else 1) if gate == "passed" else None,
            "ranking_metric": "core_coverage (synthetic test ordering)", "structure_id": "toy", "approval": "Synthetic fixture — not an approved drug",
            "coverage": {"status": "available", "value": covered / 4, "reason": "Synthetic UI test only"}, "engaged": core[:covered], "missed": core[covered:],
            "confidence": {"name": "fixture score", "value": 0.9 if gate == "passed" else 0.2, "scale": "0–1, synthetic", "gate": gate, "rule": "Synthetic threshold 0.5", "reason": "Not model output"},
            "decoys": {"status": "not_evaluated", "reason": "No scientific null was computed", "percentile": None, "scores": [], "n": 0, "metric": "core_coverage", "context": f"toy-target:toy-signature:{modality}", "tie_rule": "not evaluated"},
            "caveats": ["All candidates intentionally share a toy pose; differences in coverage are fabricated test cases, not structure-derived measurements."], "source": "synthetic_fixture"})
    manifest = {"schema_version": "1.0.0", "run_id": "ui-fixture", "data_origin": "synthetic_fixture",
        "disease": {"id": "fixture-disease", "name": "A research question"},
        "target": {"id": "toy-target", "name": "Illustrative target", "structure_id": "toy", "site_kind": "pocket", "rationale": "A synthetic interface for testing the presentation. This is not a disease mechanism or a real protein."},
        "signature": {"id": "toy-signature", "source": "known_ligand", "status": "available", "reason": "Synthetic interface, not a scientific result", "threshold": 0.6, "n_generated": 0, "n_surviving": 0, "addressability": {"status": "addressable", "reason": "Synthetic state used to exercise the small-molecule board"}, "residues": residues},
        "structures": [{"id": "toy", "path": "toy.pdb", "sha256": hashlib.sha256(pdb.encode()).hexdigest(), "format": "pdb", "origin": "synthetic_fixture", "frame_id": "toy-frame", "source": "Procedural toy coordinates; no biological meaning", "license": "MIT", "residue_map": mapping, "ligand_residues": [{"chain": "L", "number": 1, "insertion": ""}]}],
        "candidates": candidates, "evidence": [{"id": "fixture-not-a-citation", "title": "Evidence will appear here", "claim": "Real pipeline exports must supply reviewed claims and source identifiers. No paper or PMID is invented for this fixture.", "url": None}],
        "validation": [{"label": "Hotspot recovery", "status": "not_evaluated", "reason": "Requires scientific pipeline results", "value": None, "n": 0, "source": "none"}, {"label": "BoltzGen vs pocket geometry", "status": "not_evaluated", "reason": "No ablation has been run", "value": None, "n": 0, "source": "none"}],
        "provenance": {"producer": "visualization.fixture", "purpose": "UI CONTRACT TEST ONLY", "geometry": "Procedurally constructed toy atoms", "scientific_status": "Not evaluated"}}
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return destination / "manifest.json"
