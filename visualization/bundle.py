import copy
import hashlib
import io
import json
import math
import shutil
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBException
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parent
COLORS = {"target": "#cbd5e1", "engaged": "#007f73", "missed": "#b45309", "ligand": "#2563eb", "selected": "#7c3aed", "design": "#7c3aed"}
DISCLAIMER = "Computationally ranked hypotheses requiring experimental validation."


from .errors import BundleError  # re-exported: callers have always imported it from here


def load_json(path):
    def reject(value):
        raise BundleError(f"Non-finite JSON value: {value}")
    return json.loads(Path(path).read_text(), parse_constant=reject)


def safe_asset(root, path):
    root = Path(root).resolve()
    candidate = Path(path)
    if candidate.is_absolute() or ".." in candidate.parts or ":" in path or "\\" in path:
        raise BundleError(f"Asset must be a relative bundle path: {path}")
    resolved = (root / candidate).resolve()
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise BundleError(f"Missing or unsafe asset: {path}")
    if resolved.stat().st_size > 30_000_000:
        raise BundleError(f"Asset exceeds 30 MB limit: {path}")
    return resolved


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unique(records, name):
    values = [r[name] for r in records]
    if len(values) != len(set(values)):
        raise BundleError(f"Duplicate {name}")
    return set(values)


def validate(manifest, root):
    schema = load_json(ROOT / "schema.json")
    errors = sorted(Draft202012Validator(schema).iter_errors(manifest), key=lambda e: str(e.path))
    if errors:
        raise BundleError("; ".join(f"{list(e.path)}: {e.message}" for e in errors[:8]))
    for value in walk(manifest):
        if isinstance(value, float) and not math.isfinite(value):
            raise BundleError("Non-finite number")
    assets = unique(manifest["structures"], "id")
    unique(manifest["candidates"], "id")
    signature = manifest["signature"]
    unique(signature["residues"], "key")
    core = {r["key"] for r in signature["residues"] if r["core"]}
    if signature["n_surviving"] > signature["n_generated"]:
        raise BundleError("Surviving count exceeds generated count")
    if signature["source"] == "boltzgen_consensus" and signature["status"] == "available" and not signature["n_surviving"]:
        raise BundleError("Available BoltzGen consensus requires surviving designs")
    if manifest["target"]["structure_id"] is not None and manifest["target"]["structure_id"] not in assets:
        raise BundleError("Unknown target structure")
    design = manifest.get("design_reference")
    if design and design["structure_id"] not in assets:
        raise BundleError("Unknown design-reference structure")
    for asset in manifest["structures"]:
        path = safe_asset(root, asset["path"])
        if digest(path) != asset["sha256"]:
            raise BundleError(f"Checksum mismatch: {asset['id']}")
        if manifest["data_origin"] == "computed" and asset["origin"] == "synthetic_fixture":
            raise BundleError("Synthetic structure cannot appear in a computed bundle")
        alignment = asset.get("alignment")
        if alignment and (alignment["reference_id"] not in assets or alignment["reference_id"] == asset["id"]):
            raise BundleError("Invalid alignment reference")
    ranks = set()
    for c in manifest["candidates"]:
        if c["structure_id"] is not None and c["structure_id"] not in assets:
            raise BundleError(f"Unknown candidate structure: {c['id']}")
        engaged, missed = set(c["engaged"]), set(c["missed"])
        available = c["coverage"]["status"] == "available"
        if available:
            if signature["status"] != "available" or not core:
                raise BundleError("Coverage requires an available, nonempty core")
            if engaged & missed or engaged | missed != core:
                raise BundleError("Engaged/missed sets must partition the core")
            value = c["coverage"]["value"]
            if value is None or abs(value - len(engaged) / len(core)) > 1e-6:
                raise BundleError("Coverage disagrees with residue sets")
        elif c["coverage"]["value"] is not None or engaged or missed:
            raise BundleError("Unavailable coverage must have null value and no explanation sets")
        if c["rank"] is not None:
            key = (c["modality"], c["rank"])
            if key in ranks:
                raise BundleError("Duplicate rank within modality")
            ranks.add(key)
            if not available or c["confidence"]["gate"] != "passed":
                raise BundleError("Ranked candidate must pass confidence and have coverage")
            if c["modality"] == "small_molecule" and signature["addressability"]["status"] != "addressable":
                raise BundleError("Small-molecule rank requires assessed addressability")
        if c["confidence"]["gate"] != "not_evaluated" and c["confidence"]["value"] is None:
            raise BundleError("Evaluated confidence gate requires a native value")
        d = c["decoys"]
        if d["status"] == "available":
            expected = f"{manifest['target']['id']}:{signature['id']}:{c['modality']}"
            if d["context"] != expected or d["n"] <= 0 or d["n"] != len(d["scores"]) or d["percentile"] is None or not d["tie_rule"]:
                raise BundleError("Invalid decoy calibration context/count/percentile/tie rule")
        elif d["percentile"] is not None or d["scores"] or d["n"]:
            raise BundleError("Unavailable calibration must have null percentile, zero n and no scores")
    for record in manifest["validation"]:
        if (record["status"] == "available") != (record["value"] is not None):
            raise BundleError("Validation value/status mismatch")
    return manifest


def walk(value):
    yield value
    if isinstance(value, dict):
        for child in value.values():
            yield from walk(child)
    if isinstance(value, list):
        for child in value:
            yield from walk(child)


def residue_id(residue):
    return residue["chain"], residue["number"], residue["insertion"]


def atoms_for(asset, root):
    parser = PDBParser(PERMISSIVE=False, QUIET=True) if asset["format"] == "pdb" else MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(asset["id"], str(safe_asset(root, asset["path"])))
    except (PDBConstructionException, PDBException, ValueError, KeyError) as exc:
        raise BundleError(f"Structure {asset['id']} could not be parsed: {exc}") from exc
    models = list(structure.get_models())
    if len(models) != 1:
        raise BundleError("Select a single pose/model upstream before exporting")
    positions = [(r.get_parent().id, r.id[1], r.id[2].strip()) for r in models[0].get_residues()]
    if len(positions) != len(set(positions)):
        raise BundleError("Ambiguous polymer/non-polymer residue numbering; normalize upstream")
    mapped = [residue_id(r) for r in asset["residue_map"].values()]
    if len(mapped) != len(set(mapped)):
        raise BundleError("Canonical residue mapping must be one-to-one")
    if set(mapped) & {residue_id(r) for r in asset["ligand_residues"]}:
        raise BundleError("Target and ligand residue selectors must not overlap")
    atoms = []
    for atom in models[0].get_atoms():
        residue = atom.get_parent()
        atoms.append({"source_serial": int(atom.serial_number), "xyz": np.array(atom.coord, dtype=float), "name": atom.name,
                      "element": atom.element, "resname": residue.resname,
                      "residue": (residue.get_parent().id, int(residue.id[1]), residue.id[2].strip()),
                      "bfactor": float(atom.bfactor), "hetero": residue.id[0] != " "})
    if not atoms or len(atoms) > 99999:
        raise BundleError("Structure must contain 1–99999 atoms for normalized display PDB")
    if len({a["source_serial"] for a in atoms}) != len(atoms):
        raise BundleError("Source atom serials must be unique")
    if not np.isfinite([a["xyz"] for a in atoms]).all():
        raise BundleError("Non-finite structure coordinates")
    for a in atoms:
        if not a["element"] or not a["element"].isalpha() or len(a["element"]) > 2:
            raise BundleError("Unknown atomic element; normalize the structure upstream")
    return atoms


def fit(mobile, reference, max_rmsd):
    mobile, reference = np.asarray(mobile), np.asarray(reference)
    if np.linalg.matrix_rank(mobile - mobile.mean(0)) < 2 or np.linalg.matrix_rank(reference - reference.mean(0)) < 2:
        raise BundleError("Alignment atoms are collinear or coincident")
    left, _, right = np.linalg.svd((mobile - mobile.mean(0)).T @ (reference - reference.mean(0)))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    rotation = (left @ correction @ right).T
    translation = reference.mean(0) - rotation @ mobile.mean(0)
    aligned = mobile @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return rotation, translation, rmsd, rmsd <= max_rmsd


def normalize(asset, atoms):
    residues = list(dict.fromkeys(a["residue"] for a in atoms))
    if len(residues) > 9999:
        raise BundleError("Display normalization supports at most 9999 residues")
    numbers = {r: i + 1 for i, r in enumerate(residues)}
    chains = list(dict.fromkeys(r[0] for r in residues))
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(chains) > len(alphabet):
        raise BundleError("Display normalization supports at most 62 chains")
    chain_names = dict(zip(chains, alphabet))
    selection = {r: [] for r in residues}
    lines = []
    for i, atom in enumerate(atoms, 1):
        selection[atom["residue"]].append(i)
        x, y, z = atom["xyz"]
        if max(abs(x), abs(y), abs(z)) > 999.0:
            raise BundleError("Coordinates outside normalized PDB range")
        record = "HETATM" if atom["hetero"] else "ATOM"
        chain = chain_names[atom["residue"][0]]
        lines.append(f"{record:<6}{i:5d} {atom['name'][:4]:^4} {atom['resname'][:3]:>3} {chain}{numbers[atom['residue']]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{atom['bfactor']:6.2f}          {atom['element']:>2}")
    mapped = {}
    for key, residue in asset["residue_map"].items():
        if residue_id(residue) not in selection:
            raise BundleError(f"Residue mapping {key} not present in {asset['id']}")
        mapped[key] = selection[residue_id(residue)]
    ligand = []
    for residue in asset["ligand_residues"]:
        if residue_id(residue) not in selection:
            raise BundleError(f"Ligand selector absent in {asset['id']}")
        ligand.extend(selection[residue_id(residue)])
    return "\n".join(lines) + "\nEND\n", mapped, ligand


def build(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    manifest = validate(load_json(source), source.parent)
    if output.exists():
        raise BundleError("Output already exists; choose a new export directory")
    if output == source.parent or source.parent.is_relative_to(output):
        raise BundleError("Output must not contain the source directory")
    loaded = {a["id"]: atoms_for(a, source.parent) for a in manifest["structures"]}
    assets = {a["id"]: a for a in manifest["structures"]}
    display = {}
    for asset in manifest["structures"]:
        atoms = copy.deepcopy(loaded[asset["id"]])
        record = {"status": "native", "reason": "Input coordinate frame; no display alignment requested"}
        frame = asset["frame_id"]
        if alignment := asset.get("alignment"):
            ref = assets[alignment["reference_id"]]
            if ref.get("alignment"):
                raise BundleError("Chained alignments are unsupported; use the fixed reference")
            mobile_map = {a["source_serial"]: a for a in atoms}
            ref_map = {a["source_serial"]: a for a in loaded[ref["id"]]}
            mobile_pairs, ref_pairs = zip(*alignment["atom_pairs"])
            if len(set(mobile_pairs)) != len(mobile_pairs) or len(set(ref_pairs)) != len(ref_pairs):
                raise BundleError("Alignment atom correspondences must be one-to-one")
            try:
                mobile_atoms = [mobile_map[n] for n in mobile_pairs]
                ref_atoms = [ref_map[n] for n in ref_pairs]
            except KeyError as exc:
                raise BundleError(f"Unknown alignment atom {exc}") from exc
            mobile_target = {residue_id(r) for r in asset["residue_map"].values()}
            ref_target = {residue_id(r) for r in ref["residue_map"].values()}
            if any(a["residue"] not in mobile_target for a in mobile_atoms) or any(a["residue"] not in ref_target for a in ref_atoms):
                raise BundleError("Align only mapped target atoms, never ligand atoms")
            mobile_keys = {residue_id(r): key for key, r in asset["residue_map"].items()}
            reference_keys = {residue_id(r): key for key, r in ref["residue_map"].items()}
            for a, b in zip(mobile_atoms, ref_atoms):
                if a["name"] not in {"N", "CA", "C", "O"} or a["name"] != b["name"] or mobile_keys[a["residue"]] != reference_keys[b["residue"]]:
                    raise BundleError("Alignment pairs must match canonical target residues and backbone atom names")
            rotation, translation, rmsd, passed = fit([a["xyz"] for a in mobile_atoms], [a["xyz"] for a in ref_atoms], alignment["max_rmsd"])
            record = {"status": "aligned" if passed else "failed", "rmsd": rmsd, "n_atoms": len(mobile_atoms), "max_rmsd": alignment["max_rmsd"], "rotation": rotation.tolist(), "translation": translation.tolist(), "reason": "Target-only rigid fit; entire complex transformed" if passed else "Alignment exceeded configured RMSD; views remain independent"}
            if passed:
                for atom in atoms:
                    atom["xyz"] = rotation @ atom["xyz"] + translation
                frame = ref["frame_id"]
            else:
                frame = "unaligned-" + asset["id"]
        pdb, mapped, ligand = normalize(asset, atoms)
        display[asset["id"]] = {"pdb": pdb, "residue_atoms": mapped, "ligand_atoms": ligand, "frame_id": frame, "alignment": record}
    scenes = [{"id": key, "title": title, "mode": mode} for key, title, mode in [
        ("target", "Target context", "context"), ("site", "Binding site", "site"),
        ("consensus", "Consensus interface", "frequency"), ("candidate", "Candidate coverage", "coverage"),
        ("missed", "Missed core residues", "missed")]]
    if manifest.get("design_reference"):
        scenes.insert(3, {"id": "design", "title": "Representative design (not the consensus)", "mode": "design"})
    warnings = []
    keys = {r["key"] for r in manifest["signature"]["residues"]}
    for aid, asset in display.items():
        missing = sorted(keys - set(asset["residue_atoms"]))
        if missing:
            warnings.append(f"{aid}: unmapped signature residues: {', '.join(missing)}")
    bundle = {"manifest": manifest, "display": display, "scenes": scenes, "colors": COLORS, "warnings": warnings, "disclaimer": "Experimental structure reference; no repurposing prediction or ranking." if manifest["data_origin"] == "reference_example" else DISCLAIMER}
    output.mkdir(parents=True)
    (output / "bundle.json").write_text(json.dumps(bundle, indent=2, allow_nan=False))
    (output / "structures").mkdir()
    for aid, asset in display.items():
        (output / "structures" / f"{aid}.pdb").write_text(asset["pdb"])
    for path in (ROOT / "web").iterdir():
        if path.is_file():
            shutil.copy2(path, output / path.name)
    vendor = ROOT / "vendor"
    if vendor.exists():
        shutil.copytree(vendor, output / "vendor")
    # Co-folded poses for the recorded worked examples travel with every export.
    poses = ROOT / "web" / "poses"
    if poses.exists():
        shutil.copytree(poses, output / "poses")
    from .report import render_report
    (output / "report.html").write_text(render_report(bundle))
    shutil.copy2(ROOT / "pymol_render.py", output / "pymol_render.py")
    checksums = {str(p.relative_to(output)): digest(p) for p in sorted(output.rglob("*")) if p.is_file()}
    (output / "checksums.json").write_text(json.dumps(checksums, indent=2))
    return bundle
