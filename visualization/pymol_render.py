import json
import os
from pathlib import Path


def render(bundle_path, output):
    from pymol import cmd

    bundle_path, output = Path(bundle_path).resolve(), Path(output).resolve()
    bundle = json.loads(bundle_path.read_text())
    manifest, colors = bundle["manifest"], bundle["colors"]
    output.mkdir(parents=True, exist_ok=False)
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", 1)
    cmd.set("label_color", "black")
    cmd.set("label_size", 16)
    for name, color in colors.items():
        cmd.set_color("viz_" + name, [int(color[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    objects = {}
    for i, (aid, asset) in enumerate(bundle["display"].items()):
        name = f"structure_{i}"
        cmd.read_pdbstr(asset["pdb"], name)
        objects[aid] = name
    target = manifest["target"]["structure_id"]
    candidates = manifest["candidates"][:5]
    design = manifest.get("design_reference")
    choices = [None] + ([dict(design, _design=True, engaged=[], missed=[])] if design else []) + candidates
    produced = []
    for candidate in choices:
        aid = candidate["structure_id"] if candidate else target
        if not aid or aid not in objects:
            continue
        name, asset = objects[aid], bundle["display"][aid]
        is_design = bool(candidate and candidate.get("_design"))
        modes = ("design",) if is_design else ("coverage", "missed") if candidate else ("context", "site", "frequency")
        selected_scenes = [s for s in bundle["scenes"] if s["mode"] in modes]
        for scene in selected_scenes:
            cmd.disable("all")
            cmd.enable(name)
            cmd.hide("everything", "all")
            cmd.color("viz_target", name)
            cmd.show("cartoon", name)
            cmd.show("lines", name)
            cmd.set("transparency", 0.4, name)
            core_atoms = []
            for residue in manifest["signature"]["residues"]:
                atoms = asset["residue_atoms"].get(residue["key"], [])
                if not atoms:
                    continue
                selection = f"{name} and id " + "+".join(map(str, atoms))
                color = "viz_target"
                if scene["mode"] in ("frequency", "design"):
                    frequency = residue["frequency"]
                    color = f"freq_{len(produced)}_{atoms[0]}"
                    cmd.set_color(color, [0.92 * (1 - frequency), 0.96 - 0.46 * frequency, 0.96 - 0.51 * frequency])
                elif candidate:
                    color = "viz_engaged" if residue["key"] in candidate["engaged"] else "viz_missed" if residue["key"] in candidate["missed"] else "viz_target"
                elif scene["mode"] == "site":
                    color = "viz_engaged"
                cmd.color(color, selection)
                if scene["mode"] != "context":
                    cmd.show("sticks", selection)
                    cmd.show("surface", selection)
                if residue["core"]:
                    core_atoms.extend(atoms)
            ligand = asset["ligand_atoms"]
            if ligand:
                selection = f"{name} and id " + "+".join(map(str, ligand))
                cmd.hide("everything", selection)
                if candidate:
                    cmd.show("sticks", selection)
                    if is_design:
                        cmd.show("cartoon", selection)
                    cmd.color("viz_design" if is_design else "viz_ligand", selection)
            cmd.orient(name)
            if core_atoms and scene["mode"] != "context":
                cmd.zoom(f"{name} and id " + "+".join(map(str, core_atoms)), buffer=5)
            cmd.delete("caption")
            center = cmd.get_position()
            origin = manifest["data_origin"].upper().replace("_", " ")
            label = f"{origin} | {scene['title']}" + (f" | {candidate['name']}" if candidate else "")
            cmd.pseudoatom("caption", pos=center, label=label)
            cmd.hide("everything", "caption")
            cmd.show("labels", "caption")
            cmd.set("label_position", (0, 5, 0), "caption")
            key = scene["id"] + ("_" + candidate["id"] if candidate else "")
            cmd.scene(key, "store")
            path = output / (key + ".png")
            cmd.png(str(path), width=1920, height=1080, ray=1)
            produced.append({"scene": key, "candidate_id": candidate["id"] if candidate else None, "file": path.name,
                             "engaged": candidate["engaged"] if candidate else [], "missed": candidate["missed"] if candidate else [], "caption": label, "view": list(cmd.get_view())})
    cmd.save(str(output / "presentation.pse"))
    (output / "render_manifest.json").write_text(json.dumps({"pymol_version": str(cmd.get_version()), "run_id": manifest["run_id"], "origin": manifest["data_origin"], "disclaimer": bundle["disclaimer"], "scenes": produced}, indent=2))


if os.environ.get("AUTOREPURPOSE_VIZ_BUNDLE"):
    render(os.environ["AUTOREPURPOSE_VIZ_BUNDLE"], os.environ["AUTOREPURPOSE_VIZ_RENDER_OUT"])
