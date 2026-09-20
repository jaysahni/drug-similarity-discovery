"""Contact extraction from a co-folded complex.

NovaKit's `protein_interactions.analyze_interface` is marked `experimental` and
needs administrator-provisioned Modal engines, which this demo does not have.
Rather than let the toolkit silently fall back to something weaker, we do the one
piece it can't: a plain distance cutoff over the downloaded PDB.

This is deliberately the same definition the heavy pipeline uses in
scripts/interfaces.py -- heavy-atom distance <= 4.5 A -- so the demo's numbers are
comparable to results/m2_gate_CDK2.json rather than being on their own scale.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

CUTOFF = 4.5  # angstroms, heavy atom to heavy atom


def parse_pdb(path: Path) -> tuple[dict, np.ndarray, list]:
    """Return (protein_atoms, ligand_coords, ligand_names) from a PDB file.

    protein_atoms maps (chain, resseq) -> (resname, coords array).
    Hydrogens are dropped: the cutoff is defined on heavy atoms.
    """
    protein: dict[tuple[str, int], list] = {}
    resnames: dict[tuple[str, int], str] = {}
    ligand: list[list[float]] = []
    ligand_names: list[str] = []

    for line in path.read_text().splitlines():
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        element = line[76:78].strip()
        if element == "H":
            continue
        # Keep one alternate conformation. 6Q4G's ligand is modelled at two
        # altlocs; counting both double-counts its atoms and lets a conformer
        # the crystallographer ranked second contribute contacts.
        if line[16] not in (" ", "A"):
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        resname = line[17:20].strip()
        chain = line[21].strip() or "A"
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue

        if record == "HETATM" and resname not in _SOLVENT:
            ligand.append(xyz)
            ligand_names.append(resname)
        elif record == "ATOM  ":
            key = (chain, resseq)
            protein.setdefault(key, []).append(xyz)
            resnames[key] = resname

    packed = {k: (resnames[k], np.asarray(v, dtype=float)) for k, v in protein.items()}
    return packed, np.asarray(ligand, dtype=float), ligand_names


# Crystallisation additives and solvent, never counted as the ligand.
_SOLVENT = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "DMS", "ACT",
    "CL", "NA", "K", "MG", "CA", "ZN", "MN", "IOD", "BR", "NO3", "TRS",
}


def contacts_to_ligand(
    pdb: Path, *, chain: str = "A", cutoff: float = CUTOFF
) -> dict:
    """Residues of `chain` within `cutoff` of any ligand heavy atom."""
    protein, ligand, ligand_names = parse_pdb(pdb)
    if ligand.size == 0:
        return {"residues": [], "n_ligand_atoms": 0, "ligand_names": [], "reason": "no ligand atoms found"}

    engaged = {}
    for (ch, resseq), (resname, coords) in protein.items():
        if ch != chain:
            continue
        dmin = float(np.sqrt(((coords[:, None, :] - ligand[None, :, :]) ** 2).sum(-1)).min())
        if dmin <= cutoff:
            engaged[resseq] = {"resname": resname, "min_dist": round(dmin, 2)}

    return {
        "residues": sorted(engaged),
        "detail": {str(k): v for k, v in sorted(engaged.items())},
        "n_ligand_atoms": int(ligand.shape[0]),
        "ligand_names": sorted(set(ligand_names)),
        "cutoff": cutoff,
    }


def contacts_between_chains(
    pdb: Path, *, target: str = "A", binder: str = "B", cutoff: float = CUTOFF
) -> dict:
    """Target-side residues within `cutoff` of the binder chain.

    This is the design-side equivalent of contacts_to_ligand: it is what turns a
    BoltzGen design into target-side coordinates (PROJECT_GOAL.md 1.1).
    """
    protein, _, _ = parse_pdb(pdb)
    binder_coords = [c for (ch, _), (_, c) in protein.items() if ch == binder]
    if not binder_coords:
        return {"residues": [], "reason": f"no chain {binder} in {pdb.name}"}
    binder_xyz = np.vstack(binder_coords)

    engaged = {}
    for (ch, resseq), (resname, coords) in protein.items():
        if ch != target:
            continue
        dmin = float(np.sqrt(((coords[:, None, :] - binder_xyz[None, :, :]) ** 2).sum(-1)).min())
        if dmin <= cutoff:
            engaged[resseq] = {"resname": resname, "min_dist": round(dmin, 2)}

    return {
        "residues": sorted(engaged),
        "detail": {str(k): v for k, v in sorted(engaged.items())},
        "n_binder_atoms": int(binder_xyz.shape[0]),
        "cutoff": cutoff,
    }


def consensus(per_design: list[list[int]], *, core_threshold: float = 0.6) -> dict:
    """Aggregate per-design contact sets into an interface signature.

    Mirrors InterfaceSignature in PROJECT_GOAL.md 4.3: a per-residue engagement
    frequency, with `core` the residues engaged by at least `core_threshold` of
    the designs.
    """
    n = len(per_design)
    if n == 0:
        return {"n_designs": 0, "core": [], "frequency": {}, "mean_pairwise_jaccard": None}

    counts: dict[int, int] = {}
    for residues in per_design:
        for r in set(residues):
            counts[r] = counts.get(r, 0) + 1

    frequency = {r: c / n for r, c in counts.items()}
    core = sorted(r for r, f in frequency.items() if f >= core_threshold)

    # How much do the designs agree with each other? A consensus over designs
    # that disagree is not a signature of anything.
    jaccards = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = set(per_design[i]), set(per_design[j])
            if a or b:
                jaccards.append(len(a & b) / len(a | b))

    return {
        "n_designs": n,
        "core": core,
        "core_threshold": core_threshold,
        "frequency": {str(r): round(f, 3) for r, f in sorted(frequency.items())},
        "mean_pairwise_jaccard": round(float(np.mean(jaccards)), 4) if jaccards else None,
    }


def coverage(engaged: list[int], core: list[int]) -> float:
    """Fraction of core hotspot residues a molecule engages. The primary score."""
    if not core:
        return 0.0
    return len(set(engaged) & set(core)) / len(core)


def jaccard(a: list[int], b: list[int]) -> float:
    sa, sb = set(a), set(b)
    if not (sa or sb):
        return 0.0
    return len(sa & sb) / len(sa | sb)
