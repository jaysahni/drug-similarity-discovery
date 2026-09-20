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
        # Loud, not empty. An empty contact set silently becomes an empty
        # signature four steps later, which is how this was missed once already.
        present = sorted({ch for ch, _ in protein})
        raise ValueError(
            f"{pdb.name}: no polymer chain {binder!r}; chains present: {present}. "
            "Refusing to return an empty contact set."
        )
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


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def chain_sequence(pdb: Path, chain: str) -> tuple[list[int], str]:
    """Residue numbers of `chain`, in file order, and their one-letter sequence."""
    protein, _, _ = parse_pdb(pdb)
    items = sorted(
        ((resseq, resname) for (ch, resseq), (resname, _) in protein.items() if ch == chain)
    )
    if not items:
        raise ValueError(f"{pdb.name}: no polymer chain {chain!r}")
    numbers = [resseq for resseq, _ in items]
    sequence = "".join(THREE_TO_ONE.get(resname, "X") for _, resname in items)
    return numbers, sequence


def align_to_reference(pdb: Path, chain: str, reference: str) -> dict[int, int]:
    """Map this file's residue numbers to author numbering, by alignment.

    Needed because the numbering is not stable across this pipeline and no
    constant offset recovers it. `Protein.prepare` renumbers a structure
    contiguously from 1, discarding author numbers, and 6Q4G is missing 16 of
    CDK2's 298 residues, so the true mapping steps at every gap.

    Global alignment (Needleman-Wunsch, match +1 / mismatch -1 / gap -2) of the
    observed sequence against the canonical one. Refuses to return a mapping it
    cannot corroborate.

    Returns {resseq_in_file: author_number}.
    """
    numbers, observed = chain_sequence(pdb, chain)
    n, m = len(observed), len(reference)

    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    score[:, 0] = np.arange(0, -2 * (n + 1), -2)[: n + 1]
    score[0, :] = np.arange(0, -2 * (m + 1), -2)[: m + 1]
    for i in range(1, n + 1):
        row_prev, row = score[i - 1], score[i]
        for j in range(1, m + 1):
            diagonal = row_prev[j - 1] + (1 if observed[i - 1] == reference[j - 1] else -1)
            row[j] = max(diagonal, row_prev[j] - 2, row[j - 1] - 2)

    mapping: dict[int, int] = {}
    matches = 0
    i, j = n, m
    while i > 0 and j > 0:
        diagonal = score[i - 1][j - 1] + (1 if observed[i - 1] == reference[j - 1] else -1)
        if score[i][j] == diagonal:
            mapping[numbers[i - 1]] = j
            matches += observed[i - 1] == reference[j - 1]
            i, j = i - 1, j - 1
        elif score[i][j] == score[i - 1][j] - 2:
            i -= 1
        else:
            j -= 1

    identity = matches / len(observed)
    if identity < 0.95:
        raise ValueError(
            f"{pdb.name} chain {chain}: alignment to the reference sequence is only "
            f"{identity:.0%} identical over {len(observed)} residues. Refusing to map "
            "residues on an alignment this poor."
        )
    return mapping


def author_offset(pdb: Path, chain: str, reference: str) -> int:
    """Find the constant offset from this file's residue numbers to author numbering.

    Rowan hands back three different numberings across this pipeline, and guessing
    wrong produces a signature that is plausible and wrong. So it is measured:
    the offset is the one that makes the observed residue names agree with the
    canonical sequence, and it is only accepted if agreement is near total.

    Returns `offset` such that `author = resseq + offset`.
    """
    protein, _, _ = parse_pdb(pdb)
    observed = {
        resseq: THREE_TO_ONE.get(resname, "X")
        for (ch, resseq), (resname, _) in protein.items()
        if ch == chain
    }
    if not observed:
        raise ValueError(f"{pdb.name}: no polymer chain {chain!r}")

    best, best_hits = None, -1
    span = max(observed) - min(observed) + len(reference)
    for offset in range(-span, span + 1):
        hits = sum(
            1
            for resseq, aa in observed.items()
            if 1 <= resseq + offset <= len(reference)
            and reference[resseq + offset - 1] == aa
        )
        if hits > best_hits:
            best, best_hits = offset, hits

    fraction = best_hits / len(observed)
    if fraction < 0.95:
        raise ValueError(
            f"{pdb.name} chain {chain}: best offset {best:+d} matches only "
            f"{best_hits}/{len(observed)} residues ({fraction:.0%}) against the "
            "reference sequence. Refusing to map residues on a guess."
        )
    return best


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
