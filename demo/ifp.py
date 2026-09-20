"""Pose-based interaction fingerprints, and the head-to-head against ECFP4.

WHY THIS EXISTS
---------------
The small-molecule arm ranks approved drugs by ECFP4 Tanimoto on SMILES. That
measures shared substructure, and substructure is not what a target sees. The
failure is visible in results/demo/thrombin/09_smallmol_match.json: argatroban's
nearest approved drugs come out as peptidomimetics (angiotensin II, icatibant,
lisinopril) because argatroban *looks* like a peptide, and not one of them is a
thrombin drug even though argatroban is one. ECFP4 groups by chemotype.

A docked pose carries the thing ECFP4 cannot see: which residues of the target
the molecule actually touches, and how. An interaction fingerprint (IFP) turns a
pose into a fixed-length bit vector over (site residue x interaction type), so
two molecules that engage Asp189 and the 60-loop the same way score as similar
even when their scaffolds share nothing. This is the Deng/Chuaqui/Singh (J Med
Chem 2004, 47:337) SIFt idea with the interaction typing of IChem/ProLIF.

WHAT IS AND IS NOT MEASURED HERE
--------------------------------
Six channels per residue are implemented and each is a real geometric test, not
a relabelled distance: hydrophobic, ligand-as-H-bond-donor, ligand-as-acceptor,
ionic (ligand cation / ligand anion), and aromatic (pi-stacking). Ligand atom
typing uses RDKit's BaseFeatures.fdef pharmacophore definitions; protein atom
typing uses the residue/atom-name tables below.

The honest caveat: no structure in this pipeline carries explicit hydrogens --
1PPB is a 1.9 A X-ray structure deposited without them, and Vina poses come back
as heavy atoms. So the H-bond channels apply a *heavy-atom* geometric criterion
(donor-acceptor distance plus a >=90 deg angle at the donor and at the acceptor
against their own covalent neighbours), not a true D-H...A angle. Donor and
acceptor roles are assigned by atom typing, not observed from H positions.
Protonation states are likewise assigned, not measured: carboxylates are treated
as anionic and guanidinium/amidinium/primary amines as cationic, but histidine
is left neutral because its state at the thrombin active site is not known here.

SIMILARITY
----------
Tanimoto is reported because it is what the rest of this repo uses, but it
ignores the bits that are off in both fingerprints -- and for an IFP those bits
are meaningful, because the bit space is a fixed, interpretable site rather than
a hash. Racz, Bajusz & Heberger (J Cheminform 2018, 10:48, doi:10.1186/
s13321-018-0302-y) benchmarked 44 measures on binary fingerprints and found six
that outperform Tanimoto; all six are implemented here alongside it. Cosine is
deliberately absent: it was not among them.

Run:  ./env-kit/bin/python -m demo.ifp --check-poses
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import AllChem, ChemicalFeatures

from demo.contacts import residue_label, residue_order

REPO = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# geometry thresholds
#
# Taken from the IFP literature (Marcou & Rognan, J Chem Inf Model 2007, 47:195;
# ProLIF, Bouysset & Rognan, J Cheminform 2021, 13:72) rather than fitted here.
# Nothing downstream tunes them.
# ---------------------------------------------------------------------------

HYDROPHOBIC_CUTOFF = 4.5   # A, apolar heavy atom to apolar heavy atom
HBOND_CUTOFF = 3.5         # A, donor heavy atom to acceptor heavy atom
HBOND_MIN_ANGLE = 90.0     # deg, at donor and at acceptor, against covalent neighbours
IONIC_CUTOFF = 4.5         # A, charged group centroid to charged group centroid
# Pi stacking needs two cutoffs, not one: in an edge-to-face (T-shaped) pair the
# centroids sit further apart than in a stacked pair by construction, so a single
# distance would rule out exactly the geometry thrombin's aryl-binding site uses.
# These are ProLIF's published PiStacking defaults.
AROMATIC_FACE_CUTOFF = 5.5  # A, centroid to centroid, face-to-face
AROMATIC_EDGE_CUTOFF = 6.5  # A, centroid to centroid, edge-to-face
AROMATIC_FACE_MAX = 35.0    # deg between ring normals, face-to-face stacking
AROMATIC_EDGE_MIN = 50.0    # deg between ring normals, edge-to-face stacking

CONTACT_TYPES = (
    "hydrophobic",
    "hbond_donor",     # the LIGAND donates
    "hbond_acceptor",  # the LIGAND accepts
    "ionic_cation",    # the LIGAND is the cation
    "ionic_anion",     # the LIGAND is the anion
    "aromatic",
)

# ---------------------------------------------------------------------------
# protein-side atom typing
# ---------------------------------------------------------------------------

_BACKBONE_DONOR = {"N"}          # every residue but proline
_BACKBONE_ACCEPTOR = {"O", "OXT"}

_SIDECHAIN_DONOR = {
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
    "HIS": {"ND1", "NE2"}, "LYS": {"NZ"}, "SER": {"OG"}, "THR": {"OG1"},
    "TRP": {"NE1"}, "TYR": {"OH"}, "CYS": {"SG"},
}
_SIDECHAIN_ACCEPTOR = {
    "ASN": {"OD1"}, "ASP": {"OD1", "OD2"}, "GLN": {"OE1"},
    "GLU": {"OE1", "OE2"}, "HIS": {"ND1", "NE2"}, "SER": {"OG"},
    "THR": {"OG1"}, "TYR": {"OH"}, "MET": {"SD"},
}
# Charged groups, as groups: the fingerprint should not fire twice for the two
# oxygens of one carboxylate, so ionic contacts are measured centroid to
# centroid. HIS is absent on purpose -- see the module docstring.
_ANIONIC_GROUPS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_CATIONIC_GROUPS = {"ARG": ("NE", "CZ", "NH1", "NH2"), "LYS": ("NZ",)}

_AROMATIC_RINGS = {
    "PHE": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "TYR": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "HIS": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "TRP": [("CG", "CD1", "NE1", "CE2", "CD2"),
            ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")],
}

# Carbons bonded to N or O are not apolar, and without protein connectivity the
# only correct way to say that is per residue, by name. TYR CZ carries the
# hydroxyl and is excluded; PHE CZ does not and is kept.
_APOLAR_SIDECHAIN = {
    "ALA": {"CB"},
    "ARG": {"CB", "CG"},
    "ASN": {"CB"},
    "ASP": {"CB"},
    "CYS": {"CB", "SG"},
    "GLN": {"CB", "CG"},
    "GLU": {"CB", "CG"},
    "GLY": set(),
    "HIS": {"CB"},
    "ILE": {"CB", "CG1", "CG2", "CD1"},
    "LEU": {"CB", "CG", "CD1", "CD2"},
    "LYS": {"CB", "CG", "CD"},
    "MET": {"CB", "CG", "SD", "CE"},
    "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"CB", "CG", "CD"},
    "SER": set(),
    "THR": {"CG2"},
    "TRP": {"CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2"},
    "VAL": {"CB", "CG1", "CG2"},
}

_SOLVENT = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "DMS", "ACT",
    "CL", "NA", "K", "MG", "CA", "ZN", "MN", "IOD", "BR", "NO3", "TRS",
}

_COVALENT_MAX = 1.95  # A; used only to find an atom's own covalent neighbours


# ---------------------------------------------------------------------------
# receptor
# ---------------------------------------------------------------------------


@dataclass
class Receptor:
    """A parsed protein, keyed by (chain, residue label).

    demo.contacts.parse_pdb deliberately throws atom names away -- a distance
    cutoff does not need them. Interaction typing does, so this is a second,
    richer parse that follows exactly the same conventions (heavy atoms only,
    one altloc, insertion-code-aware labels via contacts.residue_label) so the
    two never disagree about which residue an atom belongs to.
    """

    residues: dict[tuple[str, int | str], dict]

    def chains(self) -> list[str]:
        return sorted({ch for ch, _ in self.residues})

    def labels(self, chain: str) -> list[int | str]:
        return sorted((lab for ch, lab in self.residues if ch == chain), key=residue_order)

    def require(self, chain: str, labels: list[int | str]) -> None:
        """Fail loudly if the site is not actually in the structure."""
        if chain not in self.chains():
            raise ValueError(
                f"chain {chain!r} is not in this receptor; chains present: {self.chains()}"
            )
        missing = [lab for lab in labels if (chain, lab) not in self.residues]
        if missing:
            raise ValueError(
                f"{len(missing)} of {len(labels)} site residues are absent from chain "
                f"{chain}: {missing[:12]}. Refusing to build a fingerprint over an axis "
                "the structure does not contain -- every missing residue would be a bit "
                "that is silently zero for every molecule."
            )


def parse_receptor(path: Path) -> Receptor:
    """Read a PDB into per-residue typed atom tables."""
    residues: dict[tuple[str, int | str], dict] = {}
    for line in Path(path).read_text().splitlines():
        if line[:6] != "ATOM  ":
            continue
        if line[76:78].strip() == "H":
            continue
        if line[16] not in (" ", "A"):
            continue
        resname = line[17:20].strip()
        if resname in _SOLVENT:
            continue
        chain = line[21].strip() or "A"
        try:
            label = residue_label(int(line[22:26]), line[26])
        except ValueError:
            continue
        key = (chain, label)
        entry = residues.setdefault(key, {"resname": resname, "names": [], "xyz": []})
        entry["names"].append(line[12:16].strip())
        entry["xyz"].append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )

    if not residues:
        raise ValueError(f"{Path(path).name}: no ATOM records parsed. Refusing to "
                         "return an empty receptor.")
    for entry in residues.values():
        entry["xyz"] = np.asarray(entry["xyz"], dtype=float)
    return Receptor(residues=residues)


def _residue_features(entry: dict) -> dict:
    """Typed atom/group coordinates for one residue."""
    resname, names, xyz = entry["resname"], entry["names"], entry["xyz"]
    index = {n: i for i, n in enumerate(names)}

    donors, acceptors, apolar = [], [], []
    for i, name in enumerate(names):
        element = name[0] if not name[0].isdigit() else name[1]
        if name in _BACKBONE_DONOR and resname != "PRO":
            donors.append(i)
        if name in _SIDECHAIN_DONOR.get(resname, ()):
            donors.append(i)
        if name in _BACKBONE_ACCEPTOR:
            acceptors.append(i)
        if name in _SIDECHAIN_ACCEPTOR.get(resname, ()):
            acceptors.append(i)
        if name in _APOLAR_SIDECHAIN.get(resname, set()) and element in ("C", "S"):
            apolar.append(i)
        elif name == "CA":
            apolar.append(i)  # the alpha carbon lines every pocket wall

    def group(names_wanted) -> np.ndarray | None:
        got = [index[n] for n in names_wanted if n in index]
        if len(got) != len(names_wanted):
            return None  # an incomplete side chain is not a charged group
        return xyz[got].mean(axis=0)

    anion = group(_ANIONIC_GROUPS[resname]) if resname in _ANIONIC_GROUPS else None
    cation = group(_CATIONIC_GROUPS[resname]) if resname in _CATIONIC_GROUPS else None

    rings = []
    for ring_names in _AROMATIC_RINGS.get(resname, ()):
        got = [index[n] for n in ring_names if n in index]
        if len(got) == len(ring_names):
            rings.append(_ring_frame(xyz[got]))

    return {
        "xyz": xyz,
        "donors": np.asarray(sorted(set(donors)), dtype=int),
        "acceptors": np.asarray(sorted(set(acceptors)), dtype=int),
        "apolar": np.asarray(sorted(set(apolar)), dtype=int),
        "anion": anion,
        "cation": cation,
        "rings": rings,
    }


def _ring_frame(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(centroid, unit normal) of a ring, the normal by SVD of the centred ring."""
    centroid = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - centroid)
    return centroid, vt[2] / np.linalg.norm(vt[2])


def _neighbours(xyz: np.ndarray, i: int) -> np.ndarray:
    """Indices of atom i's covalent neighbours within the same residue."""
    d = np.linalg.norm(xyz - xyz[i], axis=1)
    return np.where((d > 1e-6) & (d <= _COVALENT_MAX))[0]


def _angle_ok(root_xyz: np.ndarray, vertex: np.ndarray, far: np.ndarray) -> bool:
    """True if every root-vertex...far angle clears HBOND_MIN_ANGLE.

    With no explicit hydrogens this is the strongest geometric statement the
    coordinates support: an acceptor sitting directly behind the donor's own
    covalent bonds cannot be hydrogen bonded to it.
    """
    if len(root_xyz) == 0:
        return True
    v1 = root_xyz - vertex
    v2 = far - vertex
    cos = (v1 @ v2) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2) + 1e-12)
    angles = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    return bool(angles.min() >= HBOND_MIN_ANGLE)


# ---------------------------------------------------------------------------
# ligand pose
# ---------------------------------------------------------------------------

RDLogger.DisableLog("rdApp.warning")
_FACTORY = ChemicalFeatures.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)


@dataclass
class Pose:
    """A ligand conformation with its chemistry perceived."""

    mol: Chem.Mol
    name: str

    @property
    def xyz(self) -> np.ndarray:
        return self.mol.GetConformer().GetPositions()


def pose_from_pdb_block(block: str, template_smiles: str, name: str) -> Pose:
    """Build a pose from PDB HETATM lines plus the SMILES that says what it is.

    Proximity bonding alone gives every bond order 1 and no aromaticity or
    formal charge, which would silently disable four of the six channels. The
    template is therefore required, not optional.
    """
    mol = Chem.MolFromPDBBlock(block, proximityBonding=True, removeHs=False,
                               sanitize=False)
    if mol is None:
        raise ValueError(f"{name}: RDKit could not read the PDB block")
    Chem.SanitizeMol(mol)
    template = Chem.MolFromSmiles(template_smiles)
    if template is None:
        raise ValueError(f"{name}: template SMILES {template_smiles!r} did not parse")
    if template.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError(
            f"{name}: template has {template.GetNumAtoms()} heavy atoms but the pose "
            f"has {mol.GetNumAtoms()}. Refusing to type a pose against the wrong molecule."
        )
    return Pose(mol=AllChem.AssignBondOrdersFromTemplate(template, mol), name=name)


def poses_from_sdf(path: Path) -> list[Pose]:
    """Read docked poses from an SDF, which already carries bond orders."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True)
    poses = []
    for i, mol in enumerate(supplier):
        if mol is None:
            raise ValueError(f"{Path(path).name}: molecule {i} failed to parse. "
                             "Refusing to skip it silently.")
        if mol.GetNumConformers() == 0:
            raise ValueError(f"{Path(path).name}: molecule {i} has no 3D conformer")
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"pose_{i}"
        poses.append(Pose(mol=mol, name=name))
    if not poses:
        raise ValueError(f"{Path(path).name}: no poses read")
    return poses


def subpose(pose: Pose, atom_ids: list[int], name: str) -> Pose:
    """A pose restricted to some of its atoms, chemistry re-perceived.

    Used to ask what one moiety of a ligand touches -- e.g. whether the P1
    arginine of PPACK and its P3 phenyl engage different subsites of thrombin.
    """
    editable = Chem.RWMol(pose.mol)
    for i in sorted(set(range(pose.mol.GetNumAtoms())) - set(atom_ids), reverse=True):
        editable.RemoveAtom(i)
    fragment = editable.GetMol()
    for atom in fragment.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    Chem.SanitizeMol(fragment)
    return Pose(mol=fragment, name=name)


def _ligand_features(pose: Pose) -> dict:
    """Typed ligand atoms and groups, from RDKit's BaseFeatures pharmacophores."""
    mol, xyz = pose.mol, pose.xyz
    donors, acceptors, cations, anions, rings = [], [], [], [], []

    for feature in _FACTORY.GetFeaturesForMol(mol):
        family, ids = feature.GetFamily(), list(feature.GetAtomIds())
        if family == "Donor":
            # BaseFeatures marks some amide nitrogens that carry no hydrogen.
            # A donor with nothing to donate is not a donor.
            donors += [i for i in ids if mol.GetAtomWithIdx(i).GetTotalNumHs() > 0]
        elif family == "Acceptor":
            acceptors += ids
        elif family == "PosIonizable" and len(ids) > 1:
            cations.append(xyz[ids].mean(axis=0))
        elif family == "PosIonizable":
            cations.append(xyz[ids[0]])
        elif family == "NegIonizable":
            anions.append(xyz[ids].mean(axis=0) if len(ids) > 1 else xyz[ids[0]])
        elif family == "Aromatic":
            rings.append(_ring_frame(xyz[ids]))

    # Hydrophobic contact uses the IFP-literature definition (an apolar heavy
    # atom: C, S or halogen with no N/O neighbour), not RDKit's Hydrophobe
    # family, which is far narrower and would miss most of a docked scaffold.
    apolar = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() in (6, 16, 9, 17, 35, 53)
        and not any(n.GetAtomicNum() in (7, 8) for n in a.GetNeighbors())
    ]

    return {
        "xyz": xyz,
        "donors": sorted(set(donors)),
        "acceptors": sorted(set(acceptors)),
        "apolar": apolar,
        "cations": cations,
        "anions": anions,
        "rings": rings,
    }


# ---------------------------------------------------------------------------
# the fingerprint
# ---------------------------------------------------------------------------


@dataclass
class SiteAxis:
    """The fixed bit space: site residues x interaction types.

    Unlike an ECFP hash, every bit here names something -- "residue 189, ionic"
    -- which is why the both-off bits carry information and why the similarity
    measures below that use them are defensible.
    """

    chain: str
    residues: list[int | str]
    types: tuple[str, ...] = CONTACT_TYPES

    def __post_init__(self) -> None:
        self.residues = sorted(set(self.residues), key=residue_order)

    def __len__(self) -> int:
        return len(self.residues) * len(self.types)

    @property
    def names(self) -> list[str]:
        return [f"{r}:{t}" for r in self.residues for t in self.types]

    def index(self, residue: int | str, type_: str) -> int:
        return self.residues.index(residue) * len(self.types) + self.types.index(type_)


@dataclass
class Fingerprint:
    bits: np.ndarray
    axis: SiteAxis
    name: str

    def on(self) -> list[str]:
        return [n for n, b in zip(self.axis.names, self.bits) if b]

    def by_type(self) -> dict[str, list[int | str]]:
        out: dict[str, list[int | str]] = {t: [] for t in self.axis.types}
        for i, b in enumerate(self.bits):
            if b:
                out[self.axis.types[i % len(self.axis.types)]].append(
                    self.axis.residues[i // len(self.axis.types)]
                )
        return out


def fingerprint(pose: Pose, receptor: Receptor, axis: SiteAxis) -> Fingerprint:
    """Typed interaction fingerprint of one pose against one site."""
    receptor.require(axis.chain, axis.residues)
    lig = _ligand_features(pose)
    bits = np.zeros(len(axis), dtype=bool)
    n_types = len(axis.types)

    for r, label in enumerate(axis.residues):
        res = _residue_features(receptor.residues[(axis.chain, label)])
        base = r * n_types

        def set_bit(type_: str) -> None:
            if type_ in axis.types:
                bits[base + axis.types.index(type_)] = True

        # hydrophobic
        if len(res["apolar"]) and lig["apolar"]:
            d = _cdist(res["xyz"][res["apolar"]], lig["xyz"][lig["apolar"]])
            if d.min() <= HYDROPHOBIC_CUTOFF:
                set_bit("hydrophobic")

        # H bonds, both directions, with the heavy-atom angle filter
        if len(res["acceptors"]) and lig["donors"]:
            if _hbond(res["xyz"], res["acceptors"], lig["xyz"], lig["donors"]):
                set_bit("hbond_donor")
        if len(res["donors"]) and lig["acceptors"]:
            if _hbond(res["xyz"], res["donors"], lig["xyz"], lig["acceptors"]):
                set_bit("hbond_acceptor")

        # ionic, centroid to centroid so one carboxylate fires one bit
        if res["anion"] is not None and lig["cations"]:
            if min(np.linalg.norm(c - res["anion"]) for c in lig["cations"]) <= IONIC_CUTOFF:
                set_bit("ionic_cation")
        if res["cation"] is not None and lig["anions"]:
            if min(np.linalg.norm(c - res["cation"]) for c in lig["anions"]) <= IONIC_CUTOFF:
                set_bit("ionic_anion")

        # aromatic: centroid distance and an interplanar angle that is either
        # face-to-face or edge-to-face. A ring pair at 45 deg is neither.
        for p_centroid, p_normal in res["rings"]:
            for l_centroid, l_normal in lig["rings"]:
                d = float(np.linalg.norm(p_centroid - l_centroid))
                if d > AROMATIC_EDGE_CUTOFF:
                    continue
                angle = math.degrees(
                    math.acos(min(1.0, abs(float(p_normal @ l_normal))))
                )
                face = angle <= AROMATIC_FACE_MAX and d <= AROMATIC_FACE_CUTOFF
                edge = angle >= AROMATIC_EDGE_MIN and d <= AROMATIC_EDGE_CUTOFF
                if face or edge:
                    set_bit("aromatic")

    return Fingerprint(bits=bits, axis=axis, name=pose.name)


def _cdist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))


def _hbond(p_xyz, p_idx, l_xyz, l_idx) -> bool:
    """Is there a heavy-atom hydrogen bond between these two atom sets?"""
    d = _cdist(p_xyz[p_idx], l_xyz[l_idx])
    for pi, li in zip(*np.where(d <= HBOND_CUTOFF)):
        p_atom, l_atom = p_xyz[p_idx[pi]], l_xyz[l_idx[li]]
        p_roots = p_xyz[_neighbours(p_xyz, int(p_idx[pi]))]
        l_roots = l_xyz[_neighbours(l_xyz, int(l_idx[li]))]
        if _angle_ok(p_roots, p_atom, l_atom) and _angle_ok(l_roots, l_atom, p_atom):
            return True
    return False


# ---------------------------------------------------------------------------
# similarity
#
# a = on in both, b = on in A only, c = on in B only, d = off in both, n = total.
# Racz, Bajusz & Heberger 2018 (doi:10.1186/s13321-018-0302-y) benchmarked 44
# measures and reported six that beat Tanimoto; those six are the six below it.
# Every one of them uses d, which is exactly why they need a bit space where
# "off in both" means something -- an IFP over a fixed site is such a space.
# Formulae follow Todeschini et al., J Chem Inf Model 2012, 52:2884.
# ---------------------------------------------------------------------------


def counts(a_bits: np.ndarray, b_bits: np.ndarray) -> tuple[int, int, int, int]:
    if a_bits.shape != b_bits.shape:
        raise ValueError(
            f"fingerprints have different lengths ({a_bits.shape} vs {b_bits.shape}); "
            "they are not on the same site axis and must not be compared."
        )
    a = int(np.sum(a_bits & b_bits))
    b = int(np.sum(a_bits & ~b_bits))
    c = int(np.sum(~a_bits & b_bits))
    return a, b, c, int(a_bits.size) - a - b - c


def tanimoto(x, y) -> float:
    a, b, c, _ = counts(x, y)
    return 1.0 if a + b + c == 0 else a / (a + b + c)


def sokal_michener(x, y) -> float:
    a, b, c, d = counts(x, y)
    return (a + d) / (a + b + c + d)


def rogers_tanimoto(x, y) -> float:
    a, b, c, d = counts(x, y)
    return (a + d) / (a + 2 * (b + c) + d)


def sokal_sneath_2(x, y) -> float:
    a, b, c, d = counts(x, y)
    return 2 * (a + d) / (2 * (a + d) + b + c)


def consonni_todeschini_1(x, y) -> float:
    a, b, c, d = counts(x, y)
    return math.log1p(a + d) / math.log1p(a + b + c + d)


def consonni_todeschini_2(x, y) -> float:
    a, b, c, d = counts(x, y)
    n = a + b + c + d
    return (math.log1p(n) - math.log1p(b + c)) / math.log1p(n)


def austin_colwell(x, y) -> float:
    a, b, c, d = counts(x, y)
    return (2 / math.pi) * math.asin(math.sqrt((a + d) / (a + b + c + d)))


MEASURES = {
    "tanimoto": tanimoto,
    "sokal_michener": sokal_michener,
    "rogers_tanimoto": rogers_tanimoto,
    "sokal_sneath_2": sokal_sneath_2,
    "consonni_todeschini_1": consonni_todeschini_1,
    "consonni_todeschini_2": consonni_todeschini_2,
    "austin_colwell": austin_colwell,
}
BEATS_TANIMOTO = tuple(k for k in MEASURES if k != "tanimoto")


def similarity(a: Fingerprint, b: Fingerprint, measure: str = "tanimoto") -> float:
    if measure not in MEASURES:
        raise KeyError(f"unknown measure {measure!r}; have {sorted(MEASURES)}")
    return MEASURES[measure](a.bits, b.bits)


def all_similarities(a: Fingerprint, b: Fingerprint) -> dict[str, float]:
    return {name: fn(a.bits, b.bits) for name, fn in MEASURES.items()}


# ---------------------------------------------------------------------------
# the head-to-head: does IFP rank known thrombin binders better than ECFP4?
# ---------------------------------------------------------------------------


def ecfp4(smiles: str):
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES did not parse: {smiles!r}")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return gen.GetFingerprint(mol)


def ecfp4_similarity(a: str, b: str) -> float:
    from rdkit import DataStructs

    return float(DataStructs.TanimotoSimilarity(ecfp4(a), ecfp4(b)))


def auc(scores: np.ndarray, actives: np.ndarray) -> float:
    """Probability a random active outranks a random inactive (ties count half).

    This is Mann-Whitney U normalised, i.e. the ROC AUC, computed directly so
    that it is auditable at n = 2 actives where a curve would be meaningless.
    """
    scores, actives = np.asarray(scores, float), np.asarray(actives, bool)
    pos, neg = scores[actives], scores[~actives]
    if pos.size == 0 or neg.size == 0:
        raise ValueError(
            f"AUC needs both classes; got {pos.size} actives and {neg.size} inactives."
        )
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def permutation_p(scores, actives, *, n_perm: int = 20000, seed: int = 0) -> float:
    """One-sided p for AUC > 0.5, by permuting the active labels.

    Exact-in-the-limit and valid at any number of actives, unlike the normal
    approximation, which is why it is used here: the honest active set is small.
    """
    scores, actives = np.asarray(scores, float), np.asarray(actives, bool)
    observed = auc(scores, actives)
    rng = np.random.default_rng(seed)
    labels = actives.copy()
    hits = sum(auc(scores, rng.permutation(labels)) >= observed for _ in range(n_perm))
    return (hits + 1) / (n_perm + 1)


def paired_bootstrap_delta_auc(
    scores_a, scores_b, actives, *, n_boot: int = 20000, seed: int = 0
) -> dict:
    """AUC(a) - AUC(b) with a percentile CI, resampling molecules, not scores.

    The two rankings are computed on the same molecules, so their AUCs are
    correlated and an unpaired test would be wrong. Actives and inactives are
    resampled within their own strata so every bootstrap replicate still has
    both classes.
    """
    scores_a = np.asarray(scores_a, float)
    scores_b = np.asarray(scores_b, float)
    actives = np.asarray(actives, bool)
    observed = auc(scores_a, actives) - auc(scores_b, actives)

    pos_idx = np.where(actives)[0]
    neg_idx = np.where(~actives)[0]
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot)
    for k in range(n_boot):
        idx = np.concatenate([
            rng.choice(pos_idx, pos_idx.size, replace=True),
            rng.choice(neg_idx, neg_idx.size, replace=True),
        ])
        lab = np.zeros(idx.size, dtype=bool)
        lab[: pos_idx.size] = True
        deltas[k] = auc(scores_a[idx], lab) - auc(scores_b[idx], lab)
    return {
        "delta_auc": observed,
        "ci95": (float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))),
        # Floored at 1/n_boot: a bootstrap cannot resolve a p-value smaller than
        # its own resolution, and reporting 0.0 would claim it could.
        "p_two_sided": max(
            float(2 * min((deltas <= 0).mean(), (deltas >= 0).mean())), 1 / n_boot
        ),
        "p_two_sided_is_floor": bool(
            2 * min((deltas <= 0).mean(), (deltas >= 0).mean()) < 1 / n_boot
        ),
        "n_actives": int(pos_idx.size),
        "n_inactives": int(neg_idx.size),
        "n_boot": n_boot,
    }


def head_to_head(
    ifp_scores, ecfp_scores, actives, names, *, seed: int = 0
) -> dict:
    """Rank the same molecules two ways and compare where the actives land."""
    actives = np.asarray(actives, bool)
    if actives.sum() == 0:
        raise ValueError("no actives labelled; there is nothing to compare rankings on.")
    out = {
        "n": int(actives.size),
        "n_actives": int(actives.sum()),
        "auc_ifp": auc(ifp_scores, actives),
        "auc_ecfp4": auc(ecfp_scores, actives),
        "p_ifp_vs_chance": permutation_p(ifp_scores, actives, seed=seed),
        "p_ecfp4_vs_chance": permutation_p(ecfp_scores, actives, seed=seed),
        "delta": paired_bootstrap_delta_auc(ifp_scores, ecfp_scores, actives, seed=seed),
        "ranks": {},
    }
    r_ifp = _ranks(ifp_scores)
    r_ecfp = _ranks(ecfp_scores)
    for i in np.where(actives)[0]:
        out["ranks"][names[i]] = {"ifp": int(r_ifp[i]), "ecfp4": int(r_ecfp[i])}
    return out


def _ranks(scores) -> np.ndarray:
    """1 = most similar. Ties take the average rank."""
    from scipy.stats import rankdata

    return rankdata(-np.asarray(scores, float), method="average")


# ---------------------------------------------------------------------------
# poses from Rowan
# ---------------------------------------------------------------------------

DOCKING_UUID = "4b5a8570-3086-4f05-a130-8aa56416aa92"


def load_env_if_available() -> None:
    """Read .env into the environment without demanding that a key be there.

    demo.nova.load_env is the only credential path this demo uses, deliberately:
    novakit otherwise falls back to the OS keyring, which blocks on a GUI prompt
    in a non-interactive shell. This wrapper exists so a test can ask whether
    credentials are present and skip cleanly, instead of exiting.
    """
    from demo import nova

    nova.load_env()


def check_saved_poses(workflow_uuid: str = DOCKING_UUID) -> dict:
    """Did a batch-docking workflow persist its poses? Read-only; costs nothing.

    Answering this from the workflow record rather than assuming it is the point:
    a batch_docking workflow that saved nothing looks exactly like one that saved
    everything until you read num_poses_to_save.
    """
    from demo import nova

    rowan = nova.rowan_client()
    record = rowan.retrieve_workflow(workflow_uuid).model_dump()
    data = record.get("data", {})
    n_saved = data.get("num_poses_to_save")
    descendants = rowan.list_proteins(ancestor_uuid=workflow_uuid, size=50)

    return {
        "workflow_uuid": workflow_uuid,
        "workflow_type": record.get("workflow_type"),
        "num_poses_to_save": n_saved,
        "max_poses_per_ligand": data.get("docking_settings", {}).get("max_poses"),
        "n_descendant_structures": len(descendants),
        "n_ligands": len(data.get("initial_smiles_list", [])),
        "credits_charged": record.get("credits_charged"),
        "poses_available": bool(n_saved) and len(descendants) > 0,
    }


def redock_with_poses(
    smiles_list: list[str],
    *,
    n_poses: int = 1,
    name: str = "F2 arm S batch docking (poses saved)",
    cache_key: str = "dock_thrombin_poses",
    approve: bool = False,
    max_credits: int = 30,
) -> dict:
    """Dock a ligand list with pose saving switched on, memoised by cache_key.

    rowan.submit_batch_docking_workflow does not expose num_poses_to_save -- it
    is a field of stjames.BatchDockingWorkflow that defaults to 0 -- so the
    workflow payload has to be built by hand and posted through submit_workflow.
    That default is why the first run of this batch kept nothing: see
    docs/14-INTERACTION-FINGERPRINTS.md.

    Caching follows demo/nova.py: the uuid goes to results/demo/_cache/, so a
    re-run re-reads the finished workflow instead of paying for it twice.
    Spends credits, and will not submit without approve=True.
    """
    import json

    from demo import nova

    hit = nova.cached(cache_key)
    if hit and hit.get("uuid"):
        return {"submitted": False, "reason": "cached", "uuid": hit["uuid"],
                "n_ligands": hit.get("n_ligands")}

    site = json.loads((REPO / "results/demo/thrombin/02_site.json").read_text())
    docked = json.loads((REPO / "results/demo/thrombin/08_docked.json").read_text())
    payload = {
        "initial_smiles_list": smiles_list,
        "target": site["protein_uuid"],
        "protein": site["protein_uuid"],
        # The same box as the original run, so the new poses sit in the same
        # frame as the scores already on record.
        "pocket": docked["box"],
        "docking_settings": {
            "executable": "vina", "scoring_function": "vinardo", "exhaustiveness": 8,
        },
        "num_poses_to_save": n_poses,
    }
    # The first run charged 19.23 credits for 203 ligands.
    estimate = 19.23 * len(smiles_list) / 203
    if not approve:
        return {"submitted": False, "reason": "approve=False",
                "estimated_credits": round(estimate, 2),
                "n_ligands": len(smiles_list), "payload": payload}
    if estimate > max_credits:
        raise SystemExit(
            f"estimated {estimate:.1f} credits for {len(smiles_list)} ligands, over the "
            f"{max_credits} ceiling. Refusing to submit."
        )

    # rowan.submit_workflow insists on an initial_molecule, which a batch has no
    # single one of ("You must provide either `initial_smiles` or a valid
    # `initial_molecule`"), so the request goes the same way rowan's own
    # submit_batch_docking_workflow sends it: a direct POST of the stjames
    # workflow model. The only difference from that helper is the one field it
    # does not expose.
    import stjames

    rowan = nova.rowan_client()
    model = stjames.BatchDockingWorkflow(**payload)
    if model.num_poses_to_save != n_poses:
        raise RuntimeError(
            f"stjames dropped num_poses_to_save ({model.num_poses_to_save} != "
            f"{n_poses}). Submitting would repeat the run that kept nothing."
        )
    body = {
        "workflow_type": "batch_docking",
        "workflow_data": model.model_dump(serialize_as_any=True, mode="json"),
        "name": name,
        "max_credits": max_credits,
    }
    with rowan.api_client() as client:
        response = client.post("/workflow", json=body)
        response.raise_for_status()
        uuid = str(response.json()["uuid"])

    nova.cache(cache_key, {"uuid": uuid, "label": name,
                           "submitter": "POST /workflow batch_docking",
                           "n_ligands": len(smiles_list), "num_poses_to_save": n_poses})
    return {"submitted": True, "uuid": uuid,
            "estimated_credits": round(estimate, 2), "n_ligands": len(smiles_list)}


def thrombin_ligand_set() -> dict:
    """The enlarged docking set: 200 generated decoys plus measured F2 binders.

    The active set is NOT hand-picked. It is every molecule in
    data/target_annotations.csv with a quantitative ChEMBL affinity against
    P00734 (thrombin) that clears this repo's own `passes_threshold` flag --
    a criterion that already existed before this arm did, so it cannot have been
    chosen to make a result come out. SMILES come from data/approved_drugs.csv
    by struct_id, never typed from memory: a mistyped SMILES is a different
    molecule and would not announce itself.

    Molecules with a measured but sub-threshold affinity (apixaban, edoxaban --
    factor Xa drugs with weak thrombin cross-reactivity) are docked too and fall
    in the inactive class by that same pre-existing rule. Their ranks are worth
    reading separately, but the labels are not adjusted after the fact.

    Bivalirudin is the one mechanism-of-action thrombin drug left out: it is a
    20-residue peptide of 155 heavy atoms and does not belong in a 30 x 26 x 22 A
    docking box. That exclusion is by size, decided before any score was seen.
    """
    import csv
    import json

    annotations = [
        row for row in csv.DictReader(
            (REPO / "data/target_annotations.csv").open()
        )
        if row["uniprot"] == "P00734" and row["evidence"] == "chembl_pchembl"
    ]
    drugs = {
        row["struct_id"]: row
        for row in csv.DictReader((REPO / "data/approved_drugs.csv").open())
    }
    docked = json.loads((REPO / "results/demo/thrombin/08_docked.json").read_text())

    n_generated = docked["n_generated_docked"]
    generated = docked["smiles"][:n_generated]
    if len(generated) != n_generated:
        raise ValueError("08_docked.json: fewer SMILES than n_generated_docked")

    records = []
    for row in sorted(annotations, key=lambda r: -float(r["pchembl"])):
        drug = drugs.get(row["struct_id"])
        if drug is None:
            raise ValueError(
                f"{row['name']}: struct_id {row['struct_id']} has a thrombin affinity "
                "but no SMILES in approved_drugs.csv. Refusing to guess one."
            )
        records.append({
            "name": row["name"],
            "smiles": drug["smiles"],
            "pchembl": float(row["pchembl"]),
            "activity": f"{row['activity_type']}={row['activity_value']}"
                        f"{row['activity_units']}",
            "chembl_id": row["molecule_chembl_id"],
            "active": row["passes_threshold"] == "1",
        })

    smiles = generated + [r["smiles"] for r in records]
    names = [f"generated_{i:03d}" for i in range(n_generated)] + [
        r["name"] for r in records
    ]
    actives = [False] * n_generated + [r["active"] for r in records]
    return {
        "smiles": smiles,
        "names": names,
        "actives": actives,
        "n_generated": n_generated,
        "drugs": records,
        "box": docked["box"],
        "query": "argatroban",
    }


def fetch_poses(workflow_uuid: str, out_dir: Path) -> list[dict]:
    """Download every pose structure a docking workflow persisted.

    Returns one record per ligand with its downloaded PDB path, or an explicit
    note that the ligand produced no pose. Ligands that fail to dock are kept in
    the list rather than dropped, so the caller can count them.
    """
    from demo import nova

    rowan = nova.rowan_client()
    out_dir.mkdir(parents=True, exist_ok=True)
    structures = rowan.list_proteins(ancestor_uuid=workflow_uuid, size=1000)
    if not structures:
        raise RuntimeError(
            f"workflow {workflow_uuid} has no descendant structures. Either it was "
            "submitted with num_poses_to_save = 0 or it has not finished."
        )
    records = []
    for protein in structures:
        stem = str(protein.uuid)
        path = out_dir / f"{stem}.pdb"
        if not path.exists():
            protein.download_pdb_file(str(out_dir), stem)
        records.append({"uuid": stem, "name": getattr(protein, "name", None),
                        "path": str(path)})
    return records


def box_lining_residues(
    receptor: Receptor, chain: str, box: list[list[float]]
) -> list[int | str]:
    """Every residue of `chain` with a heavy atom inside the docking box.

    This is the honest bit space for poses produced by that box: it is fixed by
    the docking configuration, it cannot be influenced by an activity label, and
    it is decided before a single pose exists. `box` is Rowan's
    [[cx, cy, cz], [sx, sy, sz]] -- centre and full side lengths.
    """
    centre = np.asarray(box[0], dtype=float)
    half = np.asarray(box[1], dtype=float) / 2.0
    inside = [
        label for (ch, label), entry in receptor.residues.items()
        if ch == chain and bool((np.abs(entry["xyz"] - centre) <= half).all(axis=1).any())
    ]
    if not inside:
        raise ValueError(
            f"no residue of chain {chain!r} lies inside the docking box {box}. "
            "Refusing to return an empty fingerprint axis."
        )
    return sorted(inside, key=residue_order)


def run_head_to_head(
    *, approve: bool = False, max_credits: int = 30, wait_timeout: int = 5400
) -> dict:
    """The whole comparison, end to end. Submits only with approve=True.

    Ranks every docked molecule twice -- by IFP similarity to argatroban's pose
    and by ECFP4 Tanimoto to argatroban's SMILES -- and asks which ranking puts
    the measured thrombin binders higher. Writes
    results/demo/thrombin/ifp_head_to_head.json.
    """
    import json

    from demo import nova

    ligands = thrombin_ligand_set()
    submission = redock_with_poses(
        ligands["smiles"], n_poses=1, approve=approve, max_credits=max_credits,
        name="F2 arm S batch docking, enlarged active set, poses saved",
    )
    if not submission.get("uuid"):
        return {"stage": "not submitted", "submission": submission}

    uuid = submission["uuid"]
    record = nova.wait(uuid, label="batch docking with poses", timeout=wait_timeout)
    scores = record["data"]["best_scores"]

    pose_dir = REPO / "results/demo/thrombin/structures/poses"
    fetched = fetch_poses(uuid, pose_dir)

    receptor = parse_receptor(REPO / "results/demo/thrombin/structures/target_thrombin.pdb")
    site = json.loads((REPO / "results/demo/thrombin/02_site.json").read_text())

    # Two axes, both fixed before any pose existed and both reported.
    #
    # PRIMARY is the ensemble axis: every chain-B residue that lines the docking
    # box. It is set by the docking configuration alone -- no activity label
    # touches it -- and it is the site the ligands were actually docked into.
    #
    # SECONDARY is the 20-residue site from 02_site.json, for consistency with
    # the rest of the pipeline. It is reported as a sensitivity check rather
    # than as the headline because it was built as a 6 A shell around the
    # catalytic Ser195 and, measurably, does NOT contain Asp189 (structure
    # position 199 here; the GDSGGP motif at 203-208 fixes the offset). Asp189
    # is the base of the S1 pocket and the single residue that most
    # distinguishes a thrombin binder -- it is where PPACK's ionic bit fires in
    # the 1PPB validation. An axis without it is handicapped for this question.
    axes = {
        "box_lining": SiteAxis(
            chain=site["target_chain"],
            residues=box_lining_residues(receptor, site["target_chain"],
                                         ligands["box"]),
        ),
        "site_02": SiteAxis(
            chain=site["target_chain"],
            residues=[int(x) for x in site["chosen"]["residue_positions_1based"]],
        ),
    }
    axis = axes["box_lining"]

    # Fingerprint every pose. A ligand that failed to dock has no pose and is
    # kept with an all-zero fingerprint and flagged, not dropped.
    prints, failed = {}, []
    by_name = {rec["name"]: rec for rec in fetched}
    for i, name in enumerate(ligands["names"]):
        record_i = by_name.get(name) or (fetched[i] if i < len(fetched) else None)
        if record_i is None:
            failed.append(name)
            continue
        poses = poses_from_pose_file(Path(record_i["path"]), ligands["smiles"][i], name)
        prints[name] = fingerprint(poses, receptor, axis)

    query = ligands["query"]
    if query not in prints:
        raise RuntimeError(
            f"{query} produced no pose, so there is no query to rank against. "
            f"{len(failed)} of {len(ligands['names'])} ligands failed to dock."
        )

    q_fp = prints[query]
    q_smiles = ligands["smiles"][ligands["names"].index(query)]

    names, ifp_scores, ecfp_scores, actives = [], [], [], []
    for i, name in enumerate(ligands["names"]):
        if name == query or name not in prints:
            continue
        names.append(name)
        ifp_scores.append(similarity(prints[name], q_fp, "sokal_michener"))
        ecfp_scores.append(ecfp4_similarity(ligands["smiles"][i], q_smiles))
        actives.append(ligands["actives"][i])

    result = head_to_head(ifp_scores, ecfp_scores, actives, names)
    result["per_measure_auc"] = {
        measure: auc([MEASURES[measure](prints[n].bits, q_fp.bits) for n in names],
                     np.asarray(actives, bool))
        for measure in MEASURES
    }
    result["n_failed_to_dock"] = len(failed)
    result["failed_to_dock"] = failed
    result["workflow_uuid"] = uuid
    result["credits_charged"] = record.get("credits_charged")
    result["axis"] = {"name": "box_lining", "chain": axis.chain,
                      "n_residues": len(axis.residues), "n_bits": len(axis),
                      "contains_asp189": 199 in axis.residues}

    # The sensitivity run on the narrower pipeline axis, reported either way.
    alt = axes["site_02"]
    alt_prints = {n: fingerprint(p, receptor, alt) for n, p in prints.items()}
    alt_scores = [similarity(alt_prints[n], alt_prints[query], "sokal_michener")
                  for n in names]
    result["sensitivity_site_02"] = {
        "n_residues": len(alt.residues), "n_bits": len(alt),
        "contains_asp189": 199 in alt.residues,
        "auc_ifp": auc(alt_scores, np.asarray(actives, bool)),
        "delta": paired_bootstrap_delta_auc(alt_scores, ecfp_scores,
                                            np.asarray(actives, bool)),
    }
    result["docking_scores"] = {n: scores[ligands["names"].index(n)] for n in names}

    out = REPO / "results/demo/thrombin/ifp_head_to_head.json"
    out.write_text(json.dumps(result, indent=1, default=str))
    return result


def poses_from_pose_file(path: Path, smiles: str, name: str) -> Pose:
    """Read one docked pose, whatever container Rowan handed back.

    SDF carries bond orders; a PDB does not, so the ligand SMILES is used as the
    template. Fails loudly rather than returning a pose of the wrong molecule.
    """
    if path.suffix.lower() == ".sdf":
        poses = poses_from_sdf(path)
        return poses[0]
    text = path.read_text()
    hetatm = [
        line for line in text.splitlines()
        if line.startswith("HETATM") and line[17:20].strip() not in _SOLVENT
    ]
    if not hetatm:
        raise ValueError(
            f"{path.name}: no ligand HETATM records. Rowan returned a structure "
            "with no pose in it."
        )
    return pose_from_pdb_block("\n".join(hetatm) + "\nEND\n", smiles, name)


def main() -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-poses", action="store_true",
                        help="ask Rowan whether the docked poses were persisted")
    parser.add_argument("--ligand-set", action="store_true",
                        help="print the enlarged docking set; offline, free")
    parser.add_argument("--head-to-head", action="store_true",
                        help="re-dock with poses and run the comparison; SPENDS CREDITS")
    parser.add_argument("--approve", action="store_true",
                        help="required for --head-to-head to actually submit")
    parser.add_argument("--max-credits", type=int, default=30)
    args = parser.parse_args()

    if args.ligand_set:
        ligands = thrombin_ligand_set()
        print(f"{len(ligands['smiles'])} ligands: {ligands['n_generated']} generated "
              f"decoys + {len(ligands['drugs'])} measured F2 binders")
        print(f"actives (repo passes_threshold == 1): {sum(ligands['actives'])}")
        for record in ligands["drugs"]:
            print(f"  {'ACTIVE' if record['active'] else '  weak'}  "
                  f"{record['name']:22s} pchembl={record['pchembl']:.2f}  "
                  f"{record['activity']}  {record['chembl_id']}")
        print(f"estimated {19.23 * len(ligands['smiles']) / 203:.1f} credits")
        return 0

    if args.head_to_head:
        result = run_head_to_head(approve=args.approve, max_credits=args.max_credits)
        print(json.dumps(result, indent=2, default=str))
        return 0 if result.get("n") else 1

    if args.check_poses:
        status = check_saved_poses()
        print(json.dumps(status, indent=2))
        if not status["poses_available"]:
            print(
                "\nNo poses were saved. rowan.submit_batch_docking_workflow does not\n"
                "expose num_poses_to_save and stjames.BatchDockingWorkflow defaults it\n"
                "to 0, so the 203 Vina poses were scored and discarded. There is nothing\n"
                "to fingerprint until the batch is re-docked with pose saving on:\n"
                "  demo.ifp.redock_with_poses(n_poses=1, approve=True)  # ~20 credits\n"
            )
            return 1
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
