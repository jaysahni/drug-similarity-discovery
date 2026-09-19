"""Interface-contact toolkit: describe a binding site in target-side coordinates.

PROJECT_GOAL.md v0.2 argues that a binding site is best described by *which target
residues a ligand engages and how*, because that description is comparable across
molecules that share no chemistry. This module builds that object from real
structures.

    load_structure(path)                -> Bio.PDB structure
    het_inventory(structure)            -> every non-polymer residue + keep/drop reason
    drug_like_ligands(structure)        -> the candidates that survive the filter
    ligand_contacts(structure, resname) -> target-side contacts for one ligand
    consensus_signature([contacts])     -> InterfaceSignature (PROJECT_GOAL.md 4.3)
    weighted_jaccard(a, b)              -> signature-to-signature overlap (4.4)
    coverage(sig, observed_residues)    -> core_coverage, the primary score (4.4)

BoltzGen is unavailable here (no GPU, no Rowan key), so `consensus_signature`
is fed contacts from crystal structures rather than from designs. That is the
`source="known_ligand"` builder of PROJECT_GOAL.md C8 — the upper bound of the
8.3 ablation — and it exercises exactly the same code path a design ensemble
would.

RESIDUE NUMBERING
    Every residue id in this module is the **author sequence number**
    (mmCIF `_atom_site.auth_seq_id`, what Bio.PDB exposes as `residue.id[1]`).
    No SIFTS mapping is applied. For the VEGFR2 entries used in the self-test the
    author numbering equals UniProt P35968 numbering; `check_author_numbering`
    verifies that against the UniProt sequence rather than assuming it. Insertion
    codes are dropped from the id and reported in `warnings` if any occur, so a
    structure that has them must be checked before its ids are compared.

Usage:
    ./env/bin/python scripts/interfaces.py 4ASD 4AG8      # prints a report
"""

from __future__ import annotations

import copy
import json
import ssl
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

import certifi
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"

# The python.org 3.14 framework build ships no CA bundle (see fetch_chembl.py).
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

# --------------------------------------------------------------------------
# what counts as a ligand
# --------------------------------------------------------------------------
# Explicit exclusion list. Grouped by reason so `het_inventory` can say *why*
# something was dropped instead of silently omitting it. Codes are PDB chemical
# component IDs.
EXCLUDED_HET = {}


def _exclude(reason, codes):
    for c in codes.split():
        EXCLUDED_HET[c] = reason


_exclude("water", "HOH DOD WAT")
_exclude(
    "ion or simple salt",
    "SO4 PO4 NO3 CO3 ACT ACY FMT CL BR IOD F NA K LI CS RB MG CA SR BA MN FE FE2 "
    "CO NI CU CU1 ZN CD HG AU PT PB AG TL SM EU GD YB OH NH4 AZI CYN SCN VO4 WO4 MOO",
)
_exclude(
    "cryoprotectant, precipitant or buffer additive",
    "GOL EDO PEG PGE PG4 1PE 2PE P6G 7PE 12P XPE MPD DMS DMF TRS EPE MES BTB BIS "
    "TAM CIT FLC TLA MLA MLI SIN SUC IPA MOH ETA EOH BU3 BU1 TBU IMD DIO DOX BME "
    "DTT DTU MRD NHE CAC PIN POP UNL UNX ACE NH2 ACN SOG LDA C8E BOG BNG OCT HEZ "
    "PE4 P33 PEU PIG SPD SPM PUT",
)
_exclude(
    "glycan or covalent sugar modification",
    "NAG NDG BMA MAN BGC GLC GAL FUC FUL XYP XYS SIA NGA A2G RAM GLA",
)
# Cofactors are real binders, not crystallisation junk, but they are also not
# drug-like *candidates*: a kinase with ATP bound tells you about the natural
# substrate. Excluded by default, recoverable with include_cofactors=True.
_exclude(
    "cofactor or nucleotide (real binder, not a drug-like candidate; "
    "pass include_cofactors=True to keep)",
    "ATP ADP AMP ANP ACP AGS APC ADX GTP GDP GNP GSP GCP CTP UTP UDP UMP TTP "
    "NAD NAI NAP NDP NAJ FAD FMN FDA SAM SAH SFG COA ACO MCA HEM HEC HEA SRM "
    "PLP PMP TPP TDP B12 BTN MTA MTE MGD F43 CLA BCL CHL PQQ",
)
COFACTOR_REASON_PREFIX = "cofactor or nucleotide"

# Residues that belong to the polymer even though the parser flags them HETATM.
# Anything here is target, never ligand. Checked for, not assumed.
POLYMER_HET = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "CME", "KCX", "LLP", "MLY", "M3L",
    "ALY", "HYP", "PCA", "CGU", "SAC", "OCS", "CAS", "NEP", "HIC", "TYS", "SNN",
}

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y",
    "CSO": "C", "CME": "C", "KCX": "K", "MLY": "K", "HYP": "P", "PCA": "E",
}

# --------------------------------------------------------------------------
# interaction criteria
# --------------------------------------------------------------------------
# These are the standard geometric criteria quoted in the interaction-detection
# literature (PLIP/LigPlot-style distance cuts). This is a DISTANCE-BASED
# APPROXIMATION, not a validated reimplementation of PLIP: with the sole
# exception of the ring-centroid term, no angular criterion is applied, so
# H-bonds are not checked for D-H...A angle, halogen bonds are not checked for
# C-X...A angle, and aromatic contacts are not split into pi-stacked vs
# T-shaped. Hydrogens are absent from these crystal models, so donor/acceptor
# roles cannot be assigned from coordinates at all. Every type below should be
# read as "geometrically compatible with", not "is".
CUTOFF_HYDROPHOBIC = 4.0   # C...C
CUTOFF_HBOND = 3.5         # N/O/F ... N/O/F heavy-atom separation
CUTOFF_SALT_BRIDGE = 4.0   # charged protein group ... oppositely charged ligand atom
CUTOFF_AROMATIC = 5.5      # ring centroid ... ring centroid
CUTOFF_HALOGEN = 3.5       # Cl/Br/I ... N/O
HBOND_ELEMENTS = {"N", "O", "F"}
HALOGENS = {"CL", "BR", "I"}

# Protein side chains treated as formally charged. HIS is included as a
# *potential* cation: its protonation state is not determined by the model.
CATIONIC_N = {
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}
ANIONIC_O = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
}

AROMATIC_RINGS = {
    "PHE": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "TYR": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "HIS": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "TRP": [
        ("CG", "CD1", "NE1", "CE2", "CD2"),
        ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    ],
}


# --------------------------------------------------------------------------
# fetching (cached under data/raw/, which is gitignored)
# --------------------------------------------------------------------------
def _download(url, dest):
    """Download to a temp file and rename, so an interrupted fetch is not cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=60, context=SSL_CTX) as r:
        body = r.read()
        declared = r.headers.get("Content-Length")
    if not body:
        raise OSError(f"empty response from {url}")
    if declared is not None and len(body) != int(declared):
        raise OSError(f"truncated download from {url}: "
                      f"{len(body)} of {declared} bytes")
    tmp.write_bytes(body)
    tmp.replace(dest)
    return dest


def fetch_structure(pdb_id, cache_dir=RAW / "structures", fmt="cif"):
    """Coordinates for a PDB entry, cached. fmt is "cif" (default) or "pdb"."""
    if fmt not in ("cif", "pdb"):
        raise ValueError(f"fmt must be 'cif' or 'pdb', got {fmt!r}")
    pdb_id = pdb_id.lower()
    dest = Path(cache_dir) / f"{pdb_id}.{fmt}"
    if not dest.exists():
        _download(f"https://files.rcsb.org/download/{pdb_id}.{fmt}", dest)
    return dest


def fetch_ccd(comp_id, cache_dir=RAW / "ccd"):
    """Chemical Component Dictionary entry for a ligand code, cached."""
    comp_id = comp_id.upper()
    dest = Path(cache_dir) / f"{comp_id}.cif"
    if not dest.exists():
        _download(f"https://files.rcsb.org/ligands/download/{comp_id}.cif", dest)
    return dest


def fetch_uniprot_sequence(accession, cache_dir=RAW / "uniprot"):
    """Canonical sequence for a UniProt accession, cached. Returns the string."""
    dest = Path(cache_dir) / f"{accession}.fasta"
    if not dest.exists():
        _download(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", dest)
    return "".join(dest.read_text().splitlines()[1:])


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------
def load_structure(path, structure_id=None):
    """Parse an mmCIF (.cif/.mmcif) or PDB (.pdb/.ent) file with Bio.PDB.

    Disordered atoms keep the parser's default (highest-occupancy) altloc.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    elif suffix in (".pdb", ".ent"):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"unknown structure format {suffix!r} for {path}")
    return parser.get_structure(structure_id or path.stem, str(path))


def is_polymer_residue(residue):
    """True for target-side residues: ATOM records plus known modified residues.

    A HETATM-flagged residue is polymer only if it is on POLYMER_HET. Matching
    bare amino-acid names here would swallow a free amino acid deposited as a
    ligand: it would vanish from het_inventory and its atoms would be counted as
    target. The flag wins; POLYMER_HET is the explicit exception list.
    """
    if residue.id[0] == " ":
        return True
    return residue.get_resname().strip().upper() in POLYMER_HET


def heavy_atoms(residue):
    return [a for a in residue.get_atoms() if a.element.upper() not in ("H", "D")]


def _formula(atoms):
    counts = defaultdict(int)
    for a in atoms:
        counts[a.element.upper().capitalize()] += 1
    order = sorted(counts, key=lambda e: (e != "C", e))
    return "".join(e if counts[e] == 1 else f"{e}{counts[e]}" for e in order)


def het_inventory(structure, model=0, min_heavy_atoms=10, include_cofactors=False):
    """Every non-polymer residue in the model with a keep/drop verdict and reason.

    This is how absence gets *verified*: the caller can see that a structure
    contained nothing but waters, rather than assuming it.
    """
    out = []
    for chain in structure[model]:
        for res in chain:
            if is_polymer_residue(res):
                continue
            name = res.get_resname().strip().upper()
            atoms = heavy_atoms(res)
            rec = {
                "resname": name,
                "chain_id": chain.id,
                "resseq": res.id[1],
                "n_heavy_atoms": len(atoms),
                "formula": _formula(atoms),
            }
            reason = EXCLUDED_HET.get(name)
            if reason and reason.startswith(COFACTOR_REASON_PREFIX) and include_cofactors:
                reason = None
            if reason:
                rec["kept"], rec["reason"] = False, f"exclusion list: {reason}"
            elif len(atoms) < min_heavy_atoms:
                rec["kept"] = False
                rec["reason"] = f"{len(atoms)} heavy atoms < min_heavy_atoms={min_heavy_atoms}"
            else:
                rec["kept"], rec["reason"] = True, "drug-like candidate"
            out.append(rec)
    return sorted(out, key=lambda r: (r["chain_id"], r["resseq"]))


def drug_like_ligands(structure, model=0, min_heavy_atoms=10, include_cofactors=False):
    """Candidate ligands: >= min_heavy_atoms heavy atoms and not on the exclusion list.

    Use `het_inventory` for the full list including what was dropped and why.
    """
    inv = het_inventory(structure, model, min_heavy_atoms, include_cofactors)
    return [r for r in inv if r["kept"]]


# --------------------------------------------------------------------------
# ligand chemistry from the CCD (formal charges + aromatic rings via RDKit)
# --------------------------------------------------------------------------
_BOND_ORDER = {"SING": 1.0, "DOUB": 2.0, "TRIP": 3.0, "QUAD": 4.0, "AROM": 1.5}


def ccd_mol(comp_id):
    """RDKit mol built from the CCD definition, atoms tagged with their PDB names.

    The CCD gives element, formal charge and Kekule bond orders per *named* atom,
    so the mol maps back onto structure atoms by name. Returns None if the
    component cannot be fetched or sanitised; callers must degrade gracefully.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    try:
        d = MMCIF2Dict(str(fetch_ccd(comp_id)))
        names = [_clean_name(n) for n in d["_chem_comp_atom.atom_id"]]
        elements = [e.capitalize() for e in d["_chem_comp_atom.type_symbol"]]
        charges = d.get("_chem_comp_atom.charge", ["0"] * len(names))
    except Exception:  # noqa: BLE001 - network or malformed component
        return None

    try:
        return _build_ccd_mol(Chem, d, names, elements, charges)
    except Exception:  # noqa: BLE001 - CCD chemistry RDKit refuses to build
        return None


def _build_ccd_mol(Chem, d, names, elements, charges):
    from rdkit.Chem import BondType

    mol = Chem.RWMol()
    index = {}
    for name, elem, chg in zip(names, elements, charges):
        atom = Chem.Atom(elem)
        try:
            atom.SetFormalCharge(int(float(chg)))
        except ValueError:  # '?' or '.'
            pass
        atom.SetProp("pdb_name", name)
        index[name] = mol.AddAtom(atom)

    a1 = [_clean_name(n) for n in d.get("_chem_comp_bond.atom_id_1", [])]
    a2 = [_clean_name(n) for n in d.get("_chem_comp_bond.atom_id_2", [])]
    orders = d.get("_chem_comp_bond.value_order", [])
    for x, y, order in zip(a1, a2, orders):
        if x not in index or y not in index or mol.GetBondBetweenAtoms(index[x], index[y]):
            continue
        bt = {1.0: BondType.SINGLE, 2.0: BondType.DOUBLE, 3.0: BondType.TRIPLE,
              4.0: BondType.QUADRUPLE, 1.5: BondType.AROMATIC}[_BOND_ORDER.get(order, 1.0)]
        mol.AddBond(index[x], index[y], bt)

    out = mol.GetMol()
    Chem.SanitizeMol(out)
    return out


def _charge_groups(mol, name_of):
    """Formal charges reduced to the ones that can actually make a salt bridge.

    Two corrections to the raw CCD formal charges, both at the level of the local
    charged group (a heteroatom, its bonded neighbour, and that neighbour's
    same-element/same-degree substituents):

      * group net zero -> drop. The CCD writes a nitro group charge-separated as
        N(+1)-O(-1)-O(0). Taken atom by atom that oxygen looks anionic, and every
        nitro drug would salt-bridge to every arginine. The group is neutral, so
        no charge survives.
      * group net non-zero -> spread. A carboxylate is C(0)-O(-1)-O(0); a contact
        through the uncharged oxygen is the same interaction, so it inherits the
        charge.

    This is a heuristic on top of deposited formal charges, not a pKa model: a
    group the CCD models neutral (most drug carboxylic acids and amines) stays
    neutral here and its salt bridges are missed.
    """
    charges = {}
    for atom in mol.GetAtoms():
        q = atom.GetFormalCharge()
        if q == 0 or atom.GetIdx() not in name_of:
            continue
        # the local group: this atom, its bonded neighbours, and their substituents
        group, spread = {atom.GetIdx()}, []
        for nbr in atom.GetNeighbors():
            group.add(nbr.GetIdx())
            for x in nbr.GetNeighbors():
                if x.GetIdx() == atom.GetIdx():
                    continue
                group.add(x.GetIdx())
                if (atom.GetSymbol() in ("O", "N", "S")
                        and x.GetIdx() in name_of
                        and x.GetSymbol() == atom.GetSymbol()
                        and x.GetFormalCharge() == 0
                        and x.GetDegree() == atom.GetDegree()):
                    spread.append(x)
        net = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in group)
        if len(group) > 1 and net * q <= 0:
            continue  # locally neutralised, e.g. a nitro group
        charges[name_of[atom.GetIdx()]] = q
        for sib in spread:
            charges.setdefault(name_of[sib.GetIdx()], q)
    return charges


def _ligand_chemistry(comp_id):
    """(charge_by_atom_name, aromatic_rings_as_name_tuples, parsed_flag)."""
    mol = ccd_mol(comp_id)
    if mol is None:
        return {}, [], False
    name_of = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("pdb_name"):
            name_of[atom.GetIdx()] = atom.GetProp("pdb_name")
    charges = _charge_groups(mol, name_of)
    rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            rings.append(tuple(name_of[i] for i in ring if i in name_of))
    return charges, rings, True


# --------------------------------------------------------------------------
# buried area
# --------------------------------------------------------------------------
_SASA = ShrakeRupley()


def _sasa_by_residue(model, keep):
    """Per-residue SASA over the residues `keep(chain_id, residue)` selects."""
    work = copy.deepcopy(model)
    work.detach_parent()
    for chain in list(work):
        for res in list(chain):
            if not keep(chain.id, res):
                chain.detach_child(res.id)
                continue
            for atom in list(res):
                # ShrakeRupley KeyErrors on elements outside its radii table.
                if atom.element.upper() not in _SASA.radii_dict:
                    res.detach_child(atom.id)
        if len(chain) == 0:
            work.detach_child(chain.id)
    _SASA.compute(work, level="R")
    # keyed by (chain_id, auth_seq_id, insertion_code) - dropping the insertion
    # code here would silently collapse 100 and 100A into one number
    return {(c.id, r.id[1], r.id[2]): r.sasa for c in work for r in c}


def buried_area_by_residue(structure, ligand_key, model=0):
    """Delta-SASA per target residue on ligand binding, in square angstroms.

    Shrake-Rupley (Bio.PDB.SASA, probe 1.4 A, 100 points) over the polymer alone
    minus the polymer with this one ligand present. Waters and all other
    heteroatoms are removed from both states, and hydrogens are absent from these
    crystal models, so this is a *proxy* for buried area rather than a rigorously
    protonated ΔSASA.
    """
    m = structure[model]
    apo = _sasa_by_residue(m, lambda cid, r: is_polymer_residue(r))
    holo = _sasa_by_residue(
        m,
        lambda cid, r: is_polymer_residue(r) or (cid, r.id[1], r.get_resname().strip()) == ligand_key,
    )
    return {k: float(apo[k] - holo.get(k, apo[k])) for k in apo}


# --------------------------------------------------------------------------
# contacts
# --------------------------------------------------------------------------
def _clean_name(s):
    """Atom name with a surrounding quote pair removed, primes preserved.

    mmCIF quotes names containing a prime as "C5'". Stripping quote characters
    blindly turns C5' into C5, which collides with a genuinely different atom.
    """
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        s = s[1:-1]
    return s


def _centroid(residue, names):
    pts = [residue[n].coord for n in names if n in residue]
    return np.mean(pts, axis=0) if len(pts) == len(names) else None


def ligand_contacts(structure, ligand_resname, chain=None, cutoff=4.5,
                    resseq=None, model=0, with_buried_area=True):
    """Target-side residues within `cutoff` A of any heavy atom of one ligand copy.

    Returns
        residue_ids : sorted author sequence numbers of contacting residues
        per_residue : {residue_id: {resname, chain_id, min_distance,
                       n_atom_contacts, buried_area_proxy, interaction_types}}
        ligand_info : {resname, chain, resseq, n_heavy_atoms, formula, ...}
        warnings    : anything the caller must know before trusting the ids

    Interaction types are the distance-based approximation documented at the top
    of this module. `buried_area_proxy` is delta-SASA (see buried_area_by_residue),
    or None for every residue when with_buried_area=False -- None rather than 0.0,
    so an uncomputed area is never reported as a measured one. The returned dict
    carries `buried_area_computed` to say which it is.

    The aromatic centroid criterion is evaluated only for residues that already
    have a heavy atom within `cutoff`, so an offset ring stack whose closest
    heavy-atom approach exceeds `cutoff` is not typed. Use
    `aromatic_centroid_candidates` to measure whether that truncation bites on a
    given structure.
    """
    m = structure[model]
    want = ligand_resname.strip().upper()
    copies = [
        (c.id, r)
        for c in m
        for r in c
        if not is_polymer_residue(r) and r.get_resname().strip().upper() == want
        and (chain is None or c.id == chain)
        and (resseq is None or r.id[1] == resseq)
    ]
    if not copies:
        raise KeyError(f"no ligand {want!r} in model {model}"
                       + (f" chain {chain!r}" if chain else ""))

    warnings = []
    if len(copies) > 1:
        others = [f"{c}/{r.id[1]}" for c, r in copies[1:]]
        warnings.append(
            f"{len(copies)} copies of {want}; using chain {copies[0][0]} "
            f"resseq {copies[0][1].id[1]}, ignoring {others} "
            "(pass chain=/resseq= to select another)"
        )
    lig_chain, lig_res = copies[0]
    lig_atoms = heavy_atoms(lig_res)
    ligand_key = (lig_chain, lig_res.id[1], lig_res.get_resname().strip())

    target_atoms = [
        a
        for c in m
        for r in c
        if is_polymer_residue(r)
        for a in heavy_atoms(r)
    ]
    if not target_atoms:
        raise ValueError("structure has no polymer residues to contact")

    charges, lig_rings, parsed = _ligand_chemistry(want)
    if not parsed:
        warnings.append(
            f"CCD entry for {want} unavailable or unparseable: ligand formal "
            "charges and aromatic rings unknown, so salt_bridge and aromatic "
            "typing are not attempted for this ligand"
        )

    search = NeighborSearch(target_atoms)
    per_residue = {}
    res_obj = {}
    icodes = set()
    keys_by_rid = defaultdict(set)
    for la in lig_atoms:
        le = la.element.upper()
        lname = _clean_name(la.get_name())
        lcharge = charges.get(lname, 0)
        for ta in search.search(la.coord, cutoff):
            tres = ta.get_parent()
            tchain = tres.get_parent().id
            rid = tres.id[1]
            if tres.id[2] != " ":
                icodes.add(f"{tchain}:{rid}{tres.id[2]}")
            d = float(np.linalg.norm(la.coord - ta.coord))
            te = ta.element.upper()
            tname = ta.get_name().strip()
            tresname = tres.get_resname().strip().upper()

            rec = per_residue.setdefault(
                rid,
                {"resname": tresname, "chain_id": tchain, "min_distance": d,
                 "n_atom_contacts": 0, "buried_area_proxy": None,
                 "interaction_types": set()},
            )
            res_obj[rid] = tres
            keys_by_rid[rid].add((tchain, tres.id[2]))
            rec["n_atom_contacts"] += 1
            rec["min_distance"] = min(rec["min_distance"], d)

            types = rec["interaction_types"]
            if te == "C" and le == "C" and d <= CUTOFF_HYDROPHOBIC:
                types.add("hydrophobic")
            if te in HBOND_ELEMENTS and le in HBOND_ELEMENTS and d <= CUTOFF_HBOND:
                types.add("hbond")
            if le in HALOGENS and te in ("N", "O") and d <= CUTOFF_HALOGEN:
                types.add("halogen")
            if d <= CUTOFF_SALT_BRIDGE and (
                ((tresname, tname) in CATIONIC_N and lcharge < 0)
                or ((tresname, tname) in ANIONIC_O and lcharge > 0)
            ):
                types.add("salt_bridge")

    # aromatic: ring centroid to ring centroid, the one angular-free term that
    # cannot be decided atom-pairwise
    if lig_rings:
        lig_centroids = []
        for ring in lig_rings:
            pts = [a.coord for a in lig_atoms
                   if _clean_name(a.get_name()) in ring]
            if len(pts) == len(ring):
                lig_centroids.append(np.mean(pts, axis=0))
        for rid, rec in per_residue.items():
            if rec["resname"] not in AROMATIC_RINGS:
                continue
            tres = res_obj[rid]
            for names in AROMATIC_RINGS[rec["resname"]]:
                c = _centroid(tres, names)
                if c is None:
                    continue
                if any(np.linalg.norm(c - lc) <= CUTOFF_AROMATIC for lc in lig_centroids):
                    rec["interaction_types"].add("aromatic")

    if with_buried_area:
        buried = buried_area_by_residue(structure, ligand_key, model)
        for rid, rec in per_residue.items():
            # sum over every (chain, icode) this record merged, so nothing is lost
            rec["buried_area_proxy"] = round(sum(
                float(buried.get((ch, rid, ic), 0.0)) for ch, ic in keys_by_rid[rid]
            ), 3)

    for rec in per_residue.values():
        rec["interaction_types"] = sorted(rec["interaction_types"])
        rec["min_distance"] = round(rec["min_distance"], 3)

    chains_engaged = sorted({ch for ks in keys_by_rid.values() for ch, _ in ks})
    if len(chains_engaged) > 1:
        warnings.append(
            f"contacts span polymer chains {chains_engaged}; residue ids are "
            "author numbers and are NOT unique across chains here"
        )
    chain_clash = sorted(r for r, ks in keys_by_rid.items()
                         if len({ch for ch, _ in ks}) > 1)
    if chain_clash:
        warnings.append(
            f"residue-id collision across chains for {chain_clash}: the "
            "per_residue record for each of these ids MERGES contacts from "
            f"chains {chains_engaged}. Author numbering alone cannot separate "
            "them; split by chain before comparing such a site across structures"
        )
    icode_clash = sorted(r for r, ks in keys_by_rid.items()
                         if len({ic for _, ic in ks}) > 1)
    if icode_clash:
        warnings.append(
            f"residue-id collision across insertion codes for {icode_clash}: "
            "the per_residue record for each of these ids MERGES residues that "
            "differ only by insertion code"
        )
    if icodes:
        warnings.append(f"insertion codes dropped from residue ids: {sorted(icodes)}")

    return {
        "residue_ids": sorted(per_residue),
        "per_residue": per_residue,
        "ligand_info": {
            "resname": want,
            "chain": lig_chain,
            "resseq": lig_res.id[1],
            "n_heavy_atoms": len(lig_atoms),
            "formula": _formula(lig_atoms),
            "ccd_parsed": parsed,
            "n_aromatic_rings": len(lig_rings),
            "n_charged_atoms": len(charges),
        },
        "chains_engaged": chains_engaged,
        "cutoff": cutoff,
        "buried_area_computed": bool(with_buried_area),
        "numbering": "author (auth_seq_id); no SIFTS mapping applied",
        "structure_id": structure.id,
        "warnings": warnings,
    }


def aromatic_centroid_candidates(structure, ligand_resname, chain=None,
                                 resseq=None, centroid_cutoff=CUTOFF_AROMATIC,
                                 model=0):
    """Residue ids whose aromatic ring centroid is within centroid_cutoff of a
    ligand aromatic ring centroid, ignoring the heavy-atom contact cutoff.

    `ligand_contacts` types aromatic contacts only for residues that already
    passed its heavy-atom `cutoff`. Comparing this set against that one measures
    how much that pre-filter truncates the centroid criterion, rather than
    assuming it does not.
    """
    m = structure[model]
    want = ligand_resname.strip().upper()
    copies = [(c.id, r) for c in m for r in c
              if not is_polymer_residue(r)
              and r.get_resname().strip().upper() == want
              and (chain is None or c.id == chain)
              and (resseq is None or r.id[1] == resseq)]
    if not copies:
        raise KeyError(f"no ligand {want!r} in model {model}")
    lig_atoms = heavy_atoms(copies[0][1])
    _, lig_rings, _ = _ligand_chemistry(want)

    lig_centroids = []
    for ring in lig_rings:
        pts = [a.coord for a in lig_atoms if _clean_name(a.get_name()) in ring]
        if len(pts) == len(ring):
            lig_centroids.append(np.mean(pts, axis=0))

    out = set()
    for c in m:
        for r in c:
            if not is_polymer_residue(r):
                continue
            for names in AROMATIC_RINGS.get(r.get_resname().strip().upper(), []):
                cen = _centroid(r, names)
                if cen is None:
                    continue
                if any(np.linalg.norm(cen - lc) <= centroid_cutoff
                       for lc in lig_centroids):
                    out.add(r.id[1])
    return sorted(out)


# --------------------------------------------------------------------------
# signature
# --------------------------------------------------------------------------
def _per_residue(contacts):
    """per_residue with int keys, tolerating a dict round-tripped through JSON."""
    return {int(k): v for k, v in contacts["per_residue"].items()}


def consensus_signature(list_of_contact_dicts, core_threshold=0.6, source="known_ligand"):
    """InterfaceSignature (PROJECT_GOAL.md 4.3) over a set of contact maps.

    frequency            fraction of input complexes engaging that residue
    interaction_types    per type, fraction of ALL inputs showing it at that
                         residue (not of the inputs that engage it)
    mean_buried_area     mean delta-SASA over ALL inputs, counting 0 where the
                         residue is not engaged; mean_buried_area_when_engaged
                         is the conditional version. Both are None when any input
                         was extracted with with_buried_area=False, rather than
                         reporting an area nobody computed.
    is_core              frequency >= core_threshold

    `warnings` is non-empty when a residue id means different residues in
    different inputs, i.e. when the numbering is not actually comparable. No
    check is made that the inputs describe the same target at all -- that is the
    caller's responsibility; `members` records where each contact map came from.
    """
    if not list_of_contact_dicts:
        raise ValueError("consensus_signature needs at least one contact map")
    n = len(list_of_contact_dicts)

    have_area = all(c.get("buried_area_computed", True) for c in list_of_contact_dicts)
    seen = defaultdict(int)
    names, chains = {}, {}
    resname_seen = defaultdict(set)
    chain_seen = defaultdict(set)
    types = defaultdict(lambda: defaultdict(int))
    area = defaultdict(float)
    for contacts in list_of_contact_dicts:
        for rid, rec in _per_residue(contacts).items():
            seen[rid] += 1
            names.setdefault(rid, rec["resname"])
            chains.setdefault(rid, rec["chain_id"])
            resname_seen[rid].add(rec["resname"])
            chain_seen[rid].add(rec["chain_id"])
            if have_area:
                area[rid] += rec.get("buried_area_proxy") or 0.0
            for t in rec["interaction_types"]:
                types[rid][t] += 1

    residues = []
    for rid in sorted(seen):
        freq = seen[rid] / n
        residues.append({
            "residue_id": rid,
            "chain_id": chains[rid],
            "residue_name": names[rid],
            "frequency": freq,
            "interaction_types": {t: c / n for t, c in sorted(types[rid].items())},
            "mean_buried_area": round(area[rid] / n, 3) if have_area else None,
            "mean_buried_area_when_engaged": (round(area[rid] / seen[rid], 3)
                                              if have_area else None),
            "is_core": freq >= core_threshold,
        })

    # A residue id that resolves to different residues in different inputs means
    # the numbering is not comparable, which invalidates the whole consensus.
    warnings = []
    clash = {rid: sorted(v) for rid, v in resname_seen.items() if len(v) > 1}
    if clash:
        warnings.append(
            f"residue id disagreement across inputs for {sorted(clash)}: "
            f"{clash}. Author numbering does not line up between these "
            "structures; the consensus below is not trustworthy for those ids"
        )
    mixed = {rid: sorted(v) for rid, v in chain_seen.items() if len(v) > 1}
    if mixed:
        warnings.append(
            f"the same residue id came from different chain letters across "
            f"inputs: {mixed}. Residue names agree, so this is not a numbering "
            "failure, but check that those chains are the same entity."
        )

    members = [
        {**c["ligand_info"], "structure_id": c.get("structure_id")}
        for c in list_of_contact_dicts
    ]
    return {
        "source": source,
        "warnings": warnings,
        "n_inputs": n,
        "core_threshold": core_threshold,
        "residues": residues,
        "core_residue_ids": [r["residue_id"] for r in residues if r["is_core"]],
        "mean_buried_area": round(sum(area.values()) / n, 3) if have_area else None,
        "numbering": "author (auth_seq_id); no SIFTS mapping applied",
        "members": members,
    }


def _freqs(sig):
    return {r["residue_id"]: r["frequency"] for r in sig["residues"]}


def weighted_jaccard(sig_a, sig_b):
    """Frequency-weighted overlap: sum(min) / sum(max) over the residue union.

    1.0 for identical signatures, 0.0 for disjoint ones, 0.0 if both are empty.
    """
    fa, fb = _freqs(sig_a), _freqs(sig_b)
    keys = set(fa) | set(fb)
    if not keys:
        return 0.0
    num = sum(min(fa.get(k, 0.0), fb.get(k, 0.0)) for k in keys)
    den = sum(max(fa.get(k, 0.0), fb.get(k, 0.0)) for k in keys)
    return num / den if den else 0.0


def coverage(sig, observed_residues):
    """How much of a signature an observed contact set engages (PROJECT_GOAL.md 4.4).

    core_coverage is the primary score: fraction of core hotspots engaged. It is
    None, not 0.0, when the signature has no core residues at all -- an undefined
    score rather than a measured miss.
    `observed_residues` may be a list of ids or a contact dict from
    `ligand_contacts`.
    """
    if isinstance(observed_residues, dict) and "residue_ids" in observed_residues:
        observed_residues = observed_residues["residue_ids"]
    obs = {int(r) for r in observed_residues}

    core = [r["residue_id"] for r in sig["residues"] if r["is_core"]]
    allr = [r["residue_id"] for r in sig["residues"]]
    fa = _freqs(sig)
    tot = sum(fa.values())
    return {
        # None, not 0.0: with no core residues the primary score is undefined,
        # and writing a zero would put an uncomputed number into results/
        "core_coverage": len([r for r in core if r in obs]) / len(core) if core else None,
        "all_coverage": len([r for r in allr if r in obs]) / len(allr) if allr else 0.0,
        "weighted_coverage": (sum(f for r, f in fa.items() if r in obs) / tot) if tot else 0.0,
        "n_core": len(core),
        "n_observed": len(obs),
        "engaged_core_residues": [r for r in core if r in obs],
        "missed_core_residues": [r for r in core if r not in obs],
    }


# --------------------------------------------------------------------------
# numbering check
# --------------------------------------------------------------------------
def check_author_numbering(structure, chain_id, uniprot_seq, model=0):
    """Do author residue numbers index into the UniProt sequence correctly?

    Cross-structure comparison of residue ids is meaningless unless they mean the
    same thing, so this verifies the claim instead of assuming it.
    """
    n_cmp = n_match = 0
    mismatches = []
    for res in structure[model][chain_id]:
        if not is_polymer_residue(res):
            continue
        rid = res.id[1]
        one = AA3_TO_1.get(res.get_resname().strip().upper())
        if one is None or not (1 <= rid <= len(uniprot_seq)):
            continue
        n_cmp += 1
        if uniprot_seq[rid - 1] == one:
            n_match += 1
        elif len(mismatches) < 10:
            mismatches.append({"resid": rid, "structure": one,
                               "uniprot": uniprot_seq[rid - 1]})
    return {
        "n_compared": n_cmp,
        "n_match": n_match,
        "fraction_match": n_match / n_cmp if n_cmp else 0.0,
        "example_mismatches": mismatches,
    }


# --------------------------------------------------------------------------
# CLI - prints a report; writes nothing (scripts/test_interfaces.py owns results/)
# --------------------------------------------------------------------------
def main():
    ids = sys.argv[1:]
    if not ids:
        sys.exit("usage: ./env/bin/python scripts/interfaces.py <PDB_ID|path> ...")

    all_contacts = []
    for entry in ids:
        path = Path(entry) if Path(entry).exists() else fetch_structure(entry)
        t0 = time.perf_counter()
        st = load_structure(path)
        print(f"\n=== {entry} ({path.name}, parsed in {time.perf_counter() - t0:.2f}s) ===")

        inv = het_inventory(st)
        kept = [r for r in inv if r["kept"]]
        dropped = defaultdict(list)
        for r in inv:
            if not r["kept"]:
                dropped[r["reason"]].append(r["resname"])
        print(f"{len(inv)} non-polymer residues, {len(kept)} drug-like")
        for reason, names in sorted(dropped.items()):
            counts = {n: names.count(n) for n in sorted(set(names))}
            print(f"  excluded [{reason}]: {counts}")

        for lig in kept:
            t0 = time.perf_counter()
            c = ligand_contacts(st, lig["resname"], chain=lig["chain_id"],
                                resseq=lig["resseq"])
            dt = time.perf_counter() - t0
            all_contacts.append(c)
            info = c["ligand_info"]
            print(f"  {info['resname']} ({info['formula']}, {info['n_heavy_atoms']} "
                  f"heavy atoms) -> {len(c['residue_ids'])} residues in {dt:.2f}s")
            print(f"    {c['residue_ids']}")
            for w in c["warnings"]:
                print(f"    WARNING: {w}")

    if len(all_contacts) > 1:
        sig = consensus_signature(all_contacts)
        print(f"\nconsensus over n={sig['n_inputs']}: "
              f"{len(sig['core_residue_ids'])} core residues at "
              f"freq >= {sig['core_threshold']}")
        print(f"  core: {sig['core_residue_ids']}")
        print(json.dumps(sig["residues"][:3], indent=2))


if __name__ == "__main__":
    main()
