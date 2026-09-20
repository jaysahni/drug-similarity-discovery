"""Target structure + binding-site prep for the repurposing pipeline (task T5 / WS-B).

Given a UniProt accession it: pulls the canonical sequence, inventories every PDB
entry that maps to that accession (resolution, method, chain, bound ligands,
construct coverage), picks one primary structure by a stated rule, strips ALL
heteroatoms, runs P2Rank on the stripped structure, and writes the resulting
pocket - with grid-based geometry (volume, enclosure, concavity, max inscribed
sphere) for the addressability test in PROJECT_GOAL.md 1.4b / B10.

The known ligand is never shown to site detection (B2/B12): it is stripped before
P2Rank runs and its contact residues land in a separate, clearly namespaced file
that only validation reads. The overlap between P2Rank's top pocket and those
contacts is the measured recovery number, computed here across several holo
structures so the headline carries an n.

Usage:
    ./env/bin/python scripts/prep_target.py --uniprot P35968
    ./env/bin/python scripts/prep_target.py --uniprot P35968 --pdb 4ASD --robustness 0
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import ssl
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import warnings
from pathlib import Path

import certifi
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, cKDTree
from scipy.stats import hypergeom

warnings.filterwarnings("ignore")  # biopython's mmCIF/PDB construction warnings

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
SCRATCH = Path(
    "/private/tmp/claude-501/-Users-evanxiang-Desktop-Projects-drug-similarity-discovery"
    "/ff08da32-3cac-44f3-9083-f8fc1e02ca4d/scratchpad"
)
JAVA_HOME = SCRATCH / "tools" / "jdk-21.0.12.1+1" / "Contents" / "Home"
PRANK = SCRATCH / "tools" / "p2rank_2.5" / "prank"

UNIPROT = "https://rest.uniprot.org/uniprotkb"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
RCSB_FILES = "https://files.rcsb.org/download"

# The python.org 3.14 framework build ships no CA bundle (see fetch_chembl.py).
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

DEFAULT_SLUG = "colorectal-cancer"  # the demo disease in docs/02-PIPELINE-PLAN.md
CONTACT_CUTOFF = 4.5  # A, heavy-atom to heavy-atom, for known-ligand contacts
MIN_LIGAND_MW = 250.0  # below this a het group is a buffer/ion/cryoprotectant
# Additives that clear the MW cut but are not ligands in any useful sense.
NOT_LIGANDS = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "1PE", "MPD", "ACT", "FMT", "DMS",
    "TRS", "EPE", "MES", "BME", "DTT", "CIT", "TLA", "IMD", "NAG", "BMA", "MAN",
    "FUC", "HG", "CS", "CAC", "IOD", "AZI", "SCN", "NO3", "ACY", "BTB", "POL",
}
# van der Waals radii (A) for the elements that occur in protein heavy atoms.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90, "P": 1.80}
BACKBONE = ("N", "CA", "C")
PROBE = 1.4  # water probe radius


def is_polymer(residue):
    """A standard residue, or a modified one (PTR, SEP, TPO, CSO, CME, MSE...) that still
    carries a backbone. mmCIF flags those as heteroatoms, but they are part of the chain:
    deleting them opens a break in the activation loop the site sits in."""
    return residue.id[0] == " " or all(a in residue for a in BACKBONE)


def res_label(residue):
    """Author residue number with its insertion code, e.g. '919' or '100A'."""
    return f"{residue.id[1]}{residue.id[2].strip()}"


# --------------------------------------------------------------------------
# http with an on-disk cache under data/raw (gitignored) so reruns are offline
# --------------------------------------------------------------------------
def cached(path, fetch):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return path.read_bytes()
    blob = fetch()
    path.write_bytes(blob)
    return blob


def http(url, payload=None, retries=4):
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"Content-Type": "application/json"} if data else {}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=90, context=SSL_CTX) as r:
                blob, status = r.read(), r.status
            if not blob:
                raise ValueError(f"empty body (HTTP {status}) from {url}")
            return blob
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"unreachable: retries exhausted for {url}")


# --------------------------------------------------------------------------
# UniProt
# --------------------------------------------------------------------------
def fetch_uniprot(acc):
    fields = "accession,id,protein_name,gene_names,sequence,ft_domain"
    meta = json.loads(
        cached(
            RAW / "uniprot" / f"{acc}.json",
            lambda: http(f"{UNIPROT}/{acc}.json?{urllib.parse.urlencode({'fields': fields})}"),
        )
    )
    domains = [
        {
            "description": f.get("description", ""),
            "start": f["location"]["start"]["value"],
            "end": f["location"]["end"]["value"],
        }
        for f in meta.get("features", [])
        if f["type"] == "Domain"
    ]
    return {
        "accession": meta["primaryAccession"],
        "entry_name": meta["uniProtkbId"],
        "gene": (meta.get("genes") or [{}])[0].get("geneName", {}).get("value"),
        "protein_name": (
            meta["proteinDescription"].get("recommendedName")
            or (meta["proteinDescription"].get("submittedName") or [{}])[0]
        ).get("fullName", {}).get("value", meta["uniProtkbId"]),
        "sequence": meta["sequence"]["value"],
        "length": meta["sequence"]["length"],
        "domains": domains,
    }


def catalytic_domain(up):
    """The kinase (or single largest catalytic) domain, used as the construct."""
    kin = [d for d in up["domains"] if "kinase" in d["description"].lower()]
    if not kin:
        return None
    d = max(kin, key=lambda x: x["end"] - x["start"])
    return {"description": d["description"], "start": d["start"], "end": d["end"]}


# --------------------------------------------------------------------------
# RCSB inventory
# --------------------------------------------------------------------------
def search_entities(acc):
    q = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers."
                        "reference_sequence_identifiers.database_accession",
                        "operator": "in",
                        "value": [acc],
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers."
                        "reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 500}, "results_verbosity": "compact"},
    }
    def go():
        ids, start = [], 0
        while True:
            q["request_options"]["paginate"] = {"start": start, "rows": 500}
            page = json.loads(http(RCSB_SEARCH, q))
            ids += page.get("result_set", [])
            if len(ids) >= page.get("total_count", 0) or not page.get("result_set"):
                break
            start = len(ids)
        return json.dumps({"total_count": len(ids), "result_set": ids}).encode()

    return json.loads(cached(RAW / "rcsb" / f"search_{acc}.json", go)).get("result_set", [])


ENTITY_QUERY = """{polymer_entities(entity_ids:%s){
 rcsb_id
 entity_poly{rcsb_sample_sequence_length}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids}
 rcsb_polymer_entity_align{reference_database_accession aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}
 entry{rcsb_id
  rcsb_entry_info{resolution_combined experimental_method}
  nonpolymer_entities{rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}
   nonpolymer_comp{chem_comp{id name formula_weight}}}}}}"""


def fetch_entities(acc, entity_ids):
    def go():
        out = []
        for i in range(0, len(entity_ids), 25):
            q = ENTITY_QUERY % json.dumps(entity_ids[i : i + 25])
            out += json.loads(http(RCSB_GRAPHQL, {"query": q}))["data"]["polymer_entities"]
        return json.dumps(out).encode()

    return json.loads(cached(RAW / "rcsb" / f"entities_{acc}.json", go))


def build_inventory(entities, acc, domain):
    """One record per polymer entity: what it is, what is bound, how much it covers."""
    dom = set(range(domain["start"], domain["end"] + 1)) if domain else set()
    rows = []
    for e in entities:
        entry = e["entry"]
        info = entry["rcsb_entry_info"]
        res = (info.get("resolution_combined") or [None])[0]
        ligs = []
        for n in entry.get("nonpolymer_entities") or []:
            c = n["nonpolymer_comp"]["chem_comp"]
            ligs.append(
                {
                    "comp_id": c["id"],
                    "name": c["name"],
                    "formula_weight": c.get("formula_weight"),
                    "auth_asym_ids": n["rcsb_nonpolymer_entity_container_identifiers"]["auth_asym_ids"],
                    "drug_like": c["id"] not in NOT_LIGANDS
                    and (c.get("formula_weight") or 0) >= MIN_LIGAND_MW,
                }
            )
        covered, regions = set(), []
        for a in e.get("rcsb_polymer_entity_align") or []:
            if a.get("reference_database_accession") != acc:
                continue
            for r in a["aligned_regions"]:
                beg, ln = r["ref_beg_seq_id"], r["length"]
                regions.append([beg, beg + ln - 1])
                covered |= set(range(beg, beg + ln))
        rows.append(
            {
                "entity_id": e["rcsb_id"],
                "pdb_id": entry["rcsb_id"],
                "method": info["experimental_method"],
                "resolution_a": res,
                "auth_chain_ids": e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"],
                "entity_length": e["entity_poly"]["rcsb_sample_sequence_length"],
                "uniprot_regions": regions,
                "construct_domain_coverage": round(len(covered & dom) / len(dom), 4) if dom else None,
                "ligands": ligs,
                "drug_like_ligands": [l["comp_id"] for l in ligs if l["drug_like"]],
                "apo_or_holo": "holo" if any(l["drug_like"] for l in ligs) else "apo",
            }
        )
    return rows


def eligible_holo(inventory, bar):
    """X-ray entries with a drug-like ligand bound and enough of the catalytic domain,
    best resolution first."""
    elig = [
        r for r in inventory
        if r["apo_or_holo"] == "holo" and "X-ray" in r["method"] and r["resolution_a"]
        and (r["construct_domain_coverage"] or 0) >= bar
    ]
    return sorted(elig, key=lambda r: (r["resolution_a"], -(r["construct_domain_coverage"] or 0), r["pdb_id"]))


def rank_candidates(inventory, min_cov=0.95):
    """Stated pick rule: X-ray + a drug-like ligand + a complete catalytic domain,
    then best resolution. The coverage bar relaxes only if nothing clears it."""
    for bar in (min_cov, 0.85, 0.0):
        elig = eligible_holo(inventory, bar)
        if elig:
            return elig, bar
    return [], 0.0


# --------------------------------------------------------------------------
# structure prep
# --------------------------------------------------------------------------
def download_cif(pdb_id):
    path = RAW / "pdb" / f"{pdb_id}.cif"
    cached(path, lambda: http(f"{RCSB_FILES}/{pdb_id}.cif"))
    return path


def load_chain(pdb_id, chain_id):
    s = MMCIFParser(QUIET=True).get_structure(pdb_id, str(download_cif(pdb_id)))
    model = s[0]
    for ch in list(model):
        if ch.id != chain_id:
            model.detach_child(ch.id)
    if chain_id not in model:
        raise SystemExit(f"chain {chain_id} absent from {pdb_id}")
    return s


def write_prepared(structure, chain_id, ligand_code, holo_path, stripped_path):
    """holo = chain + the one drug-like ligand; stripped = chain alone. Modified residues
    are promoted to ATOM records so P2Rank sees an unbroken chain."""
    for model in list(structure):
        if model.id != 0:
            structure.detach_child(model.id)  # PDBIO writes every model; load_chain prunes only model 0
    chain = structure[0][chain_id]
    keep, promoted = [], []
    for r in list(chain):
        if r.id[0] == " ":
            continue
        if is_polymer(r):
            chain.detach_child(r.id)
            r.id = (" ", r.id[1], r.id[2])
            chain.add(r)
            promoted.append(r.get_resname().strip())
        elif r.get_resname().strip() == ligand_code and not keep:
            keep.append(r.id)
        else:
            chain.detach_child(r.id)
    chain.child_list.sort(key=lambda r: (r.id[1], r.id[2]))
    io = PDBIO()
    io.set_structure(structure)
    holo_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(holo_path))
    for rid in keep:
        chain.detach_child(rid)
    io.save(str(stripped_path))
    n_het = sum(1 for line in stripped_path.read_text().splitlines() if line.startswith("HETATM"))
    if n_het:
        raise SystemExit(f"ligand stripping failed: {n_het} HETATM left in {stripped_path}")
    return len(keep)


# --------------------------------------------------------------------------
# P2Rank
# --------------------------------------------------------------------------
def run_p2rank(pdb_path, outdir):
    csv_path = outdir / f"{pdb_path.name}_predictions.csv"
    if csv_path.exists():
        return csv_path
    env = dict(os.environ, JAVA_HOME=str(JAVA_HOME), PATH=f"{JAVA_HOME / 'bin'}:{os.environ['PATH']}")
    outdir.mkdir(parents=True, exist_ok=True)
    for attempt in range(2):  # the JVM start-up fails intermittently; a silent skip would shrink n
        proc = subprocess.run(
            [str(PRANK), "predict", "-f", str(pdb_path), "-o", str(outdir)],
            env=env, capture_output=True, text=True,
        )
        if proc.returncode == 0 and csv_path.exists():
            return csv_path
        print(f"  P2Rank rc={proc.returncode} on {pdb_path.name}, attempt {attempt + 1}", file=sys.stderr)
        time.sleep(2)
    raise SystemExit(f"P2Rank failed (rc={proc.returncode}):\n{proc.stdout[-3000:]}\n{proc.stderr[-2000:]}")


def parse_p2rank(csv_path):
    pockets = []
    with csv_path.open() as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = {k.strip(): v.strip() for k, v in row.items()}
            resids = []
            for tok in row["residue_ids"].split():
                ch, _, num = tok.partition("_")  # "A_919", or "A_100A" with an insertion code
                resids.append({
                    "chain": ch, "label": num,
                    "residue_id": int("".join(c for c in num if c.isdigit() or c == "-")),
                })
            pockets.append(
                {
                    "name": row["name"],
                    "rank": int(row["rank"]),
                    "score": float(row["score"]),
                    "probability": float(row["probability"]),
                    "n_sas_points": int(row["sas_points"]),
                    "n_surf_atoms": int(row["surf_atoms"]),
                    "center": [float(row["center_x"]), float(row["center_y"]), float(row["center_z"])],
                    "residues": resids,
                    "surf_atom_serials": [int(x) for x in row["surf_atom_ids"].split()],
                }
            )
    return pockets


# --------------------------------------------------------------------------
# geometry (task B10) - grid based, definitions recorded in the output json
# --------------------------------------------------------------------------
GEOMETRY_METHOD = {
    "grid_spacing_a": 0.6,
    "grid_bounds": "bounding box of the pocket-lining heavy atoms padded by 3 A; a cavity lobe "
    "reaching past that is not measured, so volume and extent are lower bounds for tunnel-shaped sites",
    "cavity_point": "grid point free of protein (dist to every heavy atom > vdW + 1.4 A probe), "
    "within 4.0 A of a pocket-lining heavy atom, and buriedness >= 0.75. The 0.75 bar is what "
    "separates the enclosed cavity from the thin solvent film that hugs the whole cleft surface; "
    "these are our definitions, so the volume is not comparable to an fpocket volume.",
    "buriedness": "fraction of 30 Fibonacci-sphere ray directions that strike a protein heavy atom "
    "within 12 A (sampled every 1.0 A)",
    "component": "cavity points are grouped into grid-connected components and only the component "
    "holding the P2Rank pocket centre is kept, so a surface groove elsewhere on the protein cannot "
    "inflate the volume or the docking box",
    "volume_a3": "n_cavity_points * spacing^3",
    "enclosure": "mean buriedness over cavity points (0-1)",
    "concavity": "fraction of cavity points lying inside the protein's convex hull (0-1)",
    "depth_a": "distance from a cavity point to the nearest convex-hull facet, mean and max",
    "max_sphere_radius_a": "max over cavity points of the distance to the nearest heavy-atom surface",
}


def fibonacci_directions(n=30):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def heavy_atoms(structure, chain_id):
    coords, radii = [], []
    for r in structure[0][chain_id]:
        if not is_polymer(r):
            continue
        for a in r:
            el = a.element.strip().upper()
            if el in ("H", "D", ""):
                continue
            coords.append(a.coord)
            radii.append(VDW.get(el, 1.7))
    return np.array(coords, dtype=float), np.array(radii, dtype=float)


def pocket_geometry(structure, chain_id, pocket, spacing=0.6, lining_dist=4.0, bur_min=0.75):
    coords, radii = heavy_atoms(structure, chain_id)
    tree = cKDTree(coords)
    lining = {(r["chain"], r["label"]) for r in pocket["residues"]}
    line_xyz = np.array(
        [a.coord for r in structure[0][chain_id]
         if is_polymer(r) and (chain_id, res_label(r)) in lining for a in r],
        dtype=float,
    )
    if len(line_xyz) == 0:
        return None
    lo, hi = line_xyz.min(0) - 3.0, line_xyz.max(0) + 3.0
    axes = [np.arange(lo[i], hi[i] + spacing, spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)

    # free of protein: no heavy atom closer than its vdW radius + probe
    d, idx = tree.query(grid, k=8)
    free = (d - radii[idx] - PROBE > 0).all(axis=1)
    near = cKDTree(line_xyz).query(grid, k=1)[0] <= lining_dist
    cand = grid[free & near]
    if len(cand) == 0:
        return None

    dirs = fibonacci_directions(30)
    steps = np.arange(1.0, 12.01, 1.0)
    hits = np.zeros(len(cand))
    for v in dirs:  # one direction at a time keeps the sample array small
        pts = cand[:, None, :] + v[None, None, :] * steps[None, :, None]
        dd, ii = tree.query(pts.reshape(-1, 3), k=1)
        blocked = (dd <= radii[ii]).reshape(len(cand), len(steps)).any(axis=1)
        hits += blocked
    buriedness = hits / len(dirs)
    keep = buriedness >= bur_min
    cav, bur = cand[keep], buriedness[keep]
    if len(cav) == 0:
        return None

    # a grid point 5 A from a lining atom can sit in an unrelated groove; keep only the
    # connected blob that actually contains the pocket P2Rank pointed at
    ctree = cKDTree(cav)
    pairs = np.array(sorted(ctree.query_pairs(spacing * 1.8)), dtype=int).reshape(-1, 2)
    n_comp, label = connected_components(
        coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(cav), len(cav))),
        directed=False,
    )
    anchor = label[ctree.query(np.array(pocket["center"], dtype=float), k=1)[1]]
    sel = label == anchor
    cav, bur = cav[sel], bur[sel]

    dc, ic = tree.query(cav, k=8)
    clearance = (dc - radii[ic]).min(axis=1)
    hull = ConvexHull(coords)
    signed = hull.equations[:, :3] @ cav.T + hull.equations[:, 3:4]  # <0 inside, |.| = distance
    depth = -signed.max(axis=0)
    extent = cav.max(0) - cav.min(0)
    return {
        "n_cavity_points": int(len(cav)),
        "n_cavity_components": int(n_comp),
        "volume_a3": round(float(len(cav) * spacing**3), 2),
        "enclosure": round(float(bur.mean()), 4),
        "concavity": round(float((depth > 0).mean()), 4),
        "mean_depth_a": round(float(depth.mean()), 3),
        "max_depth_a": round(float(depth.max()), 3),
        "max_sphere_radius_a": round(float(clearance.max()), 3),
        "cavity_centroid": [round(float(x), 3) for x in cav.mean(0)],
        "cavity_bbox_center": [round(float(x), 3) for x in (cav.min(0) + cav.max(0)) / 2],
        "cavity_extent_a": [round(float(x), 2) for x in extent],
    }


def docking_box(pocket, geom, padding=4.0, min_size=16.0, max_size=30.0):
    if geom:
        # bbox midpoint, not the centroid: for an asymmetric cavity the centroid is off-centre
        # and a box sized from the bbox but centred on it can miss part of the cavity
        center = geom["cavity_bbox_center"]
        size = [float(np.clip(e + 2 * padding, min_size, max_size)) for e in geom["cavity_extent_a"]]
        src = "cavity bounding box + 4 A padding, clipped to [16, 30] A"
    else:
        center, size, src = pocket["center"], [20.0, 20.0, 20.0], "P2Rank centre, default 20 A cube"
    return {"center": [round(float(c), 3) for c in center], "size_a": [round(s, 2) for s in size], "derived_from": src}


# --------------------------------------------------------------------------
# known-ligand contacts (validation only) and recovery
# --------------------------------------------------------------------------
def ligand_xyz(structure, chain_id, ligand_code):
    lig = [r for r in structure[0][chain_id] if r.id[0] != " " and r.get_resname().strip() == ligand_code]
    if not lig:
        return None
    return np.array([a.coord for a in lig[0] if a.element.strip().upper() not in ("H", "D", "")], dtype=float)


def ligand_contacts(structure, chain_id, ligand_code, cutoff=CONTACT_CUTOFF):
    chain = structure[0][chain_id]
    lig_xyz = ligand_xyz(structure, chain_id, ligand_code)
    if lig_xyz is None:
        return None
    tree = cKDTree(lig_xyz)
    out = []
    for r in chain:
        if not is_polymer(r):
            continue
        xyz = np.array([a.coord for a in r if a.element.strip().upper() not in ("H", "D", "")], dtype=float)
        if len(xyz) == 0:
            continue
        d = tree.query(xyz, k=1)[0].min()
        if d <= cutoff:
            out.append({"chain": chain_id, "residue_id": r.id[1], "label": res_label(r),
                        "residue_name": r.get_resname(), "min_distance_a": round(float(d), 2)})
    return {"ligand_comp_id": ligand_code, "n_ligand_heavy_atoms": int(len(lig_xyz)),
            "cutoff_a": cutoff, "contacts": sorted(out, key=lambda x: x["residue_id"])}


METALS = {"MG", "ZN", "MN", "CA", "FE", "K", "NA", "CU", "NI", "CO", "CD", "HG", "CS", "SR", "BA"}


def environment_flags(pdb_id, chain_id, pocket, ligand_code, cutoff=5.0):
    """What else occupies this site in the original crystal - metals, ordered waters, other
    het groups. Read from the unmodified mmCIF (all chains), not from the stripped file."""
    full = MMCIFParser(QUIET=True).get_structure(pdb_id, str(download_cif(pdb_id)))[0]
    lining = {r["label"] for r in pocket["residues"]}
    site_xyz = np.array(
        [a.coord for r in full[chain_id] if is_polymer(r) and res_label(r) in lining for a in r], dtype=float
    )
    tree = cKDTree(site_xyz)
    metals, waters, others, additives = set(), 0, set(), set()
    for ch in full:
        for r in ch:
            if is_polymer(r):
                continue
            xyz = np.array([a.coord for a in r], dtype=float)
            if len(xyz) == 0 or tree.query(xyz, k=1)[0].min() > cutoff:
                continue
            name = r.get_resname().strip()
            if name in ("HOH", "WAT", "DOD"):
                waters += 1
            elif len(r) == 1 and (list(r)[0].element or "").strip().upper() in METALS:
                metals.add(name)
            elif name != ligand_code and name not in NOT_LIGANDS:
                others.add(name)
            elif name != ligand_code:
                additives.add(name)
    return {
        "has_metal": bool(metals), "metal_comp_ids": sorted(metals),
        "has_ordered_water": waters > 0, "n_ordered_waters_in_site": waters,
        "has_cofactor": bool(others), "other_het_comp_ids": sorted(others),
        "additive_comp_ids": sorted(additives),
        "cutoff_a": cutoff,
    }


def n_observed_residues(structure, chain_id):
    return sum(1 for r in structure[0][chain_id] if is_polymer(r))


def recovery(pocket, contacts, n_residues):
    """Overlap between a P2Rank pocket and the known-ligand contact residues,
    with an exact hypergeometric p-value against a random residue set of the
    same size drawn from the observed chain."""
    P = {r["label"] for r in pocket["residues"]}   # label, not int: an insertion code is a distinct residue
    K = {c["label"] for c in contacts["contacts"]}
    inter = P & K
    p = float(hypergeom.sf(len(inter) - 1, n_residues, len(K), len(P))) if inter else 1.0
    return {
        "pocket_rank": pocket["rank"],
        "n_pocket_residues": len(P),
        "n_known_contact_residues": len(K),
        "n_overlap": len(inter),
        "jaccard": round(len(inter) / len(P | K), 4),
        "recall_of_known_contacts": round(len(inter) / len(K), 4),
        "precision": round(len(inter) / len(P), 4),
        "n_observed_residues_in_chain": n_residues,
        "expected_overlap_random": round(len(P) * len(K) / n_residues, 3),
        "hypergeometric_p": float(f"{p:.3e}"),
        "overlap_residue_labels": sorted(inter),
        "missed_known_contact_labels": sorted(K - P),
    }


# --------------------------------------------------------------------------
# per-structure analysis
# --------------------------------------------------------------------------
def pick_chain_and_ligand(rec):
    """First chain that actually has a drug-like ligand modelled in it - picking the chain
    and the ligand independently drops structures where the ligand sits only in chain B."""
    drug_like = [l for l in rec["ligands"] if l["drug_like"]]
    for chain in rec["auth_chain_ids"]:
        here = [l for l in drug_like if chain in l["auth_asym_ids"]]
        if here:
            return chain, max(here, key=lambda l: l["formula_weight"] or 0)["comp_id"]
    return rec["auth_chain_ids"][0], drug_like[0]["comp_id"]


def analyse(rec, workdir, with_geometry=False, n_pockets=5):
    pdb_id = rec["pdb_id"]
    chain_id, ligand = pick_chain_and_ligand(rec)
    struct = load_chain(pdb_id, chain_id)
    holo = workdir / f"{pdb_id}_{chain_id}_holo.pdb"
    stripped = workdir / f"{pdb_id}_{chain_id}_stripped.pdb"
    n_kept = write_prepared(struct, chain_id, ligand, holo, stripped)
    if n_kept == 0:
        return None  # ligand not resolved in this chain
    holo_struct = load_chain(pdb_id, chain_id)  # reload: write_prepared mutates in place
    contacts = ligand_contacts(holo_struct, chain_id, ligand)
    pockets = parse_p2rank(run_p2rank(stripped, RAW / "p2rank" / f"{pdb_id}_{chain_id}"))
    if not pockets or not contacts or not contacts["contacts"]:
        return None
    stripped_struct = load_chain(pdb_id, chain_id)
    for r in list(stripped_struct[0][chain_id]):
        if not is_polymer(r):
            stripped_struct[0][chain_id].detach_child(r.id)
    n_res = n_observed_residues(stripped_struct, chain_id)
    rows = [recovery(p, contacts, n_res) for p in pockets[:n_pockets]]
    geom = pocket_geometry(stripped_struct, chain_id, pockets[0]) if with_geometry else None
    box = docking_box(pockets[0], geom) if with_geometry else None
    box_check = None
    if box:  # validation only: does the box built without the ligand actually cover it?
        lx = ligand_xyz(holo_struct, chain_id, ligand)
        c, sz = np.array(box["center"]), np.array(box["size_a"])
        box_check = {
            "n_ligand_heavy_atoms": int(len(lx)),
            "ligand_extent_a": [round(float(x), 2) for x in (lx.max(0) - lx.min(0))],
            "fraction_inside_docking_box": round(float((np.abs(lx - c) <= sz / 2).all(axis=1).mean()), 4),
        }
    return {
        "docking_box": box, "docking_box_check": box_check,
        "environment": environment_flags(pdb_id, chain_id, pockets[0], ligand),
        "pdb_id": pdb_id, "chain_id": chain_id, "ligand": ligand,
        "resolution_a": rec["resolution_a"], "n_observed_residues": n_res,
        "structure_paths": {"holo": str(holo.relative_to(ROOT)), "stripped": str(stripped.relative_to(ROOT))},
        "pockets": pockets, "geometry": geom, "contacts": contacts,
        "recovery_top1": rows[0],
        "recovery_by_pocket": rows,
        "best_pocket_by_jaccard": max(rows, key=lambda r: (r["jaccard"], -r["pocket_rank"])),
    }


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------
def write_structures_json(out, up, dom, inventory, primary, res, rationale):
    """The whole inventory, not just the winner."""
    (out / "structures.json").write_text(json.dumps({
        "uniprot": {k: up[k] for k in ("accession", "entry_name", "gene", "protein_name", "length")},
        "catalytic_domain": dom,
        "sequence_files": {"full": "sequence.fasta", "construct": "sequence_kinase_domain.fasta" if dom else None},
        "ligand_filter": {"min_formula_weight": MIN_LIGAND_MW, "excluded_comp_ids": sorted(NOT_LIGANDS)},
        "n_entities": len(inventory),
        "n_holo": sum(1 for r in inventory if r["apo_or_holo"] == "holo"),
        "primary": {"pdb_id": res["pdb_id"], "chain_id": res["chain_id"],
                    "ligand": res["ligand"], "rationale": rationale,
                    "structure_paths": res["structure_paths"],
                    "n_observed_residues": res["n_observed_residues"]},
        "inventory": sorted(inventory, key=lambda r: (r["resolution_a"] is None, r["resolution_a"] or 0)),
    }, indent=2) + "\n")

def write_pocket_json(out, up, dom, res, n_pockets):
    """SiteSpec-shaped (PROJECT_GOAL.md 4.2), geometry included."""
    top = res["pockets"][0]
    (out / "pocket.json").write_text(json.dumps({
        "site_id": f"{res['pdb_id']}_{res['chain_id']}_pocket1",
        "uniprot_id": up["accession"],
        "structure_source": "pdb",
        "structure_id": res["pdb_id"],
        "structure_path": res["structure_paths"]["stripped"],
        "apo_or_holo": "stripped_holo",
        "chain_id": res["chain_id"],
        "construct_range": [dom["start"], dom["end"]] if dom else None,
        "site_kind": "pocket",
        "detector": "P2Rank 2.5 (default model), run on the ligand-stripped structure",
        "residue_ids": sorted(r["residue_id"] for r in top["residues"]),
        "center": [round(c, 3) for c in top["center"]],
        "p2rank_score": top["score"],
        "p2rank_probability": top["probability"],
        "environment": res["environment"],
        "geometry_method": GEOMETRY_METHOD,
        "geometry": res["geometry"],
        "docking_box": res["docking_box"],
        "pockets": [{
            "rank": p["rank"], "score": p["score"], "probability": p["probability"],
            "center": [round(c, 3) for c in p["center"]],
            "n_surf_atoms": p["n_surf_atoms"],
            "residue_ids": sorted(r["residue_id"] for r in p["residues"]),
        } for p in res["pockets"][: n_pockets]],
        "known_ligand_residue_ids": "NOT SET HERE BY DESIGN - see known_ligand_contacts.json; "
                                   "site detection must not read it (PROJECT_GOAL.md B12)",
    }, indent=2) + "\n")

def write_known_ligand_json(out, res):
    """Validation only - see the warning it carries."""
    (out / "known_ligand_contacts.json").write_text(json.dumps({
        "_WARNING": "VALIDATION USE ONLY. These residues come from the bound ligand in the holo "
                    "structure. Site detection (P2Rank) never saw them - it ran on the stripped "
                    "structure - and no signature-construction code may read this file "
                    "(PROJECT_GOAL.md task B12). Reading it into WS-C invalidates the "
                    "hotspot-recovery benchmark.",
        "structure_id": res["pdb_id"], "chain_id": res["chain_id"],
        "method": f"protein heavy atom within {CONTACT_CUTOFF} A of a ligand heavy atom",
        **res["contacts"],
        "known_ligand_residue_ids": sorted(c["residue_id"] for c in res["contacts"]["contacts"]),
        "docking_box_check": res["docking_box_check"],
    }, indent=2) + "\n")

def write_recovery_json(out, runs, skipped, jac, rc, best, hit1):
    """The measured headline."""
    (out / "pocket_recovery.json").write_text(json.dumps({
        "question": "Does P2Rank, run blind on the ligand-stripped structure, recover the residues "
                    "that the known ligand actually contacts?",
        "metric": "Jaccard / recall between the P2Rank top-1 pocket residue set and the known-ligand "
                  f"contact residue set ({CONTACT_CUTOFF} A heavy-atom cutoff)",
        "null_model": "random residue subset of the same size drawn from the observed chain; exact "
                      "hypergeometric survival function",
        "n_structures": len(runs),
        "n_structures_skipped": len(skipped),
        "skipped": skipped,
        "top1_jaccard_mean": round(float(jac.mean()), 4),
        "top1_jaccard_sd": round(float(jac.std(ddof=1)), 4) if len(jac) > 1 else None,
        "top1_jaccard_median": round(float(np.median(jac)), 4),
        "top1_jaccard_min_max": [round(float(jac.min()), 4), round(float(jac.max()), 4)],
        "top1_recall_mean": round(float(rc.mean()), 4),
        "top1_recall_sd": round(float(rc.std(ddof=1)), 4) if len(rc) > 1 else None,
        "top1_recall_median": round(float(np.median(rc)), 4),
        "top1_recall_min_max": [round(float(rc.min()), 4), round(float(rc.max()), 4)],
        "best_pocket_jaccard_mean": round(float(best.mean()), 4),
        "n_structures_where_best_pocket_is_rank1": hit1,
        "structures_where_top1_is_not_the_best_pocket": [
            {"pdb_id": r["pdb_id"], "ligand": r["ligand"], "top1_jaccard": r["recovery_top1"]["jaccard"],
             "best_pocket_rank": r["best_pocket_by_jaccard"]["pocket_rank"],
             "best_pocket_jaccard": r["best_pocket_by_jaccard"]["jaccard"]}
            for r in runs if r["best_pocket_by_jaccard"]["pocket_rank"] != 1
        ],
        "max_hypergeometric_p_top1": max(r["recovery_top1"]["hypergeometric_p"] for r in runs),
        "per_structure": [{
            "pdb_id": r["pdb_id"], "chain_id": r["chain_id"], "ligand": r["ligand"],
            "resolution_a": r["resolution_a"], "n_pockets_found": len(r["pockets"]),
            "top1": r["recovery_top1"], "best_pocket": r["best_pocket_by_jaccard"],
            "by_pocket": r["recovery_by_pocket"],
        } for r in runs],
    }, indent=2) + "\n")


def resolve_slug(arg):
    """Write next to the sibling tasks' output: the default slug if it is already there,
    otherwise the only disease dir; refuse to guess between several."""
    if arg:
        return arg
    base = ROOT / "results" / "pipeline"
    dirs = sorted(d.name for d in base.glob("*") if d.is_dir()) if base.exists() else []
    if DEFAULT_SLUG in dirs or not dirs:
        return DEFAULT_SLUG
    if len(dirs) == 1:
        return dirs[0]
    raise SystemExit(f"several disease dirs exist ({', '.join(dirs)}); pass --disease to say which")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uniprot", default="P35968")
    ap.add_argument("--disease", default=None,
                    help=f"output slug under results/pipeline/ (default {DEFAULT_SLUG})")
    ap.add_argument("--pdb", default=None, help="override the primary structure pick")
    ap.add_argument("--top-pockets", type=int, default=5)
    ap.add_argument("--robustness", type=int, default=0,
                    help="how many further holo structures to repeat the recovery test on; "
                         "0 (default) means every eligible holo entry")
    args = ap.parse_args()

    slug = resolve_slug(args.disease)
    out = ROOT / "results" / "pipeline" / slug / "target"
    (out / "structure").mkdir(parents=True, exist_ok=True)

    up = fetch_uniprot(args.uniprot)
    dom = catalytic_domain(up)
    print(f"{up['accession']} {up['entry_name']} ({up['gene']}) {up['length']} aa", file=sys.stderr)
    print(f"  catalytic domain: {dom}", file=sys.stderr)

    (out / "sequence.fasta").write_text(
        f">sp|{up['accession']}|{up['entry_name']} {up['protein_name']} len={up['length']}\n"
        + "\n".join(up["sequence"][i : i + 60] for i in range(0, up["length"], 60)) + "\n"
    )
    if dom:
        sub = up["sequence"][dom["start"] - 1 : dom["end"]]
        (out / "sequence_kinase_domain.fasta").write_text(
            f">{up['accession']}|{dom['description']}|{dom['start']}-{dom['end']} len={len(sub)} "
            f"(truncated from the {up['length']} aa canonical sequence to the UniProt "
            f"'{dom['description']}' feature; the ligand site is inside it and the crystal constructs "
            f"cover this range, so downstream co-folding does not need the rest)\n"
            + "\n".join(sub[i : i + 60] for i in range(0, len(sub), 60)) + "\n"
        )

    ids = search_entities(args.uniprot)
    inventory = build_inventory(fetch_entities(args.uniprot, ids), args.uniprot, dom)
    ranked, bar = rank_candidates(inventory)
    print(f"  {len(inventory)} PDB polymer entities, {sum(1 for r in inventory if r['apo_or_holo'] == 'holo')} holo, "
          f"{len(ranked)} eligible at coverage >= {bar}", file=sys.stderr)
    if not ranked:
        raise SystemExit("no eligible holo X-ray structure")

    if args.pdb:
        chosen = next((r for r in eligible_holo(inventory, 0.0) if r["pdb_id"] == args.pdb.upper()), None)
        if chosen is None:
            raise SystemExit(f"{args.pdb} is not an eligible holo structure for {args.uniprot}")
        ranked = [chosen] + [r for r in ranked if r is not chosen]

    primary = ranked[0]
    work = out / "structure"
    res = analyse(primary, work, with_geometry=True, n_pockets=args.top_pockets)
    if res is None:
        raise SystemExit(f"primary {primary['pdb_id']} could not be analysed (ligand unresolved?)")

    cov = (
        f"construct covers {primary['construct_domain_coverage']:.1%} of the UniProt "
        f"{dom['description']} ({dom['start']}-{dom['end']})"
        if dom and primary["construct_domain_coverage"] is not None
        else "no catalytic domain annotated in UniProt, so coverage was not used"
    )
    rule = (
        f"Rule: among X-ray entries with >=1 drug-like ligand (MW >= {MIN_LIGAND_MW:g} Da, not a known "
        f"additive) and construct coverage >= {bar}, take the best resolution; ties broken by coverage "
        f"then PDB id."
        if not args.pdb
        else f"Rule: overridden on the command line with --pdb {args.pdb.upper()}; the automatic pick "
             f"would have been {ranked[1]['pdb_id'] if len(ranked) > 1 else 'none'}."
    )
    rationale = (
        f"{res['pdb_id']} chain {res['chain_id']}: X-ray {primary['resolution_a']} A, "
        f"drug-like ligand {res['ligand']} bound, {cov}. {rule}"
    )
    print(f"  primary: {rationale}", file=sys.stderr)

    # repeat the same measurement on further holo structures so the headline has an n
    # the recovery test repeats on every eligible holo entry, not on a resolution-ranked
    # slice of them: the low-resolution tail is where the type-II / DFG-out complexes sit
    # and dropping it would flatter P2Rank
    pool, seen = [], {primary["pdb_id"]}  # one record per ENTITY; a two-entity crystal must not count twice
    for r in eligible_holo(inventory, min(bar, 0.85)):
        if r["pdb_id"] not in seen:
            seen.add(r["pdb_id"])
            pool.append(r)
    if args.robustness:
        pool = pool[: args.robustness]
    extra, skipped = [], []
    for rec in pool:
        try:
            r = analyse(rec, RAW / "prepared", with_geometry=False, n_pockets=args.top_pockets)
        except (Exception, SystemExit) as exc:  # SystemExit is BaseException: the helpers raise it
            print(f"  skip {rec['pdb_id']}: {exc}", file=sys.stderr)
            skipped.append({"pdb_id": rec["pdb_id"], "reason": str(exc)[:200]})
            continue
        if r is None:
            print(f"  skip {rec['pdb_id']}: ligand not resolved in chain", file=sys.stderr)
            skipped.append({"pdb_id": rec["pdb_id"], "reason": "drug-like ligand not modelled in the chain"})
            continue
        extra.append(r)
        print(f"  {r['pdb_id']} {r['ligand']}: top1 jaccard {r['recovery_top1']['jaccard']:.3f} "
              f"recall {r['recovery_top1']['recall_of_known_contacts']:.3f}", file=sys.stderr)

    runs = [res] + extra
    jac = np.array([r["recovery_top1"]["jaccard"] for r in runs])
    rc = np.array([r["recovery_top1"]["recall_of_known_contacts"] for r in runs])
    best = np.array([r["best_pocket_by_jaccard"]["jaccard"] for r in runs])
    hit1 = sum(1 for r in runs
               if r["best_pocket_by_jaccard"]["pocket_rank"] == 1 and r["best_pocket_by_jaccard"]["jaccard"] > 0)

    write_structures_json(out, up, dom, inventory, primary, res, rationale)
    write_pocket_json(out, up, dom, res, args.top_pockets)
    write_known_ligand_json(out, res)
    write_recovery_json(out, runs, skipped, jac, rc, best, hit1)

    t = res["recovery_top1"]
    print(f"\nP2Rank top-1 vs known ligand contacts ({res['pdb_id']}/{res['ligand']}): "
          f"overlap {t['n_overlap']}/{t['n_known_contact_residues']} known, jaccard {t['jaccard']:.3f}, "
          f"recall {t['recall_of_known_contacts']:.3f}, p={t['hypergeometric_p']:.2e}", file=sys.stderr)
    print(f"Across n={len(runs)} holo structures: jaccard {jac.mean():.3f}"
          + (f" +/- {jac.std(ddof=1):.3f}" if len(runs) > 1 else "")
          + f", recall {rc.mean():.3f}"
          + (f" +/- {rc.std(ddof=1):.3f}" if len(runs) > 1 else "")
          + f"; best pocket was rank 1 in {hit1}/{len(runs)}", file=sys.stderr)
    print(f"Wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
