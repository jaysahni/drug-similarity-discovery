"""Measure and partly fill the approved-drug hole in `data/approved_drugs.csv`.

WHY
---
`data/approved_drugs.csv` is 4,099 DrugCentral structures, 2,153 with
`approved == 1`. Every enrichment factor, base rate and positive-control list in
this project is computed against it, so a drug the corpus is missing is a
ceiling on what any method can be shown to do -- and the miss is silent.

The motivating case is thrombin (gene F2, UniProt P00734). The corpus carries 12
approved drugs with an F2 annotation and 4 with F2 as mechanism, and it is
missing the free-acid `dabigatran` (only the `dabigatran etexilate` prodrug is
there), `melagatran`, `desirudin`, `lepirudin`, `hirudin` and `heparin`.
Annotation enrichment cannot fix that: ChEMBL holds thousands of F2 activities
at pChEMBL >= 6 but they belong to investigational compounds. Only adding drugs
fixes it.

WHAT THIS DOES
--------------
1. Diffs the corpus against ChEMBL `max_phase = 4` -- every molecule ChEMBL
   records as approved by some agency, currently or historically -- matched on
   the 14-character InChIKey skeleton after the SAME standardisation the corpus
   was built with (`scripts/build_drug_corpus.py:standardise`, imported, not
   reimplemented). That measures the hole.
2. Writes the absent ones to `data/approved_drugs_extended.csv`, schema-
   compatible with `approved_drugs.csv` so the two concatenate, plus the
   provenance columns `source` / `source_url` / `withdrawn` and the ChEMBL
   fields the rows came from.
3. Adds named active moieties of approved prodrugs (dabigatran) that ChEMBL
   itself files below phase 4, with `approval_basis` recording why.
4. Carries ChEMBL mechanism-of-action annotations (gene symbols, the same
   vocabulary DrugCentral uses) so added rows work as positive controls.

WHAT IT DOES NOT DO
-------------------
It does not touch `data/approved_drugs.csv` or `data/drugs.csv`. Both are
load-bearing for already-published numbers. The additions are a separate file.

Peptides, proteins and polysaccharides with no single SMILES (lepirudin,
desirudin, heparin) are NOT forced into this file; they are reported in an
explicit "not added / why" table and, where they belong there,
`data/approved_biologics.csv` is checked for them.

STRUCTURE PROVENANCE
--------------------
No SMILES in this file was typed from memory. Every one comes from a fetched
ChEMBL record whose `source_url` is on the row, then canonicalised locally.

Usage:
    ./env-kit/bin/python demo/expand_corpus.py            # fetch (cached) + build
    ./env-kit/bin/python demo/expand_corpus.py --offline  # cache only, no network
    ./env-kit/bin/python demo/expand_corpus.py --selftest # build + assertions
    ./env-kit/bin/python demo/expand_corpus.py --spotcheck  # + PubChem round-trip
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CORPUS = DATA / "approved_drugs.csv"
BIOLOGICS = DATA / "approved_biologics.csv"
OUT = DATA / "approved_drugs_extended.csv"
CACHE = DATA / "raw" / "corpus_expansion"

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_WEB = "https://www.ebi.ac.uk/chembl/compound_report_card"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# struct_id offset. The DrugCentral corpus runs 4..5462; DrugCentral's own table
# is nowhere near 900000, so this cannot collide now or after a DrugCentral
# refresh. Asserted at the end of build().
ID_OFFSET = 900_000

# Schema of data/approved_drugs.csv, in order. The extended file starts with
# these twelve so `pd.concat` lines the two up, then adds provenance columns.
BASE_FIELDS = [
    "struct_id", "name", "smiles", "inchikey", "approved", "approval_agencies",
    "first_approval", "chembl_max_phase4", "n_targets", "targets", "moa_targets",
    "flags",
]
EXTRA_FIELDS = [
    "source", "source_url", "withdrawn", "chembl_id", "molecule_type",
    "chembl_max_phase", "approval_basis", "target_evidence", "notes",
]
FIELDS = BASE_FIELDS + EXTRA_FIELDS

# Molecules ChEMBL files below max_phase 4 that are nonetheless the circulating
# active species of an approved product. Each needs an explicit, checkable
# reason -- this list is the only place a human judgement enters the file, and
# `approval_basis` carries that judgement into every row that uses it.
ACTIVE_MOIETIES = {
    "CHEMBL48361": (
        "DABIGATRAN",
        "active moiety of dabigatran etexilate (Pradaxa, FDA 2010); ChEMBL files "
        "the free acid at max_phase 3 because the marketed entity is the prodrug. "
        "This is the species that actually binds F2.",
    ),
}

# The six drugs the task named, checked one by one at the end.
NAMED = ["DABIGATRAN", "MELAGATRAN", "LEPIRUDIN", "DESIRUDIN", "HIRUDIN", "HEPARIN"]


# ---------------------------------------------------------------------------
# reuse the corpus's own standardisation -- not a reimplementation of it
# ---------------------------------------------------------------------------
def _load_standardise():
    """Import `standardise` from scripts/build_drug_corpus.py (not a package)."""
    path = ROOT / "scripts" / "build_drug_corpus.py"
    spec = importlib.util.spec_from_file_location("build_drug_corpus", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.standardise


standardise = _load_standardise()


# ---------------------------------------------------------------------------
# cached HTTP
# ---------------------------------------------------------------------------
SESSION = requests.Session()


def _cached_get(name, url, params, offline, allow_404=False):
    """GET -> JSON, cached to disk by `name`. A re-run costs nothing.

    With `allow_404`, a 404 is a legitimate answer ("this database has no such
    record") and is cached as null, so re-runs do not re-ask.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    fp = CACHE / f"{name}.json"
    if fp.exists():
        return json.loads(fp.read_text())
    if offline:
        sys.exit(f"--offline but no cache at {fp}")
    for attempt in range(4):
        try:
            r = SESSION.get(url, params=params, timeout=120)
            if allow_404 and r.status_code == 404:
                payload = None
                break
            r.raise_for_status()
            payload = r.json()
            break
        except Exception as exc:  # noqa: BLE001 - retry any transport/5xx error
            if attempt == 3:
                raise
            print(f"    retry {attempt + 1} after {exc.__class__.__name__}")
            time.sleep(2 * (attempt + 1))
    fp.write_text(json.dumps(payload))
    return payload


MOL_ONLY = (
    "molecule_chembl_id,pref_name,max_phase,first_approval,withdrawn_flag,"
    "molecule_type,molecule_structures__canonical_smiles,"
    "molecule_structures__standard_inchi_key"
)


def fetch_phase4(offline):
    """Every ChEMBL molecule at max_phase = 4, structure or not.

    Structureless ones are kept deliberately: they are the biologic/ mixture
    part of the hole and have to be counted, not quietly skipped.
    """
    out, offset, limit, total = [], 0, 1000, None
    while True:
        page = _cached_get(
            f"phase4_{offset:06d}",
            f"{CHEMBL}/molecule.json",
            {"max_phase": 4, "limit": limit, "offset": offset, "only": MOL_ONLY},
            offline,
        )
        mols = page.get("molecules", [])
        total = page["page_meta"]["total_count"]
        if not mols:
            break
        out.extend(mols)
        offset += limit
        if offset >= total:
            break
    print(f"  ChEMBL max_phase=4 molecules fetched : {len(out)} (API total {total})")
    if len(out) != total:
        sys.exit(f"paged {len(out)} but API says {total} - incomplete fetch, refusing")
    return out


def fetch_by_name(name, offline):
    page = _cached_get(
        f"name_{name.lower().replace(' ', '_')}",
        f"{CHEMBL}/molecule.json",
        {"pref_name__iexact": name, "only": MOL_ONLY},
        offline,
    )
    return page.get("molecules", [])


def fetch_by_id(chembl_id, offline):
    page = _cached_get(
        f"id_{chembl_id}",
        f"{CHEMBL}/molecule.json",
        {"molecule_chembl_id": chembl_id, "only": MOL_ONLY},
        offline,
    )
    mols = page.get("molecules", [])
    return mols[0] if mols else None


def fetch_pubchem_by_name(name, offline):
    """PubChem name -> (CID, Title, isomeric SMILES) or None.

    Used only for ChEMBL molecules typed `Small molecule` that carry no
    structure in ChEMBL. `SMILES` is PubChem's isomeric SMILES;
    `ConnectivitySMILES` (what PubChem substitutes if you ask for the retired
    `CanonicalSMILES`) has the stereochemistry stripped and must not be used.
    """
    slug = "".join(c if c.isalnum() else "_" for c in name.lower())[:70]
    payload = _cached_get(
        f"pcname_{slug}",
        f"{PUBCHEM}/compound/name/{requests.utils.quote(name, safe='')}"
        f"/property/Title,SMILES/JSON",
        {}, offline, allow_404=True,
    )
    if not payload:
        return None
    props = (payload.get("PropertyTable") or {}).get("Properties") or []
    if not props or not props[0].get("SMILES"):
        return None
    return props[0]["CID"], props[0].get("Title", ""), props[0]["SMILES"]


# Element words that appear in drug names, and the symbol they promise. Used as
# an INDEPENDENT check on a name-resolved structure: if the name says
# technetium, the molecule must contain technetium.
#
# This is the check that matters, because PubChem's `/compound/name/` endpoint
# resolves THROUGH its synonym index, so asking afterwards whether the name is
# a synonym of the result is circular and always true. It is what caught
# 'TECHNETIUM MEDRONATE' resolving to CID 8094, heptanoic acid.
ELEMENT_WORDS = {
    "technetium": "Tc", "thall": "Tl", "ferric": "Fe", "ferrous": "Fe",
    "iron": "Fe", "cupric": "Cu", "cuprous": "Cu", "copper": "Cu",
    "strontium": "Sr", "gallium": "Ga", "indium": "In", "samarium": "Sm",
    "lutetium": "Lu", "yttrium": "Y", "radium": "Ra", "barium": "Ba",
    "bismuth": "Bi", "calcium": "Ca", "magnesium": "Mg", "sodium": "Na",
    "potassium": "K", "lithium": "Li", "zinc": "Zn", "selen": "Se",
    "chromi": "Cr", "mangan": "Mn", "cobalt": "Co", "platin": "Pt",
    "gadolin": "Gd", "alumin": "Al", "silver": "Ag", "auro": "Au",
    "gold": "Au", "mercur": "Hg", "arsen": "As", "antimon": "Sb",
    "lanthan": "La", "holmium": "Ho", "rhenium": "Re", "stann": "Sn",
    "titani": "Ti", "phosphorus": "P", "boron": "B", "nickel": "Ni",
    "molybd": "Mo", "vanad": "V", "cesium": "Cs", "rubidium": "Rb",
    "krypton": "Kr", "xenon": "Xe", "iodide": "I", "iodine": "I",
}


def element_promise_kept(name, raw_smiles):
    """(verdict, detail) - does a name-resolved structure contain what the name says?

    Checked on the RAW smiles, before salt stripping, because standardisation
    removes exactly the counter-ion the name is naming.
    """
    from rdkit import Chem

    low = name.lower()
    want = {sym for word, sym in ELEMENT_WORDS.items() if word in low}
    if not want:
        return True, ""
    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
        return False, "unparseable"
    have = {a.GetSymbol() for a in mol.GetAtoms()}
    missing = want - have
    if missing:
        return False, f"name promises {sorted(want)}, structure has none of {sorted(missing)}"
    return True, ""


def metal_stripped(raw_smiles, parent_smiles):
    """True if standardisation dropped a metal centre.

    `LargestFragmentChooser(preferOrganic=True)` keeps the largest ORGANIC
    fragment, so a coordination complex such as carboplatin loses its platinum
    and the surviving 'drug' is only the ligand. The existing corpus has the
    same behaviour; this flag makes those rows findable instead of silent.
    """
    from rdkit import Chem

    metals = set("Li Na K Rb Cs Be Mg Ca Sr Ba Al Ga In Tl Sn Pb Sb Bi Sc Ti V Cr "
                 "Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir "
                 "Pt Au Hg La Ce Gd Lu Sm Ho Y Ra Th U".split())
    a = Chem.MolFromSmiles(raw_smiles)
    b = Chem.MolFromSmiles(parent_smiles)
    if a is None or b is None:
        return False
    ma = {at.GetSymbol() for at in a.GetAtoms()} & metals
    mb = {at.GetSymbol() for at in b.GetAtoms()} & metals
    return bool(ma - mb)


def fetch_mechanisms(chembl_ids, offline):
    """molecule_chembl_id -> [(target_chembl_id, action_type, moa_text)]."""
    ids = sorted(chembl_ids)
    out = defaultdict(list)
    for i in range(0, len(ids), 40):
        chunk = ids[i : i + 40]
        offset = 0
        while True:
            page = _cached_get(
                f"mech_{i:05d}_{offset:04d}",
                f"{CHEMBL}/mechanism.json",
                {
                    "molecule_chembl_id__in": ",".join(chunk),
                    "limit": 1000,
                    "offset": offset,
                    "only": "molecule_chembl_id,target_chembl_id,action_type,"
                            "mechanism_of_action",
                },
                offline,
            )
            for m in page.get("mechanisms", []):
                if m.get("target_chembl_id"):
                    out[m["molecule_chembl_id"]].append(
                        (m["target_chembl_id"], m.get("action_type") or "",
                         m.get("mechanism_of_action") or "")
                    )
            total = page["page_meta"]["total_count"]
            offset += 1000
            if offset >= total:
                break
        if (i // 40) % 10 == 0:
            print(f"    mechanisms {i}/{len(ids)}")
    return out


def fetch_activities(chembl_ids, offline, pchembl_min=6.0, min_n=2):
    """(molecule, target) -> n measurements at pChEMBL >= `pchembl_min`, human.

    ChEMBL's curated `mechanism` table is thin for older or withdrawn drugs
    (melagatran has no mechanism row at all), so bioactivity is the fallback.
    Same evidence rule as `data/target_annotations.csv` in this repo
    (`evidence = chembl_pchembl`, pChEMBL >= 6), tightened with `min_n` so a
    single stray measurement does not become a target claim.
    """
    ids = sorted(chembl_ids)
    counts = defaultdict(Counter)
    for i in range(0, len(ids), 20):
        chunk = ids[i : i + 20]
        offset = 0
        while True:
            page = _cached_get(
                f"act_{i:05d}_{offset:05d}",
                f"{CHEMBL}/activity.json",
                {
                    "molecule_chembl_id__in": ",".join(chunk),
                    "pchembl_value__gte": pchembl_min,
                    "target_organism": "Homo sapiens",
                    "limit": 1000, "offset": offset,
                    "only": "molecule_chembl_id,target_chembl_id,pchembl_value",
                },
                offline,
            )
            for a in page.get("activities", []):
                if a.get("target_chembl_id"):
                    counts[a["molecule_chembl_id"]][a["target_chembl_id"]] += 1
            total = page["page_meta"]["total_count"]
            offset += 1000
            if offset >= total:
                break
        if (i // 20) % 5 == 0:
            print(f"    activities {i}/{len(ids)}")
    return {
        mol: {t for t, n in c.items() if n >= min_n}
        for mol, c in counts.items()
    }


def fetch_target_genes(target_ids, offline):
    """target_chembl_id -> sorted human gene symbols (DrugCentral's vocabulary).

    Non-human targets are dropped, matching `build_drug_corpus.load_targets`,
    which keeps `ORGANISM == "Homo sapiens"` rows only.
    """
    ids = sorted(target_ids)
    genes = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i : i + 40]
        page = _cached_get(
            f"target_{i:05d}",
            f"{CHEMBL}/target.json",
            {
                "target_chembl_id__in": ",".join(chunk),
                "limit": 100,
                "only": "target_chembl_id,organism,target_type,target_components",
            },
            offline,
        )
        for t in page.get("targets", []):
            if t.get("organism") != "Homo sapiens":
                genes[t["target_chembl_id"]] = []
                continue
            gs = set()
            for comp in t.get("target_components") or []:
                for syn in comp.get("target_component_synonyms") or []:
                    if syn.get("syn_type") == "GENE_SYMBOL":
                        gs.add(syn["component_synonym"])
            genes[t["target_chembl_id"]] = sorted(gs)
    return genes


# ---------------------------------------------------------------------------
def load_corpus():
    rows = list(csv.DictReader(CORPUS.open(newline="")))
    if list(rows[0]) != BASE_FIELDS:
        sys.exit(f"{CORPUS} schema changed: {list(rows[0])}")
    return rows


def report_existing_dupes(rows):
    """build_drug_corpus dedupes by struct_id, not InChIKey. Check the result."""
    by_key = defaultdict(list)
    for r in rows:
        if r["inchikey"]:
            by_key[r["inchikey"]].append((r["struct_id"], r["name"]))
    # A full-InChIKey duplicate across two DIFFERENT drug names is not a salt
    # form; it means one of the two rows carries the wrong structure. Print it
    # loudly. We do not fix approved_drugs.csv here - it is load-bearing for
    # published numbers - but it must not go unnoticed.
    dupes = {k: v for k, v in by_key.items() if len(v) > 1}
    by_skel = defaultdict(set)
    for r in rows:
        if r["inchikey"]:
            by_skel[r["inchikey"][:14]].add(r["struct_id"])
    skel_dupes = {k: v for k, v in by_skel.items() if len(v) > 1}
    print(f"  existing rows                        : {len(rows)}")
    print(f"  ...with approved == 1                : {sum(1 for r in rows if r['approved'] == '1')}")
    print(f"  distinct full InChIKeys              : {len(by_key)}")
    print(f"  full InChIKeys on >1 struct_id       : {len(dupes)}")
    for k, v in dupes.items():
        names = {n for _, n in v}
        tag = " <-- DIFFERENT DRUG NAMES ON ONE STRUCTURE: one row is mislabelled" \
            if len(names) > 1 else ""
        print(f"    ! {k}: {v}{tag}")
    print(f"  distinct 14-char skeletons           : {len(by_skel)}")
    print(f"  skeletons on >1 struct_id            : {len(skel_dupes)} "
          f"(salt/stereo variants of one parent)")
    return by_key, by_skel


# ---------------------------------------------------------------------------
def build(offline=False):
    print("=" * 74)
    print("STEP 0 - the existing corpus")
    print("=" * 74)
    corpus = load_corpus()
    corpus_keys, corpus_skels = report_existing_dupes(corpus)
    corpus_names = {r["name"].strip().lower() for r in corpus if r["name"].strip()}

    bio = list(csv.DictReader(BIOLOGICS.open(newline="")))
    bio_names = {r["name"].strip().lower() for r in bio if r["name"].strip()}
    print(f"  biologic corpus rows                 : {len(bio)}")

    print()
    print("=" * 74)
    print("STEP 1 - the hole: ChEMBL max_phase=4 vs the corpus")
    print("=" * 74)
    phase4 = fetch_phase4(offline)

    # Standardise every phase-4 molecule the same way the corpus was built, so
    # the comparison is like-for-like and not a raw-SMILES string match.
    std_ok, no_structure, unparseable = [], [], []
    for m in phase4:
        s = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if not s:
            no_structure.append(m)
            continue
        res = standardise(s)
        if res["status"] != "ok" or not res["inchikey"]:
            unparseable.append(m)
            continue
        res["raw_smiles"] = s
        std_ok.append((m, res))

    print(f"  ...with no SMILES in ChEMBL          : {len(no_structure)}")
    print(f"  ...unparseable / no InChIKey         : {len(unparseable)}")
    print(f"  ...standardised OK                   : {len(std_ok)}")

    present, absent = [], []
    for m, res in std_ok:
        if res["inchikey"][:14] in corpus_skels:
            present.append((m, res))
        else:
            absent.append((m, res))
    print(f"  present in corpus (skeleton match)   : {len(present)}")
    print(f"  ABSENT from corpus                   : {len(absent)}")

    # Characterise the hole rather than just sizing it.
    print()
    print("  the hole, by ChEMBL molecule_type:")
    for t, c in Counter(m.get("molecule_type") or "(none)" for m, _ in absent).most_common():
        print(f"    {t:<20} {c}")
    print("  the hole, by first_approval decade:")
    dec = Counter()
    for m, _ in absent:
        fa = m.get("first_approval")
        dec[f"{int(fa) // 10 * 10}s" if fa else "(not recorded)"] += 1
    for d, c in sorted(dec.items()):
        print(f"    {d:<20} {c}")
    n_wd = sum(1 for m, _ in absent if m.get("withdrawn_flag"))
    print(f"  the hole, withdrawn_flag = True      : {n_wd}")
    n_named = sum(1 for m, _ in absent if m.get("pref_name"))
    print(f"  the hole, with a pref_name           : {n_named} "
          f"({len(absent) - n_named} unnamed ChEMBL ids)")
    n_name_in_corpus = sum(
        1 for m, _ in absent
        if (m.get("pref_name") or "").strip().lower() in corpus_names
    )
    print(f"  the hole, whose NAME is in the corpus: {n_name_in_corpus} "
          f"(same name, different parent structure)")
    print("  structureless phase-4, by molecule_type:")
    for t, c in Counter(m.get("molecule_type") or "(none)" for m in no_structure).most_common():
        print(f"    {t:<20} {c}")

    # ------------------------------------------------------------------
    print()
    print("=" * 74)
    print("STEP 1b - rescuing structureless SMALL MOLECULES from PubChem")
    print("=" * 74)
    print("  ChEMBL types some phase-4 entries `Small molecule` yet carries no")
    print("  SMILES for them. Those are structures the corpus could hold, so")
    print("  resolve each by name against PubChem - a different database - then")
    print("  check the result: Title agreement where it agrees, and in every case")
    print("  the element-promise test, which is independent of the name index.")
    smallmol_nostruct = [m for m in no_structure
                         if (m.get("molecule_type") or "") == "Small molecule"]
    print(f"  structureless `Small molecule` entries : {len(smallmol_nostruct)}")
    rescued, pc_miss, pc_namemismatch, pc_badparse, pc_already = [], [], [], [], []
    pc_syn_rescued = []
    for i, m in enumerate(smallmol_nostruct):
        nm = (m.get("pref_name") or "").strip()
        if not nm:
            pc_miss.append((m["molecule_chembl_id"], "(no pref_name to search on)"))
            continue
        hit = fetch_pubchem_by_name(nm, offline)
        if hit is None:
            pc_miss.append((m["molecule_chembl_id"], nm))
            continue
        cid, title, smi = hit
        # Name agreement: PubChem Titles are often a salt or systematic form, so
        # require a shared alphabetic stem rather than string equality, and
        # record the Title on the row so the judgement stays checkable.
        a, b = nm.lower(), (title or "").lower()
        agree = a in b or b in a or (len(a) >= 6 and a.split()[0][:6] in b) \
            or (len(b) >= 6 and b.split()[0][:6] in a)
        how = "pubchem_title_agrees" if agree else "pubchem_name_index_only"
        # Independent of the name lookup: does the structure contain the
        # elements the name promises? Rejects a bad resolution outright.
        kept, why = element_promise_kept(nm, smi)
        if not kept:
            pc_namemismatch.append((nm, title, cid, why))
            continue
        if not agree:
            pc_syn_rescued.append((nm, title, cid))
        res = standardise(smi)
        if res["status"] != "ok" or not res["inchikey"]:
            pc_badparse.append((nm, smi[:60]))
            continue
        if res["inchikey"][:14] in corpus_skels:
            pc_already.append(nm)
            continue
        res["raw_smiles"] = smi
        rescued.append((m, res, cid, title, how))
        if i % 25 == 0:
            print(f"    resolved {i}/{len(smallmol_nostruct)}")
    n_title_ok = len(rescued) + len(pc_already) - len(pc_syn_rescued)
    print(f"  PubChem had no record for the name     : {len(pc_miss)}")
    print(f"  resolved, Title agreed                 : {n_title_ok}")
    print(f"  Title disagreed but element check passed: {len(pc_syn_rescued)} "
          f"(kept, flagged `pubchem_name_index_only`)")
    for nm, title, cid in pc_syn_rescued[:10]:
        print(f"    ~ '{nm}' -> CID {cid} '{title}'")
    print(f"  REJECTED by the element check          : {len(pc_namemismatch)}")
    for nm, title, cid, why in pc_namemismatch:
        print(f"    ! '{nm}' -> CID {cid} '{title}': {why}")
    print(f"  failed standardisation                 : {len(pc_badparse)}")
    print(f"  resolved but ALREADY in the corpus     : {len(pc_already)}")
    print(f"  RESCUED (new structures)               : {len(rescued)}")

    print()
    print("=" * 74)
    print("STEP 2 - named active moieties ChEMBL files below phase 4")
    print("=" * 74)
    extra = []
    for cid, (label, reason) in ACTIVE_MOIETIES.items():
        m = fetch_by_id(cid, offline)
        if m is None:
            print(f"  ! {cid} not found in ChEMBL - skipped")
            continue
        s = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if not s:
            print(f"  ! {cid} {label} has no SMILES - skipped")
            continue
        res = standardise(s)
        if res["status"] != "ok" or not res["inchikey"]:
            print(f"  ! {cid} {label} failed standardisation - skipped")
            continue
        if res["inchikey"][:14] in corpus_skels:
            print(f"  ! {cid} {label} is ALREADY in the corpus - not added")
            continue
        print(f"  + {cid} {label} max_phase={m.get('max_phase')} -> {res['inchikey']}")
        res["raw_smiles"] = s
        print(f"      basis: {reason}")
        extra.append((m, res, "active_moiety_of_approved_prodrug", reason,
                      "ChEMBL 36 REST API",
                      f"{CHEMBL}/molecule/{cid}.json"))

    # ------------------------------------------------------------------
    print()
    print("=" * 74)
    print("STEP 3 - target annotations for the added molecules")
    print("=" * 74)
    to_add = (
        [(m, r, "chembl_max_phase_4", "", "ChEMBL 36 REST API",
          f"{CHEMBL}/molecule/{m['molecule_chembl_id']}.json")
         for m, r in absent]
        + [(m, r, "chembl_max_phase_4",
            f"no structure in ChEMBL; SMILES from PubChem CID {cid} "
            f"(Title: {title}; identity confirmed by {how})",
            "PubChem PUG REST (name lookup; ChEMBL has no structure)",
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            f"/property/Title,SMILES/JSON")
           for m, r, cid, title, how in rescued]
        + extra
    )
    ids = {m["molecule_chembl_id"] for m, *_ in to_add}
    mech = fetch_mechanisms(ids, offline)
    print(f"  molecules with >=1 ChEMBL mechanism  : {len(mech)} / {len(ids)}")
    act = fetch_activities(ids, offline)
    print(f"  molecules with >=1 bioactivity target: {len(act)} / {len(ids)} "
          f"(human, pChEMBL >= 6, >= 2 measurements)")
    tgt_ids = {t for v in mech.values() for t, _, _ in v} | \
              {t for v in act.values() for t in v}
    print(f"  distinct targets to resolve          : {len(tgt_ids)}")
    genes = fetch_target_genes(tgt_ids, offline)
    n_human = sum(1 for v in genes.values() if v)
    print(f"  ...resolving to human gene symbols   : {n_human}")

    # ------------------------------------------------------------------
    print()
    print("=" * 74)
    print("STEP 4 - write, deduping within the additions by full InChIKey")
    print("=" * 74)
    # Standardisation collapses salt forms onto a shared parent, so two ChEMBL
    # drugs can land on one InChIKey (sodium vs calcium carbonate). The file
    # keeps one row per parent - but the collapsed names go into `notes` of the
    # surviving row, so nothing disappears without trace.
    rows, seen_key, internal_collisions = [], {}, []
    collapsed = defaultdict(list)
    for m, res, basis, note, source, source_url in to_add:
        key = res["inchikey"]
        label = m.get("pref_name") or m["molecule_chembl_id"]
        if key in seen_key:
            internal_collisions.append((key, seen_key[key], label))
            collapsed[key].append(f"{label} ({m['molecule_chembl_id']})")
            continue
        seen_key[key] = label

        cid = m["molecule_chembl_id"]
        g_moa = set()
        for tid, _action, _ in mech.get(cid, []):
            g_moa |= set(genes.get(tid) or [])  # a mechanism row IS an MOA claim
        g_act = set()
        for tid in act.get(cid, ()):
            g_act |= set(genes.get(tid) or [])
        g_all = g_moa | g_act
        evidence = ";".join(filter(None, [
            "chembl_mechanism" if g_moa else "",
            "chembl_pchembl>=6,n>=2" if g_act else "",
        ]))
        flags = res["flags"]
        if metal_stripped(res.get("raw_smiles", res["smiles"]), res["smiles"]):
            flags = ";".join(filter(None, [flags, "metal_stripped"]))
        if key in corpus_keys:  # cannot happen - skeleton test is weaker - but check
            flags = ";".join(filter(None, [flags, "collides_with_corpus"]))

        rows.append({
            "struct_id": "",  # assigned after the sort, so ids are stable per run
            "name": (m.get("pref_name") or cid).strip().lower(),
            "smiles": res["smiles"],
            "inchikey": key,
            "approved": 1,
            "approval_agencies": "",   # ChEMBL does not carry the agency
            "first_approval": m.get("first_approval") or "",
            "chembl_max_phase4": 1 if str(m.get("max_phase")) in ("4", "4.0") else 0,
            "n_targets": len(g_all),
            "targets": ";".join(sorted(g_all)),
            "moa_targets": ";".join(sorted(g_moa)),
            "flags": flags,
            "source": source,
            "source_url": source_url,
            "withdrawn": int(bool(m.get("withdrawn_flag"))),
            "chembl_id": cid,
            "molecule_type": m.get("molecule_type") or "",
            "chembl_max_phase": m.get("max_phase") or "",
            "approval_basis": basis,
            "target_evidence": evidence,
            "notes": note,
        })

    print(f"  internal InChIKey collisions         : {len(internal_collisions)}")
    for k, first, second in internal_collisions[:10]:
        print(f"    ! {k}: kept '{first}', dropped '{second}'")
    if len(internal_collisions) > 10:
        print(f"    ... and {len(internal_collisions) - 10} more")

    # fold the dropped names back into the survivor's notes and flags
    for r in rows:
        if r["inchikey"] in collapsed:
            r["flags"] = ";".join(filter(None, [r["flags"], "dup_parent_inchikey"]))
            r["notes"] = "; ".join(filter(None, [
                r["notes"],
                "same parent structure as: " + ", ".join(collapsed[r["inchikey"]]),
            ]))

    rows.sort(key=lambda r: (r["name"], r["chembl_id"]))
    for i, r in enumerate(rows):
        r["struct_id"] = ID_OFFSET + i

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows)} rows -> {OUT}")

    # ------------------------------------------------------------------
    print()
    print("=" * 74)
    print("STEP 5 - F2 (thrombin) coverage, before and after")
    print("=" * 74)
    def f2(rs, key, appr=True):
        return [r for r in rs
                if (not appr or str(r["approved"]) == "1")
                and "F2" in str(r[key]).split(";")]
    b_t, b_m = f2(corpus, "targets"), f2(corpus, "moa_targets")
    a_t, a_m = f2(rows, "targets"), f2(rows, "moa_targets")
    print(f"  approved rows with F2 annotated      : {len(b_t)} -> {len(b_t) + len(a_t)}")
    print(f"  approved rows with F2 as mechanism   : {len(b_m)} -> {len(b_m) + len(a_m)}")
    print(f"  existing F2 mechanism drugs          : {sorted(r['name'] for r in b_m)}")
    print(f"  added F2 mechanism drugs             : {sorted(r['name'] for r in a_m)}")
    print(f"  added F2 annotated drugs             : {sorted(r['name'] for r in a_t)}")

    # ------------------------------------------------------------------
    print()
    print("=" * 74)
    print("STEP 6 - the six named drugs, one by one")
    print("=" * 74)
    added_keys = {r["inchikey"]: r for r in rows}
    not_added = []
    for name in NAMED:
        mols = fetch_by_name(name, offline)
        in_corpus = [r for r in corpus if r["name"].strip().lower() == name.lower()]
        in_bio = name.lower() in bio_names
        if not mols:
            print(f"  {name:<12} NOT ADDED - no ChEMBL record under this pref_name")
            not_added.append((name, "no ChEMBL molecule with this preferred name; "
                                    "not a distinct approved drug entity"))
            continue
        m = mols[0]
        cid = m["molecule_chembl_id"]
        smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if not smi:
            where = "already in data/approved_biologics.csv" if in_bio else \
                    "NOT in data/approved_biologics.csv either"
            print(f"  {name:<12} NOT ADDED - {m.get('molecule_type')}, no SMILES in "
                  f"ChEMBL ({cid}); {where}")
            not_added.append((name, f"{m.get('molecule_type')} with no single small-molecule "
                                    f"structure ({cid}); belongs in the biologic corpus - {where}"))
            continue
        res = standardise(smi)
        hit = added_keys.get(res["inchikey"])
        if hit:
            print(f"  {name:<12} ADDED    struct_id={hit['struct_id']} "
                  f"{hit['inchikey']} basis={hit['approval_basis']}")
            print(f"               {hit['source_url']}")
            print(f"               targets={hit['targets'] or '(none)'} "
                  f"withdrawn={hit['withdrawn']}")
        elif in_corpus:
            print(f"  {name:<12} NOT ADDED - already in the corpus as "
                  f"struct_id={in_corpus[0]['struct_id']}")
            not_added.append((name, f"already present (struct_id {in_corpus[0]['struct_id']})"))
        else:
            print(f"  {name:<12} NOT ADDED - standardises to {res.get('inchikey')} "
                  f"which matched the corpus skeleton")
            not_added.append((name, "parent structure already in the corpus under another name"))

    print()
    print("  NOT ADDED / WHY")
    for name, why in not_added:
        print(f"    {name:<12} {why}")

    return corpus, corpus_keys, corpus_skels, rows, absent, no_structure, unparseable


# ---------------------------------------------------------------------------
def selftest(corpus, corpus_keys, corpus_skels, rows):
    print()
    print("=" * 74)
    print("SELF-TEST")
    print("=" * 74)
    from rdkit import Chem

    fails = []

    def check(label, ok, detail=""):
        print(f"  {'PASS' if ok else 'FAIL'}  {label}{(' - ' + detail) if detail else ''}")
        if not ok:
            fails.append(label)

    written = list(csv.DictReader(OUT.open(newline="")))
    check("file re-reads", len(written) == len(rows), f"n={len(written)}")
    check("schema starts with the corpus schema",
          list(written[0])[:12] == BASE_FIELDS)
    check("provenance columns present",
          all(c in written[0] for c in ("source", "source_url", "withdrawn")))

    bad = [r["name"] for r in written if Chem.MolFromSmiles(r["smiles"]) is None]
    check("every SMILES parses", not bad, f"{len(bad)} bad")

    recomputed = 0
    for r in written:
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is not None and Chem.MolToInchiKey(mol) != r["inchikey"]:
            recomputed += 1
    check("InChIKey recomputes from the stored SMILES", recomputed == 0,
          f"{recomputed} mismatches")

    keys = [r["inchikey"] for r in written]
    check("InChIKeys unique within the file", len(set(keys)) == len(keys),
          f"{len(keys) - len(set(keys))} dupes")
    check("no empty InChIKey", all(keys))

    clash = [k for k in keys if k in corpus_keys]
    check("no full-InChIKey collision with the corpus", not clash, f"{len(clash)}")
    skel_clash = [k for k in keys if k[:14] in corpus_skels]
    check("no skeleton collision with the corpus", not skel_clash, f"{len(skel_clash)}")

    existing_ids = {int(r["struct_id"]) for r in corpus}
    new_ids = [int(r["struct_id"]) for r in written]
    check("struct_ids unique", len(set(new_ids)) == len(new_ids))
    check("struct_ids do not collide with the corpus",
          not (set(new_ids) & existing_ids),
          f"corpus max {max(existing_ids)}, ours min {min(new_ids)}")
    check("struct_ids sit above the offset", min(new_ids) >= ID_OFFSET)

    check("every row has a source_url", all(r["source_url"] for r in written))
    check("withdrawn is 0/1", {r["withdrawn"] for r in written} <= {"0", "1"})
    check("n_targets matches targets",
          all(int(r["n_targets"]) == len([g for g in r["targets"].split(";") if g])
              for r in written))

    print(f"\n  {'ALL PASS' if not fails else str(len(fails)) + ' FAILURE(S): ' + ', '.join(fails)}")
    return not fails


# ---------------------------------------------------------------------------
def spotcheck(offline, names=("dabigatran", "melagatran")):
    """InChIKey -> PubChem -> name, an independent check of every structure.

    A mistyped or mis-fetched SMILES does not announce itself; this is the only
    thing that catches it. PubChem is a different database from ChEMBL, so a
    round-trip through it is genuinely independent evidence.
    """
    print()
    print("=" * 74)
    print("SPOT-CHECK - InChIKey -> PubChem name round-trip")
    print("=" * 74)
    written = {r["name"]: r for r in csv.DictReader(OUT.open(newline=""))}
    picks = [written[n] for n in names if n in written]
    # plus a deterministic spread across the file, so the check is not only the
    # rows we already believe in
    rest = sorted(written.values(), key=lambda r: r["inchikey"])
    picks += [r for r in rest[:: max(1, len(rest) // 8)] if r not in picks][:8]

    ok = bad = err = 0
    for r in picks:
        key = r["inchikey"]
        try:
            payload = _cached_get(
                f"pubchem_{key}",
                f"{PUBCHEM}/compound/inchikey/{key}/property/Title,SMILES/JSON",
                {}, offline, allow_404=True,
            )
            if not payload:
                print(f"  ?  {r['name']:<28} {key}  PubChem has no compound with "
                      f"this InChIKey")
                err += 1
                continue
            props = payload["PropertyTable"]["Properties"][0]
            title = props.get("Title", "")
            cid = props.get("CID")
        except Exception as exc:  # noqa: BLE001
            print(f"  ?  {r['name']:<28} {key}  PubChem lookup failed: "
                  f"{exc.__class__.__name__}")
            err += 1
            continue
        agree = r["name"].split()[0][:6].lower() in title.lower() or \
                title.split()[0][:6].lower() in r["name"].lower()
        print(f"  {'OK' if agree else '!!':<2} {r['name']:<28} {key}")
        print(f"       PubChem CID {cid}: {title}")
        ok += agree
        bad += not agree
    print(f"\n  agreed {ok}, disagreed {bad}, lookup failed {err} (n={len(picks)})")
    print("  'disagreed' is not automatically wrong: PubChem's Title is often a")
    print("  salt, a brand name or a systematic name. Read each one.")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    offline = "--offline" in sys.argv
    corpus, ckeys, cskels, rows, absent, nostruct, unparse = build(offline)
    if "--selftest" in sys.argv or "--spotcheck" in sys.argv:
        ok = selftest(corpus, ckeys, cskels, rows)
        if "--spotcheck" in sys.argv:
            spotcheck(offline)
        sys.exit(0 if ok else 1)
