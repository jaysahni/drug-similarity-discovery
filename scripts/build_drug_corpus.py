"""Build the unfiltered approved-drug corpus for the repurposing pipeline (T1).

`data/drugs.csv` (from prepare_data.py) keeps only drugs with >=1 annotated human
target. Repurposing needs the opposite: every DrugCentral structure with a
parseable SMILES, target-annotated or not, because a candidate may have no target
record at all. So this script keeps ALL structures and puts approval in a column
instead of filtering on it.

Approval status is NOT in structures.smiles.tsv (its columns are
SMILES/InChI/InChIKey/ID/INN/CAS_RN - checked at runtime, see check_columns()).
It is taken from DrugCentral's own `approval` table, streamed out of the release
SQL dump, and cross-checked against ChEMBL max_phase=4 matched on parent InChIKey
skeleton. Both are reported; `approved` is DrugCentral's, `chembl_max_phase4` is
the independent check.

Standardisation per structure: RDKit parse -> largest organic fragment (salt
strip) -> neutralise where unambiguous -> canonical SMILES -> InChIKey. Nothing
is dropped for being inorganic, an element or a mixture; those are flagged.

Output: data/approved_drugs.csv
    struct_id, name, smiles, inchikey, approved, approval_agencies,
    first_approval, chembl_max_phase4, n_targets, targets, moa_targets, flags

Usage:
    ./env/bin/python scripts/build_drug_corpus.py
    ./env/bin/python scripts/build_drug_corpus.py --offline   # use data/raw cache only
"""

from __future__ import annotations

import csv
import gzip
import json
import ssl
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import certifi
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RAW = DATA / "raw"
STRUCTURES = DATA / "structures.smiles.tsv"
INTERACTIONS = DATA / "drug.target.interaction.tsv"
OUT = DATA / "approved_drugs.csv"

# 2021_10_05 dump == dbversion 49 (2021-09-24), the release the local TSVs came from.
DUMP_URL = "https://unmtid-shinyapps.net/download/drugcentral.dump.010_05_2021.sql.gz"
APPROVAL_CACHE = RAW / "drugcentral_approval.tsv"
CHEMBL_CACHE = RAW / "chembl_max_phase4.tsv"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

# The python.org 3.14 build ships no CA bundle; urllib needs certifi's explicitly.
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

RDLogger.DisableLog("rdApp.*")  # we count parse failures ourselves


# --------------------------------------------------------------------------
# approval status - neither of these is in the local TSVs
# --------------------------------------------------------------------------
def fetch_drugcentral_approval(offline=False):
    """DrugCentral's `approval` table, streamed out of the release SQL dump.

    pg_dump emits table data alphabetically, so `approval` lands ~20 MB into the
    decompressed stream; we stop reading there rather than pull the full 1.3 GB.
    """
    if APPROVAL_CACHE.exists():
        return APPROVAL_CACHE
    if offline:
        sys.exit(f"--offline but no cache at {APPROVAL_CACHE}")

    print(f"streaming approval table from {DUMP_URL} ...")
    rows, cols, inside, nbytes = [], None, False, 0
    with urllib.request.urlopen(DUMP_URL, timeout=600, context=SSL_CTX) as resp:
        for raw in gzip.GzipFile(fileobj=resp):
            nbytes += len(raw)
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if not inside:
                # trailing " (" matters: approval_type is a different table
                if line.startswith("COPY public.approval ("):
                    cols = line[line.index("(") + 1 : line.index(")")].split(", ")
                    inside = True
                continue
            if line.startswith("\\."):
                break
            rows.append(line)
    if cols is None:
        sys.exit("no `COPY public.approval (` block in the dump - schema changed?")

    RAW.mkdir(parents=True, exist_ok=True)
    with APPROVAL_CACHE.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\n".join(rows) + "\n")
    print(f"  read {nbytes / 1e6:.1f} MB decompressed, cached {len(rows)} approval rows")
    return APPROVAL_CACHE


def load_approval():
    """struct_id -> (sorted agencies, earliest approval date)."""
    agencies, dates, seen = defaultdict(set), defaultdict(list), set()
    with APPROVAL_CACHE.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r["struct_id"]
            if not sid or sid == "\\N":
                continue
            seen.add(sid)  # presence is the approval fact; type/date may be NULL
            if r.get("type") and r["type"] != "\\N":
                agencies[sid].add(r["type"])
            if r.get("approval") and r["approval"] != "\\N":
                dates[sid].append(r["approval"])
    return {
        sid: (";".join(sorted(agencies.get(sid, ()))), min(dates.get(sid, [""])))
        for sid in seen
    }


def fetch_chembl_phase4(offline=False):
    """Every ChEMBL molecule with max_phase=4 and a structure (paged REST)."""
    if CHEMBL_CACHE.exists():
        return CHEMBL_CACHE
    if offline:
        sys.exit(f"--offline but no cache at {CHEMBL_CACHE}")

    print("fetching ChEMBL max_phase=4 molecules ...")
    rows, offset, limit = [], 0, 1000
    while True:
        params = {
            "max_phase": 4,
            "molecule_structures__isnull": "false",
            "limit": limit,
            "offset": offset,
            "only": "molecule_chembl_id,pref_name,first_approval,molecule_structures",
        }
        url = f"{CHEMBL_API}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(url, timeout=120, context=SSL_CTX) as resp:
            page = json.load(resp)
        mols = page["molecules"]
        if not mols:
            break
        for m in mols:
            s = m.get("molecule_structures") or {}
            if not s.get("canonical_smiles"):
                continue
            rows.append(
                {
                    "chembl_id": m["molecule_chembl_id"],
                    "pref_name": m.get("pref_name") or "",
                    "first_approval": m.get("first_approval") or "",
                    "inchikey": s.get("standard_inchi_key") or "",
                    "smiles": s["canonical_smiles"],
                }
            )
        print(f"  {len(rows)}/{page['page_meta']['total_count']}")
        offset += limit
        if offset >= page["page_meta"]["total_count"]:
            break

    if not rows:
        sys.exit("ChEMBL returned no max_phase=4 molecules with structures")
    RAW.mkdir(parents=True, exist_ok=True)
    with CHEMBL_CACHE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return CHEMBL_CACHE


# --------------------------------------------------------------------------
# structure standardisation
# --------------------------------------------------------------------------
LFC = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
UNCHARGER = rdMolStandardize.Uncharger()


def standardise(smiles):
    """Raw SMILES -> dict with a `status` of ok / parse_fail / no_parent.

    Only unparseable or empty-parent inputs fail; inorganics, single elements
    and mixtures are flagged, never dropped.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "parse_fail"}

    flags = []
    n_frags = len(Chem.GetMolFrags(mol))
    if n_frags > 1:
        flags.append("multi_fragment")

    before = Chem.MolToSmiles(mol)
    parent = LFC.choose(mol)
    if parent is None or parent.GetNumAtoms() == 0:
        return {"status": "no_parent"}
    stripped = Chem.MolToSmiles(parent) != before

    after_strip = Chem.MolToSmiles(parent)
    parent = UNCHARGER.uncharge(parent)
    neutralised = Chem.MolToSmiles(parent) != after_strip

    if not any(a.GetSymbol() == "C" for a in parent.GetAtoms()):
        flags.append("inorganic")
    if parent.GetNumHeavyAtoms() == 1:
        flags.append("element")
    if Chem.GetFormalCharge(parent) != 0:
        flags.append("charged_parent")  # e.g. quaternary ammonium - not neutralisable

    try:
        inchikey = Chem.MolToInchiKey(parent)
    except Exception:  # noqa: BLE001 - InChI rejects some valid-to-RDKit structures
        inchikey = ""
    if not inchikey:
        flags.append("no_inchikey")

    return {
        "status": "ok",
        "smiles": Chem.MolToSmiles(parent),
        "inchikey": inchikey,
        "flags": ";".join(flags),
        "stripped": stripped,
        "neutralised": neutralised,
    }


# --------------------------------------------------------------------------
# local DrugCentral TSVs
# --------------------------------------------------------------------------
def check_columns():
    """Print what the structures TSV actually carries - don't assume an approval column."""
    with STRUCTURES.open(newline="") as fh:
        cols = next(csv.reader(fh, delimiter="\t"))
    approvalish = [c for c in cols if any(k in c.lower() for k in ("approv", "phase", "status"))]
    print(f"structures.smiles.tsv columns: {cols}")
    print(f"  approval-like columns present: {approvalish or 'NONE'}")
    return cols, approvalish


def load_structures():
    with STRUCTURES.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_targets():
    """struct_id -> (human gene symbols, mechanism-of-action gene symbols).

    Same convention as prepare_data.py: human rows only, '|' splits complexes.
    """
    allt, moat, any_row = defaultdict(set), defaultdict(set), set()
    with INTERACTIONS.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid, gene = r.get("STRUCT_ID"), (r.get("GENE") or "").strip()
            if not sid:
                continue
            any_row.add(sid)
            if not gene or r.get("ORGANISM") != "Homo sapiens":
                continue
            genes = {g.strip() for g in gene.split("|") if g.strip()}
            allt[sid] |= genes
            if (r.get("MOA") or "").strip() == "1":
                moat[sid] |= genes
    return allt, moat, any_row


# --------------------------------------------------------------------------
def main():
    offline = "--offline" in sys.argv

    print("=" * 72)
    print("STEP 0 - what the local structures file actually contains")
    print("=" * 72)
    check_columns()

    fetch_drugcentral_approval(offline)
    fetch_chembl_phase4(offline)

    raw_rows = load_structures()
    approval = load_approval()
    allt, moat, any_interaction = load_targets()

    print()
    print("=" * 72)
    print("STEP 1 - standardisation")
    print("=" * 72)

    n_raw = len(raw_rows)
    n_blank_id = n_empty_smiles = n_dup_id = n_parse_fail = n_no_parent = 0
    n_stripped = n_neutralised = n_key_differs = 0
    fail_examples = []
    seen_ids = set()
    rows = []

    # Every raw row lands in exactly one bucket, so the reconciliation below is
    # an identity. Empty SMILES is tested before the duplicate guard so a blank
    # row cannot shadow a later duplicate that carries the real structure.
    for r in raw_rows:
        sid, smi = (r.get("ID") or "").strip(), (r.get("SMILES") or "").strip()
        if not sid:
            n_blank_id += 1
            continue
        if not smi:
            n_empty_smiles += 1
            continue
        if sid in seen_ids:
            n_dup_id += 1
            continue
        seen_ids.add(sid)

        std = standardise(smi)
        if std["status"] != "ok":
            if std["status"] == "parse_fail":
                n_parse_fail += 1
            else:
                n_no_parent += 1
            if len(fail_examples) < 5:
                fail_examples.append((std["status"], sid, r.get("INN", ""), smi[:60]))
            continue

        n_stripped += std["stripped"]
        n_neutralised += std["neutralised"]
        src_key = (r.get("InChIKey") or "").strip()
        if src_key and std["inchikey"] and src_key != std["inchikey"]:
            n_key_differs += 1

        targets = sorted(allt.get(sid, ()))
        agencies, first = approval.get(sid, ("", ""))
        rows.append(
            {
                "struct_id": sid,
                "name": (r.get("INN") or "").strip(),
                "smiles": std["smiles"],
                "inchikey": std["inchikey"],
                "approved": int(sid in approval),
                "approval_agencies": agencies,
                "first_approval": first,
                "chembl_max_phase4": "",  # filled below
                "n_targets": len(targets),
                "targets": ";".join(targets),
                "moa_targets": ";".join(sorted(moat.get(sid, ()))),
                "flags": std["flags"],
            }
        )

    # denominator = rows actually handed to RDKit, so the rate measures parsing
    attempted = len(rows) + n_parse_fail + n_no_parent
    parse_rate = 100.0 * len(rows) / attempted if attempted else 0.0
    print(f"raw data rows in structures.smiles.tsv : {n_raw}")
    print(f"  blank struct_id (skipped)            : {n_blank_id}")
    print(f"  empty SMILES cell                    : {n_empty_smiles}")
    print(f"  duplicate struct_id (skipped)        : {n_dup_id}")
    print(f"  SMILES handed to RDKit               : {attempted}")
    print(f"  RDKit parse failures                 : {n_parse_fail}")
    print(f"  no parent fragment survived          : {n_no_parent}")
    print(f"  standardised OK                      : {len(rows)}")
    print(f"  SMILES parse rate (of attempted)     : {parse_rate:.2f}%")
    for status, sid, name, smi in fail_examples:
        print(f"    ! {status} struct_id {sid} {name}: {smi}")
    print(f"  largest-fragment step changed        : {n_stripped}")
    print(f"  neutralisation step changed          : {n_neutralised}")
    print(f"  parent InChIKey != DrugCentral's     : {n_key_differs}")

    flag_counts = defaultdict(int)
    for r in rows:
        for f in r["flags"].split(";"):
            if f:
                flag_counts[f] += 1
    for f, c in sorted(flag_counts.items(), key=lambda kv: -kv[1]):
        print(f"  flagged {f:<16}             : {c}")

    dup_keys = defaultdict(list)
    for r in rows:
        if r["inchikey"]:
            dup_keys[r["inchikey"]].append(r["struct_id"])
    collisions = {k: v for k, v in dup_keys.items() if len(v) > 1}
    # Two struct_ids with one parent structure means an upstream duplicate, not a
    # dedupe target: flag both so a retrieval run can see it instead of silently
    # scoring the same molecule twice.
    colliding_ids = {sid for v in collisions.values() for sid in v}
    for r in rows:
        if r["struct_id"] in colliding_ids:
            r["flags"] = ";".join(filter(None, [r["flags"], "dup_parent_inchikey"]))
    print(f"  distinct parent InChIKeys            : {len(dup_keys)}")
    print(f"  keys shared by >1 struct_id          : {len(collisions)} "
          f"({sum(len(v) for v in collisions.values())} rows, flagged dup_parent_inchikey)")
    for k, v in collisions.items():
        names = [next(r["name"] for r in rows if r["struct_id"] == sid) for sid in v]
        print(f"    ! {k}: struct_ids {v} = {names}")

    print()
    print("=" * 72)
    print("STEP 2 - approval status")
    print("=" * 72)

    # ChEMBL cross-check: same standardisation both sides, matched on the 14-char
    # InChIKey skeleton so stereochemistry/salt-form differences still match.
    chembl_skeletons, n_chembl, n_chembl_unparsed = set(), 0, 0
    with CHEMBL_CACHE.open(newline="") as fh:
        for cr in csv.DictReader(fh, delimiter="\t"):
            n_chembl += 1
            std = standardise(cr["smiles"])
            if std["status"] != "ok" or not std["inchikey"]:
                n_chembl_unparsed += 1
                continue
            chembl_skeletons.add(std["inchikey"][:14])

    n_p4 = 0
    for r in rows:
        if not r["inchikey"]:
            r["chembl_max_phase4"] = ""  # unknown - no key to match on
            continue
        hit = int(r["inchikey"][:14] in chembl_skeletons)
        r["chembl_max_phase4"] = hit
        n_p4 += hit

    n_dc = sum(r["approved"] for r in rows)
    n_unknown = sum(1 for r in rows if r["chembl_max_phase4"] == "")
    both = sum(1 for r in rows if r["approved"] and r["chembl_max_phase4"] == 1)
    dc_only = sum(1 for r in rows if r["approved"] and r["chembl_max_phase4"] == 0)
    cb_only = sum(1 for r in rows if not r["approved"] and r["chembl_max_phase4"] == 1)
    neither = sum(1 for r in rows if not r["approved"] and r["chembl_max_phase4"] == 0)

    print(f"source of truth for `approved`: DrugCentral `approval` table "
          f"(dbversion 49, streamed from the release SQL dump)")
    print(f"  ChEMBL max_phase=4 molecules fetched : {n_chembl} "
          f"({n_chembl_unparsed} unparseable, {len(chembl_skeletons)} distinct skeletons)")
    print(f"  approved=1 (DrugCentral)             : {n_dc} / {len(rows)}")
    print(f"  chembl_max_phase4=1                  : {n_p4} / {len(rows)}")
    print("  agreement:")
    print(f"    both                               : {both}")
    print(f"    DrugCentral only                   : {dc_only}")
    print(f"    ChEMBL only                        : {cb_only}")
    print(f"    neither                            : {neither}")
    print(f"    not comparable (no InChIKey)       : {n_unknown}")

    n_with_targets = sum(1 for r in rows if r["n_targets"])
    n_with_moa = sum(1 for r in rows if r["moa_targets"])
    n_interaction_no_human = sum(
        1 for r in rows if r["struct_id"] in any_interaction and not r["n_targets"]
    )
    print()
    print("=" * 72)
    print("STEP 3 - targets (informational; corpus is NOT filtered on these)")
    print("=" * 72)
    print(f"  >=1 human target                     : {n_with_targets}")
    print(f"  >=1 MOA target                       : {n_with_moa}")
    print(f"  no human target                      : {len(rows) - n_with_targets}")
    print(f"    ...of which have non-human rows    : {n_interaction_no_human}")

    if not rows:
        sys.exit("no structures survived standardisation - nothing to write")
    rows.sort(key=lambda r: int(r["struct_id"]))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print()
    print("=" * 72)
    print("RECONCILIATION + ACCEPTANCE")
    print("=" * 72)
    accounted = len(rows) + n_blank_id + n_empty_smiles + n_dup_id + n_parse_fail + n_no_parent
    print(f"  {len(rows)} written + {n_blank_id} blank id + {n_empty_smiles} empty SMILES "
          f"+ {n_dup_id} duplicate ids + {n_parse_fail} parse failures + {n_no_parent} "
          f"no-parent = {accounted} (raw {n_raw}) "
          f"{'OK' if accounted == n_raw else 'MISMATCH'}")
    print(f"  acceptance >=4000 rows               : {len(rows)} "
          f"{'PASS' if len(rows) >= 4000 else 'FAIL'}")
    print(f"  acceptance >=97% parse rate          : {parse_rate:.2f}% "
          f"{'PASS' if parse_rate >= 97 else 'FAIL'}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
