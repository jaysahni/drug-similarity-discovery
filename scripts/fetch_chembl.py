"""Fetch approved drugs from ChEMBL with SMILES and ATC classification.

Ground truth for "similar": two drugs share an ATC level-4 code (chemical/
therapeutic/pharmacological subgroup). Writes data/drugs.csv.

No API key required. ChEMBL REST, paged at 1000.
"""

import json
import ssl
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import certifi

BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT = Path(__file__).resolve().parent.parent / "data" / "drugs.csv"

# The python.org 3.14 framework build ships no CA bundle, so urllib fails SSL
# verification even though curl succeeds. Use certifi's explicitly.
SSL_CTX = ssl.create_default_context(cafile=certifi.where())


def get(path, params, retries=4):
    url = f"{BASE}/{path}.json?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60, context=SSL_CTX) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return None


def fetch_molecules():
    """Approved drugs (max_phase=4) with a structure."""
    rows, offset, limit = {}, 0, 1000
    while True:
        page = get(
            "molecule",
            {
                "max_phase": 4,
                "molecule_structures__isnull": "false",
                "limit": limit,
                "offset": offset,
            },
        )
        mols = page["molecules"]
        if not mols:
            break
        for m in mols:
            struct = m.get("molecule_structures") or {}
            smiles = struct.get("canonical_smiles")
            if not smiles:
                continue
            rows[m["molecule_chembl_id"]] = {
                "chembl_id": m["molecule_chembl_id"],
                "name": (m.get("pref_name") or "").strip(),
                "smiles": smiles,
            }
        print(f"  molecules: {len(rows)} (offset {offset})", file=sys.stderr)
        offset += limit
        if offset >= page["page_meta"]["total_count"]:
            break
    return rows


def fetch_atc(rows):
    """Attach ATC level-4 codes via the atc_class endpoint on each molecule."""
    ids = list(rows)
    chunk = 50
    n_with = 0
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        page = get(
            "molecule",
            {
                "molecule_chembl_id__in": ",".join(batch),
                "limit": chunk,
                "only": "molecule_chembl_id,atc_classifications",
            },
        )
        for m in page.get("molecules", []):
            codes = m.get("atc_classifications") or []
            # ATC level 4 is the first 5 characters, e.g. N02BA from N02BA01
            lvl4 = sorted({c[:5] for c in codes if len(c) >= 5})
            rows[m["molecule_chembl_id"]]["atc4"] = ";".join(lvl4)
            if lvl4:
                n_with += 1
        print(f"  atc: {i + len(batch)}/{len(ids)}, {n_with} with codes", file=sys.stderr)
    return n_with


def main():
    print("Fetching approved molecules from ChEMBL...", file=sys.stderr)
    rows = fetch_molecules()
    print(f"Got {len(rows)} molecules with structures.", file=sys.stderr)

    print("Fetching ATC classifications...", file=sys.stderr)
    n_with = fetch_atc(rows)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    import csv

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chembl_id", "name", "smiles", "atc4"])
        w.writeheader()
        for r in rows.values():
            r.setdefault("atc4", "")
            w.writerow(r)

    print(f"\nWrote {len(rows)} drugs to {OUT}", file=sys.stderr)
    print(f"  with >=1 ATC level-4 code: {n_with}", file=sys.stderr)


if __name__ == "__main__":
    main()
