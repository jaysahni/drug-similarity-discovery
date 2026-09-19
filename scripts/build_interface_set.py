"""Assemble the benchmark set for the hotspot-recovery experiment (E2).

The experiment needs proteins that have been crystallised with SEVERAL different
drug-like ligands, so that a consensus interface signature can be built from some
of them and tested against a held-out one.

Targets are not picked by hand. They are taken from DrugCentral's
drug.target.interaction.tsv - i.e. from proteins that approved drugs actually
bind - and ranked by how many distinct approved drugs hit them. That ties this
experiment to the same corpus the ligand-side retrieval benchmark uses.

For each candidate UniProt accession this queries the RCSB search API for its PDB
entries, then pulls entry metadata in batches from the RCSB GraphQL endpoint, and
keeps entries that are X-ray, <= RESOLUTION_MAX, and carry a drug-like
non-polymer ligand.

Usage:
    ./env/bin/python scripts/build_interface_set.py            # top 200 targets
    ./env/bin/python scripts/build_interface_set.py --limit 50
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import ssl
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import certifi

ROOT = Path(__file__).resolve().parent.parent
INTERACTIONS = ROOT / "data" / "drug.target.interaction.tsv"
CACHE = ROOT / "data" / "raw" / "interface_set"
OUT_INVENTORY = CACHE / "inventory.json"
OUT_SUMMARY = ROOT / "results" / "interface_set_summary.json"

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

RESOLUTION_MAX = 3.0
MIN_LIGAND_MW = 180.0      # below this is almost always buffer, cryoprotectant or ion
MIN_LIGANDS_PER_TARGET = 3  # need >=3 distinct ligands to hold one out and still have a consensus

# Crystallisation additives, cryoprotectants, buffers and ions. These pass a naive
# molecular-weight filter but are not binding events, so they are named explicitly
# rather than left to a threshold. Counts of what each rule removed are reported.
ADDITIVES = {
    "HOH", "DOD", "SO4", "PO4", "NO3", "CL", "NA", "MG", "ZN", "CA", "K", "MN",
    "FE", "FE2", "CU", "CU1", "NI", "CD", "HG", "IOD", "BR", "F", "CO", "CS",
    "GOL", "EDO", "PEG", "PGE", "PG4", "P6G", "1PE", "2PE", "MPD", "DMS", "DMF",
    "ACT", "ACY", "FMT", "CIT", "FLC", "TRS", "MES", "EPE", "BTB", "TAR", "MLI",
    "IMD", "BME", "DTT", "DTU", "IPA", "MOH", "EOH", "ACN", "URE", "GTT", "SIN",
    "SCN", "AZI", "NH4", "TLA", "OXL", "MAL", "SRT", "BCT", "CO3", "NHE", "CXS",
    "LDA", "BOG", "C8E", "OCT", "HEZ", "12P", "15P", "33O", "PIN", "POP", "PPV",
}


def http_post(url, payload, retries=4, pause=1.5):
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90, context=SSL_CTX) as r:
                raw = r.read()
                # RCSB search answers "no hits" with 204 and an empty body, which
                # urllib reports as success - so an empty body is a real zero-hit
                # result, not a transport failure. Verified: Q969S8, Q13621 and
                # P54750 all return HTTP 204 / 0 bytes while P35968 returns 200.
                if not raw.strip():
                    return None
                return json.loads(raw)
        except urllib.error.HTTPError as exc:
            if exc.code == 204:
                return None
            if attempt == retries - 1:
                raise
        except Exception:            # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
        time.sleep(pause * (attempt + 1))
    return None


def drug_targets():
    """UniProt accession -> (gene symbol, set of distinct approved drug names)."""
    drugs = collections.defaultdict(set)
    genes = {}
    with INTERACTIONS.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("ACCESSION") or "").strip()
            org = (row.get("ORGANISM") or "").strip()
            if not acc or org != "Homo sapiens":
                continue
            # a row can carry several accessions for a complex; take each
            for a in acc.replace(",", "|").split("|"):
                a = a.strip()
                if not a:
                    continue
                drugs[a].add((row.get("DRUG_NAME") or "").strip())
                genes.setdefault(a, (row.get("GENE") or "").strip())
    return genes, drugs


def entries_for_uniprot(acc):
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 500},
                            "results_content_type": ["experimental"]},
    }
    r = http_post(SEARCH_URL, q)
    if not r:
        return []
    return [h["identifier"] for h in r.get("result_set", [])]


ENTRY_QUERY = """{ entries(entry_ids: %s) {
  rcsb_id
  exptl { method }
  rcsb_entry_info { resolution_combined }
  nonpolymer_entities {
    nonpolymer_comp { chem_comp { id name formula_weight } }
  }
} }"""


def entry_metadata(ids, batch=40, pause=0.2):
    out = {}
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        r = http_post(GRAPHQL_URL, {"query": ENTRY_QUERY % json.dumps(chunk)})
        for e in (r or {}).get("data", {}).get("entries", None) or []:
            out[e["rcsb_id"]] = e
        time.sleep(pause)
    return out


def drug_like_ligands(entry):
    """Ligand comp ids in this entry that plausibly represent a binding event."""
    keep = []
    for ne in entry.get("nonpolymer_entities") or []:
        comp = ((ne or {}).get("nonpolymer_comp") or {}).get("chem_comp") or {}
        cid, mw = comp.get("id"), comp.get("formula_weight")
        if not cid or cid in ADDITIVES or mw is None or mw < MIN_LIGAND_MW:
            continue
        keep.append({"comp_id": cid, "name": comp.get("name"), "mw": mw})
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=200,
                    help="how many top drug-bound targets to query RCSB for")
    args = ap.parse_args()

    CACHE.mkdir(parents=True, exist_ok=True)
    genes, drugs = drug_targets()
    ranked = sorted(drugs, key=lambda a: -len(drugs[a]))[:args.limit]
    print(f"{len(drugs)} human UniProt accessions bound by >=1 DrugCentral drug")
    print(f"querying RCSB for the top {len(ranked)} by distinct drug count\n")

    inventory, stats = {}, collections.Counter()
    for n, acc in enumerate(ranked, 1):
        try:
            ids = entries_for_uniprot(acc)
        except Exception as exc:                       # noqa: BLE001
            print(f"  [{n}/{len(ranked)}] {acc} SEARCH FAILED: {exc}", file=sys.stderr)
            stats["search_failed"] += 1
            continue
        stats["entries_seen"] += len(ids)
        if not ids:
            stats["no_pdb_entries"] += 1   # target has no experimental PDB entry
            continue

        meta = entry_metadata(ids)
        kept, ligands = [], {}
        for pdb_id, e in meta.items():
            methods = [m.get("method") for m in (e.get("exptl") or [])]
            if "X-RAY DIFFRACTION" not in methods:
                stats["dropped_not_xray"] += 1
                continue
            res = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            if not res or res[0] > RESOLUTION_MAX:
                stats["dropped_resolution"] += 1
                continue
            ligs = drug_like_ligands(e)
            if not ligs:
                stats["dropped_no_drug_like_ligand"] += 1
                continue
            kept.append({"pdb_id": pdb_id, "resolution": res[0],
                         "ligands": [l["comp_id"] for l in ligs]})
            for l in ligs:
                ligands.setdefault(l["comp_id"], l)

        if len(ligands) >= MIN_LIGANDS_PER_TARGET:
            inventory[acc] = {
                "uniprot": acc, "gene": genes.get(acc, ""),
                "n_approved_drugs": len(drugs[acc]),
                "n_entries_total": len(ids), "n_entries_kept": len(kept),
                "n_distinct_ligands": len(ligands),
                "entries": sorted(kept, key=lambda d: d["resolution"]),
                "ligands": ligands,
            }
            stats["targets_kept"] += 1
        else:
            stats["targets_too_few_ligands"] += 1

        if n % 10 == 0 or n == len(ranked):
            print(f"  [{n}/{len(ranked)}] kept {stats['targets_kept']} targets "
                  f"so far (last: {genes.get(acc,'?')} {acc}, "
                  f"{len(ligands)} ligands)", flush=True)

    OUT_INVENTORY.write_text(json.dumps(inventory, indent=1))
    n_entries = sum(t["n_entries_kept"] for t in inventory.values())
    n_ligs = sum(t["n_distinct_ligands"] for t in inventory.values())
    summary = {
        "targets_queried": len(ranked),
        "targets_kept": len(inventory),
        "usable_entries": n_entries,
        "distinct_ligands": n_ligs,
        "filters": {"resolution_max": RESOLUTION_MAX, "min_ligand_mw": MIN_LIGAND_MW,
                    "min_ligands_per_target": MIN_LIGANDS_PER_TARGET,
                    "n_additives_blocklisted": len(ADDITIVES)},
        "drop_counts": dict(stats),
        "top_targets": [{"gene": t["gene"], "uniprot": t["uniprot"],
                         "entries": t["n_entries_kept"], "ligands": t["n_distinct_ligands"]}
                        for t in sorted(inventory.values(),
                                        key=lambda t: -t["n_distinct_ligands"])[:25]],
    }
    OUT_SUMMARY.parent.mkdir(exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2))

    print(f"\n{len(inventory)} targets kept, {n_entries} usable entries, "
          f"{n_ligs} distinct drug-like ligands")
    print(f"drop counts: {dict(stats)}")
    print(f"wrote {OUT_INVENTORY} and {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
