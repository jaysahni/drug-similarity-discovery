"""Download the mmCIF files the hotspot-recovery experiment (E2) needs.

Reads the inventory written by build_interface_set.py and picks, per target, a
set of entries chosen to maximise DISTINCT ligands rather than to maximise
entries: the experiment holds one ligand out and builds a consensus from the
rest, so twenty structures of the same ligand carry one ligand's worth of
information. Entries are taken greedily, best resolution first, keeping an
entry only while it still contributes a ligand not already covered; the
per-target cap then fills out with the best remaining entries.

Files land in data/raw/structures/<PDBID>.cif (gitignored) and are skipped if
already present, so this is resumable.

Usage:
    ./env/bin/python scripts/fetch_structures.py --targets 60 --per-target 25
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import ssl
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import certifi

ROOT = Path(__file__).resolve().parent.parent
INVENTORY = ROOT / "data" / "raw" / "interface_set" / "inventory.json"
OUTDIR = ROOT / "data" / "raw" / "structures"
MANIFEST = ROOT / "data" / "raw" / "structures_manifest.json"
SSL_CTX = ssl.create_default_context(cafile=certifi.where())
URL = "https://files.rcsb.org/download/{}.cif"

_lock = threading.Lock()
_done = {"n": 0, "bytes": 0, "skipped": 0, "failed": 0}


def select_entries(target, cap):
    """Entries for one target, greedy on new-ligand coverage then resolution."""
    entries = sorted(target["entries"], key=lambda e: e["resolution"])
    chosen, covered = [], set()
    for e in entries:                       # pass 1: each entry must add a ligand
        new = set(e["ligands"]) - covered
        if new:
            chosen.append(e)
            covered |= set(e["ligands"])
        if len(chosen) >= cap:
            return chosen, covered
    for e in entries:                       # pass 2: fill the cap, best resolution
        if len(chosen) >= cap:
            break
        if e not in chosen:
            chosen.append(e)
    return chosen, covered


def fetch(pdb_id, retries=3):
    dest = OUTDIR / f"{pdb_id}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        with _lock:
            _done["skipped"] += 1
        return True
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL.format(pdb_id),
                                         headers={"User-Agent": "drug-similarity-discovery/0.1"})
            with urllib.request.urlopen(req, timeout=120, context=SSL_CTX) as r:
                raw = r.read()
            if not raw:
                raise ValueError("empty body")
            tmp = dest.with_suffix(".cif.part")   # atomic: a killed run leaves no
            tmp.write_bytes(raw)                  # truncated file to poison reruns
            tmp.replace(dest)
            with _lock:
                _done["n"] += 1
                _done["bytes"] += len(raw)
                if _done["n"] % 50 == 0:
                    print(f"  {_done['n']} fetched, {_done['skipped']} cached, "
                          f"{_done['bytes'] / 1e6:.0f} MB", flush=True)
            return True
        except Exception as exc:                  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                with _lock:
                    _done["failed"] += 1
                print(f"  FAILED {pdb_id}: {exc}", file=sys.stderr)
                return False
            time.sleep(1.5 * (attempt + 1))
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=60)
    ap.add_argument("--per-target", type=int, default=25)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    inv = json.loads(INVENTORY.read_text())
    OUTDIR.mkdir(parents=True, exist_ok=True)

    ranked = sorted(inv.values(), key=lambda t: -t["n_distinct_ligands"])[:args.targets]
    plan, wanted = {}, []
    for t in ranked:
        chosen, covered = select_entries(t, args.per_target)
        plan[t["uniprot"]] = {
            "gene": t["gene"], "uniprot": t["uniprot"],
            "n_distinct_ligands_available": t["n_distinct_ligands"],
            "n_distinct_ligands_selected": len(covered),
            "entries": chosen,
        }
        wanted += [e["pdb_id"] for e in chosen]

    wanted = sorted(set(wanted))
    print(f"{len(ranked)} targets, {len(wanted)} distinct entries to fetch "
          f"(cap {args.per_target}/target)\n")

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        list(ex.map(fetch, wanted))

    MANIFEST.write_text(json.dumps(
        {"targets": plan, "n_entries_requested": len(wanted),
         "fetched": _done["n"], "already_cached": _done["skipped"],
         "failed": _done["failed"], "bytes": _done["bytes"]}, indent=1))
    print(f"\nfetched {_done['n']}, cached {_done['skipped']}, failed {_done['failed']}, "
          f"{_done['bytes'] / 1e6:.0f} MB")
    print(f"wrote {MANIFEST}")


if __name__ == "__main__":
    main()
