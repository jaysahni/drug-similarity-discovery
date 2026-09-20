#!/usr/bin/env python3
"""Build data/target_annotations.csv: drug -> human protein target, from measured bioactivity.

Why this exists
---------------
`data/approved_drugs.csv` carries DrugCentral's curated `targets` / `moa_targets`
columns. 1,985 of its 4,099 rows (48%) have neither. A row with no target can
never be a true positive nor a scored negative in a retrieval benchmark, so it is
dead weight in the denominator. This script adds a second, independent evidence
layer drawn from *measured* bioactivity (ChEMBL activities carrying a pChEMBL
value) and emits both layers in one joinable table.

Design
------
* InChIKey is the primary join key (every row in approved_drugs.csv has one).
  Name matching is a labelled fallback, never silently mixed in.
* Only SINGLE PROTEIN ChEMBL targets are kept, so every row maps to exactly one
  UniProt accession.
* Human gene symbols come from ChEMBL's GENE_SYMBOL component synonyms, which are
  HGNC symbols -- the same vocabulary as the existing `targets` column.
* Non-human (pathogen) targets are kept under a separate `evidence` tier,
  `chembl_pchembl_nonhuman`. They are most of what is left once human targets are
  exhausted, because the unannotated tail is dominated by anti-infectives.
* Every raw API page is cached under data/raw/chembl/ (already gitignored), so a
  re-run costs no network. Delete that directory to force a refetch.
* No credentials. If ChEMBL is unreachable the script emits the DrugCentral
  evidence layer alone and says so rather than failing.

Threshold
---------
Rows are FETCHED at pChEMBL >= 5.0 and all of them are written out, each carrying
its `pchembl` value plus a `passes_threshold` flag for the default cut of 6.0
(1 uM). Keeping the 5.0-6.0 band in the file lets a consumer re-threshold without
refetching; `docs/12-TARGET-ANNOTATIONS.md` shows how the counts move.

Usage
-----
    ./env-kit/bin/python demo/enrich_targets.py
    ./env-kit/bin/python demo/enrich_targets.py --offline      # cache only, no network
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import statistics
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
DRUGS_CSV = ROOT / "data" / "approved_drugs.csv"
OUT_CSV = ROOT / "data" / "target_annotations.csv"
CACHE_DIR = ROOT / "data" / "raw" / "chembl"

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PAGE = 1000
FETCH_FLOOR = 5.0  # fetch everything at or above this
DEFAULT_THRESHOLD = 6.0  # 1 uM; the cut we recommend and flag

ACTIVITY_FIELDS = (
    "molecule_chembl_id,target_chembl_id,pchembl_value,standard_type,"
    "standard_value,standard_units,standard_relation,assay_type,assay_chembl_id,"
    "target_organism"
)

# Assay-quality filter. Measured justification is in docs/12-TARGET-ANNOTATIONS.md;
# the report prints how many rows each rule drops on every run.
#
#   standard_type: a concentration-response affinity/potency constant only.
#     'Potency' is excluded because it is the readout of the NCGC/Tox21 qHTS
#     panels deposited from PubChem -- sampling ChEMBL target CHEMBL1963 (TSHR)
#     at pChEMBL >= 6 returns 1000/1000 rows of standard_type=Potency,
#     assay_type=F, src_id=7. That panel is what put TSHR, THPO and SMN1 on
#     aspirin and fluorouracil. F2 (CHEMBL204), by contrast, is 779 Ki + 221
#     IC50, all src_id=1 (literature), and is untouched by this rule.
#   assay_type: B (binding) and F (functional) only. 'A' is ADMET -- the CYP
#     inhibition panels that put CYP2C9/CYP3A4 on the dihydropyridines -- and
#     'T' is toxicity. Neither is target engagement.
AFFINITY_TYPES = {"Ki", "Kd", "IC50", "EC50"}
ASSAY_TYPES = {"B", "F"}


def is_affinity(a: dict) -> tuple[bool, str]:
    """(keep?, reason-if-dropped) for one activity record."""
    if a.get("standard_type") not in AFFINITY_TYPES:
        return False, f"standard_type={a.get('standard_type')}"
    if a.get("assay_type") not in ASSAY_TYPES:
        return False, f"assay_type={a.get('assay_type')}"
    rel = a.get("standard_relation")
    if rel not in (None, "", "="):
        return False, f"relation={rel}"
    return True, ""

FIELDNAMES = [
    "struct_id",
    "name",
    "inchikey",
    "gene_symbol",
    "uniprot",
    "organism",
    "evidence",
    "activity_type",
    "activity_value",
    "activity_units",
    "pchembl",
    "passes_threshold",
    "n_activities",
    "assay_types",
    "molecule_chembl_id",
    "target_chembl_id",
    "match_method",
    "source_url",
]


# --------------------------------------------------------------------------
# plumbing: a cached, retrying, fail-soft GET
# --------------------------------------------------------------------------

class Client:
    """Cached ChEMBL client. Any hard failure flips `self.degraded` and returns None."""

    def __init__(self, cache_dir: Path, offline: bool = False):
        self.cache_dir = cache_dir
        self.offline = offline
        self.degraded = False
        self.failures: list[str] = []
        self.n_cached = 0
        self.n_fetched = 0
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "drug-similarity-discovery/enrich_targets"

    def _path(self, stage: str, key: str) -> Path:
        h = hashlib.sha1(key.encode()).hexdigest()[:20]
        return self.cache_dir / stage / f"{h}.json.gz"

    def get(self, stage: str, endpoint: str, params: dict) -> dict | None:
        key = endpoint + "?" + json.dumps(params, sort_keys=True)
        path = self._path(stage, key)
        if path.exists():
            try:
                with gzip.open(path, "rt") as fh:
                    self.n_cached += 1
                    return json.load(fh)
            except (OSError, ValueError):
                path.unlink(missing_ok=True)  # corrupt cache entry; refetch
        if self.offline:
            return None
        data = self._fetch(endpoint, params)
        if data is None:
            return None
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        with gzip.open(tmp, "wt") as fh:
            json.dump(data, fh)
        tmp.replace(path)
        self.n_fetched += 1
        return data

    def _fetch(self, endpoint: str, params: dict) -> dict | None:
        url = f"{CHEMBL}/{endpoint}.json"
        for attempt in range(4):
            try:
                r = self.session.get(url, params=params, timeout=120)
                if r.status_code == 200:
                    return r.json()
                if r.status_code in (429, 500, 502, 503, 504):
                    time.sleep(2 ** attempt)
                    continue
                self.failures.append(f"{endpoint} HTTP {r.status_code}")
                return None
            except (requests.RequestException, ValueError) as exc:
                if attempt == 3:
                    self.failures.append(f"{endpoint} {type(exc).__name__}")
                    self.degraded = True
                    return None
                time.sleep(2 ** attempt)
        return None

    def paginate(self, stage: str, endpoint: str, params: dict, collection: str) -> list | None:
        """Walk every page. Returns None if any page failed (so we never emit a partial set)."""
        out: list = []
        offset = 0
        while True:
            p = dict(params, limit=PAGE, offset=offset)
            data = self.get(stage, endpoint, p)
            if data is None:
                return None
            out.extend(data.get(collection, []))
            meta = data.get("page_meta", {})
            if not meta.get("next"):
                return out
            offset += PAGE
            if offset >= 20000:  # ChEMBL offset ceiling; caller must split the batch
                return out


def chunks(seq, n):
    seq = list(seq)
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# --------------------------------------------------------------------------
# stage 1: drugs -> ChEMBL molecule ids
# --------------------------------------------------------------------------

def match_molecules(client: Client, drugs: list[dict]) -> tuple[dict, dict]:
    """Return (struct_id -> molecule_chembl_id, struct_id -> match_method)."""
    by_ik: dict[str, list[dict]] = defaultdict(list)
    for d in drugs:
        ik = d["inchikey"].strip()
        if ik:
            by_ik[ik].append(d)

    mol_of: dict[str, str] = {}
    method: dict[str, str] = {}

    # -- exact InChIKey ----------------------------------------------------
    batches = list(chunks(sorted(by_ik), 40))
    print(f"[match] exact InChIKey: {len(by_ik)} keys in {len(batches)} batches", flush=True)

    def fetch_ik(batch):
        return batch, client.get(
            "molecule_by_inchikey", "molecule",
            {"molecule_structures__standard_inchi_key__in": ",".join(batch),
             "only": "molecule_chembl_id,molecule_structures", "limit": 100},
        )

    ik_to_mol: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for batch, data in pool.map(fetch_ik, batches):
            if data is None:
                continue
            for m in data.get("molecules", []):
                st = m.get("molecule_structures") or {}
                key = st.get("standard_inchi_key")
                if key:
                    ik_to_mol[key] = m["molecule_chembl_id"]

    for ik, ds in by_ik.items():
        if ik in ik_to_mol:
            for d in ds:
                mol_of[d["struct_id"]] = ik_to_mol[ik]
                method[d["struct_id"]] = "inchikey_exact"

    # -- name fallback, for whatever the InChIKey missed --------------------
    unmatched = [d for d in drugs if d["struct_id"] not in mol_of and d["name"].strip()]
    names = sorted({d["name"].strip().lower() for d in unmatched})
    nbatches = list(chunks(names, 25))
    print(f"[match] name fallback: {len(names)} names in {len(nbatches)} batches", flush=True)

    def fetch_name(batch):
        return client.get(
            "molecule_by_name", "molecule",
            {"pref_name__in": ",".join(b.upper() for b in batch),
             "only": "molecule_chembl_id,pref_name", "limit": 200},
        )

    name_to_mol: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for data in pool.map(fetch_name, nbatches):
            if data is None:
                continue
            for m in data.get("molecules", []):
                pn = (m.get("pref_name") or "").strip().lower()
                if pn and pn not in name_to_mol:
                    name_to_mol[pn] = m["molecule_chembl_id"]

    for d in unmatched:
        mid = name_to_mol.get(d["name"].strip().lower())
        if mid:
            mol_of[d["struct_id"]] = mid
            method[d["struct_id"]] = "name_exact"

    return mol_of, method


# --------------------------------------------------------------------------
# stage 2: molecules -> activities
# --------------------------------------------------------------------------

def fetch_activities(client: Client, molecule_ids: list[str]) -> tuple[list[dict], int]:
    """All human activities at pChEMBL >= FETCH_FLOOR. Returns (rows, n_failed_batches)."""
    batches = list(chunks(sorted(set(molecule_ids)), 10))
    print(f"[activity] {len(set(molecule_ids))} molecules in {len(batches)} batches", flush=True)
    rows: list[dict] = []
    failed = 0
    done = 0

    def fetch(batch):
        return client.paginate(
            "activity", "activity",
            {"molecule_chembl_id__in": ",".join(batch),
             "pchembl_value__gte": FETCH_FLOOR,
             "only": ACTIVITY_FIELDS},
            "activities",
        )

    def fetch_guarded(batch):
        """Split a batch that hits ChEMBL's offset ceiling, so nothing is silently lost."""
        res = fetch(batch)
        if res is not None and len(res) >= 20000 and len(batch) > 1:
            out = []
            for one in batch:
                sub = fetch([one])
                if sub is None:
                    return None
                out.extend(sub)
            return out
        return res

    with ThreadPoolExecutor(max_workers=4) as pool:
        for res in pool.map(fetch_guarded, batches):
            done += 1
            if done % 50 == 0:
                print(f"[activity]   {done}/{len(batches)} batches, {len(rows)} rows", flush=True)
            if res is None:
                failed += 1
                continue
            rows.extend(res)
    return rows, failed


# --------------------------------------------------------------------------
# stage 3: targets -> uniprot + gene symbol
# --------------------------------------------------------------------------

def fetch_targets(client: Client, target_ids: list[str]) -> dict[str, dict]:
    """target_chembl_id -> {gene_symbol, uniprot, organism} for SINGLE PROTEIN targets.

    Both human and non-human are kept. Nearly every approved drug left without a
    target after the human-only pass is an anti-infective whose target is a
    pathogen protein, so restricting to Homo sapiens throws away exactly the
    drugs that most need annotating. Non-human rows are written under a distinct
    `evidence` value so a consumer can take the human tier alone if it wants.
    `uniprot` is always populated and is the safe join key; `gene_symbol` is
    best-effort (many non-human targets have no GENE_SYMBOL synonym).
    """
    ids = sorted(set(target_ids))
    batches = list(chunks(ids, 20))
    print(f"[target] {len(ids)} targets in {len(batches)} batches", flush=True)
    out: dict[str, dict] = {}

    def fetch(batch):
        return client.get("target", "target",
                          {"target_chembl_id__in": ",".join(batch), "limit": 100})

    with ThreadPoolExecutor(max_workers=4) as pool:
        for data in pool.map(fetch, batches):
            if data is None:
                continue
            for t in data.get("targets", []):
                if t.get("target_type") != "SINGLE PROTEIN":
                    continue
                comps = t.get("target_components") or []
                if len(comps) != 1:
                    continue
                c = comps[0]
                if not c.get("accession"):
                    continue
                genes = [s["component_synonym"] for s in c.get("target_component_synonyms", [])
                         if s.get("syn_type") == "GENE_SYMBOL"]
                organism = t.get("organism") or ""
                if organism == "Homo sapiens" and not genes:
                    continue  # human rows must carry an HGNC symbol to join cleanly
                out[t["target_chembl_id"]] = {
                    "gene_symbol": sorted(genes)[0] if genes else "",
                    "uniprot": c["accession"],
                    "organism": organism,
                }
    return out


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

def build_rows(drugs, mol_of, method, activities, targets) -> tuple[list[dict], Counter]:
    by_struct = {d["struct_id"]: d for d in drugs}
    # a molecule id can back more than one struct_id (salt forms share a parent)
    structs_of_mol: dict[str, list[str]] = defaultdict(list)
    for sid, mid in mol_of.items():
        structs_of_mol[mid].append(sid)

    # aggregate: (struct, gene) -> strongest measurement + provenance
    agg: dict[tuple, dict] = {}
    drops = Counter()
    for a in activities:
        tinfo = targets.get(a.get("target_chembl_id"))
        if not tinfo:
            drops["target not a single protein (cell line / complex / no UniProt)"] += 1
            continue
        keep, why = is_affinity(a)
        if not keep:
            drops[why] += 1
            continue
        try:
            pch = float(a["pchembl_value"])
        except (TypeError, ValueError, KeyError):
            continue
        for sid in structs_of_mol.get(a["molecule_chembl_id"], []):
            key = (sid, tinfo["organism"], tinfo["gene_symbol"] or tinfo["uniprot"])
            cur = agg.get(key)
            if cur is None:
                cur = {
                    "best": pch, "act": a, "n": 0, "assays": set(),
                    "tid": a["target_chembl_id"], "uniprot": tinfo["uniprot"],
                    "organism": tinfo["organism"], "mid": a["molecule_chembl_id"],
                }
                agg[key] = cur
            cur["n"] += 1
            if a.get("assay_type"):
                cur["assays"].add(a["assay_type"])
            if pch > cur["best"]:
                cur["best"] = pch
                cur["act"] = a

    rows: list[dict] = []
    for (sid, _organism, gene), v in sorted(agg.items()):
        d = by_struct[sid]
        a = v["act"]
        rows.append({
            "struct_id": sid,
            "name": d["name"],
            "inchikey": d["inchikey"],
            "gene_symbol": gene,
            "uniprot": v["uniprot"],
            "organism": v["organism"],
            "evidence": ("chembl_pchembl" if v["organism"] == "Homo sapiens"
                         else "chembl_pchembl_nonhuman"),
            "activity_type": a.get("standard_type") or "",
            "activity_value": a.get("standard_value") or "",
            "activity_units": a.get("standard_units") or "",
            "pchembl": f"{v['best']:.2f}",
            "passes_threshold": 1 if v["best"] >= DEFAULT_THRESHOLD else 0,
            "n_activities": v["n"],
            "assay_types": "".join(sorted(v["assays"])),
            "molecule_chembl_id": v["mid"],
            "target_chembl_id": v["tid"],
            "match_method": method.get(sid, ""),
            "source_url": (f"{CHEMBL}/activity.json?molecule_chembl_id={v['mid']}"
                           f"&target_chembl_id={v['tid']}&pchembl_value__isnull=false"),
        })

    # DrugCentral curated layer, carried through so the file is self-sufficient
    uniprot_of_gene = {r["gene_symbol"]: r["uniprot"] for r in rows
                       if r["uniprot"] and r["organism"] == "Homo sapiens"}
    for d in drugs:
        moa = {t for t in d["moa_targets"].split(";") if t}
        for t in [x for x in d["targets"].split(";") if x]:
            rows.append({
                "struct_id": d["struct_id"],
                "name": d["name"],
                "inchikey": d["inchikey"],
                "gene_symbol": t,
                "uniprot": uniprot_of_gene.get(t, ""),
                "organism": "Homo sapiens",
                "evidence": "drugcentral_moa" if t in moa else "drugcentral_target",
                "activity_type": "", "activity_value": "", "activity_units": "",
                "pchembl": "", "passes_threshold": 1,
                "n_activities": "", "assay_types": "",
                "molecule_chembl_id": mol_of.get(d["struct_id"], ""),
                "target_chembl_id": "",
                "match_method": "drugcentral_curated",
                "source_url": f"https://drugcentral.org/drugcard/{d['struct_id']}",
            })
    return rows, drops


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def report(drugs, rows, mol_of, method, client, failed_batches, activities, drops):
    approved = [d for d in drugs if d["approved"] == "1"]
    out = []

    def p(s=""):
        print(s, flush=True)
        out.append(s)

    n_all, n_appr = len(drugs), len(approved)
    p(f"\n{'='*72}\nTARGET ANNOTATION REPORT   (all rows n={n_all}, approved==1 n={n_appr})\n{'='*72}")

    p("\n-- drug -> ChEMBL match rate (n=%d rows) --" % n_all)
    mc = Counter(method.values())
    for k in ("inchikey_exact", "name_exact"):
        p(f"   {k:<16} {mc.get(k,0):>5}  ({100*mc.get(k,0)/n_all:.1f}%)")
    p(f"   {'unmatched':<16} {n_all-len(mol_of):>5}  ({100*(n_all-len(mol_of))/n_all:.1f}%)")
    am = Counter(method[d["struct_id"]] for d in approved if d["struct_id"] in method)
    p(f"   approved==1 matched: {sum(am.values())}/{n_appr} "
      f"({100*sum(am.values())/n_appr:.1f}%)  [inchikey {am.get('inchikey_exact',0)}, "
      f"name {am.get('name_exact',0)}]")

    p("\n-- assay-quality filter, effect on %d fetched activity rows --" % len(activities))
    kept = len([r for r in rows if r["evidence"].startswith("chembl_")])
    for k, v in drops.most_common(8):
        p(f"   dropped {v:>7}  {k}")
    p(f"   dropped {sum(drops.values()):>7}  TOTAL "
      f"({100*sum(drops.values())/max(len(activities),1):.1f}% of fetched rows)")
    p(f"   kept    {len(activities)-sum(drops.values()):>7}  activity rows -> "
      f"{kept} distinct (drug, target) pairs "
      f"[human {len([r for r in rows if r['evidence']=='chembl_pchembl'])}, "
      f"non-human {len([r for r in rows if r['evidence']=='chembl_pchembl_nonhuman'])}]")

    # coverage at several thresholds
    p("\n-- coverage vs pChEMBL threshold --")
    p("   'human' = human SINGLE PROTEIN only; '+nonhuman' also counts pathogen targets.")
    p(f"   {'thr':>5} | {'drugs':>6} {'resc.all':>9} {'resc.appr':>10} {'rows':>7} {'targets':>8}"
      f" | {'drugs':>6} {'resc.all':>9} {'resc.appr':>10} {'rows':>7} {'targets':>8}")
    p(f"   {'':>5} | {'--- human only ---':^45} | {'--- human + nonhuman ---':^45}")
    base_unann_all = {d["struct_id"] for d in drugs if not d["targets"].strip()}
    base_unann_appr = {d["struct_id"] for d in approved if not d["targets"].strip()}
    ch = [r for r in rows if r["evidence"] == "chembl_pchembl"]
    ch_all = [r for r in rows if r["evidence"].startswith("chembl_")]
    for thr in (5.0, 5.5, 6.0, 6.5, 7.0, 8.0):
        cells = []
        for pool in (ch, ch_all):
            sel = [r for r in pool if float(r["pchembl"]) >= thr]
            dw = {r["struct_id"] for r in sel}
            tg = {r["uniprot"] or r["gene_symbol"] for r in sel}
            cells.append(f"{len(dw):>6} {len(dw & base_unann_all):>9} "
                         f"{len(dw & base_unann_appr):>10} {len(sel):>7} {len(tg):>8}")
        p(f"   {thr:>5.1f} | {cells[0]} | {cells[1]}")

    thr = DEFAULT_THRESHOLD
    sel = [r for r in ch if float(r["pchembl"]) >= thr]
    sel_all = [r for r in ch_all if float(r["pchembl"]) >= thr]
    p(f"\n-- at the chosen threshold pChEMBL >= {thr} (1 uM) --")

    for label, pool, unann in (("all rows, human targets", drugs, base_unann_all),
                               ("approved==1, human targets", approved, base_unann_appr)):
        n = len(pool)
        ids = {d["struct_id"] for d in pool}
        ch_ids = {r["struct_id"] for r in sel} & ids
        old = ids - unann
        new = old | ch_ids
        d2t = defaultdict(set)
        per = defaultdict(set)
        for d in pool:
            for t in [x for x in d["targets"].split(";") if x]:
                d2t[t].add(d["struct_id"]); per[d["struct_id"]].add(t)
        old_tgts = set(d2t)
        for r in sel:
            if r["struct_id"] in ids:
                d2t[r["gene_symbol"]].add(r["struct_id"]); per[r["struct_id"]].add(r["gene_symbol"])
        p(f"\n   [{label}] n={n}")
        p(f"     annotated  before {len(old):>5} ({100*len(old)/n:.1f}%)"
          f"   after {len(new):>5} ({100*len(new)/n:.1f}%)")
        p(f"     rescued (had nothing, now annotated): {len(ch_ids & unann)}")
        p(f"     distinct targets   before {len(old_tgts):>5}   after {len(d2t):>5}")
        pv = [len(v) for v in per.values() if v]
        p(f"     median targets/drug after {statistics.median(pv):.0f}"
          f"   median drugs/target after {statistics.median([len(v) for v in d2t.values()]):.0f}")
        big = sum(1 for v in d2t.values() if len(v) >= 10)
        p(f"     targets with >=10 drugs after: {big}")

    for label, pool, unann in (("all rows, human + nonhuman", drugs, base_unann_all),
                               ("approved==1, human + nonhuman", approved, base_unann_appr)):
        n = len(pool)
        ids = {d["struct_id"] for d in pool}
        ch_ids = {r["struct_id"] for r in sel_all} & ids
        old = ids - unann
        new = old | ch_ids
        d2t = defaultdict(set)
        per = defaultdict(set)
        for d in pool:
            for t in [x for x in d["targets"].split(";") if x]:
                d2t[t].add(d["struct_id"]); per[d["struct_id"]].add(t)
        old_tgts = set(d2t)
        for r in sel_all:
            if r["struct_id"] in ids:
                k = r["uniprot"] or r["gene_symbol"]
                d2t[k].add(r["struct_id"]); per[r["struct_id"]].add(k)
        p(f"\n   [{label}] n={n}")
        p(f"     annotated  before {len(old):>5} ({100*len(old)/n:.1f}%)"
          f"   after {len(new):>5} ({100*len(new)/n:.1f}%)")
        p(f"     rescued (had nothing, now annotated): {len(ch_ids & unann)}")
        p(f"     distinct targets   before {len(old_tgts):>5}   after {len(d2t):>5}")
        pv = [len(v) for v in per.values() if v]
        p(f"     median targets/drug after {statistics.median(pv):.0f}"
          f"   median drugs/target after {statistics.median([len(v) for v in d2t.values()]):.0f}")
        p(f"     targets with >=10 drugs after: {sum(1 for v in d2t.values() if len(v) >= 10)}")

    # F2 specifically
    p("\n-- F2 (thrombin) --")
    appr_ids = {d["struct_id"] for d in approved}
    f2_old = {d["struct_id"] for d in approved if "F2" in d["targets"].split(";")}
    f2_moa = {d["struct_id"] for d in approved if "F2" in d["moa_targets"].split(";")}
    for t in (5.0, 6.0, 7.0):
        f2_new = {r["struct_id"] for r in ch
                  if r["gene_symbol"] == "F2" and float(r["pchembl"]) >= t and r["struct_id"] in appr_ids}
        p(f"   thr {t}: curated {len(f2_old)} (MOA {len(f2_moa)})  chembl {len(f2_new)}  "
          f"union {len(f2_old | f2_new)}  net new {len(f2_new - f2_old)}")

    # agreement with the curated column
    p("\n-- agreement with DrugCentral's curated annotations --")
    cur = {d["struct_id"]: {t for t in d["targets"].split(";") if t} for d in drugs if d["targets"].strip()}
    curmoa = {d["struct_id"]: {t for t in d["moa_targets"].split(";") if t} for d in drugs if d["moa_targets"].strip()}
    new_by = defaultdict(set)
    for r in sel:
        new_by[r["struct_id"]].add(r["gene_symbol"])
    both = [s for s in cur if s in new_by]
    p(f"   drugs with BOTH curated and ChEMBL targets: n={len(both)}")
    if both:
        rec = [len(cur[s] & new_by[s]) / len(cur[s]) for s in both]
        prec = [len(cur[s] & new_by[s]) / len(new_by[s]) for s in both]
        agree_any = sum(1 for s in both if cur[s] & new_by[s])
        p(f"   drugs where the two sets share >=1 target: {agree_any}/{len(both)} "
          f"({100*agree_any/len(both):.1f}%)")
        p(f"   mean recall  of curated targets by ChEMBL: {100*statistics.mean(rec):.1f}%")
        p(f"   mean fraction of ChEMBL targets that are curated: {100*statistics.mean(prec):.1f}%")
        bm = [s for s in curmoa if s in new_by]
        if bm:
            mrec = [len(curmoa[s] & new_by[s]) / len(curmoa[s]) for s in bm]
            hit = sum(1 for s in bm if curmoa[s] & new_by[s])
            p(f"   MOA targets (n={len(bm)} drugs): >=1 MOA target recovered "
              f"{hit}/{len(bm)} ({100*hit/len(bm):.1f}%), mean MOA recall "
              f"{100*statistics.mean(mrec):.1f}%")
        worst = sorted(both, key=lambda s: len(cur[s] & new_by[s]) / len(cur[s]))[:5]
        p("   five drugs with the lowest recall (for manual inspection):")
        name_of = {d["struct_id"]: d["name"] for d in drugs}
        for s in worst:
            p(f"     {name_of[s][:28]:<28} curated={sorted(cur[s])[:6]} chembl={sorted(new_by[s])[:6]}")

    # not-annotated accounting
    p("\n-- still not annotated, and why (approved==1, n=%d) --" % n_appr)
    ch_ids = {r["struct_id"] for r in sel_all}
    buckets = Counter()
    for d in approved:
        sid = d["struct_id"]
        if d["targets"].strip() or sid in ch_ids:
            continue
        if sid not in mol_of:
            buckets["no ChEMBL molecule matched (InChIKey and name both missed)"] += 1
        elif sid in {r["struct_id"] for r in ch_all}:
            buckets[f"ChEMBL activities exist but all below pChEMBL {thr}"] += 1
        else:
            buckets["matched ChEMBL, but zero single-protein pChEMBL activities (any organism)"] += 1
    for k, v in buckets.most_common():
        p(f"   {v:>5}  {k}")
    p(f"   {sum(buckets.values()):>5}  TOTAL still unannotated "
      f"({100*sum(buckets.values())/n_appr:.1f}% of approved)")

    p("\n-- provenance --")
    p(f"   cache hits {client.n_cached}, network fetches {client.n_fetched}, "
      f"failed activity batches {failed_batches}")
    if client.failures:
        p(f"   API failures ({len(client.failures)}): {Counter(client.failures).most_common(5)}")
    if client.degraded:
        p("   WARNING: ChEMBL was degraded during this run; numbers may be incomplete.")
    return "\n".join(out)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drugs", type=Path, default=DRUGS_CSV)
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    ap.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    ap.add_argument("--offline", action="store_true",
                    help="use only what is already cached; never hit the network")
    args = ap.parse_args()

    drugs = list(csv.DictReader(args.drugs.open()))
    print(f"[load] {len(drugs)} rows from {args.drugs}", flush=True)

    client = Client(args.cache_dir, offline=args.offline)
    mol_of, method = match_molecules(client, drugs)
    print(f"[match] {len(mol_of)}/{len(drugs)} matched to a ChEMBL molecule", flush=True)

    activities, failed = fetch_activities(client, list(mol_of.values()))
    print(f"[activity] {len(activities)} rows, {failed} failed batches", flush=True)

    targets = fetch_targets(client, [a["target_chembl_id"] for a in activities])
    n_hu = sum(1 for t in targets.values() if t["organism"] == "Homo sapiens")
    print(f"[target] {len(targets)} single-protein targets resolved "
          f"({n_hu} human, {len(targets)-n_hu} non-human)", flush=True)

    rows, drops = build_rows(drugs, mol_of, method, activities, targets)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    print(f"[write] {len(rows)} rows -> {args.out}", flush=True)

    report(drugs, rows, mol_of, method, client, failed, activities, drops)
    return 0


if __name__ == "__main__":
    sys.exit(main())
