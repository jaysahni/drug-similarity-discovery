"""Pipeline task T4 - the "ideal binder" for a target, from ChEMBL bioactivities.

Resolves a UniProt accession to its human SINGLE PROTEIN ChEMBL target, pulls every
Ki/IC50/EC50/Kd activity on it, applies the potency / assay-quality / data-validity
filters, and writes a potency-ranked shortlist to
results/pipeline/<disease_slug>/ideal_binder.csv with exactly one row flagged
selected_as_ideal_binder.

The "ideal binder" here is the most potent ligand in the *literature*, not a
theoretical optimum - every row carries that caveat in a column, and downstream
affinity drop-off must be read as "vs. best literature ligand".

Every filter step prints its surviving count so the funnel is auditable. Raw API
responses are cached under data/raw/ so reruns are offline-repeatable.

Exits 2 when no ligand meets the <=10 nM + real-assay-id acceptance bar. The CSV is
still written in that case, carrying the finding - a target with no potent ligand is
a result, not a crash.

Usage:
    ./env/bin/python scripts/fetch_ligands.py --target KDR --uniprot P35968
    ./env/bin/python scripts/fetch_ligands.py --target KDR --uniprot P35968 \
        --disease colorectal-cancer --refresh
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import certifi

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
# T1's corpus is every approved structure; drugs.csv keeps only the ones with an
# annotated target, so it under-reports. Prefer T1, fall back, report which.
DRUG_SOURCES = (ROOT / "data" / "approved_drugs.csv", ROOT / "data" / "drugs.csv")
BASE = "https://www.ebi.ac.uk/chembl/api/data"

POTENCY_TYPES = ("Ki", "IC50", "EC50", "Kd")
MAX_NM = 100.0          # hard potency filter
ACCEPT_NM = 10.0        # acceptance threshold from the task spec
MIN_CONFIDENCE = 8      # ChEMBL assay confidence score, verified to exist below
MAX_MW = 600.0          # drug-likeness preference, not a hard filter
SHORTLIST = 50

CAVEAT = (
    "best literature ligand for this target, not a theoretical optimum; "
    "potency is the single most potent qualifying ChEMBL activity"
)

# The python.org 3.14 framework build ships no CA bundle, so urllib fails SSL
# verification even though curl succeeds. Use certifi's explicitly.
SSL_CTX = ssl.create_default_context(cafile=certifi.where())


def get(path, params, retries=4):
    url = f"{BASE}/{path}.json?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90, context=SSL_CTX) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return None


def fetch_all(path, key, params, cache, refresh=False):
    """Page a ChEMBL collection endpoint to exhaustion, caching the result."""
    if cache.exists() and not refresh:
        rows = json.loads(cache.read_text())
        print(f"  {key}: {len(rows)} from cache {cache.name}", file=sys.stderr)
        return rows
    rows, offset, limit, total = [], 0, 1000, None
    while True:
        page = get(path, {**params, "limit": limit, "offset": offset})
        got = page[key]
        rows.extend(got)
        total = page["page_meta"]["total_count"]
        print(f"  {key}: {len(rows)}/{total}", file=sys.stderr)
        offset += limit
        if not got or offset >= total:
            break
    # Only cache a complete page-through. A truncated fetch (HTTP 200 with a short
    # or empty page) would otherwise be replayed forever as a real empty result.
    if len(rows) != total:
        print(f"  WARNING: got {len(rows)} of {total} {key} - not caching",
              file=sys.stderr)
        return rows
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(rows))
    return rows


def resolve_target(uniprot, refresh=False):
    """UniProt accession -> the human SINGLE PROTEIN ChEMBL target id.

    Returns (target_id, all_candidates). Non-SINGLE_PROTEIN hits (selectivity
    groups, complexes, chimeras) are reported and excluded: their activities are
    not attributable to this protein alone.
    """
    cache = RAW / f"chembl_targets_{uniprot}.json"
    cands = fetch_all(
        "target",
        "targets",
        {"target_components__accession": uniprot},
        cache,
        refresh,
    )
    single = [
        t for t in cands
        if t["target_type"] == "SINGLE PROTEIN" and t["organism"] == "Homo sapiens"
    ]
    if not single:
        sys.exit(f"no human SINGLE PROTEIN ChEMBL target for {uniprot}")
    # Deterministic pick: lowest numeric ChEMBL id, so the result never depends on
    # API result ordering. More than one is unexpected, so say so.
    single.sort(key=lambda t: int(t["target_chembl_id"][6:]))
    if len(single) > 1:
        print(f"  WARNING: {len(single)} human SINGLE PROTEIN targets for {uniprot} "
              f"({[t['target_chembl_id'] for t in single]}); taking the lowest id",
              file=sys.stderr)
    return single[0]["target_chembl_id"], cands


def approved_index():
    """InChIKey -> name for the approved-drug reference table.

    Returns (exact, skeleton, source). The skeleton map is keyed on the first
    InChIKey block (connectivity only) so a salt or stereoisomer of an approved drug
    still matches - reported as a separate, weaker match type.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    src = next((p for p in DRUG_SOURCES if p.exists()), None)
    if src is None:
        sys.exit(f"no approved-drug reference table; looked for {DRUG_SOURCES}")

    exact, skeleton, n, bad = {}, {}, 0, 0
    with src.open() as fh:
        for r in csv.DictReader(fh):
            n += 1
            key = (r.get("inchikey") or "").strip()  # T1's table precomputes it
            if not key and r.get("smiles"):
                m = Chem.MolFromSmiles(r["smiles"])
                key = Chem.MolToInchiKey(m) if m is not None else ""
            if not key:
                bad += 1
                continue
            exact.setdefault(key, r["name"])
            skeleton.setdefault(key.split("-")[0], r["name"])
    print(f"  approved reference: {src.name}, {n} rows, {len(exact)} InChIKeys, "
          f"{bad} without a usable InChIKey", file=sys.stderr)
    return exact, skeleton, src.name


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="KDR", help="target symbol, for labelling")
    ap.add_argument("--uniprot", default="P35968", help="UniProt accession")
    ap.add_argument("--disease", default="colorectal-cancer", help="output slug")
    ap.add_argument("--refresh", action="store_true", help="ignore data/raw cache")
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")

    print(f"target {args.target} / UniProt {args.uniprot}", file=sys.stderr)
    tid, cands = resolve_target(args.uniprot, args.refresh)
    dropped = [c for c in cands if c["target_chembl_id"] != tid]
    print(f"  ChEMBL target: {tid} "
          f"({len(dropped)} other ChEMBL targets map to {args.uniprot}, excluded: "
          f"{sorted({c['target_type'] for c in dropped})})", file=sys.stderr)

    acts = fetch_all(
        "activity",
        "activities",
        {
            "target_chembl_id": tid,
            "standard_type__in": ",".join(POTENCY_TYPES),
            "only": ",".join([
                "activity_id", "molecule_chembl_id", "molecule_pref_name",
                "canonical_smiles", "standard_type", "standard_value",
                "standard_units", "standard_relation", "pchembl_value",
                "assay_chembl_id", "document_chembl_id", "data_validity_comment",
                "target_organism", "potential_duplicate",
            ]),
        },
        RAW / f"chembl_activities_{tid}.json",
        args.refresh,
    )

    assays = fetch_all(
        "assay",
        "assays",
        {"target_chembl_id": tid,
         "only": "assay_chembl_id,confidence_score,assay_type,description"},
        RAW / f"chembl_assays_{tid}.json",
        args.refresh,
    )
    conf = {a["assay_chembl_id"]: a.get("confidence_score") for a in assays}
    desc = {a["assay_chembl_id"]: a.get("description") or "" for a in assays}

    # Filter on assay confidence only if the field exists. Probe an unfiltered
    # /activity record for it rather than inspecting `acts`, whose fields are
    # restricted by the `only` list above and so could not show it either way.
    # Cached like every other call, so a rerun needs no network (see module docstring).
    probe_cache = RAW / f"chembl_activity_fieldprobe_{tid}.json"
    if probe_cache.exists() and not args.refresh:
        probe = json.loads(probe_cache.read_text())
    else:
        probe = get("activity", {"target_chembl_id": tid, "limit": 1})["activities"]
        probe_cache.parent.mkdir(parents=True, exist_ok=True)
        probe_cache.write_text(json.dumps(probe))
    if not probe:
        conf_on_activity = "no activity to probe"
    else:
        conf_on_activity = "confidence_score" in probe[0]
    n_conf = sum(1 for v in conf.values() if v is not None)
    print(f"\nconfidence_score on an unfiltered /activity record: "
          f"{conf_on_activity}", file=sys.stderr)
    print(f"confidence_score on /assay records: {n_conf}/{len(assays)}",
          file=sys.stderr)
    if n_conf == 0:
        print("confidence_score absent everywhere - SKIPPING that filter",
              file=sys.stderr)
    else:
        print("-> filtering on it via the assay join", file=sys.stderr)

    # ---------------------------------------------------------------- funnel
    funnel = [("activities returned (Ki/IC50/EC50/Kd, target %s)" % tid, len(acts))]

    step = [a for a in acts if a.get("target_organism") == "Homo sapiens"]
    funnel.append(("target_organism == Homo sapiens", len(step)))

    step = [a for a in step if a.get("standard_value") and a.get("standard_units")]
    funnel.append(("has standard_value and standard_units", len(step)))

    step = [a for a in step if a["standard_units"] == "nM"]
    funnel.append(("standard_units == nM", len(step)))

    step = [a for a in step if not a.get("data_validity_comment")]
    funnel.append(("data_validity_comment empty", len(step)))

    step = [a for a in step if a.get("standard_relation") in ("=", "<", "<=")]
    funnel.append(("standard_relation in (=, <, <=)", len(step)))

    # A curated standard_value of 0 is not a real potency - it would otherwise sort
    # to rank 1 and become the ideal binder.
    step = [a for a in step if 0 < float(a["standard_value"]) <= MAX_NM]
    funnel.append((f"potency in (0, {MAX_NM:g}] nM", len(step)))

    if n_conf:
        # Split the drop: a null score or an assay missing from the assay table is
        # unknown quality, not a low score. Both are dropped, counted separately.
        missing = [a for a in step if a["assay_chembl_id"] not in conf]
        null = [a for a in step if conf.get(a["assay_chembl_id"], 0) is None]
        step = [a for a in step
                if isinstance(conf.get(a["assay_chembl_id"]), int)
                and conf[a["assay_chembl_id"]] >= MIN_CONFIDENCE]
        funnel.append((f"assay confidence_score >= {MIN_CONFIDENCE} "
                       f"({len(missing)} dropped for no assay-table join, "
                       f"{len(null)} for a null score)", len(step)))

    step = [a for a in step if a.get("canonical_smiles")]
    funnel.append(("has canonical_smiles", len(step)))

    parsed, no_key = [], 0
    for a in step:
        m = Chem.MolFromSmiles(a["canonical_smiles"])
        if m is None:
            continue
        a["_mw"] = Descriptors.MolWt(m)
        a["_inchikey"] = Chem.MolToInchiKey(m) or ""
        if not a["_inchikey"]:
            no_key += 1  # kept, but it can never match an approved drug
        parsed.append(a)
    funnel.append(("SMILES parses under RDKit", len(parsed)))
    if no_key:
        funnel.append(("of those, yield no InChIKey (kept, never approved-matched)",
                       no_key))

    druglike = [a for a in parsed if a["_mw"] < MAX_MW]
    funnel.append((f"MW < {MAX_MW:g} (preference, not a hard filter)", len(druglike)))

    print("\nfilter funnel", file=sys.stderr)
    for label, count in funnel:
        print(f"  {count:7d}  {label}", file=sys.stderr)

    # Reported, not filtered: a ChEMBL potential_duplicate flag marks an activity
    # re-extracted from another paper, which is redundant evidence rather than bad
    # data. It inflates n_potent_activities, so the count is stated here.
    n_dup = sum(1 for a in parsed if a.get("potential_duplicate"))
    print(f"  ({n_dup} of the {len(parsed)} surviving activities carry ChEMBL's "
          f"potential_duplicate flag - reported, not filtered)", file=sys.stderr)

    if not parsed:
        print("\nFINDING: no activity on this target survives the filters.",
              file=sys.stderr)

    # ------------------------------------------------- one row per molecule
    per_mol = {}
    for a in parsed:
        v = float(a["standard_value"])
        mid = a["molecule_chembl_id"]
        cur = per_mol.get(mid)
        if cur is None:
            per_mol[mid] = {"best": a, "values": [v], "docs": {a["document_chembl_id"]},
                            "assays": {a["assay_chembl_id"]}}
        else:
            cur["values"].append(v)
            cur["docs"].add(a["document_chembl_id"])
            cur["assays"].add(a["assay_chembl_id"])
            if v < float(cur["best"]["standard_value"]):
                cur["best"] = a
    print(f"\n{len(per_mol)} distinct molecules survive all hard filters",
          file=sys.stderr)

    # Rank by potency as specified, drug-like first, then by supporting evidence.
    ranked = sorted(
        per_mol.values(),
        key=lambda e: (
            e["best"]["_mw"] >= MAX_MW,
            float(e["best"]["standard_value"]),
            -len(e["values"]),
        ),
    )
    shortlist = ranked[:SHORTLIST]

    # PubMed ids for the shortlist's documents only - keeps the extra calls cheap.
    doc_ids = sorted({e["best"]["document_chembl_id"] for e in shortlist
                      if e["best"].get("document_chembl_id")})
    pmid = {}
    if doc_ids:
        # Cache by document id, not by target: changing SHORTLIST or a filter changes
        # which documents are needed, and only the missing ones are fetched.
        cache = RAW / f"chembl_documents_{tid}.json"
        docs = {}
        if cache.exists() and not args.refresh:
            docs = {d["document_chembl_id"]: d for d in json.loads(cache.read_text())}
        want = [d for d in doc_ids if d not in docs]
        for i in range(0, len(want), 50):
            page = get("document", {
                "document_chembl_id__in": ",".join(want[i:i + 50]),
                "limit": 50,
                "only": "document_chembl_id,pubmed_id,doi,year,journal",
            })
            for d in page["documents"]:
                docs[d["document_chembl_id"]] = d
        if want:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(list(docs.values())))
        pmid = {k: d.get("pubmed_id") for k, d in docs.items()}
        hit = [d for d in doc_ids if pmid.get(d)]
        print(f"  documents: {len(doc_ids)} needed, {len(want)} fetched now, "
              f"{len(hit)} with a PubMed id", file=sys.stderr)

    exact, skeleton, approved_src = approved_index()

    out_dir = ROOT / "results" / "pipeline" / slugify(args.disease)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "ideal_binder.csv"

    cols = ["rank", "selected_as_ideal_binder", "chembl_id", "pref_name", "smiles",
            "standard_type", "standard_value_nm", "pchembl", "assay_chembl_id",
            "assay_confidence_score", "assay_description", "document_chembl_id",
            "pmid", "mw", "is_approved_drug", "approved_drug_name",
            "approved_match_type", "approved_reference_table",
            "n_potent_activities", "n_distinct_assays",
            "n_distinct_documents", "best_of_molecule_nm", "median_nm",
            "target_symbol", "target_uniprot", "target_chembl_id", "caveat"]

    rows = []
    for i, e in enumerate(shortlist, 1):
        a = e["best"]
        vals = sorted(e["values"])
        med = vals[len(vals) // 2] if len(vals) % 2 else (
            vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2
        key = a["_inchikey"]
        if key and key in exact:
            match, name = "inchikey_exact", exact[key]
        elif key and key.split("-")[0] in skeleton:
            match, name = "inchikey_skeleton", skeleton[key.split("-")[0]]
        else:
            match, name = "none", ""
        rows.append({
            "rank": i,
            "selected_as_ideal_binder": i == 1,
            "chembl_id": a["molecule_chembl_id"],
            "pref_name": a.get("molecule_pref_name") or "",
            "smiles": a["canonical_smiles"],
            "standard_type": a["standard_type"],
            "standard_value_nm": float(a["standard_value"]),
            "pchembl": a.get("pchembl_value") or "",
            "assay_chembl_id": a["assay_chembl_id"],
            "assay_confidence_score": conf.get(a["assay_chembl_id"], ""),
            "assay_description": desc.get(a["assay_chembl_id"], ""),
            "document_chembl_id": a.get("document_chembl_id") or "",
            "pmid": pmid.get(a.get("document_chembl_id")) or "",
            "mw": round(a["_mw"], 2),
            "is_approved_drug": match != "none",
            "approved_drug_name": name,
            "approved_match_type": match,
            "approved_reference_table": approved_src,
            "n_potent_activities": len(e["values"]),
            "n_distinct_assays": len(e["assays"]),
            "n_distinct_documents": len(e["docs"]),
            "best_of_molecule_nm": vals[0],
            "median_nm": med,
            "target_symbol": args.target,
            "target_uniprot": args.uniprot,
            "target_chembl_id": tid,
            "caveat": CAVEAT,
        })

    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        if rows:
            w.writerows(rows)
        else:
            # A target with no qualifying ligand is a finding, not an empty file.
            w.writerow({c: "" for c in cols} | {
                "rank": 0,
                "selected_as_ideal_binder": False,
                "target_symbol": args.target,
                "target_uniprot": args.uniprot,
                "target_chembl_id": tid,
                "caveat": f"FINDING: no ChEMBL activity on {tid} passed the filters "
                          f"({len(acts)} Ki/IC50/EC50/Kd activities examined). "
                          f"No ideal binder exists for this target under these "
                          f"criteria - none was invented.",
            })

    # ---------------------------------------------------------------- report
    passing = [r for r in rows if r["standard_value_nm"] <= ACCEPT_NM
               and r["assay_chembl_id"]]
    print(f"\nACCEPTANCE (>=1 ligand <= {ACCEPT_NM:g} nM with a real assay id): "
          f"{'PASS' if passing else 'FAIL'} - {len(passing)}/{len(rows)} shortlisted "
          f"rows qualify")

    print(f"\ntop 10 of {len(per_mol)} molecules by potency")
    head = f"{'#':>3} {'chembl_id':<14} {'type':<5} {'nM':>8} {'pchembl':>7} " \
           f"{'assay':<15} {'conf':>4} {'mw':>7} {'n':>3}  approved"
    print(head)
    print("-" * len(head))
    for r in rows[:10]:
        print(f"{r['rank']:>3} {r['chembl_id']:<14} {r['standard_type']:<5} "
              f"{r['standard_value_nm']:>8.4g} {str(r['pchembl']):>7} "
              f"{r['assay_chembl_id']:<15} {str(r['assay_confidence_score']):>4} "
              f"{r['mw']:>7.1f} {r['n_potent_activities']:>3}  "
              f"{r['approved_drug_name'] or '-'}")

    n_appr = sum(1 for r in rows if r["is_approved_drug"])
    print(f"\n{n_appr}/{len(rows)} shortlisted ligands match an approved drug in "
          f"{approved_src} by InChIKey")
    print(f"caveat: {CAVEAT}")
    print(f"\nwrote {out} ({len(rows)} rows)")
    if not passing:
        sys.exit(2)  # finding, not a crash - the CSV above records it


if __name__ == "__main__":
    main()
