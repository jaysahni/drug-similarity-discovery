"""Pipeline task T2 - disease name to a ranked list of associated targets.

Resolves a free-text disease name to an ontology id with the Open Targets search
endpoint (never hardcoded), pulls `disease.associatedTargets` from the Open
Targets GraphQL API, and writes the top N to
results/pipeline/<disease_slug>/targets.json with the overall association score,
every datatype sub-score, the UniProt accession, and per-target drug counts.

The drug counts exist for the repo's circularity guard: an association driven by
the known-drug datatype is partly "a drug already hits this target", so a
downstream step that wants a drug-naive target list can filter on
`has_known_drug_evidence` / `known_drug_score`.

The known-drug datatype id is NOT assumed - Open Targets renamed it between data
releases, so the script looks for every known spelling and records which one the
live API actually returned (`known_drug_datatype_ids` in the meta block).

Raw API responses are cached under data/raw/opentargets/ so reruns are cheap and
offline-repeatable.

Usage:
    ./env/bin/python scripts/disease_targets.py --disease "colorectal cancer"
    ./env/bin/python scripts/disease_targets.py --disease "pancreatic cancer" --top 25
    ./env/bin/python scripts/disease_targets.py --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
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
RAW = ROOT / "data" / "raw" / "opentargets"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"

# The python.org 3.14 framework build ships no CA bundle, so urllib fails SSL
# verification even though curl succeeds. Use certifi's explicitly.
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

# Datatypes the pipeline plan asks to break out. Open Targets renamed
# `known_drug` -> `clinical` at some point before data release 26.06; both
# spellings are probed and whichever the API returns is recorded in the meta.
KNOWN_DRUG_IDS = ("known_drug", "clinical")
REQUESTED_DATATYPES = (
    "genetic_association",
    "known_drug",
    "literature",
    "rna_expression",
    "somatic_mutation",
    "animal_model",
    "affected_pathway",
)


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def gql(query, cache_key, refresh=False, retries=4):
    """POST a GraphQL query, caching the parsed response under data/raw/.

    The key carries a hash of the query text: editing a selection set must not
    replay a cached payload that lacks the newly requested fields.
    """
    cache = RAW / f"{cache_key}-{hashlib.sha1(query.encode()).hexdigest()[:10]}.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())

    body = json.dumps({"query": query}).encode()
    req = urllib.request.Request(
        OT, data=body, headers={"Content-Type": "application/json"}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120, context=SSL_CTX) as r:
                payload = json.load(r)
            break
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    if payload.get("errors"):
        raise RuntimeError(f"GraphQL errors for {cache_key}: {payload['errors']}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(payload, indent=1))
    return payload


def api_meta(refresh):
    q = "{ meta { name apiVersion { x y z } dataVersion { year month iteration } } }"
    m = gql(q, "meta", refresh)["data"]["meta"]
    av, dv = m["apiVersion"], m["dataVersion"]
    return {
        "endpoint": OT,
        "api_name": m["name"],
        "api_version": f"{av['x']}.{av['y']}.{av['z']}",
        "data_version": f"{dv['year']}.{dv['month']}"
        + (f".{dv['iteration']}" if dv.get("iteration") is not None else ""),
    }


def resolve_disease(name, refresh):
    """Free text -> ontology id via the search endpoint. Never hardcoded."""
    esc = json.dumps(name)
    q = (
        f"{{ search(queryString: {esc}, entityNames: [\"disease\"], "
        "page: {index:0, size:10}) { total hits { id name entity score } } }"
    )
    # gql() hashes the query text into the key, so "non-small cell" and
    # "non small cell" cannot share a cache entry despite slugifying the same.
    s = gql(q, f"search-{slugify(name)}", refresh)["data"]["search"]
    hits = s["hits"]
    if not hits:
        raise SystemExit(f"No disease hit for {name!r}")

    # Prefer an exact (case-insensitive) name match; otherwise the top-scoring
    # hit. Record which rule fired so the choice is auditable.
    exact = [h for h in hits if h["name"].lower() == name.lower()]
    chosen, rule = (exact[0], "exact_name_match") if exact else (hits[0], "top_score")
    return chosen, rule, hits, s["total"]


def associated_targets(efo_id, size, refresh):
    q = (
        f"{{ disease(efoId: \"{efo_id}\") {{ id name "
        f"associatedTargets(page: {{index:0, size:{size}}}) {{ count rows {{ "
        "score datatypeScores { id score } "
        "target { id approvedSymbol approvedName biotype "
        "symbolSynonyms { label source } proteinIds { id source } } } } } }"
    )
    d = gql(q, f"assoc-{efo_id}-{size}", refresh)["data"]["disease"]
    if d is None:
        raise SystemExit(f"Open Targets has no disease {efo_id}")
    return d


def drug_candidates(ensembl_id, refresh):
    """Per-target drug + clinical-candidate rows (the v26 replacement for knownDrugs)."""
    q = (
        f"{{ target(ensemblId: \"{ensembl_id}\") {{ approvedSymbol "
        "drugAndClinicalCandidates { count rows { maxClinicalStage "
        "drug { id name drugType } } } } }"
    )
    t = gql(q, f"drugs-{ensembl_id}", refresh)["data"]["target"]
    c = (t or {}).get("drugAndClinicalCandidates") or {}
    count = c.get("count") or 0
    rows = c.get("rows") or []
    approved = [r for r in rows if r.get("maxClinicalStage") == "APPROVAL"]
    return {
        "n_drug_and_clinical_candidates": count,
        # The field takes no page argument, so rows may in principle be a page of
        # count. Record both rather than let n_approved_drugs quietly undercount.
        "n_drug_rows_returned": len(rows),
        "drug_rows_complete": len(rows) == count,
        "n_approved_drugs": len(approved),
        "approved_drug_examples": sorted(
            {r["drug"]["name"] for r in approved if r.get("drug")}
        )[:8],
    }


def uniprot_from_proteinids(protein_ids):
    """Swiss-Prot accession if Open Targets carries one."""
    for p in protein_ids or []:
        if p["source"] == "uniprot_swissprot":
            return p["id"], "opentargets_proteinIds:uniprot_swissprot"
    return None, None


def uniprot_rest(symbol, refresh):
    """Fallback: reviewed human entry for a gene symbol, via UniProt REST."""
    cache = RAW / f"uniprot-{symbol}.json"
    if cache.exists() and not refresh:
        payload = json.loads(cache.read_text())
    else:
        params = {
            "query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
            "fields": "accession",
            "format": "json",
            "size": 1,
        }
        url = f"{UNIPROT}?{urllib.parse.urlencode(params)}"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60, context=SSL_CTX) as r:
                    payload = json.load(r)
                break
            except Exception as exc:  # noqa: BLE001 - transient network
                if attempt == 3:
                    # A missing accession is a gap to record, not a dead run.
                    print(f"  uniprot lookup failed for {symbol}: {exc}", file=sys.stderr)
                    return None, f"uniprot_rest_failed:{type(exc).__name__}"
                print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
                time.sleep(2 * (attempt + 1))
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=1))
    res = payload.get("results") or []
    return (res[0]["primaryAccession"], "uniprot_rest") if res else (None, None)


def build_target(rank, row, refresh):
    tgt = row["target"]
    scores = {d["id"]: d["score"] for d in row["datatypeScores"]}

    known_id = next((k for k in KNOWN_DRUG_IDS if k in scores), None)
    known_score = scores.get(known_id) if known_id else None

    uni, uni_src = uniprot_from_proteinids(tgt.get("proteinIds"))
    if uni is None:
        uni, uni_src = uniprot_rest(tgt["approvedSymbol"], refresh)

    drugs = drug_candidates(tgt["id"], refresh)

    # Requested breakout, null where this data release carries no such datatype.
    breakout = {k: scores.get(k) for k in REQUESTED_DATATYPES}
    if breakout["known_drug"] is None:
        breakout["known_drug"] = known_score

    syns = sorted({s["label"] for s in (tgt.get("symbolSynonyms") or [])})
    return {
        "rank": rank,
        "ensembl_id": tgt["id"],
        "approved_symbol": tgt["approvedSymbol"],
        "approved_name": tgt["approvedName"],
        "biotype": tgt["biotype"],
        "symbol_synonyms": syns,
        "uniprot": uni,
        "uniprot_source": uni_src,
        "overall_association_score": row["score"],
        "datatype_scores": breakout,
        "all_datatype_scores": scores,
        "known_drug_datatype_id": known_id,
        "known_drug_score": known_score,
        # A real 0.0 is "scored, no evidence"; None is "this release has no such
        # datatype". Only the first is a known-drug score, so test for presence.
        "has_known_drug_evidence": known_score is not None and known_score > 0,
        **drugs,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease", default="colorectal cancer")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--refresh", action="store_true", help="ignore data/raw cache")
    args = ap.parse_args()

    meta = api_meta(args.refresh)
    print(
        f"Open Targets {meta['api_name']} api {meta['api_version']} "
        f"data {meta['data_version']}",
        file=sys.stderr,
    )

    hit, rule, hits, n_hits = resolve_disease(args.disease, args.refresh)
    print(
        f"Resolved {args.disease!r} -> {hit['id']} {hit['name']!r} "
        f"(score {hit['score']:.1f}, rule {rule}, {n_hits} disease hits)",
        file=sys.stderr,
    )

    assoc = associated_targets(hit["id"], args.top, args.refresh)
    rows = assoc["associatedTargets"]["rows"]
    print(
        f"{assoc['associatedTargets']['count']} associated targets; "
        f"taking top {len(rows)}",
        file=sys.stderr,
    )

    targets = []
    for i, row in enumerate(rows, 1):
        t = build_target(i, row, args.refresh)
        targets.append(t)
        print(
            f"  {i:2d} {t['approved_symbol']:<9} {t['overall_association_score']:.4f} "
            f"uniprot={t['uniprot']} drugs={t['n_drug_and_clinical_candidates']}"
            f"/approved={t['n_approved_drugs']}",
            file=sys.stderr,
        )

    observed = sorted({k for t in targets for k in t["all_datatype_scores"]})
    known_ids = sorted({t["known_drug_datatype_id"] for t in targets} - {None})
    # "known_drug" counts as present if any of its accepted spellings turned up.
    missing = [
        k
        for k in REQUESTED_DATATYPES
        if k not in observed and not (k == "known_drug" and known_ids)
    ]
    if known_ids == ["known_drug"]:
        kd_note = (
            "The known-drug datatype is returned under its historical id "
            "'known_drug' in this response."
        )
    elif "known_drug" in known_ids:
        kd_note = (
            f"This response mixes known-drug datatype ids {known_ids!r} across "
            "targets. known_drug_score takes whichever id a given target carried; "
            "the per-target known_drug_datatype_id says which."
        )
    elif known_ids:
        kd_note = (
            f"Open Targets data release {meta['data_version']} returns the "
            f"known-drug datatype under id(s) {known_ids!r}, not 'known_drug'. "
            "Established by inspecting every datatype id in this response rather "
            "than assuming a name; see datatype_ids_observed_in_these_targets."
        )
    else:
        kd_note = (
            "No known-drug datatype under any accepted spelling "
            f"{list(KNOWN_DRUG_IDS)!r} appears among these "
            f"{len(targets)} targets, so known_drug_score is null throughout. "
            "This is an observation about these targets, not proof that the "
            "release lacks the datatype."
        )

    # Acceptance: a top-10 target with literature evidence AND at least one drug.
    top10 = targets[:10]
    passing = [
        t
        for t in top10
        if (t["datatype_scores"]["literature"] or 0) > 0
        and t["n_drug_and_clinical_candidates"] >= 1
    ]

    out = {
        "meta": {
            "script": "scripts/disease_targets.py",
            "generated_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "disease_query": args.disease,
            "efo_id": hit["id"],
            "disease_name": hit["name"],
            "search_hit_chosen": hit,
            "search_hit_rule": rule,
            "search_hits_considered": hits,
            "search_total_disease_hits": n_hits,
            "n_associated_targets_total": assoc["associatedTargets"]["count"],
            "n_targets_written": len(targets),
            "known_drug_datatype_ids": known_ids,
            "known_drug_datatype_note": kd_note,
            # Scoped to the targets actually written, not the whole release.
            "datatype_ids_observed_in_these_targets": observed,
            "requested_datatypes_absent_from_these_targets": missing,
            "acceptance_top10_literature_and_known_drug": [
                t["approved_symbol"] for t in passing
            ],
            **meta,
        },
        "targets": targets,
    }

    out_dir = ROOT / "results" / "pipeline" / slugify(args.disease)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "targets.json"
    path.write_text(json.dumps(out, indent=1) + "\n")

    print(f"\nWrote {len(targets)} targets to {path}", file=sys.stderr)
    print(f"datatype ids observed in these targets: {observed}", file=sys.stderr)
    if missing:
        print(f"requested datatypes absent from these targets: {missing}", file=sys.stderr)
    print(
        f"ACCEPTANCE top-10 with literature>0 and >=1 drug: "
        f"{len(passing)}/{len(top10)} "
        f"{[t['approved_symbol'] for t in passing]}",
        file=sys.stderr,
    )

    print(f"\nTop 10 targets for {hit['name']} ({hit['id']}):")
    hdr = f"{'#':>2}  {'symbol':<9} {'uniprot':<8} {'overall':>7} {'lit':>6} {'known_drug':>10} {'drugs':>6} {'appr':>5}"
    print(hdr)
    print("-" * len(hdr))
    for t in top10:
        ds = t["datatype_scores"]
        lit = f"{ds['literature']:.3f}" if ds["literature"] is not None else "-"
        kd = f"{ds['known_drug']:.3f}" if ds["known_drug"] is not None else "-"
        print(
            f"{t['rank']:>2}  {t['approved_symbol']:<9} {t['uniprot'] or '-':<8} "
            f"{t['overall_association_score']:>7.4f} {lit:>6} {kd:>10} "
            f"{t['n_drug_and_clinical_candidates']:>6} {t['n_approved_drugs']:>5}"
        )

    if not passing:
        raise SystemExit(
            "ACCEPTANCE FAILED: no top-10 target has literature>0 and >=1 known drug"
        )


if __name__ == "__main__":
    main()
