"""Stage 0 of the demo pipeline: disease -> ranked, evidenced, *usable* targets.

`demo/pipeline.py` starts from a hardcoded `TARGETS` dict (CDK2, thrombin) that a
human picked. This module is the step that was missing in front of it: given a
disease name or an EFO/MONDO id it returns a ranked list of candidate targets,
each carrying the evidence for the pick and - the part that actually decides
whether the downstream arms are readable at all - a concrete PDB id and chain and
a count of approved drugs that bind it.

Nothing here is new science; it is the wiring of three chains that already
existed in this repo but were never connected to the demo:

  scripts/disease_targets.py  Open Targets: disease -> EFO id -> associated
                              targets with per-datatype scores and drug counts.
                              Imported wholesale (resolve_disease,
                              associated_targets, build_target, gql).
  scripts/target_evidence.py  Europe PMC: verbatim claim sentences, every PMID
                              re-fetched and title-verified. Imported wholesale
                              (search_target, verify, aliases).
  scripts/prep_target.py      RCSB: UniProt accession -> polymer entities, entry
                              resolution, bound heteroatoms. Imported for its
                              HTTP/cache helpers and its ligand-vs-additive
                              rules (NOT_LIGANDS, MIN_LIGAND_MW).

What is written fresh here is the selection logic, which is where the two
non-obvious criteria live:

1. APPROVED BINDERS. A target is only usable in this pipeline if approved drugs
   already bind it, because the matching arms need a positive control; without
   one a ranking of 2,153 approved drugs has no known-correct answer and cannot
   be scored. This is the entire difference between thrombin (4 drugs whose
   mechanism-of-action target is F2, plus bivalirudin and lepirudin as approved
   peptides) and CDK2 (0 MoA drugs). Counted here from three local corpora
   rather than from Open Targets' own drug counts, because the arms are matched
   against these exact files:
     data/approved_drugs.csv       filtered to `approved == 1` - the file has
                                   4,099 rows but only 2,153 approved ones
     data/target_annotations.csv   ChEMBL pchembl-backed binding, joined to the
                                   approved subset on struct_id
     data/approved_biologics.csv   262 approved peptides/biologics, for the
                                   peptide arm's positive control

2. CIRCULARITY. PROJECT_GOAL.md 1.3: Open Targets association evidence already
   contains ChEMBL known-drug data, so a target chosen because drugs exist for
   it and then "discovered" to have those drugs is a database lookup in a lab
   coat. Every target here reports its full `score_by_datatype` and a
   `known_drug_share_of_datatype_sum`, and is flagged when the known-drug
   datatype is its single largest line of evidence. The flag is a warning
   printed next to the result, not a filter - criterion 1 and this one pull in
   opposite directions and the honest thing is to show both.

STRUCTURE SELECTION is the other thing done fresh. Downstream needs a PDB id and
a chain, and `demo/pipeline.py` strips heteroatoms before design. A co-crystallised
*peptide* is a polymer chain, not a heteroatom, so `remove_heterogens` leaves it
sitting in the site - the mistake documented in `demo/pipeline.py`'s thrombin
comment (4UD9/4UE7/5AFY were rejected by hand for exactly this). The rule here
rejects such entries programmatically: any polymer entity in the entry that does
not map to the target's UniProt accession and is <= 50 aa disqualifies the entry.
Targets with no entry surviving are reported as such, not given an invented id.

Usage:
    ./env-kit/bin/python demo/autoresearch.py --disease "colorectal cancer"
    ./env-kit/bin/python demo/autoresearch.py --disease "venous thromboembolism"
    ./env-kit/bin/python demo/autoresearch.py --disease MONDO_0005575 --top 15
    ./env-kit/bin/python demo/autoresearch.py --disease "..." --emit-config F2
    ./env-kit/bin/python -m demo.test_autoresearch        # self-test

Free: Open Targets, RCSB, Europe PMC and UniProt are open APIs, no credentials,
no Rowan credits. Every response is cached under data/raw/ (gitignored) so a
second run is offline and instant.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sys
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import disease_targets as DT  # noqa: E402
import prep_target as PT  # noqa: E402
import target_evidence as TE  # noqa: E402

OUT_DIR = ROOT / "results" / "demo" / "autoresearch"
RAW = ROOT / "data" / "raw"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb"

APPROVED_DRUGS = ROOT / "data" / "approved_drugs.csv"
TARGET_ANNOTATIONS = ROOT / "data" / "target_annotations.csv"
APPROVED_BIOLOGICS = ROOT / "data" / "approved_biologics.csv"

# A polymer entity this short, from a different protein than the target, is a
# co-crystallised peptide sitting in the site. `remove_heterogens` does not
# remove it (it is ATOM records in a polymer chain), so an entry carrying one is
# unusable for design even though it may be the highest-resolution entry there is.
PEPTIDE_POLYMER_MAX_AA = 50

# How many polymer entities to pull detail for per target, best resolution first.
# Well-studied targets have hundreds (P00734: 900, P00533: >1000); the top slice
# by resolution is both what the selection rule wants and what keeps the run to a
# few requests per target. `n_entities_examined` vs `n_entities_total` is written
# into the output so the cap is visible rather than implied.
ENTITY_CAP = 60
ENTITY_BATCH = 25

# ACT_SITE residues that are the catalytic nucleophile - the atom a ligand
# actually sits on. Serine-protease triads list three active sites (His, Asp,
# Ser); a sphere around all three reaches behind the serine into the protein
# core, which is why demo/pipeline.py anchors thrombin on Ser195 alone. Selecting
# by residue identity reproduces that choice without hardcoding it, and returns
# nothing for a kinase whose ACT_SITE is the Asp proton acceptor.
NUCLEOPHILE_AA = {"S", "C", "T"}

ENTITY_QUERY = """{polymer_entities(entity_ids:%s){
 rcsb_id
 entity_poly{rcsb_sample_sequence_length type}
 rcsb_polymer_entity{pdbx_description}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids}
 rcsb_polymer_entity_align{reference_database_accession aligned_regions{ref_beg_seq_id length}}
 entry{rcsb_id
  rcsb_entry_info{resolution_combined experimental_method}
  polymer_entities{rcsb_id
   entity_poly{rcsb_sample_sequence_length type}
   rcsb_polymer_entity{pdbx_description}
   rcsb_polymer_entity_container_identifiers{auth_asym_ids}
   rcsb_polymer_entity_align{reference_database_accession}}
  nonpolymer_entities{
   rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}
   nonpolymer_comp{chem_comp{id name formula_weight}}}}}}"""

TARGET_EXTRAS_QUERY = """{target(ensemblId:"%s"){
 approvedSymbol
 targetClass{label level}
 tractability{label modality value}
 safetyLiabilities{event eventId datasource}}}"""

MODALITY_NAMES = {"SM": "small_molecule", "AB": "antibody", "PR": "protac", "OC": "other"}


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def split_symbols(cell):
    """`targets` / `moa_targets` are ';'-joined HGNC symbols, or empty."""
    return [s.strip() for s in (cell or "").split(";") if s.strip()]


# --------------------------------------------------------------------------
# local corpora - the arms' own positive controls
# --------------------------------------------------------------------------
def load_corpora():
    """The three local files, indexed by gene symbol and by UniProt accession.

    Only `approved == 1` rows of approved_drugs.csv are used. The file name
    promises approval but the column is what decides it: 4,099 rows, 2,153 with
    approved == 1. Both counts are returned so the filter is auditable.
    """
    drugs = list(csv.DictReader(APPROVED_DRUGS.open(newline="")))
    approved = [r for r in drugs if r["approved"] == "1"]
    approved_ids = {r["struct_id"] for r in approved}

    by_symbol = defaultdict(set)  # any annotated target, incl. off-targets
    by_symbol_moa = defaultdict(set)  # mechanism-of-action target only
    for r in approved:
        for s in split_symbols(r["targets"]):
            by_symbol[s].add(r["name"])
        for s in split_symbols(r["moa_targets"]):
            by_symbol_moa[s].add(r["name"])

    # ChEMBL bioactivity, restricted to approved drugs and to rows the corpus
    # build already judged to pass its potency threshold.
    by_uniprot_act = defaultdict(set)
    n_ann_rows = n_ann_approved = 0
    with TARGET_ANNOTATIONS.open(newline="") as fh:
        for r in csv.DictReader(fh):
            n_ann_rows += 1
            if r["struct_id"] not in approved_ids or r["passes_threshold"] != "1":
                continue
            n_ann_approved += 1
            if r["uniprot"]:
                by_uniprot_act[r["uniprot"]].add(r["name"])

    bio = list(csv.DictReader(APPROVED_BIOLOGICS.open(newline="")))
    bio_by_symbol = defaultdict(list)
    bio_by_uniprot = defaultdict(list)
    for r in bio:
        rec = {
            "name": r["name"],
            "modality": r["modality"],
            "length_aa": r["length_aa"],
            "action": r["target_action_type"],
        }
        if r["target_gene_symbol"]:
            bio_by_symbol[r["target_gene_symbol"]].append(rec)
        if r["target_uniprot"]:
            bio_by_uniprot[r["target_uniprot"]].append(rec)

    return {
        "n_drug_rows": len(drugs),
        "n_approved": len(approved),
        "n_annotation_rows": n_ann_rows,
        "n_annotation_rows_approved_and_potent": n_ann_approved,
        "n_biologics": len(bio),
        "by_symbol": by_symbol,
        "by_symbol_moa": by_symbol_moa,
        "by_uniprot_act": by_uniprot_act,
        "bio_by_symbol": bio_by_symbol,
        "bio_by_uniprot": bio_by_uniprot,
    }


def approved_binders(symbol, uniprot, synonyms, corpora):
    """Approved drugs that bind this target, from the corpora the arms use.

    Three counts, because they mean different things and collapsing them would
    hide the CDK2 case: `n_moa` is the number whose *mechanism* target is this
    gene (a real positive control), `n_annotated` includes off-target
    annotations, `n_bioactivity` is pchembl-backed binding at the UniProt level.
    Symbol lookups try the approved symbol and every Open Targets synonym, since
    the corpus is keyed on whatever spelling DrugCentral/ChEMBL carried.
    """
    names = [symbol] + [s for s in synonyms if s != symbol]
    ann, moa = set(), set()
    matched = []
    for n in names:
        if corpora["by_symbol"].get(n) or corpora["by_symbol_moa"].get(n):
            matched.append(n)
        ann |= corpora["by_symbol"].get(n, set())
        moa |= corpora["by_symbol_moa"].get(n, set())

    act = corpora["by_uniprot_act"].get(uniprot, set()) if uniprot else set()

    bio = list(corpora["bio_by_uniprot"].get(uniprot, [])) if uniprot else []
    seen = {b["name"] for b in bio}
    for n in names:
        for b in corpora["bio_by_symbol"].get(n, []):
            if b["name"] not in seen:
                seen.add(b["name"])
                bio.append(b)
    peptides = [b for b in bio if b["modality"] == "peptide"]

    return {
        "n_approved_binders_moa": len(moa),
        "n_approved_binders_annotated": len(ann),
        "n_approved_binders_bioactivity": len(act),
        "n_approved_biologic_binders": len(bio),
        "n_approved_peptide_binders": len(peptides),
        "known_approved_drugs": sorted(moa),
        "approved_binders_offtarget_only": sorted(ann - moa)[:12],
        "approved_biologic_binders": sorted(
            ({"name": b["name"], "modality": b["modality"], "action": b["action"]}
             for b in bio),
            key=lambda b: b["name"],
        ),
        "corpus_symbols_matched": matched,
    }


# --------------------------------------------------------------------------
# UniProt: mature chains and catalytic anchors
# --------------------------------------------------------------------------
def uniprot_sites(acc):
    """Sequence length, mature CHAIN features, and catalytic nucleophile anchors.

    The anchors are the residue numbers `demo/pipeline.py` calls `site_anchors`.
    An ACT_SITE is kept only if the residue at that position is S/C/T, i.e. the
    nucleophile rather than the general base or acid that sits behind it.
    """
    fields = "accession,sequence,ft_act_site,ft_chain"
    url = f"{UNIPROT_REST}/{acc}.json?{urllib.parse.urlencode({'fields': fields})}"
    try:
        blob = PT.cached(RAW / "uniprot" / f"sites_{acc}.json", lambda: PT.http(url))
        meta = json.loads(blob)
    except Exception as exc:  # noqa: BLE001 - a missing entry is a gap, not a crash
        return {"uniprot_lookup_error": f"{type(exc).__name__}: {exc}"}

    seq = meta["sequence"]["value"]
    act, anchors, chains = [], [], []
    for f in meta.get("features", []):
        loc = f["location"]
        if f["type"] == "Active site":
            pos = loc["start"]["value"]
            aa = seq[pos - 1] if 0 < pos <= len(seq) else "?"
            act.append({"position": pos, "residue": aa, "description": f.get("description")})
            if aa in NUCLEOPHILE_AA:
                anchors.append(pos)
        elif f["type"] == "Chain":
            chains.append(
                {
                    "description": f.get("description"),
                    "start": loc["start"]["value"],
                    "end": loc["end"]["value"],
                }
            )
    return {
        "sequence_length": meta["sequence"]["length"],
        "active_sites": act,
        "site_anchors": anchors,
        "site_anchor_rule": (
            "UniProt ACT_SITE positions whose residue is S/C/T (the catalytic "
            "nucleophile). His/Asp members of a charge-relay triad are excluded: "
            "they sit behind the nucleophile and a sphere over all three reaches "
            "into the core rather than the ligand site."
        ),
        "mature_chains": chains,
    }


# --------------------------------------------------------------------------
# RCSB: a usable PDB id and chain, or an honest "none"
# --------------------------------------------------------------------------
def _term(attribute, operator, value):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def _entity_search(acc, extra_nodes, cap, tag):
    nodes = [
        _term(
            "rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_accession",
            "in",
            [acc],
        ),
        _term(
            "rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_name",
            "exact_match",
            "UniProt",
        ),
        _term("rcsb_entry_info.experimental_method", "exact_match", "X-ray"),
    ] + list(extra_nodes)
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": cap},
            "results_verbosity": "compact",
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
            ],
        },
    }

    def go():
        # A search that matches nothing returns HTTP 204 with an empty body,
        # which PT.http treats as a failure. Try once so an empty result set does
        # not burn three retries, and fall back to the retrying call for anything
        # that is a real transport error rather than "no hits".
        try:
            return PT.http(PT.RCSB_SEARCH, q, retries=1)
        except Exception as exc:  # noqa: BLE001
            if "empty body" in str(exc):
                return json.dumps({"total_count": 0, "result_set": []}).encode()
            return PT.http(PT.RCSB_SEARCH, q)

    p = json.loads(PT.cached(RAW / "rcsb" / f"autoresearch_{tag}_{acc}_{cap}.json", go))
    return p.get("result_set", []), p.get("total_count", 0)


def search_entities_by_resolution(acc, cap=ENTITY_CAP):
    """X-ray polymer entities mapping to `acc`, best resolution first.

    Two searches, unioned in order, because resolution alone is a bad candidate
    generator for a well-studied target: P00734's 60 highest-resolution entities
    are almost all hirudin co-complexes, and every one of them is disqualified by
    the peptide rule, leaving the survivor apo. So the first search is biased
    toward entries that can be holo and cannot be a large complex -
    `nonpolymer_entity_count > 0` and at most two protein entities (thrombin's
    heavy and light chain are both P00734, so a bar of one would be wrong) - and
    the unrestricted search follows as the fallback for targets with no ligand-
    bound structure at all. The bias only orders the candidates; the selection
    rule in choose_structure() is what decides.
    """
    holo_nodes = [
        _term("rcsb_entry_info.nonpolymer_entity_count", "greater", 0),
        _term("rcsb_entry_info.polymer_entity_count_protein", "less_or_equal", 2),
    ]
    ids, n_holo = _entity_search(acc, holo_nodes, cap, "holo")
    plain, n_all = _entity_search(acc, [], cap, "any")
    seen = set(ids)
    ids = list(ids) + [i for i in plain if not (i in seen or seen.add(i))]
    return ids, n_all, n_holo


def fetch_entity_details(acc, entity_ids):
    def go():
        out = []
        for i in range(0, len(entity_ids), ENTITY_BATCH):
            q = ENTITY_QUERY % json.dumps(entity_ids[i : i + ENTITY_BATCH])
            out += json.loads(PT.http(PT.RCSB_GRAPHQL, {"query": q}))["data"][
                "polymer_entities"
            ]
        return json.dumps(out).encode()

    n = len(entity_ids)
    return json.loads(PT.cached(RAW / "rcsb" / f"autoresearch_ents_{acc}_{n}.json", go))


def _accessions(entity):
    return {
        a.get("reference_database_accession")
        for a in (entity.get("rcsb_polymer_entity_align") or [])
    }


def entry_records(entities, acc):
    """One record per PDB entry, with everything the selection rule needs."""
    entries = {}
    for e in entities:
        entry = e["entry"]
        pdb_id = entry["rcsb_id"]
        if pdb_id in entries:
            continue
        info = entry["rcsb_entry_info"]
        res = (info.get("resolution_combined") or [None])[0]

        target_chains, foreign_peptides, foreign_polymers = [], [], []
        for pe in entry.get("polymer_entities") or []:
            length = (pe.get("entity_poly") or {}).get("rcsb_sample_sequence_length")
            ptype = (pe.get("entity_poly") or {}).get("type")
            chains = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get(
                "auth_asym_ids"
            ) or []
            desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")
            rec = {
                "entity_id": pe["rcsb_id"],
                "length_aa": length,
                "polymer_type": ptype,
                "auth_chain_ids": chains,
                "description": desc,
                "accessions": sorted(a for a in _accessions(pe) if a),
            }
            if acc in _accessions(pe):
                target_chains.append(rec)
            elif ptype == "polypeptide(L)" and (length or 0) <= PEPTIDE_POLYMER_MAX_AA:
                foreign_peptides.append(rec)
            else:
                foreign_polymers.append(rec)

        ligands = []
        for n in entry.get("nonpolymer_entities") or []:
            c = n["nonpolymer_comp"]["chem_comp"]
            mw = c.get("formula_weight") or 0
            ligands.append(
                {
                    "comp_id": c["id"],
                    "name": c["name"],
                    "formula_weight": mw,
                    "auth_asym_ids": (
                        n["rcsb_nonpolymer_entity_container_identifiers"]["auth_asym_ids"]
                    ),
                    # Same rule scripts/prep_target.py uses, so a ligand called
                    # drug-like here is one that pipeline would also strip as a
                    # real ligand rather than a buffer component.
                    "drug_like": c["id"] not in PT.NOT_LIGANDS and mw >= PT.MIN_LIGAND_MW,
                }
            )
        drug_like = [l["comp_id"] for l in ligands if l["drug_like"]]

        # The chain to hand downstream: the longest entity that maps to the
        # target accession (thrombin's 259-aa heavy chain H, not its 36-aa
        # light chain L), first auth chain id of that entity.
        primary = max(target_chains, key=lambda r: r["length_aa"] or 0, default=None)
        entries[pdb_id] = {
            "pdb_id": pdb_id,
            "method": info.get("experimental_method"),
            "resolution_a": res,
            "chain": (primary or {}).get("auth_chain_ids", [None])[0],
            "chain_entity_length_aa": (primary or {}).get("length_aa"),
            "target_entities": target_chains,
            "foreign_peptide_chains": foreign_peptides,
            "foreign_polymer_chains": foreign_polymers,
            "drug_like_ligands": drug_like,
            "all_het_comp_ids": [l["comp_id"] for l in ligands],
            "apo_or_holo": "holo" if drug_like else "apo",
        }
    return list(entries.values())


def choose_structure(entries):
    """Pick one entry by a stated rule, and keep the rejects with their reason.

    Hard reject: a polymer entity from another protein, <= 50 aa. That is a
    co-crystallised peptide in the site; it is ATOM records, so
    `remove_heterogens` leaves it there and any binder designed against the
    entry is designed against an occupied site. This repo already made that
    mistake once and 4UD9 / 4UE7 / 5AFY were rejected by hand for it.

    Soft penalty: a larger foreign polymer (a Fab, a partner protein, a nucleic
    acid). Also occupying, but the pocket of interest may be elsewhere, so it is
    deprioritised and flagged rather than dropped.

    Order among survivors: holo before apo, no foreign polymer before one, then
    best resolution, then PDB id for determinism.
    """
    usable, rejected = [], []
    for e in entries:
        if not e["chain"]:
            rejected.append({**e, "reject_reason": "no_auth_chain_for_target_accession"})
        elif e["foreign_peptide_chains"]:
            names = [
                f"{p['entity_id']} {p['description']} ({p['length_aa']} aa)"
                for p in e["foreign_peptide_chains"]
            ]
            rejected.append(
                {
                    **e,
                    "reject_reason": "co_crystallised_peptide_polymer_chain",
                    "reject_detail": names,
                }
            )
        elif e["resolution_a"] is None:
            rejected.append({**e, "reject_reason": "no_resolution_reported"})
        else:
            usable.append(e)

    usable.sort(
        key=lambda e: (
            e["apo_or_holo"] != "holo",
            bool(e["foreign_polymer_chains"]),
            e["resolution_a"],
            e["pdb_id"],
        )
    )
    return (usable[0] if usable else None), usable, rejected


def structure_for(acc):
    """UniProt accession -> chosen entry + the audit trail behind the choice."""
    if not acc:
        return {"structure_status": "no_uniprot_accession", "structure": None}
    ids, total, n_holo_candidates = search_entities_by_resolution(acc)
    if not ids:
        return {
            "structure_status": "no_xray_entity_maps_to_this_accession",
            "n_entities_total": total,
            "n_entities_ligand_bound_total": n_holo_candidates,
            "n_entities_examined": 0,
            "structure": None,
        }
    ents = fetch_entity_details(acc, ids)
    entries = entry_records(ents, acc)
    chosen, usable, rejected = choose_structure(entries)
    reasons = Counter(r["reject_reason"] for r in rejected)
    return {
        "structure_status": "ok" if chosen else "no_entry_survived_selection",
        "n_entities_total": total,
        "n_entities_ligand_bound_total": n_holo_candidates,
        "n_entities_examined": len(ids),
        "n_entries_examined": len(entries),
        "n_entries_usable": len(usable),
        "n_entries_rejected": len(rejected),
        "reject_reason_counts": dict(reasons),
        "structure": chosen,
        # The rejects that matter for the demo's headline: the entries that
        # would have been picked on resolution alone, and were not.
        "rejected_better_resolution": [
            {
                "pdb_id": r["pdb_id"],
                "resolution_a": r["resolution_a"],
                "reject_reason": r["reject_reason"],
                "reject_detail": r.get("reject_detail"),
            }
            for r in rejected
            if chosen
            and r["resolution_a"] is not None
            and r["resolution_a"] < chosen["resolution_a"]
        ][:10],
        "runner_up_entries": [
            {
                "pdb_id": e["pdb_id"],
                "chain": e["chain"],
                "resolution_a": e["resolution_a"],
                "apo_or_holo": e["apo_or_holo"],
            }
            for e in usable[1:6]
        ],
    }


# --------------------------------------------------------------------------
# Open Targets extras: tractability, safety, protein class
# --------------------------------------------------------------------------
def target_extras(ensembl_id, refresh=False):
    try:
        d = DT.gql(
            TARGET_EXTRAS_QUERY % ensembl_id, f"autoresearch-extras-{ensembl_id}", refresh
        )["data"]["target"]
    except Exception as exc:  # noqa: BLE001 - record the gap, keep the target
        return {"tractability_error": f"{type(exc).__name__}: {exc}"}
    if d is None:
        return {"tractability_error": "open targets returned no target record"}

    buckets = defaultdict(list)
    for t in d.get("tractability") or []:
        if t["value"]:
            buckets[MODALITY_NAMES.get(t["modality"], t["modality"])].append(t["label"])
    safety = d.get("safetyLiabilities") or []
    return {
        "target_class": sorted(
            {c["label"] for c in (d.get("targetClass") or []) if c["level"] in ("l1", "l2")}
        ),
        "target_class_detail": d.get("targetClass") or [],
        "tractability": {k: sorted(v) for k, v in sorted(buckets.items())},
        "n_safety_liabilities": len(safety),
        "safety_flags": sorted({s["event"] for s in safety if s.get("event")}),
        "safety_sources": sorted({s["datasource"] for s in safety if s.get("datasource")}),
    }


# --------------------------------------------------------------------------
# circularity
# --------------------------------------------------------------------------
def circularity(all_datatype_scores, known_drug_id):
    """How much of this association is 'a drug already hits this target'.

    The share is a DESCRIPTIVE statistic over the datatype scores, not a
    decomposition of the overall score - Open Targets aggregates datatypes by a
    weighted harmonic sum, so these do not add up to `overall_association_score`
    and must not be read as if they did. It answers one question only: among the
    lines of evidence, how large is the known-drug one.
    """
    scores = {k: v for k, v in (all_datatype_scores or {}).items() if v is not None}
    kd = scores.get(known_drug_id) if known_drug_id else None
    total = sum(scores.values())
    if kd is None:
        return {
            "known_drug_datatype_id": known_drug_id,
            "known_drug_score": None,
            "known_drug_share_of_datatype_sum": None,
            "known_drug_is_largest_datatype": False,
            "known_drug_dominated": False,
            "circularity_note": "no known-drug datatype on this target",
        }
    largest = max(scores, key=scores.get)
    share = kd / total if total else None
    dominated = largest == known_drug_id and (share or 0) >= 0.25
    return {
        "known_drug_datatype_id": known_drug_id,
        "known_drug_score": kd,
        "known_drug_share_of_datatype_sum": round(share, 4) if share is not None else None,
        "known_drug_is_largest_datatype": largest == known_drug_id,
        "largest_datatype": largest,
        "known_drug_dominated": dominated,
        "circularity_note": (
            "known-drug evidence is this target's single largest datatype "
            f"({kd:.3f}, {share:.0%} of the datatype-score sum); a 'discovery' of "
            "its approved drugs is partly a lookup of the evidence that ranked it"
            if dominated
            else f"largest datatype is {largest} ({scores[largest]:.3f}); "
            f"known-drug contributes {kd:.3f}"
        ),
    }


# --------------------------------------------------------------------------
# literature
# --------------------------------------------------------------------------
def literature_for(target, disease_name, per_target, refresh=False):
    """Europe PMC via scripts/target_evidence.py, every PMID title-verified."""
    if per_target <= 0:
        return {"literature": [], "literature_status": "skipped (--lit 0)"}
    try:
        query, hits, papers = TE.search_target(target, disease_name, per_target, refresh)
    except Exception as exc:  # noqa: BLE001
        return {"literature": [], "literature_status": f"error: {type(exc).__name__}: {exc}"}

    kept, dropped = [], []
    for p in papers:
        ok, verdict, back = TE.verify(p["pmid"], p["title"], refresh)
        if ok:
            kept.append({**p, "title_verified": True, "verification": verdict})
        else:
            dropped.append({"pmid": p["pmid"], "reason": verdict, "title_returned": back})
    return {
        "literature_query": query,
        "literature_hit_count": hits,
        "literature": kept,
        "literature_dropped_failed_verification": dropped,
        "literature_status": "ok" if kept else "no verified paper",
    }


# --------------------------------------------------------------------------
# site type / modality / pipeline config
# --------------------------------------------------------------------------
def site_and_modalities(struct, binders, extras):
    """Derived from what was measured, never asserted.

    site_type_hint:
      pocket        the chosen structure holds a drug-like heteroatom ligand -
                    direct evidence of a small-molecule site in this construct
      ppi_interface no such ligand, but the entry carries a foreign polymer
                    partner
      unknown       neither, or no structure
    preferred_modalities: a modality is listed only where an approved drug of
    that kind binds the target, or the structure shows the corresponding site.
    """
    mods = []
    if binders["n_approved_binders_moa"] > 0:
        mods.append("small_molecule")
    if binders["n_approved_peptide_binders"] > 0:
        mods.append("peptide")
    if binders["n_approved_biologic_binders"] > binders["n_approved_peptide_binders"]:
        mods.append("biologic")

    s = struct.get("structure")
    if s and s["drug_like_ligands"]:
        hint = "pocket"
        why = f"chosen structure {s['pdb_id']} binds {'/'.join(s['drug_like_ligands'])} as HETATM"
        if "small_molecule" not in mods:
            mods.append("small_molecule")
    elif s and s["foreign_polymer_chains"]:
        hint = "ppi_interface"
        why = f"{s['pdb_id']} has a partner polymer chain and no drug-like heteroatom"
    elif s:
        hint = "unknown"
        why = f"{s['pdb_id']} is apo and has no partner chain"
    else:
        hint = "unknown"
        why = "no structure selected"

    # Tractability is Open Targets' opinion and is reported, not used to gate -
    # it is recorded next to the modality list so the two can be compared.
    return {
        "site_type_hint": hint,
        "site_type_evidence": why,
        "preferred_modalities": mods or ["unknown"],
        "preferred_modalities_rule": (
            "a modality is listed only where an approved drug of that kind binds "
            "this target in the local corpora, or (small_molecule) the chosen "
            "structure carries a drug-like HETATM ligand"
        ),
        "opentargets_tractability_modalities": sorted(extras.get("tractability", {})),
    }


def pipeline_config(h):
    """A `demo/pipeline.py` TARGETS entry, in exactly that dict's shape."""
    s = h["structure"]["structure"]
    if s is None:
        return None
    peptide = "peptide" in h["preferred_modalities"]
    cfg = {
        "symbol": h["approved_symbol"],
        "uniprot": h["uniprot"],
        "pdb_id": s["pdb_id"],
        "chain": s["chain"],
        "out": h["approved_symbol"].lower(),
        # An approved peptide already binds this target, so the peptide arm has
        # a positive control and an 8..16mer is the length class it will be
        # matched against. Otherwise the protein protocol, as on CDK2.
        "protocol": "peptide-anything" if peptide else "protein-anything",
        "binder_length": "8..16" if peptide else "60..90",
    }
    anchors = h.get("uniprot_sites", {}).get("site_anchors") or []
    if anchors:
        cfg["site_anchors"] = anchors
    return cfg


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
def build_hypothesis(row, corpora, disease_name, lit_per_target, refresh, structures=True):
    base = DT.build_target(row["rank"], row["row"], refresh)
    sym, uni = base["approved_symbol"], base["uniprot"]

    binders = approved_binders(sym, uni, base.get("symbol_synonyms") or [], corpora)
    extras = target_extras(base["ensembl_id"], refresh)
    circ = circularity(base["all_datatype_scores"], base["known_drug_datatype_id"])
    sites = uniprot_sites(uni) if uni else {"uniprot_lookup_error": "no accession"}
    struct = (
        structure_for(uni) if structures else {"structure_status": "skipped", "structure": None}
    )
    lit = literature_for(base, disease_name, lit_per_target, refresh)

    h = {
        **base,
        "association_score": base["overall_association_score"],
        "score_by_datatype": base["all_datatype_scores"],
        **binders,
        **extras,
        "circularity": circ,
        "uniprot_sites": sites,
        "structure": struct,
        **lit,
    }
    h.update(site_and_modalities(struct, binders, extras))

    # The two gates, applied and recorded separately so a failure says which.
    has_control = binders["n_approved_binders_moa"] >= 1
    has_struct = struct.get("structure") is not None
    reasons = []
    if not has_control:
        reasons.append(
            f"no approved drug has {sym} as its mechanism-of-action target "
            f"({binders['n_approved_binders_annotated']} off-target annotations, "
            f"{binders['n_approved_binders_bioactivity']} pchembl-backed) - the "
            "matching arms would have no positive control"
        )
    if not has_struct:
        reasons.append(f"no usable structure ({struct['structure_status']})")
    h["gate_has_approved_positive_control"] = has_control
    h["gate_has_usable_structure"] = has_struct
    h["usable_for_pipeline"] = has_control and has_struct
    h["unusable_reasons"] = reasons
    h["pipeline_config"] = pipeline_config(h) if h["usable_for_pipeline"] else None
    return h


def run(disease, top, lit_per_target, refresh, structures=True, write=True, verbose=True):
    corpora = load_corpora()
    meta = DT.api_meta(refresh)

    if re.fullmatch(r"(EFO|MONDO|HP|Orphanet|MP|OTAR)_\d+", disease):
        efo_id, disease_name, rule, hits, n_hits = disease, None, "id_given_verbatim", [], None
    else:
        hit, rule, hits, n_hits = DT.resolve_disease(disease, refresh)
        efo_id, disease_name = hit["id"], hit["name"]

    assoc = DT.associated_targets(efo_id, top, refresh)
    disease_name = disease_name or assoc["name"]
    rows = assoc["associatedTargets"]["rows"]
    if verbose:
        print(
            f"{disease!r} -> {efo_id} {disease_name!r} ({rule}); "
            f"{assoc['associatedTargets']['count']} associated targets, taking top {len(rows)}",
            file=sys.stderr,
        )

    hyps = []
    for i, row in enumerate(rows, 1):
        h = build_hypothesis(
            {"rank": i, "row": row}, corpora, disease_name, lit_per_target, refresh, structures
        )
        hyps.append(h)
        s = h["structure"].get("structure")
        if verbose:
            print(
                f"  {i:2d} {h['approved_symbol']:<9} assoc={h['association_score']:.4f} "
                f"moa={h['n_approved_binders_moa']:<3} bio={h['n_approved_biologic_binders']:<2} "
                f"pdb={(s or {}).get('pdb_id', '-'):<5} chain={(s or {}).get('chain') or '-'} "
                f"{'USABLE' if h['usable_for_pipeline'] else 'unusable'}",
                file=sys.stderr,
            )

    usable = [h for h in hyps if h["usable_for_pipeline"]]
    # Stated ranking: gates first, then Open Targets' own association order. No
    # invented composite weighting - nothing here has been validated well enough
    # to justify one.
    ranked = sorted(usable, key=lambda h: -h["association_score"])
    for k, h in enumerate(ranked, 1):
        h["pipeline_rank"] = k

    dominated = [h["approved_symbol"] for h in ranked if h["circularity"]["known_drug_dominated"]]
    no_control = [h["approved_symbol"] for h in hyps if not h["gate_has_approved_positive_control"]]
    no_struct = [h["approved_symbol"] for h in hyps if not h["gate_has_usable_structure"]]
    peptide_arm = [h["approved_symbol"] for h in ranked if h["n_approved_peptide_binders"] > 0]

    out = {
        "meta": {
            "script": "demo/autoresearch.py",
            "generated_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "disease_query": disease,
            "efo_id": efo_id,
            "disease_name": disease_name,
            "disease_resolution_rule": rule,
            "disease_search_hits_considered": hits,
            "n_disease_hits": n_hits,
            "n_associated_targets_total": assoc["associatedTargets"]["count"],
            "n_targets_examined": len(hyps),
            "n_usable": len(ranked),
            "n_failing_positive_control_gate": len(no_control),
            "n_failing_structure_gate": len(no_struct),
            "n_usable_with_peptide_control": len(peptide_arm),
            "targets_without_positive_control": no_control,
            "targets_without_usable_structure": no_struct,
            "known_drug_dominated_among_usable": dominated,
            "corpora": {
                "approved_drugs_csv_rows": corpora["n_drug_rows"],
                "approved_drugs_csv_approved_rows": corpora["n_approved"],
                "target_annotations_rows": corpora["n_annotation_rows"],
                "target_annotations_rows_used": corpora[
                    "n_annotation_rows_approved_and_potent"
                ],
                "approved_biologics_rows": corpora["n_biologics"],
            },
            "gates": {
                "positive_control": "n_approved_binders_moa >= 1",
                "structure": (
                    "an X-ray entry maps to the accession, has an auth chain for it, "
                    "reports a resolution, and carries no foreign polypeptide of "
                    f"<= {PEPTIDE_POLYMER_MAX_AA} aa"
                ),
                "ranking": "gates first, then Open Targets association score descending",
            },
            "circularity_caveat": (
                "Open Targets association evidence includes ChEMBL known-drug data "
                "(PROJECT_GOAL.md 1.3). The positive-control gate selects FOR targets "
                "with approved drugs, so it selects for exactly the evidence that is "
                "circular. known_drug_dominated_among_usable is the size of that "
                "overlap; it is reported, not filtered."
            ),
            "opentargets": meta,
            "reused_from": [
                "scripts/disease_targets.py: resolve_disease, associated_targets, "
                "build_target, gql, api_meta",
                "scripts/target_evidence.py: search_target, verify, aliases, epmc",
                "scripts/prep_target.py: http, cached, NOT_LIGANDS, MIN_LIGAND_MW, "
                "RCSB endpoints",
            ],
        },
        "targets": hyps,
    }

    path = OUT_DIR / f"{slugify(disease)}.json"
    if write:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(out, indent=1) + "\n")
    return out, path


def print_table(out):
    hyps = out["targets"]
    m = out["meta"]
    print(f"\n{m['disease_name']} ({m['efo_id']}) - top {m['n_targets_examined']} "
          f"of {m['n_associated_targets_total']} associated targets\n")
    hdr = (f"{'#':>2} {'symbol':<9}{'uniprot':<9}{'assoc':>7} {'kd%':>5} "
           f"{'moa':>4}{'bio':>4}{'pep':>4}  {'pdb':<5}{'ch':<3}{'res':>5} {'site':<14}verdict")
    print(hdr)
    print("-" * len(hdr))
    for h in hyps:
        s = h["structure"].get("structure") or {}
        c = h["circularity"]
        kd = c["known_drug_share_of_datatype_sum"]
        verdict = "USABLE" if h["usable_for_pipeline"] else (
            "no control" if not h["gate_has_approved_positive_control"] else "no structure")
        if h["usable_for_pipeline"] and c["known_drug_dominated"]:
            verdict += " (known-drug dominated)"
        print(f"{h['rank']:>2} {h['approved_symbol']:<9}{h['uniprot'] or '-':<9}"
              f"{h['association_score']:>7.4f} {f'{kd:.0%}' if kd is not None else '-':>5} "
              f"{h['n_approved_binders_moa']:>4}{h['n_approved_biologic_binders']:>4}"
              f"{h['n_approved_peptide_binders']:>4}  {s.get('pdb_id', '-'):<5}"
              f"{s.get('chain') or '-':<3}{s.get('resolution_a') or 0:>5.2f} "
              f"{h['site_type_hint']:<14}{verdict}")
    print(f"\nn = {m['n_targets_examined']} examined, {m['n_usable']} usable "
          f"({m['n_failing_positive_control_gate']} lack an approved positive control, "
          f"{m['n_failing_structure_gate']} lack a usable structure)")
    if m["known_drug_dominated_among_usable"]:
        print(f"circularity: {len(m['known_drug_dominated_among_usable'])}/{m['n_usable']} "
              f"usable targets are known-drug dominated: "
              f"{', '.join(m['known_drug_dominated_among_usable'])}")
    ranked = sorted([h for h in hyps if h["usable_for_pipeline"]], key=lambda h: h["pipeline_rank"])
    if ranked:
        top = ranked[0]
        print(f"\npipeline would pick: {top['approved_symbol']} "
              f"(association rank {top['rank']}) ->")
        print(json.dumps(top["pipeline_config"], indent=1))
    else:
        print("\npipeline would pick: nothing - no target passes both gates")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--disease", default="colorectal cancer",
                    help="free-text name, or an EFO/MONDO id")
    ap.add_argument("--top", type=int, default=25, help="associated targets to examine")
    ap.add_argument("--lit", type=int, default=2, help="verified papers per target (0 to skip)")
    ap.add_argument("--no-structures", action="store_true", help="skip RCSB (debug)")
    ap.add_argument("--refresh", action="store_true", help="ignore the data/raw cache")
    ap.add_argument("--emit-config", metavar="SYMBOL",
                    help="print only the demo/pipeline.py TARGETS entry for this symbol")
    args = ap.parse_args()

    out, path = run(args.disease, args.top, args.lit, args.refresh,
                    structures=not args.no_structures)

    if args.emit_config:
        want = args.emit_config.upper()
        hit = next((h for h in out["targets"] if h["approved_symbol"].upper() == want), None)
        if hit is None:
            raise SystemExit(f"{want} is not in the top {args.top} for {args.disease!r}")
        if hit["pipeline_config"] is None:
            raise SystemExit(
                f"{want} is not usable: " + "; ".join(hit["unusable_reasons"])
            )
        print(json.dumps({hit["pipeline_config"]["out"]: hit["pipeline_config"]}, indent=4))
        return

    print_table(out)
    print(f"\nwrote {path.relative_to(ROOT)}", file=sys.stderr)


if __name__ == "__main__":
    main()
