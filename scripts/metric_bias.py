"""E3 - is Jaccard the wrong metric for cross-modality binding-site comparison?

THE CONCERN IS ARITHMETIC, NOT TASTE. A 12-residue small-molecule footprint
lying entirely INSIDE a 40-residue epitope scores Jaccard 12/40 = 0.30, the same
number a genuinely partial 26-vs-26-residue overlap scores. If the two are
indistinguishable, then the repo's 0.6355 small-molecule Jaccard
(results/m2_gate_KDR.json) and its 0.0255-0.1836 biologic Jaccards
(results/epitope_gate.json) are not on one scale and must never be pooled.

This script answers that in three parts:

  1. ANALYTICALLY. For nested sets (|A|=a inside |B|=b) it tabulates Jaccard,
     coverage |A&B|/|A|, the overlap coefficient |A&B|/min(|A|,|B|) and F1
     across size ratios, and enumerates, for a grid of Jaccard values, the
     (perfect containment) and (partial overlap) configurations that produce
     exactly that Jaccard. Where each metric breaks is read off the table.

  2. ON REAL STRUCTURES, which is what makes this evidence rather than an
     argument. A cross-modality gold set is built from RCSB: targets with BOTH a
     drug-like small-molecule co-crystal AND a peptide/protein/antibody
     co-crystal, both footprints extracted in the SAME target-side coordinates
     by scripts/interfaces.py (ligand_contacts for the HET ligand,
     chain_contacts for the polymer binder; one 4.5 A heavy-atom cutoff, one
     numbering convention, by construction). Residue ids from two different
     entries are only comparable if both index the same UniProt sequence, so
     that is CHECKED per entry (interfaces.check_author_numbering) and a pair
     that fails is reported in not_evaluated rather than scored.

  3. AS A PIPELINE DECISION. The gate numbers are prediction-vs-truth, so the
     size-penalty question there is: what is the LARGEST Jaccard a perfect
     prediction of the size P2Rank actually emits could have scored against the
     truth set it was scored on? That ceiling is computed per class from the
     n_pred / interface-size columns already in results/epitope_gate.json and
     results/m2_gate_KDR.json, which says how much of the small-molecule/biologic
     gap is arithmetic and how much is real.

Writes results/metric_bias.json. Spends no Rowan credits; network use is RCSB
search/GraphQL/coordinates only, cached under data/raw/metric_bias/.

Usage:
    ./env/bin/python scripts/metric_bias.py            # analytic + real set
    ./env/bin/python scripts/metric_bias.py --analytic-only
    ./env/bin/python scripts/metric_bias.py --max-pairs 40 --refresh
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

import epitope_gate as E  # noqa: E402  (RCSB http_post/GraphQL helpers, one copy)
import interfaces as I  # noqa: E402
import metrics as M  # noqa: E402

CACHE = ROOT / "data" / "raw" / "metric_bias"
STRUCTURES = CACHE / "structures"
EPITOPE_CACHE = ROOT / "data" / "raw" / "epitope_gate"
PREPARED = EPITOPE_CACHE / "prepared.json"
EPITOPE_RESULT = ROOT / "results" / "epitope_gate.json"
KDR_RESULT = ROOT / "results" / "m2_gate_KDR.json"
OUT = ROOT / "results" / "metric_bias.json"

SEED = 0
CUTOFF = 4.5             # identical to interfaces.py / m2_gate / epitope_gate truth
SM_RESOLUTION_MAX = 2.5  # the small-molecule entry; the polymer side is already <= 3.0
MIN_HEAVY_ATOMS = 14     # "drug-like": above fragment size (repo default is 10)
MIN_CONTACTS = 5         # a footprint smaller than this is a crystal contact
NUMBERING_MIN_MATCH = 0.95   # author numbering must index the UniProt sequence
NUMBERING_MIN_COMPARED = 30
MAX_SM_TRIES = 3         # entries downloaded per target before giving up on it
MAX_LIGAND_COMPONENTS = 8  # distinct drug-like components tried per entry
N_BOOT = 10000

# Verified starter from the task brief: CXCR4 with the small molecule IT1t
# (3ODU) and with the 16-mer cyclic peptide CVX15 (3OE0). It is NOT in the
# epitope-gate inventory (3OE0 is 2.9 A with a fusion partner), so it is seeded
# by hand and goes through exactly the same extraction and numbering checks as
# every automatically discovered pair.
SEED_PAIRS = [{"uniprot": "P61073", "poly_pdb": "3OE0", "klass": "peptide",
               "note": "task-brief starter: CXCR4 IT1t (3ODU) vs CVX15 (3OE0)"}]


# --------------------------------------------------------------------------
# the four metrics, defined once
# --------------------------------------------------------------------------
def metric_set(a_set, b_set):
    """All four overlap metrics for two residue sets.

    A is the SMALL-MOLECULE (or, in the gate framing, the PREDICTED) footprint
    and B the biologic (or TRUTH) footprint. The asymmetry matters: coverage is
    defined against |A|, so it answers "how much of A is inside B", which is the
    question containment poses.
    """
    a, b = set(a_set), set(b_set)
    inter = len(a & b)
    union = len(a | b)
    return {
        "n_a": len(a),
        "n_b": len(b),
        "n_intersection": inter,
        "n_union": union,
        "size_ratio": (min(len(a), len(b)) / max(len(a), len(b))
                       if a and b else float("nan")),
        "jaccard": inter / union if union else float("nan"),
        "coverage_a": inter / len(a) if a else float("nan"),
        "coverage_b": inter / len(b) if b else float("nan"),
        "overlap_coefficient": (inter / min(len(a), len(b))
                                if a and b else float("nan")),
        "f1": 2 * inter / (len(a) + len(b)) if (a or b) else float("nan"),
    }


# --------------------------------------------------------------------------
# 1. analytic
# --------------------------------------------------------------------------
def analytic_nested_table(b=40, sizes=(2, 4, 6, 8, 10, 12, 16, 20, 24, 30, 36, 40)):
    """PERFECT containment: A (size a) entirely inside B (size b).

    Coverage and the overlap coefficient are 1.0 by construction at every ratio
    - they say "this footprint is completely inside that one", which is true.
    Jaccard and F1 fall away with the size ratio while nothing about the
    containment changed.
    """
    rows = []
    for a in sizes:
        A = set(range(a))
        B = set(range(b))
        rows.append({"a": a, "b": b, **metric_set(A, B)})
    return rows


def analytic_ambiguity(targets=(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), b=40):
    """For each Jaccard value, a nested pair and a partial-overlap pair that
    both score it exactly.

    Nested: |A|=a inside |B|=b gives J = a/b, so a = round(J*b), coverage 1.0.
    Partial: two sets of equal size n sharing i residues gives J = i/(2n-i);
    fixing n = b and solving, i = 2nJ/(1+J) - coverage i/n < 1.
    Two pairs with the same Jaccard and different coverage means Jaccard alone
    cannot tell containment from partial overlap.
    """
    rows = []
    for j in targets:
        a = int(round(j * b))
        nested = metric_set(set(range(a)), set(range(b)))
        n = b
        i = int(round(2 * n * j / (1 + j)))
        partial = metric_set(set(range(n)), set(range(n - i, 2 * n - i)))
        rows.append({
            "jaccard_target": j,
            "nested": {"description": f"|A|={a} entirely inside |B|={b}", **nested},
            "partial": {"description": f"|A|=|B|={n} sharing {i} residues", **partial},
            "jaccard_difference": abs(nested["jaccard"] - partial["jaccard"]),
            "coverage_a_difference": abs(nested["coverage_a"] - partial["coverage_a"]),
        })
    return rows


def analytic_breakpoints():
    """Where each metric stops being informative, as failing examples."""
    cases = [
        ("jaccard", "perfect containment of a small site in a large one reads as "
                    "partial overlap", set(range(12)), set(range(40))),
        ("f1", "same failure as Jaccard: F1 = 2J/(1+J) is a strictly increasing "
               "function of J, so it carries no information Jaccard does not",
         set(range(12)), set(range(40))),
        ("overlap_coefficient", "saturates at 1.0 whenever the smaller set is "
                                "contained, however small it is - a 2-residue "
                                "footprint inside a 40-residue epitope is a "
                                "perfect score",
         set(range(2)), set(range(40))),
        ("coverage_a", "one-sided: a prediction that engulfs the truth covers it "
                       "completely and is not penalised for the 60 residues it "
                       "invented", set(range(80)), set(range(20))),
    ]
    out = []
    for metric, breaks, A, B in cases:
        out.append({"metric": metric, "breaks_when": breaks,
                    "example": f"|A|={len(A)}, |B|={len(B)}, "
                               f"|A&B|={len(A & B)}",
                    "values": metric_set(A, B)})
    return out


# --------------------------------------------------------------------------
# 2. the real cross-modality gold set
# --------------------------------------------------------------------------
def _fetch(pdb_id, cache_dir):
    return I.fetch_structure(pdb_id, cache_dir=cache_dir)


def load_polymer_side():
    """Candidate polymer-binder complexes: the epitope-gate prepared set.

    Reused rather than re-derived: those 200 entries already passed the gate's
    X-ray/resolution/entity filters, and their receptor and binder chains were
    chosen from the coordinates (one copy of the complex, antibody H/L chains
    verified to pack together). Their cached coordinates are read, never
    written. One entry per target, chosen by best resolution then pdb id - a
    rule that cannot know anything about the overlap being measured.
    """
    rows = json.loads(PREPARED.read_text())["prepared"]
    by_uniprot = defaultdict(list)
    for r in rows:
        if r.get("receptor_uniprot"):
            by_uniprot[r["receptor_uniprot"]].append(r)
    picked = {}
    for acc, rs in by_uniprot.items():
        picked[acc] = sorted(rs, key=lambda r: (r["resolution_a"], r["pdb_id"]))[0]
    return picked


def sm_entry_candidates(accession, rows=12):
    """Entry ids with a nonpolymer entity and a chain of `accession`, best first."""
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "exptl.method", "operator": "exact_match",
            "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": SM_RESOLUTION_MAX}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.nonpolymer_entity_count",
            "operator": "greater_or_equal", "value": 1}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": accession}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}},
    ]
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "entry",
         "request_options": {
             "paginate": {"start": 0, "rows": rows},
             "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                       "direction": "asc"}],
             "results_content_type": ["experimental"]}}
    r = E.http_post(E.SEARCH_URL, q)
    if not r:
        return [], 0
    return [h["identifier"] for h in r["result_set"]], int(r["total_count"])


def _uniprot_chains(entry, accession):
    """auth chain ids in `entry` (a GraphQL record) whose entity cites `accession`.

    EVERY reference accession on the entity is checked, not just the first.
    A fusion construct - CXCR4 with T4 lysozyme spliced into ICL3 (3ODU), the
    standard GPCR crystallisation trick - is ONE entity citing two accessions,
    and reading only the first hid the task brief's own starter pair.
    """
    out = []
    for ent in entry.get("polymer_entities") or []:
        ids = (ent["rcsb_polymer_entity_container_identifiers"] or {})
        refs = ids.get("reference_sequence_identifiers") or []
        accs = {r.get("database_accession") for r in refs if r
                and r.get("database_name") == "UniProt"}
        if accession in accs:
            out.extend(ids.get("auth_asym_ids") or [])
    return out


def numbering_ok(structure, chain_id, accession):
    """Does author numbering on this chain index the UniProt canonical sequence?

    Cross-entry residue ids mean nothing unless both entries agree with the same
    reference. Verified, not assumed; the numbers are returned either way so a
    rejection is readable.
    """
    seq = I.fetch_uniprot_sequence(accession)
    chk = I.check_author_numbering(structure, chain_id, seq)
    chk["accession"] = accession
    chk["chain_id"] = chain_id
    chk["passed"] = (chk["n_compared"] >= NUMBERING_MIN_COMPARED
                     and chk["fraction_match"] >= NUMBERING_MIN_MATCH)
    return chk


def sm_footprint(structure, receptor_chains, accession):
    """Best drug-like HET ligand on one of `receptor_chains`, and its footprint.

    "Best" = most contacting residues on the receptor chain, tie-broken by heavy
    atom count, so the deposit's main ligand wins over a stray fragment. Returns
    (record, None) or (None, reason) - never a silent skip.
    """
    inv = I.het_inventory(structure, min_heavy_atoms=MIN_HEAVY_ATOMS)
    kept = [r for r in inv if r["kept"]]
    # One copy per distinct component: ligand_contacts rebuilds a NeighborSearch
    # over the whole target for every call, and a deposit with 40 copies of one
    # ligand would pay that 40 times for the same answer. Copies are selected
    # inside ligand_contacts by chain/resseq, so the first copy of each component
    # is representative; MAX_LIGAND_COMPONENTS bounds the work on a deposit that
    # really does carry many different drug-like components.
    seen_components = set()
    deduped = []
    for r in kept:
        if r["resname"] in seen_components:
            continue
        seen_components.add(r["resname"])
        deduped.append(r)
    n_components = len(deduped)
    kept = deduped[:MAX_LIGAND_COMPONENTS]
    if not kept:
        reasons = Counter(r["reason"].split(":")[0] for r in inv if not r["kept"])
        return None, (f"no HET residue with >= {MIN_HEAVY_ATOMS} heavy atoms "
                      f"survived the ligand filter ({len(inv)} HET residues, "
                      f"{dict(reasons)})")

    best = None
    for lig in kept:
        try:
            c = I.ligand_contacts(structure, lig["resname"], chain=lig["chain_id"],
                                  resseq=lig["resseq"], cutoff=CUTOFF,
                                  with_buried_area=False)
        except (KeyError, ValueError):
            continue
        collision = [w for w in c["warnings"] if "residue-id collision" in w]
        for rc in receptor_chains:
            ids = sorted(rid for rid, rec in c["per_residue"].items()
                         if rec["chain_id"] == rc)
            if len(ids) < MIN_CONTACTS:
                continue
            if collision:
                continue  # ids would merge two chains' residues; not comparable
            key = (len(ids), lig["n_heavy_atoms"])
            if best is None or key > best[0]:
                best = (key, {
                    "ligand": lig["resname"],
                    "ligand_chain": lig["chain_id"],
                    "ligand_resseq": lig["resseq"],
                    "ligand_n_heavy_atoms": lig["n_heavy_atoms"],
                    "ligand_formula": lig["formula"],
                    "receptor_chain": rc,
                    "residue_ids": ids,
                    "n_chains_engaged": len(c["chains_engaged"]),
                    "contacts_restricted_to_receptor_chain":
                        c["chains_engaged"] != [rc],
                })
    if best is None:
        return None, (f"{n_components} drug-like HET component(s) present but none "
                      f"made >= {MIN_CONTACTS} contacts on a chain of "
                      f"{accession} without a cross-chain residue-id collision")
    return best[1], None


def polymer_footprint(structure, receptor_chain, binder_chains):
    """Union of the target-side footprints of every binder chain.

    An antibody binder is a heavy AND a light chain; the epitope is what the Fab
    touches, so the two chains' target-side contacts are unioned. Same
    chain_contacts call, same cutoff, same numbering as the ligand path.
    """
    ids, per_chain = set(), {}
    for bc in binder_chains:
        c = I.chain_contacts(structure, receptor_chain, bc, cutoff=CUTOFF,
                             with_buried_area=False)
        per_chain[bc] = len(c["residue_ids"])
        ids.update(c["residue_ids"])
    return sorted(ids), per_chain


def auto_polymer_side(pdb_id, accession, cache_dir):
    """Receptor and binder chains for a hand-seeded entry, read off the coordinates."""
    st = I.load_structure(_fetch(pdb_id, cache_dir), pdb_id)
    cands = [c for c in I.candidate_binder_chains(st, cutoff=CUTOFF)
             if c["is_candidate"]]
    if not cands:
        return None, f"{pdb_id}: no candidate binder chain in the coordinates"
    best = None
    for c in cands:
        for rc, n in sorted(c["contacts_by_chain"].items(), key=lambda kv: -kv[1]):
            chk = numbering_ok(st, rc, accession)
            if not chk["passed"]:
                continue
            if best is None or n > best[0]:
                best = (n, {"pdb_id": pdb_id, "receptor_chain": rc,
                            "binder_chains": [c["chain_id"]],
                            "binder_seq_len": c["n_residues"],
                            "binder_description": f"chain {c['chain_id']} "
                                                  f"({c['modality']}, "
                                                  f"{c['n_residues']} residues)",
                            "receptor_uniprot": accession,
                            "resolution_a": None})
    if best is None:
        return None, (f"{pdb_id}: no chain pairs a binder with a receptor chain "
                      f"whose author numbering matches {accession}")
    return best[1], None


def build_pairs(max_pairs=30, refresh=False, verbose=True):
    """The cross-modality gold set: one (small molecule, biologic) pair per target."""
    dest = CACHE / f"pairs_max{max_pairs}.json"
    if dest.exists() and not refresh:
        cached = json.loads(dest.read_text())
        if cached.get("complete"):
            return cached
        print(f"  (cache {dest.name} is a partial run; rebuilding)", flush=True)

    STRUCTURES.mkdir(parents=True, exist_ok=True)
    poly_by_acc = load_polymer_side()
    accessions = sorted(poly_by_acc)

    pairs, rejected = [], []
    searched = 0

    def reject(**kw):
        rejected.append(kw)
        if verbose:
            print(f"  -- {kw.get('uniprot')} [{kw.get('stage')}] "
                  f"{str(kw.get('reason'))[:120]}", flush=True)

    work = [(s["uniprot"], s) for s in SEED_PAIRS] + [(a, None) for a in accessions]
    for accession, seed in work:
        if len(pairs) >= max_pairs:
            break
        # ---- polymer side
        if seed is not None:
            poly, why = auto_polymer_side(seed["poly_pdb"], accession, STRUCTURES)
            if poly is None:
                reject(**{"uniprot": accession, "stage": "polymer_side",
                                 "reason": why})
                continue
            poly = {**poly, "klass": seed["klass"], "note": seed["note"]}
            poly_path = _fetch(poly["pdb_id"], STRUCTURES)
        else:
            poly = poly_by_acc[accession]
            poly_path = EPITOPE_CACHE / "structures" / f"{poly['pdb_id'].lower()}.cif"
            if not poly_path.exists():
                poly_path = _fetch(poly["pdb_id"], STRUCTURES)

        # ---- small-molecule side: search, then verify in the coordinates
        t_acc = time.perf_counter()
        ids, total = sm_entry_candidates(accession)
        searched += 1
        ids = [i for i in ids][:MAX_SM_TRIES]
        if not ids:
            reject(**{"uniprot": accession, "stage": "sm_search",
                             "reason": "no X-ray entry <= "
                                       f"{SM_RESOLUTION_MAX} A with a nonpolymer "
                                       "entity and a chain of this accession"})
            continue
        entries = {e["rcsb_id"]: e for e in E.graphql_entries(ids)}

        sm_rec, sm_path, sm_why = None, None, None
        for pid in ids:
            meta = entries.get(pid.upper())
            rchains = _uniprot_chains(meta, accession) if meta else []
            if not rchains:
                sm_why = f"{pid}: GraphQL lists no auth chain for {accession}"
                continue
            try:
                sm_path = _fetch(pid, STRUCTURES)
                st_sm = I.load_structure(sm_path, pid)
            except Exception as exc:                        # noqa: BLE001
                sm_why = f"{pid}: {type(exc).__name__}: {exc}"
                continue
            ok_chains, checks = [], []
            for rc in rchains:
                if rc not in {c.id for c in st_sm[0]}:
                    continue
                chk = numbering_ok(st_sm, rc, accession)
                checks.append(chk)
                if chk["passed"]:
                    ok_chains.append(rc)
            if not ok_chains:
                sm_why = (f"{pid}: author numbering does not index {accession} on "
                          f"any chain ({[(c['chain_id'], round(c['fraction_match'], 3), c['n_compared']) for c in checks]})")
                continue
            rec, why = sm_footprint(st_sm, ok_chains, accession)
            if rec is None:
                sm_why = f"{pid}: {why}"
                continue
            res = (meta.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            sm_rec = {**rec, "pdb_id": pid.upper(),
                      "resolution_a": res[0] if res else None,
                      "title": (meta.get("struct") or {}).get("title"),
                      "numbering_check": [c for c in checks
                                          if c["chain_id"] == rec["receptor_chain"]][0],
                      "n_entries_available": total}
            break
        if sm_rec is None:
            reject(**{"uniprot": accession, "stage": "sm_side",
                             "reason": sm_why or "no candidate entry usable",
                             "tried": ids})
            continue

        # ---- polymer side footprint, same coordinates convention
        st_poly = I.load_structure(poly_path, poly["pdb_id"])
        chk_poly = numbering_ok(st_poly, poly["receptor_chain"], accession)
        if not chk_poly["passed"]:
            reject(**{
                "uniprot": accession, "stage": "polymer_numbering",
                "reason": f"{poly['pdb_id']} chain {poly['receptor_chain']}: author "
                          f"numbering matches {accession} on only "
                          f"{chk_poly['fraction_match']:.3f} of "
                          f"{chk_poly['n_compared']} residues; ids are not "
                          "comparable with the small-molecule entry"})
            continue
        try:
            poly_ids, per_chain = polymer_footprint(
                st_poly, poly["receptor_chain"], poly["binder_chains"])
        except Exception as exc:                             # noqa: BLE001
            reject(**{"uniprot": accession, "stage": "polymer_contacts",
                             "reason": f"{type(exc).__name__}: {exc}"})
            continue
        if len(poly_ids) < MIN_CONTACTS:
            reject(**{"uniprot": accession, "stage": "polymer_contacts",
                             "reason": f"{len(poly_ids)} target residues < "
                                       f"{MIN_CONTACTS}"})
            continue

        row = {
            "uniprot": accession,
            "klass": poly["klass"],
            "receptor_description": poly.get("receptor_description"),
            "same_entry": sm_rec["pdb_id"].upper() == poly["pdb_id"].upper(),
            "small_molecule": {
                "pdb_id": sm_rec["pdb_id"], "resolution_a": sm_rec["resolution_a"],
                "ligand": sm_rec["ligand"],
                "ligand_n_heavy_atoms": sm_rec["ligand_n_heavy_atoms"],
                "ligand_formula": sm_rec["ligand_formula"],
                "receptor_chain": sm_rec["receptor_chain"],
                "residue_ids": sm_rec["residue_ids"],
                "numbering_fraction_match": sm_rec["numbering_check"]["fraction_match"],
                "numbering_n_compared": sm_rec["numbering_check"]["n_compared"],
                "contacts_restricted_to_receptor_chain":
                    sm_rec["contacts_restricted_to_receptor_chain"],
                "title": sm_rec["title"],
            },
            "biologic": {
                "pdb_id": poly["pdb_id"], "resolution_a": poly.get("resolution_a"),
                "binder_description": poly.get("binder_description"),
                "binder_seq_len": poly.get("binder_seq_len"),
                "receptor_chain": poly["receptor_chain"],
                "binder_chains": poly["binder_chains"],
                "residues_per_binder_chain": per_chain,
                "residue_ids": poly_ids,
                "numbering_fraction_match": chk_poly["fraction_match"],
                "numbering_n_compared": chk_poly["n_compared"],
            },
            "metrics": metric_set(sm_rec["residue_ids"], poly_ids),
        }
        if poly.get("note"):
            row["note"] = poly["note"]
        pairs.append(row)
        CACHE.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(
            {"pairs": pairs, "rejected": rejected, "n_targets_searched": searched,
             "n_targets_offered": len(accessions) + len(SEED_PAIRS),
             "complete": False}, indent=1))
        if verbose:
            m = row["metrics"]
            print(f"  [{len(pairs):2d}] {accession} {row['klass']:<8} "
                  f"SM {row['small_molecule']['pdb_id']}/{row['small_molecule']['ligand']:<4} "
                  f"vs {row['biologic']['pdb_id']}  "
                  f"|A|={m['n_a']:3d} |B|={m['n_b']:3d} "
                  f"J={m['jaccard']:.3f} cov={m['coverage_a']:.3f} "
                  f"ovl={m['overlap_coefficient']:.3f} "
                  f"({time.perf_counter() - t_acc:.0f}s)", flush=True)

    out = {"pairs": pairs, "rejected": rejected,
           "n_targets_searched": searched,
           "n_targets_offered": len(accessions) + len(SEED_PAIRS),
           "complete": True}
    CACHE.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    return out


# --------------------------------------------------------------------------
# 3. statistics on the real set
# --------------------------------------------------------------------------
def _spearman(x, y, n_boot=N_BOOT, seed=SEED):
    """Spearman rho with n, p and a percentile bootstrap CI over the pairs.

    metrics.py owns the bootstrap machinery this repo uses, but it resamples
    MEANS of a vector; a correlation needs the two vectors resampled together,
    so scipy's paired bootstrap is used (same percentile method, same n_boot).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    n = int(x.size)
    if n < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return {"rho": float("nan"), "p_value": float("nan"), "n": n,
                "ci_lo": float("nan"), "ci_hi": float("nan"),
                "note": "constant or too-small vector; correlation undefined"}
    res = stats.spearmanr(x, y)
    try:
        boot = stats.bootstrap(
            (x, y), lambda a, b: stats.spearmanr(a, b).statistic,
            paired=True, vectorized=False, n_resamples=n_boot,
            random_state=np.random.default_rng(seed), method="percentile")
        lo, hi = float(boot.confidence_interval.low), float(boot.confidence_interval.high)
    except Exception:                                        # noqa: BLE001
        lo = hi = float("nan")
    return {"rho": float(res.statistic), "p_value": float(res.pvalue),
            "n": n, "ci_lo": lo, "ci_hi": hi}


def _partial_spearman(x, y, z, n_boot=N_BOOT, seed=SEED):
    """Spearman correlation of x and y with z partialled out, on rank residuals.

    The confound this controls: pairs whose footprints are more similar in SIZE
    may also genuinely overlap MORE, in which case Jaccard tracking the size
    ratio would be biology, not metric bias. Holding coverage (the containment
    that actually happened) fixed separates the two.

    Computed as the Pearson correlation of the residuals of rank(x) on rank(z)
    and rank(y) on rank(z) - the standard partial Spearman, written out rather
    than taken from the three-rho closed form so the SAME residuals can carry
    the test. The p value is a permutation test on those residuals (shuffling
    one of them breaks any association while keeping both marginals), which
    needs no normality assumption, and the CI is a percentile bootstrap over
    pairs with non-finite resamples dropped and counted.
    """
    x, y, z = (np.asarray(v, dtype=np.float64) for v in (x, y, z))
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[keep], y[keep], z[keep]
    n = int(x.size)
    if n < 5:
        return {"partial_rho": float("nan"), "n": n,
                "note": "fewer than 5 complete rows"}

    def resid_corr(a, b, c):
        ra, rb, rc = (stats.rankdata(v) for v in (a, b, c))
        if np.ptp(rc) == 0:
            return float("nan")
        ea = ra - np.polyval(np.polyfit(rc, ra, 1), rc)
        eb = rb - np.polyval(np.polyfit(rc, rb, 1), rc)
        if np.std(ea) == 0 or np.std(eb) == 0:
            return float("nan")
        return float(np.corrcoef(ea, eb)[0, 1])

    rho = resid_corr(x, y, z)

    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rng = np.random.default_rng(seed)
    obs = abs(float(np.corrcoef(ex, ey)[0, 1]))
    hits = 0
    for _ in range(n_boot):
        if abs(float(np.corrcoef(ex, rng.permutation(ey))[0, 1])) >= obs:
            hits += 1
    p_perm = (hits + 1) / (n_boot + 1)

    rng = np.random.default_rng(seed + 1)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        v = resid_corr(x[idx], y[idx], z[idx])
        if math.isfinite(v):
            draws.append(v)
    if draws:
        lo, hi = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    else:
        lo = hi = float("nan")
    return {"partial_rho": rho, "n": n, "ci_lo": lo, "ci_hi": hi,
            "p_permutation": float(p_perm),
            "n_bootstrap_resamples_finite": len(draws),
            "n_bootstrap_resamples_requested": n_boot,
            "controlled_for": "coverage_a",
            "method": "Pearson correlation of rank residuals; permutation p; "
                      "percentile bootstrap CI over pairs"}


def analyse_pairs(pairs):
    """Does Jaccard track the size ratio where coverage and overlap do not?"""
    if len(pairs) < 3:
        return {"n": len(pairs), "note": "too few pairs to correlate"}
    ratio = [p["metrics"]["size_ratio"] for p in pairs]
    cols = {k: [p["metrics"][k] for p in pairs]
            for k in ("jaccard", "coverage_a", "coverage_b",
                      "overlap_coefficient", "f1")}

    corr = {k: _spearman(cols[k], ratio) for k in cols}
    finite = {k: v["p_value"] for k, v in corr.items() if math.isfinite(v["p_value"])}
    holm = M.holm_bonferroni(finite) if finite else {}
    for k, v in corr.items():
        v["p_holm"] = holm.get(k)

    overlapping = [p for p in pairs if p["metrics"]["n_intersection"] > 0]
    corr_overlapping = {}
    if len(overlapping) >= 3:
        r2 = [p["metrics"]["size_ratio"] for p in overlapping]
        for k in cols:
            corr_overlapping[k] = _spearman(
                [p["metrics"][k] for p in overlapping], r2)

    # What a PERFECTLY contained footprint can score on this set: the size
    # ratio is the Jaccard ceiling, so its mean is the mean best-case Jaccard a
    # cross-modality comparison could reach even with containment perfect.
    ratio_mean, ratio_lo, ratio_hi = M.bootstrap_ci(np.array(ratio, float),
                                                    n_boot=N_BOOT, seed=SEED)
    def contained_at(thr):
        sel = [p for p in pairs if p["metrics"]["coverage_a"] >= thr]
        return {
            "coverage_threshold": thr,
            "n": len(sel),
            "pairs": [f"{p['small_molecule']['pdb_id']}/"
                      f"{p['small_molecule']['ligand']} vs "
                      f"{p['biologic']['pdb_id']}" for p in sel],
            "jaccard_values": sorted(round(p["metrics"]["jaccard"], 4)
                                     for p in sel),
            "mean_jaccard": (float(np.mean([p["metrics"]["jaccard"]
                                            for p in sel])) if sel else None),
            "mean_coverage_a": (float(np.mean([p["metrics"]["coverage_a"]
                                               for p in sel])) if sel else None),
        }
    contained_block = {
        "definition": "pairs where at least this fraction of the small-molecule "
                      "footprint lies inside the biologic footprint; coverage "
                      "calls them containments, Jaccard scores them as partial "
                      "overlaps",
        "by_threshold": [contained_at(t) for t in (1.0, 0.9, 0.8)],
    }
    contained_block["headline"] = next(
        (b for b in contained_block["by_threshold"] if b["n"] >= 3),
        contained_block["by_threshold"][1])
    # What this n could have detected. A null correlation is only informative
    # against the effects the sample had the power to see, so the smallest rho
    # detectable at 80% power (Fisher z, two-sided alpha=0.05) is computed and
    # reported next to every null below rather than left for the reader to
    # assume.
    def min_detectable_rho(n):
        return (float(np.tanh((1.959963985 + 0.8416212336) / math.sqrt(n - 3)))
                if n > 3 else float("nan"))

    return {
        "n_pairs": len(pairs),
        "n_pairs_with_any_overlap": len(overlapping),
        "power": {
            "min_detectable_rho_all_pairs": min_detectable_rho(len(pairs)),
            "min_detectable_rho_overlapping_pairs": min_detectable_rho(
                len(overlapping)),
            "definition": "smallest |rho| detectable at 80% power, two-sided "
                          "alpha=0.05, Fisher z approximation; a null below "
                          "this is 'not detected at this n', not 'absent'",
        },
        "size_ratio_is_the_jaccard_ceiling": {
            "mean_size_ratio": ratio_mean, "ci": [ratio_lo, ratio_hi],
            "n": len(ratio),
            "meaning": "the largest Jaccard these pairs could score if the "
                       "smaller footprint were entirely inside the larger",
        },
        "perfectly_contained_pairs": contained_block,
        "size_ratio_definition": "min(|A|,|B|) / max(|A|,|B|), A = small-molecule "
                                 "footprint, B = biologic footprint",
        "spearman_vs_size_ratio_all_pairs": corr,
        "spearman_vs_size_ratio_overlapping_pairs_only": corr_overlapping,
        "jaccard_vs_coverage": _spearman(cols["jaccard"], cols["coverage_a"]),
        "partial_jaccard_vs_ratio_given_coverage": _partial_spearman(
            cols["jaccard"], ratio, cols["coverage_a"]),
        "f1_is_a_monotone_function_of_jaccard": {
            "max_abs_deviation_from_2J_over_1_plus_J": max(
                abs(p["metrics"]["f1"] - 2 * p["metrics"]["jaccard"]
                    / (1 + p["metrics"]["jaccard"]))
                for p in pairs if math.isfinite(p["metrics"]["jaccard"])),
            "meaning": "F1 and Jaccard have identical rank order, so F1 inherits "
                       "every Jaccard failure and adds no information",
        },
    }


def containment_examples(pairs, coverage_hi=0.7, coverage_lo=0.45):
    """Real pairs where Jaccard and coverage disagree, named so they can be checked.

    `contained_but_low_jaccard`: the small-molecule site is mostly INSIDE the
    biologic footprint (coverage >= coverage_hi) yet Jaccard is low. If such
    pairs exist and their Jaccards sit in the same range as the genuinely
    partial pairs below, Jaccard is provably ambiguous on real structures.
    """
    hi = sorted((p for p in pairs if p["metrics"]["coverage_a"] >= coverage_hi),
                key=lambda p: p["metrics"]["jaccard"])
    lo = sorted((p for p in pairs if p["metrics"]["n_intersection"] > 0
                 and p["metrics"]["coverage_a"] <= coverage_lo),
                key=lambda p: -p["metrics"]["jaccard"])

    def brief(p):
        m = p["metrics"]
        return {"uniprot": p["uniprot"], "klass": p["klass"],
                "sm": f"{p['small_molecule']['pdb_id']}/{p['small_molecule']['ligand']}",
                "biologic": p["biologic"]["pdb_id"],
                "n_a": m["n_a"], "n_b": m["n_b"],
                "n_intersection": m["n_intersection"],
                "jaccard": m["jaccard"], "coverage_a": m["coverage_a"],
                "overlap_coefficient": m["overlap_coefficient"],
                "size_ratio": m["size_ratio"]}

    contained = [brief(p) for p in hi]
    partial = [brief(p) for p in lo]
    collisions = []
    for c in contained:
        for q in partial:
            if abs(c["jaccard"] - q["jaccard"]) <= 0.05:
                collisions.append({
                    "jaccard_gap": abs(c["jaccard"] - q["jaccard"]),
                    "contained": c, "partial": q,
                    "coverage_gap": c["coverage_a"] - q["coverage_a"]})
    collisions.sort(key=lambda c: c["jaccard_gap"])
    gaps = [abs(c["jaccard"] - q["jaccard"]) for c in contained for q in partial]
    return {
        "nearest_contained_to_partial_jaccard_gap": min(gaps) if gaps else None,
        "n_contained": len(contained), "n_partial": len(partial),
        "coverage_thresholds": {"contained_at_least": coverage_hi,
                                "partial_at_most": coverage_lo},
        "contained_but_low_jaccard": contained,
        "genuinely_partial": partial,
        "jaccard_collisions": collisions[:10],
        "n_jaccard_collisions": len(collisions),
    }


# --------------------------------------------------------------------------
# 4. what this means for the gate numbers already in results/
# --------------------------------------------------------------------------
def jaccard_ceiling():
    """The best Jaccard a perfectly-placed prediction of the emitted SIZE could get.

    For a prediction of n_pred residues scored against a truth set of n_truth,
    the maximum possible Jaccard is min(n_pred, n_truth) / max(n_pred, n_truth) -
    attained when the smaller set is entirely inside the larger. That is an
    arithmetic ceiling imposed by SIZE alone. Comparing it to the Jaccard
    actually scored says how much of the small-molecule/biologic gap is metric
    artefact and how much is a real placement failure.

    Inputs are the n_pred and interface-size columns of results/epitope_gate.json
    and results/m2_gate_KDR.json; the ceiling itself is computed here.
    """
    out = {}
    if EPITOPE_RESULT.exists():
        rows = json.loads(EPITOPE_RESULT.read_text())["per_structure"]
        by_class = defaultdict(list)
        dropped = Counter()
        total = Counter()
        for r in rows:
            n_pred = r["p2rank_top1"]["n_pred"]
            n_truth = r["n_interface_residues"]
            total[r["klass"]] += 1
            if not n_pred or not n_truth:
                # P2Rank returned no top-1 pocket at all: there is no predicted
                # SIZE, so a size-imposed ceiling is undefined. These rows are
                # counted, not silently dropped - they score Jaccard 0 in the
                # published class means, which is why the observed mean below is
                # higher than the published one.
                dropped[r["klass"]] += 1
                continue
            by_class[r["klass"]].append({
                "ceiling": min(n_pred, n_truth) / max(n_pred, n_truth),
                "observed": r["p2rank_top1"]["jaccard"],
                "n_pred": n_pred, "n_truth": n_truth})
        for k, rs in sorted(by_class.items()):
            ceil = np.array([x["ceiling"] for x in rs])
            obs = np.array([x["observed"] for x in rs])
            mean, lo, hi = M.bootstrap_ci(ceil, n_boot=N_BOOT, seed=SEED)
            out[k] = {
                "n": len(rs),
                "mean_jaccard_ceiling": mean,
                "ceiling_ci": [lo, hi],
                "mean_jaccard_observed": float(obs.mean()),
                "mean_fraction_of_ceiling_attained": float(
                    np.mean(obs / ceil)),
                "mean_n_pred": float(np.mean([x["n_pred"] for x in rs])),
                "mean_n_truth": float(np.mean([x["n_truth"] for x in rs])),
                "n_rows_in_class": int(total[k]),
                "n_rows_without_a_top1_prediction": int(dropped[k]),
                "observed_mean_is_over": "the rows with a top-1 pocket only; the "
                                         "published class mean includes the "
                                         "no-pocket rows as Jaccard 0",
                "source": "results/epitope_gate.json per_structure "
                          "(p2rank_top1.n_pred, n_interface_residues)",
            }
    if KDR_RESULT.exists():
        rows = json.loads(KDR_RESULT.read_text())["per_structure"]
        rs = [{"ceiling": min(r["p2rank_geometry"]["n_pred"], r["n_known_contacts"])
                          / max(r["p2rank_geometry"]["n_pred"], r["n_known_contacts"]),
               "observed": r["p2rank_geometry"]["jaccard"],
               "n_pred": r["p2rank_geometry"]["n_pred"],
               "n_truth": r["n_known_contacts"]}
              for r in rows
              if r["p2rank_geometry"]["n_pred"] and r["n_known_contacts"]]
        ceil = np.array([x["ceiling"] for x in rs])
        obs = np.array([x["observed"] for x in rs])
        mean, lo, hi = M.bootstrap_ci(ceil, n_boot=N_BOOT, seed=SEED)
        out["small_molecule_KDR"] = {
            "n": len(rs), "mean_jaccard_ceiling": mean, "ceiling_ci": [lo, hi],
            "mean_jaccard_observed": float(obs.mean()),
            "mean_fraction_of_ceiling_attained": float(np.mean(obs / ceil)),
            "mean_n_pred": float(np.mean([x["n_pred"] for x in rs])),
            "mean_n_truth": float(np.mean([x["n_truth"] for x in rs])),
            "source": "results/m2_gate_KDR.json per_structure "
                      "(p2rank_geometry.n_pred, n_known_contacts)",
        }
    return out


def metric_swap_on_existing_gates():
    """Re-score the gates already run under each metric, and re-test.

    The E3 question is not only "is Jaccard biased" but "does anything the repo
    has concluded change if the metric changes". The gate scripts already stored
    precision, recall, F1 and Jaccard for every structure, so the class means and
    the tests can be recomputed here under each metric without re-running
    P2Rank. Three things are asked of each metric:

      * cross-benchmark: the KDR small-molecule arm against each biologic class
        (independent samples -> epitope_gate.unpaired, bootstrap + Mann-Whitney);
      * within-class: P2Rank top-1 against the size-matched random control on
        the SAME structures (paired -> metrics.paired_bootstrap and
        metrics.wilcoxon), Holm-corrected over the whole class x metric family;
      * the ordering of the classes, which is the concavity claim.

    Note precision here is the coverage of the PREDICTION (|pred & truth|/|pred|)
    and recall the coverage of the TRUTH (|pred & truth|/|truth|) - the same two
    one-sided coverages this experiment measures on the cross-modality pairs.
    """
    if not (EPITOPE_RESULT.exists() and KDR_RESULT.exists()):
        return {"status": "skipped: gate results not present"}
    bio = json.loads(EPITOPE_RESULT.read_text())["per_structure"]
    kdr = json.loads(KDR_RESULT.read_text())["per_structure"]
    metrics = ("jaccard", "precision", "recall", "f1")

    by_class = defaultdict(list)
    for r in bio:
        by_class[r["klass"]].append(r)

    out = {"metric_definitions": {
        "precision": "|pred & truth| / |pred|  (coverage of the prediction)",
        "recall": "|pred & truth| / |truth|  (coverage of the truth)",
        "jaccard": "|pred & truth| / |pred | truth|",
        "f1": "harmonic mean of precision and recall = 2J/(1+J)"},
        "source": "per-structure scores stored by scripts/epitope_gate.py and "
                  "scripts/m2_gate.py; the means and tests below are computed here",
        "classes": {}, "cross_benchmark_vs_KDR": {}, "class_ordering": {}}

    within_p = {}
    for klass, rows in sorted(by_class.items()):
        rec = {"n": len(rows)}
        for m in metrics:
            a = np.array([r["p2rank_top1"][m] for r in rows], float)
            b = np.array([r["random"][m] for r in rows], float)
            mean, lo, hi = M.bootstrap_ci(a, n_boot=N_BOOT, seed=SEED)
            pb = M.paired_bootstrap(a, b, n_boot=N_BOOT, seed=SEED)
            wx = M.wilcoxon(a, b)
            rec[m] = {"p2rank_top1_mean": mean, "ci": [lo, hi],
                      "random_mean": float(b.mean()),
                      "delta_vs_random": pb["delta"],
                      "delta_ci": [pb["ci_lo"], pb["ci_hi"]],
                      "p_paired_bootstrap": pb["p_value"],
                      "p_wilcoxon": wx["p_value"],
                      "n_effective_wilcoxon": wx["n_effective"]}
            within_p[f"{klass}/{m}"] = wx["p_value"]
        out["classes"][klass] = rec
    holm = M.holm_bonferroni(within_p)
    for key, adj in holm.items():
        klass, m = key.split("/")
        out["classes"][klass][m]["p_wilcoxon_holm"] = adj

    cross_p = {}
    for klass, rows in sorted(by_class.items()):
        rec = {}
        for m in metrics:
            a = np.array([r["p2rank_geometry"][m] for r in kdr], float)
            b = np.array([r["p2rank_top1"][m] for r in rows], float)
            u = E.unpaired(a, b, seed=SEED)
            rec[m] = u
            cross_p[f"{klass}/{m}"] = u["p_mannwhitney"]
        out["cross_benchmark_vs_KDR"][klass] = rec
    holm_cross = M.holm_bonferroni(cross_p)
    for key, adj in holm_cross.items():
        klass, m = key.split("/")
        out["cross_benchmark_vs_KDR"][klass][m]["p_mannwhitney_holm"] = adj

    for m in metrics:
        order = sorted(by_class, key=lambda k: -float(np.mean(
            [r["p2rank_top1"][m] for r in by_class[k]])))
        out["class_ordering"][m] = order
    out["ordering_is_metric_invariant"] = len(
        {tuple(v) for v in out["class_ordering"].values()}) == 1
    return out


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-pairs", type=int, default=30)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--analytic-only", action="store_true")
    args = ap.parse_args()

    t0 = time.perf_counter()
    print("=" * 78)
    print("E3  metric bias: is Jaccard the wrong metric across modalities?")
    print("=" * 78)

    nested = analytic_nested_table()
    print("\n1. ANALYTIC - A (size a) entirely INSIDE B (size b=40). Containment "
          "is perfect\n   in every row; only the size ratio changes.\n")
    print(f"   {'a':>4}{'b':>5}{'ratio':>8}{'jaccard':>9}{'coverage':>10}"
          f"{'overlap':>9}{'f1':>8}")
    for r in nested:
        print(f"   {r['a']:>4}{r['b']:>5}{r['size_ratio']:>8.3f}"
              f"{r['jaccard']:>9.3f}{r['coverage_a']:>10.3f}"
              f"{r['overlap_coefficient']:>9.3f}{r['f1']:>8.3f}")

    ambiguity = analytic_ambiguity()
    print("\n   Same Jaccard, opposite geometry:\n")
    print(f"   {'jaccard':>8}   {'nested (coverage)':<34}{'partial (coverage)':<34}")
    for r in ambiguity:
        nd = f"{r['nested']['description']}  cov={r['nested']['coverage_a']:.2f}"
        pd = f"{r['partial']['description']}  cov={r['partial']['coverage_a']:.2f}"
        print(f"   {r['jaccard_target']:>8.2f}   {nd:<36}{pd:<36}")

    breaks = analytic_breakpoints()
    print("\n   Where each metric breaks:")
    for b in breaks:
        print(f"     {b['metric']:<21} {b['breaks_when']}")
        print(f"     {'':<21} {b['example']} -> J={b['values']['jaccard']:.3f} "
              f"cov={b['values']['coverage_a']:.3f} "
              f"ovl={b['values']['overlap_coefficient']:.3f} "
              f"F1={b['values']['f1']:.3f}")

    payload = {
        "experiment": "E3 metric bias - Jaccard vs coverage vs overlap "
                      "coefficient vs F1 for cross-modality site comparison",
        "question": "does Jaccard confound containment with partial overlap, so "
                    "that a small-molecule footprint number and a biologic "
                    "epitope number cannot be pooled?",
        "definitions": {
            "A": "small-molecule footprint: receptor residues within 4.5 A "
                 "(heavy atom) of the HET ligand (interfaces.ligand_contacts)",
            "B": "biologic footprint: receptor residues within 4.5 A (heavy "
                 "atom) of the polymer binder (interfaces.chain_contacts), "
                 "unioned over the binder's chains",
            "jaccard": "|A & B| / |A | B|",
            "coverage_a": "|A & B| / |A|",
            "overlap_coefficient": "|A & B| / min(|A|, |B|)",
            "f1": "2|A & B| / (|A| + |B|)",
        },
        "analytic": {
            "nested_containment_table": nested,
            "same_jaccard_opposite_geometry": ambiguity,
            "breakpoints": breaks,
        },
        "credits_spent": 0,
    }

    if args.analytic_only:
        payload["real_structures"] = {"status": "skipped (--analytic-only)"}
    else:
        print("\n2. REAL STRUCTURES - building the cross-modality gold set "
              "(RCSB, cached)\n")
        built = build_pairs(max_pairs=args.max_pairs, refresh=args.refresh)
        pairs = built["pairs"]
        stats_block = analyse_pairs(pairs)
        examples = containment_examples(pairs)

        print(f"\n   {len(pairs)} pairs over {built['n_targets_searched']} targets "
              f"searched of {built['n_targets_offered']} offered; "
              f"{len(built['rejected'])} rejected")
        print(f"   classes: {dict(Counter(p['klass'] for p in pairs))}")

        if pairs:
            print(f"\n   {'uniprot':<9}{'class':<9}{'SM':<12}{'biologic':<9}"
                  f"{'|A|':>5}{'|B|':>5}{'ratio':>7}{'J':>7}{'cov_a':>7}"
                  f"{'cov_b':>7}{'ovl':>7}{'F1':>7}")
            for q in sorted(pairs, key=lambda q: -q["metrics"]["jaccard"]):
                m = q["metrics"]
                sm = f"{q['small_molecule']['pdb_id']}/{q['small_molecule']['ligand']}"
                print(f"   {q['uniprot']:<9}{q['klass']:<9}{sm:<12}"
                      f"{q['biologic']['pdb_id']:<9}{m['n_a']:>5}{m['n_b']:>5}"
                      f"{m['size_ratio']:>7.3f}{m['jaccard']:>7.3f}"
                      f"{m['coverage_a']:>7.3f}{m['coverage_b']:>7.3f}"
                      f"{m['overlap_coefficient']:>7.3f}{m['f1']:>7.3f}")

        if len(pairs) >= 3:
            print("\n   Spearman against the size ratio "
                  f"(n={stats_block['n_pairs']}):\n")
            print(f"   {'metric':<22}{'rho':>8}{'p':>12}{'p_holm':>10}"
                  f"{'95% CI':>22}")
            for k, v in list(stats_block["spearman_vs_size_ratio_all_pairs"].items()):
                ci = f"[{v['ci_lo']:+.3f}, {v['ci_hi']:+.3f}]"
                ph = "  n/a" if v.get("p_holm") is None else f"{v['p_holm']:.3g}"
                print(f"   {k:<22}{v['rho']:>+8.3f}{v['p_value']:>12.3g}"
                      f"{ph:>10}{ci:>22}")
            pw = stats_block["power"]
            print(f"\n   smallest |rho| this n could detect at 80% power: "
                  f"{pw['min_detectable_rho_all_pairs']:.3f} (n="
                  f"{stats_block['n_pairs']})")
            ov = stats_block["spearman_vs_size_ratio_overlapping_pairs_only"]
            if ov:
                n_ov = stats_block["n_pairs_with_any_overlap"]
                print(f"\n   Same, on the {n_ov} pairs whose footprints overlap "
                      f"at all:\n")
                for k, v in ov.items():
                    ci = f"[{v['ci_lo']:+.3f}, {v['ci_hi']:+.3f}]"
                    print(f"   {k:<22}{v['rho']:>+8.3f}{v['p_value']:>12.3g}"
                          f"{'':>10}{ci:>22}")
            pj = stats_block["partial_jaccard_vs_ratio_given_coverage"]
            print(f"\n   Jaccard vs size ratio, coverage partialled out: "
                  f"rho={pj['partial_rho']:+.3f} "
                  f"[{pj['ci_lo']:+.3f}, {pj['ci_hi']:+.3f}] "
                  f"permutation p={pj.get('p_permutation', float('nan')):.3g} "
                  f"n={pj['n']}")
            jc = stats_block["jaccard_vs_coverage"]
            print(f"   Jaccard vs coverage:                           "
                  f"rho={jc['rho']:+.3f} "
                  f"[{jc['ci_lo']:+.3f}, {jc['ci_hi']:+.3f}] "
                  f"p={jc['p_value']:.3g} n={jc['n']}")

            print(f"\n   Jaccard collisions (a contained pair and a partial pair "
                  f"within 0.05 Jaccard): {examples['n_jaccard_collisions']}")
            for c in examples["jaccard_collisions"][:3]:
                a, b = c["contained"], c["partial"]
                print(f"     J={a['jaccard']:.3f} {a['sm']} vs {a['biologic']} "
                      f"cov={a['coverage_a']:.2f} (|A|={a['n_a']},|B|={a['n_b']})"
                      f"   ==   J={b['jaccard']:.3f} {b['sm']} vs {b['biologic']} "
                      f"cov={b['coverage_a']:.2f} (|A|={b['n_a']},|B|={b['n_b']})")

        payload["real_structures"] = {
            "gold_set": {
                "construction": "targets in the epitope-gate prepared set (one "
                                "entry per target, best resolution) that also "
                                "have an X-ray entry <= "
                                f"{SM_RESOLUTION_MAX} A carrying a HET ligand of "
                                f">= {MIN_HEAVY_ATOMS} heavy atoms on a chain of "
                                "the same UniProt accession",
                "numbering_check": "author numbering on BOTH receptor chains must "
                                   f"match the UniProt sequence on >= "
                                   f"{NUMBERING_MIN_MATCH} of >= "
                                   f"{NUMBERING_MIN_COMPARED} compared residues "
                                   "(interfaces.check_author_numbering)",
                "n_pairs": len(pairs),
                "n_targets_searched": built["n_targets_searched"],
                "n_targets_offered": built["n_targets_offered"],
                "class_counts": dict(Counter(p["klass"] for p in pairs)),
                "cutoff_a": CUTOFF,
            },
            "statistics": stats_block,
            "containment_examples": examples,
            "pairs": pairs,
            "not_evaluated": built["rejected"],
            "not_evaluated_by_stage": dict(Counter(
                r["stage"] for r in built["rejected"])),
        }

    print("\n3. WHAT THIS MEANS FOR THE GATE NUMBERS ALREADY MEASURED\n")
    ceiling = jaccard_ceiling()
    payload["jaccard_size_ceiling_on_existing_gates"] = ceiling
    print(f"   {'set':<22}{'n':>5}{'ceiling':>9}{'observed':>10}"
          f"{'attained':>10}{'n_pred':>8}{'n_truth':>9}")
    for k, v in ceiling.items():
        print(f"   {k:<22}{v['n']:>5}{v['mean_jaccard_ceiling']:>9.3f}"
              f"{v['mean_jaccard_observed']:>10.4f}"
              f"{v['mean_fraction_of_ceiling_attained']:>10.3f}"
              f"{v['mean_n_pred']:>8.1f}{v['mean_n_truth']:>9.1f}")

    swap = metric_swap_on_existing_gates()
    payload["metric_swap_on_existing_gates"] = swap
    if swap.get("classes"):
        print("\n   The same gate rows re-scored under each metric "
              "(P2Rank top-1 vs its size-matched random control):\n")
        print(f"   {'class':<10}{'n':>4}{'metric':>11}{'p2rank':>9}{'random':>9}"
              f"{'delta':>9}{'wilcoxon p':>12}{'p_holm':>10}")
        for klass, rec in swap["classes"].items():
            for m in ("jaccard", "precision", "recall", "f1"):
                v = rec[m]
                print(f"   {klass:<10}{rec['n']:>4}{m:>11}"
                      f"{v['p2rank_top1_mean']:>9.4f}{v['random_mean']:>9.4f}"
                      f"{v['delta_vs_random']:>+9.4f}{v['p_wilcoxon']:>12.3g}"
                      f"{v['p_wilcoxon_holm']:>10.3g}")
        print(f"\n   class ordering by metric: {swap['class_ordering']}")
        print(f"   ordering is metric-invariant: "
              f"{swap['ordering_is_metric_invariant']}")

    payload["verdict"] = build_verdict(payload)
    payload["not_evaluated"] = NOT_EVALUATED
    payload["caveats"] = CAVEATS
    payload["runtime_seconds"] = round(time.perf_counter() - t0, 1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=1))
    import textwrap
    print("\nVERDICT")
    print(textwrap.fill(payload["verdict"]["headline"], 76,
                        initial_indent="   ", subsequent_indent="   "))
    print("\n   findings")
    for line in payload["verdict"]["findings"]:
        print(textwrap.fill(line, 76, initial_indent="     - ",
                            subsequent_indent="       "))
    print("\n   recommendation")
    for line in payload["verdict"]["recommendation"]:
        print(textwrap.fill(line, 76, initial_indent="     - ",
                            subsequent_indent="       "))
    print(f"\nwrote {OUT.relative_to(ROOT)} in {payload['runtime_seconds']}s")


CAVEATS = [
    "The gold set is 30 targets, one pair each. A null correlation at n=30 "
    "rules out only effects above the minimum detectable rho reported in "
    "real_structures.statistics.power; it is not evidence of exactly zero.",
    "Half the pairs have footprints that do not overlap at all: the "
    "small-molecule site and the biologic site are simply different sites on "
    "the same protein (a nuclear-receptor ligand pocket versus its coactivator "
    "groove, for instance). They are kept, counted and reported, and the "
    "correlations are also given on the overlapping subset alone.",
    "The two footprints of a pair come from two crystals, so an induced-fit "
    "difference between them is inside the measured overlap. Same-entry pairs "
    "(flagged same_entry) are the exception and are not separated out.",
    "Targets are taken in accession order until the pair cap is reached, so "
    "the set is not a random sample of the PDB; the selection rule is "
    "independent of any overlap it could produce, but it is not a sample from "
    "which prevalence can be quoted.",
    "The Holm family for the gate re-scoring is the 12 class x metric "
    "comparisons made here. It is a different family from the one "
    "scripts/epitope_gate.py corrected over, so an adjusted p there and here "
    "are not the same number.",
]

NOT_EVALUATED = [
    {"item": "SIFTS residue-level mapping",
     "why": "pairs are compared in author numbering, verified per entry against "
            "the UniProt sequence; entries whose author numbering does not index "
            "UniProt are rejected rather than remapped, so the gold set is "
            "smaller than it could be but every id in it means one thing"},
    {"item": "conformational change between the two entries",
     "why": "the small-molecule and biologic footprints come from different "
            "crystals of the same protein; an induced-fit difference is part of "
            "the measured overlap and is not separated out here"},
    {"item": "weighted variants (interfaces.weighted_jaccard, coverage())",
     "why": "this experiment is about the set-level arithmetic that the gates "
            "report; the frequency-weighted signature metrics inherit the same "
            "size behaviour but were not measured on this set"},
    {"item": "Dice/Tversky and distance-based site similarity",
     "why": "Dice is F1, already shown to be a monotone function of Jaccard; "
            "Tversky needs two free parameters that would have to be fitted, "
            "which this set is too small to do honestly"},
]


def build_verdict(payload):
    """The recommendation, assembled from the numbers actually computed above.

    Every branch below is decided by a measured quantity, so if Jaccard turns
    out to be fine on this set the verdict says so instead of arguing.
    """
    real = payload.get("real_structures") or {}
    st = real.get("statistics") or {}
    ex = real.get("containment_examples") or {}
    corr_all = st.get("spearman_vs_size_ratio_all_pairs") or {}
    corr_ovl = st.get("spearman_vs_size_ratio_overlapping_pairs_only") or {}
    partial = st.get("partial_jaccard_vs_ratio_given_coverage") or {}
    ceiling = payload.get("jaccard_size_ceiling_on_existing_gates") or {}
    swap = payload.get("metric_swap_on_existing_gates") or {}

    def fmt(c, key):
        v = c.get(key) or {}
        if not v or not math.isfinite(v.get("rho", float("nan"))):
            return f"{key}: not estimable"
        holm = v.get("p_holm")
        tail = f", Holm p={holm:.3g}" if holm is not None else ""
        return (f"{key} rho={v['rho']:+.3f} [{v['ci_lo']:+.3f}, {v['ci_hi']:+.3f}] "
                f"(p={v['p_value']:.3g}{tail}, n={v['n']})")

    findings = []
    findings.append(
        "ANALYTIC: with A entirely inside B, coverage and the overlap "
        "coefficient are 1.000 at every size ratio while Jaccard equals the "
        "size ratio itself (0.050 at |A|=2/|B|=40, 0.300 at 12/40, 0.600 at "
        "24/40). Containment does not change down that column; only size does.")
    amb = payload["analytic"]["same_jaccard_opposite_geometry"]
    worst = max(amb, key=lambda r: r["coverage_a_difference"])
    findings.append(
        f"ANALYTIC: Jaccard {worst['jaccard_target']:.2f} is scored both by "
        f"{worst['nested']['description']} (coverage "
        f"{worst['nested']['coverage_a']:.2f}) and by "
        f"{worst['partial']['description']} (coverage "
        f"{worst['partial']['coverage_a']:.2f}) - identical Jaccard, coverage "
        f"differing by {worst['coverage_a_difference']:.2f}.")

    if corr_all:
        findings.append("REAL PAIRS (all): " + "; ".join(
            fmt(corr_all, k) for k in ("jaccard", "coverage_a",
                                       "overlap_coefficient")) + ".")
    if corr_ovl:
        findings.append("REAL PAIRS (pairs whose footprints overlap at all): "
                        + "; ".join(fmt(corr_ovl, k) for k in
                                    ("jaccard", "coverage_a",
                                     "overlap_coefficient")) + ".")
    if partial and math.isfinite(partial.get("partial_rho", float("nan"))):
        findings.append(
            f"REAL PAIRS: Jaccard vs size ratio with coverage partialled out, "
            f"rho={partial['partial_rho']:+.3f} "
            f"[{partial['ci_lo']:+.3f}, {partial['ci_hi']:+.3f}], "
            f"n={partial['n']} - the size dependence that remains after the "
            "containment actually observed is held fixed.")
    if ex:
        gap = ex.get("nearest_contained_to_partial_jaccard_gap")
        findings.append(
            f"REAL PAIRS: {ex.get('n_jaccard_collisions', 0)} collisions where a "
            f"contained pair (coverage >= "
            f"{ex['coverage_thresholds']['contained_at_least']}, n="
            f"{ex.get('n_contained')}) and a partial pair (coverage <= "
            f"{ex['coverage_thresholds']['partial_at_most']}, n="
            f"{ex.get('n_partial')}) score within 0.05 Jaccard of each other"
            + (f"; the closest such pair differs by {gap:.3f} Jaccard"
               if gap is not None else "") + ".")
    if ceiling:
        findings.append("GATES: size-imposed Jaccard ceiling " + "; ".join(
            f"{k} {v['mean_jaccard_ceiling']:.3f} (observed "
            f"{v['mean_jaccard_observed']:.4f}, {v['mean_fraction_of_ceiling_attained']:.3f} "
            f"of ceiling, n={v['n']})" for k, v in ceiling.items()) + ".")
    if swap.get("class_ordering"):
        findings.append(
            f"GATES: class ordering by metric {swap['class_ordering']}; "
            f"metric-invariant: {swap['ordering_is_metric_invariant']}.")

    # ---- decisions, each from a measured quantity
    def sig(c, key):
        v = c.get(key) or {}
        pv = v.get("p_holm", v.get("p_value"))
        return (pv is not None and math.isfinite(pv) and pv < 0.05
                and math.isfinite(v.get("rho", float("nan"))))

    jaccard_tracks_ratio = sig(corr_all, "jaccard") or sig(corr_ovl, "jaccard")
    coverage_tracks_ratio = sig(corr_all, "coverage_a") or sig(corr_ovl, "coverage_a")
    overlap_tracks_ratio = (sig(corr_all, "overlap_coefficient")
                            or sig(corr_ovl, "overlap_coefficient"))
    p_perm = partial.get("p_permutation")
    partial_holds = (p_perm is not None and math.isfinite(p_perm) and p_perm < 0.05)
    collisions = (ex.get("n_jaccard_collisions", 0) or 0) > 0

    contained = ((st.get("perfectly_contained_pairs") or {}).get("headline")
                 or {})
    penalty = (1.0 - contained["mean_jaccard"]
               if contained.get("mean_jaccard") is not None else None)
    jc = st.get("jaccard_vs_coverage") or {}
    rank_agrees = (math.isfinite(jc.get("rho", float("nan")))
                   and abs(jc["rho"]) >= 0.9 and jc.get("p_value", 1) < 0.05)
    ceil_vals = [v["mean_jaccard_ceiling"] for v in ceiling.values()]
    ceiling_spread = (max(ceil_vals) - min(ceil_vals)) if len(ceil_vals) > 1 else None

    if penalty is not None:
        findings.append(
            f"REAL PAIRS: the n={contained['n']} pairs whose small-molecule "
            f"footprint is >= {contained['coverage_threshold']:.0%} inside the "
            f"biologic footprint score Jaccard {contained['jaccard_values']} "
            f"(mean {contained['mean_jaccard']:.3f}, coverage "
            f"{contained['mean_coverage_a']:.3f}) - coverage calls these "
            f"containments, Jaccard marks them {penalty:.0%} wrong. "
            f"{contained['pairs']}")
    ratio_block = st.get("size_ratio_is_the_jaccard_ceiling") or {}
    if ratio_block:
        findings.append(
            f"REAL PAIRS: mean size ratio {ratio_block['mean_size_ratio']:.3f} "
            f"[{ratio_block['ci'][0]:.3f}, {ratio_block['ci'][1]:.3f}], n="
            f"{ratio_block['n']} - that IS the mean Jaccard ceiling for these "
            "cross-modality comparisons, reached only if containment is perfect.")
    if rank_agrees:
        findings.append(
            f"REAL PAIRS: Jaccard and coverage rank the pairs almost "
            f"identically (rho={jc['rho']:+.3f}, p={jc['p_value']:.3g}, "
            f"n={jc['n']}), so the problem is the SCALE Jaccard puts them on, "
            "not the order it puts them in.")

    scale_problem = ((penalty is not None and penalty >= 0.2)
                     or (ceiling_spread is not None and ceiling_spread >= 0.1)
                     or jaccard_tracks_ratio or collisions)

    if jaccard_tracks_ratio or collisions:
        headline = ("Jaccard is size-confounded on this set: it tracks the "
                    "footprint size ratio, and containment and partial overlap "
                    "land on the same value. Small-molecule and biologic "
                    "Jaccards are not on one scale.")
    elif scale_problem:
        headline = (
            "Jaccard's ORDER is fine, its SCALE is not. It did NOT track the "
            "size ratio on these pairs"
            + (f" (rho={corr_all['jaccard']['rho']:+.3f}, "
               f"p={corr_all['jaccard']['p_value']:.3g}, "
               f"n={corr_all['jaccard']['n']})" if corr_all.get("jaccard") else "")
            + ", and no contained/partial Jaccard collision occurred, so it is "
              "not mis-ranking sites. What it does do is cap a near-complete "
              "containment"
            + (f" at a mean Jaccard of {contained['mean_jaccard']:.3f} "
               f"(n={contained['n']} pairs at coverage >= "
               f"{contained['coverage_threshold']:.0%})"
               if contained.get("mean_jaccard") is not None else "")
            + (f" and give the repo's benchmark sets ceilings that differ by "
               f"{ceiling_spread:.3f} Jaccard" if ceiling_spread else "")
            + " - so a Jaccard measured against small-molecule footprints and "
              "one measured against epitopes are not comparable numbers, even "
              "though each is a valid ranking within its own set.")
    else:
        headline = ("Jaccard is fine on the evidence measured here: it does not "
                    "track the size ratio, no contained/partial collision "
                    "occurred, contained pairs are not materially penalised and "
                    "the benchmark sets' ceilings are level. No churn is "
                    "warranted.")

    primary, alongside, never = [], [], []
    if scale_problem:
        primary.append(
            "PRIMARY: report the two one-sided coverages as a pair - recall "
            "|pred & truth|/|truth| and precision |pred & truth|/|pred|. "
            "Neither works alone (recall is gamed by predicting the whole "
            "surface, precision by predicting one residue), and both are "
            "already stored per structure by epitope_gate.py and m2_gate.py, so "
            "this costs no re-run.")
        alongside.append(
            "ALONGSIDE: keep Jaccard for continuity with the published numbers, "
            "but never bare - print |pred|, |truth| and the size-imposed "
            "ceiling min(|pred|,|truth|)/max(|pred|,|truth|) beside it, so a "
            "reader sees how much of the score was reachable.")
        alongside.append(
            "DROP: F1 adds nothing. It is exactly 2J/(1+J) (verified on this "
            "set to within "
            + (f"{st['f1_is_a_monotone_function_of_jaccard']['max_abs_deviation_from_2J_over_1_plus_J']:.1e}"
               if st.get("f1_is_a_monotone_function_of_jaccard") else "machine precision")
            + "), so it has the same rank order and the same ceiling as "
              "Jaccard and is a second column saying one thing.")
        alongside.append(
            "USE WITH SIZES: the overlap coefficient |A & B|/min(|A|,|B|) is the "
            "right metric for 'is the smaller site inside the larger', but it "
            "saturates at 1.0 for any contained set however small, so it is "
            "only readable next to |A| and |B|.")
        never.append(
            "NEVER COMPARE ACROSS MODALITIES: a Jaccard (or F1) measured "
            "against small-molecule footprints against one measured against "
            "epitopes. Their size-imposed ceilings differ"
            + (f" by {ceiling_spread:.3f}" if ceiling_spread else "")
            + ", so the numbers must not be pooled, averaged or ranked "
              "together - quote each with its own ceiling or not at all.")
    else:
        primary.append("PRIMARY: keep Jaccard - it is not behaving differently "
                       "from coverage on the measured set.")
        alongside.append("ALONGSIDE: report |pred| and |truth| with it so the "
                         "arithmetic ceiling stays visible.")

    if not jaccard_tracks_ratio and corr_all.get("jaccard"):
        never.append(
            "NEGATIVE RESULT KEPT: the predicted size-ratio correlation did NOT "
            "appear on real structures"
            + (f" (Jaccard rho={corr_all['jaccard']['rho']:+.3f}, "
               f"p={corr_all['jaccard']['p_value']:.3g}, Holm p="
               f"{corr_all['jaccard'].get('p_holm', float('nan')):.3g}, "
               f"n={corr_all['jaccard']['n']}"
               + (f"; overlapping pairs only rho={corr_ovl['jaccard']['rho']:+.3f}, "
                  f"p={corr_ovl['jaccard']['p_value']:.3g}, "
                  f"n={corr_ovl['jaccard']['n']}" if corr_ovl.get("jaccard") else "")
               + ")" if corr_all.get("jaccard") else "")
            + ". The real footprints are close in size - small-molecule and "
              "biologic footprints in this set differ far less than the 12-vs-40 "
              "case the concern was posed with - so the ratio has little room "
              "to act. The arithmetic is still real; it is the SPREAD of sizes "
              "that is small here.")
    if partial_holds:
        straddles = (math.isfinite(partial.get("ci_lo", float("nan")))
                     and partial["ci_lo"] * partial["ci_hi"] <= 0)
        never.append(
            f"WITH CONTAINMENT HELD FIXED the size dependence does show: partial "
            f"rho={partial['partial_rho']:+.3f} "
            f"[{partial['ci_lo']:+.3f}, {partial['ci_hi']:+.3f}], permutation "
            f"p={p_perm:.3g}, n={partial['n']}."
            + (" THE TWO TESTS DISAGREE: the permutation p is below 0.05 while "
               "the bootstrap CI includes zero. Both are recorded; at this n "
               "the partial correlation is suggestive, not established."
               if straddles else ""))
    if coverage_tracks_ratio:
        alongside.append("CAUTION: coverage itself correlated with the size "
                         "ratio here, so read it with the sizes too.")
    if overlap_tracks_ratio:
        alongside.append("CAUTION: the overlap coefficient also correlated with "
                         "the size ratio here.")
    if swap.get("ordering_is_metric_invariant"):
        never.append(
            "NOTHING ALREADY CONCLUDED NEEDS RE-RUNNING: the gate's class "
            "ordering is identical under Jaccard, precision, recall and F1"
            + (f" {swap['class_ordering']['jaccard']}" if swap.get("class_ordering") else "")
            + ", so the concavity result is metric-invariant and changing the "
              "metric does not rescue P2Rank on epitopes.")
    elif swap.get("class_ordering"):
        never.append("RE-EXAMINE: the gate's class ordering is NOT the same "
                     "under every metric - " + str(swap["class_ordering"]))

    return {"headline": headline, "findings": findings,
            "recommendation": primary + alongside + never,
            "decision_inputs": {
                "jaccard_correlates_with_size_ratio": bool(jaccard_tracks_ratio),
                "coverage_correlates_with_size_ratio": bool(coverage_tracks_ratio),
                "overlap_coefficient_correlates_with_size_ratio":
                    bool(overlap_tracks_ratio),
                "partial_correlation_significant": bool(partial_holds),
                "n_jaccard_collisions_on_real_structures":
                    ex.get("n_jaccard_collisions", 0),
                "containment_penalty_1_minus_mean_jaccard_of_contained_pairs":
                    penalty,
                "jaccard_and_coverage_rank_agreement_rho": jc.get("rho"),
                "gate_ceiling_spread": ceiling_spread,
                "scale_problem": bool(scale_problem)}}


if __name__ == "__main__":
    main()
