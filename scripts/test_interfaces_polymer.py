"""Self-test for the POLYMER-BINDER path of scripts/interfaces.py.

`ligand_contacts` finds binders among the HET residues. A peptide or protein
binder is deposited as a polymer chain, so it is invisible there -- not an error,
a silent zero. `chain_contacts` is the polymer twin of that extractor, and this
file is the evidence that it works and that it agrees with the small-molecule
path it has to be comparable with.

No synthetic geometry: every number comes from a PDB entry that is downloaded,
verified to exist in RCSB, and parsed.

    1YCR  MDM2 / p53 transactivation-domain peptide   flat PPI epitope, 13-mer
    1CKA  Crk N-SH3 / proline-rich peptide            9-mer, PxxPxK class II
    1BRS  barnase / barstar                           87-residue PROTEIN binder
    6M0J  ACE2 / SARS-CoV-2 spike RBD                 194-residue protein binder,
                                                      the flat epitope this whole
                                                      formulation exists for
    1T4E  MDM2 / benzodiazepinedione (DIZ)            small molecule, SAME site
                                                      as the 1YCR peptide

Checks, in order:
   1. each entry exists in RCSB and the named chains are real polymer chains
   2. the old path is genuinely blind here: 1YCR/1CKA hold no drug-like HET
      ligand at all, and ligand_contacts raises rather than returning a zero
   3. candidate_binder_chains picks the short bound chain and rejects the target;
      a long protein binder is reachable only by naming it
   4. author numbering of 1YCR chain A is UniProt Q00987 numbering, so the ids
      are comparable with the small-molecule entry
   5. the 1YCR interface contains the known MDM2 p53-binding cleft residues; the
      two most-buried peptide residues are Phe19 and Trp23, and by fraction of
      free surface buried the top three are the Phe19/Trp23/Leu26 hotspot triad
      (absolute delta-SASA over-ranks a chain terminus; measured and kept)
   6. the 1CKA interface contains the SH3 aromatic cradle and the acidic
      specificity pocket that salt-bridges the peptide lysines
   7. the 1BRS interface contains the barnase active site and barstar Asp39
   8. symmetry: swapping target and binder mirrors residues, atom pairs,
      interaction types and buried areas exactly
   8b. wrong input (unknown chain, target as binder, a HET-only glycan chain)
      raises instead of returning an empty interface
   9. self-consistency: distances within cutoff, atom-contact counts add up, no
      binder-chain residue leaks into the target side
  10. cutoff monotonicity: contacts at 4.0 A are a subset of 4.5 A of 5.0 A
  11. with_buried_area=False reports None, never 0.0, and changes nothing else
  12. REGRESSION: every committed small-molecule contact map in
      results/interfaces_selftest.json is reproduced byte-for-byte, so the shared
      helpers this change introduced did not move the old path
  13. cross-modality: the 1YCR peptide epitope and the 1T4E small-molecule pocket
      are compared in one numbering with weighted_jaccard and coverage, and a
      mixed-modality consensus_signature is built

Usage:
    ./env/bin/python scripts/test_interfaces_polymer.py

Writes results/interfaces_polymer_selftest.json. Exits non-zero if any check fails.
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import interfaces as I  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "interfaces_polymer_selftest.json"
GOLDEN = ROOT / "results" / "interfaces_selftest.json"

CUTOFF = 4.5

# (pdb_id, target_chain, binder_chain, force_binder_by_name)
CASES = [
    ("1YCR", "A", "B", False),
    ("1CKA", "A", "B", False),
    ("1BRS", "A", "D", True),   # 87-residue protein binder: too long for the length rule
    ("6M0J", "A", "E", True),   # ACE2 / SARS-CoV-2 RBD: a biologic-scale flat epitope
]

MDM2_UNIPROT = "Q00987"
SM_ENTRY = ("1T4E", "DIZ", "A", 112)   # MDM2 + benzodiazepinedione, p53 cleft
# (pdb_id, target_chain, het_only_chain): NAG glycans deposited in their own
# chains, i.e. a chain that is NOT a polymer binder
HET_ONLY = ("3C45", "A", "C")

# Expected contacts. Each set is a published description of the site, and each is
# a SUBSET assertion: extra residues are reported, never required away.
#
# 1YCR: the MDM2 hydrophobic cleft that buries p53 Phe19/Trp23/Leu26
# (Kussie et al., Science 1996). 1CKA: the Crk N-SH3 PxxP groove -- the aromatic
# cradle Phe141/Trp169/Pro183/Tyr186 plus the acidic specificity pocket
# Asp147/Glu149/Asp150 that pairs with the peptide's C-terminal lysines
# (Wu et al., Structure 1995). 1BRS: the barnase active site Lys27/Arg59/Glu73/
# Arg87/His102 occluded by barstar, whose Asp39 is the interface hotspot
# (Buckle, Schreiber & Fersht, Biochemistry 1994).
# 6M0J: the ACE2 residues the RBD engages -- Gln24, Asp30, His34, Tyr41, Gln42,
# Lys353, Asp355 (Lan et al., Nature 2020), against RBD Lys417/Phe486/Gln493/
# Asn501/Tyr505, the positions the variants moved.
EXPECTED_TARGET = {
    "1YCR": {54, 58, 61, 62, 67, 72, 93, 96, 99, 100},
    "6M0J": {24, 30, 34, 41, 42, 353, 355},
    "1CKA": {141, 147, 149, 150, 169, 183, 186},
    "1BRS": {27, 59, 73, 87, 102},
}
EXPECTED_BINDER = {
    "1YCR": {19, 23, 26},     # Phe19, Trp23, Leu26
    "6M0J": {417, 486, 493, 501, 505},
    "1CKA": {2, 3, 6, 8},     # the PxxP prolines and the specificity lysine
    "1BRS": {35, 39},         # Asp35, Asp39
}
# The three peptide residues that should bury the most surface in 1YCR: the
# hotspot triad. Asserted as the measured top three, not merely as present.
MDM2_TRIAD = {19, 23, 26}

CHECKS = []
NOT_EVALUATED = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "passed": bool(ok), "detail": str(detail)})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))
    return bool(ok)


def rcsb_entry(pdb_id, cache_dir=I.RAW / "rcsb"):
    """Entry metadata from the RCSB data API, cached. Confirms the id exists."""
    dest = Path(cache_dir) / f"{pdb_id.upper()}.json"
    if not dest.exists():
        I._download(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}", dest)
    return json.loads(dest.read_text())


def load(pdb_id):
    rec = {"pdb_id": pdb_id.upper()}
    meta = rcsb_entry(pdb_id)
    rec["title"] = meta.get("struct", {}).get("title")
    rec["resolution_a"] = (meta.get("rcsb_entry_info", {})
                           .get("resolution_combined", [None])[0])
    rec["method"] = meta.get("rcsb_entry_info", {}).get("experimental_method")
    t0 = time.perf_counter()
    st = I.load_structure(I.fetch_structure(pdb_id), structure_id=pdb_id.upper())
    rec["parse_seconds"] = round(time.perf_counter() - t0, 3)
    return st, rec


def main():
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "contact_cutoff_a": CUTOFF,
        "numbering": "author (auth_seq_id); no SIFTS mapping applied",
        "cases": {},
        "timings_seconds": {},
    }

    structures = {}

    # -- 1..11: the polymer cases ------------------------------------------
    for pdb_id, tchain, bchain, forced in CASES:
        print(f"\n=== {pdb_id} target chain {tchain}, binder chain {bchain} ===")
        st, rec = load(pdb_id)
        structures[pdb_id] = st
        chains = {c["chain_id"]: c for c in I.polymer_chain_summary(st)}
        rec["chains"] = chains
        check(f"{pdb_id}: entry exists, {len(chains)} chains parsed",
              tchain in chains and bchain in chains,
              f"{rec['method']} {rec['resolution_a']} A; "
              f"target {tchain} n={chains.get(tchain, {}).get('n_residues')} res, "
              f"binder {bchain} n={chains.get(bchain, {}).get('n_residues')} res")
        check(f"{pdb_id}: binder chain {bchain} is a polymer chain, not a HET group",
              chains[bchain]["n_residues"] > 1,
              f"sequence {chains[bchain]['sequence']}")

        # 2. the old path is blind here
        kept = I.drug_like_ligands(st)
        rec["drug_like_het_ligands"] = kept
        if pdb_id in ("1YCR", "1CKA"):
            check(f"{pdb_id}: no drug-like HET ligand exists (absence verified)",
                  not kept, f"het_inventory kept {len(kept)} of "
                            f"{len(I.het_inventory(st))} non-polymer residues")
        # ask the HET path for the binder's own first residue: it is a polymer
        # residue, so is_polymer_residue excludes it and nothing is found
        first_res = next(r.get_resname().strip().upper()
                         for r in st[0][bchain] if I.is_polymer_residue(r))
        try:
            I.ligand_contacts(st, first_res)
            blind = False
            detail = f"ligand_contacts returned a site for {first_res}"
        except KeyError as exc:
            blind, detail = True, f"{first_res}: KeyError as expected ({str(exc)[:50]})"
        check(f"{pdb_id}: ligand_contacts cannot see the polymer binder (raises, "
              "does not return an empty site)", blind, detail)

        # 3. candidate_binder_chains
        t0 = time.perf_counter()
        cands = I.candidate_binder_chains(st, chain_ids=[bchain] if forced else None)
        dt_c = round(time.perf_counter() - t0, 3)
        rec["candidate_binder_chains"] = cands
        picked = [c["chain_id"] for c in cands if c["is_candidate"]]
        if forced:
            unforced = [c["chain_id"] for c in I.candidate_binder_chains(st)
                        if c["is_candidate"]]
            check(f"{pdb_id}: a {chains[bchain]['n_residues']}-residue protein binder "
                  "is NOT auto-picked and must be named",
                  unforced == [] and picked == [bchain],
                  f"auto={unforced}, named={picked}, {dt_c}s")
        else:
            check(f"{pdb_id}: candidate_binder_chains picks exactly [{bchain}]",
                  picked == [bchain],
                  f"{picked}; reason: "
                  f"{[c['reason'] for c in cands if c['is_candidate']][:1]}; {dt_c}s")

        # the contact map
        t0 = time.perf_counter()
        c = I.chain_contacts(st, tchain, bchain, cutoff=CUTOFF)
        dt = round(time.perf_counter() - t0, 3)
        t0 = time.perf_counter()
        c_nosasa = I.chain_contacts(st, tchain, bchain, cutoff=CUTOFF,
                                    with_buried_area=False)
        dt_nosasa = round(time.perf_counter() - t0, 3)
        report["timings_seconds"][f"{pdb_id}_chain_contacts"] = dt
        report["timings_seconds"][f"{pdb_id}_chain_contacts_no_sasa"] = dt_nosasa
        report["timings_seconds"][f"{pdb_id}_parse"] = rec["parse_seconds"]
        rec["contacts"] = {
            "binder_info": c["binder_info"],
            "n_target_residues": len(c["residue_ids"]),
            "residue_ids": c["residue_ids"],
            "per_residue": c["per_residue"],
            "binder_residue_ids": c["binder_residue_ids"],
            "binder_per_residue": c["binder_per_residue"],
            "n_atom_pairs": c["n_atom_pairs"],
            "warnings": c["warnings"],
        }

        # 5/6/7: the known site
        want_t, got_t = EXPECTED_TARGET[pdb_id], set(c["residue_ids"])
        check(f"{pdb_id}: target interface contains the published site residues",
              want_t <= got_t,
              f"missing {sorted(want_t - got_t)}; extracted n={len(got_t)}: "
              f"{sorted(got_t)}")
        want_b, got_b = EXPECTED_BINDER[pdb_id], set(c["binder_residue_ids"])
        check(f"{pdb_id}: binder side contains its published hotspots",
              want_b <= got_b,
              f"missing {sorted(want_b - got_b)}; engaged n={len(got_b)}: "
              f"{sorted(got_b)}")

        if pdb_id == "1YCR":
            bpr, free = c["binder_per_residue"], c["binder_sasa_free"]
            by_abs = sorted(bpr, key=lambda r: -bpr[r]["buried_area_proxy"])
            by_frac = sorted(bpr, key=lambda r: -bpr[r]["buried_area_proxy"] / free[r])
            res_b = [r for r in st[0][bchain] if I.is_polymer_residue(r)]
            termini = {res_b[0].id[1], res_b[-1].id[1]}
            # MEASURED NEGATIVE, kept rather than tuned away: by absolute
            # delta-SASA the third-ranked peptide residue is not Leu26 but the
            # C-terminal Asn29, because an isolated peptide's terminus has no
            # neighbour capping it and so is over-exposed in the reference state.
            # By fraction of its own free surface buried -- which that artifact
            # does not touch -- the top three ARE the triad.
            rank3_is_terminus = by_abs[2] in termini
            check("1YCR: the two most-buried peptide residues are Phe19 and Trp23",
                  set(by_abs[:2]) == {19, 23},
                  f"top2={sorted(by_abs[:2])} dSASA "
                  f"{[round(bpr[r]['buried_area_proxy'], 1) for r in by_abs[:2]]}")
            check("1YCR: by FRACTION of free surface buried, the top three are the "
                  "Phe19/Trp23/Leu26 triad",
                  set(by_frac[:3]) == MDM2_TRIAD,
                  f"top3={sorted(by_frac[:3])} at fractions "
                  f"{[round(bpr[r]['buried_area_proxy'] / free[r], 2) for r in by_frac[:3]]}")
            check("1YCR: absolute dSASA ranks a chain TERMINUS third, and the "
                  "reference-state SASA says why (measured caveat, not a bug)",
                  rank3_is_terminus,
                  f"rank3={by_abs[2]} ({bpr[by_abs[2]]['resname']}, terminus of "
                  f"{sorted(termini)}): dSASA "
                  f"{round(bpr[by_abs[2]]['buried_area_proxy'], 1)} but only "
                  f"{round(bpr[by_abs[2]]['buried_area_proxy'] / free[by_abs[2]], 2)} "
                  f"of its free surface, vs {round(bpr[26]['buried_area_proxy'] / free[26], 2)} "
                  "for Leu26")
            rec["terminal_dsasa_artifact"] = {
                "chain_termini": sorted(termini),
                "rank_by_absolute_dsasa": [
                    {"residue_id": r, "resname": bpr[r]["resname"],
                     "dsasa": bpr[r]["buried_area_proxy"],
                     "free_chain_sasa": free[r],
                     "buried_fraction": round(bpr[r]["buried_area_proxy"] / free[r], 3),
                     "is_terminus": r in termini}
                    for r in by_abs],
                "note": "absolute delta-SASA against an isolated-chain reference "
                        "over-ranks chain termini; buried_fraction does not",
            }
        if pdb_id in ("1CKA", "1BRS", "6M0J"):
            sb = sorted(r for r in c["residue_ids"]
                        if "salt_bridge" in c["per_residue"][r]["interaction_types"])
            check(f"{pdb_id}: the symmetric salt-bridge rule fires on the target side",
                  bool(sb), f"salt-bridged target residues {sb}")

        # 8. symmetry
        t0 = time.perf_counter()
        rev = I.chain_contacts(st, bchain, tchain, cutoff=CUTOFF)
        report["timings_seconds"][f"{pdb_id}_chain_contacts_reversed"] = round(
            time.perf_counter() - t0, 3)
        sym = {
            "target_ids_mirror": c["residue_ids"] == rev["binder_residue_ids"],
            "binder_ids_mirror": c["binder_residue_ids"] == rev["residue_ids"],
            "atom_pairs_equal": c["n_atom_pairs"] == rev["n_atom_pairs"],
            "types_mirror": all(
                c["binder_per_residue"][r]["interaction_types"]
                == rev["per_residue"][r]["interaction_types"]
                for r in c["binder_residue_ids"]),
            "buried_area_mirror": all(
                abs(c["binder_per_residue"][r]["buried_area_proxy"]
                    - rev["per_residue"][r]["buried_area_proxy"]) < 1e-6
                for r in c["binder_residue_ids"]),
        }
        rec["symmetry"] = sym
        check(f"{pdb_id}: swapping target and binder mirrors the interface exactly",
              all(sym.values()), json.dumps(sym))

        # 9. self-consistency
        cons = {
            "ids_are_sorted_keys": c["residue_ids"] == sorted(c["per_residue"]),
            "all_within_cutoff": all(v["min_distance"] <= CUTOFF
                                     for v in c["per_residue"].values()),
            "atom_counts_sum": (sum(v["n_atom_contacts"]
                                    for v in c["per_residue"].values())
                                == c["n_atom_pairs"]
                                == sum(v["n_atom_contacts"]
                                       for v in c["binder_per_residue"].values())),
            "no_binder_chain_on_target_side": all(
                v["chain_id"] != bchain for v in c["per_residue"].values()),
            "buried_area_non_negative": min(
                [v["buried_area_proxy"] for v in c["per_residue"].values()],
                default=0.0) >= -0.01,
        }
        rec["self_consistency"] = cons
        check(f"{pdb_id}: contact map is internally consistent", all(cons.values()),
              json.dumps(cons))

        # 10. cutoff monotonicity
        tight = I.chain_contacts(st, tchain, bchain, cutoff=4.0,
                                 with_buried_area=False)
        loose = I.chain_contacts(st, tchain, bchain, cutoff=5.0,
                                 with_buried_area=False)
        mono = (set(tight["residue_ids"]) <= set(c["residue_ids"])
                <= set(loose["residue_ids"]))
        rec["cutoff_series"] = {"4.0": len(tight["residue_ids"]),
                                "4.5": len(c["residue_ids"]),
                                "5.0": len(loose["residue_ids"])}
        check(f"{pdb_id}: contacts grow monotonically with the cutoff",
              mono, json.dumps(rec["cutoff_series"]))

        # 11. no uncomputed area reported as a measured one
        none_area = all(v["buried_area_proxy"] is None
                        for v in c_nosasa["per_residue"].values())
        check(f"{pdb_id}: with_buried_area=False reports None, not 0.0",
              none_area and c_nosasa["residue_ids"] == c["residue_ids"]
              and c_nosasa["buried_area_computed"] is False,
              f"same {len(c['residue_ids'])} residues, dSASA costs "
              f"{round(dt - dt_nosasa, 3)}s of the {dt}s")

        report["cases"][pdb_id] = rec

    # -- guard rails: wrong input must raise, never score zero -------------
    print("\n=== guard rails ===")
    st_y = structures["1YCR"]
    guards = {}
    for label, fn, want in [
        ("unknown binder chain", lambda: I.chain_contacts(st_y, "A", "Z"), KeyError),
        ("unknown target chain", lambda: I.chain_contacts(st_y, "Z", "B"), KeyError),
        ("target chain given as binder", lambda: I.chain_contacts(st_y, "A", "A"),
         ValueError),
    ]:
        try:
            fn()
            guards[label] = "NO ERROR"
        except want as exc:
            guards[label] = f"{type(exc).__name__}: {str(exc)[:60]}"
        except Exception as exc:  # noqa: BLE001
            guards[label] = f"WRONG TYPE {type(exc).__name__}: {str(exc)[:60]}"
    # a HET-only chain (glycans deposited in their own chain) is not a polymer
    # binder: it must be refused, not silently scored as an empty interface
    st_g = I.load_structure(I.fetch_structure(HET_ONLY[0]), structure_id=HET_ONLY[0])
    het_only_chains = [c["chain_id"] for c in I.polymer_chain_summary(st_g)
                       if c["n_residues"] == 0]
    try:
        I.chain_contacts(st_g, HET_ONLY[1], HET_ONLY[2])
        guards["HET-only chain as binder"] = "NO ERROR"
    except ValueError as exc:
        guards["HET-only chain as binder"] = f"ValueError: {str(exc)[:70]}"
    except Exception as exc:  # noqa: BLE001
        guards["HET-only chain as binder"] = f"WRONG TYPE {type(exc).__name__}"
    report["guard_rails"] = guards
    check("wrong input raises instead of returning an empty interface",
          all(not v.startswith(("NO ERROR", "WRONG TYPE")) for v in guards.values()),
          json.dumps(guards))
    check(f"{HET_ONLY[0]}: chain {HET_ONLY[2]} really is HET-only "
          "(absence of polymer residues verified, not assumed)",
          HET_ONLY[2] in het_only_chains,
          f"chains with zero polymer residues: {het_only_chains}")
    default_targets = I.chain_contacts(st_y, None, "B",
                                       with_buried_area=False)["target_chains"]
    check("target_chain=None means every other polymer chain "
          "(the ligand_contacts convention)",
          default_targets == ["A"], f"resolved to {default_targets}")

    # -- 4. numbering: 1YCR chain A must BE UniProt Q00987 numbering -------
    print("\n=== numbering ===")
    seq = I.fetch_uniprot_sequence(MDM2_UNIPROT)
    numbering = {"1YCR": I.check_author_numbering(structures["1YCR"], "A", seq)}
    check("1YCR: author numbering of chain A is UniProt Q00987 numbering",
          numbering["1YCR"]["fraction_match"] == 1.0,
          f"n_compared={numbering['1YCR']['n_compared']}, "
          f"fraction={numbering['1YCR']['fraction_match']}")

    # -- 12. REGRESSION against the committed small-molecule maps ----------
    print("\n=== small-molecule path unchanged (regression) ===")
    if not GOLDEN.exists():
        NOT_EVALUATED.append({
            "item": "small-molecule regression",
            "why": f"{GOLDEN} absent; run scripts/test_interfaces.py first",
        })
    else:
        golden = json.loads(GOLDEN.read_text())
        diffs = {}
        for pid, entry in golden["entries"].items():
            g = entry.get("contacts")
            if not g:
                NOT_EVALUATED.append({"item": f"regression {pid}",
                                      "why": "no contact map in the committed artifact"})
                continue
            info = g["ligand_info"]
            t0 = time.perf_counter()
            st = I.load_structure(I.fetch_structure(pid), structure_id=pid.upper())
            now = I.ligand_contacts(st, info["resname"], chain=info["chain"],
                                    resseq=info["resseq"], cutoff=CUTOFF)
            report["timings_seconds"][f"{pid}_ligand_contacts"] = round(
                time.perf_counter() - t0, 3)
            d = []
            if now["residue_ids"] != g["residue_ids"]:
                d.append(f"residue_ids differ: "
                         f"+{sorted(set(now['residue_ids']) - set(g['residue_ids']))} "
                         f"-{sorted(set(g['residue_ids']) - set(now['residue_ids']))}")
            for rid, want in g["per_residue"].items():
                have = now["per_residue"].get(int(rid))
                if have is None:
                    d.append(f"{rid} missing")
                    continue
                for k in ("resname", "chain_id", "min_distance", "n_atom_contacts",
                          "buried_area_proxy", "interaction_types"):
                    if have[k] != want[k]:
                        d.append(f"{rid}.{k}: {want[k]} -> {have[k]}")
            if now["ligand_info"] != info:
                d.append(f"ligand_info: {info} -> {now['ligand_info']}")
            if now["warnings"] != g["warnings"]:
                d.append(f"warnings: {g['warnings']} -> {now['warnings']}")
            diffs[pid] = d
        report["small_molecule_regression"] = {
            "source": str(GOLDEN.relative_to(ROOT)),
            "generated_utc": golden.get("generated_utc"),
            "entries_compared": sorted(diffs),
            "differences": {k: v for k, v in diffs.items() if v},
        }
        check(f"ligand_contacts reproduces all {len(diffs)} committed contact maps "
              "exactly (ids, distances, atom counts, dSASA, types, warnings)",
              all(not v for v in diffs.values()),
              json.dumps({k: v for k, v in diffs.items() if v})[:400] or
              f"n={len(diffs)} entries: {sorted(diffs)}")

    # -- 13. cross-modality: same target, same numbering, two modalities ----
    print("\n=== cross-modality (MDM2 p53 cleft: peptide vs small molecule) ===")
    pid, comp, lch, lseq = SM_ENTRY
    st_sm, rec_sm = load(pid)
    sm_num = I.check_author_numbering(st_sm, lch, seq)
    lig = [r for r in I.drug_like_ligands(st_sm) if r["resname"] == comp]
    ok_lig = check(f"{pid}: holds the small-molecule ligand {comp} (verified)",
                   bool(lig), f"drug-like ligands: {[r['resname'] for r in lig]}")
    if ok_lig:
        t0 = time.perf_counter()
        sm = I.ligand_contacts(st_sm, comp, chain=lch, resseq=lseq, cutoff=CUTOFF)
        report["timings_seconds"][f"{pid}_ligand_contacts"] = round(
            time.perf_counter() - t0, 3)
        pep = I.chain_contacts(structures["1YCR"], "A", "B", cutoff=CUTOFF)
        a, b = set(sm["residue_ids"]), set(pep["residue_ids"])
        sig_pep = I.consensus_signature([pep])
        sig_sm = I.consensus_signature([sm])
        wj = I.weighted_jaccard(sig_pep, sig_sm)
        cov = I.coverage(sig_pep, sm)
        mixed = I.consensus_signature([pep, sm])
        report["cross_modality"] = {
            "target": "MDM2 (UniProt Q00987), author numbering",
            "numbering_check_1T4E_chainA": sm_num,
            "peptide": {"pdb_id": "1YCR", "chain": "B",
                        "sequence": pep["binder_info"]["sequence"],
                        "n_target_residues": len(b), "residue_ids": sorted(b)},
            "small_molecule": {"pdb_id": pid, "comp_id": comp,
                               "n_heavy_atoms": sm["ligand_info"]["n_heavy_atoms"],
                               "n_target_residues": len(a), "residue_ids": sorted(a)},
            "shared_residues": sorted(a & b),
            "peptide_only": sorted(b - a),
            "small_molecule_only": sorted(a - b),
            "jaccard_unweighted": round(len(a & b) / len(a | b), 4),
            "weighted_jaccard": round(wj, 4),
            "coverage_of_peptide_epitope_by_small_molecule": cov,
            "mixed_consensus": {
                "n_inputs": mixed["n_inputs"],
                "n_core_residues": len(mixed["core_residue_ids"]),
                "core_residue_ids": mixed["core_residue_ids"],
                "modalities": [m.get("binder_kind", "het_ligand")
                               for m in mixed["members"]],
                "warnings": mixed["warnings"],
            },
        }
        check("cross-modality: the two maps are in one numbering and overlap",
              (a & b) and not mixed["warnings"],
              f"n_shared={len(a & b)} of {len(a | b)} union, "
              f"jaccard={round(len(a & b) / len(a | b), 3)}, "
              f"weighted_jaccard={round(wj, 3)}, "
              f"core_coverage of the epitope by {comp}="
              f"{round(cov['core_coverage'], 3)} (n_core={cov['n_core']})")
        check("cross-modality: consensus_signature accepts a polymer contact map "
              "alongside a HET one",
              mixed["n_inputs"] == 2 and len(mixed["core_residue_ids"]) == len(a & b),
              f"n_core={len(mixed['core_residue_ids'])} at frequency 1.0 "
              f"= the {len(a & b)} shared residues")
        check("cross-modality: the peptide epitope is the LARGER footprint "
              "(measured, the flat-epitope expectation)",
              len(b) > len(a), f"peptide n={len(b)} vs small molecule n={len(a)} "
                               f"target residues")

    # -- write ------------------------------------------------------------
    report["numbering_check"] = numbering
    report["not_evaluated"] = NOT_EVALUATED + [
        {"item": "SIFTS-mapped numbering",
         "why": "no SIFTS mapping is applied anywhere in interfaces.py; ids are "
                "author numbers and 1YCR chain A is verified to equal Q00987 "
                "numbering (check 4). Other entries are not."},
        {"item": "binder-side CCD chemistry",
         "why": "a polymer binder's charges and rings come from the amino-acid "
                "tables, so free termini are uncharged and non-standard residues "
                "contribute no salt_bridge/aromatic term. chain_contacts warns."},
        {"item": "sm_addressable verdict for a flat epitope (PROJECT_GOAL 1.4b)",
         "why": "out of scope for the contact extractor: it needs P2Rank pockets "
                "scored against these epitopes, which is a different script."},
    ]
    report["checks"] = CHECKS
    report["n_checks"] = len(CHECKS)
    report["n_failed"] = sum(1 for c in CHECKS if not c["passed"])
    report["environment"] = {
        "python": ".".join(str(v) for v in sys.version_info[:3]),
    }
    try:
        import Bio, numpy  # noqa: N813
        report["environment"]["biopython"] = Bio.__version__
        report["environment"]["numpy"] = numpy.__version__
    except Exception:  # noqa: BLE001
        pass
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, default=str))

    print(f"\n{report['n_checks'] - report['n_failed']}/{report['n_checks']} checks passed")
    print(f"wrote {OUT.relative_to(ROOT)}")
    for c in CHECKS:
        if not c["passed"]:
            print(f"  FAILED: {c['name']} - {c['detail']}")
    return 1 if report["n_failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
