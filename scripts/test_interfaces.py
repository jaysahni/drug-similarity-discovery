"""Self-test for scripts/interfaces.py against real crystal structures.

No synthetic geometry anywhere: every number here comes from a PDB entry that is
downloaded, verified to exist, and parsed. Three VEGFR2/KDR holo complexes carry
the main test (two sorafenib, one axitinib); two further entries act as positive
controls for the two interaction types the VEGFR2 ligands cannot exercise
(salt bridge, halogen bond), because a criterion that never fires is untested.

Checks, in order:
  1. each PDB id exists in RCSB and holds a drug-like ligand (verified, not assumed)
  2. author residue numbering agrees with UniProt P35968, so ids are comparable
  3. the hinge Cys919 is among the extracted contacts of every type-II complex
  4. consensus_signature of one complex with itself -> frequency 1.0 everywhere
  5. consensus_signature of two disjoint contact sets -> frequency 0.5 everywhere
  6. weighted_jaccard(sig, sig) == 1.0 and disjoint signatures give 0.0
  7. coverage() is 1.0 against its own residues and 0.0 against none
  8. every interaction type fires on at least one real structure
  8b. a multi-chain site reports its residue-id collisions instead of merging silently
  8c. primed atom names (C5') survive the structure <-> CCD round trip
  8d. the PDB-format parser gives the same contacts as the mmCIF parser
  9. leave-one-out: consensus of two ligands vs the held-out third's contacts

Usage:
    ./env/bin/python scripts/test_interfaces.py

Writes results/interfaces_selftest.json. Exits non-zero if any check fails.
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
OUT = ROOT / "results" / "interfaces_selftest.json"

# VEGFR2 / KDR kinase domain. UniProt P35968. Author numbering in these entries
# is expected to be UniProt numbering; check 2 verifies that rather than trusting it.
VEGFR2 = ["4ASD", "3WZE", "4AG8"]
VEGFR2_UNIPROT = "P35968"

# Residues a type-II VEGFR2 inhibitor is expected to touch. Cys919 is the hinge;
# Glu885 is the alphaC glutamate and Asp1046/Phe1047 the DFG motif, both of which
# only a DFG-out (type-II) binder reaches. Reported per structure either way.
HINGE = 919
TYPE_II_MARKERS = {885: "GLU", 1046: "ASP", 1047: "PHE", 919: "CYS"}

# Positive controls for the two criteria the VEGFR2 ligands cannot trigger.
# 1N6A: S-adenosylmethionine, whose CCD carboxylate carries a formal -1.
# 2H79: triiodothyronine in the thyroid hormone receptor; iodine halogen bonds.
CONTROLS = [("1N6A", "SAM", "salt_bridge"), ("2H79", "T3", "halogen")]

# A homodimer whose ligand sits on the two-fold axis, so the same author residue
# number occurs in both chains. Author numbering cannot separate those residues,
# and the toolkit has to say so rather than merge them silently.
MULTICHAIN = ("1HXW", "RIT")

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "passed": bool(ok), "detail": str(detail)})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))
    return ok


def rcsb_entry(pdb_id, cache_dir=I.RAW / "rcsb"):
    """Entry metadata from the RCSB data API, cached. Proves the id exists."""
    dest = Path(cache_dir) / f"{pdb_id.upper()}.json"
    if not dest.exists():
        I._download(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}", dest)
    return json.loads(dest.read_text())


def load_entry(pdb_id, include_cofactors=False):
    """Fetch + parse one entry, timed. Returns a record; 'error' set if it failed."""
    rec = {"pdb_id": pdb_id.upper()}
    try:
        meta = rcsb_entry(pdb_id)
        rec["title"] = meta.get("struct", {}).get("title")
        rec["resolution_a"] = (meta.get("rcsb_entry_info", {})
                               .get("resolution_combined", [None])[0])
        rec["method"] = meta.get("rcsb_entry_info", {}).get("experimental_method")
    except Exception as exc:  # noqa: BLE001 - network or unknown id
        rec["error"] = f"RCSB metadata fetch failed: {type(exc).__name__}: {exc}"
        return rec
    try:
        path = I.fetch_structure(pdb_id)
        t0 = time.perf_counter()
        rec["structure"] = I.load_structure(path, structure_id=pdb_id.upper())
        rec["parse_seconds"] = round(time.perf_counter() - t0, 3)
    except Exception as exc:  # noqa: BLE001 - download or parse
        rec["error"] = f"parse failed: {type(exc).__name__}: {exc}"
        return rec

    st = rec["structure"]
    rec["n_polymer_residues"] = sum(
        1 for c in st[0] for r in c if I.is_polymer_residue(r)
    )
    inv = I.het_inventory(st, include_cofactors=include_cofactors)
    rec["n_het_residues"] = len(inv)
    excluded = {}
    for r in inv:
        if not r["kept"]:
            excluded.setdefault(r["reason"], {}).setdefault(r["resname"], 0)
            excluded[r["reason"]][r["resname"]] += 1
    rec["excluded"] = excluded
    rec["ligands"] = [
        {k: v for k, v in r.items() if k != "kept"}
        for r in inv if r["kept"]
    ]
    return rec


def extract(rec, resname, with_buried_area=True):
    """Timed contact extraction for one ligand of a loaded entry."""
    lig = next(l for l in rec["ligands"] if l["resname"] == resname)
    t0 = time.perf_counter()
    c = I.ligand_contacts(rec["structure"], resname, chain=lig["chain_id"],
                          resseq=lig["resseq"], with_buried_area=with_buried_area)
    return c, round(time.perf_counter() - t0, 3)


def main():
    OUT.parent.mkdir(exist_ok=True)
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "environment": {
            "python": sys.version.split()[0],
            "biopython": __import__("Bio").__version__,
            "rdkit": __import__("rdkit").__version__,
            "numpy": __import__("numpy").__version__,
        },
        "numbering": "author (auth_seq_id); no SIFTS mapping applied",
        "contact_cutoff_a": 4.5,
        "core_threshold": 0.6,
    }

    # --- 1. entries exist and hold a drug-like ligand -----------------------
    print("\n1. VEGFR2 entries: existence and drug-like ligand")
    entries, contacts, failures = {}, {}, []
    for pdb in VEGFR2:
        rec = load_entry(pdb)
        entries[pdb] = rec
        if "error" in rec:
            failures.append({"pdb_id": pdb, "reason": rec["error"]})
            check(f"{pdb} loads", False, rec["error"])
            continue
        ok = check(f"{pdb} has >=1 drug-like ligand", bool(rec["ligands"]),
                   f"{[l['resname'] for l in rec['ligands']]} "
                   f"from {rec['n_het_residues']} het residues; "
                   f"{rec['resolution_a']} A")
        if ok:
            c, secs = extract(rec, rec["ligands"][0]["resname"])
            contacts[pdb] = c
            rec["contact_seconds"] = secs
            rec["contact_seconds_no_sasa"] = extract(
                rec, rec["ligands"][0]["resname"], with_buried_area=False)[1]
            print(f"       {c['ligand_info']['resname']} "
                  f"({c['ligand_info']['formula']}) -> "
                  f"{len(c['residue_ids'])} residues in {secs}s")
            print(f"       {c['residue_ids']}")

    if len(contacts) < 2:
        print("\nABORT: fewer than 2 VEGFR2 complexes usable", file=sys.stderr)
        report["aborted"] = "fewer than 2 VEGFR2 complexes usable"
        report["entries"] = {p: {k: v for k, v in r.items() if k != "structure"}
                             for p, r in entries.items()}
        report["failures"] = failures
        report["checks"] = CHECKS
        report["n_checks"] = len(CHECKS)
        report["n_failed"] = sum(1 for c in CHECKS if not c["passed"])
        OUT.write_text(json.dumps(report, indent=2))
        print(f"wrote {OUT}")
        return 1

    # --- 2. numbering is comparable across structures ----------------------
    print("\n2. author numbering vs UniProt", VEGFR2_UNIPROT)
    seq = I.fetch_uniprot_sequence(VEGFR2_UNIPROT)
    report["uniprot"] = {"accession": VEGFR2_UNIPROT, "length": len(seq)}
    numbering = {}
    for pdb, rec in entries.items():
        if "error" in rec or pdb not in contacts:
            continue
        engaged = contacts[pdb]["chains_engaged"]
        if not engaged:
            check(f"{pdb} ligand contacts at least one polymer chain", False,
                  "no target residue within the cutoff")
            continue
        chain = engaged[0]
        res = I.check_author_numbering(rec["structure"], chain, seq)
        numbering[pdb] = res
        check(f"{pdb} chain {chain} author numbering == UniProt numbering",
              res["fraction_match"] > 0.99,
              f"{res['n_match']}/{res['n_compared']} residues match "
              f"({res['fraction_match']:.4f}); mismatches {res['example_mismatches']}")
    report["numbering_check"] = numbering

    # --- 3. hinge recovery -------------------------------------------------
    print(f"\n3. hinge Cys{HINGE} and type-II markers in the extracted contacts")
    markers = {}
    for pdb, c in contacts.items():
        per = {int(k): v for k, v in c["per_residue"].items()}
        present = {rid: (per[rid]["resname"] if rid in per else None)
                   for rid in sorted(TYPE_II_MARKERS)}
        markers[pdb] = {
            "ligand": c["ligand_info"]["resname"],
            "engaged": {str(k): v for k, v in present.items()},
            "hinge_types": per.get(HINGE, {}).get("interaction_types", []),
            "hinge_min_distance": per.get(HINGE, {}).get("min_distance"),
        }
        check(f"{pdb} contacts include hinge Cys{HINGE}",
              present[HINGE] == TYPE_II_MARKERS[HINGE],
              f"residue {HINGE} = {present[HINGE]}, "
              f"min_distance {per.get(HINGE, {}).get('min_distance')} A, "
              f"types {per.get(HINGE, {}).get('interaction_types')}")
        missing = [f"{r}{TYPE_II_MARKERS[r]}" for r in TYPE_II_MARKERS
                   if present[r] != TYPE_II_MARKERS[r]]
        check(f"{pdb} reaches the DFG-out back pocket (type-II markers)",
              not missing, f"missing: {missing}" if missing else "all of "
              f"{[f'{r}{n}' for r, n in sorted(TYPE_II_MARKERS.items())]} engaged")
    report["type_ii_markers"] = markers

    # --- 4/5. consensus frequencies ----------------------------------------
    print("\n4. consensus_signature(same complex twice) -> frequency 1.0")
    first_id = next(iter(contacts))
    first = contacts[first_id]
    dup = I.consensus_signature([first, first])
    freqs = [r["frequency"] for r in dup["residues"]]
    # (4) min()/max() are computed only when there is something to compute, so an
    # empty signature records a FAIL instead of raising out of the run
    check(f"duplicate consensus ({first_id}): every frequency == 1.0",
          bool(freqs) and all(f == 1.0 for f in freqs),
          f"n_residues={len(freqs)}"
          + (f", min={min(freqs)}, max={max(freqs)}" if freqs else ""))
    check("duplicate consensus: every residue is core",
          len(dup["core_residue_ids"]) == len(dup["residues"]),
          f"{len(dup['core_residue_ids'])}/{len(dup['residues'])}")

    print("\n5. consensus_signature(two disjoint contact sets) -> frequency 0.5")
    ctrl = {}
    for pdb, resname, _ in CONTROLS:
        rec = load_entry(pdb, include_cofactors=True)
        entries[pdb] = rec
        if "error" in rec or not any(l["resname"] == resname for l in rec["ligands"]):
            failures.append({"pdb_id": pdb,
                             "reason": rec.get("error", f"no {resname} ligand kept")})
            check(f"{pdb} control loads with {resname}", False,
                  rec.get("error", "ligand absent"))
            continue
        c, secs = extract(rec, resname)
        rec["contact_seconds"] = secs
        ctrl[pdb] = c
        check(f"{pdb} control loads with {resname}", True,
              f"{len(c['residue_ids'])} contact residues in {secs}s")

    disjoint_pair = None
    for pdb, c in ctrl.items():
        overlap = set(first["residue_ids"]) & set(c["residue_ids"])
        if not overlap:  # verified empty, not assumed
            disjoint_pair = (first_id, pdb, c)
            break
    if disjoint_pair is None:
        check("a genuinely disjoint real contact pair exists", False,
              "every control shares residue numbers with 4ASD")
        sig_disjoint = None
    else:
        a_id, b_id, cb = disjoint_pair
        sig_disjoint = I.consensus_signature([first, cb])
        f2 = [r["frequency"] for r in sig_disjoint["residues"]]
        check(f"disjoint consensus ({a_id} + {b_id}): every frequency == 0.5",
              bool(f2) and all(f == 0.5 for f in f2),
              f"n_residues={len(f2)} = {len(first['residue_ids'])}"
              f"+{len(cb['residue_ids'])}"
              + (f", min={min(f2)}, max={max(f2)}" if f2 else ""))
        check("disjoint consensus: no residue is core at threshold 0.6",
              not sig_disjoint["core_residue_ids"],
              f"core={sig_disjoint['core_residue_ids']}")

    # --- 6. weighted_jaccard ------------------------------------------------
    print("\n6. weighted_jaccard")
    sigs = {p: I.consensus_signature([c]) for p, c in contacts.items()}
    consensus = I.consensus_signature(list(contacts.values()))
    check("weighted_jaccard(sig, sig) == 1.0",
          I.weighted_jaccard(consensus, consensus) == 1.0,
          f"{I.weighted_jaccard(consensus, consensus)}")
    if sig_disjoint is not None:
        a_id, b_id, cb = disjoint_pair
        j = I.weighted_jaccard(sigs[a_id], I.consensus_signature([cb]))
        check(f"weighted_jaccard of disjoint signatures ({a_id} vs {b_id}) == 0.0",
              j == 0.0, f"{j}")

    cross = {}
    for i, a in enumerate(contacts):
        for b in list(contacts)[i + 1:]:
            cross[f"{a}_vs_{b}"] = round(I.weighted_jaccard(sigs[a], sigs[b]), 4)
    print(f"       pairwise VEGFR2 overlap: {cross}")
    report["pairwise_weighted_jaccard"] = cross

    # --- 7. coverage --------------------------------------------------------
    print("\n7. coverage")
    self_cov = I.coverage(consensus, consensus["core_residue_ids"])
    no_core = I.coverage({"residues": [{"residue_id": 1, "chain_id": "A",
                                        "residue_name": "ALA", "frequency": 0.2,
                                        "is_core": False}]}, [1])
    check("coverage() of a signature with no core residues is None, not 0.0",
          no_core["core_coverage"] is None, f"{no_core['core_coverage']}")
    check("coverage(sig, its own core residues).core_coverage == 1.0",
          self_cov["core_coverage"] == 1.0, f"{self_cov['core_coverage']}")
    empty_cov = I.coverage(consensus, [])
    check("coverage(sig, []).core_coverage == 0.0",
          empty_cov["core_coverage"] == 0.0, f"{empty_cov['core_coverage']}")

    # --- 8. every interaction type fires somewhere --------------------------
    print("\n8. interaction-type coverage across all real structures")
    fired = {}
    for pdb, c in list(contacts.items()) + list(ctrl.items()):
        for rec_r in c["per_residue"].values():
            for t in rec_r["interaction_types"]:
                fired.setdefault(t, []).append(pdb)
    for t in ("hydrophobic", "hbond", "aromatic", "salt_bridge", "halogen"):
        where = sorted(set(fired.get(t, [])))
        check(f"interaction type '{t}' fires on a real structure",
              bool(where), f"seen in {where}" if where else "never triggered")
    report["interaction_types_observed"] = {
        t: sorted(set(v)) for t, v in sorted(fired.items())
    }

    # VEGFR2-specific negatives, measured rather than assumed
    neg = {
        "vegfr2_ligand_formal_charges": {
            p: c["ligand_info"]["n_charged_atoms"] for p, c in contacts.items()},
        "note": ("no salt bridge is detectable for these ligands because the CCD "
                 "models all three as formally neutral; halogen bonds are absent "
                 "because the nearest sorafenib Cl...O separation exceeds the "
                 "3.5 A cutoff (measured below)"),
    }
    import numpy as np

    for pdb, c in contacts.items():
        info = c["ligand_info"]
        m = entries[pdb]["structure"][0]
        lig = [r for ch in m for r in ch
               if r.get_resname().strip().upper() == info["resname"]
               and ch.id == info["chain"] and r.id[1] == info["resseq"]]
        hal = [a for a in I.heavy_atoms(lig[0]) if a.element.upper() in I.HALOGENS] if lig else []
        tgt = [a for ch in m for r in ch if I.is_polymer_residue(r)
               for a in I.heavy_atoms(r) if a.element.upper() in ("N", "O")]
        key = f"{pdb}_{info['resname']}_min_halogen_to_NO_a"
        if not hal:
            neg[key] = "no Cl/Br/I in this ligand"
        elif not tgt:
            neg[key] = "no N/O atoms in the polymer"
        else:
            neg[key] = round(min(float(np.linalg.norm(a.coord - t.coord))
                                 for a in hal for t in tgt), 3)
    report["measured_negatives"] = neg
    print(f"       {neg}")

    # does the 4.5 A heavy-atom pre-filter truncate the 5.5 A aromatic centroid
    # criterion on any of these structures? measured, not assumed
    trunc = {}
    for pdb, c in list(contacts.items()) + list(ctrl.items()):
        info = c["ligand_info"]
        st = entries[pdb]["structure"]
        cand = I.aromatic_centroid_candidates(st, info["resname"],
                                              chain=info["chain"],
                                              resseq=info["resseq"])
        trunc[pdb] = sorted(set(cand) - set(c["residue_ids"]))
    report["aromatic_centroid_truncation"] = trunc
    check("no aromatic ring is within 5.5 A by centroid yet outside the 4.5 A "
          "heavy-atom cutoff", all(not v for v in trunc.values()),
          f"missed per structure: {trunc}")

    # --- 8b. multi-chain site ----------------------------------------------
    print("\n8b. multi-chain control: residue-id collisions are reported")
    mc_pdb, mc_lig = MULTICHAIN
    mc_rec = load_entry(mc_pdb)
    entries[mc_pdb] = mc_rec
    mc = None
    if "error" in mc_rec or not any(l["resname"] == mc_lig for l in mc_rec["ligands"]):
        failures.append({"pdb_id": mc_pdb,
                         "reason": mc_rec.get("error", f"no {mc_lig} ligand kept")})
        check(f"{mc_pdb} multi-chain control loads", False,
              mc_rec.get("error", "ligand absent"))
    else:
        mc, mc_secs = extract(mc_rec, mc_lig)
        mc_rec["contact_seconds"] = mc_secs
        collide = [w for w in mc["warnings"] if "collision" in w]
        check(f"{mc_pdb} contacts span >1 polymer chain",
              len(mc["chains_engaged"]) > 1, f"chains {mc['chains_engaged']}, "
              f"{len(mc['residue_ids'])} residue ids in {mc_secs}s")
        check(f"{mc_pdb} reports the residue-id collision instead of merging silently",
              bool(collide), collide[0] if collide else "no collision warning emitted")
        report["multichain_control"] = {
            "pdb_id": mc_pdb, "ligand": mc_lig,
            "chains_engaged": mc["chains_engaged"],
            "warnings": mc["warnings"],
        }

    # all VEGFR2 structures are single-chain: verified here, not assumed, because
    # the leave-one-out comparison below would be meaningless if they were not
    single = {p: c["chains_engaged"] for p, c in contacts.items()}
    check("every VEGFR2 test structure engages exactly one polymer chain",
          all(len(v) == 1 for v in single.values()), f"{single}")

    # --- 8c. atom-name round trip ------------------------------------------
    # mmCIF quotes names containing a prime. Stripping quote characters blindly
    # turns C5' into C5, which is a different atom; the CCD lookup then misses
    # and charge/ring typing silently degrades. SAM carries eight primed heavy
    # atoms, so it is the regression test for that.
    print("\n8c. primed atom names survive the structure <-> CCD round trip")
    if "1N6A" in ctrl:
        mol = I.ccd_mol("SAM")
        ccd_names = {a.GetProp("pdb_name") for a in mol.GetAtoms()} if mol else set()
        lig = [r for c in entries["1N6A"]["structure"][0] for r in c
               if r.get_resname().strip().upper() == "SAM"][0]
        obs = [I._clean_name(a.get_name()) for a in I.heavy_atoms(lig)]
        primed = [n for n in obs if "'" in n]
        unmapped = [n for n in obs if n not in ccd_names]
        check("SAM has primed heavy-atom names to test with", bool(primed),
              f"{sorted(primed)}")
        check("every SAM heavy atom in the structure maps to a CCD atom",
              not unmapped, f"unmapped: {unmapped}" if unmapped
              else f"{len(obs)}/{len(obs)} mapped")

    # --- 8d. format agnosticism --------------------------------------------
    # load_structure claims to read both formats; that claim gets checked rather
    # than trusted, because author numbering is the only thing tying them together
    print("\n8d. PDB format vs mmCIF for the same entry")
    try:
        pdb_path = I.fetch_structure(first_id, fmt="pdb")
        st_pdb = I.load_structure(pdb_path, structure_id=first_id)
        lig0 = first["ligand_info"]
        c_pdb = I.ligand_contacts(st_pdb, lig0["resname"], chain=lig0["chain"],
                                  resseq=lig0["resseq"], with_buried_area=False)
        same = c_pdb["residue_ids"] == first["residue_ids"]
        check(f"{first_id}: .pdb and .cif give identical contact residue ids", same,
              f"{len(c_pdb['residue_ids'])} vs {len(first['residue_ids'])} residues"
              + ("" if same else
                 f"; cif-only {sorted(set(first['residue_ids']) - set(c_pdb['residue_ids']))}, "
                 f"pdb-only {sorted(set(c_pdb['residue_ids']) - set(first['residue_ids']))}"))
        report["format_check"] = {
            "pdb_id": first_id,
            "n_residues_cif": len(first["residue_ids"]),
            "n_residues_pdb": len(c_pdb["residue_ids"]),
            "identical": same,
        }
    except Exception as exc:  # noqa: BLE001 - legacy PDB format may not exist
        check(f"{first_id}: .pdb and .cif give identical contact residue ids", False,
              f"{type(exc).__name__}: {exc}")

    # the consensus is only meaningful if the ids mean the same thing everywhere
    check("VEGFR2 consensus reports no residue-id disagreement",
          not consensus.get("warnings"), f"{consensus.get('warnings')}")

    # --- 9. leave-one-out: the central experiment, at n=3 -------------------
    print("\n9. leave-one-out consensus -> held-out ligand coverage (n=3)")
    loo = {}
    for held in contacts:
        rest = [c for p, c in contacts.items() if p != held]
        sig = I.consensus_signature(rest)
        cov = I.coverage(sig, contacts[held])
        loo[held] = {
            "held_out_ligand": contacts[held]["ligand_info"]["resname"],
            "trained_on": [p for p in contacts if p != held],
            "n_core_in_consensus": cov["n_core"],
            "core_coverage": (round(cov["core_coverage"], 4)
                              if cov["core_coverage"] is not None else None),
            "weighted_coverage": round(cov["weighted_coverage"], 4),
            "missed_core_residues": cov["missed_core_residues"],
        }
        print(f"       hold out {held} ({loo[held]['held_out_ligand']}): "
              f"core_coverage {loo[held]['core_coverage']} over "
              f"{cov['n_core']} core residues, missed {cov['missed_core_residues']}")
    report["leave_one_out"] = loo

    # --- assemble ----------------------------------------------------------
    report["entries"] = {}
    for pdb, rec in entries.items():
        out = {k: v for k, v in rec.items() if k not in ("structure",)}
        if pdb in contacts or pdb in ctrl or (mc is not None and pdb == mc_pdb):
            c = contacts.get(pdb) or ctrl.get(pdb) or mc
            out["contacts"] = {
                "ligand_info": c["ligand_info"],
                "n_contact_residues": len(c["residue_ids"]),
                "residue_ids": c["residue_ids"],
                "per_residue": {str(k): v for k, v in c["per_residue"].items()},
                "warnings": c["warnings"],
            }
        report["entries"][pdb] = out

    report["timings_seconds"] = {
        pdb: {"parse": rec.get("parse_seconds"),
              "contacts_with_sasa": rec.get("contact_seconds"),
              "contacts_without_sasa": rec.get("contact_seconds_no_sasa")}
        for pdb, rec in entries.items() if "error" not in rec
    }
    report["consensus_vegfr2"] = consensus
    report["failures"] = failures
    report["not_evaluated"] = [
        {"item": "source='boltzgen_consensus' signature",
         "why": "BoltzGen needs a GPU or a Rowan API key; neither exists on this "
                "machine. The same consensus_signature code path is exercised "
                "with source='known_ligand' from crystal contacts instead."},
        {"item": "Boltz-2 co-folding of approved molecules into the site",
         "why": "no GPU, no Rowan API key"},
        {"item": "SIFTS residue mapping to UniProt numbering",
         "why": "not applied. Author numbering is used throughout and verified "
                "against the UniProt P35968 sequence in check 2; a structure "
                "whose author numbering is not UniProt numbering would fail that "
                "check rather than silently produce incomparable ids."},
        {"item": "validation of interaction typing against PLIP",
         "why": "PLIP is not installed and no ground-truth contact-type set was "
                "available offline. Types are distance-based geometric "
                "compatibility only; no angular term is applied except the "
                "ring-centroid criterion, and hydrogens are absent from these "
                "crystal models so donor/acceptor roles are not assigned."},
        {"item": "salt bridges through groups the CCD models as neutral",
         "why": "ligand formal charges come from the deposited CCD definition and "
                "there is no pKa model available offline. Verified: BAX, AXI, T3, "
                "RIT and ATP all return zero charged atoms, so a carboxylic acid "
                "or amine that is charged at pH 7.4 is invisible to the "
                "salt_bridge criterion. The criterion itself is exercised on "
                "1N6A, where SAM's CCD carboxylate does carry a formal -1."},
        {"item": "cation-pi and pi-stacking geometry split",
         "why": "not implemented; 'aromatic' is a single centroid-distance type"},
        {"item": "sm_addressable subsignature (PROJECT_GOAL.md 1.4b)",
         "why": "needs pocket geometry from P2Rank, which is a different task's "
                "scope; the signature's source discriminator is in place for it"},
    ]

    n_fail = sum(1 for c in CHECKS if not c["passed"])
    report["checks"] = CHECKS
    report["n_checks"] = len(CHECKS)
    report["n_failed"] = n_fail

    OUT.write_text(json.dumps(report, indent=2))
    print(f"\n{len(CHECKS) - n_fail}/{len(CHECKS)} checks passed")
    if failures:
        print(f"structures that failed to load: {failures}")
    print(f"wrote {OUT}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
