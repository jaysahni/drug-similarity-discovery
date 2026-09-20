"""Leave-one-ligand-out hotspot recovery: does a consensus interface signature
carry information beyond the pocket P2Rank finds in under a second?

READ THIS FIRST -- THE DECOMPOSITION, NOT THE HEADLINE
    The headline number of this experiment is the consensus's margin over
    P2Rank. On its own it is misleading, because a null that costs nothing --
    take ONE other ligand of the same target, at random, and predict the
    held-out ligand's contacts with that single ligand's contacts -- already
    beats P2Rank by most of that margin. Most of what ANY arm here "predicts"
    is that the pocket is the pocket, which is exactly what P2Rank returns for
    free. The share of the margin that AGGREGATING many ligands actually buys
    is reported explicitly in `aggregation_decomposition`, and that is the
    number this experiment exists to produce. Quote the decomposition, not the
    headline, and never quote the headline against the random-surface floor.

WHICH ABLATION THIS IS -- READ THIS BEFORE QUOTING THE NUMBER
    This is ablation I6.2 of PROJECT_GOAL.md, NOT I6.1, and it does NOT settle
    the M2 gate. PROJECT_GOAL.md I6 lists them verbatim: I6.1 is "BoltzGen
    signature vs `p2rank_geometry` signature (S8.3). *The* experiment"; I6.2 is
    "Signature vs `known_ligand` signature -- the ceiling". Section 8.3 is
    equally explicit: `source="known_ligand"` GIVES THE UPPER BOUND. Everything
    below is built with source="known_ligand", so it measures the ceiling the
    absent BoltzGen arm would have had to reach, and nothing about BoltzGen.

    The asymmetry that makes this worth running without a GPU: a NEGATIVE result
    here would have foreclosed the BoltzGen arm a fortiori, because no generator
    can beat P2Rank through a ceiling that does not itself beat P2Rank. The
    result is positive, so it forecloses nothing and licenses nothing about
    BoltzGen -- it establishes only that headroom over P2Rank exists.
    The M2 gate is neither passed nor failed by this result.

Section 7 is blunt about the stakes of I6.1: if the consensus adds nothing over
a cheap pocket predictor, the honest pipeline is the cheap one and that is what
should ship. This run cannot answer that question; it only bounds it.

WHAT IS AND IS NOT BEING TESTED
    PROJECT_GOAL.md 4.3/6.I build the consensus from thousands of BoltzGen
    designs. That arm is ABSENT here: no design ensemble exists to aggregate.
    The recorded reason is docs/04-BOLTZGEN-MODAL.md (verified 2026-09-19): the
    Modal account authenticates but every GPU tier tried (T4, L4, A10G, A100)
    was refused pending a payment method. That is a DATED OBSERVATION, not a
    standing fact -- this script never probes Modal, so if the account has since
    been given a payment method the sentence is stale and only re-running
    scripts/modal_boltzgen.py can say so. What is load-bearing here and is
    verified by this script's own output is only that no BoltzGen arm was
    scored: `generated_arm.present` is False. PROJECT_GOAL.md 4.3 deliberately
    makes signature `source` a discriminator so the builder is swappable, and
    8.1 says the primary metric is "just geometry against crystal structures",
    so the same experiment is run with the consensus built from REAL
    CO-CRYSTALS: the `source="known_ligand"` builder of task C8. A result here
    is evidence about the *interface-signature representation*, not about
    BoltzGen.

THE TRIAL
    One trial = one held-out ligand of one target.
        O = the held-out ligand's own contact residues (ground truth)
    The arms, as a NULL LADDER from the cheapest floor to the arm of interest:
        rung 1  random                uniform draw from the binding chain's
                                      solvent-exposed residues. The trivial floor
        rung 2  buried                the most BURIED residues of the binding
                                      chain by isolated-chain SASA. A shape-only
                                      null that costs nothing and needs no ligand.
                                      Scored over the whole chain (`buried*`,
                                      which selects the deep core) and over the
                                      solvent-exposed residues only
                                      (`buried_surface_size_matched`, the "any
                                      concave spot" reading), because which one
                                      "most buried" means is itself a choice
        rung 3  null1_one_random_other_ligand
                                      ONE other ligand of the same target, drawn
                                      at random, averaged over --null-resamples
                                      draws. THE KEY NULL: the consensus beats it
                                      only by whatever AGGREGATION buys
        rung 4  null2_most_dissimilar_other_ligand
                                      the other ligand with the LOWEST ECFP4
                                      Tanimoto to the held-out one. The hardest
                                      single-ligand null, and the one that speaks
                                      to the cross-chemistry claim
        A       consensus             core_residue_ids of consensus_signature
                                      over every OTHER distinct ligand
        B       p2rank_top1           rank-1 P2Rank pocket of the held-out
                                      ligand's OWN structure, restricted to the
                                      binding chain. A vs B is the decisive
                                      comparison (PROJECT_GOAL.md 8.3)
        C       p2rank_best_of_3      whichever of the top 3 pockets overlaps O
                                      best. GENEROUS TO THE BASELINE AND NOT A
                                      FAIR ARM: it picks the pocket using the
                                      answer. An upper bound on B, never B
        E       union_all             union of every other ligand's contacts with
                                      no frequency threshold. The ceiling on what
                                      any consensus over these ligands can cover

TWO SET-CONSTRUCTION RULES, BOTH REPORTED
    A consensus is a FREQUENCY MAP over residues. Turning it into a predicted
    set needs a rule, and the rule is a researcher degree of freedom that can
    change the answer. Both rules are therefore computed for every arm where
    they apply and reported side by side in `set_construction_rules`:
        core@threshold   keep residues with frequency >= --core-threshold
                         (PROJECT_GOAL.md 4.3's default, 0.6)
        size_matched     keep the |O| highest-frequency residues, so the
                         prediction is the same size as the truth and the
                         precision/recall trade is removed by construction
    Ties in frequency at the size-matched cut are broken by the mean minimum
    heavy-atom distance across the ligands that engage the residue (tighter
    first), then by residue key. How often a tie straddles the cut is reported
    (`notes.trials_with_size_matched_tie_at_cut`) and the whole arm is re-scored
    with the tie-break removed (`set_construction_rules.size_matched_tie_break`),
    so the tie-break's weight is measured rather than assumed to be small.
    A threshold sweep over core@{0.2..0.8} is reported as well.

RESIDUE KEYS ARE CHAIN- AND INSERTION-CODE-AWARE
    A residue is keyed by (chain_id, auth_seq_id, insertion_code), not by bare
    auth_seq_id. Two things follow, and both are measured here:
      - Contacts outside the binding chain are DROPPED, never merged into it.
        The previous bare-id keying merged a residue id seen in two chains into
        one record and then filtered by that record's chain, which could both
        attribute another chain's contacts to the binding chain and lose a real
        binding-chain residue.
      - Thrombin's 60A..60I no longer collapse onto 60. `insertion_code_audit`
        reports how many trials carry insertion-coded contacts and by how much
        the keying changes them.
    Within a trial every scored set is restricted to ONE chain (the binding
    chain), so the chain component of the key is constant and the key is written
    "<auth_seq_id><insertion_code>" (e.g. 60A). Across entries of one target the
    binding chain's LETTER is a per-deposition label, not an identity -- entry
    1ABC's chain A and entry 2DEF's chain B are the same entity -- so the letter
    is normalised away by that restriction rather than being compared literally;
    `site_chain_letter_audit` reports how often it differs, which is how much a
    literal chain-letter key would have destroyed. Each ligand's own letter is
    kept in the per-trial record.
    P2Rank residues carry a chain but NO insertion code (verified at load:
    `inputs.p2rank_residue_fields`), so an insertion-coded truth residue is
    unreachable for that arm; the count is reported and a sensitivity analysis
    drops the affected targets.

PRIMARY METRIC
    PROJECT_GOAL.md 4.4 makes core coverage primary, so the primary metric here
    is `precision` = |predicted hotspots INTERSECT O| / |predicted hotspots|:
    the fraction of the predicted hotspot set that the held-out ligand actually
    engages. For arm A this is exactly interfaces.coverage(sig, O)["core_coverage"]
    -- asserted numerically at runtime, see `core_coverage_identity_max_deviation`
    in the output. Defining it this way is what makes the number comparable to
    the P2Rank arms, which have no notion of a "core".
    Precision alone rewards predicting few residues, so recall, F1 and Jaccard
    are computed and reported for every arm alongside it, and every arm's
    predicted set size is recorded. If the arms disagree between precision and
    F1, that disagreement is a finding and is stated in the output.

CONVENTIONS THAT AFFECT THE HEADLINE (all recorded in the output)
    empty prediction   precision 0.0, not undefined. A method that predicts no
                       residue in the binding chain earns no credit; scoring it
                       NaN would silently drop its worst cases from the paired
                       test and flatter it. `primary_defined_only` re-runs the
                       headline on trials where both A and B are non-empty.
    binding chain      the polymer chain contributing the most contact residues
                       to that ligand, not the ligand's own chain letter (in a
                       few entries the ligand is deposited under its own chain).
                       Contacts, P2Rank pockets and the null arms are all
                       restricted to that one chain.
    one copy           if a ligand appears in several chains of one entry, the
                       first copy by (chain, resseq) is used.
    one entry          if a ligand appears in several entries of one target, the
                       best-resolution entry is used and the rest are counted as
                       collapsed duplicates. Twenty structures of one ligand are
                       one ligand's worth of evidence.
    numbering          author numbering (auth_seq_id + insertion code), no SIFTS
                       mapping. Whether it lines up across a target's entries is
                       measured per trial (`n_residue_id_disagreements`) and a
                       sensitivity analysis on the clean trials is reported.
    lipids out         membrane lipids, detergents and organomercurials are
                       excluded from BOTH the held-out set and the consensus
                       (interfaces.LIPID_DETERGENT_REASON). The headline is
                       re-reported with them left in, as it was before.
    buried area        NOT computed (with_buried_area=False). Every metric here
                       is set-membership, and delta-SASA costs two
                       Shrake-Rupley passes per ligand. `mean_buried_area` is
                       therefore None throughout, never 0.0.

Usage:
    ./env/bin/python scripts/hotspot_recovery.py                  # full cached set
    ./env/bin/python scripts/hotspot_recovery.py --limit-targets 3 --workers 4
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
import time
import zlib
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import interfaces as I  # noqa: E402
import metrics as M  # noqa: E402

STRUCTURES = ROOT / "data" / "raw" / "structures"
MANIFEST = ROOT / "data" / "raw" / "structures_manifest.json"
P2RANK_JSON = ROOT / "results" / "p2rank_pockets.json"
OUT_JSON = ROOT / "results" / "hotspot_recovery.json"
OUT_CSV = ROOT / "results" / "hotspot_recovery_by_target.csv"

# Intermediate contact extraction is cached OUTSIDE the repo: this script owns
# only its two results/ files. Override with HOTSPOT_CACHE.
_SCRATCH = Path("/private/tmp/claude-501/-Users-evanxiang-Desktop-Projects-"
                "drug-similarity-discovery/ff08da32-3cac-44f3-9083-f8fc1e02ca4d/"
                "scratchpad")
CACHE = Path(os.environ.get("HOTSPOT_CACHE")
             or (_SCRATCH if _SCRATCH.exists() else Path("/tmp")) / "hotspot_cache")
# v3: chain/insertion-code-aware keys, per-residue min distance kept so the
# contact cutoff can be swept without re-parsing, per-residue SASA kept for the
# buried-null arm, lipid-class ligands kept but flagged.
# v4: min distances kept to 6 decimals instead of 3 (see _full_key_contacts).
CACHE_VERSION = "contacts_v4"

# Rules for turning a consensus frequency map into a predicted set.
RULE_CORE = "core@threshold"
RULE_SIZE = "size_matched_top_O"

ARMS = [
    "consensus",                          # rule core@threshold
    "consensus_size_matched",             # rule size_matched_top_O
    "p2rank_top1",
    "p2rank_best_of_3",
    "null1_one_random_other_ligand",      # rung 3 -- the null that matters
    "null2_most_dissimilar_other_ligand",  # rung 4
    "buried",                             # rung 2, |core|-sized, whole chain
    "buried_size_matched",                # rung 2, |O|-sized, whole chain
    "buried_surface_size_matched",        # rung 2, |O|-sized, exposed residues only
    "random",                             # rung 1, |core|-sized
    "random_size_matched",                # rung 1, |O|-sized
    "union_all",
]
# The ladder, weakest first, as it is reported and printed.
LADDER = ["random_size_matched", "buried_size_matched",
          "buried_surface_size_matched",
          "null1_one_random_other_ligand", "null2_most_dissimilar_other_ligand",
          "p2rank_top1"]
METRICS = ["precision", "recall", "f1", "jaccard"]
PRIMARY = "precision"
SURFACE_MIN_SASA = 10.0  # square angstroms, on the isolated binding chain
THRESHOLD_SWEEP = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
CUTOFF_SWEEP = (4.0, 4.5, 5.0)   # contact cutoffs in angstrom
EXTRACT_CUTOFF = max(CUTOFF_SWEEP)


# --------------------------------------------------------------------------
# residue keys
# --------------------------------------------------------------------------
def res_key(residue):
    """Insertion-code-aware key of a residue WITHIN its chain: '60', '60A'."""
    return f"{residue.id[1]}{residue.id[2].strip()}"


def key_seq(key):
    """The bare auth_seq_id of a residue key -- the LEGACY key."""
    i = 0
    while i < len(key) and (key[i].isdigit() or (i == 0 and key[i] == "-")):
        i += 1
    return key[:i]


# --------------------------------------------------------------------------
# stage 1: contacts per structure (cached, one JSON per entry)
# --------------------------------------------------------------------------
def _sasa_by_residue(model, chain_id):
    """SASA of every polymer residue of one chain, computed on that chain ALONE.

    Feeds two arms: the random arm draws from residues at or above
    SURFACE_MIN_SASA (drawing from the buried core would let it sample residues
    no ligand can reach and would understate the floor), and the buried arm
    takes the LOWEST-SASA residues.
    """
    from Bio.PDB.SASA import ShrakeRupley

    sr = ShrakeRupley()
    work = copy.deepcopy(model)
    work.detach_parent()
    for ch in list(work):
        if ch.id != chain_id:
            work.detach_child(ch.id)
            continue
        for res in list(ch):
            if not I.is_polymer_residue(res):
                ch.detach_child(res.id)
                continue
            for atom in list(res):
                if atom.element.upper() not in sr.radii_dict:
                    res.detach_child(atom.id)
            if len(res) == 0:
                ch.detach_child(res.id)
    if chain_id not in [c.id for c in work] or len(work[chain_id]) == 0:
        return {}
    sr.compute(work, level="R")
    return {res_key(r): round(float(r.sasa), 2) for r in work[chain_id]}


def _polymer_search(model):
    from Bio.PDB import NeighborSearch

    atoms = [a for ch in model for r in ch if I.is_polymer_residue(r)
             for a in I.heavy_atoms(r)]
    return (NeighborSearch(atoms) if atoms else None)


def _full_key_contacts(search, lig_res, cutoff):
    """{('CHAIN', 'KEY'): {...}} for every polymer residue within `cutoff`.

    Keyed by chain AND insertion code, so nothing is merged. `min_distance` is
    kept so any cutoff <= `cutoff` can be applied later without re-parsing. It
    is kept to 6 decimals, not 3: at 3 a contact at 4.5004 A rounds to 4.500 and
    passes a `<= 4.5` test that interfaces.ligand_contacts' radius search fails,
    which is exactly the disagreement the crosscheck caught on three ligands.
    """
    out = {}
    for la in I.heavy_atoms(lig_res):
        for ta in search.search(la.coord, cutoff):
            tres = ta.get_parent()
            ch = tres.get_parent().id
            k = f"{ch}|{res_key(tres)}"
            d = float(np.linalg.norm(la.coord - ta.coord))
            rec = out.get(k)
            if rec is None:
                out[k] = {"resname": tres.get_resname().strip().upper(),
                          "chain_id": ch, "min_distance": round(d, 6),
                          "n_atom_contacts": 1}
            else:
                rec["n_atom_contacts"] += 1
                rec["min_distance"] = min(rec["min_distance"], round(d, 6))
    return out


def extract_one(job):
    """One cached structure -> its ligands' target-side contacts, full keys."""
    pdb_id, path, min_heavy, primary_cutoff, crosscheck = job
    # the cache key carries every parameter that changes what is stored
    dest = (CACHE / f"{CACHE_VERSION}_c{primary_cutoff}_h{min_heavy}"
            / f"{pdb_id}.json")
    if dest.exists():
        try:
            return json.loads(dest.read_text())
        except Exception:  # noqa: BLE001 - corrupt cache entry, recompute
            pass

    t0 = time.perf_counter()
    rec = {"pdb_id": pdb_id, "error": None, "ligands": [], "sasa": {},
           "n_het": 0, "n_drug_like_copies": 0, "n_lipid_class_copies": 0}
    try:
        st = I.load_structure(path, structure_id=pdb_id)
        model = st[0]
        inv = I.het_inventory(st, min_heavy_atoms=min_heavy)
        rec["n_het"] = len(inv)
        lipid_tag = "exclusion list: " + I.LIPID_DETERGENT_REASON_PREFIX
        kept = []
        for r in inv:
            if r["kept"]:
                kept.append({**r, "lipid_class": False})
            elif r.get("reason", "").startswith(lipid_tag):
                # kept in the cache but FLAGGED, so the run can report the
                # headline with and without them without re-extracting
                kept.append({**r, "lipid_class": True})
        rec["n_drug_like_copies"] = sum(1 for r in kept if not r["lipid_class"])
        rec["n_lipid_class_copies"] = sum(1 for r in kept if r["lipid_class"])

        search = _polymer_search(model)
        if search is None:
            raise ValueError("structure has no polymer residues to contact")

        seen_comp = set()
        for lig in sorted(kept, key=lambda r: (r["chain_id"], r["resseq"])):
            comp = lig["resname"]
            if comp in seen_comp:          # one copy per component id per entry
                continue
            seen_comp.add(comp)
            try:
                lig_res = st[0][lig["chain_id"]][
                    [r.id for r in st[0][lig["chain_id"]]
                     if r.id[1] == lig["resseq"]
                     and r.get_resname().strip().upper() == comp][0]]
                full = _full_key_contacts(search, lig_res, EXTRACT_CUTOFF)
            except Exception as exc:  # noqa: BLE001
                rec["ligands"].append({"comp_id": comp, "error": repr(exc)})
                continue
            at_primary = {k: v for k, v in full.items()
                          if v["min_distance"] <= primary_cutoff}
            if not at_primary:
                rec["ligands"].append({"comp_id": comp, "error": "no contacts"})
                continue
            by_chain = defaultdict(list)
            for k, v in at_primary.items():
                by_chain[v["chain_id"]].append(k)
            site_chain = max(sorted(by_chain), key=lambda ch: len(by_chain[ch]))

            cross = None
            if crosscheck:
                # Independent check that this module's own contact pass agrees
                # with interfaces.ligand_contacts, whose residue sets are the
                # published definition. Compared on BARE seq ids over ALL
                # chains, because that is what the library returns.
                try:
                    c = I.ligand_contacts(st, comp, chain=lig["chain_id"],
                                          resseq=lig["resseq"],
                                          cutoff=primary_cutoff,
                                          with_buried_area=False)
                    lib = {str(r) for r in c["per_residue"]}
                    mine = {key_seq(k.split("|", 1)[1]) for k in at_primary}
                    cross = {"n_only_library": len(lib - mine),
                             "n_only_here": len(mine - lib),
                             "n_shared": len(lib & mine),
                             "warnings": c["warnings"]}
                except Exception as exc:  # noqa: BLE001
                    cross = {"error": repr(exc)}

            warn = " | ".join(cross.get("warnings", [])) if cross else ""
            rec["ligands"].append({
                "comp_id": comp,
                "lipid_class": lig["lipid_class"],
                "ligand_chain": lig["chain_id"],
                "resseq": lig["resseq"],
                "n_heavy_atoms": lig["n_heavy_atoms"],
                "site_chain": site_chain,
                "site_chain_is_ligand_chain": site_chain == lig["chain_id"],
                "contacts": full,
                "n_contacts_at_primary_all_chains": len(at_primary),
                "chains_engaged": sorted(by_chain),
                "multi_chain_site": len(by_chain) > 1,
                "chain_id_collision": "residue-id collision across chains" in warn,
                "insertion_codes_dropped": "insertion codes dropped" in warn,
                "crosscheck_vs_interfaces": cross,
                "error": None,
            })
        for ch in sorted({l["site_chain"] for l in rec["ligands"]
                          if not l.get("error")}):
            rec["sasa"][ch] = _sasa_by_residue(model, ch)
    except Exception as exc:  # noqa: BLE001 - unparseable entry
        rec["error"] = repr(exc)

    rec["seconds"] = round(time.perf_counter() - t0, 2)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    tmp.write_text(json.dumps(rec))
    tmp.replace(dest)
    return rec


# --------------------------------------------------------------------------
# stage 2: fingerprints (ECFP4 from the PDB chemical component dictionary)
# --------------------------------------------------------------------------
def fingerprint(comp_id):
    """ECFP4 (Morgan r=2, 2048 bits) for a PDB component, or None."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    mol = I.ccd_mol(comp_id)
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return gen.GetFingerprint(mol)
    except Exception:  # noqa: BLE001
        return None


def tanimoto(a, b):
    from rdkit import DataStructs

    return float(DataStructs.TanimotoSimilarity(a, b))


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def set_metrics(pred, obs):
    """precision / recall / F1 / Jaccard of a predicted residue set against O.

    An empty prediction scores 0.0 precision by convention (see module
    docstring), not NaN: a method that names no residue earns no credit.
    """
    pred, obs = set(pred), set(obs)
    inter = len(pred & obs)
    prec = inter / len(pred) if pred else 0.0
    rec = inter / len(obs) if obs else float("nan")
    denom = prec + rec
    f1 = (2 * prec * rec / denom) if denom > 0 and math.isfinite(denom) else 0.0
    union = len(pred | obs)
    return {"precision": prec, "recall": rec, "f1": f1,
            "jaccard": inter / union if union else float("nan"),
            "n_pred": len(pred), "n_obs": len(obs), "n_intersect": inter,
            "empty_prediction": not pred}


def _nan_arm(k, obs, note):
    return {"precision": float("nan"), "recall": float("nan"),
            "f1": float("nan"), "jaccard": float("nan"),
            "n_pred": k, "n_obs": len(set(obs)), "n_intersect": float("nan"),
            "empty_prediction": True, "note": note}


def _mean_arm(sets, obs, extra=None):
    """Mean set_metrics over a list of predicted sets (a resampled arm)."""
    if not sets:
        return None
    acc = defaultdict(float)
    for s in sets:
        m = set_metrics(s, obs)
        for key in METRICS + ["n_intersect", "n_pred"]:
            acc[key] += m[key]
    out = {k: acc[k] / len(sets) for k in METRICS + ["n_intersect", "n_pred"]}
    out.update({"n_obs": len(set(obs)), "empty_prediction": False,
                "n_samples": len(sets)})
    if extra:
        out.update(extra)
    return out


def random_arm(pool, k, obs, n_draws, seed):
    """Mean set_metrics over n_draws uniform draws of k residues from `pool`."""
    if k <= 0 or not pool:
        return _nan_arm(k, obs, "no core residues to match in size" if k <= 0
                        else "no surface residues found")
    rng = np.random.default_rng(seed)
    arr = np.asarray(pool, dtype=object)
    k_eff = min(k, arr.size)
    sets = [rng.choice(arr, size=k_eff, replace=False).tolist()
            for _ in range(n_draws)]
    out = _mean_arm(sets, obs)
    out["n_draws"] = n_draws
    out["note"] = (f"pool of {arr.size} residues < k={k}" if k_eff < k else None)
    return out


def buried_arm(sasa, k, obs, min_sasa=None):
    """The k MOST BURIED residues of the binding chain, by isolated-chain SASA.

    A shape-only null: it needs no ligand and no pocket finder, only the fold.
    With min_sasa set, the ranking runs over SOLVENT-EXPOSED residues only, so
    it returns the least-exposed residues a ligand could actually reach -- the
    "any concave spot" floor. Without it, the ranking runs over every residue
    and returns the deep core, which no ligand can touch; both are reported,
    because which one "most buried" means is itself a choice.
    """
    if k <= 0:
        return _nan_arm(k, obs, "no core residues to match in size")
    pool = ({r: v for r, v in sasa.items() if v >= min_sasa}
            if min_sasa is not None else sasa)
    if not pool:
        return _nan_arm(k, obs, "no SASA computed for the binding chain"
                        if not sasa else f"no residue with SASA >= {min_sasa}")
    order = sorted(pool, key=lambda r: (pool[r], r))
    m = set_metrics(order[:k], obs)
    m["note"] = f"{len(pool)} residues in the ranking pool"
    return m


def pocket_residues(pockets, site_chain, rank):
    """Residue KEYS of the rank-N P2Rank pocket that lie in the binding chain.

    P2Rank carries no insertion code, so its keys are bare auth_seq_ids and an
    insertion-coded truth residue can never be matched. That is a real property
    of the baseline, counted in `insertion_code_audit`, not papered over.
    """
    for p in pockets:
        if p["rank"] == rank:
            return [str(r[1]) for r in p["residues"]
                    if r[0] == site_chain and isinstance(r[1], int)]
    return None


# --------------------------------------------------------------------------
# ligand table
# --------------------------------------------------------------------------
def build_ligand_table(manifest, contacts, min_ligands, drop_lipid_class):
    """Per target: one contact record per distinct ligand, best resolution wins."""
    targets, exclusions = {}, []
    for uniprot, t in manifest["targets"].items():
        res_of = {e["pdb_id"].upper(): e.get("resolution") for e in t["entries"]}
        best, dupes, failed = {}, defaultdict(list), []
        n_lipid_dropped, lipid_codes = 0, defaultdict(int)
        for pdb_id in sorted(res_of):
            rec = contacts.get(pdb_id)
            if rec is None or rec.get("error"):
                failed.append({"pdb_id": pdb_id,
                               "why": "not cached" if rec is None else rec["error"]})
                continue
            for lig in rec["ligands"]:
                if lig.get("error"):
                    failed.append({"pdb_id": pdb_id, "comp_id": lig["comp_id"],
                                   "why": lig["error"]})
                    continue
                if lig.get("lipid_class"):
                    lipid_codes[lig["comp_id"]] += 1
                    if drop_lipid_class:
                        n_lipid_dropped += 1
                        continue
                comp = lig["comp_id"]
                r = res_of.get(pdb_id)
                score = (r if r is not None else 99.0, pdb_id)
                if comp not in best or score < best[comp][0]:
                    if comp in best:
                        dupes[comp].append(best[comp][1]["pdb_id"])
                    best[comp] = (score, {**lig, "pdb_id": pdb_id, "resolution": r})
                else:
                    dupes[comp].append(pdb_id)
        ligands = {c: v[1] for c, v in best.items()}
        entry = {
            "uniprot": uniprot, "gene": t.get("gene"),
            "n_distinct_ligands_available_in_pdb": t.get("n_distinct_ligands_available"),
            "n_entries_cached": len(res_of),
            "n_distinct_ligands_with_contacts": len(ligands),
            "n_duplicate_entries_collapsed": sum(len(v) for v in dupes.values()),
            "n_ligand_extraction_failures": len(failed),
            "ligand_extraction_failures": failed[:50],
            "n_lipid_class_instances_seen": sum(lipid_codes.values()),
            "lipid_class_codes_seen": dict(sorted(lipid_codes.items())),
            "n_lipid_class_instances_dropped": n_lipid_dropped,
            "median_resolution": float(np.median([v for v in res_of.values()
                                                  if v is not None]))
                                 if any(v is not None for v in res_of.values())
                                 else None,
            "ligands": ligands,
        }
        if len(ligands) < min_ligands:
            exclusions.append({"uniprot": uniprot, "gene": t.get("gene"),
                               "n_distinct_ligands_with_contacts": len(ligands),
                               "why": f"< min_ligands={min_ligands}"})
            entry["included"] = False
        else:
            entry["included"] = True
        targets[uniprot] = entry
    return targets, exclusions


# --------------------------------------------------------------------------
# one trial
# --------------------------------------------------------------------------
def site_keys(lig, cutoff, legacy=False):
    """Residue keys of one ligand's contacts, restricted to its binding chain.

    legacy=True reconstructs the OLD bare-auth_seq_id behaviour: residue ids are
    merged across chains first (first chain by sorted key wins, which is the
    deterministic stand-in for the original atom-order-dependent first-wins) and
    only then filtered to the binding chain. That is what makes chain A 350 and
    chain B 350 one element and collapses 60A..60I onto 60.
    """
    site = lig["site_chain"]
    hit = {k: v for k, v in lig["contacts"].items() if v["min_distance"] <= cutoff}
    if not legacy:
        return {k.split("|", 1)[1] for k in hit if v_chain(k) == site}
    owner = {}
    for k in sorted(hit):
        seq = key_seq(k.split("|", 1)[1])
        owner.setdefault(seq, v_chain(k))
    return {seq for seq, ch in owner.items() if ch == site}


def v_chain(k):
    return k.split("|", 1)[0]


def per_residue_for_consensus(lig, cutoff):
    """interfaces.consensus_signature input for one ligand, binding chain only."""
    site = lig["site_chain"]
    return {k.split("|", 1)[1]: {"resname": v["resname"], "chain_id": site,
                                 "interaction_types": [],
                                 "buried_area_proxy": None,
                                 "min_distance": v["min_distance"]}
            for k, v in lig["contacts"].items()
            if v["min_distance"] <= cutoff and v["chain_id"] == site}


def rank_consensus(freqs, per_res_list):
    """Residue keys ranked for the size-matched rule.

    Frequency first; ties broken by the mean minimum heavy-atom distance over
    the ligands that engage the residue (tighter contacts first), then by key.
    """
    dsum, dn = defaultdict(float), defaultdict(int)
    for pr in per_res_list:
        for k, v in pr.items():
            dsum[k] += v["min_distance"]
            dn[k] += 1
    return sorted(freqs, key=lambda k: (-freqs[k], dsum[k] / max(dn[k], 1), k))


def trial_arms(h, others, pockets, fps, comps, uniprot, held, args, cutoff,
               legacy=False, arms="all"):
    """Score every arm for one held-out ligand. Returns (arms, diagnostics)."""
    obs = site_keys(h, cutoff, legacy=legacy)
    per_res = [per_residue_for_consensus(o, cutoff) for o in others]
    if legacy:
        per_res = [{key_seq(k): v for k, v in pr.items()} for pr in per_res]
    # interfaces._per_residue casts residue ids to int, so a key that carries an
    # insertion code ("60A") cannot be passed through the library directly. Every
    # key of this trial is therefore encoded as an integer id for the call and
    # decoded straight back; `enc` is a bijection, so nothing is merged. This is
    # also what lets the core_coverage identity be asserted against the library
    # rather than reimplemented.
    enc = {k: i for i, k in enumerate(sorted({k for pr in per_res for k in pr} | obs))}
    dec = {i: k for k, i in enc.items()}
    contact_dicts = [{
        "per_residue": {enc[k]: v for k, v in pr.items()},
        "ligand_info": {"resname": o["comp_id"], "chain": o["site_chain"],
                        "resseq": o["resseq"], "n_heavy_atoms": o["n_heavy_atoms"]},
        "buried_area_computed": False,
        "structure_id": o["pdb_id"],
    } for pr, o in zip(per_res, others)]
    sig = I.consensus_signature(contact_dicts, core_threshold=args.core_threshold,
                                source="known_ligand")
    core = [dec[r] for r in sig["core_residue_ids"]]
    freqs = {dec[r["residue_id"]]: r["frequency"] for r in sig["residues"]}
    order = rank_consensus(freqs, per_res)
    union = sorted(freqs)

    out, diag = {}, {}
    out["consensus"] = set_metrics(core, obs)
    k_o = len(obs)
    out["consensus_size_matched"] = set_metrics(order[:k_o], obs)
    diag["size_matched_tie_at_cut"] = bool(
        0 < k_o < len(order) and freqs[order[k_o - 1]] == freqs[order[k_o]])
    # The size-matched rule needs a tie-break, and ties at the cut are common,
    # so the rule is also scored with the tie-break REMOVED (frequency, then
    # residue key). If the two agree the tie-break is not carrying the result.
    alt = sorted(freqs, key=lambda k: (-freqs[k], k))
    out["consensus_size_matched_alt_tiebreak"] = set_metrics(alt[:k_o], obs)
    out["union_all"] = set_metrics(union, obs)

    # --- P2Rank -------------------------------------------------------------
    pk = pockets.get(h["pdb_id"])
    if pk is None:
        diag["no_p2rank_data"] = True
        out["p2rank_top1"] = None
        out["p2rank_best_of_3"] = None
    else:
        chains = {r[0] for pp in pk for r in pp["residues"]}
        diag["binding_chain_in_no_p2rank_pocket"] = bool(pk and h["site_chain"]
                                                         not in chains)
        p1 = pocket_residues(pk, h["site_chain"], 1)
        diag["p2rank_predicted_no_pocket"] = p1 is None
        if p1 is None:
            p1 = []
        diag["p2rank_top1_empty_in_binding_chain"] = (not p1) and not diag[
            "p2rank_predicted_no_pocket"]
        out["p2rank_top1"] = set_metrics(p1, obs)
        cands = []
        for rank in (1, 2, 3):
            rr = pocket_residues(pk, h["site_chain"], rank)
            if rr is None:
                continue
            m = set_metrics(rr, obs)
            m["rank"] = rank
            cands.append(m)
        if cands:
            best = max(cands, key=lambda m: (m["jaccard"]
                                             if math.isfinite(m["jaccard"]) else -1.0))
        else:
            best = set_metrics([], obs)
            best["rank"] = None
        out["p2rank_best_of_3"] = best

    if arms == "core":
        return out, diag, freqs, core, obs, per_res, enc, sig

    # --- rung 3: ONE random other ligand ------------------------------------
    other_sets = [site_keys(o, cutoff, legacy=legacy) for o in others]
    base = zlib.crc32(f"{uniprot}:{held}".encode())
    rng = np.random.default_rng(base ^ (args.seed + 1))
    idx = rng.integers(0, len(other_sets), size=args.null_resamples)
    out["null1_one_random_other_ligand"] = _mean_arm(
        [other_sets[i] for i in idx], obs,
        {"n_candidate_ligands": len(other_sets),
         "n_resamples": args.null_resamples})
    # the exact expectation over all other ligands, as a check on the sampling
    exact = _mean_arm(other_sets, obs)
    diag["null1_sampling_minus_exact"] = {
        m: out["null1_one_random_other_ligand"][m] - exact[m] for m in METRICS}

    # --- rung 4: the most chemically dissimilar other ligand ----------------
    fp_h = fps.get(held)
    sims = {}
    if fp_h is not None:
        sims = {c: tanimoto(fp_h, fps[c]) for c in comps
                if c != held and fps.get(c) is not None}
    if sims:
        worst = min(sorted(sims), key=lambda c: sims[c])
        pos = {o["comp_id"]: i for i, o in enumerate(others)}
        m = set_metrics(other_sets[pos[worst]], obs)
        m["ligand"] = worst
        m["tanimoto"] = sims[worst]
        out["null2_most_dissimilar_other_ligand"] = m
    else:
        out["null2_most_dissimilar_other_ligand"] = None
        diag["null2_unavailable"] = True

    # --- rungs 1-2: shape-only floors ---------------------------------------
    sasa = h.get("_sasa") or {}
    surface = sorted(r for r, s in sasa.items() if s >= SURFACE_MIN_SASA)
    out["random"] = random_arm(surface, len(core), obs, args.random_draws,
                               base ^ args.seed)
    out["random_size_matched"] = random_arm(surface, k_o, obs, args.random_draws,
                                            base ^ (args.seed + 2))
    out["buried"] = buried_arm(sasa, len(core), obs)
    out["buried_size_matched"] = buried_arm(sasa, k_o, obs)
    out["buried_surface_size_matched"] = buried_arm(sasa, k_o, obs,
                                                    min_sasa=SURFACE_MIN_SASA)
    diag["max_tanimoto"] = max(sims.values()) if sims else None
    diag["mean_tanimoto"] = float(np.mean(list(sims.values()))) if sims else None
    diag["n_tanimoto_pairs"] = len(sims)
    return out, diag, freqs, core, obs, per_res, enc, sig


def run_trials(targets, pockets, fps, args, cutoff, legacy=False, arms="all"):
    trials, max_dev = [], 0.0
    notes = defaultdict(int, {
        "trials_without_p2rank_data": 0,
        "trials_p2rank_predicted_no_pocket": 0,
        "trials_p2rank_top1_empty_in_binding_chain": 0,
        "trials_binding_chain_in_no_p2rank_pocket": 0,
        "trials_with_size_matched_tie_at_cut": 0,
        "trials_without_null2": 0,
    })
    for uniprot, t in targets.items():
        if not t["included"]:
            continue
        ligands = t["ligands"]
        comps = sorted(ligands)
        for held in comps:
            h = ligands[held]
            others = [ligands[c] for c in comps if c != held]
            a, diag, freqs, core, obs, per_res, enc, sig_raw = trial_arms(
                h, others, pockets, fps, comps, uniprot, held, args, cutoff,
                legacy=legacy, arms=arms)
            notes["trials_without_p2rank_data"] += bool(diag.get("no_p2rank_data"))
            notes["trials_p2rank_predicted_no_pocket"] += bool(
                diag.get("p2rank_predicted_no_pocket"))
            notes["trials_p2rank_top1_empty_in_binding_chain"] += bool(
                diag.get("p2rank_top1_empty_in_binding_chain"))
            notes["trials_binding_chain_in_no_p2rank_pocket"] += bool(
                diag.get("binding_chain_in_no_p2rank_pocket"))
            notes["trials_with_size_matched_tie_at_cut"] += bool(
                diag.get("size_matched_tie_at_cut"))
            notes["trials_without_null2"] += bool(diag.get("null2_unavailable"))

            # does the residue key mean the same residue across entries?
            names = defaultdict(set)
            held_pr = per_residue_for_consensus(h, cutoff)
            if legacy:
                held_pr = {key_seq(k): v for k, v in held_pr.items()}
            for pr in per_res + [held_pr]:
                for k, v in pr.items():
                    names[k].add(v["resname"])
            n_disagree = sum(1 for v in names.values() if len(v) > 1)

            sweep = {}
            for th in THRESHOLD_SWEEP:
                m = set_metrics([r for r, f in freqs.items() if f >= th], obs)
                sweep[f"{th:.1f}"] = {k: m[k] for k in METRICS + ["n_pred"]}

            if arms == "all" and not legacy:
                # the primary metric of the consensus arm IS interfaces.coverage's
                # core_coverage -- asserted at runtime, not asserted by comment.
                cov = I.coverage(sig_raw, [enc[k] for k in obs])
                if cov["core_coverage"] is not None:
                    max_dev = max(max_dev,
                                  abs(cov["core_coverage"] - a["consensus"]["precision"]))

            icoded = sorted(k for k in obs if not k.isdigit() and not k.lstrip("-").isdigit())
            trials.append({
                "uniprot": uniprot, "gene": t["gene"], "held_out_ligand": held,
                "pdb_id": h["pdb_id"], "resolution": h["resolution"],
                "site_chain": h["site_chain"],
                "site_chain_is_ligand_chain": h["site_chain_is_ligand_chain"],
                "multi_chain_site": h["multi_chain_site"],
                "held_out_is_lipid_class": bool(h.get("lipid_class")),
                "n_consensus_ligands": len(others),
                "n_consensus_ligands_lipid_class":
                    sum(1 for o in others if o.get("lipid_class")),
                "n_consensus_ligands_from_same_entry":
                    sum(1 for o in others if o["pdb_id"] == h["pdb_id"]),
                "n_consensus_ligands_other_site_chain_letter":
                    sum(1 for o in others if o["site_chain"] != h["site_chain"]),
                "n_observed_residues": len(obs),
                "n_observed_residues_with_insertion_code": len(icoded),
                "observed_insertion_coded_residues": icoded[:20],
                "n_core_residues": len(core),
                "consensus_threshold_sweep": sweep,
                "n_union_residues": len(freqs),
                "n_residue_id_disagreements": n_disagree,
                "size_matched_tie_at_cut": bool(diag.get("size_matched_tie_at_cut")),
                "max_tanimoto": diag.get("max_tanimoto"),
                "mean_tanimoto": diag.get("mean_tanimoto"),
                "n_tanimoto_pairs": diag.get("n_tanimoto_pairs", 0),
                "null1_sampling_minus_exact": diag.get("null1_sampling_minus_exact"),
                "arms": a,
            })
    return trials, dict(notes), max_dev


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def vec(trials, arm, metric):
    return np.array([
        (t["arms"][arm][metric] if t["arms"].get(arm) is not None else float("nan"))
        for t in trials], dtype=float)


def compare(trials, a, b, metric, n_boot, seed):
    va, vb = vec(trials, a, metric), vec(trials, b, metric)
    keep = np.isfinite(va) & np.isfinite(vb)
    if keep.sum() == 0:
        return None
    out = {"arm_a": a, "arm_b": b, "metric": metric}
    out.update(M.paired_bootstrap(va, vb, n_boot=n_boot, seed=seed))
    w = M.wilcoxon(va, vb)
    out["wilcoxon_p"] = w["p_value"]
    out["wilcoxon_n_effective"] = w["n_effective"]
    out["mean_a"] = float(np.mean(va[keep]))
    out["mean_b"] = float(np.mean(vb[keep]))
    return out


def compare_vectors(a, b, n_boot, seed):
    a, b = np.asarray(a, float), np.asarray(b, float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() == 0:
        return None
    out = M.paired_bootstrap(a, b, n_boot=n_boot, seed=seed)
    out["wilcoxon"] = M.wilcoxon(a, b)
    out["mean_a"] = float(np.mean(a[keep]))
    out["mean_b"] = float(np.mean(b[keep]))
    return out


def summarise_arms(trials, arms=ARMS):
    out = {}
    for arm in arms:
        row = {}
        for metric in METRICS:
            v = vec(trials, arm, metric)
            v = v[np.isfinite(v)]
            mean, lo, hi = M.bootstrap_ci(v) if v.size else (float("nan"),) * 3
            row[metric] = {"mean": mean, "ci_lo": lo, "ci_hi": hi, "n": int(v.size)}
        sizes = np.array([t["arms"][arm]["n_pred"] for t in trials
                          if t["arms"].get(arm) is not None], dtype=float)
        row["mean_n_predicted"] = float(np.mean(sizes)) if sizes.size else float("nan")
        out[arm] = row
    return out


def per_target_means(trials, arms=ARMS):
    """{uniprot: {arm: {metric: mean}}} -- the cluster-aware level."""
    by = defaultdict(list)
    for t in trials:
        by[t["uniprot"]].append(t)
    out = {}
    for u, ts in by.items():
        row = {}
        for arm in arms:
            row[arm] = {}
            for m in METRICS:
                v = vec(ts, arm, m)
                v = v[np.isfinite(v)]
                row[arm][m] = float(np.mean(v)) if v.size else float("nan")
        out[u] = row
    return out


def decomposition(trials, n_boot, seed, consensus_arm, metric, level="per_trial"):
    """How much of the consensus's margin over P2Rank one ligand already buys.

    total       = consensus - p2rank_top1
    single      = null1     - p2rank_top1   (what knowing ONE binding event buys)
    aggregation = consensus - null1         (what AGGREGATING buys on top)

    All three are computed on the SAME trials (those where all three arms are
    finite), so total = single + aggregation exactly; the residual is reported.
    """
    null = "null1_one_random_other_ligand"
    if level == "per_target":
        pt = per_target_means(trials)
        keys = sorted(pt)
        va = np.array([pt[u][consensus_arm][metric] for u in keys])
        vn = np.array([pt[u][null][metric] for u in keys])
        vb = np.array([pt[u]["p2rank_top1"][metric] for u in keys])
    else:
        va, vn, vb = (vec(trials, consensus_arm, metric), vec(trials, null, metric),
                      vec(trials, "p2rank_top1", metric))
    keep = np.isfinite(va) & np.isfinite(vn) & np.isfinite(vb)
    if keep.sum() == 0:
        return None
    va, vn, vb = va[keep], vn[keep], vb[keep]
    total = compare_vectors(va, vb, n_boot, seed)
    single = compare_vectors(vn, vb, n_boot, seed)
    agg = compare_vectors(va, vn, n_boot, seed)
    frac_single = (single["delta"] / total["delta"]) if total["delta"] else None
    return {
        "level": level, "metric": metric, "consensus_arm": consensus_arm,
        "n": int(keep.sum()),
        "mean_p2rank_top1": float(vb.mean()),
        "mean_one_random_other_ligand": float(vn.mean()),
        "mean_consensus": float(va.mean()),
        "delta_total_consensus_minus_p2rank": total,
        "delta_single_ligand_null_minus_p2rank": single,
        "delta_aggregation_consensus_minus_single_ligand": agg,
        "additivity_residual": float(total["delta"] - single["delta"] - agg["delta"]),
        "fraction_of_total_already_delivered_by_one_ligand": frac_single,
        "fraction_of_total_added_by_aggregation":
            (1.0 - frac_single) if frac_single is not None else None,
        "reading": "the one-random-other-ligand null needs no aggregation, no "
                   "generator and no pocket finder; whatever it already scores "
                   "is the part of the consensus's margin that is NOT about "
                   "aggregating many ligands. Quote "
                   "delta_aggregation_consensus_minus_single_ligand as the "
                   "margin aggregation buys.",
    }


def _jsonable(obj):
    """NaN/Inf -> null, so the file is strict JSON and an unmeasured value reads
    as absent rather than as a number."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if math.isfinite(obj) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--min-ligands", type=int, default=5)
    ap.add_argument("--core-threshold", type=float, default=0.6)
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--min-heavy-atoms", type=int, default=10)
    ap.add_argument("--random-draws", type=int, default=100)
    ap.add_argument("--null-resamples", type=int, default=40)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit-targets", type=int, default=0)
    ap.add_argument("--no-crosscheck", action="store_true",
                    help="skip the per-ligand agreement check against "
                         "interfaces.ligand_contacts (affects the cache)")
    ap.add_argument("--contacts-only", action="store_true",
                    help="build the contact cache and stop")
    args = ap.parse_args()

    t_start = time.time()
    manifest = json.loads(MANIFEST.read_text())
    if args.limit_targets:
        keys = list(manifest["targets"])[:args.limit_targets]
        manifest = {**manifest, "targets": {k: manifest["targets"][k] for k in keys}}

    wanted = sorted({e["pdb_id"].upper()
                     for t in manifest["targets"].values() for e in t["entries"]})
    jobs, missing_files = [], []
    for pdb_id in wanted:
        hits = [p for p in (STRUCTURES / f"{pdb_id}.cif",
                            STRUCTURES / f"{pdb_id.lower()}.cif") if p.exists()]
        if not hits:
            missing_files.append(pdb_id)
            continue
        jobs.append((pdb_id, str(hits[0]), args.min_heavy_atoms, args.cutoff,
                     not args.no_crosscheck))
    print(f"{len(wanted)} entries in manifest, {len(jobs)} cached on disk, "
          f"{len(missing_files)} missing", flush=True)

    contacts, t0 = {}, time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(extract_one, j): j[0] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            contacts[rec["pdb_id"]] = rec
            if i % 100 == 0 or i == len(jobs):
                print(f"  contacts {i}/{len(jobs)}  ({time.time() - t0:.0f}s)",
                      flush=True)
    contact_seconds = time.time() - t0
    n_parse_fail = sum(1 for r in contacts.values() if r.get("error"))
    print(f"contacts done in {contact_seconds:.0f}s, {n_parse_fail} structures failed",
          flush=True)
    if args.contacts_only:
        return

    # agreement of this module's contact pass with interfaces.ligand_contacts
    cross = {"n_ligands_checked": 0, "n_ligands_identical": 0,
             "n_residues_only_in_library": 0, "n_residues_only_here": 0,
             "n_ligands_not_checked": 0}
    for rec in contacts.values():
        for lig in rec.get("ligands", []):
            c = lig.get("crosscheck_vs_interfaces")
            if not c or "error" in c:
                cross["n_ligands_not_checked"] += 1
                continue
            cross["n_ligands_checked"] += 1
            cross["n_residues_only_in_library"] += c["n_only_library"]
            cross["n_residues_only_here"] += c["n_only_here"]
            cross["n_ligands_identical"] += int(c["n_only_library"] == 0
                                                and c["n_only_here"] == 0)

    # ---------------- ligand tables, lipids out (primary) and in ------------
    targets, exclusions = build_ligand_table(manifest, contacts, args.min_ligands,
                                             drop_lipid_class=True)
    targets_lip, _ = build_ligand_table(manifest, contacts, args.min_ligands,
                                        drop_lipid_class=False)
    for tt in (targets, targets_lip):
        for t in tt.values():
            for lig in t["ligands"].values():
                lig["_sasa"] = contacts[lig["pdb_id"]]["sasa"].get(lig["site_chain"], {})
    included = [u for u, t in targets.items() if t["included"]]
    print(f"{len(included)}/{len(targets)} targets included "
          f"(>= {args.min_ligands} distinct ligands), {len(exclusions)} excluded",
          flush=True)

    p2 = json.loads(P2RANK_JSON.read_text()) if P2RANK_JSON.exists() else {"pockets": {}}
    pockets = p2.get("pockets", {})
    # verified, not assumed: what fields does a P2Rank residue actually carry?
    widths = {len(r) for pp in pockets.values() for p in pp for r in p["residues"]}
    p2_fields = {"residue_tuple_widths": sorted(widths),
                 "carries_chain": True, "carries_insertion_code": False,
                 "carries_per_residue_score": False,
                 "checked": "every residue entry in results/p2rank_pockets.json is "
                            "[chain, auth_seq_id]; there is no third element, so "
                            "neither an insertion code nor a per-residue score is "
                            "available and the size-matched rule CANNOT be applied "
                            "to the P2Rank arms"}
    print(f"P2Rank pockets for {len(pockets)} structures, residue tuple widths "
          f"{sorted(widths)}", flush=True)

    comps = sorted({c for u in included for c in targets_lip[u]["ligands"]})
    print(f"fetching ECFP4 for {len(comps)} distinct ligands ...", flush=True)
    fps, t0 = {}, time.time()
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(fingerprint, c): c for c in comps}
        for i, fut in enumerate(as_completed(futs), 1):
            fps[futs[fut]] = fut.result()
            if i % 250 == 0 or i == len(comps):
                print(f"  ccd {i}/{len(comps)} ({time.time() - t0:.0f}s)", flush=True)
    n_fp_fail = sum(1 for c in comps if fps.get(c) is None)
    print(f"{n_fp_fail}/{len(comps)} ligands without a usable CCD structure",
          flush=True)

    # ---------------- the primary run ---------------------------------------
    trials, notes, max_dev = run_trials(targets, pockets, fps, args, args.cutoff)
    print(f"{len(trials)} trials (lipid/detergent class excluded)", flush=True)

    overall = summarise_arms(trials)
    family, comparisons = {}, []
    for i, a in enumerate(ARMS):
        for b in ARMS[i + 1:]:
            for metric in METRICS:
                c = compare(trials, a, b, metric, args.n_boot, args.seed)
                if c is None:
                    continue
                comparisons.append(c)
                family[f"{a}_vs_{b}::{metric}"] = c["p_value"]
    holm_boot = M.holm_bonferroni(family)
    holm_wil = M.holm_bonferroni({f"{c['arm_a']}_vs_{c['arm_b']}::{c['metric']}":
                                  c["wilcoxon_p"] for c in comparisons})
    for c in comparisons:
        k = f"{c['arm_a']}_vs_{c['arm_b']}::{c['metric']}"
        c["p_holm"] = holm_boot[k]
        c["wilcoxon_p_holm"] = holm_wil[k]
    print(f"{len(comparisons)} arm comparisons, Holm-corrected", flush=True)

    headline = next((c for c in comparisons
                     if c["arm_a"] == "consensus" and c["arm_b"] == "p2rank_top1"
                     and c["metric"] == PRIMARY), None)
    if headline is None:
        sys.exit("no trial has both a consensus and a P2Rank arm; the decisive "
                 "comparison cannot be made. Check results/p2rank_pockets.json.")

    # ---------------- the null ladder ---------------------------------------
    pt_means = per_target_means(trials)
    ladder = {}
    for arm in LADDER + ["union_all", "p2rank_best_of_3"]:
        row = {"arm": arm,
               "per_trial": {m: overall[arm][m] for m in METRICS},
               "mean_n_predicted": overall[arm]["mean_n_predicted"],
               "consensus_vs_this": {}, "consensus_size_matched_vs_this": {}}
        for cons, dest in (("consensus", "consensus_vs_this"),
                           ("consensus_size_matched", "consensus_size_matched_vs_this")):
            for metric in METRICS:
                row[dest][metric] = {
                    "per_trial": compare(trials, cons, arm, metric, args.n_boot,
                                         args.seed),
                    "per_target": compare_vectors(
                        [pt_means[u][cons][metric] for u in sorted(pt_means)],
                        [pt_means[u][arm][metric] for u in sorted(pt_means)],
                        args.n_boot, args.seed),
                }
        ladder[arm] = row
    # Holm across the ladder family: consensus vs every rung, both set rules,
    # every metric, at BOTH levels, and on the Wilcoxon p as well as the
    # bootstrap p. The bootstrap p has a resolution floor of 2/(n_boot+1)
    # (metrics.paired_bootstrap), so a whole family can sit on that floor and
    # read as non-significant for arithmetic rather than statistical reasons --
    # which is why the Wilcoxon column is corrected alongside it.
    for level, pkey in (("per_trial", "p_value"), ("per_target", "p_value")):
        fam, famw = {}, {}
        for arm, row in ladder.items():
            for dest in ("consensus_vs_this", "consensus_size_matched_vs_this"):
                for metric, cell in row[dest].items():
                    c = cell.get(level)
                    if not c:
                        continue
                    k = f"{dest}:{arm}::{metric}"
                    fam[k] = c[pkey]
                    famw[k] = (c["wilcoxon_p"] if "wilcoxon_p" in c
                               else c["wilcoxon"]["p_value"])
        h, hw = M.holm_bonferroni(fam), M.holm_bonferroni(famw)
        for arm, row in ladder.items():
            for dest in ("consensus_vs_this", "consensus_size_matched_vs_this"):
                for metric, cell in row[dest].items():
                    c = cell.get(level)
                    k = f"{dest}:{arm}::{metric}"
                    if c and k in h:
                        c["p_holm_within_ladder"] = h[k]
                        c["wilcoxon_p_holm_within_ladder"] = hw[k]

    # ---------------- the decomposition -------------------------------------
    decomp = {}
    for cons in ("consensus", "consensus_size_matched"):
        for metric in (PRIMARY, "jaccard"):
            for level in ("per_trial", "per_target"):
                d = decomposition(trials, args.n_boot, args.seed, cons, metric, level)
                if d:
                    decomp[f"{cons}::{metric}::{level}"] = d
    novel = [t for t in trials
             if t["max_tanimoto"] is not None and t["max_tanimoto"] < 0.3]
    decomp["chemically_novel_stratum_max_tanimoto_lt_0.3"] = {
        "n_trials": len(novel),
        "why": "the stratum where the held-out ligand is chemically unlike every "
               "ligand in the consensus, so a single-ligand null has the least to "
               "go on",
        "consensus_vs_one_random_other_ligand": (
            compare(novel, "consensus", "null1_one_random_other_ligand", PRIMARY,
                    args.n_boot, args.seed) if novel else None),
        "consensus_size_matched_vs_one_random_other_ligand": (
            compare(novel, "consensus_size_matched", "null1_one_random_other_ligand",
                    PRIMARY, args.n_boot, args.seed) if novel else None),
        "consensus_vs_p2rank_top1": (
            compare(novel, "consensus", "p2rank_top1", PRIMARY, args.n_boot,
                    args.seed) if novel else None),
    }

    # ---------------- both set-construction rules, side by side -------------
    rules = {
        "researcher_degree_of_freedom": (
            "A consensus is a frequency map; turning it into a SET needs a rule, "
            "and the rule is a free parameter that can change the answer. Both "
            "rules are reported here for every arm where they apply, and the "
            "core@threshold sweep is reported as well, precisely so no reader "
            "inherits an unregistered choice. Neither rule is privileged by this "
            "artifact beyond core@threshold being PROJECT_GOAL.md 4.3's default, "
            "which is why it carries the name `consensus`."),
        "rules": {
            RULE_CORE: "keep residues with consensus frequency >= core_threshold "
                       f"(reported at {args.core_threshold})",
            RULE_SIZE: "keep the |O| highest-frequency residues, |O| = the size of "
                       "the held-out ligand's own contact set; ties at the cut "
                       "broken by mean minimum heavy-atom distance, then key",
        },
        "arms": {},
        "not_applicable": {
            "p2rank_top1 / p2rank_best_of_3":
                "P2Rank residues in results/p2rank_pockets.json are [chain, "
                "auth_seq_id] with NO per-residue score (verified: "
                "inputs.p2rank_residue_fields), so the pocket cannot be truncated "
                "to its |O| best residues. The size-matched rule is therefore NOT "
                "applied to the baseline and the size asymmetry of the P2Rank arm "
                "is left standing and reported (mean_n_predicted).",
            "null1 / null2":
                "a single ligand's contact set has no ranking to truncate; its "
                "size is what it is.",
            "union_all":
                "the union IS the frequency>0 set; its size-matched truncation is "
                "exactly the consensus_size_matched arm.",
        },
    }
    for base, pair in (("consensus", ("consensus", "consensus_size_matched")),
                       ("random", ("random", "random_size_matched")),
                       ("buried", ("buried", "buried_size_matched")),
                       ("buried_surface", (None, "buried_surface_size_matched"))):
        rules["arms"][base] = {}
        for rule, arm in zip((RULE_CORE, RULE_SIZE), pair):
            if arm is None:     # no |core|-sized variant of this arm was scored
                rules["arms"][base][rule] = None
                continue
            rules["arms"][base][rule] = {
                "arm": arm,
                "per_trial": {m: overall[arm][m] for m in METRICS},
                "mean_n_predicted": overall[arm]["mean_n_predicted"],
                "vs_p2rank_top1_primary": compare(trials, arm, "p2rank_top1",
                                                  PRIMARY, args.n_boot, args.seed),
                "vs_p2rank_top1_jaccard": compare(trials, arm, "p2rank_top1",
                                                  "jaccard", args.n_boot, args.seed),
            }
    tb = compare(trials, "consensus_size_matched",
                 "consensus_size_matched_alt_tiebreak", PRIMARY, args.n_boot,
                 args.seed)
    rules["size_matched_tie_break"] = {
        "primary_tie_break": "mean minimum heavy-atom distance across the "
                             "engaging ligands (tighter first), then residue key",
        "alternative_tie_break": "residue key alone (no distance term)",
        "n_trials_with_a_frequency_tie_straddling_the_cut":
            notes["trials_with_size_matched_tie_at_cut"],
        "of_n_trials": len(trials),
        "primary_minus_alternative": tb,
        "reading": "ties at the cut are common, so this says how much of the "
                   "size-matched arm is the ranking and how much is the "
                   "tie-break.",
    }
    a_core = rules["arms"]["consensus"][RULE_CORE]["vs_p2rank_top1_primary"]
    a_size = rules["arms"]["consensus"][RULE_SIZE]["vs_p2rank_top1_primary"]
    rules["two_rules_agree_in_sign"] = bool(
        math.copysign(1, a_core["delta"]) == math.copysign(1, a_size["delta"]))
    # "excludes zero" means the whole interval is on one side of it. Testing
    # (ci_lo > 0) == (ci_lo > 0) would also return True when NEITHER interval
    # excluded zero, which is the opposite reading, so each rule is tested on
    # its own interval and the weaker "same side" statement is kept separately.
    rules["two_rules_both_exclude_zero"] = bool(
        (a_core["ci_lo"] > 0 or a_core["ci_hi"] < 0)
        and (a_size["ci_lo"] > 0 or a_size["ci_hi"] < 0))
    rules["two_rules_ci_on_the_same_side_of_zero"] = bool(
        (a_core["ci_lo"] > 0) == (a_size["ci_lo"] > 0))

    # ---------------- sensitivity subsets -----------------------------------
    defined = [t for t in trials
               if t["arms"].get("p2rank_top1") is not None
               and not t["arms"]["p2rank_top1"]["empty_prediction"]
               and not t["arms"]["consensus"]["empty_prediction"]]
    icode_genes = sorted({t["gene"] for t in trials
                          if t["n_observed_residues_with_insertion_code"]})
    subsets = {
        "primary_defined_only": [t for t in defined],
        "numbering_clean_only": [t for t in trials
                                 if t["n_residue_id_disagreements"] == 0],
        "single_chain_sites_only": [t for t in trials if not t["multi_chain_site"]],
        "no_same_entry_consensus_ligand":
            [t for t in trials if t["n_consensus_ligands_from_same_entry"] == 0],
        "drop_ADORA2A": [t for t in trials if t["gene"] != "ADORA2A"],
        "drop_targets_with_insertion_coded_contacts":
            [t for t in trials if t["gene"] not in icode_genes],
        "drop_trials_whose_consensus_mixes_site_chain_letters":
            [t for t in trials
             if t["n_consensus_ligands_other_site_chain_letter"] == 0],
    }
    sens = {}
    for name, subset in subsets.items():
        sens[name] = {
            "n_trials": len(subset),
            "n_targets": len({t["uniprot"] for t in subset}),
            "consensus_vs_p2rank_top1": (
                compare(subset, "consensus", "p2rank_top1", PRIMARY,
                        args.n_boot, args.seed) if subset else None),
            "consensus_size_matched_vs_p2rank_top1": (
                compare(subset, "consensus_size_matched", "p2rank_top1", PRIMARY,
                        args.n_boot, args.seed) if subset else None),
            "consensus_vs_one_random_other_ligand": (
                compare(subset, "consensus", "null1_one_random_other_ligand",
                        PRIMARY, args.n_boot, args.seed) if subset else None),
        }
    sens["drop_targets_with_insertion_coded_contacts"]["genes_dropped"] = icode_genes

    # ---------------- contact-cutoff sweep ----------------------------------
    # The ground truth O is "residues within `cutoff` of a ligand heavy atom",
    # which is itself a choice. Contacts are extracted once at EXTRACT_CUTOFF
    # with each residue's minimum distance kept, so a tighter cutoff is a filter,
    # not a re-parse. The binding chain and the SASA pools stay fixed at the
    # primary cutoff so that the only thing varying across the sweep is the
    # cutoff itself.
    cutoffs = {"_note": (
        "consensus, union and the single-ligand nulls all move with the cutoff "
        "because both the prediction and the truth are contact sets; the P2Rank "
        "arm does not move, because a pocket is not a contact set. The binding "
        "chain and the surface/buried pools are held at the primary cutoff "
        f"({args.cutoff} A) so the sweep varies one thing.")}
    for cut in CUTOFF_SWEEP:
        if cut == args.cutoff:
            tr = trials
        else:
            tr, _, _ = run_trials(targets, pockets, fps, args, cut, arms="core")
        cutoffs[f"{cut:.1f}"] = {
            "cutoff_angstrom": cut,
            "n_trials": len(tr),
            "mean_n_observed_residues": float(np.mean([t["n_observed_residues"]
                                                       for t in tr])),
            "consensus_vs_p2rank_top1": compare(tr, "consensus", "p2rank_top1",
                                                PRIMARY, args.n_boot, args.seed),
            "consensus_size_matched_vs_p2rank_top1": compare(
                tr, "consensus_size_matched", "p2rank_top1", PRIMARY,
                args.n_boot, args.seed),
        }
        print(f"  cutoff {cut} delta="
              f"{cutoffs[f'{cut:.1f}']['consensus_vs_p2rank_top1']['delta']:+.4f}",
              flush=True)

    # ---------------- legacy residue keys -----------------------------------
    leg_trials, _, _ = run_trials(targets, pockets, fps, args, args.cutoff,
                                  legacy=True, arms="core")
    leg_by = {(t["uniprot"], t["held_out_ligand"]): t for t in leg_trials}
    changed, deltas, changed_genes = 0, [], defaultdict(int)
    for t in trials:
        o = leg_by.get((t["uniprot"], t["held_out_ligand"]))
        if o is None:
            continue
        d = [abs(t["arms"][a][PRIMARY] - o["arms"][a][PRIMARY])
             for a in ("consensus", "p2rank_top1")
             if t["arms"].get(a) is not None and o["arms"].get(a) is not None]
        d.append(abs(t["n_observed_residues"] - o["n_observed_residues"]))
        if any(x > 1e-12 for x in d):
            changed += 1
            changed_genes[t["gene"]] += 1
            deltas.append(max(d[:-1]) if len(d) > 1 else 0.0)
    key_audit = {
        "what_changed": "residues are keyed by (chain, auth_seq_id, insertion "
                        "code) instead of bare auth_seq_id, on BOTH arms",
        "legacy_reconstruction": "bare auth_seq_id, merged across chains first "
                                 "(deterministic first-by-sorted-key wins) and "
                                 "filtered to the binding chain afterwards, which "
                                 "is what the previous version did",
        "n_trials_compared": len(leg_trials),
        "n_trials_changed": changed,
        "n_trials_changed_per_gene": dict(sorted(changed_genes.items(),
                                                 key=lambda kv: -kv[1])),
        "mean_abs_primary_change_on_changed_trials":
            float(np.mean(deltas)) if deltas else 0.0,
        "max_abs_primary_change_on_changed_trials":
            float(np.max(deltas)) if deltas else 0.0,
        "headline_under_legacy_keys": compare(leg_trials, "consensus",
                                              "p2rank_top1", PRIMARY,
                                              args.n_boot, args.seed),
        "headline_under_chain_and_icode_keys": headline,
    }
    ic_trials = [t for t in trials if t["n_observed_residues_with_insertion_code"]]
    insertion_audit = {
        "n_trials_with_insertion_coded_contact_residues": len(ic_trials),
        "genes": sorted({t["gene"] for t in ic_trials}),
        "n_insertion_coded_residues_total":
            int(sum(t["n_observed_residues_with_insertion_code"] for t in ic_trials)),
        "mean_share_of_truth_set": float(np.mean(
            [t["n_observed_residues_with_insertion_code"] / t["n_observed_residues"]
             for t in ic_trials])) if ic_trials else 0.0,
        "p2rank_cannot_match_them": "P2Rank residues carry no insertion code, so "
                                    "for these trials part of the truth set is "
                                    "unreachable for the baseline; the "
                                    "drop_targets_with_insertion_coded_contacts "
                                    "sensitivity removes them entirely",
        "consensus_vs_p2rank_top1_on_these_trials": (
            compare(ic_trials, "consensus", "p2rank_top1", PRIMARY, args.n_boot,
                    args.seed) if ic_trials else None),
    }
    mixed = [t for t in trials if t["n_consensus_ligands_other_site_chain_letter"]]
    chain_audit = {
        "n_trials_whose_consensus_mixes_binding_chain_letters": len(mixed),
        "n_targets_affected": len({t["uniprot"] for t in mixed}),
        "why_the_letter_is_normalised":
            "the binding chain letter is a per-deposition label: entry 1ABC chain "
            "A and entry 2DEF chain B are the same entity. Every scored set is "
            "restricted to one chain, so the letter carries no within-trial "
            "information; keying by the literal letter instead would have made "
            "these trials' consensus and truth sets DISJOINT by construction.",
        "sensitivity_dropping_them":
            sens["drop_trials_whose_consensus_mixes_site_chain_letters"],
    }

    # ---------------- lipid / detergent exclusion ---------------------------
    lip_trials, lip_notes, _ = run_trials(targets_lip, pockets, fps, args,
                                          args.cutoff, arms="core")
    per_target_lipid = {}
    for u, t in targets_lip.items():
        if not t["included"]:
            continue
        seen = t["lipid_class_codes_seen"]
        if not seen:
            continue
        kept_as_ligand = sorted(c for c in t["ligands"]
                                if t["ligands"][c].get("lipid_class"))
        per_target_lipid[t["gene"] or u] = {
            "n_lipid_class_instances_seen": sum(seen.values()),
            "codes_seen": seen,
            "n_distinct_lipid_class_ligands_removed": len(kept_as_ligand),
            "codes_removed_as_trials_and_consensus_members": kept_as_ligand,
            "n_trials_before": sum(1 for x in lip_trials if x["uniprot"] == u),
            "n_trials_after": sum(1 for x in trials if x["uniprot"] == u),
        }
    lipid_block = {
        "class": I.LIPID_DETERGENT_REASON,
        "where": "scripts/interfaces.py EXCLUDED_HET, named class "
                 "LIPID_DETERGENT_REASON",
        "n_codes_in_class": sum(1 for r in I.EXCLUDED_HET.values()
                                if r.startswith(I.LIPID_DETERGENT_REASON_PREFIX)),
        "codes_in_class": sorted(c for c, r in I.EXCLUDED_HET.items()
                                 if r.startswith(I.LIPID_DETERGENT_REASON_PREFIX)),
        "codes_present_in_this_dataset": sorted(
            {c for t in targets_lip.values() for c in t["lipid_class_codes_seen"]}),
        "codes_in_class_absent_from_this_dataset": sorted(
            {c for c, r in I.EXCLUDED_HET.items()
             if r.startswith(I.LIPID_DETERGENT_REASON_PREFIX)}
            - {c for t in targets_lip.values() for c in t["lipid_class_codes_seen"]}),
        "verified_how": "every code was matched by name against the local PDB "
                        "chemical component dictionary (data/raw/ccd/*.cif); "
                        "presence/absence in this dataset is counted from the "
                        "contact cache, not assumed",
        "removed_from_consensus_as_well_as_from_trials": True,
        "per_target": per_target_lipid,
        "n_trials_with_lipids_in": len(lip_trials),
        "n_trials_with_lipids_out": len(trials),
        "headline_with_lipids_in": compare(lip_trials, "consensus", "p2rank_top1",
                                           PRIMARY, args.n_boot, args.seed),
        "headline_with_lipids_out": headline,
        "empty_core_trials_with_lipids_in":
            int(sum(1 for t in lip_trials if t["n_core_residues"] == 0)),
        "empty_core_trials_with_lipids_out":
            int(sum(1 for t in trials if t["n_core_residues"] == 0)),
        "component_level_exclusion_is_target_blind": {
            "what": "the class is a list of component ids, so it cannot know that "
                    "a fatty acid is a membrane lipid in a GPCR and a genuine "
                    "orthosteric ligand in a fatty-acid receptor",
            "collateral_removals": {
                g: v["codes_removed_as_trials_and_consensus_members"]
                for g, v in per_target_lipid.items()
                if g in ("PPARA", "PPARD", "PPARG", "RXRA", "NR1H4")},
        },
    }

    # ---------------- core-threshold sweep ----------------------------------
    sweep_summary = {}
    for th in THRESHOLD_SWEEP:
        key = f"{th:.1f}"
        va = np.array([t["consensus_threshold_sweep"][key][PRIMARY] for t in trials])
        vb = vec(trials, "p2rank_top1", PRIMARY)
        keep = np.isfinite(va) & np.isfinite(vb)
        row = {
            "core_threshold": th,
            "consensus_mean_primary": float(np.mean(va[np.isfinite(va)])),
            "consensus_mean_f1": float(np.mean(
                [t["consensus_threshold_sweep"][key]["f1"] for t in trials])),
            "consensus_mean_n_predicted": float(np.mean(
                [t["consensus_threshold_sweep"][key]["n_pred"] for t in trials])),
            "n_trials_with_empty_core": int(sum(
                1 for t in trials if t["consensus_threshold_sweep"][key]["n_pred"] == 0)),
        }
        if keep.sum():
            row["vs_p2rank_top1"] = M.paired_bootstrap(va[keep], vb[keep],
                                                       n_boot=args.n_boot,
                                                       seed=args.seed)
            row["vs_p2rank_top1"]["wilcoxon"] = M.wilcoxon(va[keep], vb[keep])
        sweep_summary[key] = row

    empty_core = [t for t in trials if t["n_core_residues"] == 0]
    empty_core_by_target = defaultdict(int)
    for t in empty_core:
        empty_core_by_target[t["gene"]] += 1
    empty_core_lip = defaultdict(int)
    for t in lip_trials:
        if t["n_core_residues"] == 0:
            empty_core_lip[t["gene"]] += 1

    # ---------------- per target --------------------------------------------
    by_target = {}
    per_target_trials = defaultdict(list)
    for t in trials:
        per_target_trials[t["uniprot"]].append(t)
    for u, ts in per_target_trials.items():
        row = {"uniprot": u, "gene": targets[u]["gene"], "n_trials": len(ts),
               "n_duplicate_entries_collapsed": targets[u]["n_duplicate_entries_collapsed"],
               "n_distinct_ligands_available_in_pdb":
                   targets[u]["n_distinct_ligands_available_in_pdb"],
               "median_resolution": targets[u]["median_resolution"],
               "mean_n_observed_residues": float(np.mean([x["n_observed_residues"]
                                                          for x in ts])),
               "mean_n_core_residues": float(np.mean([x["n_core_residues"]
                                                      for x in ts])),
               "n_trials_with_empty_core": int(sum(1 for x in ts
                                                   if x["n_core_residues"] == 0)),
               "n_lipid_class_instances_seen":
                   targets[u]["n_lipid_class_instances_seen"],
               "arms": {}}
        for arm in ARMS:
            row["arms"][arm] = {}
            for metric in METRICS:
                v = vec(ts, arm, metric)
                v = v[np.isfinite(v)]
                row["arms"][arm][metric] = float(np.mean(v)) if v.size else float("nan")
        row["delta_consensus_minus_p2rank_primary"] = (
            row["arms"]["consensus"][PRIMARY] - row["arms"]["p2rank_top1"][PRIMARY])
        row["delta_consensus_minus_one_random_other_ligand_primary"] = (
            row["arms"]["consensus"][PRIMARY]
            - row["arms"]["null1_one_random_other_ligand"][PRIMARY])
        mt = [x["max_tanimoto"] for x in ts if x["max_tanimoto"] is not None]
        row["mean_max_tanimoto"] = float(np.mean(mt)) if mt else float("nan")
        by_target[u] = row

    def target_weighted(arm, metric):
        v = np.array([by_target[u]["arms"][arm][metric] for u in by_target],
                     dtype=float)
        v = v[np.isfinite(v)]
        mean, lo, hi = M.bootstrap_ci(v) if v.size else (float("nan"),) * 3
        return {"mean": mean, "ci_lo": lo, "ci_hi": hi, "n_targets": int(v.size)}

    tw = {arm: {m: target_weighted(arm, m) for m in METRICS} for arm in ARMS}
    keys = sorted(by_target)
    tw_delta = compare_vectors([by_target[u]["arms"]["consensus"][PRIMARY] for u in keys],
                               [by_target[u]["arms"]["p2rank_top1"][PRIMARY] for u in keys],
                               args.n_boot, args.seed)
    agree = (tw_delta is not None
             and math.copysign(1, tw_delta["delta"]) == math.copysign(1, headline["delta"]))

    # resolution confound: is the per-target margin explained by crystal quality?
    from scipy import stats as _st

    rr = [(by_target[u]["median_resolution"],
           by_target[u]["delta_consensus_minus_p2rank_primary"]) for u in keys
          if by_target[u]["median_resolution"] is not None
          and math.isfinite(by_target[u]["delta_consensus_minus_p2rank_primary"])]
    sp = _st.spearmanr([x for x, _ in rr], [y for _, y in rr]) if len(rr) > 2 else None
    resolution_confound = {
        "n_targets": len(rr),
        "spearman_rho": float(sp.statistic) if sp is not None else None,
        "p_value": float(sp.pvalue) if sp is not None else None,
        "x": "per-target median resolution of the cached entries",
        "y": "per-target mean (consensus - p2rank_top1) precision",
        "reading": "a strong positive rho would mean the margin is a crystal-"
                   "quality artifact rather than an information difference",
    }

    # ---------------- negative results, recomputed --------------------------
    losing = sorted(
        [{"gene": by_target[u]["gene"], "uniprot": u,
          "n_trials": by_target[u]["n_trials"],
          "n_trials_with_empty_core": by_target[u]["n_trials_with_empty_core"],
          "consensus": by_target[u]["arms"]["consensus"][PRIMARY],
          "p2rank_top1": by_target[u]["arms"]["p2rank_top1"][PRIMARY],
          "one_random_other_ligand":
              by_target[u]["arms"]["null1_one_random_other_ligand"][PRIMARY],
          "consensus_size_matched":
              by_target[u]["arms"]["consensus_size_matched"][PRIMARY],
          "delta": by_target[u]["delta_consensus_minus_p2rank_primary"]}
         for u in keys
         if by_target[u]["delta_consensus_minus_p2rank_primary"] < 0],
        key=lambda r: r["delta"])
    losing_null1 = sorted(
        [{"gene": by_target[u]["gene"],
          "delta": by_target[u]["delta_consensus_minus_one_random_other_ligand_primary"]}
         for u in keys
         if by_target[u]["delta_consensus_minus_one_random_other_ligand_primary"] < 0],
        key=lambda r: r["delta"])
    recall_cmp = next(c for c in comparisons if c["arm_a"] == "consensus"
                      and c["arm_b"] == "p2rank_top1" and c["metric"] == "recall")
    f1_cmp = next(c for c in comparisons if c["arm_a"] == "consensus"
                  and c["arm_b"] == "p2rank_top1" and c["metric"] == "f1")
    negatives = {
        "consensus_loses_on_recall_to_p2rank_top1": {
            "point_estimate_is_a_loss": bool(recall_cmp["delta"] < 0),
            "ci_excludes_zero": bool(recall_cmp["ci_lo"] > 0
                                     or recall_cmp["ci_hi"] < 0),
            "delta": recall_cmp["delta"], "ci": [recall_cmp["ci_lo"], recall_cmp["ci_hi"]],
            "p_value": recall_cmp["p_value"], "p_holm": recall_cmp["p_holm"],
            "wilcoxon_p": recall_cmp["wilcoxon_p"],
            "wilcoxon_p_holm": recall_cmp["wilcoxon_p_holm"],
            "n": recall_cmp["n"],
            "under_the_size_matched_rule": next(
                (c for c in comparisons
                 if c["arm_a"] == "consensus_size_matched"
                 and c["arm_b"] == "p2rank_top1" and c["metric"] == "recall"),
                None),
            "note": "the core@threshold rule under-predicts by construction; it "
                    "buys precision with recall. Read point_estimate_is_a_loss "
                    "together with ci_excludes_zero: a negative point estimate "
                    "whose interval spans zero is not a demonstrated loss, and "
                    "saying so is the whole point of reporting both."},
        "consensus_and_p2rank_disagree_between_precision_and_f1": {
            "is_true": bool(math.copysign(1, headline["delta"])
                            != math.copysign(1, f1_cmp["delta"])),
            "primary_delta": headline["delta"], "f1_delta": f1_cmp["delta"]},
        "targets_where_p2rank_beats_the_consensus": {
            "n_targets": len(losing), "of_n_targets": len(keys), "targets": losing,
            "all_of_them_are_empty_core_targets": bool(losing) and all(
                r["n_trials_with_empty_core"] > 0.5 * r["n_trials"] for r in losing),
            "empty_core_target_criterion":
                "a target counts as an empty-core target when MORE THAN HALF its "
                "trials predict no residue at core_threshold; it is not the claim "
                "that every trial of every such target is empty. Each row's "
                "n_trials_with_empty_core / n_trials is printed beside it so the "
                "share can be read directly.",
            "note": "where this is true the consensus does not lose by naming the "
                    "WRONG residues, it loses by naming NONE: at core_threshold "
                    "no residue clears the bar, the empty prediction scores 0 by "
                    "convention, and the target's whole column is zero. The "
                    "consensus_size_matched column of each row shows what the "
                    "same evidence gives when the set rule cannot return nothing, "
                    "which separates 'the signature is wrong here' from 'the "
                    "threshold rule fires never here'."},
        "targets_where_one_random_other_ligand_beats_the_consensus": {
            "n_targets": len(losing_null1), "of_n_targets": len(keys),
            "targets": losing_null1},
        "trials_scoring_zero_because_the_core_is_empty": {
            "n_trials": len(empty_core), "per_gene": dict(sorted(
                empty_core_by_target.items(), key=lambda kv: -kv[1]))},
        "p2rank_oracle_best_of_3_is_about_the_single_ligand_null": {
            "p2rank_top1": overall["p2rank_top1"][PRIMARY]["mean"],
            "p2rank_best_of_3_oracle_selected":
                overall["p2rank_best_of_3"][PRIMARY]["mean"],
            "one_random_other_ligand":
                overall["null1_one_random_other_ligand"][PRIMARY]["mean"],
            "best_of_3_minus_top1": compare(trials, "p2rank_best_of_3",
                                            "p2rank_top1", PRIMARY, args.n_boot,
                                            args.seed),
            "best_of_3_minus_one_random_other_ligand": compare(
                trials, "p2rank_best_of_3", "null1_one_random_other_ligand",
                PRIMARY, args.n_boot, args.seed),
            "note": "p2rank_best_of_3 picks the pocket USING the answer, so it is "
                    "not a fair arm and is never the baseline. Read as a "
                    "diagnostic it is a negative result for the consensus: most "
                    "of the gap between p2rank_top1 and the arms above it is "
                    "recovered by letting P2Rank's own top 3 be oracle-ranked, so "
                    "much of what the consensus beats p2rank_top1 by is pocket "
                    "RANKING rather than pocket FINDING."},
        "buried_null_is_near_the_floor": {
            "mean_primary_whole_chain": overall["buried_size_matched"][PRIMARY]["mean"],
            "mean_primary_exposed_residues_only":
                overall["buried_surface_size_matched"][PRIMARY]["mean"],
            "mean_primary_uniform_random":
                overall["random_size_matched"][PRIMARY]["mean"],
            "note": "ranking by absolute isolated-chain SASA over the WHOLE chain "
                    "selects deep-core residues, which no ligand can reach, so "
                    "that version of rung 2 sits at or below the uniform-random "
                    "floor rather than above it. Restricting the same ranking to "
                    "solvent-exposed residues (the 'any concave spot' reading) is "
                    "reported beside it. Either way the finding is the same: "
                    "shape alone, without a cavity detector, buys little, which "
                    "is why rung 3 and not rung 2 is the null that matters."},
    }

    # ---------------- stratified by chemical novelty -------------------------
    have = [t for t in trials if t["max_tanimoto"] is not None]
    strata = {}
    if have:
        q = np.quantile([t["max_tanimoto"] for t in have], [0.25, 0.5, 0.75])
        edges = [(-0.01, q[0]), (q[0], q[1]), (q[1], q[2]), (q[2], 1.01)]
        for i, (lo, hi) in enumerate(edges, 1):
            sub = [t for t in have if lo < t["max_tanimoto"] <= hi]
            if not sub:
                continue
            strata[f"Q{i}"] = {
                "max_tanimoto_range": [float(lo), float(hi)],
                "n_trials": len(sub),
                "arms": {arm: {m: float(np.nanmean(vec(sub, arm, m))) for m in METRICS}
                         for arm in ARMS},
                "consensus_vs_p2rank_top1": compare(sub, "consensus", "p2rank_top1",
                                                    PRIMARY, args.n_boot, args.seed),
                "consensus_vs_one_random_other_ligand": compare(
                    sub, "consensus", "null1_one_random_other_ligand", PRIMARY,
                    args.n_boot, args.seed),
            }

    null1_signed = np.array([t["null1_sampling_minus_exact"][PRIMARY]
                             for t in trials
                             if t.get("null1_sampling_minus_exact")], dtype=float)
    null1_dev = np.abs(null1_signed).tolist()
    # Absolute deviation alone cannot say whether the seeded sampling moved the
    # headline decomposition, because it hides the sign. The rung-3 arm mean and
    # the aggregation margin are therefore also recomputed against the EXACT
    # expectation over all other ligands (sampled minus its own deviation).
    _n1 = vec(trials, "null1_one_random_other_ligand", PRIMARY)
    _cs = vec(trials, "consensus", PRIMARY)
    _ok = np.isfinite(_n1) & np.isfinite(_cs)
    null1_exact_effect = {
        "rung3_mean_sampled": float(_n1[_ok].mean()),
        "rung3_mean_exact": float((_n1 - null1_signed)[_ok].mean()),
        "aggregation_delta_vs_sampled_rung3": float((_cs - _n1)[_ok].mean()),
        "aggregation_delta_vs_exact_rung3":
            float((_cs - (_n1 - null1_signed))[_ok].mean()),
        "n": int(_ok.sum()),
    }

    out = {
        "experiment": "leave-one-ligand-out hotspot recovery; PROJECT_GOAL.md "
                      "ablation I6.2 (known_ligand signature -- THE CEILING), "
                      "NOT I6.1 and NOT the M2 gate",
        "read_this_first": {
            "claim": "the consensus beats P2Rank on the primary metric, but a "
                     "null that costs nothing -- ONE other ligand of the same "
                     "target, chosen at random -- already delivers most of that "
                     "margin. Most of what any arm here predicts is that the "
                     "pocket is the pocket.",
            "where": "aggregation_decomposition",
            "never_quote": "the consensus's margin over the random-surface floor, "
                           "or its margin over P2Rank, without the "
                           "one-random-other-ligand rung beside it",
        },
        "what_this_is_not": {
            "is_ablation": "I6.2",
            "is_not_ablation": "I6.1",
            "m2_gate_status": "neither passed nor failed",
            "why": "PROJECT_GOAL.md I6 defines I6.1 as 'BoltzGen signature vs "
                   "p2rank_geometry signature -- *the* experiment' and I6.2 as "
                   "'Signature vs known_ligand signature -- the ceiling'; 8.3 "
                   "states that source='known_ligand' gives the UPPER BOUND. "
                   "The boltzgen_consensus arm was NOT run, so the M2 gate is "
                   "not passed and not failed by this result. What the result "
                   "bounds is the ceiling that absent arm would have had to "
                   "reach. A negative result here would have foreclosed the "
                   "BoltzGen arm a fortiori (no generator beats P2Rank through "
                   "a ceiling that does not itself beat P2Rank); this positive "
                   "result forecloses nothing and licenses nothing about "
                   "BoltzGen.",
        },
        "primary_metric": PRIMARY,
        "primary_metric_definition":
            "|predicted hotspot residues INTERSECT observed contact residues of the "
            "held-out ligand| / |predicted hotspot residues|; for the consensus arm "
            "this is exactly interfaces.coverage(sig, O)['core_coverage'] "
            "(PROJECT_GOAL.md 4.4)",
        "core_coverage_identity_max_deviation": max_dev,
        "generated_arm": {
            "present": False,
            "why": "No BoltzGen design ensemble was available to this script, so no "
                   "generated arm was scored. The consensus is built from real "
                   "co-crystals instead: the source='known_ligand' builder of "
                   "PROJECT_GOAL.md C8.",
            "recorded_reason_not_reverified_here": {
                "claim": "the Modal account authenticates but every GPU tier tried "
                         "(T4, L4, A10G, A100) was refused pending a payment method",
                "source": "docs/04-BOLTZGEN-MODAL.md",
                "verified_on": "2026-09-19",
                "caveat": "this script does not probe Modal; the GPU status may have "
                          "changed since. Re-verify with scripts/modal_boltzgen.py "
                          "before quoting it.",
            },
        },
        "parameters": vars(args),
        "conventions": {
            "empty_prediction_precision": 0.0,
            "binding_chain": "polymer chain contributing the most contact residues",
            "one_copy_per_component_per_entry": True,
            "one_entry_per_component_per_target": "best resolution",
            "buried_area_computed": False,
            "p2rank_best_of_3_is_oracle_selected": True,
            "residue_key":
                "(chain_id, auth_seq_id, insertion_code). Every scored set is "
                "restricted to the binding chain, so the key is written "
                "'<auth_seq_id><insertion_code>' within that chain; contacts in "
                "other chains are DROPPED, never merged in. See "
                "residue_key_audit and insertion_code_audit for what this "
                "changed, and site_chain_letter_audit for why the chain LETTER "
                "is not compared literally across entries.",
            "numbering": "author (auth_seq_id + insertion code), no SIFTS mapping",
            "lipid_detergent_class_excluded": True,
            "benchmark_set_selection":
                "Entries were chosen greedily to MAXIMISE distinct ligands per "
                "target (<= 25 entries/target, <= 3.0 A) by "
                "scripts/fetch_structures.py. by_target."
                "n_distinct_ligands_available_in_pdb shows what was left on the "
                "table (ESR1: 434 available, 26 used). The set is therefore "
                "DELIBERATELY diversity-maximising and NOT representative of the "
                "PDB: any near-duplicate rate measured on it is a LOWER bound and "
                "must not be quoted as a property of the PDB. Diversity "
                "maximisation makes each trial HARDER for the consensus, so it "
                "makes a positive result conservative and a null partly an "
                "artifact of deliberate hardness.",
            "random_arm_denominator":
                "The `random` and `buried` arms are size-matched to |core|, so on "
                "the trials where the core is empty they have no set to draw and "
                "score NaN rather than 0. Their n is therefore SMALLER than every "
                "other arm's (see overall_per_trial.random.*.n), and every "
                "comparison involving them is pairwise-complete on that easier "
                "subset, which FLATTERS the consensus in those rows: see "
                "all_arm_comparisons_holm_corrected[consensus vs random].mean_a, "
                "which is the consensus mean on that subset, not its overall "
                "mean. Read those rows as a floor, not as a paired result. The "
                "`random_size_matched` and `buried_size_matched` arms are sized "
                "to |O|, which is never empty, so they are defined on every trial "
                "and are the versions to quote.",
            "null1_sampling":
                f"rung 3 averages {args.null_resamples} seeded draws of one other "
                "ligand per trial. The exact expectation over ALL other ligands is "
                "computed too and the sampling error is reported "
                "(null1_sampling_vs_exact).",
        },
        "inputs": {
            "n_entries_in_manifest": len(wanted),
            "n_entries_cached_on_disk": len(jobs),
            "n_entries_missing_from_disk": len(missing_files),
            "entries_missing_from_disk": missing_files[:50],
            "n_structures_failed_to_parse": n_parse_fail,
            "n_structures_with_p2rank_pockets": len(pockets),
            "p2rank_source": str(P2RANK_JSON),
            "p2rank_n_structures_parsed": p2.get("n_structures_parsed"),
            "p2rank_residue_fields": p2_fields,
            "contact_pass_crosscheck_vs_interfaces_ligand_contacts": cross,
        },
        "targets": {
            "n_in_manifest": len(targets),
            "n_included": len(included),
            "n_excluded": len(exclusions),
            "exclusions": exclusions,
        },
        "n_trials": len(trials),
        "notes": notes,
        "overall_per_trial": overall,
        "overall_per_target_weighted": tw,
        "headline_consensus_vs_p2rank_top1": headline,
        "per_target_weighted_consensus_vs_p2rank_top1": tw_delta,
        "per_trial_and_per_target_agree_in_sign": agree,
        "aggregation_decomposition": decomp,
        "null_ladder": ladder,
        "set_construction_rules": rules,
        "all_arm_comparisons_holm_corrected": comparisons,
        "bootstrap_p_resolution_floor": {
            "floor": 2.0 / (args.n_boot + 1),
            "n_comparisons_in_the_family": len(comparisons),
            "n_comparisons_sitting_on_the_floor":
                int(sum(1 for c in comparisons
                        if c["p_value"] <= 2.0 / (args.n_boot + 1) + 1e-12)),
            "holm_of_the_floor": (2.0 / (args.n_boot + 1)) * len(comparisons),
            "what": "metrics.paired_bootstrap cannot resolve a p below "
                    "2/(n_boot+1), so a comparison that is overwhelmingly "
                    "significant reports exactly that floor and Holm multiplies it "
                    "by the family size. A p_holm near or above 0.05 in this family "
                    "is therefore a RESOLUTION limit, not weak evidence -- read "
                    "wilcoxon_p_holm, which has no floor, alongside it. The "
                    "headline's wilcoxon p is ~1e-112.",
        },
        "sensitivity": sens,
        "contact_cutoff_sweep": cutoffs,
        "residue_key_audit": key_audit,
        "insertion_code_audit": insertion_audit,
        "site_chain_letter_audit": chain_audit,
        "lipid_detergent_exclusion": lipid_block,
        "resolution_confound": resolution_confound,
        "negative_results": negatives,
        "core_threshold_sweep": sweep_summary,
        "null1_sampling_vs_exact": {
            "max_abs_deviation_primary": float(np.max(null1_dev)) if null1_dev else None,
            "mean_abs_deviation_primary": float(np.mean(null1_dev)) if null1_dev else None,
            "mean_signed_deviation_primary":
                float(null1_signed.mean()) if null1_signed.size else None,
            "n_trials": len(null1_dev),
            "effect_on_the_decomposition": null1_exact_effect,
            "what": "the seeded resampled rung-3 estimate minus the exact mean over "
                    "all other ligands of that target",
            "reading": "per trial the deviation can be large (one target's other "
                       "ligands are a small, heterogeneous pool), so the SIGNED "
                       "mean and effect_on_the_decomposition are what say whether "
                       "the sampling moved the reported aggregation margin; "
                       "compare aggregation_delta_vs_sampled_rung3 with "
                       "aggregation_delta_vs_exact_rung3 against the CI width in "
                       "aggregation_decomposition.",
        },
        "empty_core_diagnostics": {
            "n_trials_with_empty_core_at_reported_threshold": len(empty_core),
            "per_gene": dict(sorted(empty_core_by_target.items(),
                                    key=lambda kv: -kv[1])),
            "with_lipid_class_left_in": {
                "n_trials": int(sum(1 for t in lip_trials if t["n_core_residues"] == 0)),
                "per_gene": dict(sorted(empty_core_lip.items(),
                                        key=lambda kv: -kv[1]))},
            "reduction_per_gene_from_the_lipid_exclusion": {
                g: empty_core_lip[g] - empty_core_by_target.get(g, 0)
                for g in sorted(empty_core_lip)},
            "note": "An empty core means no residue is engaged by >= core_threshold "
                    "of the target's other ligands, so the arm predicts nothing and "
                    "is scored 0 rather than dropped. Removing the membrane "
                    "lipids, detergents and organomercurials was expected to be the "
                    "cure; MEASURED, it is not -- see "
                    "reduction_per_gene_from_the_lipid_exclusion, which is a "
                    "single-gene effect. The rest is the one-site assumption in the "
                    "known_ligand consensus builder failing on targets whose cached "
                    "ligands genuinely occupy more than one site, and it is a "
                    "standing negative result, not a fixed one.",
        },
        "not_evaluated": [
            {"what": "BoltzGen-generated interface signature (PROJECT_GOAL.md 4.3, "
                     "6.I, the arm the M2 gate was written for)",
             "why": "no design ensemble was available to this script. The recorded "
                    "reason (docs/04-BOLTZGEN-MODAL.md, verified 2026-09-19, NOT "
                    "re-verified here) is that the Modal account authenticates but "
                    "every GPU tier tried was refused pending a payment method. "
                    "Everything here uses the source='known_ligand' builder, so it "
                    "bounds what a design ensemble would have to beat, and says "
                    "nothing about BoltzGen itself. This is why the run is ablation "
                    "I6.2 and not I6.1, and why the M2 gate is neither passed nor "
                    "failed."},
            {"what": "size-matched truncation of the P2Rank arms (the symmetric "
                     "version of the size_matched rule)",
             "why": "results/p2rank_pockets.json stores each pocket residue as "
                    "[chain, auth_seq_id] with no per-residue score -- verified at "
                    "load, inputs.p2rank_residue_fields -- so there is no ranking "
                    "to truncate. Re-running P2Rank to keep its *_residues.csv "
                    "scores would be needed. The consequence is left standing and "
                    "visible: the P2Rank arms predict ~2x more residues than the "
                    "consensus (overall_per_trial.*.mean_n_predicted), which costs "
                    "them precision and earns them recall."},
            {"what": "rank-based, cut-free comparison (AUC / average precision of "
                     "each method's residue ranking against the binary contact "
                     "label)",
             "why": "same cause: the P2Rank arm has no per-residue ranking in the "
                    "cached artifact, so only the consensus could be scored that "
                    "way and the comparison would not be like for like."},
            {"what": "PEG oligomers and buffer components not on any exclusion "
                     "class (PE8, 1PG, 15P, P4C, BCN, CXS and similar)",
             "why": "they belong to the EXISTING 'cryoprotectant, precipitant or "
                    "buffer additive' class rather than to the lipid/detergent "
                    "class added here, and editing a published class's membership "
                    "was outside this change. They are a handful of instances; "
                    "residual_non_drug_het counts them exactly so the omission is "
                    "visible rather than silent."},
            {"what": "cross-target transfer (does a signature built on target X "
                     "predict hotspots on a homologue Y)",
             "why": "out of scope for this ablation; author numbering is not "
                    "mapped through SIFTS here, so keys are not comparable across "
                    "targets."},
            {"what": "buried-area-weighted hotspot definitions (PROJECT_GOAL.md 4.3 "
                     "mean_buried_area)",
             "why": "contacts were extracted with with_buried_area=False for "
                    "speed; every metric reported is set-membership, so no "
                    "buried-area number is computed and none is reported."},
            {"what": "whether the held-out ligand is chemically novel to the "
                     "*retrieval* half of the project",
             "why": "measured only against the same target's other co-crystal "
                    "ligands (max/mean ECFP4 Tanimoto), not against a training "
                    "corpus, because no retrieval model is trained here."},
        ],
        "stratified_by_max_tanimoto_quartile": strata,
        "fingerprints": {"n_ligands": len(comps), "n_without_structure": n_fp_fail},
        "by_target": by_target,
        "trials": trials,
        "runtime_seconds": round(time.time() - t_start, 1),
        "contact_extraction_seconds": round(contact_seconds, 1),
    }

    # residual non-drug het still in the trial set, counted not assumed
    residual = defaultdict(int)
    for u in included:
        for c in targets[u]["ligands"]:
            if c in ("PE8", "1PG", "15P", "P4C", "BCN", "CXS", "PEG", "P6G"):
                residual[c] += 1
    out["residual_non_drug_het"] = {
        "codes_counted": ["PE8", "1PG", "15P", "P4C", "BCN", "CXS", "PEG", "P6G"],
        "found_as_trials": dict(sorted(residual.items())),
        "n_trials": sum(residual.values()),
        "why_not_excluded": "see not_evaluated; they are precipitant/buffer class, "
                            "not lipid/detergent class",
    }

    OUT_JSON.parent.mkdir(exist_ok=True)
    OUT_JSON.write_text(json.dumps(_jsonable(out), indent=1, default=float))

    cols = (["uniprot", "gene", "n_trials", "n_distinct_ligands_available_in_pdb",
             "n_duplicate_entries_collapsed", "median_resolution",
             "mean_n_observed_residues", "mean_n_core_residues",
             "n_trials_with_empty_core", "n_lipid_class_instances_seen",
             "mean_max_tanimoto", "delta_consensus_minus_p2rank_primary",
             "delta_consensus_minus_one_random_other_ligand_primary"]
            + [f"{arm}::{m}" for arm in ARMS for m in METRICS])
    with OUT_CSV.open("w", newline="") as fh:
        import csv

        w = csv.writer(fh)
        w.writerow(cols)
        for u in sorted(by_target, key=lambda k: by_target[k]["gene"] or k):
            r = by_target[u]
            w.writerow([r[c] if "::" not in c
                        else r["arms"][c.split("::")[0]][c.split("::")[1]]
                        for c in cols])

    # ---------------- console ------------------------------------------------
    print("\n=== PRIMARY (precision = core coverage), per trial ===")
    for arm in ARMS:
        s = overall[arm][PRIMARY]
        print(f"  {arm:36s} {s['mean']:.4f}  [{s['ci_lo']:.4f}, {s['ci_hi']:.4f}]  "
              f"n={s['n']:4d}  mean |pred|={overall[arm]['mean_n_predicted']:.1f}")
    print(f"\nHEADLINE consensus - p2rank_top1 = {headline['delta']:+.4f} "
          f"[{headline['ci_lo']:+.4f}, {headline['ci_hi']:+.4f}]  "
          f"p={headline['p_value']:.2e} (holm {headline['p_holm']:.2e}), "
          f"n={headline['n']}; wilcoxon p={headline['wilcoxon_p']:.2e} "
          f"(n_eff={headline['wilcoxon_n_effective']})")
    if tw_delta:
        print(f"per-target-weighted delta = {tw_delta['delta']:+.4f} "
              f"[{tw_delta['ci_lo']:+.4f}, {tw_delta['ci_hi']:+.4f}] "
              f"p={tw_delta['p_value']:.2e} n_targets={tw_delta['n']}  "
              f"agrees in sign with per-trial: {agree}")

    print("\n=== DECOMPOSITION: what does AGGREGATION actually buy? ===")
    for key in (f"consensus::{PRIMARY}::per_trial",
                f"consensus::{PRIMARY}::per_target",
                f"consensus_size_matched::{PRIMARY}::per_trial"):
        d = decomp.get(key)
        if not d:
            continue
        print(f"  {key}  (n={d['n']})")
        print(f"    p2rank_top1                      {d['mean_p2rank_top1']:.4f}")
        print(f"    one random OTHER ligand (rung 3) {d['mean_one_random_other_ligand']:.4f}"
              f"   -> already +{d['delta_single_ligand_null_minus_p2rank']['delta']:.4f}"
              f" [{d['delta_single_ligand_null_minus_p2rank']['ci_lo']:+.4f},"
              f"{d['delta_single_ligand_null_minus_p2rank']['ci_hi']:+.4f}] over P2Rank")
        print(f"    consensus                        {d['mean_consensus']:.4f}"
              f"   -> aggregation adds "
              f"{d['delta_aggregation_consensus_minus_single_ligand']['delta']:+.4f}"
              f" [{d['delta_aggregation_consensus_minus_single_ligand']['ci_lo']:+.4f},"
              f"{d['delta_aggregation_consensus_minus_single_ligand']['ci_hi']:+.4f}]"
              f" p={d['delta_aggregation_consensus_minus_single_ligand']['p_value']:.1e}")
        print(f"    total {d['delta_total_consensus_minus_p2rank']['delta']:+.4f}"
              f" = single {d['delta_single_ligand_null_minus_p2rank']['delta']:+.4f}"
              f" + aggregation "
              f"{d['delta_aggregation_consensus_minus_single_ligand']['delta']:+.4f}"
              + ("  (one ligand already delivers "
                 f"{100 * d['fraction_of_total_already_delivered_by_one_ligand']:.0f}%)"
                 if d["fraction_of_total_already_delivered_by_one_ligand"] is not None
                 else "  (total margin is zero; no share to attribute)"))

    print("\n=== NULL LADDER (primary, per trial; consensus vs each rung) ===")
    for arm in LADDER:
        row = ladder[arm]
        c = row["consensus_vs_this"][PRIMARY]["per_trial"]
        ct = row["consensus_vs_this"][PRIMARY]["per_target"]
        s = row["consensus_size_matched_vs_this"][PRIMARY]["per_trial"]
        print(f"  {arm:36s} {row['per_trial'][PRIMARY]['mean']:.4f}  "
              f"core@{args.core_threshold} {c['delta']:+.4f} "
              f"[{c['ci_lo']:+.4f},{c['ci_hi']:+.4f}] p={c['p_value']:.1e} "
              f"holm={c.get('p_holm_within_ladder', float('nan')):.1e} | "
              f"per-target {ct['delta']:+.4f} "
              f"[{ct['ci_lo']:+.4f},{ct['ci_hi']:+.4f}] | "
              f"size-matched {s['delta']:+.4f}")

    print("\n=== SET-CONSTRUCTION RULES side by side (primary) ===")
    for base, cell in rules["arms"].items():
        for rule, r in cell.items():
            if r is None:
                continue
            d = r["vs_p2rank_top1_primary"]
            j = r["vs_p2rank_top1_jaccard"]
            print(f"  {base:10s} {rule:20s} mean={r['per_trial'][PRIMARY]['mean']:.4f} "
                  f"|pred|={r['mean_n_predicted']:5.1f}  vs p2rank: "
                  f"precision {d['delta']:+.4f} [{d['ci_lo']:+.4f},{d['ci_hi']:+.4f}]  "
                  f"jaccard {j['delta']:+.4f} [{j['ci_lo']:+.4f},{j['ci_hi']:+.4f}]")
    print(f"  the two rules agree in sign: {rules['two_rules_agree_in_sign']}; "
          f"both CIs exclude zero: {rules['two_rules_both_exclude_zero']}; "
          f"same side of zero: {rules['two_rules_ci_on_the_same_side_of_zero']}")
    tbk = rules["size_matched_tie_break"]
    print(f"  size-matched tie-break: ties at the cut in "
          f"{tbk['n_trials_with_a_frequency_tie_straddling_the_cut']}/"
          f"{tbk['of_n_trials']} trials; distance tie-break minus key-only "
          f"tie-break = {tbk['primary_minus_alternative']['delta']:+.4f} "
          f"[{tbk['primary_minus_alternative']['ci_lo']:+.4f},"
          f"{tbk['primary_minus_alternative']['ci_hi']:+.4f}]")

    print("\n=== LIPID / DETERGENT / HEAVY-ATOM EXCLUSION ===")
    print(f"  codes in class: {lipid_block['n_codes_in_class']}, present here: "
          f"{lipid_block['codes_present_in_this_dataset']}")
    for g, v in sorted(per_target_lipid.items(),
                       key=lambda kv: -kv[1]["n_lipid_class_instances_seen"]):
        print(f"    {g:10s} instances={v['n_lipid_class_instances_seen']:3d} "
              f"distinct ligands removed={v['n_distinct_lipid_class_ligands_removed']:2d} "
              f"trials {v['n_trials_before']} -> {v['n_trials_after']}  "
              f"{v['codes_removed_as_trials_and_consensus_members']}")
    hi_, ho_ = lipid_block["headline_with_lipids_in"], lipid_block["headline_with_lipids_out"]
    print(f"  headline WITH lipids in : {hi_['delta']:+.4f} "
          f"[{hi_['ci_lo']:+.4f},{hi_['ci_hi']:+.4f}] n={hi_['n']}")
    print(f"  headline WITH lipids out: {ho_['delta']:+.4f} "
          f"[{ho_['ci_lo']:+.4f},{ho_['ci_hi']:+.4f}] n={ho_['n']}")
    print(f"  empty cores {lipid_block['empty_core_trials_with_lipids_in']} -> "
          f"{lipid_block['empty_core_trials_with_lipids_out']}")

    print("\n=== RESIDUE KEYS (chain + insertion code) ===")
    print(f"  {key_audit['n_trials_changed']}/{key_audit['n_trials_compared']} trials "
          f"change; mean |delta primary| on those = "
          f"{key_audit['mean_abs_primary_change_on_changed_trials']:.4f}, max "
          f"{key_audit['max_abs_primary_change_on_changed_trials']:.4f}")
    print(f"  per gene: {key_audit['n_trials_changed_per_gene']}")
    print(f"  headline legacy keys {key_audit['headline_under_legacy_keys']['delta']:+.4f} "
          f"-> new keys {headline['delta']:+.4f}")
    print(f"  insertion-coded contacts in {insertion_audit['n_trials_with_insertion_coded_contact_residues']} "
          f"trials, genes {insertion_audit['genes']}")

    print("\n=== CONTACT CUTOFF SWEEP ===")
    for k, v in cutoffs.items():
        if k.startswith("_"):
            continue
        c = v["consensus_vs_p2rank_top1"]
        print(f"  {k} A  mean|O|={v['mean_n_observed_residues']:.1f}  "
              f"delta={c['delta']:+.4f} [{c['ci_lo']:+.4f},{c['ci_hi']:+.4f}]")

    print("\n=== core-threshold sweep (consensus arm) ===")
    for key, row in sweep_summary.items():
        c = row.get("vs_p2rank_top1")
        print(f"  th={key}  primary={row['consensus_mean_primary']:.4f}  "
              f"f1={row['consensus_mean_f1']:.4f}  "
              f"|pred|={row['consensus_mean_n_predicted']:5.1f}  "
              f"empty_core={row['n_trials_with_empty_core']:4d}  "
              + (f"delta_vs_p2rank_top1={c['delta']:+.4f} p={c['p_value']:.1e}" if c else ""))

    print("\n=== NEGATIVE RESULTS ===")
    print(f"  consensus loses on RECALL (point estimate): "
          f"{negatives['consensus_loses_on_recall_to_p2rank_top1']['point_estimate_is_a_loss']}"
          f", CI excludes zero: "
          f"{negatives['consensus_loses_on_recall_to_p2rank_top1']['ci_excludes_zero']} "
          f"({recall_cmp['delta']:+.4f} [{recall_cmp['ci_lo']:+.4f},"
          f"{recall_cmp['ci_hi']:+.4f}], n={recall_cmp['n']})")
    print(f"  targets where P2Rank wins: {len(losing)}/{len(keys)} -> "
          f"{[r['gene'] for r in losing]}")
    print(f"  targets where ONE random other ligand beats the consensus: "
          f"{len(losing_null1)}/{len(keys)} -> {[r['gene'] for r in losing_null1]}")
    print(f"  trials scoring 0 on an empty core: {len(empty_core)} "
          f"{dict(sorted(empty_core_by_target.items(), key=lambda kv: -kv[1]))}")
    print(f"  buried (shape-only) null primary: whole chain "
          f"{overall['buried_size_matched'][PRIMARY]['mean']:.4f}, exposed only "
          f"{overall['buried_surface_size_matched'][PRIMARY]['mean']:.4f}, "
          f"uniform random {overall['random_size_matched'][PRIMARY]['mean']:.4f}")
    if resolution_confound["spearman_rho"] is not None:
        print(f"  resolution confound: spearman rho="
              f"{resolution_confound['spearman_rho']:.3f} "
              f"p={resolution_confound['p_value']:.2f} "
              f"n={resolution_confound['n_targets']}")
    else:
        print("  resolution confound: not computable (< 3 targets with a "
              "recorded resolution)")
    print(f"  p2rank oracle best-of-3 {overall['p2rank_best_of_3'][PRIMARY]['mean']:.4f} "
          f"vs one random other ligand "
          f"{overall['null1_one_random_other_ligand'][PRIMARY]['mean']:.4f} "
          f"(top1 {overall['p2rank_top1'][PRIMARY]['mean']:.4f}): much of the "
          f"margin over top1 is pocket RANKING")
    e_ = out["null1_sampling_vs_exact"]["effect_on_the_decomposition"]
    print(f"  rung-3 sampling: mean signed deviation "
          f"{null1_signed.mean():+.5f}; aggregation margin "
          f"{e_['aggregation_delta_vs_sampled_rung3']:+.5f} sampled vs "
          f"{e_['aggregation_delta_vs_exact_rung3']:+.5f} exact")
    print(f"  core-coverage identity max deviation = {max_dev:.2e}")
    print(f"  contact-pass crosscheck: {cross['n_ligands_identical']}/"
          f"{cross['n_ligands_checked']} ligands identical to "
          f"interfaces.ligand_contacts")

    print(f"\nwrote {OUT_JSON}\nwrote {OUT_CSV}\n"
          f"{out['runtime_seconds']:.0f}s total")


if __name__ == "__main__":
    main()
