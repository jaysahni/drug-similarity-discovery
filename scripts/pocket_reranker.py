"""E2 - P2Rank finds the peptide site and then picks the wrong one. Fix the picking.

results/epitope_gate.json measured, on 70 peptide-binder complexes, that P2Rank's
rank-1 pocket recovers the true interface at Jaccard 0.1836, while the BEST pocket
it offered at any rank recovers it at 0.3112, and at least one offered pocket
touches the true interface in 78.6% of cases. The detector is not the problem. The
ranker is. That +0.128 is the largest measured headroom anywhere in this pipeline.

This script asks whether the right pocket can be picked out of P2Rank's own
candidate list WITHOUT the answer, by re-ranking the offered pockets on
truth-free descriptors:

    p2rank_rank       do nothing - P2Rank's own order. The baseline to beat.
    <rule>            one rule, CHOSEN on a development half of the peptide
                      complexes and REPORTED on the held-out half. A rule chosen
                      and reported on the same data is not a result, so the
                      headline number here is the held-out one and the
                      development number is printed beside it to show the gap.
    oracle_bestany    the best pocket at any rank, chosen BY THE TRUTH. An
                      ORACLE. It is the ceiling, it is not achievable, and it is
                      labelled as such everywhere it appears.

GROUND TRUTH is not recomputed here. It is read from the cached preparation that
epitope_gate.py wrote (data/raw/epitope_gate/prepared.json, field
interface_residues): receptor residues with any heavy atom within 4.5 A of any
binder heavy atom, in the receptor's renumbered 1..N space. Scoring is the same
Jaccard over residue-number sets that m2_gate.py and epitope_gate.py use, so
every number here is on the same axis as the numbers they report.

FEATURES are computed from the binder-stripped receptor PDB and from P2Rank's own
prediction and residue tables. None of them can see the binder:

    p2rank        rank, score, probability, sas_points, surf_atoms, and the
                  per-residue P2Rank probability mass inside the pocket
    size          number of lining residues, fraction of the receptor they are
    shape         PCA elongation of the lining CA cloud, maximum CA-CA extent
                  (a peptide site is a GROOVE, a long shallow trench; P2Rank
                  scores enclosed concave volume, which is a different object)
    burial        heavy atoms within 8 A and 12 A of the pocket centre, and a
                  ray-cast enclosure fraction over 92 directions
    exposure      Shrake-Rupley solvent accessibility of the lining residues,
                  absolute and relative to Tien et al. (2013) maxima
    position      distance of the pocket centre from the protein centroid, raw
                  and divided by the radius of gyration
    composition   hydrophobic / aromatic / charged / polar lining fraction
    secondary     helix / strand / coil lining fraction from a CA-only P-SEA
                  style assignment (approximate - no DSSP binary on this machine)
    contiguity    number of contiguous sequence runs in the lining set

The same locked rule is then run on the ppi and antibody classes. On antibodies it
is expected NOT to help, because no pocket touching the epitope is offered at all
in 69% of cases - there is nothing there to re-rank. Confirming that is the point.

Usage:
    ./env/bin/python scripts/pocket_reranker.py
    ./env/bin/python scripts/pocket_reranker.py --refresh      # recompute features
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interfaces as I  # noqa: E402
import metrics as M  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "raw" / "epitope_gate"
P2RANK_OUT = CACHE / "p2rank" / "out"
OUT = ROOT / "results" / "pocket_reranker.json"
# Derived cache, not a result: it lives under data/raw/ with every other cached
# derivation in this repo (gitignored), so the script reproduces on any machine.
FEATURE_CACHE = ROOT / "data" / "raw" / "pocket_reranker" / "pocket_features.json"

SEED = 0
N_BOOT = 10000

# Tien et al. 2013 theoretical maximum accessible surface area, A^2.
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
    "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174, "MSE": 224,
}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "MSE"}
AROMATIC = {"PHE", "TRP", "TYR", "HIS"}
CHARGED = {"ASP", "GLU", "LYS", "ARG", "HIS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}


# --------------------------------------------------------------------------
# 1. scoring - identical to epitope_gate.score / m2_gate.score
# --------------------------------------------------------------------------
def score(pred, truth):
    if not pred:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "jaccard": 0.0, "n_pred": 0}
    tp = len(pred & truth)
    p = tp / len(pred)
    r = tp / len(truth) if truth else float("nan")
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"precision": p, "recall": r, "f1": f1,
            "jaccard": tp / len(pred | truth), "n_pred": len(pred)}


def pocket_residues(pocket):
    """Residue numbers of one P2Rank pocket, in the receptor's renumbered space.

    Same rule as epitope_gate.pocket_residues: integer residue numbers only, so
    an insertion-coded token is dropped rather than guessed at. The receptors
    this reads were rewritten 1..N with no insertion codes, so nothing is lost -
    verified below by n_noninteger_residue_tokens in the output.
    """
    return {num for _chain, num in pocket["residues"] if isinstance(num, int)}


# --------------------------------------------------------------------------
# 2. P2Rank tables - the cached CSVs, including columns pockets.json drops
# --------------------------------------------------------------------------
def read_predictions(path):
    """predictions.csv -> pockets, keeping sas_points and surf_atoms.

    run_p2rank.parse_predictions drops those two columns; they are pocket-size
    descriptors this experiment needs, so the same file is read again here
    rather than changing a module another experiment is running against.
    """
    pockets, bad_tokens = [], 0
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            if not row.get("name"):
                continue
            residues = []
            for tok in (row.get("residue_ids") or "").split():
                chain, _, num = tok.rpartition("_")
                if not num:
                    continue
                try:
                    residues.append([chain, int(num)])
                except ValueError:
                    residues.append([chain, num])
                    bad_tokens += 1
            pockets.append({
                "rank": int(row["rank"]), "score": float(row["score"]),
                "probability": float(row["probability"]),
                "sas_points": float(row.get("sas_points") or "nan"),
                "surf_atoms": float(row.get("surf_atoms") or "nan"),
                "center": [float(row["center_x"]), float(row["center_y"]),
                           float(row["center_z"])],
                "residues": residues,
            })
    return sorted(pockets, key=lambda p: p["rank"]), bad_tokens


def read_residue_scores(path):
    """residues.csv -> {residue number: P2Rank per-residue probability}."""
    out = {}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            try:
                out[int(row["residue_label"])] = float(row["probability"])
            except (KeyError, ValueError):
                continue
    return out


# --------------------------------------------------------------------------
# 3. structure descriptors
# --------------------------------------------------------------------------
_SPHERE = None


def sphere_directions(n=92):
    """Roughly uniform unit vectors, for the ray-cast enclosure feature."""
    global _SPHERE
    if _SPHERE is None:
        i = np.arange(n) + 0.5
        phi = np.arccos(1 - 2 * i / n)
        theta = np.pi * (1 + 5 ** 0.5) * i
        _SPHERE = np.stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
    return _SPHERE


def psea_ss(ca_by_num):
    """CA-only secondary structure, P-SEA style. Approximate; no DSSP here.

    Helix and strand are called from the CA(i)-CA(i+3) and CA(i)-CA(i+4)
    distances, which separate the two cleanly in ideal geometry (helix ~5.3 and
    ~6.4 A, strand ~9.9 and ~12.4 A). Everything else is coil.
    """
    nums = sorted(ca_by_num)
    xyz = np.array([ca_by_num[n] for n in nums])
    ss = {n: "C" for n in nums}
    for k, n in enumerate(nums):
        # only call on a contiguous stretch of the renumbered chain
        if k + 4 >= len(nums) or nums[k + 4] != n + 4:
            continue
        d3 = np.linalg.norm(xyz[k] - xyz[k + 3])
        d4 = np.linalg.norm(xyz[k] - xyz[k + 4])
        if 4.6 <= d3 <= 6.2 and 5.0 <= d4 <= 7.6:
            for j in range(k, k + 5):
                ss[nums[j]] = "H"
        elif 8.8 <= d3 <= 11.0 and 11.0 <= d4 <= 13.8:
            for j in range(k, k + 5):
                if ss[nums[j]] == "C":
                    ss[nums[j]] = "E"
    return ss


def receptor_descriptors(pdb_path, pdb_id):
    """Everything about the receptor that the pocket features are built from."""
    from Bio.PDB.SASA import ShrakeRupley

    model = I.load_structure(pdb_path, pdb_id)[0]
    residues = [r for r in model.get_residues() if I.is_polymer_residue(r)]
    ShrakeRupley().compute(model, level="R")

    heavy, ca, sasa, rsa, name = [], {}, {}, {}, {}
    for r in residues:
        num = r.id[1]
        nm = r.get_resname().strip().upper()
        name[num] = nm
        sasa[num] = float(getattr(r, "sasa", 0.0))
        rsa[num] = sasa[num] / MAX_ASA.get(nm, 200.0)
        for a in I.heavy_atoms(r):
            heavy.append(a.get_coord())
        if "CA" in r:
            ca[num] = np.asarray(r["CA"].get_coord(), dtype=float)
    heavy = np.asarray(heavy, dtype=float)
    centroid = heavy.mean(axis=0)
    rg = float(np.sqrt(((heavy - centroid) ** 2).sum(axis=1).mean()))
    return {"heavy": heavy, "centroid": centroid, "rg": rg, "ca": ca,
            "sasa": sasa, "rsa": rsa, "name": name, "ss": psea_ss(ca),
            "n_residues": len(residues)}


def pocket_features(pk, rec, res_prob, n_offered, score_sum):
    """Truth-free descriptors of ONE candidate pocket."""
    lining = sorted(pocket_residues(pk))
    f = {
        "p2rank_rank": float(pk["rank"]),
        "p2rank_score": float(pk["score"]),
        "p2rank_probability": float(pk["probability"]),
        "p2rank_score_share": float(pk["score"] / score_sum) if score_sum else 0.0,
        "p2rank_sas_points": float(pk["sas_points"]),
        "p2rank_surf_atoms": float(pk["surf_atoms"]),
        "n_pockets_offered": float(n_offered),
        "n_lining_residues": float(len(lining)),
        "lining_fraction_of_receptor": len(lining) / max(1, rec["n_residues"]),
    }
    center = np.asarray(pk["center"], dtype=float)

    # position
    d = float(np.linalg.norm(center - rec["centroid"]))
    f["center_dist_to_centroid"] = d
    f["center_dist_over_rg"] = d / rec["rg"] if rec["rg"] else 0.0

    # burial and enclosure, from the pocket centre
    delta = rec["heavy"] - center
    dist = np.linalg.norm(delta, axis=1)
    f["heavy_atoms_within_8a"] = float((dist <= 8.0).sum())
    f["heavy_atoms_within_12a"] = float((dist <= 12.0).sum())
    f["burial_ratio_8_over_12"] = (f["heavy_atoms_within_8a"] /
                                   max(1.0, f["heavy_atoms_within_12a"]))
    near = delta[dist <= 15.0]
    if len(near):
        unit = near / np.linalg.norm(near, axis=1, keepdims=True)
        cos = sphere_directions() @ unit.T
        f["enclosure_fraction"] = float((cos.max(axis=1) > 0.95).mean())
    else:
        f["enclosure_fraction"] = 0.0

    # lining-residue properties
    if lining:
        cas = np.array([rec["ca"][n] for n in lining if n in rec["ca"]])
        f["mean_rsa"] = float(np.mean([rec["rsa"].get(n, 0.0) for n in lining]))
        f["sum_sasa"] = float(np.sum([rec["sasa"].get(n, 0.0) for n in lining]))
        f["p2rank_residue_prob_sum"] = float(np.sum([res_prob.get(n, 0.0)
                                                     for n in lining]))
        f["p2rank_residue_prob_mean"] = float(np.mean([res_prob.get(n, 0.0)
                                                       for n in lining]))
        nm = [rec["name"].get(n, "") for n in lining]
        f["hydrophobic_fraction"] = float(np.mean([x in HYDROPHOBIC for x in nm]))
        f["aromatic_fraction"] = float(np.mean([x in AROMATIC for x in nm]))
        f["charged_fraction"] = float(np.mean([x in CHARGED for x in nm]))
        f["polar_fraction"] = float(np.mean([x in POLAR for x in nm]))
        ss = [rec["ss"].get(n, "C") for n in lining]
        f["helix_fraction"] = float(np.mean([x == "H" for x in ss]))
        f["strand_fraction"] = float(np.mean([x == "E" for x in ss]))
        f["coil_fraction"] = float(np.mean([x == "C" for x in ss]))
        runs = 1 + sum(1 for a, b in zip(lining, lining[1:]) if b - a > 1)
        f["n_sequence_runs"] = float(runs)
        f["mean_run_length"] = len(lining) / runs
        if len(cas) >= 3:
            c = cas - cas.mean(axis=0)
            ev = np.sort(np.linalg.eigvalsh(np.cov(c.T)))[::-1]
            ev = np.maximum(ev, 1e-9)
            f["elongation_l1_over_l2"] = float(ev[0] / ev[1])
            f["anisotropy_l1_share"] = float(ev[0] / ev.sum())
            f["max_ca_extent"] = float(np.max(
                np.linalg.norm(cas[:, None] - cas[None], axis=-1)))
        else:
            f["elongation_l1_over_l2"] = 1.0
            f["anisotropy_l1_share"] = 1 / 3
            f["max_ca_extent"] = 0.0
    else:
        for k in ("mean_rsa", "sum_sasa", "p2rank_residue_prob_sum",
                  "p2rank_residue_prob_mean", "hydrophobic_fraction",
                  "aromatic_fraction", "charged_fraction", "polar_fraction",
                  "helix_fraction", "strand_fraction", "coil_fraction",
                  "n_sequence_runs", "mean_run_length", "max_ca_extent"):
            f[k] = 0.0
        f["elongation_l1_over_l2"] = 1.0
        f["anisotropy_l1_share"] = 1 / 3
    return f, lining


FEATURES = None  # filled from the first pocket built, so the list cannot drift


def build_features(prepared, refresh=False):
    """One row per (complex, offered pocket), with truth attached for scoring only."""
    global FEATURES
    if FEATURE_CACHE.exists() and not refresh:
        d = json.loads(FEATURE_CACHE.read_text())
        FEATURES = d["feature_names"]
        return d
    t0 = time.time()
    cases, bad_tokens, missing = [], 0, []
    for i, r in enumerate(prepared):
        pid = r["pdb_id"]
        pred = P2RANK_OUT / f"{pid}.pdb_predictions.csv"
        if not pred.exists():
            missing.append(pid)
            continue
        pks, bad = read_predictions(pred)
        bad_tokens += bad
        res_prob = read_residue_scores(P2RANK_OUT / f"{pid}.pdb_residues.csv")
        rec = receptor_descriptors(Path(r["receptor_pdb"]), pid)
        ssum = sum(p["score"] for p in pks)
        truth = set(r["interface_residues"])
        pockets = []
        for pk in pks:
            f, lining = pocket_features(pk, rec, res_prob, len(pks), ssum)
            if FEATURES is None:
                FEATURES = sorted(f)
            pockets.append({"rank": pk["rank"], "features": f,
                            "residues": lining, **score(set(lining), truth)})
        cases.append({
            "pdb_id": pid, "klass": r["klass"],
            "n_receptor_residues": rec["n_residues"],
            "n_interface_residues": len(truth),
            "n_pockets_offered": len(pks), "pockets": pockets,
        })
        if (i + 1) % 50 == 0:
            print(f"  features {i + 1}/{len(prepared)}", flush=True)
    d = {"n_cases": len(cases), "missing_predictions": missing,
         "n_noninteger_residue_tokens": bad_tokens,
         "feature_names": FEATURES, "wall_seconds": round(time.time() - t0, 1),
         "cases": cases}
    FEATURE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    FEATURE_CACHE.write_text(json.dumps(d))
    print(f"[features] {len(cases)} complexes, "
          f"{sum(len(c['pockets']) for c in cases)} pockets, "
          f"{d['wall_seconds']:.0f}s")
    return d


# --------------------------------------------------------------------------
# 4. the rule family - declared in full before any of it is scored
# --------------------------------------------------------------------------
def zstats(cases, names):
    """Mean and sd of each feature over the DEVELOPMENT pockets only."""
    X = np.array([[p["features"][n] for n in names]
                  for c in cases for p in c["pockets"]], dtype=float)
    if not len(X):
        return np.zeros(len(names)), np.ones(len(names))
    mu, sd = X.mean(axis=0), X.std(axis=0)
    return mu, np.where(sd > 1e-9, sd, 1.0)


def make_rules(names):
    """Every candidate rule. Each maps a pocket's features to a score; the
    highest-scoring pocket in a complex is the one picked."""
    rules = {"baseline_p2rank_top1": ("single", "p2rank_rank", -1.0)}
    for n in names:
        if n == "p2rank_rank":
            continue
        rules[f"high_{n}"] = ("single", n, +1.0)
        rules[f"low_{n}"] = ("single", n, -1.0)
    # a small number of declared composites, from the groove hypothesis: a
    # peptide site is long, shallow and less enclosed than the deep pocket
    # P2Rank prefers, but P2Rank's own confidence still carries information.
    rules["composite_score_x_extent"] = ("z", {"p2rank_score": 1.0,
                                               "max_ca_extent": 1.0})
    rules["composite_score_x_elongation"] = ("z", {"p2rank_score": 1.0,
                                                   "anisotropy_l1_share": 1.0})
    rules["composite_prob_minus_enclosure"] = ("z", {"p2rank_probability": 1.0,
                                                     "enclosure_fraction": -1.0})
    rules["composite_size_x_exposure"] = ("z", {"n_lining_residues": 1.0,
                                                "mean_rsa": 1.0})
    rules["learned_logistic"] = ("logistic", None)
    return rules


def rule_scores(rule, case, names, mu, sd, model):
    kind = rule[0]
    if kind == "single":
        _, name, sign = rule
        return [sign * p["features"][name] for p in case["pockets"]]
    X = np.array([[p["features"][n] for n in names] for p in case["pockets"]],
                 dtype=float)
    if not len(X):
        return []
    Z = (X - mu) / sd
    if kind == "z":
        w = np.array([rule[1].get(n, 0.0) for n in names])
        return list(Z @ w)
    return list(model.decision_function(Z))


def apply_rule(rule, cases, names, mu, sd, model):
    """Pick one pocket per complex; return the per-complex Jaccard and the picks."""
    js, picks = [], []
    for c in cases:
        if not c["pockets"]:
            js.append(0.0)
            picks.append(None)
            continue
        s = rule_scores(rule, c, names, mu, sd, model)
        # ties broken by P2Rank rank, so a rule with no signal degrades exactly
        # to the baseline rather than to an arbitrary pocket
        k = min(range(len(s)), key=lambda i: (-s[i], c["pockets"][i]["rank"]))
        js.append(c["pockets"][k]["jaccard"])
        picks.append(c["pockets"][k]["rank"])
    return np.array(js), picks


def fit_logistic(cases, names, mu, sd):
    """Label a DEVELOPMENT pocket 1 if it is the best pocket its complex offers
    and it actually touches the interface. Trained on development only."""
    from sklearn.linear_model import LogisticRegression

    X, y = [], []
    for c in cases:
        if not c["pockets"]:
            continue
        best = max(p["jaccard"] for p in c["pockets"])
        for p in c["pockets"]:
            X.append([p["features"][n] for n in names])
            y.append(1 if (p["jaccard"] == best and p["jaccard"] > 0) else 0)
    X = np.array(X, dtype=float)
    y = np.array(y)
    if len(set(y)) < 2:
        return None
    Z = (X - mu) / sd
    m = LogisticRegression(C=0.1, max_iter=5000, random_state=SEED)
    m.fit(Z, y)
    return m


# --------------------------------------------------------------------------
# 5. the split
# --------------------------------------------------------------------------
def fold_of(pdb_id):
    """Deterministic, outcome-blind half: parity of the md5 of the PDB ID."""
    return int(hashlib.md5(pdb_id.encode()).hexdigest(), 16) % 2


def paired(a, b, label_a, label_b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    pb = M.paired_bootstrap(a, b, n_boot=N_BOOT, seed=SEED)
    w = M.wilcoxon(a, b)
    return {
        "arm_a": label_a, "arm_b": label_b, "n": int(len(a)),
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
        "delta": float(a.mean() - b.mean()),
        "ci_lo": float(pb["ci_lo"]), "ci_hi": float(pb["ci_hi"]),
        "p_wilcoxon": float(w["p_value"]),
        "wilcoxon_statistic": float(w["statistic"]),
        "p_paired_bootstrap": float(pb["p_value"]),
        "n_pairs_entering_wilcoxon": int(w["n_effective"]),
        "n_a_better": int((a > b).sum()), "n_b_better": int((b > a).sum()),
        "n_tied": int((a == b).sum()),
    }


def arms(cases, rule, names, mu, sd, model):
    """Every arm on one set of complexes, on the same pairing."""
    top1 = np.array([c["pockets"][0]["jaccard"] if c["pockets"] else 0.0
                     for c in cases])
    oracle = np.array([max((p["jaccard"] for p in c["pockets"]), default=0.0)
                       for c in cases])
    rr, picks = apply_rule(rule, cases, names, mu, sd, model)
    return top1, rr, oracle, picks


def block(cases, rule, rule_name, names, mu, sd, model):
    top1, rr, oracle, picks = arms(cases, rule, names, mu, sd, model)
    gap = float(oracle.mean() - top1.mean())
    closed = float(rr.mean() - top1.mean())
    return {
        "n": len(cases),
        "rule": rule_name,
        "p2rank_top1_jaccard": float(top1.mean()),
        "reranked_jaccard": float(rr.mean()),
        "ORACLE_bestany_jaccard": float(oracle.mean()),
        "oracle_note": ("ORACLE - the pocket is chosen by the truth. It is the "
                        "ceiling for what re-ranking could ever reach on these "
                        "complexes. It is not a method and not achievable."),
        "headroom_oracle_minus_top1": gap,
        "headroom_closed": closed,
        "fraction_of_headroom_closed": closed / gap if gap > 1e-12 else float("nan"),
        "reranked_vs_top1": paired(rr, top1, "reranked", "p2rank_top1"),
        "reranked_ci": [float(x) for x in
                        M.bootstrap_ci(rr, n_boot=N_BOOT, seed=SEED)[1:]],
        "p2rank_top1_ci": [float(x) for x in
                           M.bootstrap_ci(top1, n_boot=N_BOOT, seed=SEED)[1:]],
        "n_complexes_with_no_pocket": int(sum(1 for c in cases
                                              if not c["pockets"])),
        "n_complexes_with_one_pocket": int(sum(1 for c in cases
                                               if len(c["pockets"]) == 1)),
        "n_complexes_rerankable": int(sum(1 for c in cases
                                          if len(c["pockets"]) > 1)),
        "n_picks_changed": int(sum(1 for c, p in zip(cases, picks)
                                   if p is not None and p != 1)),
        "n_picks_changed_with_no_effect_on_jaccard":
            int(sum(1 for c, p, a, b in zip(cases, picks, rr, top1)
                    if p is not None and p != 1 and a == b)),
        "n_complexes_where_no_offered_pocket_touches_the_site":
            int(sum(1 for c in cases
                    if not any(pk["jaccard"] > 0 for pk in c["pockets"]))),
        "mean_picked_rank": float(np.mean([p for p in picks if p is not None]))
        if any(p is not None for p in picks) else float("nan"),
    }


# --------------------------------------------------------------------------
# 6. mechanism - can any descriptor tell the right pocket from the wrong one?
# --------------------------------------------------------------------------
def within_complex_auc(cases, names):
    """Pooled WITHIN-complex pairwise AUC for each feature.

    For every complex offering two or more pockets, every (best pocket, other
    pocket) pair contributes: the feature scores a win if it ranks the best
    pocket higher, half a win on a tie. 0.5 is chance. This asks the question
    the re-ranker needs answered - inside one protein, does this descriptor
    separate the pocket that is the epitope from the ones that are not - and it
    is free of the between-protein variance that a pooled AUC would confound.

    Computed on ALL complexes of the class, so it is a diagnostic and NOT a
    held-out number. It is here to say why the re-ranker does what it does.
    """
    usable = []
    for c in cases:
        pk = c["pockets"]
        if len(pk) < 2:
            continue
        best = max(p["jaccard"] for p in pk)
        if best <= 0:
            continue                         # nothing offered touches the site
        pos = [p for p in pk if p["jaccard"] == best]
        neg = [p for p in pk if p["jaccard"] < best]
        if pos and neg:
            usable.append((pos, neg))

    def tally(subset, n):
        wins = total = 0.0
        for pos, neg in subset:
            for a in pos:
                for b in neg:
                    va, vb = a["features"][n], b["features"][n]
                    wins += 1.0 if va > vb else (0.5 if va == vb else 0.0)
                    total += 1
        return (wins / total, int(total)) if total else (float("nan"), 0)

    # the CI resamples COMPLEXES, not pairs: pairs inside one protein are not
    # independent draws and a pair-level interval would be far too narrow
    rng = np.random.default_rng(SEED)
    idx = [rng.integers(0, len(usable), len(usable)) for _ in range(2000)] \
        if usable else []
    out = {}
    for n in names:
        auc, npairs = tally(usable, n)
        boots = [tally([usable[i] for i in ix], n)[0] for ix in idx]
        boots = [b for b in boots if np.isfinite(b)]
        out[n] = {"auc": float(auc), "n_pairs": npairs,
                  "n_complexes": len(usable),
                  "ci_lo": float(np.percentile(boots, 2.5)) if boots else float("nan"),
                  "ci_hi": float(np.percentile(boots, 97.5)) if boots else float("nan")}
    return out


def headroom_by_pocket_count(cases):
    """Where the oracle headroom actually sits, by how many pockets were offered."""
    buckets = {"0_pockets": [], "1_pocket": [], "2_or_more": []}
    for c in cases:
        k = ("0_pockets" if not c["pockets"]
             else "1_pocket" if len(c["pockets"]) == 1 else "2_or_more")
        buckets[k].append(c)
    out = {}
    for k, cs in buckets.items():
        if not cs:
            out[k] = {"n": 0}
            continue
        t = np.array([c["pockets"][0]["jaccard"] if c["pockets"] else 0.0 for c in cs])
        o = np.array([max((p["jaccard"] for p in c["pockets"]), default=0.0) for c in cs])
        out[k] = {"n": len(cs), "p2rank_top1_jaccard": float(t.mean()),
                  "ORACLE_bestany_jaccard": float(o.mean()),
                  "headroom": float(o.mean() - t.mean()),
                  "note": ("no pocket to re-rank - top-1 and the oracle are the "
                           "same object here" if k != "2_or_more" else
                           "the only complexes a re-ranker can act on")}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    prep_path = CACHE / "prepared.json"
    if not prep_path.exists():
        sys.exit(f"missing {prep_path}; run scripts/epitope_gate.py first")
    prepared = json.loads(prep_path.read_text())["prepared"]

    feats = build_features(prepared, refresh=args.refresh)
    names = feats["feature_names"]
    by_class = {}
    for c in feats["cases"]:
        by_class.setdefault(c["klass"], []).append(c)
    pep = sorted(by_class.get("peptide", []), key=lambda c: c["pdb_id"])

    # ---- split, declared before anything is scored -----------------------
    dev = [c for c in pep if fold_of(c["pdb_id"]) == 0]
    test = [c for c in pep if fold_of(c["pdb_id"]) == 1]
    print(f"\npeptide n={len(pep)}  development n={len(dev)}  held-out n={len(test)}")

    mu, sd = zstats(dev, names)
    model = fit_logistic(dev, names, mu, sd)
    rules = make_rules(names)
    if model is None:
        rules.pop("learned_logistic", None)

    # ---- choose on development ONLY --------------------------------------
    dev_top1 = np.array([c["pockets"][0]["jaccard"] if c["pockets"] else 0.0
                         for c in dev])
    leaderboard = []
    for rn, rule in rules.items():
        j, _ = apply_rule(rule, dev, names, mu, sd, model)
        leaderboard.append({"rule": rn, "dev_jaccard": float(j.mean()),
                            "dev_delta_vs_top1": float(j.mean() - dev_top1.mean())})
    leaderboard.sort(key=lambda r: -r["dev_jaccard"])
    chosen_name = leaderboard[0]["rule"]
    chosen = rules[chosen_name]
    print(f"\n{len(rules)} candidate rules searched on the development half.")
    print("top 8 on development (these numbers are OPTIMISTIC - the rule was "
          "chosen on this data):")
    for r in leaderboard[:8]:
        print(f"  {r['rule']:<38} dev jaccard {r['dev_jaccard']:.4f} "
              f"({r['dev_delta_vs_top1']:+.4f})")
    print(f"\nchosen: {chosen_name}")

    # ---- report on the held-out half -------------------------------------
    dev_block = block(dev, chosen, chosen_name, names, mu, sd, model)
    test_block = block(test, chosen, chosen_name, names, mu, sd, model)
    full_block = block(pep, chosen, chosen_name, names, mu, sd, model)

    # ---- the swap, as a robustness check ---------------------------------
    mu2, sd2 = zstats(test, names)
    model2 = fit_logistic(test, names, mu2, sd2)
    rules2 = make_rules(names)
    if model2 is None:
        rules2.pop("learned_logistic", None)
    t_top1 = np.array([c["pockets"][0]["jaccard"] if c["pockets"] else 0.0
                       for c in test])
    lb2 = []
    for rn, rule in rules2.items():
        j, _ = apply_rule(rule, test, names, mu2, sd2, model2)
        lb2.append({"rule": rn, "dev_jaccard": float(j.mean()),
                    "dev_delta_vs_top1": float(j.mean() - t_top1.mean())})
    lb2.sort(key=lambda r: -r["dev_jaccard"])
    chosen2_name = lb2[0]["rule"]
    swap_block = block(dev, rules2[chosen2_name], chosen2_name, names,
                       mu2, sd2, model2)

    # ---- the locked rule on the other classes ----------------------------
    others = {}
    for k in ("ppi", "antibody"):
        cs = sorted(by_class.get(k, []), key=lambda c: c["pdb_id"])
        if cs:
            others[k] = block(cs, chosen, chosen_name, names, mu, sd, model)

    # ---- how much can this feature set express AT ALL? --------------------
    # Fit and scored on the same complexes, so it is an upper bound and is
    # labelled one. If even this cannot close the gap, the descriptors are the
    # limit and no amount of honest fitting will change that.
    in_sample = {}
    for k, cs in (("peptide", pep), ("ppi", by_class.get("ppi", [])),
                  ("antibody", by_class.get("antibody", []))):
        if not cs:
            continue
        mu_a, sd_a = zstats(cs, names)
        m_a = fit_logistic(cs, names, mu_a, sd_a)
        if m_a is None:
            in_sample[k] = {"n": len(cs), "ran": False,
                            "why": "no complex offers a pocket that touches the site"}
            continue
        b = block(cs, ("logistic", None), "learned_logistic_FIT_ON_THIS_DATA",
                  names, mu_a, sd_a, m_a)
        in_sample[k] = {kk: b[kk] for kk in
                        ("n", "p2rank_top1_jaccard", "reranked_jaccard",
                         "ORACLE_bestany_jaccard", "headroom_oracle_minus_top1",
                         "headroom_closed", "fraction_of_headroom_closed",
                         "n_picks_changed")}

    # ---- PASS gate --------------------------------------------------------
    need = test_block["p2rank_top1_jaccard"] + 0.5 * test_block["headroom_oracle_minus_top1"]
    passed = (test_block["reranked_jaccard"] >= need and
              test_block["reranked_vs_top1"]["p_wilcoxon"] < 0.05)

    # ---- Holm over the family of reported comparisons ---------------------
    fam = [("peptide_heldout", test_block["reranked_vs_top1"]["p_wilcoxon"]),
           ("peptide_development", dev_block["reranked_vs_top1"]["p_wilcoxon"]),
           ("ppi", others.get("ppi", {}).get("reranked_vs_top1", {}).get("p_wilcoxon")),
           ("antibody", others.get("antibody", {}).get("reranked_vs_top1", {}).get("p_wilcoxon"))]
    fam = [(k, p) for k, p in fam if p is not None]
    adj = M.holm_bonferroni(dict(fam))
    holm = {k: {"p_raw": float(p), "p_holm": float(adj[k]),
                "reject_at_0.05": bool(adj[k] < 0.05)} for k, p in fam}

    # ---- what could not be evaluated, and why -----------------------------
    pesto = sorted(ROOT.glob("results/*pesto*")) + sorted(ROOT.glob("results/*PeSTo*"))
    not_evaluated = [
        {"feature": "PeSTo PPI-head score mass inside the pocket",
         "why": (("this feature set was frozen before E1 existed and does not "
                  "use PeSTo. Verified: "
                  f"{len(pesto)} file(s) match results/*pesto* at run time, so "
                  "a PeSTo output is present now but is NOT read here.")
                 if pesto else
                 ("E1 has not produced a PeSTo output. Verified: "
                  f"{len(pesto)} files match results/*pesto* at run time.")),
         "checked": [str(x.relative_to(ROOT)) for x in pesto]},
        {"feature": "sequence conservation of the lining residues",
         "why": ("needs an MSA per receptor (192 distinct UniProts). No local "
                 "database and no network budget inside this experiment, so it "
                 "is not cheap here and was not attempted.")},
        {"feature": "DSSP secondary structure",
         "why": ("no mkdssp/dssp binary on this machine (verified with "
                 "shutil.which). Helix/strand/coil are a CA-only P-SEA style "
                 "approximation instead, and are labelled as such.")},
        {"feature": "re-running P2Rank",
         "why": ("not needed - the cached predictions under "
                 "data/raw/epitope_gate/p2rank/out are reused unchanged, so "
                 "top-1 here reproduces epitope_gate.json exactly.")},
    ]

    if passed:
        verdict = (
            f"RE-RANK. A rule chosen on {len(dev)} development peptide complexes "
            f"lifts held-out Jaccard from {test_block['p2rank_top1_jaccard']:.4f} "
            f"to {test_block['reranked_jaccard']:.4f} on the {len(test)} held out "
            f"(p={test_block['reranked_vs_top1']['p_wilcoxon']:.3g}), closing "
            f"{test_block['fraction_of_headroom_closed']:.0%} of the oracle headroom.")
    else:
        best_auc = max(
            ((n, v) for n, v in within_complex_auc(pep, names).items()),
            key=lambda kv: abs(kv[1]["auc"] - 0.5))
        verdict = (
            f"DO NOT RE-RANK - the headroom is real but these descriptors cannot "
            f"reach it. The rule chosen on {len(dev)} development complexes moves "
            f"held-out peptide Jaccard {test_block['headroom_closed']:+.4f} "
            f"({test_block['p2rank_top1_jaccard']:.4f} -> "
            f"{test_block['reranked_jaccard']:.4f}, n={len(test)}, Wilcoxon "
            f"p={test_block['reranked_vs_top1']['p_wilcoxon']:.3g}), which is "
            f"{test_block['fraction_of_headroom_closed']:.0%} of the "
            f"{test_block['headroom_oracle_minus_top1']:.4f} oracle headroom on "
            f"those same complexes, against a {need - test_block['p2rank_top1_jaccard']:+.4f} "
            f"bar. Running the folds the other way round moves it "
            f"{swap_block['headroom_closed']:+.4f}, i.e. the wrong way. The "
            f"mechanism is measurable: inside one protein the best descriptor of "
            f"the {len(names)} tried separates the right pocket from the wrong "
            f"ones at AUC {best_auc[1]['auc']:.3f} "
            f"[{best_auc[1]['ci_lo']:.3f}, {best_auc[1]['ci_hi']:.3f}] "
            f"({best_auc[0]}, {best_auc[1]['n_pairs']} pairs over "
            f"{best_auc[1]['n_complexes']} complexes), and even a rule FIT AND "
            f"SCORED on all 70 peptides - an upper bound, not a result - closes "
            f"only {in_sample['peptide']['fraction_of_headroom_closed']:.0%}. "
            f"The +0.128 gap is a fact about P2Rank's ranking; it is not "
            f"recoverable from pocket geometry, P2Rank's own scores, or "
            f"lining-residue composition.")

    # ---- does top-1 here reproduce epitope_gate.json exactly? -------------
    eg = ROOT / "results" / "epitope_gate.json"
    repro = {"checked": False, "why": f"{eg} not present"}
    if eg.exists():
        prev = {r["pdb_id"]: r for r in json.loads(eg.read_text())["per_structure"]}
        diffs = []
        for c in feats["cases"]:
            a = c["pockets"][0]["jaccard"] if c["pockets"] else 0.0
            b = prev.get(c["pdb_id"], {}).get("p2rank_top1", {}).get("jaccard")
            if b is None or abs(a - b) > 1e-9:
                diffs.append({"pdb_id": c["pdb_id"], "here": a, "epitope_gate": b})
        repro = {"checked": True, "n_compared": len(feats["cases"]),
                 "n_disagreeing": len(diffs), "disagreements": diffs[:10],
                 "note": ("p2rank_top1 is re-derived here from the same cached "
                          "CSVs; a non-zero disagreement would mean the two "
                          "scripts are not measuring the same thing")}

    out = {
        "experiment": "e2 - re-rank P2Rank's own pockets to pick the peptide site",
        "question": ("P2Rank offers a pocket touching the true interface in 78.6% "
                     "of peptide complexes but ranks the wrong one first. Can the "
                     "right one be picked without the answer?"),
        "ground_truth": ("receptor residues within 4.5 A (heavy atom) of any binder "
                         "heavy atom, read from data/raw/epitope_gate/prepared.json "
                         "as written by scripts/epitope_gate.py - not recomputed here"),
        "scoring": "Jaccard over residue-number sets, identical to m2_gate.score",
        "verdict": verdict,
        "pass_criterion": {
            "stated": ("held-out peptide Jaccard >= P2Rank top-1 + half the "
                       "oracle headroom on the SAME held-out complexes, and "
                       "paired Wilcoxon p < 0.05"),
            "threshold_jaccard": float(need),
            "achieved_jaccard": test_block["reranked_jaccard"],
            "p_wilcoxon": test_block["reranked_vs_top1"]["p_wilcoxon"],
            "passed": bool(passed),
        },
        "split": {
            "rule": "md5(pdb_id) parity, fixed before any scoring",
            "n_peptide": len(pep), "n_development": len(dev), "n_held_out": len(test),
            "development_pdb_ids": [c["pdb_id"] for c in dev],
            "held_out_pdb_ids": [c["pdb_id"] for c in test],
        },
        "rule_search": {
            "n_candidate_rules": len(rules),
            "chosen": chosen_name,
            "chosen_on": "development half only",
            "note": ("development numbers are optimistic by construction: the "
                     "rule is the winner of a search over this data. The held-out "
                     "block is the result."),
            "leaderboard_development": leaderboard,
        },
        "peptide_held_out": test_block,
        "peptide_development": dev_block,
        "peptide_all_70_not_a_held_out_number": full_block,
        "swapped_folds": {
            "note": ("the same procedure run the other way round - chosen on the "
                     "held-out half, reported on the development half. A "
                     "robustness check, not a second independent result."),
            "chosen": chosen2_name,
            **swap_block,
        },
        "other_classes_same_locked_rule": others,
        "holm_bonferroni": holm,
        "mechanism_within_complex_discrimination": {
            "what": ("pooled within-complex pairwise AUC of each descriptor at "
                     "ranking the best offered pocket above a worse one. 0.5 is "
                     "chance."),
            "not_held_out": ("computed on every complex of the class; a "
                             "diagnostic of why the re-ranker behaves as it "
                             "does, not a result"),
            "peptide": within_complex_auc(pep, names),
            "ppi": within_complex_auc(by_class.get("ppi", []), names),
            "antibody": within_complex_auc(by_class.get("antibody", []), names),
        },
        "in_sample_ceiling_of_this_feature_set": {
            "what": ("the learned rule FIT AND SCORED on the same complexes. "
                     "Maximally optimistic; it is an upper bound on what these "
                     f"{len(names)} descriptors can express, not a result."),
            **in_sample,
        },
        "headroom_by_pockets_offered": {
            "peptide": headroom_by_pocket_count(pep),
            "ppi": headroom_by_pocket_count(by_class.get("ppi", [])),
            "antibody": headroom_by_pocket_count(by_class.get("antibody", [])),
        },
        "features": {
            "n": len(names), "names": names,
            "none_can_see_the_binder": True,
            "secondary_structure": ("CA-only P-SEA style approximation; no DSSP "
                                    "binary on this machine"),
        },
        "reproduces_epitope_gate_top1": repro,
        "p2rank": {"version": "2.5", "predictions_reused_from": str(P2RANK_OUT),
                   "n_complexes": feats["n_cases"],
                   "missing_predictions": feats["missing_predictions"],
                   "n_noninteger_residue_tokens": feats["n_noninteger_residue_tokens"]},
        "credits_spent": {"rowan": 0, "note": "local CPU only"},
        "not_evaluated": not_evaluated,
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))

    # ---- console ---------------------------------------------------------
    def line(tag, b):
        print(f"  {tag:<26}{b['n']:>5}{b['p2rank_top1_jaccard']:>10.4f}"
              f"{b['reranked_jaccard']:>11.4f}{b['ORACLE_bestany_jaccard']:>10.4f}"
              f"{b['headroom_closed']:>+10.4f}"
              f"{b['reranked_vs_top1']['p_wilcoxon']:>10.3g}")

    print("\n" + "=" * 86)
    print(f"  {'set':<26}{'n':>5}{'top-1':>10}{'reranked':>11}"
          f"{'ORACLE':>10}{'delta':>10}{'p':>10}")
    print("-" * 86)
    line("peptide HELD-OUT", test_block)
    line("peptide development", dev_block)
    line("peptide all 70", full_block)
    line("peptide swapped folds", out["swapped_folds"])
    for k, b in others.items():
        line(k, b)
    print("-" * 86)
    print(f"  ORACLE is the best pocket at ANY rank, chosen BY THE TRUTH. It is a "
          f"ceiling,\n  not a method, and it is not achievable.")
    print(f"\n  PASS needed held-out Jaccard >= {need:.4f} "
          f"(top-1 {test_block['p2rank_top1_jaccard']:.4f} + half the "
          f"{test_block['headroom_oracle_minus_top1']:.4f} oracle headroom) "
          f"at p < 0.05")
    print(f"  got {test_block['reranked_jaccard']:.4f} at p="
          f"{test_block['reranked_vs_top1']['p_wilcoxon']:.3g}  ->  "
          f"{'PASS' if passed else 'FAIL'}")
    print(f"\n  {verdict}")
    print(f"\nwrote {OUT}")
    return out


if __name__ == "__main__":
    main()
