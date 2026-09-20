"""E1 - can ONE model replace the site-definition stage for every modality?

P2Rank 2.5 defines the binding site today, and results/epitope_gate.json measured
that it fails outside concave pockets: top-1 Jaccard 0.6355 on KDR small molecules
but 0.0255 on antibody epitopes, which is chance. The axis is concavity, not
modality. This script asks whether PeSTo can do the whole job instead.

PeSTo (Krapp, Gainza, Schwede, Correia, "PeSTo: parameter-free geometric deep
learning for accurate prediction of protein binding interfaces", Nat Commun 14,
2175, 2023; https://github.com/LBM-EPFL/PeSTo) is a structure-only transformer
that emits FIVE per-residue interface heads from one forward pass - protein,
nucleic acid, ion, ligand, lipid. Two of them are exactly the two site-definition
problems this project has: the protein head is the epitope/PPI problem, the ligand
head is the pocket problem. One model, one pass, both modalities.

LICENCE - PeSTo's weights and source are CC BY-NC-SA 4.0. This repository is MIT.
Nothing from PeSTo is vendored here. This script LOCATES a PeSTo checkout at
runtime (--pesto-root, or $PESTO_ROOT, or --fetch to git-clone one) and refuses
to run against a checkout that lives inside this repository, so no CC BY-NC-SA
file can be committed by accident. Only numbers computed from it are stored.

    ./env/bin/python scripts/pesto_site.py --fetch /tmp/pesto
    ./env/bin/python scripts/pesto_site.py --pesto-root /tmp/pesto/PeSTo

=============================================================================
PRE-REGISTRATION - frozen before the first forward pass on either benchmark
=============================================================================
Copied verbatim into results/pesto_vs_p2rank.json under "preregistration" so the
rule and the numbers cannot drift apart.

MODEL
  Release i_v4_1_2021-09-07_11-21, checkpoint model_ckpt.pt - the release and the
  file that PeSTo's own apply_model.ipynb selects. CPU, torch.no_grad, one
  forward pass per structure.

HEADS
  That release's config_data['r_types'] orders the five output channels
  [protein, dna+rna, ion, ligand, lipid]. The channel indices are not hard-coded:
  head_indices() re-derives them by matching each r_types entry against
  src.data_encoding.categ_to_resnames and asserts PPI == 0 and LIGAND == 3.
  p_i = sigmoid(z[i, head]), one probability per residue.

DECISION RULE (primary)
  1. S = { residues with p >= 0.50 }.
  2. Graph on S: residues u, v adjacent iff their minimum heavy-atom distance
     is <= 5.0 A.
  3. Predicted site = the largest connected component of that graph; ties broken
     by the component with the larger summed p.
  4. S empty -> the prediction is the empty set and scores 0 on every metric,
     which is how epitope_gate.py already treats a structure P2Rank offers no
     pocket for.

SCORING
  m2_gate.score() imported, not re-implemented - precision, recall, F1, Jaccard.
  Ground truth is the 4.5 A heavy-atom contact set ALREADY in this repo
  (scripts/interfaces.py wrote it; results/pipeline/*/target/pocket_recovery.json
  for small molecules, data/raw/epitope_gate/prepared.json for biologics). No
  third contact routine is written here.

FLOOR
  Per structure, 100 draws of |predicted site| residues uniformly from the
  receptor's observed residues, averaged. Same construction as epitope_gate.py's
  random arm, but size-matched to PeSTo instead of to P2Rank, so the floor tracks
  whatever PeSTo predicts.

STATISTICS - scripts/metrics.py only
  paired_bootstrap (20000, seed 0) and wilcoxon on the per-structure vectors;
  bootstrap_ci (10000, seed 0) on each arm's Jaccard; holm_bonferroni over ONE
  family containing every primary paired test in the run (the convention
  m2_gate.py uses). Everything is paired: every arm is scored on the same
  structures, and the P2Rank vectors are read back from the published result
  files rather than recomputed, so they are the same numbers as the headline.

PASS CRITERIA (set with the experiment, before the run)
  (a) small molecule - the LIGAND head passes if mean Jaccard reaches 0.6355
      (KDR) / 0.4938 (CDK2), or the paired CI on (PeSTo - P2Rank) contains 0.
  (b) biologic - the PPI head passes a class if it beats its size-matched random
      floor with Holm p < 0.05 on Jaccard.

LABELLED SWEEP (secondary; the list is fixed here so no cell can be picked after
the fact) threshold in {0.30, 0.50, 0.70, 0.90} x {largest component, no
component filter}, plus a top-k variant with k = |P2Rank top-1 pocket| for that
same structure. The 0.50 + largest-component cell is the primary and is the one
quoted anywhere outside the sweep block.
=============================================================================

ADDED AFTER THE FREEZE - two things, both recorded in the output JSON under
preregistration.added_after_the_freeze_and_why, neither of them a new primary:

  * pesto_bestcc_ORACLE, the best-scoring connected component. A 4-structure
    smoke test showed the failure mode the pre-registered rule has: on 7RKE the
    PPI head puts p >= 0.97 on the true epitope and also >= 0.5 on a second,
    LARGER protein-interface patch, so "largest component" returns the wrong
    blob and scores 0. The oracle arm separates "cannot see the site" from
    "sees it, picks the wrong one". It is the exact counterpart of
    p2rank_bestany, an oracle results/epitope_gate.json already publishes, and
    it is only ever compared against that oracle.
  * per-residue ROC AUC for BOTH methods, from PeSTo's probabilities and from
    P2Rank's own *_residues.csv score column. No threshold and no component
    rule on either side, so it measures the two scorers rather than this
    script's reducer. It is a diagnostic, not a pass criterion.

The pre-registered primary above is unchanged and is reported first everywhere.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M            # noqa: E402  bootstrap/wilcoxon/holm - not hand-rolled
import m2_gate as G            # noqa: E402  G.score is THE scoring function

ROOT = Path(__file__).resolve().parent.parent
SEED = 0

PESTO_REPO = "https://github.com/LBM-EPFL/PeSTo.git"
RELEASE = "i_v4_1_2021-09-07_11-21"
CHECKPOINT = "model_ckpt.pt"

# --- the pre-registered constants, in one place ---------------------------
THRESHOLD = 0.50
CONTACT_A = 5.0            # residue-residue adjacency for the component filter
N_FLOOR_DRAWS = 100
SWEEP_THRESHOLDS = (0.30, 0.50, 0.70, 0.90)
METRICS = ("precision", "recall", "f1", "jaccard")

SMALL_MOLECULE = (
    # label, results/pipeline dir, results/m2_gate_*.json, published P2Rank top-1 Jaccard
    ("KDR", "colorectal-cancer", "m2_gate_KDR.json", 0.6355),
    ("CDK2", "cdk2-second-target", "m2_gate_CDK2.json", 0.4938),
)


# --------------------------------------------------------------------------
# 0. locate PeSTo at runtime - never vendored, never committed
# --------------------------------------------------------------------------
def resolve_pesto(root=None, fetch=None):
    """Return a usable PeSTo checkout, cloning into `fetch` if asked.

    Refuses any checkout inside this repository: PeSTo is CC BY-NC-SA 4.0 and
    this repo is MIT, so the weights must never sit where git can stage them.
    """
    if fetch:
        dest = Path(fetch).expanduser().resolve()
        dest.mkdir(parents=True, exist_ok=True)
        root = dest / "PeSTo"
        if not (root / "src" / "structure_io.py").exists():
            print(f"[pesto] git clone {PESTO_REPO} -> {root}", flush=True)
            subprocess.run(["git", "clone", "--depth", "1", PESTO_REPO, str(root)],
                           check=True)
    if root is None:
        root = os.environ.get("PESTO_ROOT")
    if root is None:
        raise SystemExit(
            "PeSTo not located. Pass --pesto-root /path/to/PeSTo, set $PESTO_ROOT, "
            "or pass --fetch <dir> to clone it. It must live OUTSIDE this repo: "
            "its weights are CC BY-NC-SA 4.0 and this repo is MIT.")
    root = Path(root).expanduser().resolve()
    if ROOT == root or ROOT in root.parents:
        raise SystemExit(
            f"refusing to use a PeSTo checkout inside this repository ({root}). "
            "PeSTo is CC BY-NC-SA 4.0; this repo is MIT-licensed and must not "
            "contain it. Clone it elsewhere.")
    ckpt = root / "model" / "save" / RELEASE / CHECKPOINT
    if not ckpt.exists():
        raise SystemExit(f"no PeSTo checkpoint at {ckpt}")
    return root


def head_indices(pesto_root):
    """Re-derive the output channel of the PPI and ligand heads from the config.

    Verified rather than assumed: each r_types entry is matched back to the
    category in src.data_encoding.categ_to_resnames whose residue list it is.
    """
    _import_pesto(pesto_root)
    from src.data_encoding import categ_to_resnames        # noqa: PLC0415
    from config import config_data                         # noqa: PLC0415
    named = []
    for entry in config_data["r_types"]:
        s = set(entry)
        hit = [c for c, rn in categ_to_resnames.items() if set(rn) <= s]
        named.append("+".join(sorted(hit)))
    ppi = named.index("protein")
    lig = named.index("ligand")
    assert (ppi, lig) == (0, 3), f"unexpected head order {named}"
    return {"order": named, "ppi": ppi, "ligand": lig}


_PESTO_LOADED = None


def _import_pesto(pesto_root):
    """Put PeSTo on sys.path the way its own apply_model.ipynb does.

    Order matters and is not cosmetic: the checkout ships TWO copies of src/,
    the top-level one and an archived one under model/save/<release>/src, and
    they differ. apply_model.ipynb imports src.* first, from the top level, and
    only then puts the release directory on the path for config/model. Getting
    this backwards silently runs the archived parser. So the release directory
    goes on first, the repo root goes on top of it, and src.* is imported at once
    to pin it in sys.modules before anything else can claim the name.
    """
    global _PESTO_LOADED
    if _PESTO_LOADED == str(pesto_root):
        return
    for p in (str(Path(pesto_root) / "model" / "save" / RELEASE), str(pesto_root)):
        while p in sys.path:
            sys.path.remove(p)
        sys.path.insert(0, p)
    import src.data_encoding, src.structure, src.structure_io, src.dataset   # noqa: F401,PLC0415
    top = str(Path(pesto_root) / "src")
    got = str(Path(src.structure.__file__).parent)
    if got != top:
        raise RuntimeError(f"PeSTo src resolved to {got}, expected {top}")
    _PESTO_LOADED = str(pesto_root)


# --------------------------------------------------------------------------
# 1. one forward pass -> per-residue probabilities + a residue adjacency graph
# --------------------------------------------------------------------------
_MODEL = None


def _model(pesto_root):
    global _MODEL
    if _MODEL is None:
        import torch as pt                                  # noqa: PLC0415
        _import_pesto(pesto_root)
        from config import config_model                     # noqa: PLC0415
        from model import Model                             # noqa: PLC0415
        m = Model(config_model)
        m.load_state_dict(pt.load(
            Path(pesto_root) / "model" / "save" / RELEASE / CHECKPOINT,
            map_location="cpu"))
        _MODEL = m.eval()
    return _MODEL


def run_one(pesto_root, pdb_path):
    """PeSTo on one receptor PDB.

    Returns {"resnums": [author/renumbered residue number per output row],
             "probs": [[5 sigmoid probabilities] per residue],
             "edges": [[i, j] residue-index pairs within CONTACT_A heavy-atom],
             "seconds": wall time of the forward pass}.

    PeSTo's clean_structure() renumbers residues 1..N internally, so the original
    residue number is carried through as an extra per-atom key ('origid'). Every
    helper in PeSTo's src filters keys generically, so it survives untouched to
    the other end, and the mapping back is exact rather than positional.
    """
    import torch as pt                                      # noqa: PLC0415
    from scipy.spatial import cKDTree                       # noqa: PLC0415
    _import_pesto(pesto_root)
    from src.structure_io import read_pdb                                  # noqa: PLC0415
    from src.structure import (clean_structure, tag_hetatm_chains,         # noqa: PLC0415
                               split_by_chain, filter_non_atomic_subunits,
                               remove_duplicate_tagged_subunits,
                               concatenate_chains)
    from src.data_encoding import (encode_structure, encode_features,      # noqa: PLC0415
                                   extract_topology)
    from src.dataset import collate_batch_features                         # noqa: PLC0415

    st = read_pdb(str(pdb_path))
    chains = {c.split(":")[0] for c in st["chain_name"]}
    if len(chains) != 1:
        raise ValueError(f"{pdb_path}: expected one chain, found {sorted(chains)}")
    if any(i.strip() for i in st["icode"]):
        raise ValueError(f"{pdb_path}: insertion codes present, residue mapping unsafe")
    st["origid"] = np.array(st["resid"], dtype=np.int64)

    st = clean_structure(st)
    st = tag_hetatm_chains(st)
    sub = remove_duplicate_tagged_subunits(filter_non_atomic_subunits(split_by_chain(st)))
    st = concatenate_chains(sub)

    X, Mres = encode_structure(st)
    q = encode_features(st)[0]
    ids_topk, _, _, _, _ = extract_topology(X, 64)
    Xb, idsb, qb, Mb = collate_batch_features([[X, ids_topk, q, Mres]])
    t0 = time.time()
    with pt.no_grad():
        z = _model(pesto_root)(Xb, idsb, qb, Mb.float())
    secs = time.time() - t0
    probs = pt.sigmoid(z).numpy().astype(float)

    ures = np.unique(st["resid"])
    if probs.shape[0] != ures.size:
        raise ValueError(f"{pdb_path}: {probs.shape[0]} outputs for {ures.size} residues")
    resnums, xyz_by_res = [], []
    for r in ures:
        m = st["resid"] == r
        orig = np.unique(st["origid"][m])
        if orig.size != 1:
            raise ValueError(f"{pdb_path}: residue {r} maps to {orig.size} originals")
        resnums.append(int(orig[0]))
        xyz_by_res.append(st["xyz"][m])
    if len(set(resnums)) != len(resnums):
        raise ValueError(f"{pdb_path}: original residue numbers are not unique")

    # residue adjacency: minimum heavy-atom distance <= CONTACT_A. Computed once
    # here because it does not depend on the threshold, so the sweep reuses it.
    flat = np.concatenate(xyz_by_res)
    owner = np.concatenate([np.full(len(a), i) for i, a in enumerate(xyz_by_res)])
    pairs = cKDTree(flat).query_pairs(CONTACT_A, output_type="ndarray")
    e = np.unique(np.sort(owner[pairs], axis=1), axis=0) if len(pairs) else np.zeros((0, 2), int)
    edges = [[int(a), int(b)] for a, b in e if a != b]
    return {"resnums": resnums, "probs": probs.tolist(), "edges": edges,
            "seconds": round(secs, 2)}


def _worker(args):
    pesto_root, key, path = args
    try:
        return key, run_one(pesto_root, path), None
    except Exception as exc:                                 # noqa: BLE001
        return key, None, f"{type(exc).__name__}: {exc}"


def score_structures(pesto_root, jobs, cache_dir, workers, refresh=False):
    """{key: per-residue record} for every receptor, cached outside the repo."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out, todo = {}, []
    for key, path in jobs:
        c = cache_dir / f"{key}.json"
        if c.exists() and not refresh:
            out[key] = json.loads(c.read_text())
        else:
            todo.append((str(pesto_root), key, str(path)))
    if todo:
        print(f"[pesto] {len(todo)} structures to run, {len(out)} cached, "
              f"{workers} workers", flush=True)
        import multiprocessing as mp                         # noqa: PLC0415
        done, failed = 0, {}
        ctx = mp.get_context("spawn")
        with ctx.Pool(workers, initializer=_init_worker) as pool:
            for key, rec, err in pool.imap_unordered(_worker, todo):
                done += 1
                if err:
                    failed[key] = err
                    print(f"[pesto] FAIL {key}: {err}", flush=True)
                else:
                    (cache_dir / f"{key}.json").write_text(json.dumps(rec))
                    out[key] = rec
                if done % 20 == 0:
                    print(f"[pesto] {done}/{len(todo)}", flush=True)
        if failed:
            (cache_dir / "_failures.json").write_text(json.dumps(failed, indent=1))
    return out


def _init_worker():
    import torch as pt                                       # noqa: PLC0415
    pt.set_num_threads(3)


# --------------------------------------------------------------------------
# 2. the pre-registered decision rule
# --------------------------------------------------------------------------
def components(sel, edges):
    """Connected components of `sel` under `edges`, as a list of sets."""
    if not sel:
        return []
    parent = {i: i for i in sel}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    comps = {}
    for i in sel:
        comps.setdefault(find(i), set()).add(i)
    return list(comps.values())


def largest_component(sel, edges, p_of):
    """Largest connected component of `sel`; ties -> larger summed p."""
    cs = components(sel, edges)
    if not cs:
        return set()
    return max(cs, key=lambda c: (len(c), sum(p_of[i] for i in c)))


def predict(rec, head, threshold=THRESHOLD, component=True, top_k=None):
    """Per-residue probabilities -> a predicted site, as residue NUMBERS."""
    p = np.asarray(rec["probs"])[:, head]
    if top_k is not None:
        k = min(int(top_k), p.size)
        idx = set(int(i) for i in np.argsort(-p)[:k]) if k else set()
    else:
        idx = set(int(i) for i in np.where(p >= threshold)[0])
        if component:
            idx = largest_component(idx, rec["edges"], {i: float(p[i]) for i in idx})
    return {rec["resnums"][i] for i in idx}


def best_component(rec, head, truth, threshold=THRESHOLD):
    """ORACLE: of the components of S, the one that scores best against truth.

    This arm cheats - it is told the answer. It exists because
    results/epitope_gate.json already publishes the same cheat for P2Rank
    (p2rank_bestany, the best of all offered pockets), and the two oracles have
    to be compared to each other, not one oracle against one blind arm. It
    separates "the model cannot see the site" from "the model sees the site but
    the reducer picks the wrong blob".
    """
    p = np.asarray(rec["probs"])[:, head]
    sel = set(int(i) for i in np.where(p >= threshold)[0])
    cs = components(sel, rec["edges"])
    if not cs:
        return set(), None
    scored = [({rec["resnums"][i] for i in c}, c) for c in cs]
    best = max(scored, key=lambda t: G.score(t[0], truth)["jaccard"])
    ranked = sorted(cs, key=lambda c: (-len(c), -sum(float(p[i]) for i in c)))
    return best[0], ranked.index(best[1]) + 1


def p2rank_residue_scores(path):
    """{residue number: P2Rank ligandability score} from a *_residues.csv.

    The raw `score` column is used, not `probability`: AUC only needs a ranking,
    and `probability` is P2Rank's calibration of `score` (checked monotone on 19
    of 20 sampled files), so the ranking is what carries the information.
    """
    import csv                                              # noqa: PLC0415
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            try:
                out[int(row["residue_label"])] = float(row["score"])
            except (KeyError, ValueError):
                continue
    return out


def ligand_vocabulary_check(pesto_root, sm_rows):
    """Do the ligands we score against even belong to PeSTo's ligand class?

    Verified rather than assumed, because it decides how the small-molecule
    result should be read. PeSTo's ligand head is not trained on "a small
    molecule": categ_to_resnames['ligand'] is a fixed list of 31 common PDB HET
    codes - buffers, sugars, cofactors, nucleotides - and a residue outside that
    list was not a positive for any of the five heads during training.
    """
    _import_pesto(pesto_root)
    from src.data_encoding import categ_to_resnames            # noqa: PLC0415
    vocab = set(categ_to_resnames["ligand"])
    ours = sorted({r["ligand"] for r in sm_rows})
    inside = sorted(set(ours) & vocab)
    by_set = {}
    for label, *_ in SMALL_MOLECULE:
        s = sorted({r["ligand"] for r in sm_rows if r["set"] == label})
        by_set[label] = {"n_distinct_ligands": len(s),
                         "n_in_pesto_ligand_vocabulary": len(set(s) & vocab),
                         "ligands": s}
    return {
        "pesto_ligand_vocabulary": sorted(vocab),
        "n_pesto_ligand_vocabulary": len(vocab),
        "our_distinct_ligands": ours,
        "n_our_distinct_ligands": len(ours),
        "n_ours_in_pesto_vocabulary": len(inside),
        "ours_in_pesto_vocabulary": inside,
        "by_set": by_set,
        "reading": (
            "PeSTo's ligand head predicts where a member of THAT 31-code list "
            "binds. Every ligand in this benchmark is a drug-like inhibitor and "
            "none of them is in the list, so the ligand head is being asked a "
            "question adjacent to, but not the same as, the one P2Rank is asked. "
            "The ATP-site cofactors ATP/ADP/ANP/AMP/ANP ARE in the list and share "
            "the kinase pocket with these inhibitors, so the head is not blind to "
            "the site by construction - but a loss here is a statement about "
            "transfer to drug-like chemistry, not about PeSTo failing its own "
            "benchmark."),
        "no_such_mismatch_on_the_ppi_head": (
            "the protein head's training positives are the 20 standard amino "
            "acids, which is exactly what every binder in the biologic benchmark "
            "is made of, so the biologic arm carries no equivalent caveat."),
    }


def residue_auc(scores_by_res, truth):
    """Per-residue ROC AUC of a score over one receptor. None if degenerate."""
    from sklearn.metrics import roc_auc_score                # noqa: PLC0415
    res = sorted(scores_by_res)
    y = np.array([1 if r in truth else 0 for r in res])
    if y.sum() == 0 or y.sum() == y.size:
        return None
    return float(roc_auc_score(y, np.array([scores_by_res[r] for r in res])))


def random_floor(rng, pool, k, truth):
    """|pred|-sized uniform draws from the receptor - epitope_gate.py's floor."""
    k = min(int(k), len(pool))
    if not k:
        return {**G.score(set(), truth), "n_pred": 0}
    rows = [G.score(set(int(x) for x in rng.choice(pool, size=k, replace=False)), truth)
            for _ in range(N_FLOOR_DRAWS)]
    return {**{m: float(np.mean([r[m] for r in rows])) for m in METRICS}, "n_pred": k}


# --------------------------------------------------------------------------
# 3. the two benchmarks, both read from ground truth this repo already owns
# --------------------------------------------------------------------------
def small_molecule_rows():
    """The 38 KDR + 31 CDK2 apo receptors, their 4.5 A truth and P2Rank's score.

    Truth is reconstructed exactly the way m2_gate.py reconstructs it, and the
    P2Rank vector is READ BACK from results/m2_gate_*.json rather than recomputed,
    so it is the published number and the comparison is genuinely paired.
    """
    out = []
    for label, pipe, res_file, published in SMALL_MOLECULE:
        rec = json.loads((ROOT / "results" / "pipeline" / pipe / "target" /
                          "pocket_recovery.json").read_text())
        m2 = json.loads((ROOT / "results" / res_file).read_text())
        p2 = {r["pdb_id"]: r for r in m2["per_structure"]}
        for s in rec["per_structure"]:
            pid, ch, t1 = s["pdb_id"], s["chain_id"], s["top1"]
            if pid not in p2:
                continue
            truth = ({int(x) for x in t1["overlap_residue_labels"]} |
                     {int(x) for x in t1["missed_known_contact_labels"]})
            assert len(truth) == t1["n_known_contact_residues"], pid
            assert len(truth) == p2[pid]["n_known_contacts"], pid
            pdb = ROOT / "data" / "raw" / "prepared" / f"{pid}_{ch}_stripped.pdb"
            if not pdb.exists():                      # same file, P2Rank's own copy
                pdb = (ROOT / "data" / "raw" / "p2rank" / f"{pid}_{ch}" /
                       "visualizations" / "data" / f"{pid}_{ch}_stripped.pdb")
            bp = s["best_pocket"]
            pr, rc = float(bp["precision"]), float(bp["recall_of_known_contacts"])
            out.append({"set": label, "published_p2rank_jaccard": published,
                        "key": f"sm_{label}_{pid}_{ch}", "pdb_id": pid, "chain": ch,
                        "ligand": s["ligand"], "resolution_a": s["resolution_a"],
                        "pdb": pdb, "truth": truth,
                        "p2rank_residues_csv": (ROOT / "data" / "raw" / "p2rank" /
                                                f"{pid}_{ch}" /
                                                f"{pid}_{ch}_stripped.pdb_residues.csv"),
                        "p2rank_top1": p2[pid]["p2rank_geometry"],
                        "p2rank_bestany": {
                            "precision": pr, "recall": rc,
                            "f1": 2 * pr * rc / (pr + rc) if (pr + rc) else 0.0,
                            "jaccard": float(bp["jaccard"]),
                            "n_pred": int(bp["n_pocket_residues"])}})
    return out


def biologic_rows():
    """The 200 cached biologic complexes, their 4.5 A truth and P2Rank's score."""
    prep = {r["pdb_id"]: r for r in json.loads(
        (ROOT / "data" / "raw" / "epitope_gate" / "prepared.json").read_text())["prepared"]}
    eg = json.loads((ROOT / "results" / "epitope_gate.json").read_text())
    out = []
    for r in eg["per_structure"]:
        pid = r["pdb_id"]
        p = prep[pid]
        truth = set(p["interface_residues"])
        assert len(truth) == r["n_interface_residues"], pid
        out.append({"set": r["klass"], "key": f"bio_{pid}", "pdb_id": pid,
                    "resolution_a": r["resolution_a"],
                    "receptor_uniprot": r["receptor_uniprot"],
                    "binder_description": r["binder_description"],
                    "pdb": Path(p["receptor_pdb"]), "truth": truth,
                    "p2rank_residues_csv": (ROOT / "data" / "raw" / "epitope_gate" /
                                            "p2rank" / "out" / f"{pid}.pdb_residues.csv"),
                    "p2rank_top1": r["p2rank_top1"],
                    "p2rank_bestany": r["p2rank_bestany"]})
    return out


# --------------------------------------------------------------------------
# 4. arms, summaries, paired statistics
# --------------------------------------------------------------------------
def build_arms(rows, scores, head, sweep=True):
    """Score every arm for every row. Returns per-structure records."""
    rng = np.random.default_rng(SEED)
    per = []
    for r in rows:
        rec = scores[r["key"]]
        truth, pool = r["truth"], np.array(rec["resnums"])
        pred = predict(rec, head)
        bcc, bcc_rank = best_component(rec, head, truth)
        p2res = p2rank_residue_scores(r["p2rank_residues_csv"])
        pesto_res = {rn: float(p) for rn, p in
                     zip(rec["resnums"], np.asarray(rec["probs"])[:, head])}
        shared = set(p2res) & set(pesto_res)
        row = {k: r[k] for k in r if k not in ("pdb", "truth", "p2rank_residues_csv")}
        row.update({
            "n_receptor_residues": len(rec["resnums"]),
            "n_truth_residues": len(truth),
            "pesto_seconds": rec["seconds"],
            "pesto": {**G.score(pred, truth), "n_pred": len(pred)},
            "pesto_bestcc_ORACLE": {**G.score(bcc, truth), "n_pred": len(bcc)},
            "pesto_bestcc_rank_by_size": bcc_rank,
            "random_pesto_sized": random_floor(rng, pool, len(pred), truth),
            "n_above_threshold": int((np.asarray(rec["probs"])[:, head] >= THRESHOLD).sum()),
            "n_components_above_threshold":
                len(components(set(int(i) for i in
                                   np.where(np.asarray(rec["probs"])[:, head] >= THRESHOLD)[0]),
                               rec["edges"])),
            "max_p": float(np.asarray(rec["probs"])[:, head].max()),
            # reducer-free: rank every residue, ask how well the ranking separates
            # the true interface. No threshold, no component rule, either side.
            "auc_pesto": residue_auc({k: pesto_res[k] for k in shared}, truth),
            "auc_p2rank": residue_auc({k: p2res[k] for k in shared}, truth),
            "n_residues_scored_by_both": len(shared),
        })
        if sweep:
            sw = {}
            for t in SWEEP_THRESHOLDS:
                for comp in (True, False):
                    s = predict(rec, head, threshold=t, component=comp)
                    sw[f"p>={t:.2f}|{'largest_cc' if comp else 'all'}"] = {
                        **G.score(s, truth), "n_pred": len(s)}
            k = r["p2rank_top1"]["n_pred"]
            s = predict(rec, head, top_k=k)
            sw["top_k=|p2rank_top1|"] = {**G.score(s, truth), "n_pred": len(s)}
            row["sweep"] = sw
        per.append(row)
    return per


def summarise(per, arms):
    out = {}
    for a in arms:
        v = {m: float(np.mean([r[a][m] for r in per])) for m in METRICS}
        lo_hi = M.bootstrap_ci(np.array([r[a]["jaccard"] for r in per]),
                               n_boot=10000, seed=SEED)
        v["jaccard_ci"] = [lo_hi[1], lo_hi[2]]
        v["mean_n_pred"] = float(np.mean([r[a]["n_pred"] for r in per]))
        v["n"] = len(per)
        out[a] = v
    return out


def paired(per, a, b):
    """Paired bootstrap + Wilcoxon on every metric, plus sign counts."""
    o = {}
    for m in METRICS:
        va = np.array([r[a][m] for r in per])
        vb = np.array([r[b][m] for r in per])
        st = M.paired_bootstrap(va, vb, n_boot=20000, seed=SEED)
        st["wilcoxon"] = M.wilcoxon(va, vb)
        d = va - vb
        nz = d[d != 0]
        st["sign_counts"] = {"a_greater": int((d > 0).sum()),
                             "b_greater": int((d < 0).sum()),
                             "tied": int((d == 0).sum()),
                             "median_nonzero_delta": float(np.median(nz)) if nz.size else 0.0}
        st["n"] = len(per)
        o[m] = st
    return o


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pesto-root", default=None)
    ap.add_argument("--fetch", default=None,
                    help="git clone PeSTo into this directory (must be outside the repo)")
    ap.add_argument("--cache", default=None,
                    help="where per-residue scores are cached (default: beside PeSTo)")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "results" / "pesto_vs_p2rank.json"))
    args = ap.parse_args()

    pesto = resolve_pesto(args.pesto_root, args.fetch)
    heads = head_indices(pesto)
    cache = Path(args.cache) if args.cache else pesto.parent / "pesto_scores_cache"
    print(f"[pesto] {pesto}\n[pesto] heads {heads['order']} -> "
          f"ppi={heads['ppi']} ligand={heads['ligand']}\n[pesto] cache {cache}")

    sm_rows, bio_rows = small_molecule_rows(), biologic_rows()
    vocab = ligand_vocabulary_check(pesto, sm_rows)
    print(f"[vocab] PeSTo's ligand class is {vocab['n_pesto_ligand_vocabulary']} fixed "
          f"HET codes; {vocab['n_ours_in_pesto_vocabulary']} of our "
          f"{vocab['n_our_distinct_ligands']} distinct ligands are in it "
          f"{vocab['ours_in_pesto_vocabulary']}")
    missing = [str(r["pdb"]) for r in sm_rows + bio_rows if not Path(r["pdb"]).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} receptor files missing, first: {missing[:3]}")
    print(f"[data] {len(sm_rows)} small-molecule receptors, {len(bio_rows)} biologic")

    t0 = time.time()
    scores = score_structures(pesto, [(r["key"], r["pdb"]) for r in sm_rows + bio_rows],
                              cache, args.workers, args.refresh)
    wall = time.time() - t0
    lost = [r["key"] for r in sm_rows + bio_rows if r["key"] not in scores]
    sm_rows = [r for r in sm_rows if r["key"] in scores]
    bio_rows = [r for r in bio_rows if r["key"] in scores]

    # ---- arms ------------------------------------------------------------
    sm_per = build_arms(sm_rows, scores, heads["ligand"])
    bio_per = build_arms(bio_rows, scores, heads["ppi"])
    ARMS = ("pesto", "pesto_bestcc_ORACLE", "p2rank_top1", "p2rank_bestany",
            "random_pesto_sized")

    sets = {}
    for label, _, _, published in SMALL_MOLECULE:
        rs = [r for r in sm_per if r["set"] == label]
        if rs:
            sets[f"small_molecule_{label}"] = {"head": "ligand", "rows": rs,
                                               "published_p2rank_jaccard": published}
    for k in ("peptide", "ppi", "antibody"):
        rs = [r for r in bio_per if r["set"] == k]
        if rs:
            sets[f"biologic_{k}"] = {"head": "ppi", "rows": rs}
    sets["biologic_ALL"] = {"head": "ppi", "rows": bio_per}

    blocks, family = {}, {}
    for name, s in sets.items():
        rs = s["rows"]
        b = {"head": s["head"], "n": len(rs), "arms": summarise(rs, ARMS),
             "paired": {}}
        if "published_p2rank_jaccard" in s:
            b["published_p2rank_jaccard"] = s["published_p2rank_jaccard"]
        for a, bb in (("pesto", "p2rank_top1"), ("pesto", "random_pesto_sized"),
                      ("p2rank_top1", "random_pesto_sized"),
                      ("pesto_bestcc_ORACLE", "p2rank_bestany"),
                      ("pesto_bestcc_ORACLE", "random_pesto_sized")):
            blk = paired(rs, a, bb)
            b["paired"][f"{a} - {bb}"] = blk
            for m in METRICS:
                family[f"[{name}] {a} - {bb} [{m}]"] = blk[m]["p_value"]
        # reducer-free per-residue AUC, paired on the same structures
        va = np.array([r["auc_pesto"] for r in rs if r["auc_pesto"] is not None
                       and r["auc_p2rank"] is not None])
        vb = np.array([r["auc_p2rank"] for r in rs if r["auc_pesto"] is not None
                       and r["auc_p2rank"] is not None])
        st = M.paired_bootstrap(va, vb, n_boot=20000, seed=SEED)
        st["wilcoxon"] = M.wilcoxon(va, vb)
        st["n"] = int(va.size)
        st["n_dropped_degenerate"] = len(rs) - int(va.size)
        for arm, v in (("pesto", va), ("p2rank", vb)):
            lo_hi = M.bootstrap_ci(v, n_boot=10000, seed=SEED)
            st[f"mean_auc_{arm}"] = float(v.mean())
            st[f"ci_auc_{arm}"] = [lo_hi[1], lo_hi[2]]
            st[f"frac_above_0.5_{arm}"] = float((v > 0.5).mean())
        b["residue_auc"] = st
        family[f"[{name}] residue_auc pesto - p2rank"] = st["p_value"]
        blocks[name] = b
    for k, adj in M.holm_bonferroni(family).items():
        name = k[1:k.index("]")]
        rest = k[k.index("] ") + 2:]
        if rest == "residue_auc pesto - p2rank":
            blocks[name]["residue_auc"]["p_holm"] = adj
            continue
        pair, met = rest[:rest.rindex(" [")], rest[rest.rindex("[") + 1:-1]
        blocks[name]["paired"][pair][met]["p_holm"] = adj

    # ---- sweep (secondary) ----------------------------------------------
    sweep = {}
    for name, s in sets.items():
        cells = sorted(s["rows"][0]["sweep"])
        sweep[name] = {c: {"jaccard": float(np.mean([r["sweep"][c]["jaccard"] for r in s["rows"]])),
                           "precision": float(np.mean([r["sweep"][c]["precision"] for r in s["rows"]])),
                           "recall": float(np.mean([r["sweep"][c]["recall"] for r in s["rows"]])),
                           "mean_n_pred": float(np.mean([r["sweep"][c]["n_pred"] for r in s["rows"]])),
                           "n": len(s["rows"])}
                       for c in cells}

    # ---- console ---------------------------------------------------------
    print(f"\n=== E1: PeSTo vs P2Rank, identical 4.5 A truth, paired per structure ===")
    print(f"{'set':<22}{'head':>8}{'n':>5}{'PeSTo J':>10}{'P2Rank J':>10}"
          f"{'rand J':>9}{'delta':>9}{'Holm p':>11}{'Wilcox p':>11}{'|PeSTo|':>9}")
    print("-" * 104)
    order = [f"small_molecule_{l}" for l, *_ in SMALL_MOLECULE] + \
            ["biologic_peptide", "biologic_ppi", "biologic_antibody", "biologic_ALL"]
    for name in order:
        if name not in blocks:
            continue
        b = blocks[name]
        st = b["paired"]["pesto - p2rank_top1"]["jaccard"]
        print(f"{name:<22}{b['head']:>8}{b['n']:>5}"
              f"{b['arms']['pesto']['jaccard']:>10.4f}"
              f"{b['arms']['p2rank_top1']['jaccard']:>10.4f}"
              f"{b['arms']['random_pesto_sized']['jaccard']:>9.4f}"
              f"{st['delta']:>+9.4f}{st['p_holm']:>11.3g}"
              f"{st['wilcoxon']['p_value']:>11.3g}"
              f"{b['arms']['pesto']['mean_n_pred']:>9.1f}")
    print(f"  Holm p is adjusted over a family of {len(family)} paired tests. The "
          f"bootstrap p has a\n  resolution floor of 2/(20000+1) = 1.0e-4, so "
          f"{len(family)} x that floor = {len(family) * 2 / 20001:.4f} is the "
          f"smallest\n  Holm value reachable - a value sitting exactly there means "
          f"'below the bootstrap's\n  resolution', not 'p = 0.0126'. The Wilcoxon "
          f"column has no such floor and is the\n  cross-check; both are stored per "
          f"comparison in the JSON.")

    print(f"\n--- is each arm above its own size-matched random floor? (Jaccard) ---")
    for name in order:
        if name not in blocks:
            continue
        b = blocks[name]
        for arm in ("pesto", "pesto_bestcc_ORACLE", "p2rank_top1"):
            st = b["paired"][f"{arm} - random_pesto_sized"]["jaccard"]
            tag = "ABOVE CHANCE" if st["delta"] > 0 and st["p_holm"] < 0.05 else "at chance"
            print(f"  {name:<22}{arm:<20} delta {st['delta']:+.4f} "
                  f"[{st['ci_lo']:+.4f}, {st['ci_hi']:+.4f}]  Holm p={st['p_holm']:.3g}  {tag}")

    print(f"\n--- ORACLE vs ORACLE: best PeSTo component vs best P2Rank pocket (Jaccard) ---")
    print("    both arms are told the answer; this asks whether the SITE is visible "
          "at all,\n    separately from whether the blind reducer picks it.")
    for name in order:
        if name not in blocks:
            continue
        b = blocks[name]
        st = b["paired"]["pesto_bestcc_ORACLE - p2rank_bestany"]["jaccard"]
        print(f"  {name:<22} PeSTo-bestCC {b['arms']['pesto_bestcc_ORACLE']['jaccard']:.4f}"
              f"   P2Rank-bestany {b['arms']['p2rank_bestany']['jaccard']:.4f}"
              f"   delta {st['delta']:+.4f}  Holm p={st['p_holm']:.3g}")

    print(f"\n--- REDUCER-FREE: per-residue ROC AUC of the raw score (n = structures) ---")
    print("    no threshold and no component rule on either side - does the ranking "
          "separate\n    the true interface from the rest of the surface?")
    for name in order:
        if name not in blocks:
            continue
        a = blocks[name]["residue_auc"]
        print(f"  {name:<22} n={a['n']:<4} PeSTo {a['mean_auc_pesto']:.4f} "
              f"[{a['ci_auc_pesto'][0]:.4f}, {a['ci_auc_pesto'][1]:.4f}]   "
              f"P2Rank {a['mean_auc_p2rank']:.4f} "
              f"[{a['ci_auc_p2rank'][0]:.4f}, {a['ci_auc_p2rank'][1]:.4f}]   "
              f"delta {a['delta']:+.4f}  Holm p={a['p_holm']:.3g}")

    print(f"\n--- labelled sweep, Jaccard (primary cell = p>=0.50|largest_cc) ---")
    cells = sorted(next(iter(sweep.values())))
    print(f"{'set':<22}" + "".join(f"{c:>22}" for c in cells))
    for name in order:
        if name in sweep:
            print(f"{name:<22}" + "".join(f"{sweep[name][c]['jaccard']:>22.4f}" for c in cells))

    # ---- verdict ---------------------------------------------------------
    verdict = {}
    for label, _, _, published in SMALL_MOLECULE:
        name = f"small_molecule_{label}"
        if name not in blocks:
            continue
        b = blocks[name]
        st = b["paired"]["pesto - p2rank_top1"]["jaccard"]
        verdict[name] = {
            "criterion": f"mean Jaccard >= {published} OR paired CI on (PeSTo - P2Rank) contains 0",
            "pesto_jaccard": b["arms"]["pesto"]["jaccard"],
            "reaches_published": b["arms"]["pesto"]["jaccard"] >= published,
            "ci_contains_zero": bool(st["ci_lo"] <= 0 <= st["ci_hi"]),
            "pass": bool(b["arms"]["pesto"]["jaccard"] >= published or
                         (st["ci_lo"] <= 0 <= st["ci_hi"])),
        }
    for k in ("peptide", "ppi", "antibody"):
        name = f"biologic_{k}"
        if name not in blocks:
            continue
        st = blocks[name]["paired"]["pesto - random_pesto_sized"]["jaccard"]
        p2 = blocks[name]["paired"]["p2rank_top1 - random_pesto_sized"]["jaccard"]
        au = blocks[name]["residue_auc"]
        verdict[name] = {
            "criterion": "PeSTo PPI head beats its size-matched random floor, Holm p < 0.05 on Jaccard",
            "pesto_jaccard": blocks[name]["arms"]["pesto"]["jaccard"],
            "delta_vs_floor": st["delta"], "p_holm": st["p_holm"],
            "pass": bool(st["delta"] > 0 and st["p_holm"] < 0.05),
            "p2rank_above_same_floor": bool(p2["delta"] > 0 and p2["p_holm"] < 0.05),
            "reducer_free_auc_pesto": au["mean_auc_pesto"],
            "reducer_free_auc_p2rank": au["mean_auc_p2rank"],
        }
    print("\n=== PASS / FAIL against the pre-registered criteria ===")
    for k, v in verdict.items():
        print(f"  {k:<22}{'PASS' if v['pass'] else 'FAIL'}   {v['criterion']}")

    out = {
        "experiment": "e1-pesto - can one model replace the site-definition stage "
                      "for every modality?",
        "question": "PeSTo emits a protein-interface head and a ligand-interface head "
                    "from one forward pass. Does either beat P2Rank 2.5 on the ground "
                    "truth this repo already owns - and does the PPI head beat chance "
                    "on antibody epitopes, where P2Rank does not?",
        "ground_truth": "receptor residues with any heavy atom within 4.5 A of any "
                        "binder heavy atom. Read back from "
                        "results/pipeline/*/target/pocket_recovery.json (small molecule) "
                        "and data/raw/epitope_gate/prepared.json (biologic); both were "
                        "written by scripts/interfaces.py. No contact routine is "
                        "implemented in this script.",
        "preregistration": {
            "frozen": "before the first forward pass on either benchmark",
            "model": {"repo": PESTO_REPO, "release": RELEASE, "checkpoint": CHECKPOINT,
                      "device": "cpu", "licence": "CC BY-NC-SA 4.0 - located at runtime, "
                                                  "never vendored into this MIT repo"},
            "heads": heads,
            "decision_rule": {
                "1_threshold": f"S = residues with sigmoid(logit) >= {THRESHOLD}",
                "2_adjacency": f"u,v adjacent iff min heavy-atom distance <= {CONTACT_A} A",
                "3_site": "largest connected component of S; ties -> larger summed p",
                "4_empty": "S empty -> empty prediction, scores 0 (same treatment "
                           "epitope_gate.py gives a structure with no P2Rank pocket)",
            },
            "scoring": "m2_gate.score() imported verbatim",
            "floor": f"{N_FLOOR_DRAWS} uniform draws of |pred| residues from the "
                     f"receptor, averaged; size-matched to PeSTo",
            "statistics": "metrics.paired_bootstrap(20000, seed 0), metrics.wilcoxon, "
                          "metrics.bootstrap_ci(10000, seed 0), one "
                          "metrics.holm_bonferroni family over every paired test below",
            "pass_criteria": {
                "small_molecule": "mean Jaccard >= 0.6355 (KDR) / 0.4938 (CDK2), or the "
                                  "paired CI on (PeSTo - P2Rank) contains 0",
                "biologic": "PPI head beats its size-matched random floor, Holm p < 0.05",
            },
            "sweep_declared_in_advance": {
                "thresholds": list(SWEEP_THRESHOLDS),
                "component_filter": [True, False],
                "extra": "top-k with k = |P2Rank top-1 pocket| for the same structure",
                "primary_cell": f"p>={THRESHOLD:.2f}|largest_cc",
            },
            "added_after_the_freeze_and_why": [
                {"item": "pesto_bestcc_ORACLE arm",
                 "added": "after a 4-structure smoke test, before the full run",
                 "why": "on 7RKE the PPI head put p >= 0.97 on the true epitope AND "
                        "p >= 0.5 on a second, larger protein-interface patch, so the "
                        "pre-registered largest-component rule returned the wrong blob "
                        "and scored 0. Without an oracle arm the table cannot tell "
                        "'the model cannot see the site' from 'the reducer picks the "
                        "wrong site'. This is not a new primary: it is the exact "
                        "counterpart of p2rank_bestany, an oracle results/epitope_gate.json "
                        "already publishes, and it is compared only against that oracle.",
                 "not_a_tuning_knob": "it is labelled ORACLE everywhere and can never "
                                      "be run blind; the pre-registered primary is "
                                      "unchanged and is reported first."},
                {"item": "reducer-free per-residue ROC AUC (both methods)",
                 "added": "after the same smoke test, before the full run",
                 "why": "AUC is the metric PeSTo's own paper reports and it uses no "
                        "threshold and no component rule, so it measures the heads "
                        "rather than this script's reducer. P2Rank's *_residues.csv "
                        "gives the same kind of per-residue score, so the comparison "
                        "is symmetric.",
                 "not_a_tuning_knob": "no decision rule was changed to improve it; it "
                                      "is a diagnostic reported alongside, not a "
                                      "pass criterion."},
            ],
        },
        "verdict": verdict,
        "p_value_floor": {
            "family_size": len(family),
            "bootstrap_resolution_floor": 2 / 20001,
            "smallest_reachable_holm_p": len(family) * 2 / 20001,
            "note": "metrics.paired_bootstrap's p is an achieved significance level "
                    "with add-one smoothing, so it cannot go below 2/(n_boot+1). Any "
                    "Holm value sitting exactly at smallest_reachable_holm_p means the "
                    "bootstrap could not resolve a smaller p, not that the p IS that "
                    "number. Every comparison also carries an unfloored Wilcoxon "
                    "signed-rank p under ['wilcoxon']['p_value']; read that one for "
                    "how strong the evidence actually is.",
            "family_note": "biologic_ALL re-tests the same 200 structures already "
                           "tested as peptide/ppi/antibody, so the family is larger "
                           "than the number of independent questions and the "
                           "correction is conservative.",
        },
        "ligand_vocabulary_check": vocab,
        "sets": {k: {kk: vv for kk, vv in v.items()} for k, v in blocks.items()},
        "sweep_jaccard": sweep,
        "runtime": {
            "n_structures": len(scores),
            "n_run_this_invocation": len(scores) - len(
                [k for k in scores if (cache / f"{k}.json").stat().st_mtime < t0]),
            "wall_seconds_this_invocation": round(wall, 1),
            "mean_forward_seconds": float(np.mean([s["seconds"] for s in scores.values()])),
            "sum_forward_seconds": float(np.sum([s["seconds"] for s in scores.values()])),
            "workers": args.workers, "device": "cpu",
            "note": "wall_seconds_this_invocation is ~0 on a cached re-run; the cost "
                    "that matters is sum_forward_seconds, which is CPU time across all "
                    f"{len(scores)} structures, divided by however many workers are "
                    "used. mean_forward_seconds is higher than the ~5 s PeSTo reports "
                    "because 5 workers share this machine's cores; it is not a "
                    "like-for-like single-structure latency.",
        },
        "per_structure": {"small_molecule": sm_per, "biologic": bio_per},
        "not_evaluated": [
            {"item": "PeSTo ion / nucleic-acid / lipid heads",
             "why": "this repo has no 4.5 A ground truth for ion, nucleic-acid or lipid "
                    "interfaces, so they cannot be scored here. They are emitted by the "
                    "same forward pass and are recorded in the cached per-residue "
                    "probabilities if a benchmark is ever built."},
            {"item": "PeSTo releases i_v3_0, i_v3_1, i_v4_0",
             "why": "the pre-registration fixed i_v4_1 (the release apply_model.ipynb "
                    "selects) before any run. Sweeping releases after seeing a number "
                    "would be tuning on the answer."},
            {"item": "a PeSTo-conditioned co-folding arm",
             "why": "E1 measures site definition only. Whether a better site actually "
                    "improves the downstream ranking is a separate experiment."},
            {"item": "antibody-specific epitope predictors (Discotope, SEMA, epitope3D)",
             "why": "out of scope for E1, which asks whether ONE model covers every "
                    "modality. They are single-modality by construction."},
            {"item": "whether PeSTo saw these structures in training",
             "why": "PeSTo's training split is a subunit list over the whole PDB and "
                    "these 269 receptors were not selected to avoid it. Measuring the "
                    "overlap needs the authors' subunits_train_set.txt matched to "
                    "these entries, which was not done here. Stated as a caveat "
                    "instead of asserted either way."},
            {"item": "a blind reducer that picks the right component",
             "why": "the gap between the pre-registered blind arm and the oracle arm "
                    "below says how much such a reducer could be worth, but inventing "
                    "one and scoring it on the same 200 complexes that motivated it "
                    "would be fitting to this benchmark. It belongs in a follow-up "
                    "with a held-out set."},
        ],
        "credits_spent": {"rowan": 0,
                          "note": "local CPU only - this script makes no Rowan call "
                                  "and imports no Rowan client; its only network use "
                                  "is the optional --fetch git clone of PeSTo."},
        "licence_note": "PeSTo is CC BY-NC-SA 4.0. This script locates a checkout at "
                        "runtime (--pesto-root / $PESTO_ROOT / --fetch) and refuses one "
                        "inside this repository. No PeSTo weights or source are stored "
                        "here; only numbers computed from them.",
        "caveats": [
            "PeSTo was trained on protein interfaces from the PDB; some of these 200 "
            "biologic complexes and 69 small-molecule receptors may sit in its training "
            "set. This is a ceiling on PeSTo, not on P2Rank, and any PeSTo win here is "
            "an optimistic estimate. Nothing in this script removes that overlap.",
            "the receptors are the ligand- and binder-STRIPPED forms, so both methods "
            "are blind in the same way; but they are holo conformations with the binder "
            "deleted, not true apo structures.",
            "the biologic P2Rank vectors are read back from results/epitope_gate.json "
            "and the small-molecule ones from results/m2_gate_*.json, so they are the "
            "published numbers; nothing is re-run for P2Rank.",
        ],
    }
    if lost:
        out["structures_lost_to_pesto_errors"] = lost
    dest = Path(args.out)
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
