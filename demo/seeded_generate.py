"""Arm S, seeded: build molecules FROM the known binders instead of from all of chemical space.

    seed        the approved drugs whose *mechanism* target is F2        [data/approved_drugs.csv]
      -> BRICS  decompose them into a fragment pool                      [RDKit]
      -> build  recombine the pool (scrambleReagents=True)               [RDKit]
      -> gate   the SAME validity / MW / QED / not-already-approved gate [demo/generate.py]
      -> match  ECFP4-2048 Tanimoto vs the 2,153 approved drugs          [scripts/representations.py]
      -> test   seeded vs unseeded, plus a decoy-seeded circularity control

**Why this exists.** `demo/generate.py` samples `entropy/gpt2_zinc_87m`, a SMILES
language model that is not conditioned on the target. It samples all of drug-like
space, so almost nothing lands near thrombin chemistry: of its 591 kept molecules,
the best nearest neighbour in the approved library was dithranol (an anthralin for
psoriasis) at Tanimoto 0.579, and nothing in its top ten had a thrombin
annotation. Seeding puts the generator in the right neighbourhood by construction.

**THE CIRCULARITY, STATED UP FRONT.** This module draws its fragments from four
approved drugs and then matches the products back against the library those four
drugs are in. That is partly circular *by construction*, and no amount of
statistics undoes it:

  * Every heavy atom of every molecule produced here came from a seed fragment --
    BRICS cannot do anything else. `seed_atom_fraction` puts a number on it
    against the unseeded model's own background.
  * A high nearest-neighbour Tanimoto is therefore the expected result, not a
    discovery. Recombining pieces of approved drugs produces things that look
    like approved drugs. That is arithmetic.

So this module does not report the seeded-vs-unseeded Tanimoto lift on its own.
It reports it next to three controls designed to take it away:

  1. **Decoy-seeded panels** -- the identical BRICS pipeline seeded from randomly
     drawn, MW-matched approved drugs with no F2 annotation. Any lift that also
     shows up here is generic "recombining approved drugs", not target chemistry.
  2. **Leave-the-seeds-out** -- the four seed drugs deleted from the corpus before
     matching. A thrombin hit that survives this one is at least not the molecule
     recognising its own parents.
  3. **LibInvent** (step 2, REINVENT4) -- scaffold fixed, but the R-groups come
     from a ChEMBL prior rather than from the approved library.

What those controls did to the result, on the run in results/demo/seeded/ (see
docs/13-SEEDED-GENERATION.md for the full table):

  * Seeding raises median NN Tanimoto 0.2639 -> 0.3066 (p = 1.25e-45) -- and the
    DECOY-seeded arm reaches 0.3623, beating it. None of the Tanimoto lift is
    target-specific.
  * Seeding raises "nearest neighbour is F2-annotated" from 5/591 to 151/378, but
    135 of those 151 are the seed drugs themselves, and all ten of the top ten are
    argatroban. Delete the seeds and 151 becomes 23, all of them captopril, an
    off-target annotation.
  * The one non-circular hit is LibInvent finding nafamostat -- n = 1, p = 0.54.

Seeding fixes the symptom, not the problem. It is kept because the controls are
the reusable part.

Run:  ./env-kit/bin/python -m demo.seeded_generate
      ./env-kit/bin/python -m demo.test_seeded_generate
No network, no credentials, no credits. RDKit and scipy only.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from demo import generate as unseeded_arm   # noqa: E402  the gate and the corpus loader
from demo import match_direct as md         # noqa: E402

SYMBOL = "F2"
UNSEEDED_BASELINE = REPO / "results" / "demo" / "thrombin" / "07_generated.json"
OUT_DIR = REPO / "results" / "demo" / "seeded"

N_TARGET = 600          # match demo/generate.py's DEMO_N_GENERATE so n is comparable
BUILD_BUDGET_S = 180.0  # BRICSBuild is a generator over a combinatorial space; cap it
MAX_FRAGMENTS = 40      # the build is superexponential in pool size; cap it too
MAX_DEPTH = 3
MIN_FRAGMENT_ATOMS = 5   # size floor for the circularity meter; see fragment_cores
SEED = 0
N_DECOY_PANELS = 3
DECOY_MAX_MW = 900.0    # BRICS on a peptide is a hang, not a control; see mw_matched_decoys

# REINVENT4 pins torch==2.12; env-kit is on 2.14, so it gets its own interpreter.
#   uv venv --python 3.12 --seed env-reinvent
#   env-reinvent/bin/python -m pip install /path/to/REINVENT4   # NOT uv pip: see docs/13
#   uv pip install --python env-reinvent/bin/python scipy
#   curl -L -o priors/libinvent.prior \
#     https://zenodo.org/api/records/20701824/files/libinvent.prior/content
REINVENT_PYTHON = Path(os.environ.get("REINVENT_PYTHON", REPO / "env-reinvent/bin/python"))
LIBINVENT_PRIOR = Path(os.environ.get("LIBINVENT_PRIOR", REPO / "priors/libinvent.prior"))

# Thrombin's S1 pocket ends in Asp189, so its small-molecule inhibitors are built
# around a basic arginine mimetic. These are the three groups that shows up as:
# amidine (dabigatran's benzamidine), guanidine (argatroban's arginine side
# chain), amidoxime (ximelagatran's prodrug-masked benzamidine).
ARG_MIMETIC = {
    "amidine": "[CX3](=[NX2])[NX3]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "benzamidine": "c1ccccc1[CX3](=[NX2])[NX3]",
}


# ---------------------------------------------------------------------------
# seeds and fragments
# ---------------------------------------------------------------------------


def seed_drugs(approved: list[dict], symbol: str = SYMBOL) -> list[dict]:
    """The approved drugs whose *mechanism* target is `symbol`.

    Mechanism, not any-target: the 12 approved drugs with F2 anywhere in
    `targets` include captopril, sitosterol and cianidanol, which are F2
    annotations of a very different kind. The four with F2 as `moa_targets` are
    argatroban, bivalirudin, dabigatran etexilate and ximelagatran.
    """
    return [r for r in approved if symbol in (r["moa_targets"] or "").split(";")]


def fragment_pool(seeds: list[dict]) -> dict[str, list[str]]:
    """BRICS fragments of the seed drugs, mapped back to the drugs they came from."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import BRICS

    RDLogger.DisableLog("rdApp.*")
    pool: dict[str, list[str]] = {}
    for drug in seeds:
        mol = Chem.MolFromSmiles(drug["smiles"])
        if mol is None:
            continue
        for frag in BRICS.BRICSDecompose(mol):
            pool.setdefault(frag, []).append(drug["name"])
    return dict(sorted(pool.items()))


def fragment_cores(frag_smiles: list[str], *, min_atoms: int = MIN_FRAGMENT_ATOMS):
    """Fragments with their BRICS dummy atoms stripped, usable as substructure queries.

    `min_atoms` is not a tuning knob, it is what makes the measurement mean
    anything. Without a size floor the meter is degenerate: the F2 pool contains
    fragments like `OC` and `NC(=N)N`, which occur in most drug-like molecules,
    so median `seed_atom_fraction` comes out at exactly 1.000 for BRICS products,
    for LibInvent products, for the unseeded GPT-2 molecules AND for aspirin.
    Measured, all four. At a floor of 5 heavy atoms the same four read 0.700,
    0.682, 0.286 and 0.538, which is a measurement rather than a tautology.
    """
    from rdkit import Chem

    cores = []
    for smi in frag_smiles:
        frag = Chem.MolFromSmiles(smi)
        if frag is None:
            continue
        core = Chem.DeleteSubstructs(frag, Chem.MolFromSmarts("[#0]"))
        if core.GetNumHeavyAtoms() < max(1, min_atoms):
            continue
        Chem.SanitizeMol(core, catchErrors=True)
        cores.append((smi, core))
    return cores


def _brics_worker(frag_smiles: list[str], n: int, seed: int, max_depth: int, queue) -> None:
    """BRICSBuild in a child process, streaming products back as it finds them.

    In-process this was not safe. `BRICSBuild` with `uniquify=True` can spend
    unbounded wall-clock between yields -- a pool of peptide fragments
    (a decoy panel matched to bivalirudin's size) generates duplicate after
    duplicate and emits nothing -- so a budget check inside the `for` body never
    runs and the whole thing hangs. Measured: two runs pinned a core for 3-5
    minutes with a 20s budget set. The budget has to be enforced by someone who
    is not blocked on the generator, hence a child process the parent can kill.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import BRICS

    RDLogger.DisableLog("rdApp.*")
    random.seed(seed)                       # scrambleReagents uses the `random` module
    mols = [m for m in (Chem.MolFromSmiles(s) for s in sorted(frag_smiles)) if m is not None]
    made = 0
    for product in BRICS.BRICSBuild(mols, scrambleReagents=True, maxDepth=max_depth):
        product.UpdatePropertyCache(strict=False)
        queue.put(Chem.MolToSmiles(product))
        made += 1
        if made >= n:
            break
    queue.put(None)                         # sentinel: the generator finished on its own


def brics_enumerate(frag_smiles: list[str], n: int, *, seed: int = SEED,
                    budget_s: float = BUILD_BUDGET_S, max_depth: int = MAX_DEPTH
                    ) -> tuple[list[str], dict]:
    """Recombine the fragment pool. Returns raw SMILES and how the build terminated.

    `BRICSBuild` enumerates a combinatorial space that does not finish, so the
    honest thing is to record *why* enumeration stopped -- quota, wall-clock
    budget, or genuine exhaustion -- rather than silently truncate.
    """
    import multiprocessing as mp
    import queue as queue_mod

    pool = sorted(frag_smiles)[:MAX_FRAGMENTS]
    started = time.time()
    ctx = mp.get_context("spawn")
    channel = ctx.Queue()
    child = ctx.Process(target=_brics_worker, args=(pool, n, seed, max_depth, channel),
                        daemon=True)
    child.start()

    # Polled with get_nowait() rather than get(timeout=...) on purpose. macOS has
    # no sem_timedwait, so CPython emulates a timed semaphore acquire with its own
    # poll loop, and a run of this module was caught wedged inside exactly that
    # (`sample` showed the parent parked in select_poll_poll for 10+ minutes with a
    # 240s budget set, while its worker kept a core at 99%). A plain non-blocking
    # get plus an explicit sleep keeps the deadline in this function's own hands.
    out: list[str] = []
    stop = "exhausted"
    deadline = started + budget_s
    while True:
        if time.time() >= deadline:
            stop = "time_budget"
            break
        try:
            item = channel.get_nowait()
        except queue_mod.Empty:
            if not child.is_alive():
                stop = "child_died"
                break
            time.sleep(0.02)
            continue
        if item is None:
            stop = "exhausted"
            break
        out.append(item)
        if len(out) >= n:
            stop = "quota"
            break

    channel.close()
    channel.cancel_join_thread()        # never block this process on the dead queue
    child.terminate()
    child.join(timeout=5)
    if child.is_alive():
        child.kill()
        child.join(timeout=5)

    return out, {"stop_reason": stop, "n_raw": len(out),
                 "seconds": round(time.time() - started, 1),
                 "n_fragments": len(pool), "n_fragments_available": len(frag_smiles),
                 "max_fragments": MAX_FRAGMENTS, "max_depth": max_depth}


# ---------------------------------------------------------------------------
# what a molecule inherited from the seeds
# ---------------------------------------------------------------------------


def motif_hits(smiles: str) -> list[str]:
    """Which arginine-mimetic groups this molecule carries. Empty list is a real answer."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [name for name, smarts in ARG_MIMETIC.items()
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))]


def seed_atom_fraction(smiles: str, cores) -> float:
    """Fraction of this molecule's heavy atoms covered by a seed fragment of >=5 atoms.

    This is the circularity meter, and it is only interpretable against its own
    background. The unseeded GPT-2 molecules -- which never saw a seed drug --
    still score a median 0.292, because a five-atom fragment is not a rare thing.
    So the number to read is the gap between a generator and that 0.292, not the
    number itself. See `fragment_cores` for why the size floor exists.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return float("nan")
    covered: set[int] = set()
    for _, core in cores:
        for match in mol.GetSubstructMatches(core, useChirality=False):
            covered.update(match)
    return len(covered) / mol.GetNumHeavyAtoms()


# ---------------------------------------------------------------------------
# generation
# ---------------------------------------------------------------------------


def generate_seeded(seeds: list[dict], approved: list[dict], *, n: int = N_TARGET,
                    seed: int = SEED, budget_s: float = BUILD_BUDGET_S) -> dict:
    """Seed -> fragments -> BRICS build -> the same drug-likeness gate demo/generate.py uses."""
    pool = fragment_pool(seeds)
    cores = fragment_cores(list(pool))
    raw, build = brics_enumerate(list(pool), n, seed=seed, budget_s=budget_s)

    keys = {r["inchikey"] for r in approved if r["inchikey"]}
    kept, reasons = unseeded_arm.drug_like(raw, keys)

    for mol in kept:
        mol["motifs"] = motif_hits(mol["smiles"])
        mol["seed_atom_fraction"] = round(seed_atom_fraction(mol["smiles"], cores), 4)

    return {
        "generator": "RDKit BRICS recombination of the F2-mechanism approved drugs",
        "seeds": [{"name": s["name"], "smiles": s["smiles"]} for s in seeds],
        "fragment_pool": {"n": len(pool), "fragments": pool},
        "build": build,
        "filters": {"mw_range": unseeded_arm.MW_RANGE, "min_qed": unseeded_arm.MIN_QED,
                    "source": "demo/generate.py drug_like(), imported not reimplemented"},
        "dropped": reasons,
        "n_kept": len(kept),
        "molecules": kept,
    }


# ---------------------------------------------------------------------------
# step 2: a LEARNED decoration model (REINVENT4 LibInvent)
#
# BRICS can only re-shuffle atoms the seed drugs already had, which is why its
# circularity meter reads 1.0. LibInvent is a different proposition: the scaffold
# is fixed (so the pharmacophore is guaranteed) but the R-groups come from a
# prior trained on ChEMBL, not from the approved library. If seeding is going to
# produce anything that is not a rearrangement of approved drugs, this is where
# it comes from.
#
# REINVENT4 pins torch==2.12 and this repo's env-kit is on torch 2.14, so it gets
# its own interpreter and is driven as a subprocess. Point --reinvent-python at
# it; if it is absent the arm records itself as not evaluated rather than being
# silently skipped.
# ---------------------------------------------------------------------------

# Each scaffold below is checked against the seed drugs at run time by
# `verify_scaffolds` and dropped if it is not a substructure of one of them, so
# "derived from the known binders" is a computed claim rather than an assertion.
LIBINVENT_SCAFFOLDS = {
    "arginine_guanidine": "[*:0]C(CCCNC(=N)N)[*:1]",        # argatroban / bivalirudin
    "benzamidine_para_amine": "[*:0]NC(=N)c1ccc(N[*:1])cc1",  # dabigatran etexilate
    "benzylamine_amidine": "[*:0]Cc1ccc(C(=N)N[*:1])cc1",     # ximelagatran
}


def verify_scaffolds(scaffolds: dict[str, str], seeds: list[dict]) -> dict[str, list[str]]:
    """Which seed drugs each scaffold is actually a substructure of. Empty list = drop it."""
    import re

    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    parsed = [(s["name"], Chem.MolFromSmiles(s["smiles"])) for s in seeds]
    out: dict[str, list[str]] = {}
    for label, smi in scaffolds.items():
        query = Chem.MolFromSmarts(re.sub(r"\[\*:\d+\]", "[*]", smi))
        out[label] = [name for name, mol in parsed
                      if mol is not None and query is not None and mol.HasSubstructMatch(query)]
    return out


def libinvent_generate(seeds: list[dict], approved: list[dict], *, python: Path, prior: Path,
                       n: int = N_TARGET, device: str = "mps", workdir: Path | None = None,
                       cores=None) -> dict:
    """Decorate the verified thrombin scaffolds with REINVENT4's LibInvent prior."""
    import subprocess
    import tempfile

    provenance = verify_scaffolds(LIBINVENT_SCAFFOLDS, seeds)
    used = {k: v for k, v in LIBINVENT_SCAFFOLDS.items() if provenance[k]}
    dropped = {k: "not a substructure of any seed drug" for k in LIBINVENT_SCAFFOLDS
               if not provenance[k]}
    if not used:
        return {"status": "not_evaluated",
                "reason": "no candidate scaffold is a substructure of a seed drug",
                "scaffold_provenance": provenance}

    binary = python.parent / "reinvent"
    if not binary.exists():
        return {"status": "not_evaluated",
                "reason": f"REINVENT4 console script not found at {binary}"}
    if not prior.exists():
        return {"status": "not_evaluated", "reason": f"LibInvent prior not found at {prior}"}

    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="libinvent-"))
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "scaffolds.smi").write_text("\n".join(used.values()) + "\n")
    per_scaffold = max(1, -(-n // len(used)))
    (tmp / "sample.toml").write_text(
        f'run_type = "sampling"\ndevice = "{device}"\n\n[parameters]\n'
        f'model_file = "{prior}"\nsmiles_file = "scaffolds.smi"\n'
        f'output_file = "sampled.csv"\nnum_smiles = {per_scaffold}\n'
        f"unique_molecules = true\nrandomize_smiles = true\n")

    started = time.time()
    proc = subprocess.run([str(binary), "sample.toml"], cwd=tmp,
                          capture_output=True, text=True, timeout=1800)
    log = (proc.stdout + proc.stderr)[-4000:]
    csv_path = tmp / "sampled.csv"
    if proc.returncode != 0 or not csv_path.exists():
        return {"status": "not_evaluated",
                "reason": f"reinvent exited {proc.returncode} or wrote no CSV",
                "log_tail": log}

    import csv as csv_mod

    raw, scaffold_of = [], {}
    with csv_path.open() as handle:
        for row in csv_mod.DictReader(handle):
            if row.get("SMILES"):
                raw.append(row["SMILES"])
                scaffold_of[row["SMILES"]] = row.get("Scaffold", "")

    keys = {r["inchikey"] for r in approved if r["inchikey"]}
    kept, reasons = unseeded_arm.drug_like(raw, keys)
    for mol in kept:
        mol["motifs"] = motif_hits(mol["smiles"])
        mol["seed_atom_fraction"] = (
            round(seed_atom_fraction(mol["smiles"], cores), 4) if cores else None)

    return {
        "status": "evaluated",
        "generator": "REINVENT4 LibInvent prior, scaffold decoration",
        "prior": str(prior),
        "device": device,
        "scaffolds": used,
        "scaffold_provenance": {k: provenance[k] for k in used},
        "scaffolds_dropped": dropped,
        "num_smiles_per_scaffold": per_scaffold,
        "n_raw": len(raw),
        "seconds": round(time.time() - started, 1),
        "filters": {"mw_range": unseeded_arm.MW_RANGE, "min_qed": unseeded_arm.MIN_QED,
                    "source": "demo/generate.py drug_like(), imported not reimplemented"},
        "dropped": reasons,
        "n_kept": len(kept),
        "molecules": kept,
        "determinism": ("REINVENT4's sampling run mode exposes no RNG seed, so this arm is "
                        "NOT reproducible bit-for-bit; n and the statistics are."),
        "workdir": str(tmp),
        "log_tail": log[-1500:],
    }


def mw_matched_decoys(approved: list[dict], seeds: list[dict], *, rng: random.Random,
                      tolerance: float = 0.25) -> list[dict]:
    """One decoy panel: an F2-free approved drug of similar size for each seed drug.

    Size-matched because fragment count scales with molecular size, and a panel
    of small drugs would produce a smaller, less productive fragment pool --
    which would make the control easier to beat for a reason that has nothing to
    do with thrombin.

    Both sides are capped at `DECOY_MAX_MW`, which in practice drops bivalirudin
    (MW 2,180, a 20-residue peptide) from the matching. That is a stated
    limitation, not a silent one: BRICS on a peptide yields a pool of
    interchangeable amide fragments whose recombination generates duplicates
    indefinitely, so a bivalirudin-sized decoy is a hang rather than a control.
    The panel is therefore matched to the three small-molecule seeds.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    banned = {s["name"] for s in seeds}
    eligible = []
    for row in approved:
        if row["name"] in banned or SYMBOL in (row["targets"] or "").split(";"):
            continue
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        if mw > DECOY_MAX_MW:
            continue
        eligible.append((mw, row))

    panel, used = [], set()
    for drug in seeds:
        target_mw = Descriptors.MolWt(Chem.MolFromSmiles(drug["smiles"]))
        if target_mw > DECOY_MAX_MW:
            continue
        band = [r for mw, r in eligible
                if abs(mw - target_mw) <= tolerance * target_mw and r["name"] not in used]
        if not band:                       # nothing that size: fall back to nearest by MW
            band = [r for _, r in sorted(eligible, key=lambda t: abs(t[0] - target_mw))
                    if r["name"] not in used][:20]
        pick = rng.choice(band)
        used.add(pick["name"])
        panel.append(pick)
    return panel


# ---------------------------------------------------------------------------
# matching and the comparison
# ---------------------------------------------------------------------------


def nearest_neighbours(queries: list[str], approved: list[dict], symbol: str = SYMBOL,
                       *, drop_names: set[str] | None = None) -> list[dict]:
    """Nearest approved drug for each query by ECFP4-2048 Tanimoto.

    One fingerprint call for corpus and queries together, via
    scripts/representations.py REPRESENTATIONS['morgan'], so this arm is measured
    with the same representation as the rest of the repo.
    """
    keep = [r for r in approved if not drop_names or r["name"] not in drop_names]
    mat = unseeded_arm.fingerprints([r["smiles"] for r in keep] + queries)
    corpus_fp, query_fp = mat[: len(keep)], mat[len(keep):]

    annotated = {i for i, r in enumerate(keep)
                 if md.classify_novelty(r["targets"], r["moa_targets"], symbol) != "novel_pairing"}
    rows = []
    for smi, fp in zip(queries, query_fp):
        sims = md.tanimoto(fp, corpus_fp)
        best = int(np.argmax(sims))
        rows.append({
            "smiles": smi,
            "nearest": keep[best]["name"],
            "tanimoto": round(float(sims[best]), 4),
            "nn_is_target_annotated": best in annotated,
            "novelty": md.classify_novelty(keep[best]["targets"], keep[best]["moa_targets"], symbol),
            "best_similarity_to_a_known_binder": round(
                float(max((sims[j] for j in annotated), default=0.0)), 4),
        })
    return rows


def summarise(rows: list[dict], label: str) -> dict:
    """n, the Tanimoto distribution, and the count that actually matters."""
    sims = np.array([r["tanimoto"] for r in rows], dtype=float)
    top10 = sorted(rows, key=lambda r: -r["tanimoto"])[:10]
    return {
        "label": label,
        "n": len(rows),
        "best_nn_tanimoto": round(float(sims.max()), 4) if len(sims) else None,
        "median_nn_tanimoto": round(float(np.median(sims)), 4) if len(sims) else None,
        "mean_nn_tanimoto": round(float(sims.mean()), 4) if len(sims) else None,
        "n_nn_target_annotated": int(sum(r["nn_is_target_annotated"] for r in rows)),
        "frac_nn_target_annotated": round(
            float(np.mean([r["nn_is_target_annotated"] for r in rows])), 4) if rows else None,
        "n_target_annotated_in_top10": int(sum(r["nn_is_target_annotated"] for r in top10)),
        # WHICH annotated drugs. If they are all seed drugs, the "enrichment" is
        # the generator recognising its own parents and nothing more.
        "target_annotated_nn_names": dict(sorted(Counter(
            r["nearest"] for r in rows if r["nn_is_target_annotated"]).items(),
            key=lambda kv: -kv[1])),
        "top10": [{"tanimoto": r["tanimoto"], "nearest": r["nearest"],
                   "novelty": r["novelty"], "smiles": r["smiles"]} for r in top10],
    }


def compare(seeded: list[dict], baseline: list[dict], *, label: str) -> dict:
    """Mann-Whitney on the Tanimoto distributions, Fisher on the annotated-NN counts."""
    from scipy.stats import fisher_exact, mannwhitneyu

    a = np.array([r["tanimoto"] for r in seeded], dtype=float)
    b = np.array([r["tanimoto"] for r in baseline], dtype=float)
    u, p = mannwhitneyu(a, b, alternative="greater")
    rank_biserial = 2.0 * u / (len(a) * len(b)) - 1.0      # 0 = no difference, 1 = total

    a_hit = int(sum(r["nn_is_target_annotated"] for r in seeded))
    b_hit = int(sum(r["nn_is_target_annotated"] for r in baseline))
    odds, p_fisher = fisher_exact([[a_hit, len(a) - a_hit], [b_hit, len(b) - b_hit]],
                                  alternative="greater")
    return {
        "comparison": label,
        "n_seeded": len(a),
        "n_reference": len(b),
        "tanimoto": {
            "median_seeded": round(float(np.median(a)), 4),
            "median_reference": round(float(np.median(b)), 4),
            "mannwhitney_U": float(u),
            "p_seeded_greater": float(f"{p:.4g}"),
            "rank_biserial": round(float(rank_biserial), 4),
        },
        "nn_target_annotated": {
            "seeded": f"{a_hit}/{len(a)}",
            "reference": f"{b_hit}/{len(b)}",
            "fisher_odds_ratio": (None if not np.isfinite(odds) else round(float(odds), 4)),
            "fisher_odds_ratio_infinite": bool(not np.isfinite(odds)),
            "p_seeded_greater": float(f"{p_fisher:.4g}"),
        },
    }


def load_unseeded(path: Path = UNSEEDED_BASELINE) -> list[str]:
    """The unseeded arm's kept molecules, read from its own output. Not regenerated."""
    if not path.exists():
        raise SystemExit(
            f"unseeded baseline not found at {path}. Run demo/generate.py --stage generate "
            f"first; this module will not invent a baseline it did not read.")
    return [m["smiles"] for m in json.loads(path.read_text())["molecules"]]


# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbol", default=SYMBOL)
    parser.add_argument("--n", type=int, default=N_TARGET)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--budget", type=float, default=BUILD_BUDGET_S)
    parser.add_argument("--decoy-panels", type=int, default=N_DECOY_PANELS)
    parser.add_argument("--unseeded", type=Path, default=UNSEEDED_BASELINE)
    parser.add_argument("--out", type=Path, default=OUT_DIR)
    parser.add_argument("--reinvent-python", type=Path, default=REINVENT_PYTHON,
                        help="interpreter of the REINVENT4 venv (it pins torch==2.12, "
                             "so it cannot share env-kit)")
    parser.add_argument("--reinvent-prior", type=Path, default=LIBINVENT_PRIOR,
                        help="libinvent.prior from Zenodo record 20701824")
    parser.add_argument("--reinvent-device", default="mps")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    approved = unseeded_arm.load_approved()          # approved == 1 only: 2,153 rows
    seeds = seed_drugs(approved, args.symbol)
    if not seeds:
        raise SystemExit(f"no approved drug has {args.symbol} as a mechanism target")
    n_annot = sum(1 for r in approved
                  if md.classify_novelty(r["targets"], r["moa_targets"], args.symbol)
                  != "novel_pairing")
    print(f"[1] corpus {len(approved)} approved; {n_annot} annotated to {args.symbol} "
          f"(base rate {n_annot / len(approved):.4f})")
    print(f"    seeds: {', '.join(s['name'] for s in seeds)}")

    print(f"[2] BRICS build, quota {args.n}, budget {args.budget:.0f}s")
    gen = generate_seeded(seeds, approved, n=args.n, seed=args.seed, budget_s=args.budget)
    print(f"    pool {gen['fragment_pool']['n']} fragments -> {gen['build']['n_raw']} raw "
          f"({gen['build']['stop_reason']}, {gen['build']['seconds']}s) -> {gen['n_kept']} kept; "
          f"dropped {gen['dropped']}")
    (args.out / "10_seeded_generated.json").write_text(json.dumps(gen, indent=1))

    seeded_smiles = [m["smiles"] for m in gen["molecules"]]
    if not seeded_smiles:
        raise SystemExit("seeded generation kept nothing; there is nothing to compare")
    motif_frac = float(np.mean([bool(m["motifs"]) for m in gen["molecules"]]))
    atom_frac = [m["seed_atom_fraction"] for m in gen["molecules"]]
    print(f"    {motif_frac:.1%} carry an arginine mimetic")

    print("[3] matching seeded and unseeded against the same corpus")
    unseeded_smiles = load_unseeded(args.unseeded)
    # The circularity meter's own background: what the unseeded model scores on a
    # measure it has no way to be good at.
    cores = fragment_cores(list(gen["fragment_pool"]["fragments"]))
    background = [seed_atom_fraction(s, cores) for s in unseeded_smiles]
    print(f"    seed_atom_fraction: seeded median {np.median(atom_frac):.3f} vs "
          f"unseeded background {np.median(background):.3f} "
          f"(fragments >= {MIN_FRAGMENT_ATOMS} heavy atoms)")
    seeded_rows = nearest_neighbours(seeded_smiles, approved, args.symbol)
    unseeded_rows = nearest_neighbours(unseeded_smiles, approved, args.symbol)

    seed_names = {s["name"] for s in seeds}
    seeded_heldout = nearest_neighbours(seeded_smiles, approved, args.symbol,
                                        drop_names=seed_names)

    print(f"[4] decoy-seeded control, {args.decoy_panels} panels")
    rng = random.Random(args.seed + 1)
    panels = []
    for p in range(args.decoy_panels):
        panel = mw_matched_decoys(approved, seeds, rng=rng)
        pgen = generate_seeded(panel, approved, n=args.n, seed=args.seed + 100 + p,
                               budget_s=args.budget)
        prows = nearest_neighbours([m["smiles"] for m in pgen["molecules"]],
                                   approved, args.symbol) if pgen["molecules"] else []
        panels.append({"panel": [d["name"] for d in panel],
                       "build": pgen["build"], "n_kept": pgen["n_kept"],
                       "summary": summarise(prows, f"decoy panel {p}") if prows else
                       {"label": f"decoy panel {p}", "n": 0,
                        "status": "not_evaluated",
                        "reason": "BRICS build produced no molecule passing the gate"},
                       "_rows": prows})
        s = panels[-1]["summary"]
        print(f"    panel {p}: {', '.join(d['name'] for d in panel)} -> n={s['n']} "
              f"best={s.get('best_nn_tanimoto')} annotated_nn={s.get('n_nn_target_annotated')}")

    pooled_decoy = [r for p in panels for r in p["_rows"]]
    for p in panels:
        p.pop("_rows")

    # Step 2: the learned decoration model. Scaffold fixed, R-groups from a
    # ChEMBL prior rather than from approved drugs, so its seed_atom_fraction is
    # the direct comparison to BRICS's 1.0.
    print(f"[5] LibInvent (REINVENT4) via {args.reinvent_python}")
    lib = libinvent_generate(seeds, approved, python=args.reinvent_python,
                             prior=args.reinvent_prior, n=args.n,
                             device=args.reinvent_device,
                             workdir=args.out / "libinvent", cores=cores)
    lib_rows: list[dict] = []
    if lib.get("status") == "evaluated":
        lib_rows = nearest_neighbours([m["smiles"] for m in lib["molecules"]],
                                      approved, args.symbol)
        lib_motif = float(np.mean([bool(m["motifs"]) for m in lib["molecules"]]))
        lib_atom = [m["seed_atom_fraction"] for m in lib["molecules"]]
        lib["frac_with_arg_mimetic"] = round(lib_motif, 4)
        lib["median_seed_atom_fraction"] = round(float(np.median(lib_atom)), 4)
        print(f'    {lib["n_raw"]} raw -> {lib["n_kept"]} kept in {lib["seconds"]}s; '
              f'{lib_motif:.1%} carry an arginine mimetic; '
              f'median seed_atom_fraction {np.median(lib_atom):.3f}')
        (args.out / "12_libinvent_generated.json").write_text(json.dumps(lib, indent=1))
    else:
        print(f'    not evaluated: {lib.get("reason")}')

    payload = {
        "arm": "S-seeded (scaffold-seeded small molecule, direct matching)",
        "target": args.symbol,
        "representation": "ECFP4-2048 (scripts/representations.py REPRESENTATIONS['morgan'])",
        "corpus": {"file": "data/approved_drugs.csv", "n": len(approved),
                   "n_annotated_to_target": n_annot,
                   "base_rate": round(n_annot / len(approved), 4)},
        "seeds": [s["name"] for s in seeds],
        "circularity": {
            "statement": (
                "Fragments were taken from approved drugs and the products were matched back "
                "against approved drugs. A Tanimoto lift over the unseeded model is expected "
                "by construction and is NOT a result on its own."),
            "seed_atom_fraction": {
                "min_fragment_atoms": MIN_FRAGMENT_ATOMS,
                "seeded_median": round(float(np.median(atom_frac)), 4),
                "unseeded_background_median": round(float(np.median(background)), 4),
                "note": ("The unseeded molecules never saw a seed drug, so their median is "
                         "the floor this meter can reach. Read the gap, not the value."),
            },
            "frac_with_arg_mimetic": round(motif_frac, 4),
            "controls": ["decoy_seeded_panels", "leave_seeds_out", "libinvent"],
        },
        "seeded": summarise(seeded_rows, "seeded (F2 mechanism drugs)"),
        "unseeded": summarise(unseeded_rows, f"unseeded {json.loads(args.unseeded.read_text())['model']}"),
        "seeded_leave_seeds_out": summarise(
            seeded_heldout, "seeded, the 4 seed drugs deleted from the corpus"),
        "decoy_seeded": {
            "panels": panels,
            "pooled": summarise(pooled_decoy, "decoy-seeded, pooled") if pooled_decoy else
                      {"status": "not_evaluated", "reason": "no decoy panel produced molecules"},
        },
        "libinvent": (
            {**{k: v for k, v in lib.items() if k != "molecules"},
             "summary": summarise(lib_rows, "LibInvent scaffold decoration"),
             "leave_seeds_out": summarise(
                 nearest_neighbours([m["smiles"] for m in lib["molecules"]], approved,
                                    args.symbol, drop_names=seed_names),
                 "LibInvent, the 4 seed drugs deleted from the corpus")}
            if lib_rows else lib),
        "tests": {
            "seeded_vs_unseeded": compare(seeded_rows, unseeded_rows,
                                          label="seeded vs unseeded GPT-2"),
            "seeded_vs_decoy_seeded": (
                compare(seeded_rows, pooled_decoy,
                        label="seeded vs decoy-seeded (the circularity control)")
                if pooled_decoy else
                {"status": "not_evaluated", "reason": "no decoy panel produced molecules"}),
            "decoy_vs_unseeded": (
                compare(pooled_decoy, unseeded_rows,
                        label="decoy-seeded vs unseeded: how much lift is generic BRICS")
                if pooled_decoy else
                {"status": "not_evaluated", "reason": "no decoy panel produced molecules"}),
            "libinvent_vs_unseeded": (
                compare(lib_rows, unseeded_rows, label="LibInvent vs unseeded GPT-2")
                if lib_rows else
                {"status": "not_evaluated", "reason": lib.get("reason")}),
            "libinvent_vs_decoy_seeded": (
                compare(lib_rows, pooled_decoy,
                        label="LibInvent vs decoy-seeded (the circularity control)")
                if lib_rows and pooled_decoy else
                {"status": "not_evaluated", "reason": lib.get("reason")}),
        },
        "caveats": [
            "Tanimoto similarity is not affinity and implies nothing about binding.",
            "Nothing here is docked. The pocket enters only through the choice of seeds.",
            "BRICS recombination is enumeration, not learning: it cannot propose a "
            "substructure that was not already in a seed drug.",
            "The products are novel only in the weak sense of not being in the approved "
            "library by InChIKey; they are rearrangements of approved-drug fragments.",
        ],
        "results": sorted(seeded_rows, key=lambda r: -r["tanimoto"]),
    }
    path = args.out / "11_seeded_vs_unseeded.json"
    path.write_text(json.dumps(payload, indent=1))

    print(f"\n{'set':34s} {'n':>5s} {'best':>7s} {'median':>7s} {'annot NN':>9s} {'top10':>6s}")
    blocks = [payload["seeded"], payload["seeded_leave_seeds_out"]]
    if lib_rows:
        blocks += [payload["libinvent"]["summary"], payload["libinvent"]["leave_seeds_out"]]
    blocks += [payload["decoy_seeded"]["pooled"], payload["unseeded"]]
    for block in blocks:
        if block.get("n"):
            print(f'{block["label"][:34]:34s} {block["n"]:5d} {block["best_nn_tanimoto"]:7.4f} '
                  f'{block["median_nn_tanimoto"]:7.4f} {block["n_nn_target_annotated"]:9d} '
                  f'{block["n_target_annotated_in_top10"]:6d}')
    for key, test in payload["tests"].items():
        if "tanimoto" in test:
            print(f'  {key:26s} tanimoto p={test["tanimoto"]["p_seeded_greater"]:<10} '
                  f'rb={test["tanimoto"]["rank_biserial"]:+.3f}   '
                  f'annotated-NN p={test["nn_target_annotated"]["p_seeded_greater"]}')
    print(f"  -> {path.relative_to(REPO) if path.is_relative_to(REPO) else path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
