"""Arm S: generate small molecules, dock them, find the approved drug they resemble.

    generate           GPT-2 trained on ZINC                  [local, MIT]
      -> filter        RDKit validity + drug-likeness          [local]
      -> dock          Rowan batch docking into the site       [Rowan]
      -> match         ECFP4 Tanimoto vs 4,099 approved drugs  [local]
      -> null, novelty, three leaderboards                     [demo/match_direct.py]

**What this is, stated plainly.** The generator is `entropy/gpt2_zinc_87m`, a
SMILES language model. It is NOT conditioned on the binding site -- it samples
drug-like chemical space, and the pocket only enters afterwards, as a docking
filter. So this arm is virtual screening with a generative front end, not
structure-based de novo design, and it must not be described as the latter.

That was a forced choice, and the reason is worth recording: every
pocket-conditioned generator (DiffSBDD, Pocket2Mol, TargetDiff, PocketFlow,
DecompDiff, LiGAN) depends on `torch-scatter`/`torch-sparse`/`torch-cluster`,
which have no build for torch 2.14 and have never shipped a native macOS arm64
wheel. Rowan hosts no de novo generator either -- its 41 workflows all take
molecules as input, and `analogue_docking` requires you to supply the analogues.
See docs/10-DIRECT-MATCHING.md.

Unlike arm P, the comparison here is small-molecule to small-molecule, so ECFP4
Tanimoto applies -- the representation this repo's own benchmark found best of
16 (README result #1).

Run:  DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.generate --stage all
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

from demo import match_direct as md, nova

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

APPROVED = REPO / "data" / "approved_drugs.csv"
MODEL = "entropy/gpt2_zinc_87m"          # MIT, GPT-2 87M trained on ZINC drug-like
N_GENERATE = int(os.environ.get("DEMO_N_GENERATE", "600"))
N_DOCK = int(os.environ.get("DEMO_N_DOCK", "24"))
SEED = 0

# Drug-likeness gate. Deliberately loose -- this is a plausibility filter, not a
# selection criterion, and a tight one would just reshape the library into
# whatever it was tuned on.
MW_RANGE = (200.0, 600.0)
MIN_QED = 0.4


def generate(n: int, seed: int = SEED) -> list[str]:
    """Sample SMILES from the ZINC-trained language model."""
    import torch
    from transformers import AutoTokenizer, GPT2LMHeadModel

    torch.manual_seed(seed)
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = GPT2LMHeadModel.from_pretrained(MODEL).eval()

    out: list[str] = []
    batch = 64
    while len(out) < n:
        with torch.no_grad():
            ids = model.generate(
                do_sample=True, max_length=64, num_return_sequences=batch,
                top_k=30, temperature=1.0,
                bos_token_id=tok.bos_token_id, eos_token_id=tok.eos_token_id,
                pad_token_id=tok.pad_token_id or tok.eos_token_id,
            )
        out.extend(s.strip() for s in tok.batch_decode(ids, skip_special_tokens=True))
    return out[:n]


def drug_like(smiles: list[str], approved_inchikeys: set[str]) -> tuple[list[dict], dict]:
    """Keep valid, drug-like, novel molecules. Report every reason for dropping one."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
    kept: list[dict] = []
    seen: set[str] = set()
    reasons = {"invalid": 0, "duplicate": 0, "mw": 0, "qed": 0, "already_approved": 0}

    for raw in smiles:
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            reasons["invalid"] += 1
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen:
            reasons["duplicate"] += 1
            continue
        mw = Descriptors.MolWt(mol)
        if not (MW_RANGE[0] <= mw <= MW_RANGE[1]):
            reasons["mw"] += 1
            continue
        qed = QED.qed(mol)
        if qed < MIN_QED:
            reasons["qed"] += 1
            continue
        key = rdMolDescriptors.CalcInchiKey(mol) if hasattr(rdMolDescriptors, "CalcInchiKey") else Chem.MolToInchiKey(mol)
        if key in approved_inchikeys:
            # Regenerating an approved drug outright is not a repurposing hit,
            # it is the model reciting its training set.
            reasons["already_approved"] += 1
            continue
        seen.add(canonical)
        kept.append({"smiles": canonical, "mw": round(mw, 1), "qed": round(qed, 3), "inchikey": key})

    return kept, reasons


def docking_box(pdb: Path, chain: str, positions: list[int], pad: float = 4.0) -> list[list[float]]:
    """Centre and size of a box around the site residues, for Vina."""
    from demo import contacts

    protein, _, _ = contacts.parse_pdb(pdb)
    ordered = [lbl for (ch, lbl) in sorted(
        (k for k in protein if k[0] == chain), key=lambda k: contacts.residue_order(k[1])
    )]
    wanted = {ordered[i - 1] for i in positions if 1 <= i <= len(ordered)}
    coords = np.vstack([c for (ch, lbl), (_, c) in protein.items() if ch == chain and lbl in wanted])
    lo, hi = coords.min(0), coords.max(0)
    centre = (lo + hi) / 2.0
    size = (hi - lo) + 2 * pad
    return [[round(float(x), 3) for x in centre], [round(float(x), 3) for x in size]]


def load_approved() -> list[dict]:
    """Approved small molecules with a SMILES, plus their target annotations."""
    rows = []
    for row in csv.DictReader(APPROVED.open()):
        if not row.get("smiles"):
            continue
        rows.append(
            {
                "name": row["name"],
                "smiles": row["smiles"],
                "inchikey": row.get("inchikey") or "",
                "targets": row.get("targets") or "",
                "moa_targets": row.get("moa_targets") or "",
            }
        )
    return rows


def fingerprints(smiles: list[str]) -> np.ndarray:
    """ECFP4-2048 via the repo's own registry, so the arm matches the benchmark."""
    import representations  # scripts/representations.py

    return representations.REPRESENTATIONS["morgan"](smiles)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="all", choices=["generate", "dock", "match", "all"])
    args = parser.parse_args()

    from demo import pipeline  # after DEMO_TARGET is set

    symbol = pipeline.TARGET["symbol"]
    out = pipeline.OUT
    out.mkdir(parents=True, exist_ok=True)

    if args.stage in ("generate", "all"):
        approved = load_approved()
        keys = {r["inchikey"] for r in approved if r["inchikey"]}
        print(f"[S1] generating {N_GENERATE} molecules with {MODEL}")
        raw = generate(N_GENERATE)
        kept, reasons = drug_like(raw, keys)
        print(f"  {len(kept)}/{len(raw)} kept; dropped {reasons}")
        (out / "07_generated.json").write_text(json.dumps(
            {"model": MODEL, "n_requested": N_GENERATE, "n_kept": len(kept),
             "filters": {"mw_range": MW_RANGE, "min_qed": MIN_QED},
             "dropped": reasons, "molecules": kept}, indent=1))
        print(f"  -> {(out / '07_generated.json').relative_to(REPO)}")

    if args.stage in ("dock", "all"):
        gen = json.loads((out / "07_generated.json").read_text())
        site = json.loads((out / "02_site.json").read_text())
        pdb = pipeline._download_structure(site["protein_uuid"], f'target_{pipeline.TARGET["out"]}')
        box = docking_box(pdb, site["target_chain"], site["chosen"]["residue_positions_1based"])
        print(f"[S2] docking {N_DOCK} molecules; box centre {box[0]} size {box[1]}")

        # Dock the generated molecules AND the known binders, in one batch, so the
        # positive control is scored by exactly the same run.
        approved = {r["name"]: r for r in load_approved()}
        controls = [n for n in ("argatroban", "ximelagatran", "dabigatran etexilate") if n in approved]
        smiles = [m["smiles"] for m in gen["molecules"][:N_DOCK]] + [approved[n]["smiles"] for n in controls]

        record = nova.run_workflow(
            f'dock_{pipeline.TARGET["out"]}', "submit_batch_docking_workflow",
            label=f"batch dock x{len(smiles)}", max_credits=120,
            smiles_list=smiles, protein=site["protein_uuid"], pocket=box,
            name=f"{symbol} arm S batch docking",
        )
        scores = record["data"].get("best_scores") or record["data"].get("scores")
        (out / "08_docked.json").write_text(json.dumps(
            {"workflow_uuid": record["uuid"], "box": box, "controls": controls,
             "n_generated_docked": len(smiles) - len(controls),
             "smiles": smiles, "scores": scores,
             "credits_charged": record.get("credits_charged")}, indent=1, default=str))
        print(f"  -> {(out / '08_docked.json').relative_to(REPO)}  credits {record.get('credits_charged')}")

    if args.stage in ("match", "all"):
        run_match(out, symbol)


def run_match(out: Path, symbol: str) -> None:
    """ECFP4 nearest approved drug for each generated molecule, with the null."""
    gen = json.loads((out / "07_generated.json").read_text())
    docked_path = out / "08_docked.json"
    docked = json.loads(docked_path.read_text()) if docked_path.exists() else None

    approved = load_approved()
    for row in approved:
        row["novelty"] = md.classify_novelty(row["targets"], row["moa_targets"], symbol)
    n_known = sum(1 for r in approved if r["novelty"] != "novel_pairing")
    print(f"[S3] corpus {len(approved)} approved molecules; {n_known} annotated to {symbol}")

    molecules = gen["molecules"]
    queries = [m["smiles"] for m in molecules]
    mat = fingerprints([r["smiles"] for r in approved] + queries)
    corpus_fp, query_fp = mat[: len(approved)], mat[len(approved):]

    # The null must be INDEPENDENT of the molecules being scored. Scoring every
    # generated molecule against a null built from those same molecules only
    # recovers its rank within its own set, which is not a null at all.
    #
    # So the split is by pocket evidence: molecules the docking filter selected
    # are the queries, the rest are the null. Both come from the same generator
    # and the same drug-likeness filter, so the only difference between the two
    # groups is whether the pocket chose them. The percentile then answers a
    # real question -- do pocket-selected molecules resemble approved drugs any
    # more closely than pocket-agnostic ones?
    docked_smiles = set((docked or {}).get("smiles") or [])
    controls = set((docked or {}).get("controls") or [])
    is_query = [m["smiles"] in docked_smiles for m in molecules]
    if not any(is_query):
        null_idx = list(range(len(molecules)))
        query_idx = list(range(len(molecules)))
        null_note = (
            "DEGENERATE: docking has not run, so queries and null are the same set "
            "and the percentiles below are ranks within that set, not a null test."
        )
    else:
        query_idx = [i for i, q in enumerate(is_query) if q]
        null_idx = [i for i, q in enumerate(is_query) if not q]
        null_note = (
            f"{len(null_idx)} generated molecules not in the docked subset, scored "
            f"against the same corpus. The docked subset was taken by list position "
            f"rather than by fit, so this null asks whether a query beats a TYPICAL "
            f"generated molecule -- not whether the pocket chose it. See docking_signal."
        )

    null = md.null_distribution(query_fp[null_idx], corpus_fp, "tanimoto")

    rows = []
    for i in query_idx:
        mol, fp = molecules[i], query_fp[i]
        sims = md.tanimoto(fp, corpus_fp)
        best = int(np.argmax(sims))
        rows.append({
            "smiles": mol["smiles"], "mw": mol["mw"], "qed": mol["qed"],
            "nearest": approved[best]["name"],
            "tanimoto": round(float(sims[best]), 4),
            "novelty": approved[best]["novelty"],
            "nearest_targets": approved[best]["targets"][:60],
            "null_percentile": round(md.percentile_of(float(sims[best]), null), 4),
            "best_similarity_to_a_known_binder": round(float(
                max((sims[j] for j, r in enumerate(approved) if r["novelty"] != "novel_pairing"),
                    default=0.0)), 4),
        })
    rows.sort(key=lambda r: -r["tanimoto"])

    # Does the POCKET add anything? The docked subset was taken by list position,
    # not by fit, so "docked vs undocked" is not a pocket-evidence split and the
    # percentile above only asks whether a molecule beats a typical generated one.
    # The pocket question is separate and is answered here: among the docked
    # molecules, does a better Vina score go with greater similarity to the drugs
    # already known to bind this target? If it does not, docking is selecting for
    # something unrelated to the known chemistry of this site.
    docking_signal = {"status": "not_evaluated", "reason": "docking has not run"}
    if docked and docked.get("scores"):
        from scipy.stats import spearmanr

        n_gen = docked["n_generated_docked"]
        by_smiles = {m["smiles"]: i for i, m in enumerate(molecules)}
        pairs = []
        for smi, score in list(zip(docked["smiles"], docked["scores"]))[:n_gen]:
            if score is None or smi not in by_smiles:
                continue
            sims = md.tanimoto(query_fp[by_smiles[smi]], corpus_fp)
            best_known = max(
                (sims[j] for j, r in enumerate(approved) if r["novelty"] != "novel_pairing"),
                default=0.0,
            )
            pairs.append((float(score), float(best_known)))

        controls_scores = [s for s in docked["scores"][n_gen:] if s is not None]

        # The pocket-evidence test proper: split the docked molecules into the
        # ones that dock WELL and the ones that dock badly, and ask whether the
        # good dockers resemble approved drugs -- and known binders in
        # particular -- any more than the bad ones. Unlike a split by list
        # position, this one is actually about the pocket.
        if len(pairs) >= 20:
            from scipy.stats import mannwhitneyu

            ordered = sorted(pairs, key=lambda p: p[0])          # most negative = best
            cut = max(5, len(ordered) // 4)
            good, poor = ordered[:cut], ordered[-cut:]
            good_known = [k for _, k in good]
            poor_known = [k for _, k in poor]
            u, pv = mannwhitneyu(good_known, poor_known, alternative="greater")
            quartile_test = {
                "n_per_group": cut,
                "mean_similarity_to_known_binders": {
                    "good_dockers": round(float(np.mean(good_known)), 4),
                    "poor_dockers": round(float(np.mean(poor_known)), 4),
                },
                "mannwhitney_U": float(u),
                "p_good_gt_poor": float(f"{pv:.4g}"),
                "score_ranges": {
                    "good_dockers": [round(good[0][0], 3), round(good[-1][0], 3)],
                    "poor_dockers": [round(poor[0][0], 3), round(poor[-1][0], 3)],
                },
            }
        else:
            quartile_test = {
                "status": "not_evaluated",
                "reason": f"only {len(pairs)} molecules docked; need 20 for quartiles",
            }

        if len(pairs) >= 5:
            scores, sims_known = zip(*pairs)
            rho, pval = spearmanr(scores, sims_known)
            docking_signal = {
                "status": "evaluated",
                "question": "does a better Vina score go with similarity to known binders?",
                "n": len(pairs),
                "spearman_rho": round(float(rho), 4),
                "p_value": float(f"{pval:.4g}"),
                "note": "Vina scores are negative-is-better, so a NEGATIVE rho means "
                        "better docking goes with greater similarity to known binders",
                "generated_score_range": [round(min(scores), 3), round(max(scores), 3)],
                "control_scores": {
                    name: score for name, score in zip(docked["controls"], controls_scores)
                },
                "controls_beat_generated_median": (
                    bool(np.median(controls_scores) < np.median(scores))
                    if controls_scores else None
                ),
                "good_vs_poor_dockers": quartile_test,
            }

    # Positive control: a known binder as the query must retrieve its analogues.
    # Positive control: where do the OTHER drugs with this target as their
    # mechanism rank when one of them is the query?
    #
    # The top-5 list alone is misleading here and nearly produced a false
    # negative: argatroban's five nearest approved drugs are all peptidomimetics
    # (angiotensin II, icatibant, lisinopril ...) and none is a thrombin drug,
    # which reads as failure. The rank of the co-mechanism drugs is the
    # informative number, and it says the opposite.
    control = {}
    by_name = {r["name"]: i for i, r in enumerate(approved)}
    moa_idx = [i for i, r in enumerate(approved) if symbol in r["moa_targets"].split(";")]
    for i in moa_idx:
        sims = md.tanimoto(corpus_fp[i], corpus_fp)
        sims[i] = -np.inf
        order = list(np.argsort(-sims))
        control[approved[i]["name"]] = {
            "top_5_any": [
                {"name": approved[j]["name"], "tanimoto": round(float(sims[j]), 4)}
                for j in order[:5]
            ],
            "ranks_of_other_moa_drugs": sorted(
                (
                    {
                        "name": approved[j]["name"],
                        "rank": order.index(j) + 1,
                        "of": len(order),
                        "tanimoto": round(float(sims[j]), 4),
                    }
                    for j in moa_idx
                    if j != i
                ),
                key=lambda d: d["rank"],
            ),
        }
    control["_n_moa_drugs"] = len(moa_idx)
    control["_expected_rank_by_chance"] = round(len(approved) / 2)

    payload = {
        "arm": "S (small molecule, direct matching)",
        "target": symbol,
        "framing": (
            "virtual screening with a generative front end: the generator is NOT "
            "conditioned on the binding site; the pocket enters only as a docking filter"
        ),
        "generator": gen["model"],
        "n_generated_kept": len(molecules),
        "corpus": {"file": str(APPROVED.relative_to(REPO)), "n": len(approved),
                   "n_annotated_to_target": n_known,
                   "base_rate": round(n_known / len(approved), 4)},
        "representation": "ECFP4-2048 (scripts/representations.py REPRESENTATIONS['morgan'])",
        "null": {"construction": null_note, **md.calibration(null)},
        "n_queries": len(query_idx),
        "positive_control": control,
        "docking_signal": docking_signal,
        "docking": {"ran": docked is not None,
                    "workflow_uuid": (docked or {}).get("workflow_uuid"),
                    "box": (docked or {}).get("box")},
        "results": rows,
        "caveats": [
            "The generator is unconditioned; any pocket relevance comes from the "
            "docking filter, not from generation.",
            "A nearest neighbour exists for every molecule by construction. Only "
            "the null percentile distinguishes a real resemblance from the "
            "background rate of the 4,099-molecule library.",
            "Tanimoto similarity is not affinity and implies nothing about binding.",
        ],
    }
    path = out / "09_smallmol_match.json"
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"  -> {path.relative_to(REPO)}")
    print(f'  null: mean={payload["null"]["mean"]} p95={payload["null"]["p95"]}')
    print(f'\n{"tanimoto":>9s} {"pctile":>7s} {"nearest approved":26s} {"novelty":16s} smiles')
    for r in rows[:10]:
        print(f'{r["tanimoto"]:9.4f} {r["null_percentile"]:7.3f} {r["nearest"][:26]:26s} '
              f'{r["novelty"]:16s} {r["smiles"][:44]}')


if __name__ == "__main__":
    main()
