"""THE PIPELINE. Binding site -> designed binders -> interface signature -> approved
drugs that engage the same site.

    disease/target + pocket
        -> BoltzGen designs a few binders against that pocket      (scripts/boltzgen_signature.py)
        -> consensus INTERFACE SIGNATURE: which target residues a good binder engages
        -> co-fold each approved drug WITH the target (Boltz-2 via Rowan),
           unconstrained by default, then score the pose against that site
        -> extract each drug's target-side contacts
        -> rank drugs by how much of the signature they engage
        -> check: do the target's KNOWN binders come out on top?

This is PROJECT_GOAL.md's v0.2 pipeline. The match is computed in TARGET-SIDE
coordinates - which residues are engaged - not in chemical space, because a designed
miniprotein and a small molecule share no chemistry to compare. That is the whole
point: scripts/match_candidates.py ranks the same drugs by chemical similarity and
buries sunitinib at 464 and pazopanib at 764, both real binders of this target.

NO AFFINITY IS USED ANYWHERE IN THE RANKING. Boltz-2 can emit an affinity score and
this pipeline deliberately does not read it: PROJECT_GOAL.md 4.4 removes affinity from
the ranking path and G8 bans the vocabulary. Ranking is interface overlap, gated on
structural confidence.

Contacts for designs and for drugs come from the SAME function in interfaces.py
(PROJECT_GOAL.md F5) - if the two drifted apart, every score would be meaningless.

Phases, each resumable; state in results/pipeline/<slug>/repurpose_state.json:
    shortlist   pick the drugs to co-fold, with positive controls and decoys
    submit      send one co-folding job per drug. UNCONSTRAINED by default:
                the pocket is used to SCORE the pose afterwards, not to
                condition the folding. --constrain-pocket changes that
    collect     poll, download poses
    score       extract contacts, score against the signature, rank, validate

Usage:
    ./env/bin/python scripts/repurpose.py shortlist --n-decoys 20
    ./env/bin/python scripts/repurpose.py submit --max-credits 120
    ./env/bin/python scripts/repurpose.py collect
    ./env/bin/python scripts/repurpose.py score --designs v4
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interfaces  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEED = 0
CUTOFF = 4.5
CONFIDENCE_GATE = 0.5      # documented below; records below it are kept and flagged


def load_env():
    env = ROOT / ".env"
    if env.exists():
        for line in env.read_text().splitlines():
            if "=" in line and not line.strip().startswith("#"):
                k, _, v = line.partition("=")
                os.environ.setdefault(k.strip(), v.strip())


def paths(slug):
    p = ROOT / "results" / "pipeline" / slug
    return {"dir": p, "state": p / "repurpose_state.json",
            "poses": p / "cofold_poses", "target": p / "target"}


def state_load(slug):
    f = paths(slug)["state"]
    return json.loads(f.read_text()) if f.exists() else {}


def state_save(slug, st):
    paths(slug)["state"].write_text(json.dumps(st, indent=1))


def target_info(slug):
    t = paths(slug)["target"]
    pocket = json.loads((t / "pocket.json").read_text())
    fa = next(x for x in (t / "sequence_kinase_domain.fasta", t / "sequence.fasta") if x.exists())
    header = fa.read_text().splitlines()[0]
    seq = "".join(l.strip() for l in fa.read_text().splitlines() if not l.startswith(">"))
    # the kinase-domain fasta header records the UniProt slice it was cut from
    start = 1
    for tok in header.replace("|", " ").split():
        if "-" in tok and tok.split("-")[0].isdigit():
            start = int(tok.split("-")[0])
            break
    return {"pocket": pocket, "sequence": seq, "seq_start": start,
            "symbol": pocket.get("gene") or slug, "uniprot": pocket.get("uniprot_id")}


# --------------------------------------------------------------------------
def cmd_shortlist(args):
    """Positive controls (known binders) + property-matched decoys.

    Both are needed: the known binders say whether the ranking works at all, and the
    decoys give it something to be enriched against. Decoys are matched on MW and
    cLogP so that "looks like a drug that binds kinases" is not what is being measured.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors
    RDLogger.DisableLog("rdApp.*")

    info = target_info(args.pipeline)
    sym = args.target or info["symbol"]
    rows = [r for r in csv.DictReader((ROOT / "data" / "approved_drugs.csv").open())
            if r["smiles"]]
    known = [r for r in rows if sym in (r.get("targets") or "").split(";")]
    if not known:
        sys.exit(f"no approved drug in the corpus is annotated against {sym}")

    def props(smi):
        m = Chem.MolFromSmiles(smi)
        return (Descriptors.MolWt(m), Crippen.MolLogP(m)) if m else None

    kp = [(r, props(r["smiles"])) for r in known]
    kp = [(r, p) for r, p in kp if p]
    kp.sort(key=lambda x: -float(x[0].get("n_targets") or 0))
    positives = [r for r, _ in kp[:args.n_positives]]
    pos_props = [p for _, p in kp[:args.n_positives]]

    # decoys.
    # --hard-decoys picks APPROVED KINASE INHIBITORS that are not annotated against
    # this target (INN convention: -tinib). This is a much harder null than MW/cLogP
    # matching, which happily returns steroids and perfluorocarbons - chemically
    # matched on bulk properties but not plausible ATP-site binders, so separating
    # them proves little. A -tinib decoy binds SOME kinase ATP pocket by construction.
    # Caveat worth stating in any result: DrugCentral annotation is incomplete, so
    # "not annotated against KDR" is not "does not bind KDR" - several of these
    # (crizotinib, bosutinib) are promiscuous. A hard decoy scoring well may be a
    # labelling gap rather than a false positive.
    pool = [(r, props(r["smiles"])) for r in rows
            if sym not in (r.get("targets") or "").split(";")
            and (not args.hard_decoys or r["name"].strip().endswith("tinib"))]
    pool = [(r, p) for r, p in pool if p]
    rng = np.random.default_rng(SEED)
    chosen, used = [], set()
    for mw, lp in pos_props:
        scored = sorted(pool, key=lambda rp: abs(rp[1][0] - mw) / 100 + abs(rp[1][1] - lp))
        for r, _ in scored:
            if r["struct_id"] not in used:
                used.add(r["struct_id"])
                chosen.append(r)
                break
        if len(chosen) >= args.n_decoys:
            break
    while len(chosen) < args.n_decoys and pool:
        r, _ = pool[int(rng.integers(len(pool)))]
        if r["struct_id"] not in used:
            used.add(r["struct_id"])
            chosen.append(r)

    shortlist = ([{"struct_id": r["struct_id"], "name": r["name"], "smiles": r["smiles"],
                   "role": "known_binder"} for r in positives] +
                 [{"struct_id": r["struct_id"], "name": r["name"], "smiles": r["smiles"],
                   "role": "hard_decoy" if args.hard_decoys else "decoy"}
                  for r in chosen])
    st = state_load(args.pipeline)
    if args.append and st.get("shortlist"):
        seen = {d["struct_id"] for d in st["shortlist"]}
        shortlist = st["shortlist"] + [d for d in shortlist if d["struct_id"] not in seen]
    st.update({"target_symbol": sym, "shortlist": shortlist,
               "n_known_in_corpus": len(known),
               "decoy_matching": "nearest neighbour in (MW/100, cLogP) to each positive, "
                                 "excluding anything annotated against the target"})
    state_save(args.pipeline, st)
    print(f"{sym}: {len(known)} approved drugs annotated against it in the corpus")
    print(f"shortlist: {len(positives)} known binders + {len(chosen)} property-matched decoys "
          f"= {len(shortlist)} co-folds")
    for s in shortlist[:8]:
        print(f"   {s['role']:<12} {s['name']}")


# --------------------------------------------------------------------------
def cmd_submit(args):
    import rowan
    from stjames.workflows.protein_cofolding import PocketConstraint, Token
    load_env()
    rowan.api_key = os.environ.get("ROWAN_API_KEY") or sys.exit("no ROWAN_API_KEY in .env")

    info = target_info(args.pipeline)
    st = state_load(args.pipeline)
    if not st.get("shortlist"):
        sys.exit("run `shortlist` first")
    jobs = st.setdefault("jobs", {})
    tok = [a - info["seq_start"] for a in info["pocket"]["residue_ids"]]

    spent, submitted = 0.0, 0
    for d in st["shortlist"]:
        if d["struct_id"] in jobs:
            continue
        if spent + args.per_job_credits > args.max_credits:
            print(f"budget guard: stopping at {submitted} submissions "
                  f"(~{spent:.0f} of {args.max_credits} credits)")
            break
        # UNCONSTRAINED BY DEFAULT, and that is a measured decision that departs from
        # PROJECT_GOAL.md F3. F3 argues for a pocket constraint so a promiscuous drug
        # does not dock somewhere irrelevant. Measured on this target, the constraint
        # does the opposite damage - it forces EVERY ligand into the site, so a
        # non-binder scores like a binder and the ranking says nothing:
        #
        #                      constrained            unconstrained
        #   axitinib (binder)  16/17 pocket, 13.1 cr  16/17 pocket, 5.1 cr
        #   aspirin  (control) 11/17 pocket, 9-11 cr   6/17 pocket, 5.7 cr
        #
        # Unconstrained keeps the true binder in the pocket, halves the false signal
        # from the control, doubles the separation, and costs 2.5x less. Pass
        # --constrain-pocket to restore F3's behaviour and see this for yourself.
        pcs = None
        if args.constrain_pocket:
            pcs = [PocketConstraint(input_type="ligand", input_index=0,
                                    contacts=[Token(input_type="protein", input_index=0,
                                                    token_index=t) for t in tok],
                                    max_distance=args.pocket_distance, force=False)]
        wf = rowan.submit_protein_cofolding_workflow(
            initial_protein_sequences=[info["sequence"]],
            initial_smiles_list=[d["smiles"]],
            pocket_constraints=pcs, model="boltz_2", use_msa_server=True,
            max_credits=int(args.per_job_credits * 3),
            name=f"{st['target_symbol']} cofold - {d['name'][:40]}")
        jobs[d["struct_id"]] = {"uuid": str(wf.uuid), "name": d["name"], "role": d["role"]}
        spent += args.per_job_credits
        submitted += 1
        print(f"  submitted {d['role']:<12} {d['name'][:38]}")
    state_save(args.pipeline, st)
    print(f"\n{submitted} submitted, {len(jobs)} total tracked")


# --------------------------------------------------------------------------
def cmd_collect(args):
    import rowan
    load_env()
    rowan.api_key = os.environ.get("ROWAN_API_KEY") or sys.exit("no ROWAN_API_KEY in .env")
    st = state_load(args.pipeline)
    jobs = st.get("jobs") or sys.exit("no jobs; run `submit` first")
    out = paths(args.pipeline)["poses"]
    out.mkdir(parents=True, exist_ok=True)

    pending = [k for k, v in jobs.items() if not v.get("pose_path")]
    for attempt in range(args.max_polls):
        still = []
        for sid in pending:
            j = jobs[sid]
            wf = rowan.retrieve_workflow(j["uuid"])
            if not wf.is_finished():
                still.append(sid)
                continue
            j["status"] = str(wf.status)
            j["credits"] = wf.credits_charged
            if "OK" not in str(wf.status).upper():
                j["pose_path"] = None
                print(f"  FAILED {j['name']}: {wf.status}")
                continue
            res = (wf.data or {}).get("result") or wf.data
            uid = None
            for key in ("predicted_structure_uuid", "predicted_refined_structure_uuid"):
                uid = (res or {}).get(key) or uid
            if not uid:
                j["pose_path"] = None
                j["why"] = f"no structure uuid in result (keys {list((res or {}).keys())[:8]})"
                print(f"  NO POSE {j['name']}: {j['why']}")
                continue
            d = out / sid
            if not list(d.glob("*.pdb")):
                rowan.retrieve_protein(uid).download_pdb_file(str(d))
            files = list(d.glob("*.pdb"))
            j["pose_path"] = str(files[0].relative_to(ROOT)) if files else None
            j["scores"] = (res or {}).get("scores")
            j["constrained"] = bool((res or {}).get("pocket_constraints"))
            print(f"  got {j['role']:<12} {j['name'][:36]}  credits={wf.credits_charged}")
        pending = still
        state_save(args.pipeline, st)
        if not pending:
            break
        print(f"  {len(pending)} still running; waiting...", flush=True)
        time.sleep(args.poll_seconds)

    done = sum(1 for v in jobs.values() if v.get("pose_path"))
    total_cr = sum(v.get("credits") or 0 for v in jobs.values())
    print(f"\n{done}/{len(jobs)} poses collected, {total_cr:.1f} credits charged")


# --------------------------------------------------------------------------
def cmd_score(args):
    """Rank the co-folded drugs by how much of the design signature they engage."""
    st = state_load(args.pipeline)
    jobs = st.get("jobs") or sys.exit("no jobs; run submit/collect first")
    info = target_info(args.pipeline)
    start = info["seq_start"]

    # THREE signatures, not one. m2_gate.py measured that the BoltzGen consensus is
    # WORSE than the P2Rank pocket at locating a real ligand's contacts, so ranking
    # drugs only against the design signature would hide the cheaper, better option.
    # PROJECT_GOAL.md 4.3 makes `source` a discriminator precisely so these are swappable.
    sig_file = ROOT / "results" / f"boltzgen_signature_{args.designs}.json"
    sig = json.loads(sig_file.read_text())
    design_sets = [set(d["contacts_auth"]) for d in sig["per_design"]]
    freq = Counter(r for s in design_sets for r in s)
    n_des = len(design_sets)
    signatures = {
        "boltzgen_consensus": {
            "core": {r for r, k in freq.items() if k / n_des >= args.core_threshold},
            "weights": {r: k / n_des for r, k in freq.items()},
            "provenance": f"{n_des} BoltzGen designs, ipTM "
                          f"{min(d['iptm'] or 0 for d in sig['per_design']):.3f}-"
                          f"{max(d['iptm'] or 0 for d in sig['per_design']):.3f} "
                          f"(all below the 0.85 the plan suggests filtering at)",
        },
        "p2rank_geometry": {
            "core": set(info["pocket"]["residue_ids"]),
            "weights": {r: 1.0 for r in info["pocket"]["residue_ids"]},
            "provenance": ("P2Rank 2.5 rank-1 pocket on the ligand-stripped "
                           "structure; 0.47 s/structure amortised over a 1,531-structure batch at 12 threads; ~7 s for a single structure in isolation, JVM startup included"),
        },
    }
    klc = paths(args.pipeline)["target"] / "known_ligand_contacts.json"
    if klc.exists():
        k = json.loads(klc.read_text())
        ids = k.get("residue_ids") or k.get("known_ligand_residue_ids") or []
        if ids:
            signatures["known_ligand"] = {
                "core": {int(x) for x in ids},
                "weights": {int(x): 1.0 for x in ids},
                "provenance": "contacts of the ligand in the primary holo structure; "
                              "VALIDATION REFERENCE ONLY - it is the answer for one ligand",
            }
    for name, sg in signatures.items():
        print(f"signature {name:<20} {len(sg['core']):3d} residues  ({sg['provenance']})")
    print()
    core = signatures["boltzgen_consensus"]["core"]
    weights = signatures["boltzgen_consensus"]["weights"]

    rows = []
    for sid, j in jobs.items():
        if not j.get("pose_path"):
            rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                         "status": "no_pose", "why": j.get("why") or j.get("status")})
            continue
        path = ROOT / j["pose_path"]
        try:
            struct = interfaces.load_structure(path, sid)
            ligs = interfaces.drug_like_ligands(struct)
            if not ligs:
                rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                             "status": "no_ligand_in_pose"})
                continue
            lg = ligs[0]
            con = interfaces.ligand_contacts(struct, lg["resname"], chain=lg.get("chain_id"),
                                             cutoff=CUTOFF)
            engaged = {int(r) + start - 1 for r in con["residue_ids"]}
        except Exception as exc:                        # noqa: BLE001
            rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                         "status": "extract_failed", "why": repr(exc)[:120]})
            continue
        row = {"struct_id": sid, "name": j["name"], "role": j["role"],
               "status": "scored", "n_engaged": len(engaged),
               "engaged_residues": sorted(engaged), "by_signature": {}}
        for sname, sg in signatures.items():
            c, w = sg["core"], sg["weights"]
            hit = engaged & c
            num = sum(w.get(r, 0) for r in engaged & set(w))
            den = sum(w.values()) + len(engaged - set(w))
            # FIVE scores, because a bigger ligand touches more residues and so
            # scores higher for reasons that have nothing to do with binding -
            # PROJECT_GOAL.md 1.4a. E4 asks for at least two normalisations
            # compared and the default justified. Measured on this screen
            # (correlation of the score with the number of residues engaged,
            # and where a 39-residue lipopeptide decoy lands out of 29):
            #   core_coverage        r=+0.641  decoy rank  1
            #   weighted_jaccard     r=+0.520  decoy rank  2
            #   f1_engaged_core      r=+0.413  decoy rank  6
            #   jaccard_engaged_core r=+0.409  decoy rank  6
            #   precision_in_core    r=+0.135  decoy rank 14   <- default
            # Enrichment and median known-binder rank are IDENTICAL under all
            # five, so the separation is not an artefact of this choice; the
            # metric only decides where the oversized decoy lands.
            row["by_signature"][sname] = {
                "core_coverage": len(hit) / len(c) if c else float("nan"),
                "weighted_jaccard": num / den if den else 0.0,
                "f1_engaged_core": (2 * len(hit) / (len(engaged) + len(c))
                                    if (engaged or c) else 0.0),
                "jaccard_engaged_core": (len(hit) / len(engaged | c)
                                         if (engaged | c) else 0.0),
                "precision_in_core": len(hit) / len(engaged) if engaged else 0.0,
                "engaged_core": sorted(hit), "missed_core": sorted(c - engaged)}
        row["core_coverage"] = row["by_signature"]["boltzgen_consensus"]["core_coverage"]
        row["weighted_jaccard"] = row["by_signature"]["boltzgen_consensus"]["weighted_jaccard"]
        rows.append(row)

    scored = [r for r in rows if r["status"] == "scored"]
    scored.sort(key=lambda r: -r["by_signature"][args.rank_by][args.metric])
    for i, r in enumerate(scored, 1):
        r["rank"] = i
    print(f"{len(scored)}/{len(rows)} drugs scored "
          f"({sum(1 for r in rows if r['status'] != 'scored')} without a usable pose)\n")
    print(f"{'rank':>4}  {'drug':<34}{'role':<13}"
          f"{args.metric[:9]:>9}{'core_cov':>10}{'engaged':>8}")
    print("-" * 80)
    for r in scored[:30]:
        b = r["by_signature"][args.rank_by]
        print(f"{r['rank']:>4}  {r['name'][:34]:<34}{r['role']:<13}"
              f"{b[args.metric]:9.3f}{b['core_coverage']:10.3f}{r['n_engaged']:8d}")

    # the same validation for EVERY signature, so the cheap arm is not hidden
    print(f"\n{'signature':<22}{'enrichment@25%':>16}{'median rank of a known binder':>32}")
    print("-" * 72)
    per_sig = {}
    for sname in signatures:
        order = sorted(scored, key=lambda r: -r["by_signature"][sname][args.metric])
        rk = [i for i, r in enumerate(order, 1) if r["role"] == "known_binder"]
        npos = len(rk)
        kk = max(1, round(0.25 * len(order)))
        hh = sum(1 for r in order[:kk] if r["role"] == "known_binder")
        bb = npos / len(order) if order else float("nan")
        e = (hh / kk) / bb if bb else float("nan")
        per_sig[sname] = {"enrichment_at_top_quartile": e, "known_binder_ranks": rk,
                          "median_rank": float(np.median(rk)) if rk else None}
        print(f"{sname:<22}{e:16.2f}{(np.median(rk) if rk else float('nan')):32.1f}")

    # validation: do the known binders come out on top?
    n_pos = sum(1 for r in scored if r["role"] == "known_binder")
    k = max(1, round(0.25 * len(scored)))
    hits = sum(1 for r in scored[:k] if r["role"] == "known_binder")
    base = n_pos / len(scored) if scored else float("nan")
    ef = (hits / k) / base if base else float("nan")
    ranks = [r["rank"] for r in scored if r["role"] == "known_binder"]
    print(f"\nVALIDATION: {n_pos} known binders, {len(scored) - n_pos} decoys")
    print(f"  known-binder ranks: {ranks}")
    print(f"  top {k}: {hits}/{k} are known binders (base rate {base:.3f}) "
          f"-> enrichment {ef:.2f}x")
    if ranks:
        print(f"  median rank of a known binder: {int(np.median(ranks))} of {len(scored)}")

    out = {
        "pipeline": "binding site -> BoltzGen designs -> interface signature -> "
                    "co-folded approved drugs -> ranked by interface overlap",
        "ranking_metric": args.metric,
        "ranking_metric_choice": ("size-normalised; PROJECT_GOAL.md 1.4a warns raw "
                                  "overlap ranks bigger ligands higher for non-binding "
                                  "reasons. All five scores are stored per drug so the "
                                  "choice can be checked, and enrichment is identical "
                                  "under all of them on this screen."),
        "affinity_used": False,
        "affinity_note": "Boltz-2 can emit an affinity score; it is deliberately not read. "
                         "PROJECT_GOAL.md 4.4 removes affinity from the ranking path.",
        "ranked_by": args.rank_by,
        "signatures": {k: {"n_core": len(v["core"]), "provenance": v["provenance"],
                           "core_residues_auth": sorted(v["core"])}
                       for k, v in signatures.items()},
        "validation_by_signature": per_sig,
        "target": {"symbol": st.get("target_symbol"), "uniprot": info["uniprot"],
                   "structure": info["pocket"].get("structure_id")},
        "n_scored": len(scored), "n_known_binders": n_pos,
        "enrichment_at_top_quartile": ef, "known_binder_ranks": ranks,
        "results": rows,
    }
    dest = ROOT / "results" / f"repurpose_{args.pipeline}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pipeline", default="colorectal-cancer")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("shortlist"); s.set_defaults(fn=cmd_shortlist)
    s.add_argument("--target"); s.add_argument("--n-positives", type=int, default=10)
    s.add_argument("--n-decoys", type=int, default=20)
    s.add_argument("--hard-decoys", action="store_true",
                   help="use approved kinase inhibitors (-tinib) that do not hit the "
                        "target, instead of MW/cLogP-matched drugs")
    s.add_argument("--append", action="store_true",
                   help="add to the existing shortlist rather than replacing it")

    s = sub.add_parser("submit"); s.set_defaults(fn=cmd_submit)
    s.add_argument("--max-credits", type=float, default=120)
    s.add_argument("--per-job-credits", type=float, default=6.0,
                   help="budget guard; measured ~5.4 credits per unconstrained co-fold")
    s.add_argument("--pocket-distance", type=float, default=6.0)
    s.add_argument("--constrain-pocket", action="store_true",
                   help="restore PROJECT_GOAL.md F3's pocket constraint. Measured to "
                        "destroy discrimination on this target - see the comment in "
                        "cmd_submit before using it")

    s = sub.add_parser("collect"); s.set_defaults(fn=cmd_collect)
    s.add_argument("--poll-seconds", type=int, default=30)
    s.add_argument("--max-polls", type=int, default=120)

    s = sub.add_parser("score"); s.set_defaults(fn=cmd_score)
    s.add_argument("--designs", default="v4")
    s.add_argument("--core-threshold", type=float, default=0.6)
    s.add_argument("--rank-by", default="boltzgen_consensus",
                   choices=["boltzgen_consensus", "p2rank_geometry", "known_ligand"])
    s.add_argument("--metric", default="precision_in_core",
                   choices=["precision_in_core", "core_coverage", "weighted_jaccard",
                            "f1_engaged_core", "jaccard_engaged_core"],
                   help="size-normalised by default; see the comment in cmd_score")

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
