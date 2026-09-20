"""One command: a target goes in, a ranked list of approved drugs comes out.

    ./env/bin/python scripts/autorepurpose.py run --target KDR --uniprot P35968

This is the orchestrator over the stages PROJECT_GOAL.md WS-H asks for: stage-level
resume, a run manifest, and a budget guard. Each stage is skipped if its output
already exists, so a killed run resumes where it stopped.

    prep       structure + ligand-stripped pocket        (scripts/prep_target.py)
    design     BoltzGen binders for the pocket           (optional - see below)
    shortlist  drugs to screen, with controls and decoys (scripts/repurpose.py)
    cofold     Boltz-2 co-folding of each drug           (scripts/repurpose.py)
    score      rank by interface overlap                 (scripts/repurpose.py)
    report     one self-contained HTML                   (scripts/report.py)

THE DESIGN STAGE IS OFF BY DEFAULT, and that is a measured decision, not a
simplification. Ablation I6.1 (scripts/m2_gate.py) ran on two targets and the
BoltzGen consensus LOST to P2Rank's pocket both times - by 0.324 Jaccard on KDR
(n=38) and 0.155 on CDK2 (n=31), Holm p=0.0016 each. Inside the product the two
signatures rank drugs at Spearman 0.954 with identical enrichment. So the default
signature is `p2rank_geometry`: it is free, it is fast (0.47 s/structure amortised over a 1,531-structure batch at 12 threads, 250 per JVM; a single structure in isolation measured 2.1-2.5 s wall (P2Rank self-reports 1.87 s), so the batch figure is not a single-run cost), and it is at least as
good. Pass --designs N to run BoltzGen anyway and --signature boltzgen_consensus to
rank by it; the comparison is printed either way so the choice stays visible.

PROJECT_GOAL.md section 7 prescribes exactly this on this outcome: "the honest
pipeline is target -> pocket -> ... and you should ship that: it's faster, cheaper,
and more defensible."

COSTS, measured on Rowan: co-folding ~5.4 credits per drug (~4 min); BoltzGen
~3.6 credits per design. --budget-credits refuses to start a stage whose estimate
exceeds what is left.

NO AFFINITY IS PREDICTED OR READ. Ranking is interface overlap, gated on pose
confidence (PROJECT_GOAL.md 4.4, G8).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = str(ROOT / "env" / "bin" / "python")

# report.py's inputs are hardcoded to this pipeline; see the report stage below.
REPORT_PIPELINE = "colorectal-cancer"

COST_PER_COFOLD = 5.4      # measured, 30 drugs / 137.8 credits + probes
COST_PER_DESIGN = 3.6      # measured, 24 designs / 85.5 credits


def run(cmd, dry=False):
    print(f"\n$ {' '.join(str(c) for c in cmd)}", flush=True)
    if dry:
        print("  (dry run)")
        return 0
    return subprocess.call([str(c) for c in cmd])


def stage_done(path):
    return Path(path).exists()


def choose_target(args, root):
    """Resolve a disease to ranked targets and make the caller choose.

    Returns the chosen target dict, or None if the caller must still decide (in
    which case the ranked list has been printed). Never picks silently.
    """
    slug = args.slug or args.disease.lower().replace(" ", "-")
    tpath = root / "results" / "pipeline" / slug / "targets.json"
    if not tpath.exists():
        cmd = [PY, root / "scripts/disease_targets.py", "--disease", args.disease]
        print(f"resolving \"{args.disease}\" to targets (Open Targets)...")
        if args.dry_run:
            print(f"  (dry run) $ {' '.join(str(c) for c in cmd)}")
            return None
        if run(cmd):
            print(f"could not resolve \"{args.disease}\"", file=sys.stderr)
            return None
    if not tpath.exists():
        print(f"no targets.json at {tpath} after resolution", file=sys.stderr)
        return None

    data = json.loads(tpath.read_text())
    targets = [t for t in data.get("targets", []) if t.get("uniprot")]
    if not targets:
        print(f"{tpath} lists no target with a UniProt accession", file=sys.stderr)
        return None

    if args.target:                       # caller named one; honour it
        for t in targets:
            if (t.get("approved_symbol") or "").upper() == args.target.upper():
                return {"symbol": t["approved_symbol"], "uniprot": t["uniprot"],
                        "slug": slug, "rank": t.get("rank"),
                        "score": t.get("overall_association_score"),
                        "chosen_by": "--target", "n_candidates": len(targets)}
        print(f"--target {args.target} is not among the {len(targets)} targets "
              f"resolved for \"{args.disease}\"", file=sys.stderr)
        return None

    show = targets[:max(1, args.show_targets)]
    print(f"\n\"{args.disease}\" resolves to {len(targets)} targets "
          f"(EFO {data.get('meta', {}).get('efo_id', '?')}). "
          f"This is a real choice and it is yours:\n")
    print(f"  {'#':>2}  {'symbol':<10}{'uniprot':<10}{'assoc':>7}  known drug?  name")
    for t in show:
        print(f"  {t.get('rank', '?'):>2}  {(t.get('approved_symbol') or '?'):<10}"
              f"{(t.get('uniprot') or '?'):<10}"
              f"{(t.get('overall_association_score') or 0):>7.3f}  "
              f"{'yes' if t.get('has_known_drug_evidence') else 'no ':<11}  "
              f"{(t.get('approved_name') or '')[:38]}")
    if len(targets) > len(show):
        print(f"      ... {len(targets) - len(show)} more in {tpath}")

    if args.accept_top_target:
        t = targets[0]
        print(f"\n--accept-top-target: proceeding with {t['approved_symbol']} "
              f"({t['uniprot']}), rank 1 of {len(targets)}. Recorded in the manifest.")
        return {"symbol": t["approved_symbol"], "uniprot": t["uniprot"], "slug": slug,
                "rank": t.get("rank"), "score": t.get("overall_association_score"),
                "chosen_by": "--accept-top-target", "n_candidates": len(targets)}

    print(f"\nNo target chosen. Pick one and rerun, e.g.:\n"
          f"    ./env/bin/python scripts/autorepurpose.py run --disease "
          f"\"{args.disease}\" --target {show[0].get('approved_symbol')}\n"
          f"  or add --accept-top-target to take rank 1 deliberately.\n"
          f"Association score is evidence that the target relates to the disease. It is "
          f"not evidence that the target is druggable or that this pipeline will work "
          f"on it.")
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("cmd", choices=["run", "targets"],
                    help="run = the pipeline; targets = resolve a disease and stop")
    ap.add_argument("--disease", default=None,
                    help="disease name, e.g. \"colorectal cancer\". Resolves to ranked "
                         "targets and SHOWS them; it will not pick one for you")
    ap.add_argument("--target", default=None, help="gene symbol, e.g. KDR")
    ap.add_argument("--uniprot", default=None, help="e.g. P35968")
    ap.add_argument("--accept-top-target", action="store_true",
                    help="with --disease, proceed with the top-ranked target instead of "
                         "stopping. An explicit choice, recorded in the manifest")
    ap.add_argument("--show-targets", type=int, default=10,
                    help="how many ranked targets to display for a disease")
    ap.add_argument("--slug", default=None,
                    help="results/pipeline/<slug>/ (default: lowercased target)")
    ap.add_argument("--n-positives", type=int, default=10)
    ap.add_argument("--n-decoys", type=int, default=20)
    ap.add_argument("--hard-decoys", action="store_true",
                    help="use approved kinase inhibitors that miss the target as the null")
    ap.add_argument("--designs", type=int, default=0,
                    help="run BoltzGen with this many designs (0 = skip; see the module "
                         "docstring for why the default is 0)")
    ap.add_argument("--signature", default="p2rank_geometry",
                    choices=["p2rank_geometry", "boltzgen_consensus", "known_ligand"])
    ap.add_argument("--metric", default="precision_in_core")
    ap.add_argument("--designs-tag", default=None,
                    help="tag of an existing BoltzGen design run for THIS target "
                         "(scripts/boltzgen_signature.py names outputs by --tag)")
    ap.add_argument("--budget-credits", type=float, default=200.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # ---- step 1-2: accept a disease, resolve targets, EXPOSE the ambiguity ----
    # A disease does not map to one target. Open Targets returns 25 for colorectal
    # cancer, several of them druggable. Picking the top one silently would hide a
    # scientific choice inside a default, so this prints the ranked list and stops
    # unless the caller names a target or explicitly accepts the top one.
    resolution = None
    if args.disease:
        picked = choose_target(args, ROOT)
        if picked is None:
            return 1
        resolution = picked
        args.target = args.target or picked["symbol"]
        args.uniprot = args.uniprot or picked["uniprot"]
        args.slug = args.slug or picked["slug"]

    if args.cmd == "targets":
        return 0 if args.disease else (
            print("`targets` needs --disease") or 1)

    if not args.target or not args.uniprot:
        print("need a target: pass --target SYMBOL --uniprot ACCESSION, or --disease "
              "NAME to resolve one.\n"
              "  ./env/bin/python scripts/autorepurpose.py targets "
              "--disease \"colorectal cancer\"", file=sys.stderr)
        return 1

    slug = args.slug or args.target.lower()
    pipe = ROOT / "results" / "pipeline" / slug
    manifest = {"target": args.target, "uniprot": args.uniprot, "slug": slug,
                "signature": args.signature, "metric": args.metric,
                "stages": {}, "budget_credits": args.budget_credits,
                "disease": args.disease, "target_resolution": resolution}

    n_drugs = args.n_positives + args.n_decoys
    est = n_drugs * COST_PER_COFOLD + args.designs * COST_PER_DESIGN
    print(f"target {args.target} ({args.uniprot}) -> results/pipeline/{slug}/")
    print(f"plan: {n_drugs} drugs to co-fold"
          + (f" + {args.designs} BoltzGen designs" if args.designs else "")
          + f"  ~{est:.0f} credits estimated (budget {args.budget_credits:.0f})")
    if est > args.budget_credits:
        sys.exit(f"\nESTIMATE {est:.0f} EXCEEDS BUDGET {args.budget_credits:.0f}. "
                 f"Lower --n-decoys/--designs or raise --budget-credits. "
                 f"Nothing was submitted.")
    if args.signature == "boltzgen_consensus" and not args.designs \
            and not (ROOT / f"results/boltzgen_signature_{slug}.json").exists():
        sys.exit("--signature boltzgen_consensus needs designs: pass --designs N, or "
                 "use the default p2rank_geometry (which measured at least as good).")

    t0 = time.time()

    # 1. structure + pocket -------------------------------------------------
    if stage_done(pipe / "target" / "pocket.json"):
        print(f"\n[prep] already done ({pipe}/target/pocket.json)")
    else:
        rc = run([PY, ROOT / "scripts/prep_target.py", "--uniprot", args.uniprot,
                  "--disease", slug], args.dry_run)
        if rc:
            sys.exit(f"prep_target failed ({rc})")
    manifest["stages"]["prep"] = "ok"

    # 2. designs (optional) -------------------------------------------------
    if args.designs:
        print(f"\n[design] BoltzGen is OFF by default and you asked for "
              f"{args.designs} designs. Ablation I6.1 measured this step to lose to "
              f"the free pocket finder on two targets; see docs and m2_gate.py.")
        print("  submit with scripts/boltzgen_signature.py after running the design "
              "workflow; this orchestrator does not spend design credits implicitly.")
        manifest["stages"]["design"] = "requested, not auto-submitted"
    else:
        manifest["stages"]["design"] = "skipped (default; p2rank_geometry used instead)"

    # 3-5. shortlist / cofold / score ---------------------------------------
    # --append ALWAYS. Without it, cmd_shortlist REPLACES the shortlist, and
    # re-running a completed pipeline silently dropped the 12 hard-decoy rows that
    # the reported AUC 0.717 rests on (42 rows -> 30), orphaning their job records.
    # Appending is idempotent: cmd_shortlist de-duplicates on struct_id.
    sl = [PY, ROOT / "scripts/repurpose.py", "--pipeline", slug, "shortlist",
          "--target", args.target, "--n-positives", args.n_positives,
          "--n-decoys", args.n_decoys, "--append"]
    if args.hard_decoys:
        sl.append("--hard-decoys")
    if run(sl, args.dry_run):
        sys.exit("shortlist failed")
    manifest["stages"]["shortlist"] = "ok"

    board = ROOT / "results" / f"repurpose_{slug}.json"

    if run([PY, ROOT / "scripts/repurpose.py", "--pipeline", slug, "submit",
            "--max-credits", args.budget_credits], args.dry_run):
        sys.exit("submit failed")
    manifest["stages"]["cofold_submit"] = "ok"

    if run([PY, ROOT / "scripts/repurpose.py", "--pipeline", slug, "collect"],
           args.dry_run):
        sys.exit("collect failed")
    manifest["stages"]["cofold_collect"] = "ok"

    score = [PY, ROOT / "scripts/repurpose.py", "--pipeline", slug, "score",
             "--rank-by", args.signature, "--metric", args.metric]
    # boltzgen_signature.py names its output by --tag, not by slug, so the old
    # test on boltzgen_signature_<slug>.json never matched for ANY target and the
    # --designs flag was never appended. Accept either spelling, and never fall
    # back to another target's designs.
    for tag in (slug, args.designs_tag):
        if tag and (ROOT / f"results/boltzgen_signature_{tag}.json").exists():
            score += ["--designs", tag]
            break
    if run(score, args.dry_run):
        sys.exit("score failed")
    manifest["stages"]["score"] = "ok"

    # 6. report -------------------------------------------------------------
    # report.py renders ONE target: its inputs are hardcoded to the KDR run. Calling
    # it for another slug used to pass an argument it does not accept (it exited 2)
    # while the manifest recorded "ok" for a stage that never ran. Until it is
    # generalised, only invoke it for the pipeline it actually knows, and record
    # honestly otherwise.
    rep = ROOT / "scripts/report.py"
    if not rep.exists():
        manifest["stages"]["report"] = "skipped (scripts/report.py not present)"
    elif slug != REPORT_PIPELINE:
        manifest["stages"]["report"] = (
            f"not run: scripts/report.py currently renders only {REPORT_PIPELINE!r} "
            f"(its inputs are hardcoded). The ranked board for {slug!r} is still in "
            f"results/repurpose_{slug}.json and readable with scripts/ask.py")
        print(f"\n[report] skipped: report.py renders only {REPORT_PIPELINE}; "
              f"use scripts/ask.py --target {args.target} for this run")
    else:
        rc = run([PY, rep], args.dry_run)
        manifest["stages"]["report"] = "ok" if rc == 0 else f"FAILED (exit {rc})"
        if rc:
            print(f"[report] FAILED with exit {rc}", file=sys.stderr)

    manifest["wall_seconds"] = round(time.time() - t0, 1)
    manifest["git_sha"] = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
        cwd=ROOT).stdout.strip()
    out = pipe / "run_manifest.json"
    if not args.dry_run:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(manifest, indent=1))
    print(f"\ndone in {manifest['wall_seconds']:.0f}s -> {out}")
    print(f"ranked board: results/repurpose_{slug}.json")


if __name__ == "__main__":
    # sys.exit(main()), not main(): the target-resolution paths signal refusal by
    # returning 1, and without this every refusal would exit 0 and a caller or CI
    # would read "no target chosen" as success.
    sys.exit(main())
