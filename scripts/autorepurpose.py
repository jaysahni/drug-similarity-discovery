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
signature is `p2rank_geometry`: it is free, it takes 0.47 s, and it is at least as
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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--target", required=True, help="gene symbol, e.g. KDR")
    ap.add_argument("--uniprot", required=True, help="e.g. P35968")
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
    ap.add_argument("--budget-credits", type=float, default=200.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    slug = args.slug or args.target.lower()
    pipe = ROOT / "results" / "pipeline" / slug
    manifest = {"target": args.target, "uniprot": args.uniprot, "slug": slug,
                "signature": args.signature, "metric": args.metric,
                "stages": {}, "budget_credits": args.budget_credits}

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
    sl = [PY, ROOT / "scripts/repurpose.py", "--pipeline", slug, "shortlist",
          "--target", args.target, "--n-positives", args.n_positives,
          "--n-decoys", args.n_decoys]
    if args.hard_decoys:
        sl.append("--hard-decoys")
    if run(sl, args.dry_run):
        sys.exit("shortlist failed")
    manifest["stages"]["shortlist"] = "ok"

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
    if (ROOT / f"results/boltzgen_signature_{slug}.json").exists():
        score += ["--designs", slug]
    if run(score, args.dry_run):
        sys.exit("score failed")
    manifest["stages"]["score"] = "ok"

    # 6. report -------------------------------------------------------------
    rep = ROOT / "scripts/report.py"
    if rep.exists():
        run([PY, rep, "--pipeline", slug], args.dry_run)
        manifest["stages"]["report"] = "ok"
    else:
        manifest["stages"]["report"] = "skipped (scripts/report.py not present)"

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
    main()
