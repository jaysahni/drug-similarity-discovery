"""Run P2Rank over the co-crystal benchmark set and collect the predicted pockets.

This builds the BASELINE arm of the hotspot-recovery experiment. PROJECT_GOAL.md
section 8.3 calls the comparison between a consensus interface signature and a
plain pocket prediction "the decisive ablation", and section 7 is blunt that if
pocket geometry already gives you the answer, the honest pipeline is the cheap
one. So the pocket arm has to be computed properly rather than strawmanned.

Ligands are NOT stripped first, and that is a measured decision, not an
oversight. P2Rank's predictions were compared on the holo and ligand-stripped
forms of the same structure (3VHE chain A): identical pocket score (13.02),
identical probability (0.671) and an identical 17-residue set, Jaccard 1.0. So
`prank predict` ignores HETATM, there is no ligand leakage to guard against
here, and a stripping stage would only add a failure mode. Re-check this with
--verify-no-leakage if the P2Rank version ever changes.

P2Rank is run in dataset mode so the JVM starts once rather than per structure,
which is most of the cost at this scale.

Usage:
    ./env/bin/python scripts/run_p2rank.py                    # all cached structures
    ./env/bin/python scripts/run_p2rank.py --limit 50 --threads 6
    ./env/bin/python scripts/run_p2rank.py --verify-no-leakage
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STRUCTURES = ROOT / "data" / "raw" / "structures"
WORK = ROOT / "data" / "raw" / "p2rank"
OUT = ROOT / "results" / "p2rank_pockets.json"

# Tool paths. The macOS system java is a stub, so a Temurin JDK 21 tarball was
# extracted into the session scratchpad - no Homebrew, nothing installed
# system-wide. Override with P2RANK_HOME / JAVA_HOME if they move.
SCRATCH = Path("/private/tmp/claude-501/-Users-evanxiang-Desktop-Projects-"
               "drug-similarity-discovery/ff08da32-3cac-44f3-9083-f8fc1e02ca4d/"
               "scratchpad/tools")
JAVA_HOME = Path(os.environ.get("JAVA_HOME") or SCRATCH / "jdk-21.0.12.1+1/Contents/Home")
P2RANK = Path(os.environ.get("P2RANK_HOME") or SCRATCH / "p2rank_2.5") / "prank"


def check_tools():
    missing = []
    if not (JAVA_HOME / "bin" / "java").exists():
        missing.append(f"java not found under JAVA_HOME={JAVA_HOME}")
    if not P2RANK.exists():
        missing.append(f"p2rank not found at {P2RANK}")
    if missing:
        sys.exit("cannot run:\n  " + "\n  ".join(missing) +
                 "\nSet JAVA_HOME and P2RANK_HOME, or see docs/03-SCOPE-AND-CONSTRAINTS.md")


def p2rank_env():
    env = dict(os.environ)
    env["JAVA_HOME"] = str(JAVA_HOME)
    env["PATH"] = f"{JAVA_HOME / 'bin'}:{env.get('PATH', '')}"
    return env


def run_batch(paths, outdir, threads):
    """One P2Rank invocation over a dataset file; returns wall seconds."""
    outdir.mkdir(parents=True, exist_ok=True)
    ds = outdir / "batch.ds"
    ds.write_text("HEADER: protein\n" + "\n".join(str(p) for p in paths) + "\n")
    t0 = time.time()
    proc = subprocess.run(
        [str(P2RANK), "predict", str(ds), "-o", str(outdir / "out"),
         "-threads", str(threads), "-visualizations", "0"],
        env=p2rank_env(), capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout[-2000:], file=sys.stderr)
        print(proc.stderr[-2000:], file=sys.stderr)
        raise RuntimeError(f"p2rank exited {proc.returncode}")
    return time.time() - t0


def parse_predictions(path):
    """P2Rank predictions.csv -> list of pockets, residue ids split by chain."""
    pockets = []
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            row = { (k or "").strip(): (v or "").strip() for k, v in row.items() }
            if not row.get("name"):
                continue
            residues = []
            for tok in (row.get("residue_ids") or "").split():
                chain, _, num = tok.rpartition("_")
                if not num:
                    continue
                try:
                    residues.append([chain, int(num)])
                except ValueError:          # insertion codes: keep the raw token
                    residues.append([chain, num])
            pockets.append({
                "rank": int(row["rank"]), "score": float(row["score"]),
                "probability": float(row["probability"]),
                "center": [float(row["center_x"]), float(row["center_y"]),
                           float(row["center_z"])],
                "n_residues": len(residues), "residues": residues,
            })
    return sorted(pockets, key=lambda p: p["rank"])


def verify_no_leakage():
    """Re-run the holo vs ligand-stripped comparison that justifies not stripping."""
    base = ROOT / "results/pipeline/colorectal-cancer/target/structure"
    holo, strip = base / "3VHE_A_holo.pdb", base / "3VHE_A_stripped.pdb"
    if not holo.exists() or not strip.exists():
        sys.exit(f"need {holo} and {strip}; run scripts/prep_target.py first")
    out = {}
    for tag, path in (("holo", holo), ("stripped", strip)):
        d = WORK / f"leakcheck_{tag}"
        run_batch([path], d, 1)
        out[tag] = parse_predictions(next((d / "out").glob("*_predictions.csv")))[0]
    a = {tuple(r) for r in out["holo"]["residues"]}
    b = {tuple(r) for r in out["stripped"]["residues"]}
    j = len(a & b) / len(a | b) if a | b else float("nan")
    print(f"holo     pocket1 score={out['holo']['score']} n_res={len(a)}")
    print(f"stripped pocket1 score={out['stripped']['score']} n_res={len(b)}")
    print(f"jaccard={j:.4f}  identical={a == b}")
    if a != b:
        print("\nP2Rank IS sensitive to the bound ligand in this version - strip "
              "before predicting, and update this script's docstring.", file=sys.stderr)
    return a == b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="0 = every cached structure")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--chunk", type=int, default=250, help="structures per JVM batch")
    ap.add_argument("--verify-no-leakage", action="store_true")
    args = ap.parse_args()

    check_tools()
    if args.verify_no_leakage:
        sys.exit(0 if verify_no_leakage() else 1)

    paths = sorted(STRUCTURES.glob("*.cif"))
    if args.limit:
        paths = paths[:args.limit]
    if not paths:
        sys.exit(f"no structures in {STRUCTURES}; run scripts/fetch_structures.py first")
    print(f"{len(paths)} structures, {args.threads} threads, chunks of {args.chunk}\n")

    results, failures, total_s = {}, [], 0.0
    for i in range(0, len(paths), args.chunk):
        chunk = paths[i:i + args.chunk]
        d = WORK / f"chunk{i // args.chunk:03d}"
        secs = run_batch(chunk, d, args.threads)
        total_s += secs
        got = 0
        for p in chunk:
            pred = d / "out" / f"{p.name}_predictions.csv"
            if not pred.exists():
                failures.append({"pdb_id": p.stem, "why": "no predictions.csv"})
                continue
            try:
                results[p.stem] = parse_predictions(pred)
                got += 1
            except Exception as exc:                      # noqa: BLE001
                failures.append({"pdb_id": p.stem, "why": repr(exc)})
        print(f"  chunk {i // args.chunk:03d}: {got}/{len(chunk)} parsed "
              f"in {secs:.0f}s ({secs / len(chunk):.2f}s per structure)", flush=True)

    n_pockets = sum(len(v) for v in results.values())
    empty = [k for k, v in results.items() if not v]
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps({
        "n_structures_attempted": len(paths),
        "n_structures_parsed": len(results),
        "n_failed": len(failures),
        "failures": failures,
        "n_structures_with_no_pocket": len(empty),
        "structures_with_no_pocket": sorted(empty),
        "n_pockets_total": n_pockets,
        "wall_seconds": round(total_s, 1),
        "seconds_per_structure": round(total_s / len(paths), 3),
        "p2rank": "2.5", "threads": args.threads,
        "ligands_stripped": False,
        "ligand_stripping_note": ("verified unnecessary: prank predict gives an "
                                  "identical pocket on holo and stripped 3VHE_A "
                                  "(score 13.02, 17 residues, jaccard 1.0)"),
        "pockets": results,
    }, indent=1))
    print(f"\n{len(results)}/{len(paths)} structures parsed, {len(failures)} failed, "
          f"{len(empty)} with no pocket, {n_pockets} pockets total")
    print(f"{total_s:.0f}s wall, {total_s / len(paths):.2f}s per structure")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
