"""Local driver for the BoltzGen arm: plan it, price it, run it, record it.

    ./env/bin/python scripts/run_design.py --dry-run            # plan + cost, submits nothing
    ./env/bin/python scripts/run_design.py smoke                # ~cents: does the image build?
    ./env/bin/python scripts/run_design.py smoke --gpu A100-40GB
    ./env/bin/python scripts/run_design.py weights              # one-off 6.35 GB download
    ./env/bin/python scripts/run_design.py check                # boltzgen check, CPU, no GPU spend
    ./env/bin/python scripts/run_design.py design --num-designs 50

Everything that costs money goes through the cost guard in `estimate()` and is
refused above --budget-usd without --force (PROJECT_GOAL.md rule 10). `--dry-run`
works with no credentials at all and submits nothing.

THE COST NUMBERS IN THIS FILE ARE ESTIMATES AND ARE LABELLED AS SUCH.

CLAUDE.md forbids reporting a number that was not computed. So, precisely: the
arithmetic here is computed, but two of its inputs are not ours and one is a
guess:

  * Per-stage seconds/design come from BoltzGen's own Figure 12
    (assets/fig_seconds_per_design.png), measured on a single A100. There is no
    published numeric table -- the values were read off the plot, so +/-10%.
  * `design_folding` (refold the binder alone, which `protein-anything` runs and
    Figure 12 does not plot) is EXTRAPOLATED from the refold curve evaluated at
    the binder length. Call it +/-20 s/design of extra uncertainty.
  * Non-A100 GPUs use a slowdown factor that NOBODY has measured for BoltzGen.
    Those factors are guesses from bandwidth/throughput ratios and are the
    softest numbers here. A GPU with no factor is refused rather than guessed at.

Nothing in `results/` may quote these as measured. After a real run the manifest
carries the actual wall time and the actual cost, and those are the numbers that
belong in the README.

RESIDUE NUMBERING. `--pocket` is read in author numbering (auth_seq_id), which
is what pocket.json, P2Rank and scripts/interfaces.py all use. BoltzGen wants
1-based indices over observed residues. The conversion happens in
scripts/modal_boltzgen.py, is refused if any residue is unmapped, and is echoed
here by --dry-run so it can be checked by eye before anything is submitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
DEFAULT_DISEASE = "colorectal-cancer"
ENV_FILE = ROOT / ".env"

EXIT_OK, EXIT_USAGE, EXIT_NO_CREDS, EXIT_BUDGET, EXIT_REMOTE = 0, 1, 2, 3, 4


# ==========================================================================
# credentials -- never crash, never pretend
# ==========================================================================
def load_dotenv(path: Path = ENV_FILE) -> dict:
    """Minimal KEY=VALUE reader. No dependency added for four lines of parsing.

    Does not overwrite variables already set in the environment, so an explicit
    `export` always wins over the file.
    """
    found = {}
    if not path.exists():
        return found
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        if key.startswith("export "):
            key = key[len("export "):].strip()
        val = val.strip().strip('"').strip("'")
        if not key:
            continue
        found[key] = val
        os.environ.setdefault(key, val)
    return found


def _modal_toml_profiles() -> list:
    """Profile names in ~/.modal.toml that carry both a token id and secret.

    Parsed by hand (tomllib exists in 3.11+, but this stays dependency- and
    version-free and never reads the secret values into anything we print).
    """
    cfg = Path(os.environ.get("MODAL_CONFIG_PATH") or Path.home() / ".modal.toml")
    if not cfg.exists():
        return []
    profiles, current, keys = [], None, set()
    for raw in cfg.read_text().splitlines():
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            if current and {"token_id", "token_secret"} <= keys:
                profiles.append(current)
            current, keys = line[1:-1], set()
        elif "=" in line:
            keys.add(line.split("=")[0].strip())
    if current and {"token_id", "token_secret"} <= keys:
        profiles.append(current)
    return profiles


def resolve_credentials() -> dict:
    """Where Modal credentials would come from, without validating them remotely."""
    dotenv = load_dotenv()
    if os.environ.get("MODAL_TOKEN_ID") and os.environ.get("MODAL_TOKEN_SECRET"):
        src = ".env" if "MODAL_TOKEN_ID" in dotenv else "environment"
        return {"ok": True, "source": src,
                "detail": f"MODAL_TOKEN_ID / MODAL_TOKEN_SECRET from {src}"}
    profiles = _modal_toml_profiles()
    if profiles:
        cfg = os.environ.get("MODAL_CONFIG_PATH") or "~/.modal.toml"
        return {"ok": True, "source": cfg,
                "detail": f"{cfg} profile(s): {', '.join(profiles)}"}
    return {"ok": False, "source": None, "detail": "no Modal credentials found"}


NO_CREDS_MESSAGE = f"""\
No Modal credentials found, so nothing was submitted and nothing was charged.

Put both lines in {ENV_FILE} (already gitignored -- .env is in .gitignore):

    MODAL_TOKEN_ID=ak-xxxxxxxxxxxxxxxxxxxxxx
    MODAL_TOKEN_SECRET=as-xxxxxxxxxxxxxxxxxxxxxx

Get them from https://modal.com/settings/tokens (New token), or run
`./env/bin/modal token set --token-id ak-... --token-secret as-...`, which
writes ~/.modal.toml instead; this driver accepts either source.

`--dry-run` needs no credentials and will print the plan and the cost estimate.
"""


# ==========================================================================
# cost model -- ESTIMATES, with every input's provenance attached
# ==========================================================================
# BoltzGen Figure 12, one A100, seconds per design, by TOTAL residues
# (design + target). Read off the PNG by eye: +/-10%. There is no published
# numeric table.
FIG12_A100 = {
    100: {"design": 4.9, "inverse_folding": 0.8, "folding": 21.1, "analysis": 1.0},
    200: {"design": 6.2, "inverse_folding": 1.3, "folding": 36.0, "analysis": 1.9},
    300: {"design": 6.5, "inverse_folding": 1.5, "folding": 46.0, "analysis": 1.8},
    400: {"design": 9.6, "inverse_folding": 2.0, "folding": 60.5, "analysis": 1.8},
    500: {"design": 12.9, "inverse_folding": 4.1, "folding": 81.3, "analysis": 1.8},
}
# Filtering is a one-off, not per design: README "~15 sec", preprint "around 20
# seconds, independent of the number of designs".
FILTERING_SECONDS = 20.0
# Steps that consume GPU time. `analysis` is CPU-side; `affinity` only runs for
# protein-small_molecule, which we do not use (PROJECT_GOAL.md 6.C2).
GPU_STEPS = {"design", "inverse_folding", "folding", "design_folding"}
PROTOCOL_STEPS = {
    "protein-anything": ["design", "inverse_folding", "folding", "design_folding",
                         "analysis", "filtering"],
    "peptide-anything": ["design", "inverse_folding", "folding", "analysis",
                         "filtering"],
    "nanobody-anything": ["design", "inverse_folding", "folding", "analysis",
                          "filtering"],
    "antibody-anything": ["design", "inverse_folding", "folding", "analysis",
                          "filtering"],
    "protein-redesign": ["design", "inverse_folding", "folding", "analysis",
                         "filtering"],
}
# Modal list prices, $/second, modal.com/pricing as fetched 2026-09-19.
GPU_USD_PER_SEC = {
    "T4": 0.000164, "L4": 0.000222, "A10": 0.000306, "L40S": 0.000542,
    "A100": 0.000583, "A100-40GB": 0.000583, "A100-80GB": 0.000694,
    "RTX-PRO-6000": 0.000842, "H100": 0.001097, "H200": 0.001261,
    "B200": 0.001736, "B300": 0.001972,
}
CPU_USD_PER_CORE_SEC = 0.0000131
MEM_USD_PER_GIB_SEC = 0.00000222
# Slowdown vs A100. GUESSES from bandwidth/throughput ratios -- no one has
# published BoltzGen timings on any GPU but A100. None = refuse to estimate.
GPU_SLOWDOWN_VS_A100 = {
    "A100": 1.0, "A100-40GB": 1.0, "A100-80GB": 1.0, "H100": 0.6,
    "L40S": 1.6, "A10": 2.75, "L4": 3.5,
    "T4": None, "H200": None, "B200": None, "B300": None, "RTX-PRO-6000": None,
}
# One-off container costs, ESTIMATED wall-clock (neither has been measured).
IMAGE_BUILD_SECONDS = 900.0      # torch-sized pip install, first build only
WEIGHTS_DOWNLOAD_SECONDS = 1200.0  # 6.35 GB onto the Volume, first time only
WEIGHTS_CPU, WEIGHTS_MEM_GIB = 4.0, 16.0
DESIGN_CPU, DESIGN_MEM_GIB = 8.0, 64.0   # 64 GiB is forced: boltzgen issue #208
                                         # OOM'd in `analysis` at 16 GiB host RAM


def _interp(table_key: str, residues: float) -> float:
    xs = sorted(FIG12_A100)
    if residues <= xs[0]:
        return FIG12_A100[xs[0]][table_key]
    if residues >= xs[-1]:
        # Clamped, not extrapolated: past 500 residues the curve is unknown and
        # a linear guess would understate it.
        return FIG12_A100[xs[-1]][table_key]
    for lo, hi in zip(xs, xs[1:]):
        if lo <= residues <= hi:
            f = (residues - lo) / (hi - lo)
            return FIG12_A100[lo][table_key] + f * (
                FIG12_A100[hi][table_key] - FIG12_A100[lo][table_key])
    raise AssertionError                         # pragma: no cover


def seconds_per_design(target_residues: int, binder_min: int, binder_max: int,
                       protocol: str, steps=None) -> dict:
    """A100-seconds per design, by stage. Inputs are Figure 12 + one extrapolation."""
    binder_mean = (binder_min + binder_max) / 2.0
    total = target_residues + binder_mean
    plan = steps or PROTOCOL_STEPS.get(protocol, PROTOCOL_STEPS["protein-anything"])
    per = {}
    for step in plan:
        if step == "design":
            per[step] = _interp("design", total)
        elif step == "inverse_folding":
            per[step] = _interp("inverse_folding", total)
        elif step == "folding":
            per[step] = _interp("folding", total)
        elif step == "design_folding":
            # EXTRAPOLATION: Figure 12 does not plot this step. Approximated as
            # a refold of the binder alone, i.e. the refold curve at binder-only
            # length. Worth ~+/-20 s/design of uncertainty on its own.
            per[step] = _interp("folding", binder_mean)
        elif step == "analysis":
            per[step] = _interp("analysis", total)
        elif step == "filtering":
            per[step] = 0.0                      # per-run, added separately
        else:
            per[step] = 0.0
    return {
        "per_stage": {k: round(v, 2) for k, v in per.items()},
        "gpu_seconds": round(sum(v for k, v in per.items() if k in GPU_STEPS), 2),
        "cpu_seconds": round(sum(v for k, v in per.items() if k not in GPU_STEPS), 2),
        "total_residues": round(total, 1),
        "binder_mean": binder_mean,
    }


def estimate(num_designs: int, gpu: str, target_residues: int, binder_min: int,
             binder_max: int, protocol: str, steps=None,
             include_overheads: bool = True) -> dict:
    """Full cost estimate. Raises ValueError for a GPU with no slowdown factor."""
    if gpu not in GPU_USD_PER_SEC:
        raise ValueError(
            f"unknown Modal GPU {gpu!r}. Known: {', '.join(sorted(GPU_USD_PER_SEC))}. "
            "Note the current Modal string is 'A10', not 'A10G'.")
    factor = GPU_SLOWDOWN_VS_A100.get(gpu)
    if factor is None:
        raise ValueError(
            f"no BoltzGen slowdown factor is known for {gpu!r}, so its cost "
            "cannot be estimated without inventing a number. Use A100-40GB "
            "(cheapest per design of the tiers with a factor) or H100.")
    per = seconds_per_design(target_residues, binder_min, binder_max, protocol, steps)
    container_seconds = (per["gpu_seconds"] * factor + per["cpu_seconds"]) * num_designs
    container_seconds += FILTERING_SECONDS
    rate = (GPU_USD_PER_SEC[gpu]
            + DESIGN_CPU * CPU_USD_PER_CORE_SEC
            + DESIGN_MEM_GIB * MEM_USD_PER_GIB_SEC)
    gpu_usd = container_seconds * rate
    overhead_rate = (WEIGHTS_CPU * CPU_USD_PER_CORE_SEC
                     + WEIGHTS_MEM_GIB * MEM_USD_PER_GIB_SEC)
    overhead_usd = ((IMAGE_BUILD_SECONDS + WEIGHTS_DOWNLOAD_SECONDS) * overhead_rate
                    if include_overheads else 0.0)
    return {
        "gpu": gpu,
        "slowdown_vs_a100": factor,
        "num_designs": num_designs,
        "protocol": protocol,
        "per_design": per,
        "container_seconds": round(container_seconds, 1),
        "container_hours": round(container_seconds / 3600, 2),
        "usd_per_second": round(rate, 8),
        "gpu_usd": round(gpu_usd, 2),
        "one_off_overhead_usd": round(overhead_usd, 2),
        "total_usd": round(gpu_usd + overhead_usd, 2),
        "usd_per_design": round(gpu_usd / num_designs, 4) if num_designs else 0.0,
        "basis": "ESTIMATE: BoltzGen Figure 12 (A100, read off a PNG, +/-10%) + "
                 "an extrapolated design_folding stage + a GUESSED per-GPU "
                 "slowdown factor + Modal list prices of 2026-09-19. Not measured.",
    }


# ==========================================================================
# target inputs
# ==========================================================================
def target_paths(disease: str, structure: str | None, pocket: str | None) -> dict:
    base = ROOT / "results" / "pipeline" / disease / "target"
    s = Path(structure) if structure else base / "structure" / "3VHE_A_stripped.pdb"
    p = Path(pocket) if pocket else base / "pocket.json"
    return {"structure": s, "pocket": p, "base": base}


def load_target(disease, structure, pocket) -> dict:
    paths = target_paths(disease, structure, pocket)
    missing = [str(v) for k, v in paths.items()
               if k in ("structure", "pocket") and not v.exists()]
    if missing:
        raise FileNotFoundError(
            "target inputs are missing: " + ", ".join(missing)
            + "\nRun scripts/prep_target.py first, or pass --structure/--pocket.")
    pdb_bytes = paths["structure"].read_bytes()
    site = json.loads(paths["pocket"].read_text())
    residue_ids = site.get("residue_ids")
    if not residue_ids:
        raise ValueError(f"{paths['pocket']} has no non-empty 'residue_ids'")
    return {
        "structure_path": paths["structure"],
        "pocket_path": paths["pocket"],
        "pdb_bytes": pdb_bytes,
        "sha256": hashlib.sha256(pdb_bytes).hexdigest(),
        "residue_ids": [int(r) for r in residue_ids],
        "site_id": site.get("site_id", paths["structure"].stem),
        "chain_id": site.get("chain_id"),
        "uniprot_id": site.get("uniprot_id"),
        "name": paths["structure"].stem,
    }


def import_app():
    """Import the Modal app module. Returns (module, error-string)."""
    sys.path.insert(0, str(SCRIPTS))
    try:
        import modal_boltzgen

        return modal_boltzgen, None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def spec_hash(target, args, mod) -> str:
    """Content address for this exact design batch (parallel-execution rule 6)."""
    payload = json.dumps({
        "target_sha256": target["sha256"],
        "residues": sorted(target["residue_ids"]),
        "chain": target["chain_id"],
        "protocol": args.protocol,
        "num_designs": args.num_designs,
        "seed": args.seed,
        "binder": [args.binder_min, args.binder_max],
        "steps": args.steps,
        "boltzgen": getattr(mod, "BOLTZGEN_VERSION", "unknown"),
        "python": getattr(mod, "PYTHON_VERSION", "unknown"),
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def run_dir_for(disease, target, args, digest) -> Path:
    return (ROOT / "results" / "pipeline" / disease / "designs"
            / f"{target['site_id']}__{args.protocol}__n{args.num_designs}"
              f"__seed{args.seed}__{digest}")


# ==========================================================================
# printing
# ==========================================================================
def hr(title):
    print(f"\n=== {title} " + "=" * max(0, 64 - len(title)))


def print_plan(target, args, est, mapping, creds, digest, out_dir, mod):
    hr("plan")
    print(f"disease           {args.disease}")
    print(f"target            {target['site_id']}  ({target['uniprot_id'] or 'no uniprot id'})")
    print(f"structure         {target['structure_path'].relative_to(ROOT)}")
    print(f"                  sha256 {target['sha256'][:16]}...")
    print(f"pocket            {target['pocket_path'].relative_to(ROOT)}  "
          f"({len(target['residue_ids'])} residues, author numbering)")
    print(f"protocol          {args.protocol}")
    print(f"num_designs       {args.num_designs}   (batches of {args.designs_per_batch})")
    print(f"binder length     {args.binder_min}..{args.binder_max}")
    print(f"steps             {args.steps or 'full protocol pipeline'}")
    print(f"seed              {args.seed}  (applied only if this boltzgen build "
          f"advertises --seed; the manifest records which)")
    print(f"gpu               {args.gpu}")
    print(f"boltzgen          {getattr(mod, 'BOLTZGEN_VERSION', '?')} on python "
          f"{getattr(mod, 'PYTHON_VERSION', '?')} (Modal image)")
    print(f"spec hash         {digest}")
    print(f"output            {out_dir.relative_to(ROOT)}")
    print(f"credentials       {creds['detail']}")

    if mapping:
        hr("residue mapping (author -> BoltzGen 1-based index)")
        print(f"{mapping['n_residues']} observed residues, "
              f"author {mapping['author_first']}-{mapping['author_last']}, "
              f"gaps {mapping['gaps']}")
        pairs = [(r, mapping['map'].get(r)) for r in sorted(target['residue_ids'])]
        print("  " + ", ".join(f"{a}->{b}" for a, b in pairs))
        unmapped = [a for a, b in pairs if b is None]
        if unmapped:
            print(f"  !! UNMAPPED: {unmapped} -- design would be refused")
        for w in mapping["warnings"]:
            print(f"  !! {w}")
        print("  This mapping is DERIVED, not confirmed by boltzgen. Run the "
              "`check` subcommand before spending GPU time.")

    if est:
        hr("cost estimate (NOT measured)")
        per = est["per_design"]
        print(f"system size       {per['total_residues']} residues "
              f"({target['pdb_residues']} target + {per['binder_mean']} mean binder)")
        print(f"per design        " + ", ".join(
            f"{k} {v}s" for k, v in per["per_stage"].items()))
        print(f"                  {per['gpu_seconds']}s GPU + {per['cpu_seconds']}s CPU "
              f"per design on A100")
        print(f"slowdown factor   x{est['slowdown_vs_a100']} for {est['gpu']} (GUESS)")
        print(f"container time    {est['container_seconds']}s = {est['container_hours']}h")
        print(f"rate              ${est['usd_per_second']}/s  "
              f"(GPU + {DESIGN_CPU} cores + {DESIGN_MEM_GIB} GiB)")
        print(f"compute           ${est['gpu_usd']}  (${est['usd_per_design']}/design)")
        print(f"one-off overhead  ${est['one_off_overhead_usd']}  "
              f"(image build + 6.35 GB weight download, first run only)")
        print(f"TOTAL ESTIMATE    ${est['total_usd']}   vs budget ${args.budget_usd}")
        print(f"basis             {est['basis']}")


def print_ladder(args, target):
    hr("cost ladder (same basis; --reuse makes the ladder cost the top row once)")
    print(f"{'N':>6}  {'GPU-h':>7}  {'$compute':>9}  {'$/design':>9}")
    for n in (50, 100, 250, 500, 1000):
        e = estimate(n, args.gpu, target["pdb_residues"], args.binder_min,
                     args.binder_max, args.protocol, args.steps,
                     include_overheads=False)
        print(f"{n:>6}  {e['container_hours']:>7}  {e['gpu_usd']:>9}  "
              f"{e['usd_per_design']:>9}")


# ==========================================================================
# subcommands
# ==========================================================================
def require_creds(creds):
    if not creds["ok"]:
        print(NO_CREDS_MESSAGE)
        return False
    return True


def cmd_smoke(args, mod, creds):
    import modal

    with modal.enable_output():
        with mod.app.run():
            fn = mod.smoke_gpu if args.gpu_probe else mod.smoke
            if args.gpu_probe and args.gpu != mod.DEFAULT_GPU:
                fn = fn.with_options(gpu=args.gpu)
            result = fn.remote()
    out = ROOT / "results" / "pipeline" / args.disease / "designs" / "smoke.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    prior = json.loads(out.read_text()) if out.exists() else []
    prior.append({"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  **result})
    out.write_text(json.dumps(prior, indent=2))
    print(json.dumps(result, indent=2)[:6000])
    print(f"\nwrote {out.relative_to(ROOT)}")
    return EXIT_OK


def cmd_weights(args, mod, creds):
    import modal

    with modal.enable_output():
        with mod.app.run():
            result = mod.fetch_weights.remote(force=args.force)
    print(json.dumps({k: v for k, v in result.items()
                      if k != "weights_volume"}, indent=2)[:6000])
    vol = result.get("weights_volume", {})
    print(f"weights volume now holds {vol.get('gib')} GiB in "
          f"{vol.get('n_files')} files")
    return EXIT_OK if result.get("status") != "failed" else EXIT_REMOTE


def cmd_check(args, mod, creds, target, digest, out_dir):
    import modal

    out_dir.mkdir(parents=True, exist_ok=True)
    with modal.enable_output():
        with mod.app.run():
            result = mod.check_spec.remote(
                target_pdb=target["pdb_bytes"],
                target_name=target["name"],
                pocket_residue_ids=target["residue_ids"],
                chain_id=target["chain_id"],
                binder_min=args.binder_min,
                binder_max=args.binder_max,
                run_key=f"{args.disease}/{digest}",
            )
    (out_dir / "check.json").write_text(json.dumps(result, indent=2))
    print(result["yaml"])
    print(f"boltzgen check returncode={result['returncode']} "
          f"in {result['seconds']}s; emitted {len(result['emitted_files'])} files")
    if not result["ok"]:
        print(result["stderr_tail"][-2000:])
    print(f"\nwrote {(out_dir / 'check.json').relative_to(ROOT)}")
    print(f"structures written to the '{mod.RUNS_VOLUME}' Volume under "
          f"{result['runs_volume_path']}; fetch with:\n"
          f"  ./env/bin/modal volume get {mod.RUNS_VOLUME} "
          f"{result['runs_volume_path']} {out_dir.relative_to(ROOT)}/check_out")
    return EXIT_OK if result["ok"] else EXIT_REMOTE


def cmd_design(args, mod, creds, target, digest, out_dir, est):
    import modal

    out_dir.mkdir(parents=True, exist_ok=True)
    batches = []
    remaining, i = args.num_designs, 0
    while remaining > 0:
        n = min(args.designs_per_batch, remaining)
        batches.append((i, n))
        remaining -= n
        i += 1

    manifest_path = out_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "run_id": out_dir.name,
        "spec_hash": digest,
        "disease": args.disease,
        "site_id": target["site_id"],
        "uniprot_id": target["uniprot_id"],
        "structure": str(target["structure_path"].relative_to(ROOT)),
        "structure_sha256": target["sha256"],
        "pocket_residue_ids_author": target["residue_ids"],
        "protocol": args.protocol,
        "num_designs_requested": args.num_designs,
        "designs_per_batch": args.designs_per_batch,
        "seed": args.seed,
        "binder_length": [args.binder_min, args.binder_max],
        "steps": args.steps,
        "gpu": args.gpu,
        "boltzgen_version": mod.BOLTZGEN_VERSION,
        "image_python": mod.PYTHON_VERSION,
        "modal_app": mod.APP_NAME,
        "cost_estimate_before_run": est,
        "batches": {},
    }

    exit_code = EXIT_OK
    with modal.enable_output():
        with mod.app.run():
            fn = mod.design
            if args.gpu != mod.DEFAULT_GPU:
                fn = fn.with_options(gpu=args.gpu)
            for idx, n in batches:
                bdir = out_dir / f"batch-{idx:04d}"
                done = bdir / "result.json"
                if done.exists() and not args.force_rerun:
                    prev = json.loads(done.read_text())
                    if prev.get("ok"):
                        print(f"batch {idx}: already complete "
                              f"({prev.get('n_structures_returned')} structures), skipping")
                        continue
                bdir.mkdir(parents=True, exist_ok=True)
                t0 = time.time()
                print(f"batch {idx}: submitting {n} designs on {args.gpu} ...",
                      flush=True)
                result = fn.remote(
                    target_pdb=target["pdb_bytes"],
                    target_name=target["name"],
                    pocket_residue_ids=target["residue_ids"],
                    protocol=args.protocol,
                    num_designs=n,
                    seed=args.seed,
                    chain_id=target["chain_id"],
                    binder_min=args.binder_min,
                    binder_max=args.binder_max,
                    steps=args.steps,
                    diffusion_batch_size=args.diffusion_batch_size,
                    use_kernels=args.use_kernels,
                    reuse=True,
                    run_key=f"{args.disease}/{digest}/batch-{idx:04d}",
                )
                wall = time.time() - t0
                _write_batch(bdir, result)
                rate = (GPU_USD_PER_SEC[args.gpu]
                        + DESIGN_CPU * CPU_USD_PER_CORE_SEC
                        + DESIGN_MEM_GIB * MEM_USD_PER_GIB_SEC)
                manifest["batches"][f"batch-{idx:04d}"] = {
                    "n_designs": n,
                    "ok": result.get("ok"),
                    "local_wall_seconds": round(wall, 1),
                    "container_wall_seconds": result.get("wall_seconds"),
                    "billed_usd_estimate_from_container_wall": round(
                        (result.get("wall_seconds") or 0) * rate, 4),
                    "usd_rate_per_second": round(rate, 8),
                    "gpu_reported": result.get("gpu", {}).get("nvidia_smi"),
                    "boltzgen_version_reported": result.get("boltzgen", {}).get("version"),
                    "seed_applied": result.get("seed_applied"),
                    "seed_note": result.get("seed_note"),
                    "structure_source": result.get("structure_source"),
                    "n_structures": result.get("n_structures_returned"),
                    "structures_truncated": result.get("structures_truncated"),
                    "runs_volume_path": result.get("runs_volume_path"),
                    "returncode": result.get("returncode"),
                }
                manifest_path.write_text(json.dumps(manifest, indent=2))
                if not result.get("ok"):
                    print(f"batch {idx}: FAILED (returncode "
                          f"{result.get('returncode')}). Stopping so the next "
                          f"batch does not repeat the same failure on the clock.")
                    print((result.get("stderr_tail") or "")[-3000:])
                    exit_code = EXIT_REMOTE
                    break

    billed = sum(b.get("billed_usd_estimate_from_container_wall", 0.0)
                 for b in manifest["batches"].values())
    manifest["billed_usd_estimate_total"] = round(billed, 2)
    manifest["note"] = (
        "billed_* figures are container wall time x Modal list rate. They are "
        "computed from the actual wall time of this run, but they are not an "
        "invoice -- check the Modal dashboard for what was really charged.")
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest: {manifest_path.relative_to(ROOT)}")
    print(f"estimated spend from measured wall time: ${billed:.2f} "
          f"(pre-run estimate was ${est['total_usd']})")
    print(f"\nFull BoltzGen output trees stay on the '{mod.RUNS_VOLUME}' Volume. "
          f"Fetch one with:\n  ./env/bin/modal volume get {mod.RUNS_VOLUME} "
          f"{args.disease}/{digest} {out_dir.relative_to(ROOT)}/volume")
    return exit_code


def _write_batch(bdir: Path, result: dict):
    designs = bdir / "designs"
    designs.mkdir(parents=True, exist_ok=True)
    for name, text in (result.get("structures") or {}).items():
        (designs / name).write_text(text)
    for m in result.get("metrics") or []:
        if "rows" not in m:
            continue
        name = Path(m["path"]).name
        rows = m["rows"]
        if not rows:
            continue
        import csv

        with open(bdir / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    slim = {k: v for k, v in result.items() if k not in ("structures", "metrics")}
    slim["metrics_summary"] = [{k: v for k, v in m.items() if k != "rows"}
                               for m in (result.get("metrics") or [])]
    (bdir / "result.json").write_text(json.dumps(slim, indent=2))
    if result.get("design_spec_yaml"):
        (bdir / "design_spec.yaml").write_text(result["design_spec_yaml"])


# ==========================================================================
# cli
# ==========================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Estimates are estimates. Nothing here is a measured number "
               "until a run writes a manifest.")
    p.add_argument("command", nargs="?", default=None,
                   choices=["smoke", "weights", "check", "design"],
                   help="smoke: cheap validation. weights: one-off download. "
                        "check: boltzgen check on the YAML (no GPU). "
                        "design: the real run. Omitting it is only allowed with "
                        "--dry-run, so a bare invocation can never spend money.")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and the cost estimate; submit nothing. "
                        "Needs no credentials.")
    p.add_argument("--disease", default=DEFAULT_DISEASE)
    p.add_argument("--structure", default=None,
                   help="target .pdb (default: the prepared KDR/VEGFR2 3VHE_A_stripped)")
    p.add_argument("--pocket", default=None,
                   help="pocket JSON with author-numbered residue_ids "
                        "(default: the P2Rank pocket.json)")
    p.add_argument("--protocol", default="protein-anything")
    p.add_argument("--num-designs", type=int, default=50)
    p.add_argument("--designs-per-batch", type=int, default=100,
                   help="one Modal call per batch, so a 6 h function timeout "
                        "never truncates a large N")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--binder-min", type=int, default=80)
    p.add_argument("--binder-max", type=int, default=140)
    p.add_argument("--steps", nargs="+", default=None,
                   help="subset of the pipeline, e.g. --steps design "
                        "(~10x cheaper, but no confidence scores and no filtering)")
    p.add_argument("--diffusion-batch-size", type=int, default=None,
                   help="BoltzGen defaults to 1 below 100 designs and 10 at or "
                        "above; designs in a batch share a binder length, so a "
                        "large batch under-samples the length range")
    p.add_argument("--use-kernels", default="auto")
    p.add_argument("--gpu", default="A100-40GB",
                   help="Modal GPU string. Note it is 'A10', not 'A10G'.")
    p.add_argument("--gpu-probe", action="store_true",
                   help="smoke only: run the GPU probe instead of the CPU one")
    p.add_argument("--budget-usd", type=float, default=10.0,
                   help="cost guard (PROJECT_GOAL.md rule 10): refuse to submit "
                        "an estimate above this without --force")
    p.add_argument("--force", action="store_true",
                   help="submit anyway despite the cost guard (or re-download weights)")
    p.add_argument("--force-rerun", action="store_true",
                   help="re-run batches that already have output")
    p.add_argument("--ladder", action="store_true",
                   help="also print the cost ladder over N")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        # A bare `run_design.py` must never submit a GPU run by accident.
        if not args.dry_run:
            print("no subcommand given. Use one of: smoke, weights, check, "
                  "design -- or --dry-run to see the plan and the cost of a "
                  "design run without submitting anything.")
            return EXIT_USAGE
        args.command = "design"
    if args.protocol == "protein-small_molecule" and not args.force:
        print("protein-small_molecule designs a protein that BINDS a small "
              "molecule -- the reverse of what the interface signature needs "
              "(PROJECT_GOAL.md 6.C2). Pass --force if you really mean it.")
        return EXIT_USAGE
    if args.num_designs < 1 or args.designs_per_batch < 1:
        print("--num-designs and --designs-per-batch must be >= 1")
        return EXIT_USAGE

    creds = resolve_credentials()
    mod, import_err = import_app()
    if mod is None:
        print(f"could not import scripts/modal_boltzgen.py: {import_err}")
        print("Is the venv active?  ./env/bin/python scripts/run_design.py ...")
        return EXIT_USAGE

    try:
        target = load_target(args.disease, args.structure, args.pocket)
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc))
        return EXIT_USAGE

    mapping = None
    try:
        mapping = mod.residue_index_map(target["pdb_bytes"],
                                        chain_id=target["chain_id"])
        target["pdb_residues"] = mapping["n_residues"]
    except Exception as exc:
        print(f"could not read residues from {target['structure_path']}: {exc}")
        return EXIT_USAGE

    digest = spec_hash(target, args, mod)
    out_dir = run_dir_for(args.disease, target, args, digest)

    est = None
    if args.command == "design":
        try:
            est = estimate(args.num_designs, args.gpu, target["pdb_residues"],
                           args.binder_min, args.binder_max, args.protocol,
                           args.steps)
        except ValueError as exc:
            print(str(exc))
            return EXIT_USAGE

    if args.dry_run or args.command == "design":
        print_plan(target, args, est, mapping, creds, digest, out_dir, mod)
        if args.ladder:
            print_ladder(args, target)

    if args.dry_run:
        hr("dry run")
        print("Nothing was submitted and nothing was charged.")
        if args.command == "design" and est["total_usd"] > args.budget_usd:
            print(f"NOTE: the estimate ${est['total_usd']} exceeds the "
                  f"${args.budget_usd} budget; a real run would be refused "
                  f"without --force.")
        return EXIT_OK

    if not require_creds(creds):
        return EXIT_NO_CREDS

    # cost guard -- PROJECT_GOAL.md rule 10
    if args.command == "design" and est["total_usd"] > args.budget_usd:
        if not args.force:
            hr("refused by the cost guard")
            print(f"estimated ${est['total_usd']} > --budget-usd "
                  f"{args.budget_usd}. Nothing was submitted.")
            print("Options: lower --num-designs, use --steps design (roughly "
                  "10x cheaper, no confidence scores), raise --budget-usd, or "
                  "pass --force if you have decided to spend it.")
            return EXIT_BUDGET
        print(f"\ncost guard OVERRIDDEN with --force: submitting an estimated "
              f"${est['total_usd']} of work.\n")

    try:
        if args.command == "smoke":
            return cmd_smoke(args, mod, creds)
        if args.command == "weights":
            return cmd_weights(args, mod, creds)
        if args.command == "check":
            return cmd_check(args, mod, creds, target, digest, out_dir)
        return cmd_design(args, mod, creds, target, digest, out_dir, est)
    except KeyboardInterrupt:
        print("\ninterrupted; any completed batch is already on disk and will "
              "be skipped on the next run")
        return EXIT_REMOTE


if __name__ == "__main__":
    sys.exit(main())
