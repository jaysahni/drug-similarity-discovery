"""Modal app for the BoltzGen arm of the M2 ablation (PROJECT_GOAL.md 6.C, 8.3).

The question this exists to answer is I6.1: does a consensus interface signature
built from BoltzGen designs carry information beyond the pocket P2Rank finds in
7.1 s? The P2Rank arm is built and measured. This is the other arm.

WHY MODAL AND NOT THIS MACHINE (hard constraints, verified, stated before we
build on them per CLAUDE.md):

  1. No NVIDIA GPU on this Mac. BoltzGen's design/refold steps are CUDA-only.
  2. `env/` is Python 3.14.7. boltzgen pins `numpy==2.0.2`, which publishes
     wheels only for cp39-cp312, so boltzgen cannot be installed into `env/`
     at all -- its declared `requires-python = ">=3.11"` has no upper bound but
     the effective ceiling is 3.12.
  3. Stronger still: boltzgen requires cuequivariance_ops_cu12,
     cuequivariance_ops_torch_cu12 and cuequivariance_torch, none of which
     publish ANY macOS wheel. No local Python version fixes that.

On Modal we choose the image, so (2) and (3) disappear and (1) is rented. This
module is the only place that knows how to build that image and run the tool.

WHAT IT EXPOSES

    smoke()          CPU. Python + boltzgen version, weights-volume state, disk.
                     A few cents. Run this FIRST -- it is the cheapest possible
                     proof that the image builds and boltzgen imports.
    smoke_gpu()      GPU. nvidia-smi + torch CUDA capability (>= 8.0 is what
                     `--use_kernels auto` needs). Separate function because
                     nvidia-smi does not exist on a CPU container.
    fetch_weights()  CPU. ~6.35 GB of checkpoints into a persistent Volume,
                     once, so GPU seconds are never spent downloading.
    check_spec()     CPU. `boltzgen check` on the generated YAML. This is the
                     guard on the residue-index trap described below.
    design()         GPU. The real run. Returns per-design confidence scores and
                     (optionally) the designed complex structures inline; always
                     writes everything to the runs Volume.

THE RESIDUE-INDEX TRAP -- the single highest-risk item in this file.

BoltzGen's design YAML takes 1-based *canonical mmCIF* residue indices, NOT PDB
author numbering. For a .pdb with no SEQRES, boltzgen's parse_pdb builds the
sequence from the residues physically present and assigns label_seq = j + 1,
i.e. 1..N over observed residues. Our target
(results/pipeline/colorectal-cancer/target/structure/3VHE_A_stripped.pdb) has
0 SEQRES lines, 303 observed residues, author numbering 811-1169 with two gaps.
Every one of the 17 P2Rank pocket residue ids (840..1047) is therefore OUT OF
RANGE as a BoltzGen index -- pasting them in would address residues that do not
exist, or silently the wrong ones.

So `design()` takes pocket ids in AUTHOR numbering (the numbering pocket.json
and scripts/interfaces.py use) and computes the mapping itself, from the same
bytes the container parses, via `residue_index_map()`. It refuses to run if any
requested residue is unmapped. The mapping is returned in the result so it can
be audited.

UNVERIFIED (nobody has run this; do not present any of it as fact):
  * That boltzgen's label_seq assignment matches `residue_index_map()` for this
    file. It was derived by reading pdb_parser.py, not by running boltzgen.
    `check_spec()` exists to settle it before GPU money is spent.
  * That `boltzgen check` runs without a GPU (it loads no checkpoint, which
    strongly implies CPU-only, but that is inference).
  * The exact `boltzgen download` argument shape, and whether `--cache` is
    honoured by every step. `_run_download()` tries more than one form and
    reports what happened rather than assuming.
  * Whether `boltzgen run` accepts `--seed`. `_supported_flags()` reads
    `--help` at runtime and the manifest records whether the seed was actually
    applied, so a run is never described as seeded when it was not.
  * Column names in the metrics CSVs. They are read as-is; expected confidence
    keys are reported as found/missing rather than assumed.
  * Whether MSA generation is triggered for the target and whether that needs
    network access inside the container.
  * VRAM needed for a ~413-residue system. Nobody has published one; the
    BoltzGen maintainers say so explicitly.

Nothing here executes at import time, so this module imports fine with no Modal
token present.
"""

from __future__ import annotations

import modal

# --------------------------------------------------------------------------
# pins -- everything that could change a number is pinned, per CLAUDE.md
# --------------------------------------------------------------------------
APP_NAME = "drug-similarity-boltzgen"
BOLTZGEN_VERSION = "0.3.2"
PYTHON_VERSION = "3.12"          # NOT optional: debian_slim() with no
                                 # python_version copies the LOCAL interpreter
                                 # (3.14 here), where boltzgen cannot install.
WEIGHTS_VOLUME = "boltzgen-weights"
RUNS_VOLUME = "boltzgen-runs"
WEIGHTS_MOUNT = "/weights"
RUNS_MOUNT = "/runs"
BOLTZGEN_CACHE = f"{WEIGHTS_MOUNT}/boltzgen"
HF_HOME = f"{WEIGHTS_MOUNT}/hf"
WEIGHTS_MARKER = f"{BOLTZGEN_CACHE}/.download-complete.json"

# Modal GPU strings, current reference (modal.com/docs/guide/gpu): T4, L4, A10,
# L40S, A100, A100-40GB, A100-80GB, RTX-PRO-6000, H100, H200, B200, B300.
# "A10G" is NOT in that list and the client does no validation -- a wrong string
# is only rejected server-side, after the image has built.
DEFAULT_GPU = "A100-40GB"

# Host RAM, not VRAM, is the documented OOM: boltzgen issue #208 OOM'd in the
# *analysis* step on an A100 with 16 GB host RAM at only 100 designs, and
# succeeded at 64 GB. Upstream's own slurm example asks for --mem=64G.
DESIGN_CPU = 8.0
DESIGN_MEMORY_MB = 65536

PROTOCOLS = (
    "protein-anything",
    "peptide-anything",
    "protein-small_molecule",   # designs a protein that binds a SMALL MOLECULE
    "nanobody-anything",        # -- the reverse of what we want (6.C2). Listed
    "antibody-anything",        # for completeness; run_design.py blocks it.
    "protein-redesign",
)
PIPELINE_STEPS = (
    "design", "inverse_folding", "design_folding", "folding",
    "affinity", "analysis", "filtering",
)

# const.eval_keys_confidence in boltzgen 0.3.2. Read from source, not from a
# real CSV -- treated as EXPECTED, and reported as found/missing.
EXPECTED_CONFIDENCE_KEYS = (
    "design_iptm", "design_ptm", "design_iiptm", "design_to_target_iptm",
    "design_residue_iptm", "target_ptm", "iptm", "ptm", "protein_iptm",
    "ligand_iptm", "interaction_pae", "min_interaction_pae",
    "min_design_to_target_pae", "design_ipsae_min", "design_to_target_ipsae",
    "target_to_design_ipsae", "complex_plddt", "complex_iplddt",
    "complex_pde", "complex_ipde",
)

IMAGE = (
    modal.Image.debian_slim(python_version=PYTHON_VERSION)
    .apt_install("git")
    .pip_install(f"boltzgen=={BOLTZGEN_VERSION}")
    .env(
        {
            "HF_HOME": HF_HOME,
            "PYTHONUNBUFFERED": "1",
            # boltzgen needs no HF token (boltzgen/boltzgen-1 is public,
            # gated=False, MIT). The README's hardcoded hf_eOOQ... default is
            # stale for 0.3.2; the code reads os.environ.get("HF_TOKEN").
        }
    )
)

weights_volume = modal.Volume.from_name(WEIGHTS_VOLUME, create_if_missing=True)
runs_volume = modal.Volume.from_name(RUNS_VOLUME, create_if_missing=True)

app = modal.App(APP_NAME, image=IMAGE)


# ==========================================================================
# pure helpers -- no modal, no boltzgen, importable and testable anywhere
# ==========================================================================
def residue_index_map(pdb_bytes: bytes, chain_id: str | None = None) -> dict:
    """Author residue number -> 1-based BoltzGen residue index, for one chain.

    Mirrors boltzgen's parse_pdb for a SEQRES-less .pdb: the sequence is the
    residues physically present, in file order, numbered 1..N.

    Returned dict also carries the evidence needed to distrust it:
      chain_id, n_residues, n_ca, author_first/last, gaps, altlocs, warnings.

    A residue with no CA atom would still occupy an index here (we key on
    residue identity, not on CA), which is what parse_pdb does too; n_ca is
    reported separately so a mismatch is visible rather than silent.
    """
    text = pdb_bytes.decode("utf-8", errors="replace")
    order: list[tuple[str, int, str]] = []
    seen: set[tuple[str, int, str]] = set()
    chains: set[str] = set()
    altlocs: set[str] = set()
    n_ca = 0
    n_hetatm = 0
    n_seqres = 0
    for line in text.splitlines():
        tag = line[:6]
        if tag == "SEQRES":
            n_seqres += 1
            continue
        if tag == "HETATM":
            n_hetatm += 1
            continue
        if tag != "ATOM  ":
            continue
        ch = line[21]
        chains.add(ch)
        if chain_id is not None and ch != chain_id:
            continue
        altlocs.add(line[16])
        try:
            resseq = int(line[22:26])
        except ValueError:                      # pragma: no cover - malformed
            continue
        key = (ch, resseq, line[26])
        if key not in seen:
            seen.add(key)
            order.append(key)
        if line[12:16] == " CA ":
            n_ca += 1

    warnings: list[str] = []
    if chain_id is None and len(chains) > 1:
        warnings.append(
            f"structure has {len(chains)} chains {sorted(chains)} and no chain "
            "was requested; indices below run over ALL of them in file order, "
            "which is almost certainly not what the YAML's `include` selects"
        )
    if n_seqres:
        warnings.append(
            f"{n_seqres} SEQRES lines present. boltzgen's parse_pdb uses SEQRES "
            "for the canonical sequence when it exists, so indices may be over "
            "the FULL sequence including unobserved residues, not 1..N over "
            "observed ones. This mapping assumes observed-only and is unsafe here"
        )
    dirty = {a for a in altlocs if a not in (" ", "A")}
    if dirty:
        warnings.append(f"alternate location indicators present: {sorted(dirty)}")
    icodes = {k[2] for k in order if k[2] != " "}
    if icodes:
        warnings.append(
            f"insertion codes present {sorted(icodes)}; author residue numbers "
            "are not unique, so the author->index map is ambiguous"
        )
    if n_ca != len(order):
        warnings.append(
            f"{len(order)} residues but {n_ca} CA atoms; some residue lacks a CA "
            "and index alignment should be checked by hand"
        )

    nums = [k[1] for k in order]
    gaps = [(nums[i], nums[i + 1]) for i in range(len(nums) - 1)
            if nums[i + 1] != nums[i] + 1]
    mapping: dict[int, int] = {}
    for i, key in enumerate(order):
        mapping.setdefault(key[1], i + 1)
    return {
        "chain_id": chain_id if chain_id is not None else (
            sorted(chains)[0] if len(chains) == 1 else None),
        "chains_in_file": sorted(chains),
        "map": mapping,
        "n_residues": len(order),
        "n_ca": n_ca,
        "n_hetatm": n_hetatm,
        "n_seqres": n_seqres,
        "author_first": nums[0] if nums else None,
        "author_last": nums[-1] if nums else None,
        "gaps": gaps,
        "warnings": warnings,
    }


def map_pocket_residues(pdb_bytes: bytes, author_ids, chain_id: str | None = None):
    """(boltzgen indices, mapping-report). Raises if any id cannot be mapped.

    Refusing is deliberate. An unmapped or out-of-range index does not fail
    loudly inside boltzgen -- it addresses a residue that is not the one we
    meant -- so the only safe behaviour is to stop here.
    """
    report = residue_index_map(pdb_bytes, chain_id=chain_id)
    m = report["map"]
    missing = [int(r) for r in author_ids if int(r) not in m]
    if missing:
        raise ValueError(
            f"{len(missing)} of {len(list(author_ids))} pocket residues are not "
            f"present in the structure and cannot be mapped to a BoltzGen index: "
            f"{missing}. Author numbering in the file runs "
            f"{report['author_first']}-{report['author_last']} over "
            f"{report['n_residues']} observed residues with gaps {report['gaps']}. "
            "Refusing to build a design spec that points at residues that do not "
            "exist."
        )
    idx = [m[int(r)] for r in author_ids]
    bad = [i for i in idx if i < 1 or i > report["n_residues"]]
    if bad:                                     # pragma: no cover - impossible
        raise ValueError(f"mapped indices out of 1..{report['n_residues']}: {bad}")
    return idx, report


def build_design_yaml(
    target_filename: str,
    chain_id: str,
    binding_indices,
    binder_min: int,
    binder_max: int,
    designed_chain_id: str = "G",
) -> str:
    """The BoltzGen design spec, shaped like example/vanilla_peptide_with_target_binding_site.

    `sequence: 80..140` means "sample a designed length uniformly in [80,140]".
    Paths inside the YAML resolve relative to the YAML's own directory.
    """
    if designed_chain_id == chain_id:
        designed_chain_id = next(
            c for c in "GHIJKLMNOPQRSTUVWXYZ" if c != chain_id
        )
    binding = ",".join(str(i) for i in binding_indices)
    return (
        "# generated by scripts/modal_boltzgen.py -- do not hand-edit\n"
        "# residue indices below are 1-based BoltzGen/mmCIF indices over the\n"
        "# OBSERVED residues of the file, NOT PDB author numbering.\n"
        "entities:\n"
        f"  - protein:\n"
        f"      id: {designed_chain_id}\n"
        f"      sequence: {binder_min}..{binder_max}\n"
        f"  - file:\n"
        f"      path: {target_filename}\n"
        f"      include:\n"
        f"        - chain:\n"
        f"            id: {chain_id}\n"
        f"      binding_types:\n"
        f"        - chain:\n"
        f"            id: {chain_id}\n"
        f"            binding: {binding}\n"
        f'      structure_groups: "all"\n'
    )


# ==========================================================================
# container-side helpers (import stdlib only at call time)
# ==========================================================================
def _sh(cmd, cwd=None, timeout=None, env=None):
    import os
    import subprocess
    import time

    t0 = time.time()
    proc = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout,
        env={**os.environ, **(env or {})},
    )
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "seconds": round(time.time() - t0, 2),
        # trimmed: a boltzgen run's stdout is large and we return this to the
        # client. The full log is in the Modal dashboard.
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def _boltzgen_version() -> dict:
    out = {"import_ok": False, "version": None, "error": None, "cli_ok": None}
    try:
        from importlib.metadata import version

        out["version"] = version("boltzgen")
    except Exception as exc:                    # pragma: no cover
        out["error"] = f"{type(exc).__name__}: {exc}"
    try:
        import boltzgen  # noqa: F401

        out["import_ok"] = True
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    res = _sh(["boltzgen", "--help"], timeout=300)
    out["cli_ok"] = res["returncode"] == 0
    out["cli_stderr_tail"] = res["stderr_tail"][-600:]
    return out


def _supported_flags(subcommand: str) -> set:
    """Flags `boltzgen <sub> --help` actually advertises.

    Used instead of assuming --seed / --reuse / --cache exist, because every
    CLI claim available to us was read off argparse source, never from a real
    --help. If --help itself fails we return an empty set and callers degrade
    to passing nothing optional.
    """
    import re

    res = _sh(["boltzgen", subcommand, "--help"], timeout=300)
    if res["returncode"] != 0:
        return set()
    return set(re.findall(r"--[A-Za-z0-9][A-Za-z0-9_-]*", res["stdout_tail"] or ""))


def _volume_state(root: str) -> dict:
    import os

    total = 0
    files = []
    if not os.path.isdir(root):
        return {"root": root, "exists": False, "n_files": 0, "bytes": 0, "files": []}
    for dirpath, _dirs, names in os.walk(root):
        for n in names:
            p = os.path.join(dirpath, n)
            try:
                sz = os.path.getsize(p)
            except OSError:
                continue
            total += sz
            if len(files) < 200:
                files.append({"path": os.path.relpath(p, root), "bytes": sz})
    return {
        "root": root,
        "exists": True,
        "n_files": sum(len(n) for _d, _s, n in os.walk(root)),
        "bytes": total,
        "gib": round(total / 1024 ** 3, 3),
        "files": sorted(files, key=lambda f: -f["bytes"])[:40],
    }


def _gpu_report() -> dict:
    rep = {"nvidia_smi": None, "nvidia_smi_error": None,
           "torch_cuda": None, "capability": None, "kernels_eligible": None}
    res = _sh(["nvidia-smi",
               "--query-gpu=name,memory.total,driver_version",
               "--format=csv,noheader"], timeout=120)
    if res["returncode"] == 0:
        rep["nvidia_smi"] = res["stdout_tail"].strip()
    else:
        rep["nvidia_smi_error"] = (res["stderr_tail"] or "").strip()[:400]
    try:
        import torch

        rep["torch_version"] = torch.__version__
        rep["torch_cuda"] = bool(torch.cuda.is_available())
        if rep["torch_cuda"]:
            cap = torch.cuda.get_device_capability(0)
            rep["capability"] = list(cap)
            # boltzgen's `--use_kernels auto` enables cuEquivariance kernels at
            # device capability >= 8 (A100/H100 yes; T4/V100 no).
            rep["kernels_eligible"] = cap[0] >= 8
            rep["device_name"] = torch.cuda.get_device_name(0)
            rep["vram_bytes"] = torch.cuda.get_device_properties(0).total_memory
    except Exception as exc:
        rep["torch_error"] = f"{type(exc).__name__}: {exc}"
    return rep


def _run_download(cache: str) -> dict:
    """`boltzgen download` -- argument shape is NOT verified, so try forms.

    Recorded honestly: every attempt and its exit code is returned.
    """
    attempts = []
    for cmd in (
        ["boltzgen", "download", "all", "--cache", cache],
        ["boltzgen", "download", "--cache", cache],
        ["boltzgen", "download", "all"],
    ):
        res = _sh(cmd, timeout=90 * 60)
        attempts.append(res)
        if res["returncode"] == 0:
            return {"ok": True, "attempts": attempts}
    return {"ok": False, "attempts": attempts}


def _read_csv(path) -> dict:
    import csv

    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    cols = list(rows[0].keys()) if rows else []
    return {
        "path": str(path),
        "n_rows": len(rows),
        "columns": cols,
        "expected_confidence_keys_present": sorted(
            set(cols) & set(EXPECTED_CONFIDENCE_KEYS)),
        "expected_confidence_keys_missing": sorted(
            set(EXPECTED_CONFIDENCE_KEYS) - set(cols)),
        "rows": rows,
    }


def _collect_outputs(out_dir, return_structures: bool, max_return_bytes: int) -> dict:
    """Find the designed COMPLEX structures and the metrics CSVs.

    Preference order, and why: `refold_cif/` holds the refolded target+binder
    complex (one mmCIF per design) and is what final_<budget>_designs/ copies
    from. `intermediate_designs_inverse_folded/*.cif` must NOT be used for
    contacts -- designed residues there carry backbone atoms only, sidechains
    are literally 0,0,0. `intermediate_designs/` is the design-step-only output
    and does contain target+binder, which is what a `--steps design` run leaves.

    Nothing is assumed to exist: whatever is found is reported, with counts.
    """
    from pathlib import Path

    out_dir = Path(out_dir)
    candidates = [
        ("final_designs", sorted(out_dir.glob("final_ranked_designs/final_*_designs"))),
        ("refold_cif", sorted(out_dir.glob("**/refold_cif"))),
        ("intermediate_ranked",
         sorted(out_dir.glob("final_ranked_designs/intermediate_ranked_*_designs"))),
        ("intermediate_designs", sorted(out_dir.glob("intermediate_designs"))),
    ]
    inventory = {}
    chosen_kind, chosen_dir = None, None
    for kind, dirs in candidates:
        found = [d for d in dirs if d.is_dir()]
        n = sum(len(list(d.glob("*.cif"))) for d in found)
        inventory[kind] = {"dirs": [str(d.relative_to(out_dir)) for d in found],
                           "n_cif": n}
        if chosen_dir is None and n:
            chosen_kind, chosen_dir = kind, found[0]

    structures, returned_bytes, truncated = {}, 0, False
    if chosen_dir is not None and return_structures:
        for cif in sorted(chosen_dir.glob("*.cif")):
            raw = cif.read_bytes()
            if returned_bytes + len(raw) > max_return_bytes:
                truncated = True
                break
            structures[cif.name] = raw.decode("utf-8", errors="replace")
            returned_bytes += len(raw)

    metrics = []
    for csv_path in sorted(out_dir.glob("**/*metrics*.csv")):
        try:
            metrics.append(_read_csv(csv_path))
        except Exception as exc:
            metrics.append({"path": str(csv_path),
                            "error": f"{type(exc).__name__}: {exc}"})

    return {
        "structure_source": chosen_kind,
        "structure_dir": str(chosen_dir.relative_to(out_dir)) if chosen_dir else None,
        "structure_inventory": inventory,
        "n_structures_returned": len(structures),
        "structures_truncated": truncated,
        "structures": structures,
        "metrics_csvs": metrics,
        "sidechain_warning": (
            "contacts must come from the refolded COMPLEX; designed residues in "
            "intermediate_designs_inverse_folded/*.cif have backbone atoms only"
        ),
    }


def _prepare_workdir(target_pdb: bytes, target_name, pocket_residue_ids,
                     chain_id, binder_min, binder_max, workdir):
    from pathlib import Path

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    target_path = workdir / f"{target_name}.pdb"
    target_path.write_bytes(target_pdb)
    indices, report = map_pocket_residues(target_pdb, pocket_residue_ids,
                                          chain_id=chain_id)
    yaml_text = build_design_yaml(
        target_path.name, report["chain_id"] or "A", indices, binder_min, binder_max)
    yaml_path = workdir / f"{target_name}.yaml"
    yaml_path.write_text(yaml_text)
    mapping = {
        "author_to_boltzgen_index": {int(a): int(b) for a, b in
                                     zip(pocket_residue_ids, indices)},
        "n_observed_residues": report["n_residues"],
        "n_ca": report["n_ca"],
        "chain_id": report["chain_id"],
        "author_range": [report["author_first"], report["author_last"]],
        "gaps": report["gaps"],
        "warnings": report["warnings"],
        "numbering_note": (
            "pocket ids in, author numbering (auth_seq_id, what pocket.json and "
            "scripts/interfaces.py use); binding: ids out, 1-based BoltzGen "
            "index over observed residues. UNVERIFIED against a real boltzgen "
            "parse -- run check_spec() before spending GPU time."
        ),
    }
    return yaml_path, target_path, yaml_text, mapping


# ==========================================================================
# Modal functions
# ==========================================================================
@app.function(
    volumes={WEIGHTS_MOUNT: weights_volume},
    timeout=30 * 60,          # generous: first call also pays the image pull
    cpu=2.0,
    memory=8192,
    max_containers=1,
)
def smoke() -> dict:
    """Cheapest possible validation: does the image build and boltzgen import?

    CPU only, so `nvidia-smi` is absent BY DESIGN and its absence here is not a
    finding -- use smoke_gpu() for the GPU probe.
    """
    import os
    import platform
    import shutil
    import sys

    print("smoke: starting", flush=True)
    bg = _boltzgen_version()
    state = _volume_state(BOLTZGEN_CACHE)
    du = shutil.disk_usage("/")
    result = {
        "kind": "smoke-cpu",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "boltzgen": bg,
        "boltzgen_run_flags": sorted(_supported_flags("run")),
        "boltzgen_subcommands_expected": ["run", "execute", "configure",
                                          "download", "check", "merge"],
        "weights_volume": state,
        "weights_marker_present": os.path.exists(WEIGHTS_MARKER),
        "hf_home": os.environ.get("HF_HOME"),
        "hf_token_set": bool(os.environ.get("HF_TOKEN")),
        "disk_free_gib": round(du.free / 1024 ** 3, 1),
        "gpu": {"note": "CPU container: nvidia-smi is expected to be absent"},
    }
    print(f"smoke: boltzgen={bg['version']} import_ok={bg['import_ok']} "
          f"weights={state['gib'] if state['exists'] else 0} GiB", flush=True)
    return result


@app.function(
    gpu=DEFAULT_GPU,
    volumes={WEIGHTS_MOUNT: weights_volume},
    timeout=30 * 60,
    cpu=2.0,
    memory=16384,
    max_containers=1,
)
def smoke_gpu() -> dict:
    """GPU probe: nvidia-smi, VRAM, and whether `--use_kernels auto` will fire.

    Call it with .with_options(gpu=...) to probe a specific tier. This is the
    only way to settle the VRAM question -- no primary source publishes a
    memory-vs-length figure for BoltzGen, and the maintainers say so.
    """
    import sys

    rep = _gpu_report()
    print("smoke_gpu:", rep.get("nvidia_smi") or rep.get("nvidia_smi_error"),
          flush=True)
    return {
        "kind": "smoke-gpu",
        "python": sys.version.split()[0],
        "gpu": rep,
        "boltzgen": _boltzgen_version(),
        "weights_volume": _volume_state(BOLTZGEN_CACHE),
    }


@app.function(
    volumes={WEIGHTS_MOUNT: weights_volume},
    timeout=2 * 60 * 60,
    cpu=4.0,
    memory=16384,
    max_containers=1,
)
def fetch_weights(force: bool = False) -> dict:
    """Download ~6.35 GB of checkpoints into the Volume, once.

    On a CPU container on purpose: GPU seconds must never be spent downloading,
    and the default function timeout of 300 s would otherwise kill a cold run
    before any design work started.

    protein-anything needs: boltzgen1_diverse 1.93 + boltzgen1_adherence 1.93 +
    boltzgen1_ifold 0.01 + boltz2_conf_final 2.09 GB, plus mols.zip 0.39 GB.
    All from the public HF repo boltzgen/boltzgen-1 (gated=False, MIT) and the
    dataset repo boltzgen/inference-data -- no token, no licence click.
    """
    import json
    import os
    import time

    os.makedirs(BOLTZGEN_CACHE, exist_ok=True)
    before = _volume_state(BOLTZGEN_CACHE)
    if os.path.exists(WEIGHTS_MARKER) and not force:
        print("fetch_weights: already complete, skipping", flush=True)
        return {"status": "cached", "weights_volume": before,
                "marker": json.load(open(WEIGHTS_MARKER))}

    t0 = time.time()
    res = _run_download(BOLTZGEN_CACHE)
    after = _volume_state(BOLTZGEN_CACHE)
    record = {
        "status": "downloaded" if res["ok"] else "failed",
        "seconds": round(time.time() - t0, 1),
        "boltzgen_version": _boltzgen_version()["version"],
        "cache": BOLTZGEN_CACHE,
        "bytes_before": before.get("bytes", 0),
        "bytes_after": after.get("bytes", 0),
        "attempts": [{"cmd": a["cmd"], "returncode": a["returncode"],
                      "seconds": a["seconds"]} for a in res["attempts"]],
    }
    if res["ok"]:
        with open(WEIGHTS_MARKER, "w") as fh:
            json.dump(record, fh, indent=2)
    weights_volume.commit()
    record["weights_volume"] = _volume_state(BOLTZGEN_CACHE)
    record["last_stderr_tail"] = res["attempts"][-1]["stderr_tail"]
    return record


@app.function(
    volumes={WEIGHTS_MOUNT: weights_volume, RUNS_MOUNT: runs_volume},
    timeout=60 * 60,
    cpu=4.0,
    memory=16384,
    max_containers=1,
)
def check_spec(
    target_pdb: bytes,
    target_name: str,
    pocket_residue_ids: list,
    chain_id: str | None = None,
    binder_min: int = 80,
    binder_max: int = 140,
    run_key: str = "check",
) -> dict:
    """`boltzgen check` on the generated YAML -- the guard before GPU money.

    boltzgen check loads no model checkpoint (it takes only --moldir and the
    download options), which is why this is a CPU function. That it truly runs
    without a GPU is UNVERIFIED; if it turns out to need one, this call fails
    here rather than mid-design-run, which is still the cheap failure.

    The emitted mmCIF renders the binding site in a different colour; it is
    copied to the runs Volume so the author->index mapping can be eyeballed.
    """
    import json
    import shutil
    from pathlib import Path

    work = Path("/tmp/check")
    if work.exists():
        shutil.rmtree(work)
    yaml_path, _target, yaml_text, mapping = _prepare_workdir(
        target_pdb, target_name, pocket_residue_ids, chain_id,
        binder_min, binder_max, work)

    out = Path(RUNS_MOUNT) / run_key / "check_out"
    out.mkdir(parents=True, exist_ok=True)
    flags = _supported_flags("check")
    cmd = ["boltzgen", "check", str(yaml_path), "--output", str(out)]
    if "--cache" in flags:
        cmd += ["--cache", BOLTZGEN_CACHE]
    res = _sh(cmd, cwd=str(work), timeout=45 * 60)
    emitted = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
    (Path(RUNS_MOUNT) / run_key / "design_spec.yaml").write_text(yaml_text)
    (Path(RUNS_MOUNT) / run_key / "residue_mapping.json").write_text(
        json.dumps(mapping, indent=2))
    runs_volume.commit()
    return {
        "kind": "check",
        "ok": res["returncode"] == 0,
        "command": res["cmd"],
        "returncode": res["returncode"],
        "seconds": res["seconds"],
        "stdout_tail": res["stdout_tail"],
        "stderr_tail": res["stderr_tail"],
        "yaml": yaml_text,
        "residue_mapping": mapping,
        "emitted_files": emitted[:200],
        "runs_volume_path": f"{run_key}/check_out",
        "check_flags_advertised": sorted(flags),
    }


@app.function(
    gpu=DEFAULT_GPU,
    volumes={WEIGHTS_MOUNT: weights_volume, RUNS_MOUNT: runs_volume},
    timeout=6 * 60 * 60,      # a 1000-design protein-anything run is ~26 GPU-h,
                              # so large N must be split into batches by the
                              # driver; 6 h is a batch ceiling, not a run ceiling
    startup_timeout=45 * 60,  # image pull + weight staging, kept off the clock
    cpu=DESIGN_CPU,
    memory=DESIGN_MEMORY_MB,
    retries=0,                # retries restart the timeout per attempt; a
                              # retried GPU run is a doubled bill. --reuse makes
                              # a manual re-run resume instead.
    max_containers=1,
)
def design(
    target_pdb: bytes,
    target_name: str,
    pocket_residue_ids: list,
    protocol: str = "protein-anything",
    num_designs: int = 50,
    seed: int | None = 0,
    chain_id: str | None = None,
    binder_min: int = 80,
    binder_max: int = 140,
    steps: list | None = None,
    diffusion_batch_size: int | None = None,
    use_kernels: str = "auto",
    reuse: bool = True,
    run_key: str = "run",
    return_structures: bool = True,
    max_return_bytes: int = 32 * 1024 * 1024,
    extra_args: list | None = None,
) -> dict:
    """Run BoltzGen against a pocket and return designs + confidence scores.

    `pocket_residue_ids` are AUTHOR residue numbers; the mapping to BoltzGen's
    1-based observed-residue indices happens here (see module docstring).

    Everything is written to the runs Volume under `run_key` as it goes, so a
    timeout or a lost connection does not lose the designs; `--reuse` then makes
    a re-invocation continue rather than restart. Only JSON-able types cross the
    boundary (the local client is 3.14, the container is 3.12).
    """
    import json
    import time
    from pathlib import Path

    t_start = time.time()
    run_dir = Path(RUNS_MOUNT) / run_key
    work = run_dir / "spec"
    out_dir = run_dir / "boltzgen_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    yaml_path, _target, yaml_text, mapping = _prepare_workdir(
        target_pdb, target_name, pocket_residue_ids, chain_id,
        binder_min, binder_max, work)
    (run_dir / "design_spec.yaml").write_text(yaml_text)
    (run_dir / "residue_mapping.json").write_text(json.dumps(mapping, indent=2))
    runs_volume.commit()

    gpu = _gpu_report()
    bg = _boltzgen_version()
    flags = _supported_flags("run")
    print(f"design: boltzgen={bg['version']} gpu={gpu.get('nvidia_smi')} "
          f"protocol={protocol} n={num_designs}", flush=True)

    cmd = ["boltzgen", "run", str(yaml_path),
           "--output", str(out_dir),
           "--protocol", protocol,
           "--num_designs", str(num_designs)]
    unsupported = []

    def _add(flag, *values):
        # Never pass a flag this build does not advertise: an unknown flag is
        # an argparse error six minutes into a GPU container.
        if not flags or flag in flags:
            cmd.extend([flag, *values])
            return True
        unsupported.append(flag)
        return False

    seed_applied = False
    if seed is not None:
        seed_applied = _add("--seed", str(seed))
    if reuse:
        _add("--reuse")
    _add("--cache", BOLTZGEN_CACHE)
    if diffusion_batch_size:
        _add("--diffusion_batch_size", str(diffusion_batch_size))
    if use_kernels:
        _add("--use_kernels", str(use_kernels))
    if steps:
        bad = [s for s in steps if s not in PIPELINE_STEPS]
        if bad:
            raise ValueError(f"unknown pipeline steps {bad}; known: {PIPELINE_STEPS}")
        if not flags or "--steps" in flags:
            cmd += ["--steps", *steps]
        else:
            unsupported.append("--steps")
    if extra_args:
        cmd += list(extra_args)

    res = _sh(cmd, cwd=str(work), timeout=None)
    runs_volume.commit()
    collected = _collect_outputs(out_dir, return_structures, max_return_bytes)

    manifest = {
        "kind": "design",
        "ok": res["returncode"] == 0,
        "run_key": run_key,
        "command": res["cmd"],
        "returncode": res["returncode"],
        "wall_seconds": round(time.time() - t_start, 1),
        "boltzgen_seconds": res["seconds"],
        "boltzgen": bg,
        "gpu": gpu,
        "protocol": protocol,
        "num_designs_requested": num_designs,
        "seed": seed,
        "seed_applied": seed_applied,
        "seed_note": (
            "boltzgen run does not advertise --seed in this build; the run is "
            "NOT seeded and is not bit-reproducible"
            if (seed is not None and not seed_applied) else None
        ),
        "steps": steps,
        "binder_length": [binder_min, binder_max],
        "diffusion_batch_size": diffusion_batch_size,
        "flags_not_supported_by_this_build": unsupported,
        "residue_mapping": mapping,
        "design_spec_yaml": yaml_text,
        "runs_volume_path": run_key,
        "structure_source": collected["structure_source"],
        "structure_inventory": collected["structure_inventory"],
        "metrics_csvs": [{k: v for k, v in m.items() if k != "rows"}
                         for m in collected["metrics_csvs"]],
        "stdout_tail": res["stdout_tail"],
        "stderr_tail": res["stderr_tail"],
    }
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    runs_volume.commit()

    manifest["metrics"] = collected["metrics_csvs"]
    manifest["structures"] = collected["structures"]
    manifest["n_structures_returned"] = collected["n_structures_returned"]
    manifest["structures_truncated"] = collected["structures_truncated"]
    return manifest


@app.local_entrypoint()
def main(what: str = "smoke"):
    """`modal run scripts/modal_boltzgen.py --what smoke|smoke-gpu|weights`.

    The real entry point is scripts/run_design.py, which carries the cost guard,
    the manifest writing and the resume logic. This exists so the app can be
    exercised with the bare modal CLI.
    """
    import json

    if what == "smoke":
        print(json.dumps(smoke.remote(), indent=2)[:8000])
    elif what == "smoke-gpu":
        print(json.dumps(smoke_gpu.remote(), indent=2)[:8000])
    elif what == "weights":
        print(json.dumps(fetch_weights.remote(), indent=2)[:8000])
    else:
        raise SystemExit(f"unknown --what {what!r}; use smoke, smoke-gpu or weights")
