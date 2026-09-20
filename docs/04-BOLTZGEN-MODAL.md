# 04 - The BoltzGen arm on Modal

Written **before** the first paid run, per `CLAUDE.md`: hard constraints stated up
front, nothing reported as measured that was not measured. Every number below is
labelled either **computed here**, **read from a published source**, or
**ESTIMATE**. Verified 2026-09-19.

## What this is for

`PROJECT_GOAL.md` §8.3 calls one comparison decisive, and §7 makes it the M2 gate:

> does a consensus interface signature built from **BoltzGen designs** carry
> information beyond the pocket **P2Rank** finds in 7.1 s?

The P2Rank arm is built and measured (`scripts/run_p2rank.py`,
`results/p2rank_pockets.json`, `results/pipeline/colorectal-cancer/target/pocket.json`).
This document covers the other arm: generate a design ensemble against the same
site, reduce each design to target-side contacts, and aggregate them into the
`InterfaceSignature` of §4.3 with `source="boltzgen_consensus"`.

Two files implement it:

| Path | Role |
|---|---|
| `scripts/modal_boltzgen.py` | the Modal app: image, weights Volume, runs Volume, and five functions (`smoke`, `smoke_gpu`, `fetch_weights`, `check_spec`, `design`) |
| `scripts/run_design.py` | the local driver: credentials, plan, cost guard, batching, resume, manifest |

`scripts/interfaces.py` already turns a complex into target-side contacts and
aggregates a list of them into a signature. The design run's job is only to hand
it **refolded target+binder complexes**, one per design.

## Why this cannot run on this machine

Three independent blockers, each verified rather than assumed:

| Blocker | Evidence |
|---|---|
| No NVIDIA GPU | Apple silicon; `docs/03-SCOPE-AND-CONSTRAINTS.md` |
| `env/` is Python 3.14.7, boltzgen pins `numpy==2.0.2` | numpy 2.0.2 publishes wheels for cp39-cp312 only; on 3.13/3.14 pip would build it from source. boltzgen's declared `requires-python = ">=3.11"` has **no upper bound** and is misleading: the effective ceiling is **3.12** |
| boltzgen needs `cuequivariance_ops_cu12`, `cuequivariance_ops_torch_cu12`, `cuequivariance_torch` | all three publish **zero** macOS wheels (manylinux x86_64 / aarch64 only). No local Python version fixes this |

The third is the strongest: boltzgen cannot be installed on this Mac at any
Python version. On Modal we choose the image, so the runner pins
`python_version="3.12"` explicitly — `modal.Image.debian_slim()` with no
`python_version` copies the **local** interpreter (3.14 here), which would
reproduce the failure remotely and late.

Also verified and relevant: boltzgen does **not** depend on the `boltz` PyPI
package (whose `>=3.10,<3.13` constraint is therefore irrelevant here); it
vendors its model code and pulls Boltz-2 *weights* from HuggingFace. The weight
repos `boltzgen/boltzgen-1` and `boltzgen/inference-data` are public
(`gated: False`, MIT), so **no HF token and no licence click** is needed.

## The residue-index trap

This is the single highest-risk item in the whole arm and the runner is built
around it.

BoltzGen design YAMLs take **1-based canonical mmCIF indices**, not PDB author
numbering. For a `.pdb` with no SEQRES, boltzgen's `parse_pdb` builds the
sequence from the residues physically present and assigns `label_seq = j + 1`,
i.e. **1..N over observed residues**.

Computed here from the actual file
(`results/pipeline/colorectal-cancer/target/structure/3VHE_A_stripped.pdb`):
**0 SEQRES lines, one chain A, 303 observed residues = 303 CA atoms, no altlocs,
no insertion codes, no HETATM, author numbering 811-1169 with gaps after 939 and
after 1062.** So valid BoltzGen indices are **1-303**, and every one of the 17
P2Rank pocket residue ids (840 … 1047) is out of range. Pasting them into
`binding:` would not mis-target quietly — it would address residues that do not
exist.

`residue_index_map()` computes the mapping from the same bytes the container
parses, and `map_pocket_residues()` **refuses to build a spec** if any requested
residue is unmapped. For our pocket it yields:

```
840->30  848->38  866->56  868->58  885->75  889->79  899->89
914->104 916->106 918->108 919->109 922->112 923->113
1035->171 1045->181 1046->182 1047->183
```

**This mapping is derived from reading `pdb_parser.py`, not confirmed by running
BoltzGen.** `run_design.py check` exists to settle it on a CPU container before
any GPU money is spent; it writes the emitted mmCIF to the runs Volume, where the
binding site renders in a different colour.

## Cost

All cost figures are **ESTIMATES**. The arithmetic is computed here; two of its
three inputs are not ours and one is a guess:

1. Per-stage seconds/design come from BoltzGen's Figure 12
   (`assets/fig_seconds_per_design.png`, "on a single NVIDIA A100 GPU"). There is
   **no published numeric table** — the values were read off the plot, so ±10%.
2. `design_folding` (refold the binder alone; `protein-anything` runs it, Figure
   12 does not plot it) is **extrapolated** as the refold curve evaluated at the
   binder length. ±20 s/design of extra uncertainty on its own.
3. Per-GPU slowdown factors vs A100 are **guesses** from bandwidth/throughput
   ratios. Nobody has published BoltzGen timings on any GPU but A100. A GPU with
   no factor (T4, H200, B200, B300, RTX-PRO-6000) is **refused** rather than
   guessed at.

System size, computed: 303 target residues + a mean designed binder of 110
(`sequence: 80..140`) = **413 residues**. Interpolating Figure 12 there gives, per
design on A100: design 10.03 s, inverse_folding 2.27 s, folding 63.20 s,
design_folding 22.59 s (extrapolated) = **98.1 GPU-s** plus 1.8 CPU-s.

Modal list prices, fetched 2026-09-19. The billed rate for the design function is
GPU + 8 cores + 64 GiB = **$0.00082988/s** on A100-40GB. 64 GiB is not optional:
boltzgen issue #208 OOM'd in the *analysis* step on an A100 with 16 GiB **host**
RAM at only 100 designs, and upstream's own SLURM example asks for `--mem=64G`.

`protein-anything`, full pipeline, A100-40GB — printed by `--ladder`:

| N | GPU-hours | compute | $/design |
|---|---|---|---|
| 50 | 1.39 | $4.16 | $0.0832 |
| 100 | 2.78 | $8.31 | $0.0831 |
| 250 | 6.94 | $20.74 | $0.0830 |
| 500 | 13.88 | $41.47 | $0.0829 |
| 1000 | 27.76 | $82.92 | $0.0829 |

Plus a one-off **$0.18** for the first image build and the 6.35 GB weight
download (both estimated wall-clock; neither measured).

`--steps design` only — the design step emits `intermediate_designs/*.cif`, which
per the README contains "CIF and NPZ for the designed proteins **and targets**",
i.e. the full complex `scripts/interfaces.py` needs:

| N | GPU-hours | compute | $/design |
|---|---|---|---|
| 50 | 0.14 | $0.43 | $0.0087 |
| 1000 | 2.79 | $8.34 | $0.0083 |

**~10x cheaper**, at the price of no refolding, therefore no `design_iptm` /
RMSD, therefore no filtering — the consensus would be over *all* designs rather
than the few percent that survive. Both arms are worth running at the same N:
"does filtering change the consensus contact map?" is itself a cheap I6.1 result.

Cheapest GPU is not cheapest per design, computed at N=1000:

| GPU | slowdown (GUESS) | $/design | total |
|---|---|---|---|
| A100-40GB | x1.0 | $0.0829 | $82.92 |
| H100 | x0.6 | $0.0815 | $81.55 |
| L40S | x1.6 | $0.1253 | $125.26 |
| A10 | x2.75 | $0.1502 | $150.16 |
| L4 | x3.5 | $0.1618 | $161.84 |

A100-40GB is the default: cheapest per design, the GPU the published timings were
measured on, no VRAM question, and it avoids the L40S `pynvml.NVMLError_NotSupported`
crash reported in boltzgen issue #11.

## The cost guard (PROJECT_GOAL.md rule 10)

`run_design.py` estimates before it submits and **refuses** above `--budget-usd`
(default **$10**) unless `--force` is passed, exiting 3. A bare invocation with no
subcommand refuses too — only `--dry-run` works without one, so nothing can spend
money by accident. Real output:

```
$ ./env/bin/python scripts/run_design.py design --num-designs 1000
TOTAL ESTIMATE    $83.11   vs budget $10.0
=== refused by the cost guard =======================================
estimated $83.11 > --budget-usd 10.0. Nothing was submitted.
Options: lower --num-designs, use --steps design (roughly 10x cheaper, no
confidence scores), raise --budget-usd, or pass --force if you have decided to
spend it.                                                      # exit code 3
```

## Exact commands

```bash
cd /Users/evanxiang/Desktop/Projects/drug-similarity-discovery

# 0. credentials. Either works; .env wins. .env is gitignored.
printf 'MODAL_TOKEN_ID=ak-...\nMODAL_TOKEN_SECRET=as-...\n' >> .env
# or: ./env/bin/modal token set --token-id ak-... --token-secret as-...

# 1. plan and price, submits nothing, needs no credentials
./env/bin/python scripts/run_design.py --dry-run --num-designs 100 --ladder

# 2. a few cents: does the image build and does boltzgen import?
./env/bin/python scripts/run_design.py smoke

# 3. GPU probe: nvidia-smi, VRAM, and whether --use_kernels auto will fire
./env/bin/python scripts/run_design.py smoke --gpu-probe --gpu A100-40GB

# 4. one-off 6.35 GB of checkpoints into the persistent Volume (CPU container)
./env/bin/python scripts/run_design.py weights

# 5. THE GUARD ON THE INDEX TRAP: boltzgen check on the generated YAML, no GPU
./env/bin/python scripts/run_design.py check
./env/bin/modal volume get boltzgen-runs colorectal-cancer/<hash> ./check_out

# 6. the real run, smallest useful N first (README: validate the spec at 50)
./env/bin/python scripts/run_design.py design --num-designs 50 --budget-usd 6

# 7. the convergence ladder (PROJECT_GOAL C5). --reuse makes each rung cost only
#    the increment, so the whole ladder costs about what the top rung costs.
for n in 100 250 500 1000; do
  ./env/bin/python scripts/run_design.py design --num-designs $n --budget-usd 90 --force
done

# the bare modal CLI also works for the cheap probes
./env/bin/modal run scripts/modal_boltzgen.py --what smoke
```

Outputs land in `results/pipeline/<disease>/designs/<site>__<protocol>__n<N>__seed<S>__<spechash>/`:
`batch-0000/designs/*.cif`, the metrics CSVs, `batch-0000/result.json`,
`design_spec.yaml`, and a run-level `manifest.json` recording the boltzgen
version, the Modal app and image Python, the GPU string, the GPU actually
reported by `nvidia-smi`, the seed and **whether the seed was actually applied**,
the structure sha256, the pocket residue ids, the wall time and a cost computed
from that wall time. A batch whose `result.json` says `ok` is never re-run
(`--force-rerun` overrides). The full BoltzGen output tree stays on the
`boltzgen-runs` Volume and is fetched with `modal volume get`.

## Still UNVERIFIED — nothing below has been run

No Modal job has been submitted by this work. Credentials do exist locally
(`~/.modal.toml`, profile `evankxiang`), and `docs/03` records that CPU functions
run while **every GPU tier was refused pending a payment method**. So the design
run may still be blocked for a reason that has nothing to do with this code.

| Item | Why it is open | What settles it |
|---|---|---|
| The author→BoltzGen index mapping | derived by reading `pdb_parser.py` + counting CA records, not by running boltzgen | `run_design.py check`, then look at the emitted mmCIF |
| That `boltzgen check` needs no GPU | it loads no checkpoint (only `--moldir`), which is inference, not observation | running `check` |
| That boltzgen 0.3.2 installs on `debian_slim(python_version="3.12")` | never attempted | `run_design.py smoke` |
| `boltzgen download` argument shape, and whether `--cache` is honoured everywhere | the CLI was read off argparse source; no real `--help` output exists anywhere in this project | `weights`; `_run_download()` tries three forms and records every exit code |
| Whether `boltzgen run` accepts `--seed` | not advertised in any source we have | `smoke` prints the advertised flags; the manifest records `seed_applied` and refuses to describe an unseeded run as seeded |
| Metrics CSV column names | the metric dict **keys** were verified in `const.eval_keys_confidence`, but no real CSV was seen | the run; the CSVs are read as-is and expected keys reported found/missing |
| VRAM needed for a 413-residue system | **no primary source publishes one**; the maintainers say so explicitly (issue #32). A 16 GB T4 completed a ~300-residue run; a third-party skill names 24 GB as a floor to avoid | `smoke --gpu-probe`, then a real run |
| Whether MSA generation is triggered for the target, and whether it needs network access in the container | the file-level `msa` flag defaults to 0 = "automatic"; that path was never traced | the first `check` or `design` run |
| Per-GPU slowdown factors | guesses; only A100 timings are published | measuring one batch on two GPUs |
| Actual filter pass rate for KDR | the 3/100 figure is AWS's, on a different and smaller target | the first full-pipeline run |
| Modal's acceptance of the exact GPU strings | the client does **no** client-side validation — a wrong string is rejected server-side after the image builds. Current docs list `A10`, **not** `A10G` | any submission |

## Not evaluated / why

| Not evaluated | Why |
|---|---|
| `peptide-anything`, `nanobody-anything` (PROJECT_GOAL 6.C2 wants all three modalities) | the runner accepts `--protocol`, but the M2 gate only needs one modality; multi-modality is a cost multiple of an arm that has not run once |
| `protein-small_molecule` | designs a protein that *binds* a small molecule — the reverse of what a target-side signature needs. The driver refuses it without `--force` |
| The Rowan backend (`PROJECT_GOAL` C1's other option) | no Rowan key in the repo; `.env.example` is a template only |
| Affinity step / `boltz2_aff.ckpt` | only runs under `protein-small_molecule`, so it is neither downloaded nor timed |
| foldseek-based novelty and clustering | `analysis.yaml` ships `novelty_*: false` and `run_clustering: false`, so upstream's hardcoded internal foldseek paths are never reached. Left off |
| A bit-reproducible run | contingent on `--seed` existing in this build (see above). Until then the manifest says so rather than implying determinism |
| Contact extraction from the designs | `scripts/interfaces.py` already does it and is out of this task's scope; the runner's contract is to deliver refolded **complexes** (`refold_cif/`, or `intermediate_designs/` for a `--steps design` run) — never `intermediate_designs_inverse_folded/*.cif`, where designed sidechains are literally `0,0,0` |
| Anything about how many designs the consensus needs | no published measurement exists; the preprint has no convergence, epitope-adherence or contact-recovery analysis. That is PROJECT_GOAL C5's experiment, and it is what this runner exists to make affordable |

## What was validated, with no token

Run offline, on this machine, today:

- Both scripts parse (`ast.parse`) and `scripts/modal_boltzgen.py` **imports with
  no credentials present** — the `App`, `Image` and both `Volume`s construct, and
  all five functions expose `.remote` / `.with_options`. No side effects at import.
- `run_design.py --dry-run` prints the plan, the residue mapping and the cost with
  no credentials, and submits nothing.
- With credentials deliberately hidden (`MODAL_CONFIG_PATH=/nonexistent`),
  `run_design.py smoke` exits **2** with the "put these two lines in .env" message
  and no traceback.
- The cost guard refuses N=1000 at the default budget with exit **3**.
- Unknown GPU strings (`A10G`) and GPUs with no published slowdown (`T4`) are
  refused with an explanation rather than estimated.
- The residue mapping reproduces independently of the research note, and
  `map_pocket_residues` refuses an id that is not in the structure. Synthetic
  inputs confirm the multi-chain, SEQRES-present, insertion-code and empty-file
  warnings all fire.
- The cost arithmetic reproduces the published-price cross-checks: the billed rate
  is $0.00082988/s, and at the researchers' 93.2 s/design basis N=1000 is $77.34,
  matching their independently computed $77.

**Not validated, and cannot be from here:** everything that requires Modal to
actually run. No claim in this document about what BoltzGen does on a GPU is a
measurement.
