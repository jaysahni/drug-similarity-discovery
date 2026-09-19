# 00 - Environment

Verified 2026-09-19 on the primary dev machine.

## Hardware

- Apple M5 Pro, 24 GB unified memory, macOS arm64 (Darwin 25.6.0)
- **No CUDA GPU.** Anything needing fast neural inference either runs slowly
  on CPU or goes through a hosted API (Rowan).

## Python

- `env/` = python.org **Python 3.14.7** venv. No conda, no Homebrew.
- `pip install` works; cp314 arm64 wheels exist for rdkit, numpy, pandas,
  scikit-learn, torch.
- **Hard constraint: `boltz` (Boltz-1/2) requires Python >=3.10,<3.13** — it
  cannot be installed in `env/`. Local Boltz-2 would need a separate Python
  3.12 install + venv, and would run CPU-only. Preferred path is hosted
  Boltz-2 via Rowan (see `02-PIPELINE-PLAN.md`).
- `rowan-python` requires Python >=3.12 — installs fine in `env/`.

## Network / secrets

- The python.org 3.14 build ships no CA bundle; HTTP scripts must use
  certifi (see `scripts/fetch_chembl.py`).
- `ROWAN_API_KEY` goes in `.env` (gitignored). Free tier: 500 credits on
  signup + 20/week. CPU jobs 1 credit/min, GPU jobs 3 credits/min
  (7/min on H100/H200).

## Remote

- `origin` = github.com/jaysahni/drug-similarity-discovery
