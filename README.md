# drug-similarity-discovery

Benchmarking representations for drug similarity matching — and building an
auto-research pipeline that, given a disease, finds a binding target, an ideal
small-molecule binder, and the closest FDA-approved drugs ranked by predicted
affinity drop-off.

Built at HackMIT 2026. MIT license.

## Layout

| Path | Contents |
|---|---|
| `CLAUDE.md` | project rules: rigor standards, deadlines, conventions |
| `CONTRIBUTING.md` | branching + file-ownership workflow for parallel work |
| `docs/` | numbered notes, incl. `02-PIPELINE-PLAN.md` — the executable task plan |
| `scripts/` | runnable entry points; each writes under `results/` |
| `data/` | DrugCentral corpus (`drugs.csv` = 2,114 drugs w/ SMILES + human targets) |
| `results/` | computed outputs only — never hand-edited numbers |

## Quickstart

```bash
python3.14 -m venv env
./env/bin/pip install -r requirements.txt
./env/bin/python scripts/benchmark.py   # target-retrieval benchmark, ~1 min
```

Working on the auto-research pipeline? Read `docs/02-PIPELINE-PLAN.md`,
claim a task, branch as `<name>/<task-id>-<slug>` per `CONTRIBUTING.md`.
