# Contributing

This repo is worked on by multiple people and coding agents in parallel.
Follow this workflow so branches don't collide.

## Branching

- `main` is the shared trunk. Nothing lands on `main` except through a PR.
- One branch per person per task: `<name>/<task-id>-<slug>`, e.g.
  `evan/T2-disease-targets`, `jay/T6-embeddings`.
- Claim a task by editing its row in `docs/02-PIPELINE-PLAN.md` (set
  `Claimed by` + `Status`) — do this in your first commit on the branch so
  the claim itself is the merge point.
- Open a PR as soon as the branch exists (draft is fine). Small PRs, one
  logical change per commit.

## File ownership

Merge conflicts are avoided by construction, not by rebasing:

- Each task owns exactly its entry point `scripts/<task>.py` and its output
  directory `results/<task>/` (or `results/pipeline/<disease>/` — see the
  plan). Don't edit another task's files; if you need a change there, ask
  or leave a TODO.
- Shared files (`data/drugs.csv`, this file, the plan, `CLAUDE.md`) are
  edited only in a PR dedicated to that change.
- New data files go in `data/` with the source and fetch date in a comment
  or docstring at the top of the script that writes them.

## Environment

```bash
python3.14 -m venv env          # no conda, no Homebrew
./env/bin/pip install -r requirements.txt
./env/bin/python scripts/<task>.py
```

API keys go in `.env` (gitignored; see `.env.example`). Never commit keys.

## Rigor rules (from CLAUDE.md — they apply to every task)

- Never report a number that wasn't computed. No placeholder metrics.
- Keep negative results — an underperforming method goes in the table.
- Every headline number carries n; comparisons carry a significance test.
- Anything not evaluated goes in the plan's "not evaluated" table with a
  reason, rather than being silently skipped.
