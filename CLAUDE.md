# drug-similarity-discovery

Benchmarking a range of models for **drug similarity matching** — given a drug, find
the most similar drugs, and measure which representation actually does this well.

Built at HackMIT 2026. Submission targets: HackMIT **Healthcare** track, **Regeneron**,
**Voloridge**. Judging is Innovation 30% / Technical Complexity 30% / Impact 30% /
Learning & Collaboration 10%.

## Deadlines

| When | What |
|---|---|
| Midnight **Sat 2026-09-20** | Project must exist in Plume, or no judging table is assigned |
| 11am **Sun 2026-09-21** | Full project details submitted |

## Working agreement

This is a benchmark, so the results are the product. The standard is the one set by
`../PMC8136837` and `../EVEREST`:

- **Never report a number that wasn't computed.** No placeholder metrics, no
  illustrative values in the README, no "expected" results written ahead of the run.
- **Keep the negative results.** A representation that underperforms is a finding and
  belongs in the results table with its number. Do not quietly drop models that
  disappoint.
- Every headline claim carries **n** and, where it's a comparison, a significance test.
- Anything that can't be evaluated goes in a short "not evaluated / why" table rather
  than being silently skipped.
- Verify absence programmatically before asserting it (empty column, missing file,
  no overlap) — don't assume.
- State hard environment or data constraints up front, in writing, before building on
  them.

Under hackathon time pressure, cut **scope**, not rigor: fewer models and one dataset,
each measured properly, beats a wide grid of unvalidated numbers.

## Conventions

- Python venv at `env/` (Python 3.14.7, arm64). No conda, no Homebrew.
  `pip install` works — cp314 arm64 wheels exist for rdkit, numpy, pandas,
  scikit-learn and torch.
- `scripts/` holds runnable entry points; each writes its outputs under `results/`.
- `docs/` holds numbered notes (`00-ENVIRONMENT.md`, `01-DATA.md`, ...) in the style
  of `../PMC8136837/docs/`.
- **MIT license** — Regeneron's challenge prefers an MIT-licensed public repo.
- Commit as you go, one logical change per commit. Nothing is pushed without asking.

## Context

Durable project facts (sponsor requirements, deadlines, toolchain findings) live in
this session's memory directory and are indexed in its `MEMORY.md`. When a decision
or an external requirement is discovered, write it there — an earlier session's
sponsor research was lost to a `/clear` and had to be recovered from a raw transcript.
