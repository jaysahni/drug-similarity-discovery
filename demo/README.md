# demo — AutoRepurpose, lightweight

The PROJECT_GOAL.md pipeline end to end, in ~900 lines, with the scientific work
done by the [`cheminformatics-kit`](https://github.com/jaysahni/cheminformatics-kit)
toolkit (`novakit`) rather than by this repo's `scripts/` (12,305 lines).

```
disease/target → site → BoltzGen designs → interface signature → approved drugs ranked
```

## Run it

```bash
uv venv --python 3.12 env-kit                       # novakit needs >=3.12,<3.13
uv pip install --python env-kit/bin/python \
  "novakit[rowan,chem] @ git+https://github.com/jaysahni/cheminformatics-kit"
cp .env.example .env                                 # add ROWAN_API_KEY

./env-kit/bin/python -m demo.pipeline --stage all
./env-kit/bin/python -m demo.report
```

Stages also run individually: `--stage research|site|design|signature|match`.

Every Rowan workflow uuid is cached under `results/demo/_cache/`, so a re-run
re-reads the finished workflows and **costs no credits**. Delete a cache entry to
force that stage to run again.

## Files

| File | Lines | Does |
|---|---|---|
| `pipeline.py` | ~480 | the five stages, and the residue-numbering map between them |
| `nova.py` | ~170 | toolkit + Rowan wrapper: `.env` loading, submit/poll, disk cache, credit guard |
| `contacts.py` | ~160 | 4.5 Å heavy-atom contacts, consensus signature, coverage |
| `library.py` | ~100 | the 14-drug shortlist, in three tiers |
| `report.py` | ~180 | renders `results/demo/cdk2/report.md` from the computed JSON |

Outputs land in `results/demo/cdk2/`: one JSON per stage, plus `report.md`.

## Read this first

**`docs/06-LIGHTWEIGHT-DEMO.md`** covers the parts that will bite you:

- why the demo needs its own `env-kit/` virtualenv (novakit is 3.12-only, `env/`
  is 3.14)
- the three incompatible residue numberings that meet in this pipeline
- why the site is chosen by UniProt annotation and not by pocket score — the ATP
  site is Rowan's **4th-ranked** pocket of 5
- what the toolkit did not cover, and what the demo does instead
- the scope cuts: 14 drugs not 4,099, 12 designs not 2,000, one target

## What this demo does not claim

`results/m2_gate_CDK2.json` already measured the BoltzGen consensus against plain
pocket geometry on 31 held-out CDK2 co-crystals, and **BoltzGen lost** (Jaccard
0.338 vs 0.494, Holm p = 0.0016). So stage 5 scores every drug against *both*
signatures and reports both columns. The demo shows the pipeline running; it does
not claim the generated arm is the better one, because this repo has already
measured that it is not.
