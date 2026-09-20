# demo — AutoRepurpose, lightweight

The PROJECT_GOAL.md pipeline end to end, with the scientific work done by
[`cheminformatics-kit`](https://github.com/jaysahni/cheminformatics-kit)
(`novakit`) rather than by this repo's `scripts/` (12,305 lines). Three
Rowan-hosted stages; everything else runs locally on a Mac.

Two questions are answered here, and they are not the same question:

```
FOOTPRINT MATCHING   design a binder -> throw it away, keep the residues it touched
                     -> rank approved drugs by how many of those they also touch

DIRECT MATCHING      generate a candidate -> find the approved drug most similar
                     TO THE CANDIDATE ITSELF
```

## Run it

```bash
uv venv --python 3.12 env-kit                       # novakit needs >=3.12,<3.13
uv pip install --python env-kit/bin/python \
  "novakit[rowan,chem] @ git+https://github.com/jaysahni/cheminformatics-kit" \
  torch transformers biopython scipy
cp .env.example .env                                 # add ROWAN_API_KEY
```

**Footprint arm** (CDK2, 5 stages — research → site → design → signature → match):

```bash
DEMO_TARGET=cdk2 ./env-kit/bin/python -m demo.pipeline --stage all
DEMO_TARGET=cdk2 ./env-kit/bin/python -m demo.report      # -> results/demo/cdk2/report.md
DEMO_TARGET=cdk2 ./env-kit/bin/python -m demo.embed       # interface-vector similarity
```

**Direct-matching arms** (thrombin). Stages 1–3 come from the same pipeline:

```bash
DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.pipeline --stage research
DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.pipeline --stage site
DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.pipeline --stage design

DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.peptide_arm          # arm P
DEMO_TARGET=thrombin ./env-kit/bin/python -m demo.generate --stage all # arm S
```

**Tests** — none need credentials or credits:

```bash
./env-kit/bin/python -m demo.selftest           # contact extraction vs the heavy pipeline
./env-kit/bin/python -m demo.test_match_direct  # similarity + NULL CALIBRATION
./env-kit/bin/python -m demo.test_esmc          # ESM-C embeddings are meaningful
```

`demo.test_match_direct` is the one that gates everything: if the null is not
calibrated, every percentile downstream is wrong. It has already caught one real
bug (see `docs/10`).

## Knobs

| variable | default | |
|---|---|---|
| `DEMO_TARGET` | `cdk2` | `cdk2` or `thrombin` |
| `DEMO_N_DESIGNS` | 20 | BoltzGen designs; `budget` is automatically `n/5` |
| `DEMO_N_GENERATE` | 600 | molecules sampled in arm S |
| `DEMO_N_DOCK` | 24 | how many of them get docked |
| `DEMO_BUDGET_CREDITS` | 400 | hard ceiling per run |
| `DEMO_SITE_RADIUS` | 6.0 | geometric-site radius, Å |

Every Rowan workflow uuid is cached under `results/demo/_cache/`, keyed by
target, so a re-run re-reads finished workflows and **costs no credits**. Delete
a cache entry to force that stage to run again.

## Files

| file | does |
|---|---|
| `pipeline.py` | stages 1–5, target config, and the residue-numbering maps between them |
| `nova.py` | toolkit + Rowan wrapper: `.env` loading, submit/poll, disk cache, credit guard |
| `contacts.py` | 4.5 Å contacts, insertion codes, sequence alignment, consensus signature |
| `match_direct.py` | similarity, the null, novelty classification, three leaderboards |
| `peptide_arm.py` | arm P — ESM-C + BLOSUM62, with positive control and hub check |
| `generate.py` | arm S — SMILES generation, docking, ECFP4 matching |
| `esmc.py` | ESM-C embeddings (MIT, ESMC-300M) |
| `embed.py` | interface-vector similarity for the footprint arm |
| `library.py` | the 14-drug CDK2 shortlist, in three tiers |
| `build_biologics.py` | builds `data/approved_biologics.csv` |
| `report.py` | renders `results/demo/<target>/report.md` |

Outputs land in `results/demo/<target>/`: one JSON per stage.

## Read this before trusting any of it

**`docs/10-DIRECT-MATCHING.md`** — the direct-matching arms, both results, and
the controls that make them readable.
**`docs/06-LIGHTWEIGHT-DEMO.md`** — the footprint arm, the environment, and the
three residue-numbering systems that meet in this pipeline.
**`docs/07-BOLTZGEN-BEHAVIOR.md`** — what ipTM does and does not mean here.
**`docs/08`, `docs/09`** — the biologics corpus and the ESM-C setup.

Four things in this pipeline produced output that looked completely ordinary
while being wrong, and each is now guarded by a test: a ligand left in the pocket
during design, a PDB parser that dropped 28 of thrombin's 259 residues, a
structure preparation step that renames chains, and a null that scored every tie
against the candidate.

## What this demo does not claim

Every result here is negative, and deliberately reported as such. The generated
arm loses to a pocket finder; hotspot coverage does not separate known binders
from decoys; neither the peptide nor the small-molecule arm finds an approved
drug resembling its candidate. In each case a positive control passes, which is
what makes the negative a finding rather than a broken pipeline.

Similarity is not affinity. Nothing here predicts binding.
