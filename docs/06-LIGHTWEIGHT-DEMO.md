# 06 — Lightweight demo on the NovaKit toolkit

A five-stage, runnable version of the PROJECT_GOAL.md pipeline that does the
scientific work through [`cheminformatics-kit`](https://github.com/jaysahni/cheminformatics-kit)
(the `novakit` package) instead of this repo's own `scripts/`.

Target: **CDK2 (P24941)**, design structure **6Q4G** — the same target and
structure as `results/m2_gate_CDK2.json`, so the demo's numbers can be compared
against the heavy pipeline's rather than sitting on their own scale.

```
research → site → design → signature → match
UniProt    Rowan   Rowan     contacts   Rowan cofold
PubMed     pockets BoltzGen  consensus  → ranked approved drugs
```

Entry point: `./env-kit/bin/python -m demo.pipeline --stage all`

## Why a second virtualenv

**`novakit` requires Python `>=3.12,<3.13`. This project's `env/` is 3.14.7.**
They cannot be the same interpreter. The demo therefore runs in `env-kit/`
(Python 3.12.13, created with `uv`), which is gitignored:

```bash
uv venv --python 3.12 env-kit
uv pip install --python env-kit/bin/python \
  "novakit[rowan,chem] @ git+https://github.com/jaysahni/cheminformatics-kit"
```

This is the same class of constraint already recorded in `requirements.txt` for
`boltz` ("requires Python <3.13 and cannot be installed in `env/`"). Nothing in
`demo/` is importable from `scripts/` and nothing in `scripts/` changed.

`novakit` is installed from a pinned commit (`60056db`). The repo is private, so
the install needs the git credentials already configured on this machine.

## Credentials

`ROWAN_API_KEY` in `.env` (gitignored; `.env.example` shows the shape).

`demo/nova.py:load_env` deliberately reads `.env` into `os.environ` **before**
any NovaKit call. NovaKit falls back to the OS keychain when the variable is
absent, and that call blocks on a GUI prompt in a non-interactive shell — it
hung this session once before being killed.

## Cost

Rowan charges per workflow. Measured on this run:

| Stage | Tool | Credits |
|---|---|---|
| site | `rowan.detect_pockets` | 0.05 |
| design | `rowan.design_protein_binder` (12 designs) | see `03_designs.json` |
| match | `rowan.cofold` × 14 drugs | see `05_ranked.json` |

For reference, the heavy pipeline's 24-design CDK2 run cost 55.14 credits and
took 1,103 s on an A100-80GB, and a single KDR co-fold cost 13.09.

`demo/nova.py` enforces a ceiling (`DEMO_BUDGET_CREDITS`, default 400) and passes
an explicit `max_credits` on every billable call — NovaKit refuses a billable
Rowan request that does not. Every workflow uuid is cached under
`results/demo/_cache/`, so **a re-run re-reads finished workflows and costs
nothing**.

### One inefficiency, left in and flagged

PROJECT_GOAL.md §2.2 calls MSA reuse "the single biggest lever on total runtime":
compute the MSA once per target, reuse it across every ligand. The demo does not
do this. `rowan.cofold` exposes `use_msa_server: bool` and no way to hand it a
precomputed alignment, so each of the 14 co-folds recomputes the same CDK2 MSA.

The toolkit has `rowan.generate_msa`, but wiring its output into the co-fold
workflow is not something the pinned SDK's signature supports. At 14 drugs this
costs time rather than correctness; at 4,099 it would be the dominant cost and
would have to be solved first.

## Scope cuts, stated plainly

This is a demo, so scope was cut — not rigor (CLAUDE.md working agreement).

- **14 drugs, not 4,099.** Co-folding the full approved library in
  `data/approved_drugs.csv` at ~13 credits each is ~53,000 credits. The
  shortlist in `demo/library.py` is fixed in source and split into three tiers
  (positive / candidate / decoy) so the leaderboard is interpretable. A
  shortlist chosen by an upstream model would make the result a measure of that
  model.
- **12 designs, not 2,000.** PROJECT_GOAL.md §1.2 budgets 2,000; the M2 gate used
  24. Twelve is enough to build a consensus and is labelled as such. The
  convergence question (task C5) is not answered here.
- **One target, one structure.** No cross-target claims.
- **Small molecules only.** Tiers 2 and 3 of PROJECT_GOAL.md §2.2 (peptides,
  biologics) are not evaluated.

## Site selection, and what it is not

`rowan.detect_pockets` returns 5 pockets for 6Q4G. The demo picks the one
overlapping the **UniProt-curated ATP/Mg binding annotation** for P24941
(`ft_binding`, `ft_act_site` — 20 residues), not the top-scoring one.

That rule matters, and the reason is a finding in itself:

> **The ATP site is Rowan's 4th-ranked pocket of 5** (score 3.00 vs 5.85 for the
> top pocket). The top two pockets overlap the annotated site in **zero**
> residues. Taking the highest-scoring pocket would have designed against the
> wrong site entirely.

The annotation is used rather than `known_ligand_contacts.json` because that file
is marked **validation-use-only** under task B12: it is derived from the ligand
bound in the very structure the signature is later scored against. UniProt's
annotation is curated independently of 6Q4G, so it can drive design without
leaking. The known-ligand contacts are read back in stage 4 **only** to score the
signature after it is built.

## Residue numbering

Three numbering systems meet in this pipeline, and conflating them silently
produces plausible nonsense:

| Where | Numbering |
|---|---|
| UniProt, PDB author records, all reported output | author numbering (1–298) |
| `rowan.detect_pockets` `residue_numbers` | **0-based** over residues present in the file |
| BoltzGen `binding` spec | **1-based** over residues present in the file |

6Q4G chain A has 282 of CDK2's 298 residues, so the offset is not constant — it
steps at every gap. `demo/pipeline.py:_residue_index_map` builds the map from the
downloaded PDB. `scripts/boltzgen_signature.py` documents the same trap, learned
the hard way when author numbering made BoltzGen exit 1.

Co-folded complexes are the exception: they are numbered 1..N over the sequence
we supply, which is the UniProt canonical sequence, so there residue *i* is
author residue *i*.

## What the toolkit did not cover

| Needed | Toolkit status | What the demo does |
|---|---|---|
| interface contact extraction | `protein_interactions.analyze_interface` is `experimental` and needs administrator-provisioned Modal engines | `demo/contacts.py`, a 4.5 Å heavy-atom cutoff — the same definition as `scripts/interfaces.py`, so the numbers stay comparable |
| BoltzGen directly | no `boltz` code in the toolkit | `rowan.design_protein_binder`; the account's `boltzgen_design_limit_100` feature flag identifies the Rowan `protein_binder_design` workflow as BoltzGen-backed, which is also what `scripts/boltzgen_signature.py` assumes |
| disease → target ranking (Open Targets) | not a toolkit domain | out of scope for the demo; the target is fixed in `demo/pipeline.py:TARGET`. The real chain exists in `scripts/disease_targets.py` and `results/pipeline/colorectal-cancer/` |

## Embeddings: why not ESM on one side and fingerprints on the other

The natural request is "embed the top BoltzGen design, embed the FDA drugs, take
the nearest one". `demo/embed.py` does exactly that shape — but not with two
different encoders, and the reason is not model accuracy.

ESM-2 returns 1280 dimensions of learned protein-sequence features. ECFP4 returns
2048 dimensions of hashed chemical substructures. Dimension 7 is "some sequence
motif" on one side and "substructure hash 7" on the other. A cosine between them
is arithmetic over quantities that are not commensurable; it returns a float, and
the float carries no information. A better protein encoder produces a better
protein vector that is still not comparable to a molecule vector. This is the
wall PROJECT_GOAL.md §1.1 describes: *"A de novo miniprotein has no chemical
space in common with an approved drug. Similarity between them is undefined."*

So both sides are embedded on axes that do mean the same thing — **the 298
residues of CDK2**:

| Side | Vector |
|---|---|
| design | engagement frequency across surviving designs (consensus), or binary contact (one design) |
| drug | binary contact, ligand heavy atom within 4.5 Å of the residue, from the co-folded pose |

Cosine is then well-posed: same axes, same units, same length. Output is
`06_embedding_similarity.json`, ranked, with both a `cosine_consensus` and a
`cosine_top_design` column.

### This has published precedent — it is not an improvisation

The interface-overlap formulation is an established analysis, and citing it
properly also supplies the significance test this repo's working agreement
demands:

- **Davis FP & Sali A**, "The Overlap of Small Molecule and Protein Binding Sites
  within Families of Protein Structures," *PLOS Comput Biol* 6(2):e1000668, 2010.
  [doi:10.1371/journal.pcbi.1000668](https://doi.org/10.1371/journal.pcbi.1000668)
  Of 2,619 protein-binding families, 1,028 also bind small molecules and 197 show
  statistically significant overlap (p<0.01) between the protein-binding and
  ligand-binding positions. Their metric is *the fraction of interface residues
  aligned to ligand-binding-site residues*, tested with Fisher's exact test
  against a null of random independent placement. `demo/embed.py:overlap_pvalue`
  implements that null exactly (hypergeometric, via `math.comb`).
- **Rácz, Bajusz & Héberger**, "Life beyond the Tanimoto coefficient: similarity
  measures for interaction fingerprints," *J Cheminform* 10:48, 2018.
  [doi:10.1186/s13321-018-0302-y](https://doi.org/10.1186/s13321-018-0302-y)
  Benchmarks 44 similarity measures on interaction fingerprints across ten
  targets. Six beat Tanimoto; **cosine is not among the recommended measures.**
  So `demo/embed.py` reports `tanimoto_core` beside the cosine columns. Cosine is
  kept because it uses the frequency-weighted consensus vector rather than a
  binary set, but it carries no published endorsement for IFPs and the docs
  should not imply otherwise.

### Not evaluated, and why

| Approach | Why not |
|---|---|
| ESM-C embedding of the design vs fingerprints of the drugs | No protein language model embeds small molecules. ESM-C (300M/960-d, 600M/1152-d, 6B/2560-d) is protein-sequence-only in every variant. Undefined, not inaccurate. |
| ESM-2 (the toolkit's model) on the design | ESM-2's training set **explicitly excluded de novo designs** — 1,027 sequences tagged "artificial sequence" in UniProt, plus 58,462 similar to 81 known de novo designs. A BoltzGen sequence is out of distribution by construction. |
| [ConPLex](https://github.com/samsledje/ConPLex) | A genuinely shared 1024-d space (trained projections from ProtBert + Morgan FP), but trained on a *binding* objective with a triplet loss against DUD-E decoys chosen to be chemically similar yet non-binding — the objective deliberately decorrelates chemical similarity from embedding distance. It also takes the protein as the *target*, which inverts our query. Validated on 5 kinases, CDK2 not among them. Release 0.1.12 (Feb 2024) is a self-described pre-release, classifiers capped at Python 3.11. |
| [DrugCLIP](https://github.com/bowen-gao/DrugCLIP) | Requires a 3D pocket to encode; a de novo miniprotein has none. Checkpoint is Google-Drive-hosted and the README self-describes the code as "a raw version." |
| MolTrans, HyperAttentionDTI | No separable per-entity embedding — they score a drug–protein *pair* through a joint network, so there is nothing to index or retrieve against. |

## The one method that does compare the two directly

Joint protein–ligand embedding models (ConPLex, DrugCLIP and relatives) put a
protein encoder and a molecule encoder into a shared latent space by *training
projections* on known drug–target interactions. That genuinely makes a sequence
vector and a molecule vector comparable.

It is not used here, for two reasons worth stating rather than discovering later:
a BoltzGen design is a de novo sequence with no evolutionary history and is out of
distribution for a model trained on natural proteins; and such a model scores
"does this molecule bind this protein", which treats the design as a *target*.
The question here is the opposite — which approved drug engages the same hotspots
on CDK2 that the design does. Neither the toolkit nor this repo ships one.

## Carrying the M2 gate forward

`results/m2_gate_CDK2.json` measured the BoltzGen consensus against pocket
geometry on 31 held-out CDK2 co-crystals, and **BoltzGen lost**: Jaccard 0.338 vs
0.494, Δ = −0.155, Holm-corrected p = 0.0016, n = 31.

The demo does not pretend otherwise. Stage 5 scores every drug against **both**
signatures — `core_coverage_boltzgen` and `core_coverage_pocket` — and reports
both leaderboards. Ranking drugs by the BoltzGen signature alone would present an
arm this repo has already measured as the weaker one.
