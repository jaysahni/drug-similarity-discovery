# 02 - Plan: auto-research pipeline for drug repurposing

Status: **ready for task claims**. Feasibility verified 2026-09-19.

## Goal

Given a disease, produce a ranked list of FDA-approved small molecules that
are predicted to bind a disease-relevant target, quantified as the affinity
drop-off vs. an "ideal" binder.

```
disease ──(T2 Open Targets + T3 literature)──▶ ranked targets
target  ──(T5 structure prep)────────────────▶ protein input spec
target  ──(T4 known ligands)─────────────────▶ "ideal binder" (best literature ligand)
drugs   ──(T1 corpus + T6 embeddings)────────▶ embedded approved-drug corpus
ideal binder ──(T7 embedding match)──────────▶ top-k approved candidates
target + candidates ──(T8 Boltz-2 affinity)──▶ pIC50 table, Δ vs ideal
all ──(T9 report + T10 rediscovery check)────▶ REPORT.md
```

## Feasibility findings (verified 2026-09-19)

| Component | Option | Verdict | Notes |
|---|---|---|---|
| Disease → target | Open Targets GraphQL `api.platform.opentargets.org/api/v4/graphql` | ✅ | Free, no key. `disease(efoId).associatedTargets` gives overall + per-datatype scores incl. **Literature** — directly satisfies the "proven literature" requirement. Also `knownDrugs` for validation. |
| Literature evidence | Europe PMC REST / PubMed E-utilities | ✅ | Free, no key. Pull PMIDs + titles per target. |
| Target structure | UniProt → AlphaFold DB / PDB | ✅ | AlphaFold DB serves predicted structures by UniProt ID; PDB for co-crystal pocket definition. |
| **"Ideal binder"** | ~~BoltzGen / Rowan binder design~~ | ⚠️ reframed | BoltzGen and Rowan's `protein_binder_design` workflow design **proteins/peptides/nanobodies** — not small molecules. (`protein-small_molecule` designs a protein to *bind* a small molecule — the inverse of what we need.) Embedding-matching a protein binder against small-molecule drugs is meaningless. **Plan: the ideal binder is the most potent literature ligand for the target** (ChEMBL/BindingDB bioactivity). This is stronger anyway — it's literature-proven by construction. De novo small-molecule generation (DiffSBDD/TargetDiff) needs a GPU we don't have — stretch only. |
| Affinity prediction | Boltz-2 via **Rowan `submit_protein_cofolding_workflow`** | ✅ primary | `rowan-python` (PyPI, py≥3.12 → works in env/). Inputs: protein sequence + SMILES → `pred_value` (pIC50) + `probability_binary`. Runs on Rowan's GPUs; ~3–7 credits/min. |
| Affinity fallback A | Rowan docking (`submit_docking_workflow`, Vina/Gnina) | ✅ | Cheaper (CPU 1 cr/min), needs a prepared protein + pocket box. Use to pre-filter top-k before spending GPU credits. |
| Affinity fallback B | local `boltz` in a separate py3.12 env | ⚠️ | boltz requires py<3.13 (env/ is 3.14) and is CPU-only here — one cofold ≈ 10–60 min. Only if Rowan unavailable. |
| Embeddings | ChemBERTa / MolFormer via `transformers`+torch | ✅ | torch cp314 wheel exists. Morgan/Tanimoto is the already-benchmarked baseline (p@1 0.66 on shared-target retrieval, n=2069). |
| Approved-drug corpus | DrugCentral (already in `data/`) | ✅ | Needs rebuild: `drugs.csv` drops drugs with no annotated target; matching needs all ~4.6k structures. |

## Hard constraints (state before building on them)

1. `env/` is Python 3.14; `boltz` cannot be installed there. All Boltz-2 use
   goes through Rowan, or a separate py3.12 venv (task T11).
2. No GPU. Anything diffusion/generative is hosted-only.
3. Rowan free tier = 500 credits + 20/week. Budget the demo run to
   **≤ ~300 credits**: dock-first filter, Boltz-2 only on ideal binder +
   top ~10 approved drugs. Record credits spent per run in the output JSON.
4. Deadlines: repo in Plume by Sat 2026-09-20 midnight; full submission
   Sun 2026-09-21 11am.
5. ChEMBL REST was flaky earlier — every ChEMBL task needs a BindingDB or
   DrugCentral fallback path.

## How to claim and run tasks

Per `CONTRIBUTING.md`: branch `<name>/<task-id>-<slug>` off `main`, edit your
task's `Claimed by`/`Status` in the table below in your first commit, open a
draft PR. You own only your task's files. Output convention:
`results/pipeline/<disease_slug>/<artifact>` so two agents can run the whole
pipeline on different diseases without colliding.

Suggested demo disease: **angiogenesis-dependent solid tumors**
(e.g. colorectal cancer, MONDO_0005575-ish — resolve exact EFO in T2).
Expected target: KDR/VEGFR2. This makes T10 validation meaningful because
approved VEGFR inhibitors exist (sorafenib, sunitinib, pazopanib, …).

## Task board

| ID | Task | Depends on | Owns | Status | Claimed by |
|---|---|---|---|---|---|
| T1 | Approved-drug corpus | — | `scripts/build_drug_corpus.py`, `data/approved_drugs.csv` | open | |
| T2 | Disease → ranked targets | — | `scripts/disease_targets.py`, `results/pipeline/*/targets.json` | open | |
| T3 | Literature evidence per target | T2 | `scripts/target_evidence.py`, `results/pipeline/*/evidence.md` | open | |
| T4 | Ideal binder from known ligands | T2 | `scripts/fetch_ligands.py`, `results/pipeline/*/ideal_binder.csv` | open | |
| T5 | Target structure prep | T2 | `scripts/prep_target.py`, `results/pipeline/*/target/` | open | |
| T6 | Drug embeddings | T1 | `scripts/embed_drugs.py`, `results/pipeline/embeddings/` | open | |
| T7 | Embedding match → top-k | T4, T6 | `scripts/match_candidates.py`, `results/pipeline/*/candidates.csv` | open | |
| T8 | Affinity prediction + drop-off | T5, T7 | `scripts/predict_affinity.py`, `results/pipeline/*/affinities.csv` | open | |
| T9 | Final report | T8 | `scripts/pipeline_report.py`, `results/pipeline/*/REPORT.md` | open | |
| T10 | Validation: rediscovery check | T8 | `scripts/validate_rediscovery.py`, `results/pipeline/*/validation.json` | open | |
| T11 | Stretch: local boltz py3.12 env | — | `scripts/setup_boltz_env.sh`, `docs/03-BOLTZ-LOCAL.md` | open | |
| T12 | Stretch: Rowan peptide-binder demo | T5 | `scripts/binder_design_demo.py`, `results/pipeline/*/binders/` | open | |

## Task specs

### T1 — Approved-drug corpus
Build `data/approved_drugs.csv`: every DrugCentral structure with a valid
SMILES + approval status + name — **not** filtered by target annotation.
Check the structures TSV for an approval/max-phase column; if absent, join
against DrugCentral's `approval` table or mark source.
**Acceptance:** ≥4,000 rows; ≥97% SMILES parse under RDKit (report rate);
columns `struct_id,name,smiles,approved`.

### T2 — Disease → ranked targets
`scripts/disease_targets.py --disease "colorectal cancer"` → resolves EFO via
the Open Targets search endpoint, pulls `associatedTargets` (overall score +
`datatypeScores`), writes `targets.json` with top 25: ensembl id, symbol,
scores, uniprot id. Flag druggable/novel.
**Acceptance:** runs end-to-end; ≥1 top-10 target has literature-datatype
score >0 and ≥1 known drug; JSON records the EFO id used.

### T3 — Literature evidence
For each of T2's top targets, query Europe PMC for `<target> <disease>
inhibitor OR binder` papers; write `evidence.md` — per target: 1–3 PMIDs,
title, year, one-line claim. A human picks the target from this + T2 scores.
**Acceptance:** every shortlisted target has ≥1 real PMID (verify each PMID
resolves — no hallucinated citations).

### T4 — Ideal binder from known ligands
For the chosen target: pull ChEMBL bioactivities (Ki/IC50/EC50 ≤ 100 nM,
human); fallback BindingDB. Rank by potency, prefer drug-like (MW<600).
Output `ideal_binder.csv`: chembl_id, smiles, best_pIC50, assay refs.
**Acceptance:** ≥1 ligand with documented ≤10 nM potency + assay ID; if the
target has no potent ligand that is itself a finding — write it in the file.

### T5 — Target structure prep
UniProt canonical sequence for the chosen target; locate a PDB co-crystal to
define the pocket (center+size for docking) and the residues for a Boltz-2
pocket constraint. Write `target/sequence.fasta`, `target/pocket.json`,
`target/rowan_cofold_input.yaml` template.
**Acceptance:** pocket json exists and cites a real PDB id; sequence length
< ~1500 aa (truncate to kinase/binding domain if needed and document why).

### T6 — Drug embeddings
`scripts/embed_drugs.py` reads `data/approved_drugs.csv`, writes
`results/pipeline/embeddings/{morgan,chemberta}.{npy,csv}` aligned to
struct_id. Morgan reuses the generator from `benchmark.py`; ChemBERTa =
`DeepChem/ChemBERTa-77M-MLM` mean-pooled CLS. Install torch+transformers
(record versions in output meta). Keep the `list[smiles] -> ndarray`
interface so embeddings plug back into `benchmark.py`.
**Acceptance:** shape (n_drugs, d) for both; zero-vector rows counted and
reported, not dropped silently.

### T7 — Embedding match → top-k
Embed the ideal binder SMILES with each T6 representation; cosine (Tanimoto
for morgan) against the corpus; write `candidates.csv` — per representation,
top 25 with struct_id, name, similarity.
**Acceptance:** both representations produce a ranked list; report the
overlap between lists (low overlap = finding, not a bug).

### T8 — Affinity prediction + drop-off
Rowan `submit_protein_cofolding_workflow` (Boltz-2, `ligand_binding_affinity_index`)
for: ideal binder + union of T7 top-10 per representation. Parse `pred_value`
(pIC50) + `probability_binary` → `affinities.csv` with ΔpIC50 vs ideal.
Optional pre-filter: T5 pocket + Vina docking to cut GPU credits.
**Acceptance:** every submitted ligand has a pIC50 or a logged failure row;
`credits_used` recorded; **no silent skips**.

### T9 — Final report
`results/pipeline/<disease>/REPORT.md`: disease → chosen target (why, PMIDs)
→ ideal binder (potency, source) → ranked approved drugs by ΔpIC50 → top
pick + its known indications. State n and every failure.
**Acceptance:** every number in the report traces to a CSV/JSON in
`results/`; includes the "not evaluated" table.

### T10 — Validation: rediscovery check
For the chosen target, list its *approved* inhibitors (DrugCentral/ChEMBL
max_phase=4, e.g. KDR → sorafenib…). Check: do they appear in T8's ranked
output, and where? Write `validation.json`: hits found, ranks, misses with
hypotheses.
**Acceptance:** this is the pipeline's honest signal that it works; a miss
is reported with a reason, not hidden.

### T11 — Stretch: local boltz
Separate `env-boltz/` on python.org 3.12 (no brew/conda), `pip install
"boltz"` (no `[cuda]`), verify one `boltz predict` affinity run on CPU and
time it. Output: `docs/03-BOLTZ-LOCAL.md` with timing + verdict.

### T12 — Stretch: peptide-binder demo
Rowan `submit_protein_binder_design_workflow` (BoltzGen) on the T5 target —
designs a *peptide* binder as a demo artifact, clearly labeled "not the
small-molecule path". Only if credits remain after T8.

## Risks

- **Credits**: worst case T8 ≈ 11 cofolds × ~5 GPU-min × 3–7 cr ≈ 165–385 cr.
  Mitigation: docking pre-filter, shrink k, T11 fallback.
- **Embedding mismatch**: a potent novel ligand may sit far from all approved
  drugs → few/no close candidates. That's a reportable finding; T10
  quantifies it.
- **Ideal-binder caveat**: best-known-ligand ≠ theoretical ideal. Report
  drop-off as "vs. best literature ligand", not absolute.

## Not evaluated / why (running list)

| Item | Why |
|---|---|
| De novo small-molecule generation (DiffSBDD etc.) | no GPU; BoltzGen designs proteins, not small molecules |
| Boltz-2 local inference | py3.14 env incompatible; CPU-only — deferred to T11 |
| FEP-accuracy validation of Boltz-2 | out of hackathon scope; we take the published benchmark |
