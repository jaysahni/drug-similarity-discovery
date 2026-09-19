# AutoRepurpose v0.2 — Interface-Signature Matching for Multi-Modal Drug Repurposing

**Project plan + parallel task decomposition**
Version 0.2 · supersedes v0.1 · Status: pre-implementation · Owner: Jay

---

## 0. TL;DR for agents

The pipeline answers: *given a disease, which already-approved molecule — of any modality —
engages the same target hotspots that a de novo designed binder would?*

```
disease → targets → pocket/epitope → BoltzGen designs (protein/peptide/nanobody)
       → consensus INTERFACE SIGNATURE → co-fold every approved molecule
       → rank by hotspot coverage → report
```

**There is no affinity number anywhere in the ranking path.** Ranking is interface-contact
overlap, gated on structural confidence. Read §1 and §2 before writing code; find your
workstream in §6; work only inside your owned paths.

---

## 1. What changed from v0.1, and why

| v0.1 | v0.2 | Reason |
|---|---|---|
| Generator produces small molecules (FLOWR/DiffSBDD) | Generator is **BoltzGen** producing proteins, peptides, nanobodies | BoltzGen is what Rowan hosts and what you want to use |
| Match generated binder ↔ approved drugs by **chemical similarity** | Match by **target-side interface signature** | A de novo miniprotein has no chemical space in common with an approved drug. Similarity between them is undefined. |
| Rank by predicted affinity Δ vs ideal binder | Rank by **core hotspot coverage**, gated on ipTM | BoltzGen gives interface confidence, not affinity, and no metric reliably converts one to the other |
| Library = ~3k approved small molecules | Library = small molecules **+ peptides + biologics** | Scope expansion; the interface representation is modality-independent so this is now free |

### 1.1 Why interface matching is the right pivot

The binder BoltzGen designs is not the useful artifact — the **hotspot map it reveals** is.
A design's value is the information "a good binder engages target residues {L95, R91, D105,
Y120…} with these interaction types and buries this much surface." That statement is
expressed entirely in **target-side coordinates**, so it is directly comparable across a
small molecule, a cyclic peptide, and a nanobody, none of which resemble each other.

This is the only formulation where the scope expansion and the similarity-only ranking are
compatible. It also gives the pipeline a genuine failure mode: a drug that engages none of
the hotspots scores zero. Chemical-similarity ranking over a ~400-entry biologics library
could not fail, and therefore could not inform.

### 1.2 What BoltzGen actually gives you

Confirmed capabilities, so nobody designs against the wrong assumption:

- **Modalities produced:** proteins, peptides, nanobodies/VHHs, disulfide-bonded peptides.
  Rowan exposes protocols `protein-anything`, `peptide-anything`, `protein-small_molecule`,
  `nanobody-anything`. Note `protein-small_molecule` designs a *protein that binds a small
  molecule* — the reverse of repurposing; we do not use it for target engagement.
- **Design sizes:** validated campaigns used VHHs and proteins of 80–140 residues; peptide
  protocols go down to 8–16-mer cyclic peptides.
- **Pipeline steps (Rowan/local, 5–7 stages):** Design (diffusion) → Inverse folding
  (BoltzIF) → Design folding (Boltz-2 refold to validate the complex) → Folding (standalone
  stability; skipped for peptides/nanobodies) → Analysis (RMSD, SASA, contacts) → Filtering
  and ranking.
- **Scores emitted:** `design_iptm`, `design_ptm`, interface PAE, hydrogen bonds, salt
  bridges, ΔSASA, refold RMSD, combined into a heuristic Quality Score via a worst-case
  weighted metric rank.
- **What it does not emit:** any calibrated affinity. Treat all of the above as *relative*
  ranking signals valid within a target and meaningless across targets.
- **Yield:** brutal. The AWS reference run had 3 of 100 designs pass all filters; the paper
  generated on the order of 10,000–60,000 designs per target to reach ~15 wet-lab-ready
  candidates; the CLI docs suggest `--num_designs 10000` as typical.
  **But we don't need 15 wet-lab candidates — we need a converged contact map**, which
  should stabilize far earlier. Quantifying that is task C5 and it's the main cost lever.

### 1.3 Corrections that carry over from v0.1

- **Circularity via Open Targets.** Target–disease evidence includes ChEMBL known-drug data,
  so rediscovering a drug that already hits the target is a database lookup wearing a lab
  coat. Novelty flags and three leaderboards are mandatory (§6.G).
- **Model circularity, now worse.** Boltz-2 scores BoltzGen designs and the two share
  training data. As one peptide-design paper puts it plainly about this exact setup: high
  scores may reflect "native-like appearance" to a model with shared training priors rather
  than true binding. On top of that, co-folding an approved drug against its known target
  means Boltz-2 has probably seen that complex. **Decoy nulls and the leakage audit are
  not optional garnish; they are what makes any number here interpretable.**
- **Scope honesty.** Output is computationally-ranked, evidence-linked repurposing
  hypotheses requiring experimental validation. Enforced by CI lint on the report template.

### 1.4 The three hard technical problems

These are new to v0.2 and each one can sink the project if handled naively.

**(a) Size asymmetry.** A 400 Da small molecule buries ~300–400 Å² and contacts 10–15
residues. A 110-residue miniprotein buries 1,500–2,500 Å² and contacts 30–50. Raw contact
overlap will rank every biologic above every small molecule, for reasons that have nothing
to do with biology.
→ **Primary score must be coverage of the *core* hotspot subset, not total overlap**, plus
per-modality leaderboards, plus an explicit size-normalized variant. See §4.4 and §6.E.

**(b) Addressability.** A signature derived from protein binders may describe a large flat
PPI epitope that no small molecule can ever drug. Matching small molecules against it is
guaranteed to return junk.
→ Compute an **`sm_addressable` subsignature**: the subset of hotspot residues that form an
enclosed, concave subpocket (cross-referenced against P2Rank/fpocket geometry). If a target
has no addressable subsignature, the report says so and the small-molecule leaderboard is
suppressed rather than populated with noise. *This is a useful output in itself:* "this
hotspot is flat — biologics only."

**(c) The generator may add nothing.** The honest null hypothesis is that P2Rank's top
pocket, computed in under a second, gives you the same hotspot map as 2,000 BoltzGen
designs and many GPU-hours. **If that's true you need to know it early**, so it's a gate
condition at M2, not an afterthought at M4. See §8.3.

---

## 2. Feasibility and compute

### 2.1 Verdict

Feasible. Every component is open-source or API-accessible. The engineering is
straightforward; the risk is concentrated in whether the interface signature carries
information beyond pocket geometry (§1.4c) and in interpretation discipline (§1.3).

### 2.2 Sizing

**Design generation (per target).** The dominant unknown. Budget for 2,000 designs across
3 modalities as the default, with C5 to determine the real convergence point. At roughly
a few GPU-seconds per design for the diffusion step plus refolding, plan on **10–30
GPU-hours per target** until measured. Rowan absorbs this if you don't want to run GPUs.

**Library co-folding (per target).**

| Tier | Contents | Count | Co-fold cost | Est. GPU-hours |
|---|---|---|---|---|
| 1 | Approved small molecules | ~2,500–3,000 | ~20 s each | ~15–17 |
| 2 | Approved peptides | ~130 | ~60 s each | ~2 |
| 3 | Biologics (Fv/scFv domains only) | ~260 | minutes each | ~15–30 |

MSA is computed **once per target** and reused across every ligand — the single biggest
lever on total runtime.

**Tier 3 carries a specific caveat:** on antibody benchmarks Boltz-2 shows an improvement
over Boltz-1 while still lagging behind AlphaFold3. Full IgGs at ~150 kDa are impractical
to co-fold; restrict to Fv/scFv and mark Tier 3 as a stretch goal with an explicit
low-confidence flag. Tiers 1–2 are the reliable core.

**Total: roughly 40–80 GPU-hours per target.** One 8-GPU node for a day. Affordable.

---

## 3. Architecture

```
  disease (EFO)
       │
       ▼
┌──────────────────────────┐
│ WS-A  Target Selection   │  Open Targets + EuropePMC + tractability
│       + novelty flags    │  + modality-intent (is this PPI or pocket?)
└──────────┬───────────────┘
           │ TargetHypothesis[]
           ▼
┌──────────────────────────┐
│ WS-B  Structure & Site   │  PDB/AFDB → strip ligands → P2Rank
│       Prep               │  → SiteSpec + cached MSA + geometry
└──────────┬───────────────┘
           │ SiteSpec
    ┌──────┴────────────────────────────┐
    ▼                                   ▼
┌───────────────────────────┐   ┌──────────────────────────┐
│ WS-C  BoltzGen Design     │   │ WS-D  Approved Library   │
│  protein/peptide/nanobody │   │  Tier1 small molecules   │
│  → filter → CONSENSUS     │   │  Tier2 peptides          │
│     INTERFACE SIGNATURE   │   │  Tier3 biologics (Fv)    │
└──────────┬────────────────┘   └──────────┬───────────────┘
           │ InterfaceSignature            │ Molecule[]
           └──────────────┬────────────────┘
                          ▼
              ┌───────────────────────────┐
              │ WS-F  Co-Folding Service  │  Boltz-2 / Rowan / Chai
              │  every approved molecule  │  → pose + ipTM + contacts
              │  + decoy nulls            │  (NO affinity head used)
              └───────────┬───────────────┘
                          │ CofoldResult[]
                          ▼
              ┌───────────────────────────┐
              │ WS-E  Interface Matching  │  contact extraction,
              │  core coverage, weighted  │  size normalization,
              │  Jaccard, type agreement  │  addressability gating
              └───────────┬───────────────┘
                          │ InterfaceMatch[]
                          ▼
              ┌───────────────────────────┐
              │ WS-G  Rank & Report       │  per-modality boards,
              │                           │  novelty split, evidence
              └───────────────────────────┘

  WS-H Orchestration / contracts / CI      — spans everything, starts first
  WS-I Validation: hotspot recovery,       — spans everything
       P2Rank null, leakage audit
```

**Note the reordering vs v0.1:** co-folding (F) now runs *before* matching (E), because the
match is computed from poses rather than from precomputed embeddings. E is no longer a
prefilter; it's the scoring function.

**Repo layout** — one directory per workstream, exclusive ownership:

```
autorepurpose/
  contracts/      # FROZEN after M0. Owner: WS-H.
  targets/        # WS-A
  structures/     # WS-B
  design/         # WS-C   (BoltzGen + signature construction)
  library/        # WS-D
  matching/       # WS-E   (was retrieval/)
  cofold/         # WS-F   (was affinity/)
  ranking/        # WS-G
  orchestrator/   # WS-H
  validation/     # WS-I
  fixtures/       # golden I/O per stage, committed
  runs/           # gitignored
  configs/
```

---

## 4. Contracts (freeze at M0)

Pydantic v2 in `contracts/`, exported to JSON Schema. Common `Provenance` envelope as in
v0.1 (`run_id`, `stage`, `tool`, `tool_version`, `db_version`, `created_at`, `config_hash`,
`license`).

### 4.1 TargetHypothesis (WS-A → WS-B)

Unchanged from v0.1 except for two added fields:

```python
    site_type_hint: Literal["pocket","ppi_interface","allosteric","unknown"]
    preferred_modalities: list[str]   # ["small_molecule","peptide","biologic"]
```

Everything else carries over: `association_score`, `score_by_datatype`, `tractability`,
`safety_flags`, `mechanism_rationale`, `desired_effect`, `known_drug_evidence`,
`known_approved_drugs`, `literature[]`.

### 4.2 SiteSpec (WS-B → WS-C, WS-F) — renamed from PocketSpec

```python
class SiteSpec(BaseModel):
    site_id: str
    target_id: str
    uniprot_id: str
    structure_source: Literal["pdb","afdb","predicted"]
    structure_id: str
    structure_path: str
    apo_or_holo: Literal["apo","holo","stripped_holo"]
    chain_id: str
    construct_range: tuple[int,int]
    # site definition
    site_kind: Literal["pocket","epitope","interface"]
    residue_ids: list[int]
    center: tuple[float,float,float]
    # geometry — feeds addressability (§1.4b)
    volume_a3: float | None
    enclosure: float | None          # 0-1, buriedness
    concavity: float | None
    max_sphere_radius: float | None
    p2rank_score: float | None
    p2rank_probability: float | None
    # environment flags
    has_metal: bool
    has_cofactor: bool
    has_ordered_water: bool
    # known-ligand reference, when one exists — used by WS-I, NOT by WS-C
    known_ligand_residue_ids: list[int] = []
    msa_cache_key: str | None
    provenance: Provenance
```

### 4.3 InterfaceSignature (WS-C → WS-E) — **the central new object**

```python
class ResidueEngagement(BaseModel):
    residue_id: int
    chain_id: str
    residue_name: str
    frequency: float                  # fraction of surviving designs engaging it
    interaction_types: dict[str,float]  # {"hbond":0.6,"salt_bridge":0.2,"hydrophobic":0.9}
    mean_buried_area: float
    is_core: bool                     # frequency >= core_threshold (default 0.6)

class InterfaceSignature(BaseModel):
    signature_id: str
    target_id: str
    site_id: str
    source: Literal["boltzgen_consensus","p2rank_geometry","known_ligand","hybrid"]
    # provenance of the design ensemble
    modalities: list[str]             # ["protein","peptide","nanobody"]
    n_designs_generated: int
    n_designs_surviving: int
    filter_thresholds: dict[str,float]
    # the signature itself
    residues: list[ResidueEngagement]
    core_residue_ids: list[int]
    mean_buried_area: float
    # addressability (§1.4b) — gates which library tiers get scored
    sm_addressable: bool
    sm_addressable_residue_ids: list[int]
    addressability_rationale: str
    # convergence (§1.2, task C5)
    convergence: dict                 # {"n_at_stable":..., "bootstrap_jaccard_ci":[..]}
    provenance: Provenance
```

`source` is a discriminator, not decoration — WS-I builds `p2rank_geometry` and
`known_ligand` signatures with the *same* model so the ablation in §8.3 is an
apples-to-apples swap.

### 4.4 InterfaceMatch (WS-E output) — **no affinity field exists by design**

```python
class InterfaceMatch(BaseModel):
    mol_id: str
    signature_id: str
    target_id: str
    modality: Literal["small_molecule","peptide","biologic"]
    is_decoy: bool
    # pose + confidence gate
    cofold_backend: str
    pose_path: str
    iptm: float | None                # polymers
    plddt: float | None
    pose_confidence: float            # unified 0-1, backend-documented
    passes_confidence_gate: bool
    # observed contacts
    engaged_residue_ids: list[int]
    engaged_interaction_types: dict[int, list[str]]
    buried_area: float
    n_hbonds: int
    n_salt_bridges: int
    # scores — interface overlap only
    core_coverage: float              # PRIMARY: frac of core hotspots engaged
    weighted_jaccard: float           # frequency-weighted overlap, full signature
    interaction_type_agreement: float
    size_normalized_score: float      # controls for buried-area asymmetry (§1.4a)
    composite: float
    percentile_vs_decoys: float       # calibrated within modality
    low_confidence_flags: list[str]
    provenance: Provenance
```

### 4.5 Molecule (WS-C, WS-D → WS-F)

As v0.1, plus modality support:

```python
    modality: Literal["small_molecule","peptide","biologic","designed"]
    sequence: str | None              # peptides/biologics/designs
    is_cyclic: bool = False
    chains: list[str] = []            # Fv = 2 chains
    smiles: str | None                # small molecules only
```

### 4.6 RankedHit (WS-G)

```python
class RankedHit(BaseModel):
    rank: int
    rank_within_modality: int
    mol_id: str
    name: str | None
    modality: str
    target_id: str
    symbol: str
    disease_efo: str
    core_coverage: float
    composite: float
    percentile_vs_decoys: float
    pose_confidence: float
    engaged_core_residues: list[int]
    missed_core_residues: list[int]   # explains the gap — the "drop-off", restated
    novelty: Literal["known_moa","known_offtarget","novel_pairing"]
    supporting_literature: list[LiteratureEvidence]
    caveats: list[str]
    leakage_flags: list[str]
```

`missed_core_residues` is the v0.2 replacement for v0.1's affinity drop-off: instead of
"2.1 log units worse than ideal," you get "engages 7 of 9 core hotspots, misses R91 and
D105." Same intent, honest units, and actionable.

---

## 5. Rules for parallel execution

Unchanged from v0.1 and still binding:

1. **Exclusive path ownership** — write only your directories + `fixtures/<your-ws>/`.
2. **Contracts freeze at M0**; changes only via batched `contract-change` issues to WS-H.
3. **`--mock` mode ships before the real implementation.** This is what unblocks everyone
   downstream. Schema-valid synthetic output, under a second, zero dependencies.
4. **Golden fixtures committed.** Downstream agents test against your fixtures, not your code.
5. **CLI + Python entry point** for every stage:
   `python -m autorepurpose.<ws> run --config <yaml> --in <path> --out <path> [--mock]`
6. **Content-addressed caching + idempotency.** Mandatory for WS-C and WS-F.
7. **No cross-workstream imports except `contracts`.**
8. **External API calls behind a retry/backoff client with recorded cassettes.**
9. **Tests pass with zero credentials and zero GPU.** Real-backend tests are
   `@pytest.mark.integration`, opt-in.
10. **Cost guard** — WS-C and WS-F refuse batches exceeding `budget_gpu_hours` without
    `--force`.

---

## 6. Workstreams

| WS | Name | Changed from v0.1? | Can start at M0 |
|---|---|---|---|
| A | Target Selection | minor (+2 fields) | yes |
| B | Structure & Site Prep | moderate (geometry for addressability) | yes |
| C | BoltzGen Design + Signature | **rewritten** | yes |
| D | Approved Library | **expanded to 3 tiers** | yes |
| E | Interface Matching | **rewritten** | yes |
| F | Co-Folding Service | **rewritten** (affinity head dropped) | yes |
| G | Rank & Report | moderate | yes (fixtures) |
| H | Orchestration | unchanged | **starts first** |
| I | Validation | **rewritten** | yes |

---

### WS-A — Target Selection

**Owns:** `targets/`
Substantially as v0.1. Tasks A1–A7 carry over unchanged: Open Targets GraphQL client,
disease→EFO resolution, tractability filter, circularity guard with
`--exclude-known-drug-evidence`, EuropePMC literature with PMIDs attached to claims,
mechanism direction classification, mock + fixture.

**New tasks**
- A8. **Site-type hint.** Classify whether the therapeutic hypothesis is a small-molecule
  pocket, a PPI interface, or allosteric. Drives `preferred_modalities` and sets
  expectations for addressability downstream.
- A9. **Relax the small-molecule-only tractability filter.** v0.1 dropped targets not
  small-molecule tractable. Now keep antibody- and PPI-tractable targets and record the
  modality buckets — the whole point of the scope expansion.

**Done when:** ≥10 schema-valid TargetHypothesis records per disease with correct
known-drug flags, and ranking shifts sensibly under `--exclude-known-drug-evidence`.

---

### WS-B — Structure & Site Prep

**Owns:** `structures/`
Tasks B1–B9 from v0.1 carry over: UniProt resolution, structure acquisition (holo > apo >
AFDB), **mandatory ligand stripping before site detection**, P2Rank with the AlphaFold
profile, site triage, construct definition, environment flags, MSA caching, optional
ensembles.

**New / changed tasks**
- B10. **Pocket geometry extraction** — volume, enclosure, concavity, max inscribed sphere.
  These feed the `sm_addressable` determination in WS-C. Without them §1.4b can't be
  implemented. fpocket gives most of this directly.
- B11. **Epitope-mode site definition.** For PPI targets there is no pocket. Define the site
  from interface residues of a known complex, from conservation + surface exposure, or from
  a user-supplied residue list. `site_kind` discriminates.
- B12. **Known-ligand residue extraction** into `known_ligand_residue_ids` — recorded for
  WS-I validation only. **WS-C must never read this field**; assert it in code, because
  leaking it into signature construction invalidates the entire hotspot-recovery benchmark.

**Done when:** for 20 diverse targets, ≥18 produce a SiteSpec with populated geometry, and
MSA reuse is demonstrated across ≥100 co-folding calls.

---

### WS-C — BoltzGen Design + Signature Construction

**Owns:** `design/` — **fully rewritten**
**Goal:** `SiteSpec` → `InterfaceSignature`.

**Tasks**
- C1. **BoltzGen adapter** with two backends behind one interface:
  - `rowan` — the hosted `protein-binder-design` workflow. Rowan encodes your spec into
    BoltzGen's YAML, runs it on GPU, and parses designs back into a structured data model.
    Async submit + poll. Lowest ops burden; start here.
  - `local` — the `boltzgen` CLI (`boltzgen run --steps ...`). Full control, needed for
    C5's convergence sweep and for cost control at volume.
  Both must emit identical `Molecule` records with `modality="designed"`.
- C2. **Multi-modality generation.** Run `protein-anything`, `peptide-anything`, and
  `nanobody-anything` against the same site. Do **not** use `protein-small_molecule` — it
  designs a protein that binds a small molecule, which is the reverse of what we need.
  Record which protocol produced each design.
- C3. **Filtering.** Use BoltzGen's own filter/rank stage, then apply explicit thresholds.
  Sensible defaults from the community: ipTM > 0.85, pTM > 0.8, complex RMSD < 2.5 Å,
  plus refold RMSD < 2.5 Å and balanced amino-acid composition. Make all thresholds
  config-driven — expect to loosen them, since strict settings can leave you with single
  digits out of a hundred.
- C4. **Contact extraction → per-design contact map.** For each surviving design, compute
  the target-side residues engaged and the interaction type per contact (H-bond donor /
  acceptor, salt bridge, hydrophobic, aromatic stacking, cation-π), plus per-residue buried
  area. Use a standard interaction-detection library rather than hand-rolled distance cuts.
- C5. **Convergence study — do this early, it sets your entire compute budget.**
  How many designs until the consensus contact map stabilizes? Generate 4,000 for 3 test
  sites, then bootstrap-subsample at N ∈ {50, 100, 250, 500, 1000, 2000, 4000} and measure
  the weighted-Jaccard stability of the resulting signature. Pick the default N from the
  knee of that curve. **Expect the answer to be far below the 10k–60k used for wet-lab
  candidate selection** — a converged contact map is a much weaker requirement than 15
  synthesizable designs.
- C6. **Signature aggregation.** Per-residue engagement frequency across surviving designs,
  interaction-type distributions, core-residue determination at the configured threshold.
  Report per-modality signatures *and* the union, since peptides and nanobodies will
  emphasize different subregions and that difference is informative.
- C7. **Addressability determination** (§1.4b). Intersect the signature with WS-B's pocket
  geometry: is there an enclosed, concave subregion of the hotspot large enough for a
  drug-like small molecule? Emit `sm_addressable`, the residue subset, and a human-readable
  rationale. **When false, say so loudly** — it's a real finding, not a failure.
- C8. **Baseline signature builders** for the §8.3 ablation, same output type:
  `p2rank_geometry` (top pocket residues, uniform frequency) and `known_ligand` (contacts
  from a real holo complex; WS-I use only).
- C9. `--mock` returning a fixed signature + fixture.

**Done when:** three test sites produce converged signatures with documented N, the
addressability call is correct on a known-flat PPI target and a known-pocket target, and
the P2Rank baseline signature builds through the same code path.

**Pitfalls:** yield is low and the filters do enormous work — the AWS reference run passed
3 of 100. Budget accordingly and cache aggressively. Also: BoltzGen doesn't consistently
demonstrate that designs bind where predicted, so treat the *consensus across many designs*
as the signal and never trust a single design's contact map.

---

### WS-D — Approved Molecule Library

**Owns:** `library/` — expanded to three tiers.

**Tasks**
- D1. **Tier 1 — small molecules.** As v0.1: DrugCentral (~4,100 approved small organic
  molecules) + ChEMBL max_phase=4 cross-check → ~2,500–3,000 after filtering to human,
  systemic, non-salt, drug-like parents.
- D2. **Standardization** (Tier 1): parent extraction, salt stripping, neutralization,
  tautomer canonicalization, stereo handling, `parent_id` mapping.
- D3. **Tier 2 — approved peptides** (~133 in DrugCentral). Sequences, cyclization state,
  non-standard residues, disulfide connectivity. **Non-standard residues are the main
  engineering headache** — decide early whether to model, substitute, or exclude them, and
  record the decision per entry.
- D4. **Tier 3 — biologics** (~260). Extract Fv/scFv domains from mAbs; full IgG co-folding
  is impractical. Flag every Tier 3 entry low-confidence given Boltz-2's antibody-benchmark
  gap vs AlphaFold3. **This tier is a stretch goal — ship Tiers 1–2 first.**
- D5. Annotations across all tiers: known targets (for novelty), approval year, routes, ATC,
  withdrawn status, prodrug flag.
- D6. **Reproducible build** — one command, content-hashed parquet + manifest, byte-identical
  across machines, with count reconciliation against published DrugCentral totals.
- D7. **Decoy sets, per modality.** Property-matched within modality (MW/logP for small
  molecules; length/composition for peptides; domain type for biologics), drawn from the
  library itself rather than generated. WS-F scores these to build the null.

**Done when:** Tiers 1–2 built reproducibly with full annotation coverage and per-modality
decoy sets; Tier 3 scoped with a documented Fv extraction procedure.

---

### WS-F — Co-Folding Service

**Owns:** `cofold/` — rewritten from v0.1's `affinity/`.
**Goal:** (SiteSpec, Molecule) → pose + confidence + contacts. **The affinity head is not
used and its outputs are not stored.**

**Tasks**
- F1. **Backend adapters**, one interface:
  - `rowan` — the co-folding workflow, a unified interface over Chai-1r, Boltz-1, Boltz-2,
    Boltz-2.1, OpenFold-3 and Decaf-Boltz, with pocket and contact constraints settable at
    the Rowan level. Default choice; no GPU ops.
  - `boltz2_local` — open weights, full control, cheapest at volume.
  - `mock` — instant, deterministic.
- F2. **All three modalities through one path.** Small molecule via SMILES; peptide via
  sequence + cyclization; biologic via multi-chain Fv. Verify each modality's YAML encoding
  against fixtures before scaling.
- F3. **Pocket/contact constraints** so you co-fold at the intended site rather than
  wherever the model prefers. Without this, a promiscuous drug will dock somewhere
  irrelevant and score zero for the wrong reason.
- F4. **MSA reuse** from `msa_cache_key`. Mandatory; verify per-ligand wall time drops after
  the first call on a target.
- F5. **Contact extraction from poses** — identical code path to C4. Factor it into a shared
  module under `contracts/` or a tiny `interactions/` lib both import, because the two must
  agree exactly or every match score is garbage. *This is the one permitted shared
  dependency beyond contracts; WS-H owns it.*
- F6. **Decoy scoring** — every target gets its per-modality null from D7.
- F7. Batching, retries, partial-failure isolation, content-addressed caching, resumability.
- F8. Confidence gating: ipTM for polymers, composite confidence for small molecules.
  Records below the gate are retained but flagged, never silently dropped.
- F9. Cost accounting + budget enforcement.
- F10. **Measured throughput report** per modality on your hardware, replacing §2.2's
  estimates with real numbers.

**Done when:** 3,000 small molecules + 130 peptides co-folded against one target with nulls,
inside budget, fully resumable after `kill -9`, with contact extraction byte-identical to
WS-C's.

---

### WS-E — Interface Matching

**Owns:** `matching/` — rewritten from v0.1's `retrieval/`.
**Goal:** (InterfaceSignature, CofoldResult[]) → `InterfaceMatch[]`. This is the scoring
function, not a prefilter.

**Tasks**
- E1. **Core coverage** (primary score): fraction of `core_residue_ids` engaged by the
  molecule's pose. Simple, interpretable, and naturally resistant to the size asymmetry in
  §1.4a because the core set is small.
- E2. **Weighted Jaccard** over the full signature, weighting each residue by its engagement
  frequency. Secondary score.
- E3. **Interaction-type agreement** — engaging R91 via a salt bridge where the designs
  consistently formed a salt bridge is worth more than brushing past it hydrophobically.
- E4. **Size normalization** (§1.4a). Implement and compare at least two approaches:
  normalize by buried area, and rank-normalize within modality. Report both; make the
  default explicit and justified in an ADR. **Per-modality leaderboards are mandatory
  regardless**, because no normalization fully rescues cross-modality comparison.
- E5. **Addressability gating.** If `sm_addressable` is false, suppress the small-molecule
  board and state why. Do not emit a ranked list you know to be noise.
- E6. **Null calibration** — percentile of each score against the per-modality decoy
  distribution for that target. This is what converts a raw overlap number into something
  interpretable.
- E7. **Composite** score with documented, config-driven weights. Resist the urge to tune
  weights against the validation set; if you do tune, hold out a separate set.
- E8. **Explanation output**: engaged core residues, missed core residues, per-residue
  contribution — this is what WS-G renders and what makes a hit reviewable by a human.

**Done when:** scores computed for all three modalities on one target, nulls calibrated,
and the size-normalization comparison documented.

---

### WS-G — Ranking & Reporting

**Owns:** `ranking/`

**Tasks**
- G1. Rank by composite within modality; produce a cross-modality board explicitly labeled
  as indicative only.
- G2. **Novelty classification** — `known_moa` / `known_offtarget` / `novel_pairing` from
  WS-D annotations + WS-A's `known_approved_drugs`.
- G3. **Three leaderboards × per-modality**: all / novel-only / known-only. The known-only
  board is the sanity check and renders every time.
- G4. **`missed_core_residues` as the headline explanatory field** — replaces v0.1's
  affinity drop-off with something in honest units.
- G5. Caveat auto-population from `low_confidence_flags` + `leakage_flags`.
- G6. Evidence linkback: literature claims with PMIDs, the design ensemble the signature
  came from, convergence stats.
- G7. **Report artifact**: one self-contained HTML per run — target summary, site render,
  signature heatmap over the target surface (per-residue engagement frequency is
  genuinely worth visualizing), per-modality boards, null distributions, addressability
  verdict, provenance footer with every tool and DB version. Plus `results.parquet`.
- G8. **Language guard**: header states outputs are computational hypotheses requiring
  experimental validation; CI lints the template for "cures", "proves", "demonstrates
  efficacy", and — new for v0.2 — any use of "affinity", "Kd", "IC50", or "potency" in the
  results section, since none of those are computed anywhere in this pipeline.

**Done when:** one command turns a completed run into a report a collaborator can read
unaided.

---

### WS-H — Orchestration, Contracts, Infrastructure

**Owns:** `contracts/`, `orchestrator/`, `configs/`, CI, and the shared `interactions/` lib.
**Starts first; ship H1–H4 within 48 hours.**

- H1. Write `contracts/` from §4. Pydantic v2 + JSON Schema export.
- H2. Mock factories for every model.
- H3. Scaffold all nine packages with CLI entry point, `--mock`, passing placeholder test.
- H4. Conventions doc: ownership table, contract-change process, fixture conventions.
- H5. **Shared `interactions/` library** — contact and interaction-type detection, imported
  by both WS-C and WS-F (§6.F5). Single source of truth; versioned; changes are
  contract-level events because they silently move every score.
- H6. DAG runner with stage-level resume, retries, structured logging.
- H7. Run manifest: config hash, DB versions, tool versions, licenses, git SHA, hardware,
  cost, wall time. Written incrementally so crashed runs still have one.
- H8. Content-addressed cache shared across stages.
- H9. Budget enforcement with pre-submission estimation.
- H10. CI: full mock-mode pipeline smoke test per PR, contract schema-diff check, report
  language lint.
- H11. Secrets handling for Rowan / hosted backends.

**Done when:** `autorepurpose run --config configs/example.yaml --mock` runs every stage and
produces a report, in CI, under 60 seconds.

---

### WS-I — Validation Harness

**Owns:** `validation/` — rewritten. **This is what makes the project defensible.**

**Tasks**
- I1. **Hotspot recovery — the new headline metric, and the best thing about v0.2.**
  For targets with a known holo complex, does the BoltzGen consensus signature overlap the
  experimentally-observed ligand contact residues? Measure per-residue precision/recall and
  weighted Jaccard against `known_ligand_residue_ids`. **This validates the core premise
  directly and requires no affinity data at all** — it's just geometry against crystal
  structures, and it can run over hundreds of PDB entries. Ship it first.
- I2. **Retrospective drug recovery.** Ground-truth (approved molecule, target, indication)
  triples from DrugCentral/ChEMBL mechanism-of-action data. With the target's known-drug
  evidence excluded from selection, does interface matching rank the known binder in the
  top-k? Report recall@k, enrichment at 1% and 5%, BEDROC (α=20).
- I3. **Temporal holdout.** Old DB snapshot for target selection, newer for ground truth.
  The strongest available guard against the leakage critique.
- I4. **Decoy calibration** per modality, feeding F6.
- I5. **Leakage audit.** Two layers now: (a) was the approved molecule's complex with this
  target plausibly in Boltz-2's training set (PDB deposition date vs training cutoff,
  ligand similarity)? (b) does the BoltzGen design ensemble resemble known binders of this
  target more than chance — i.e. is it recalling rather than designing? Annotate every hit.
- I6. **Ablations.** In priority order:
  1. **BoltzGen signature vs `p2rank_geometry` signature** (§8.3). *The* experiment.
  2. Signature vs `known_ligand` signature — the ceiling.
  3. Per-modality design ensembles vs the union.
  4. Core coverage alone vs weighted Jaccard vs composite.
  5. Size-normalization variants.
  6. Holo vs apo vs AFDB structures.
  7. Design count N — reuses C5's convergence data.
- I7. **Negative controls.** Shuffled disease–target mapping must collapse performance. A
  target with no plausible approved binder must not confidently produce one. A signature
  from one target scored against another target's library must score at chance.

**Done when:** I1 produces numbers, I6.1 has an answer, and the negative controls behave.

---

## 7. Milestones

| ID | Milestone | Gate criteria |
|---|---|---|
| **M0** | Contracts frozen, skeletons + mocks live | All 9 packages run in `--mock`; CI green; shared `interactions/` lib exists. **48h.** |
| **M1** | Component prototypes | A: targets for 1 disease. B: sites + geometry for 20 targets. **C: convergence study done, default N chosen.** D: Tiers 1–2 built. F: 100 molecules co-folded, measured throughput. E: scores on fixtures. I: I1 benchmark set assembled. |
| **M2** | **Premise validated** | **I1 hotspot recovery reported. I6.1 (BoltzGen vs P2Rank) answered.** First end-to-end run on one target, all modalities, with nulls. |
| **M3** | Retrospective validation | I2 recall@k + I3 temporal holdout + remaining ablations. Negative controls pass. Leakage audit clean. |
| **M4** | Hardened | Tier 3 biologics. Multi-disease batch runs. Docs. |

**Gate discipline.** M2 is the real gate and it is deliberately early. Two ways to fail it:

- **I1 fails** — the designed-binder consensus doesn't recover known binding sites. The
  signature is not measuring what you think it is. Stop and diagnose.
- **I6.1 fails** — BoltzGen gives you the same signature as P2Rank's top pocket for a
  fraction of a second of compute. Then the honest pipeline is
  `target → pocket → co-fold library → rank by pocket-residue coverage`, and you should
  ship that: it's faster, cheaper, and more defensible.

**Be genuinely willing to accept either outcome.** Finding this out at M2 costs a few weeks.
Finding it out at M4, after building the whole design stack, costs the project.

---

## 8. Validation design

### 8.1 Primary metric
**Hotspot recovery (I1).** Weighted Jaccard between the BoltzGen consensus signature and the
experimentally-observed ligand contact set, over targets with holo structures. Direct,
geometric, affinity-free, and runnable at scale.

### 8.2 Secondary
Retrospective drug recovery (I2) with recall@k, EF@1%/5%, BEDROC(α=20), under temporal
holdout (I3).

### 8.3 The decisive ablation
`source="boltzgen_consensus"` vs `source="p2rank_geometry"`, everything else held fixed, on
both I1 and I2. Because `InterfaceSignature` has a `source` discriminator and both builders
emit the same type, this is a one-line config swap — which is exactly why the contract was
designed that way. `source="known_ligand"` gives the upper bound.

### 8.4 Calibration
Every reported score carries a percentile against a per-modality, property-matched decoy
null for that specific target. Raw overlap numbers never appear without their percentile.

### 8.5 Success at M3
- Hotspot recovery well above what pocket geometry alone achieves (or a clear, reported
  finding that it isn't).
- Known binders recovered in the top 5% within their modality.
- Temporal holdout doesn't collapse relative to the non-temporal benchmark.
- A handful of `novel_pairing` hits with coherent mechanistic rationale that survive
  expert review.

---

## 9. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| BoltzGen signature ≈ P2Rank pocket; generator adds nothing | **Medium-High** | High (but survivable) | I6.1 at M2; fallback pipeline is simpler and still useful |
| Size asymmetry makes biologics dominate every board | High if unhandled | High | Core coverage as primary (E1), per-modality boards (E4) |
| Signature describes a flat PPI epitope no small molecule can address | Medium-High | Medium | `sm_addressable` gating (C7, E5); report as a finding |
| Model circularity — Boltz-2 scoring BoltzGen designs, shared priors | **High** | High | Decoy nulls (I4), two-layer leakage audit (I5), hotspot recovery against crystals (I1) |
| Contact-extraction drift between WS-C and WS-F | Medium | **Fatal, silent** | Single shared `interactions/` lib owned by WS-H (H5) |
| BoltzGen yield too low / cost blowup | Medium | Medium | Convergence study (C5); Rowan-hosted fallback; budget guards |
| Peptide non-standard residues break co-folding | Medium | Medium | Decide model/substitute/exclude per entry at D3 |
| Tier 3 biologics unreliable (Boltz-2 antibody gap) | High | Low (stretch tier) | Fv-only, low-confidence flags, ship Tiers 1–2 first |
| Apo/collapsed sites degrade everything | Medium | High | Prefer holo (B2), ensembles (B9), flag in output |
| Circular rediscovery presented as novel | High | Fatal to credibility | Novelty flags (G2), three boards (G3), `--exclude-known-drug-evidence` |
| Overclaiming — affinity language creeping into outputs | Medium | High | CI lint bans affinity/Kd/IC50/potency in results (G8) |

---

## 10. Open decisions

- **D1 (M1):** Rowan-hosted vs local BoltzGen as default? Trade ops burden against per-run cost.
- **D2 (M1):** Design count N per target — C5 decides empirically.
- **D3 (M1):** Core-residue frequency threshold (default 0.6). Sensitivity-test it.
- **D4 (M1):** Peptide non-standard residues — model, substitute, or exclude?
- **D5 (M2):** **Does BoltzGen earn its place?** I6.1 decides. Be willing to ship the
  P2Rank-only pipeline.
- **D6 (M2):** Default size-normalization scheme.
- **D7 (M3):** Composite weights — fixed heuristic or learned on a held-out set?
- **D8 (M4):** Ship Tier 3 biologics, or leave the pipeline at small molecules + peptides?

---

## 11. Reference stack

**Targets:** Open Targets Platform (GraphQL, quarterly releases) · EuropePMC
**Structure/site:** UniProt · RCSB PDB · AlphaFold DB · P2Rank (v2.4+, AlphaFold profile) ·
fpocket (geometry) · BioPython / PyMOL
**Design:** BoltzGen (`boltzgen` CLI, or Rowan `protein-binder-design`) · BoltzIF (inverse
folding, internal to the pipeline)
**Co-folding:** Boltz-2 (open weights) · Rowan co-folding workflow (Chai-1r, Boltz-1,
Boltz-2, Boltz-2.1, OpenFold-3, Decaf-Boltz) · AlphaFold3 server as an independent
cross-check on top designs
**Interactions:** PLIP or ProLIF for contact/interaction typing (pick one, wrap it in H5)
**Drug data:** DrugCentral · ChEMBL · FDA Orange Book / UNII
**Infra:** Pydantic v2 · Prefect or Dagster · parquet · pytest + cassettes

---

## 12. Example run config

```yaml
run_id: auto
disease:
  efo_id: EFO_0000311
  name: "cancer"
targets:
  max_targets: 5
  min_association_score: 0.4
  exclude_known_drug_evidence: false
  modalities: [small_molecule, peptide, biologic]
site:
  prefer: [holo, apo, afdb]
  ensemble: 1
design:
  backend: rowan            # or "local"
  protocols: [protein-anything, peptide-anything, nanobody-anything]
  n_designs_per_protocol: 1000     # set from C5 convergence study
  filters:
    iptm_min: 0.85
    ptm_min: 0.80
    refold_rmsd_max: 2.5
signature:
  core_frequency_threshold: 0.6
  source: boltzgen_consensus       # swap to p2rank_geometry for the I6.1 ablation
library:
  tiers: [small_molecule, peptide]  # add "biologic" at M4
  version: "2025.4"
  include_withdrawn: false
cofold:
  backend: rowan
  use_site_constraints: true
  score_decoys: true
  n_decoys_per_modality: 300
matching:
  primary_score: core_coverage
  size_normalization: rank_within_modality
budget:
  gpu_hours: 80
  usd: 400
```

---

*v0.2 replaces chemical similarity with target-side interface overlap and removes affinity
from the ranking path entirely. The consequence worth internalizing: the pipeline's central
claim is now falsifiable by geometry alone (§8.1), and its central risk is that pocket
detection already gives you the answer for free (§8.3). Both are answered at M2. Build
toward that gate.*
