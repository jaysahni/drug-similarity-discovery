# 03 - Scope and constraints

Written **before** building on any of it, per `CLAUDE.md`. Every row below was
checked by running the command, not by assuming. Verified 2026-09-19.

## What this machine can and cannot do

| Capability | Status | How it was checked |
|---|---|---|
| Open Targets GraphQL | available | `POST api.platform.opentargets.org/api/v4/graphql` returned `KDR` for `ENSG00000128052` |
| RCSB PDB (REST + file download) | available | entry metadata for `4ASD`; `4ASD.cif` downloaded, 357,753 bytes |
| EuropePMC REST | available | search returned a non-zero `hitCount` |
| ChEMBL REST | available | `molecule/CHEMBL25.json` returned aspirin's record |
| UniProt REST | available | `P35968` sequence returned |
| RDKit / numpy / pandas / scipy / sklearn | available | `env/` venv, Python 3.14.7 |
| torch 2.14.0 + transformers 5.17.0 | available | installed into `env/`, CPU only |
| **Java** | **installed by us** | macOS ships a stub only (`Unable to locate a Java Runtime`). Temurin JDK 21.0.12.1 aarch64 tarball extracted under the session scratchpad — no Homebrew, no admin, nothing installed system-wide |
| **P2Rank 2.5** | **available** | runs on the local JDK; `prank` prints its version banner |
| **CUDA GPU** | **absent** | Apple silicon, no NVIDIA device |
| **Rowan API key** | **absent** | no `.env` in the repo; `.env.example` is a template only |
| `boltz` (Boltz-1/2) local | **not installable** | requires Python `>=3.10,<3.13`; `env/` is 3.14 (recorded in `00-ENVIRONMENT.md`) |

Two of those are load-bearing and they are the reason the plan below departs from
`PROJECT_GOAL.md`: **there is no GPU and no hosted-inference credential.**

## Consequence for `PROJECT_GOAL.md` (AutoRepurpose v0.2)

v0.2's pipeline is `disease -> target -> pocket -> BoltzGen designs -> consensus
interface signature -> co-fold every approved molecule -> rank by hotspot coverage`.
Two of those stages are neural-structure inference on GPU:

| v0.2 workstream | Reachable here? | Why |
|---|---|---|
| WS-A Target selection | yes | Open Targets + EuropePMC, free APIs |
| WS-B Structure & site prep | yes | RCSB + P2Rank 2.5 + fpocket-style geometry |
| WS-C BoltzGen design | **no** | BoltzGen is GPU diffusion; Rowan is the hosted route and there is no key |
| WS-C signature construction | yes, from a different source | the aggregation maths is source-agnostic by design (`source` is a discriminator in §4.3) |
| WS-D Approved library Tier 1 | yes | DrugCentral is already on disk |
| WS-D Tiers 2-3 (peptides, biologics) | no | only meaningful with a co-folding backend |
| WS-E Interface matching | partly | the scoring functions are pure geometry; they need poses to score |
| WS-F Co-folding | **no** | Boltz-2 needs a GPU or Rowan |
| WS-G Rank & report | yes | |
| WS-I Validation / I1 hotspot recovery | **yes** | §8.1 says it plainly: "just geometry against crystal structures" |

The blocked stages are exactly the ones that produce *poses for arbitrary
molecules*. Without poses we cannot score an approved drug against a target it has
never been crystallised with, so the v0.2 ranking path cannot run at scale. That is
a hard stop, not a scoping preference.

## What we build instead, and why it is the same question

`PROJECT_GOAL.md` §7 anticipates this. Its M2 gate is `I6.1`: does the BoltzGen
consensus signature beat the P2Rank pocket signature? And it instructs, if the
generator adds nothing, to ship the pocket-based pipeline because it is "faster,
cheaper, and more defensible". We cannot run the BoltzGen arm at all — so we run
the arms we can, and we say which arm is missing rather than implying it was tested.

The substitution that makes this work: **a consensus interface signature does not
have to come from generated designs.** §4.3 makes `source` a discriminator precisely
so that different signature builders are swappable. Where v0.2 aggregates contacts
over 2,000 BoltzGen designs, we aggregate contacts over the *real co-crystal
structures a target already has in the PDB*. The object is identical; the evidence
behind it is experimental rather than generated.

Three experiments, all CPU-feasible, all with computed numbers:

**E1 - Ligand-side retrieval (the repo's stated product).** Given a drug, rank all
others; a hit shares a protein target. A panel of representations across families
(fingerprints, pharmacophore, 3D shape, physicochemical, scaffold, learned
chemical-language embeddings) against a random floor and a perfect-ranking ceiling,
with bootstrap CIs and paired significance tests. This already has a validated
baseline at `n=2069` and it is what `CLAUDE.md` says the product is.

**E2 - Target-side hotspot recovery (the v0.2 thesis, GPU-free).** For targets with
several holo complexes: build the consensus signature from the other ligands and
predict a held-out ligand's contact residues. Compare against the P2Rank geometry
signature computed on the ligand-stripped structure. This is §8.3's decisive
ablation with the generated arm replaced by an experimental one, and it answers
whether a consensus of observed binding events carries information beyond pocket
geometry.

**E3 - The bridge.** Do drugs that share a target but are chemically dissimilar
engage the same residues? E1 measures when chemical similarity finds a target-mate;
E3 asks what the misses have in common. This is the join between the two halves.

## Not evaluated / why

| Item | Why |
|---|---|
| BoltzGen de novo binder design | GPU diffusion; no GPU, no Rowan key. The v0.2 `boltzgen_consensus` signature arm is therefore **absent, not tested** |
| Boltz-2 co-folding of the approved library | same; `boltz` also cannot install on Python 3.14 |
| Predicted binding affinity (pIC50, Kd, IC50) | no affinity model is run anywhere, so no affinity number is reported anywhere. v0.2 §G8 bans the vocabulary outright and we keep that ban |
| Approved peptide + biologic tiers | only scoreable with a co-folding backend |
| Docking as a co-folding substitute | Checked, not assumed: `vina` 1.2.7 is on PyPI but publishes no cp314 wheel, and its sdist fails to build on Python 3.14 (setuptools `sdist._add_defaults_ext` path). A separate 3.12 venv could host it; that was judged a worse use of the remaining time than E1-E3, so it is recorded here rather than silently skipped |
