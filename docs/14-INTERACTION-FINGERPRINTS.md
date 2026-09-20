# 14 — Pose-based interaction fingerprints

The small-molecule arm ranks approved drugs by ECFP4 Tanimoto on SMILES. This
note records why that is the wrong measure for a repurposing claim, what was
built to replace it, what was verified about it, and the one thing that stopped
the replacement from being evaluated.

**Headline: the fingerprints are built and validated against published
crystallographic facts. The head-to-head against ECFP4 was NOT run. Nothing in
this note is a result about which representation wins.**

Two findings here stand on their own, independent of that comparison:

1. **Rowan's documented batch-docking helper cannot save poses at all.**
   `rowan.submit_batch_docking_workflow()` does not expose `num_poses_to_save`,
   and `stjames.BatchDockingWorkflow` defaults it to `0`. Every pose from a run
   we had already paid for was scored and discarded, silently. Anyone building
   on Rowan batch docking will hit this. Section: *The toolkit limitation*.
2. **The head-to-head as originally scoped could not have answered the
   question.** With 3 known actives it has 0.23 power against a true ΔAUC of
   0.15 and 0.48 against 0.30. Spending on it would have bought a null that
   meant nothing. Section: *Power comes first*.

## The problem ECFP4 has

`results/demo/thrombin/09_smallmol_match.json` ranks approved drugs against
argatroban. The five nearest are peptidomimetics — angiotensin II, icatibant,
lisinopril — and not one of them is a thrombin drug, although argatroban is one.
ECFP4 counts shared substructure, argatroban reads as a peptide, and so the
neighbourhood it returns is a chemotype, not a target.

A docked pose carries what ECFP4 cannot see: which residues the molecule
actually touches, and how. An interaction fingerprint (IFP) is a fixed-length
bit vector over *site residue × interaction type* — the SIFt construction of
Deng, Chuaqui & Singh (*J Med Chem* 2004, 47:337) with the interaction typing of
IChem/ProLIF. Two molecules that both salt-bridge Asp189 and pack under the
60-loop come out similar even when their scaffolds share nothing.

## What was built — `demo/ifp.py`

Six channels per site residue, each a distinct geometric test rather than one
distance wearing six labels:

| Channel | Criterion | Typing source |
|---|---|---|
| `hydrophobic` | apolar heavy atom ↔ apolar heavy atom ≤ 4.5 Å | ligand: C/S/halogen with no N/O neighbour; protein: per-residue apolar atom table |
| `hbond_donor` (ligand donates) | D···A ≤ 3.5 Å, ≥ 90° at both heavy atoms against their own covalent neighbours | RDKit `BaseFeatures.fdef` Donor, filtered to atoms with H > 0; protein donor table |
| `hbond_acceptor` (ligand accepts) | as above, roles reversed | RDKit Acceptor; protein acceptor table |
| `ionic_cation` (ligand cation) | charged-group centroid ↔ Asp/Glu carboxylate centroid ≤ 4.5 Å | RDKit PosIonizable |
| `ionic_anion` (ligand anion) | centroid ↔ Arg/Lys centroid ≤ 4.5 Å | RDKit NegIonizable |
| `aromatic` | centroid ≤ 5.5 Å and normals ≤ 35° (face-to-face), or centroid ≤ 6.5 Å and normals ≥ 50° (edge-to-face) | RDKit Aromatic rings; protein ring table |

Thresholds are the published IFP defaults (Marcou & Rognan, *JCIM* 2007, 47:195;
ProLIF, Bouysset & Rognan, *J Cheminform* 2021, 13:72). None was fitted here.

Residue labelling reuses `demo/contacts.py`'s insertion-code-aware
`residue_label`, because thrombin is chymotrypsin-numbered and the 60-loop that
forms the S2 lid is *entirely* inside insertion codes 60A–60I. Keying on
`resseq` alone would merge them and delete the most informative part of the site.

### What is typed, and what is only asserted

These are real typed interactions, **not** distance-only contacts. But two
honest limits, both stated in the module docstring:

- **No explicit hydrogens exist anywhere in this pipeline.** 1PPB is a 1.9 Å
  X-ray structure deposited without them and Vina poses come back as heavy
  atoms. The H-bond channels therefore apply a heavy-atom criterion (distance
  plus a ≥ 90° angle at donor and acceptor against their own covalent
  neighbours), not a true D–H···A angle. Donor/acceptor roles are assigned by
  atom typing, not observed.
- **Protonation states are assigned, not measured.** Carboxylates are treated as
  anionic and guanidinium/amidinium/primary amines as cationic. Histidine is
  left neutral, so a His-mediated salt bridge would be missed — at thrombin's
  active site His57 is the catalytic base and its state is not known here.

### Similarity measures

Tanimoto is reported because the rest of the repo uses it, but it discards the
bits that are off in both fingerprints. For an IFP those bits mean something:
unlike a hashed ECFP, every bit names a real place ("residue 189, ionic"), so
"neither molecule touches Asp189" is evidence, not padding.

Rácz, Bajusz & Héberger (*J Cheminform* 2018, 10:48,
doi:10.1186/s13321-018-0302-y) benchmarked 44 measures on binary fingerprints
and found six that outperform Tanimoto. All six are implemented —
Sokal–Michener, Rogers–Tanimoto, Sokal–Sneath 2, Consonni–Todeschini 1 and 2,
Austin–Colwell — and all six use `d`, the both-off count. Cosine is deliberately
absent: it was not among the six. Formulae follow Todeschini et al., *JCIM*
2012, 52:2884, and are pinned in the self-test to values worked out by hand.

## What was verified — `demo/test_ifp.py`

`./env-kit/bin/python -m demo.test_ifp` — no credits, network only for one
optional read-only Rowan GET. All ten checks pass.

The validation target is 1PPB: thrombin with PPACK (D-Phe-Pro-Arg
chloromethylketone) covalently trapped in the active site, solved by Bode et al.
(*EMBO J* 1989, 8:3467). That paper says what PPACK touches, so the fingerprint
has something external to be wrong against. Site axis = the 59 chain-H residues
with an atom within 8 Å of the ligand, × 6 channels = **354 bits, 20 on**.

**Every published contact is reproduced, in the right channel:**

| Channel | Residue | What Bode et al. report |
|---|---|---|
| `ionic_cation` | Asp189 | P1 Arg guanidinium salt bridge at the base of S1 |
| `hbond_donor` | Gly216 | antiparallel β-sheet pairing with the ligand backbone |
| `hbond_acceptor` | Gly193 | carbonyl in the oxyanion hole |
| `aromatic` | Trp215 | D-Phe ring, edge-to-face, in the aryl-binding site |
| `hydrophobic` | Tyr60A, Trp60D | the 60-loop lid over S2 |
| `hydrophobic` | Leu99, Ile174 | floor and wall of the aryl-binding site |

Asp189 is also the *only* ionic bit in the whole fingerprint, which is the
correct answer: S1 has exactly one acidic residue at its base.

The full result: `hydrophobic` {42, 60A, 60D, 98, 99, 174, 191, 192, 195, 213,
215}; `hbond_donor` {189, 214, 216, 219}; `hbond_acceptor` {193, 195, 216};
`ionic_cation` {189}; `aromatic` {215}.

**Meaningfulness, not just execution.** Three further checks:

1. *Identity.* The real pose against itself is exactly 1.000000 on all seven
   measures.

2. *Subsite separation.* PPACK's three moieties were fingerprinted separately.
   They land where the literature puts them — P1 Arg on Asp189/Gly219 and the
   S1 walls (191, 192, 195, 213, 215); P2 Pro on the 60-loop (60A, 60D, 99); P3
   D-Phe in the aryl site (98, 99, 174, 215 plus the Trp215 stack). The
   prediction registered *before* the numbers were seen — that P2 and P3, which
   share the shallow upper cleft, resemble each other more than either resembles
   the buried P1 — holds on **all seven measures** (e.g. Tanimoto 0.125 vs 0.000
   and 0.077; Sokal–Michener 0.9802 vs 0.9689 and 0.9661). Only the P1 moiety
   reaches Asp189.

3. *Decay with displacement.* The real pose was slid out of the pocket in six
   axis directions, 0–10 Å in 1 Å steps (n = 66 fingerprints). Spearman
   ρ = **−0.902, p = 4.6 × 10⁻²⁵**. Zero displacement reproduces the pose
   exactly; at 10 Å similarity has fallen to 0.042–0.167. The assertion is the
   *direction* of the correlation, fixed in advance — no threshold was chosen
   after seeing the value.

Plus: the seven formulae against hand-computed a/b/c/d; symmetry and
self-similarity; a site residue absent from the structure raises rather than
zero-filling a bit; fingerprints on different axes refuse to be compared.

## The toolkit limitation

The test that matters — rank the docked molecules by IFP similarity to
argatroban's pose, rank them again by ECFP4 Tanimoto to argatroban's SMILES, see
which puts the known thrombin binders higher — **could not be run. The poses do
not exist.**

Verified, not assumed (`test_poses_are_missing_upstream`):

- `results/demo/thrombin/08_docked.json` has keys `box, controls,
  credits_charged, n_generated_docked, scores, smiles, workflow_uuid`. Nothing
  three-dimensional: 203 SMILES and 203 Vina scores.
- Rowan workflow `4b5a8570-3086-4f05-a130-8aa56416aa92` reports
  **`num_poses_to_save = 0`**, and `rowan.list_proteins(ancestor_uuid=...)`
  returns **0 descendant structures**. The Vina run generated up to 4 poses per
  ligand (`docking_settings.max_poses = 4`), scored them, and discarded all of
  them.

The cause is a default. `rowan.submit_batch_docking_workflow()` does not expose
`num_poses_to_save`; it is a field of `stjames.BatchDockingWorkflow` that
defaults to `0`, so no call through the documented helper can ever keep a pose.
The workflow has to be built by hand and posted through `rowan.submit_workflow`.
`demo.ifp.redock_with_poses()` does exactly that and refuses to submit without
`approve=True`.

**Cost: the original 203-ligand batch was charged 19.23 credits.** Re-docking is
the only fix. There is no free route: no local docking binary exists in this
environment (`vina`, `smina`, `gnina`, `obabel` all absent, nothing in
`env-kit/`), and a crystal-structure-only test cannot have a negative class
because an untested molecule has no pose in thrombin to fingerprint.

### Status of the re-dock

The submission path was found, built and **validated offline**:
`rowan.submit_workflow` refuses a batch (*"You must provide either
`initial_smiles` or a valid `initial_molecule`"*), so the request goes the way
rowan's own helper sends it — a direct POST of the `stjames` workflow model,
differing only in the one field the helper does not expose.
`stjames.BatchDockingWorkflow(**payload)` builds cleanly and
`num_poses_to_save = 1` survives `model_dump(mode="json")`, which is checked at
submit time so a silent drop cannot repeat the run that kept nothing.

**The run itself has not happened: the submission was refused by this
environment's permission system as a real-world transaction.** Nothing was
posted and **no credits were spent** — every Rowan call made here was a
read-only GET. When permission is granted it is one command:

```
./env-kit/bin/python -m demo.ifp --head-to-head --approve   # ~20.5 credits
```

which submits, waits, downloads the poses, fingerprints them, runs the
comparison on both axes and writes
`results/demo/thrombin/ifp_head_to_head.json`. The workflow uuid is cached under
`results/demo/_cache/dock_thrombin_poses.json` following `demo/nova.py`, so a
re-run is free. `redock_with_poses()` refuses to submit without `approve=True`
and refuses outright if the estimate exceeds the credit ceiling.

One step remains unverified because it cannot be tested without running: how
Rowan hands the saved poses back. `fetch_poses()` reads them as descendant
`Protein` structures via `list_proteins(ancestor_uuid=...)` and
`download_pdb_file`, and raises with what it actually found rather than
returning an empty list if that shape is wrong.

## Power comes first

This project has already been burned once by an underpowered test: a docking
comparison reported as a null at n = 23 (rho −0.084, p = 0.703) came out
significant at n = 199 (rho −0.146, p = 0.040). The effect had been there the
whole time. So the power of this comparison was computed *before* the spend, not
after the result.

The original docked set contains three known thrombin binders (argatroban,
ximelagatran, dabigatran etexilate), and argatroban is the query, so the active
set is 2–3 of 203. A simulated power analysis of the paired bootstrap used in
`demo.ifp.paired_bootstrap_delta_auc` (200 replicates × 600 bootstraps per cell,
200 inactives):

| n actives | power at true ΔAUC = 0.15 | power at ΔAUC = 0.30 |
|---|---|---|
| 3 | 0.23 | 0.48 |
| 5 | 0.26 | 0.54 |
| 10 | 0.25 | 0.77 |
| 20 | 0.39 | 0.94 |

A coin-flip chance of detecting even an enormous effect. "IFP beats chance" is
answerable (the exact one-sided permutation floor at 3 actives of 203 is
7.3 × 10⁻⁷), but "IFP beats ECFP4" is not.

So the plan is **not** to re-dock the same 203 molecules. It is to re-dock with
the active set enlarged, against the same 200 generated decoys, so the
background is unchanged and only the number of actives moves.

## The enlarged docking set — 216 ligands, 14 actives

`demo.ifp.thrombin_ligand_set()`; print it with
`./env-kit/bin/python -m demo.ifp --ligand-set`. Free and offline.

The active set is **not hand-picked, and no SMILES was typed from memory** — a
mistyped SMILES is a different molecule and would not announce itself. Actives
are every molecule in `data/target_annotations.csv` with a quantitative ChEMBL
affinity against **P00734** that clears this repo's own `passes_threshold` flag,
a criterion that existed before this arm did and so cannot have been chosen to
make a result come out. SMILES are joined from `data/approved_drugs.csv` by
`struct_id`. All three original controls keep their exact docked SMILES
(verified equal), so the comparison really is against the same background.

| | n | source |
|---|---|---|
| generated decoys | 200 | `08_docked.json`, unchanged |
| actives (`passes_threshold == 1`) | 14 | ChEMBL pchembl 6.02–8.92 |
| measured but sub-threshold | 2 | apixaban (5.51), edoxaban (5.22) |

Actives: liothyronine, ximelagatran, argatroban, dabigatran etexilate,
betrixaban, quercetin, captopril, sulfaguanidine, hexamidine, sitosterol,
nafamostat, gabexate, camostat, dibrompropamidine. Argatroban is the query, so
**13 actives** are rankable — power 0.77–0.94 at ΔAUC 0.30.

Three honest caveats, all recorded before the run:

- Some of these are almost certainly promiscuous ChEMBL hits rather than
  pharmacological thrombin drugs (liothyronine, sitosterol, quercetin,
  captopril). They are kept because the selection rule is pre-committed;
  dropping them after seeing which ranking they helped would be exactly the
  manipulation this repo's working agreement forbids. The mechanism-of-action
  subset (argatroban, ximelagatran, dabigatran etexilate — bivalirudin excluded
  as a 155-heavy-atom peptide that does not belong in a 30 × 26 × 22 Å box) is
  the n = 3 set, and is reported alongside as the underpowered secondary.
- The 200 generated molecules are **presumed** inactive, not measured inactive.
  They were generated *against thrombin*, so they are enriched for
  thrombin-likeness and are harder decoys than random. Both facts bias toward
  the null, which is the conservative direction.
- apixaban and edoxaban fall in the inactive class by that same pre-committed
  rule. Their ranks are worth reading separately; the labels are not adjusted
  afterwards.

Estimated cost: **20.5 credits** for 216 ligands (the original run was charged
19.23 for 203), inside the 30-credit ceiling.

## The fingerprint axis — and a problem with the pipeline's site

Two axes, both fixed before any pose existed, both reported:

| Axis | Residues | Bits | Contains Asp189? |
|---|---|---|---|
| `box_lining` (**primary**) | 129 | 774 | **yes** |
| `site_02` (sensitivity) | 20 | 120 | **no** |

`box_lining` is every chain-B residue with a heavy atom inside the docking box.
It is set by the docking configuration alone — no activity label touches it —
and it is literally the site the ligands were docked into.

`site_02` is the 20-residue site from `02_site.json`, built upstream as a 6 Å
shell around the catalytic Ser195. **It does not contain Asp189**, the aspartate
at the base of the S1 pocket — the single residue that most distinguishes a
thrombin binder, and the one where PPACK's ionic bit fires in the 1PPB
validation above. Asp189 is position 199 in this structure's numbering; the
offset is not assumed but pinned by the GDSGGP catalytic motif at positions
203–208, and both facts are asserted in `test_axes_are_declared_and_differ`.

An axis that omits S1 is handicapped for this question, which is why it is the
sensitivity check and not the headline. The pipeline site is a strict subset of
the box-lining set, so the two are nested and directly comparable.

## Not evaluated, and why

| Claim | Status | Why |
|---|---|---|
| IFP ranks known thrombin binders above ECFP4 | **not evaluated** | re-dock built and validated offline, blocked at submission by the permission system; ~20.5 credits |
| Which of the seven similarity measures is best for IFPs | **not evaluated** | same blocker; all seven are implemented and formula-checked, and per-measure AUCs are computed by the same run |
| Whether the `site_02` axis' missing Asp189 changes the answer | **not evaluated** | same blocker; both axes are wired and reported by the same run |
| IFP fixes argatroban's peptidomimetic neighbours in the approved-drug corpus | **not evaluated** | the 2,153-drug corpus has never been docked at all, only the 203 generated molecules |
| H-bond geometry with real hydrogens | **not evaluated** | no structure in this pipeline has explicit H; heavy-atom criterion used and labelled as such |
| Histidine-mediated ionic contacts | **not implemented** | protonation state unknown; a neutral His is the conservative choice |
| Water-mediated contacts | **not implemented** | Vina poses carry no waters |

## Files

- `demo/ifp.py` — fingerprints, seven similarity measures, the head-to-head
  driver and its statistics, the pose-availability check, the guarded re-dock.
- `demo/test_ifp.py` — the ten checks above (all passing).
- `./env-kit/bin/python -m demo.ifp --check-poses` — free, read-only; prints the
  pose status and the exact fix.
- `./env-kit/bin/python -m demo.ifp --ligand-set` — free, offline; prints the
  216-ligand set with each active's pchembl and ChEMBL id.
- `./env-kit/bin/python -m demo.ifp --head-to-head --approve` — spends ~20.5
  credits; the whole comparison end to end.
