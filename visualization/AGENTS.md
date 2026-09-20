# Isolated visualization workstream

This subtree is owned by the visualization worker. The scientific pipeline is being developed concurrently by another agent. Do not modify `scripts/`, `data/`, `results/`, `autorepurpose/`, root requirements, root ignore rules, or the shared Git branch/index while doing visualization work. Do not commit or stage another worker's files. On a shared checkout, do not switch branches; use a separate clone/worktree if Git isolation is required.

## Stable handoff to the scientific agent

The integration boundary is a JSON manifest validated against `visualization/schema.json` (version `1.0.0`), plus immutable PDB/mmCIF files referenced by bundle-relative paths and SHA-256 hashes. No scientific code imports or callbacks are required. Never read an in-progress output: the producer should publish the manifest only after all referenced files are complete, preferably by atomic rename. Build a fresh export for each run/revision.

Generate a complete, intentionally synthetic example:

```
visualization/.venv/bin/python -m visualization fixture --out visualization/out/example-source
```

`example-source/manifest.json` demonstrates every required field. It and the toy structure are synthetic UI test data, not scientific output. Do not relabel the example as a real run. Schema field names are scoped to the display boundary, not a demand to rename pipeline models.

Producer checklist:
- Stable disease, target, signature, molecule, and structure IDs.
- Actual signature source: BoltzGen consensus, pocket geometry, known-ligand reference, or hybrid. A crystal-derived signature must not be labeled generated.
- Author-numbered residue mappings explicitly include chain and insertion code. Canonical `key` strings are opaque to the viewer; the same key must denote the same target-side residue in every structure. No integer-only inference or SIFTS guessing.
- One model per structure. PDB/mmCIF parser uses author chain/residue numbering. Alternate conformers follow BioPython's selected (normally highest-occupancy) conformer; resolve disorder upstream when another choice is needed.
- `ligand_residues` must explicitly identify every binder residue. No residue-name guessing, ligand building, docking, contact extraction, or candidate scoring happens here.
- `frame_id` asserts a shared coordinate frame, not merely a shared target. Independently predicted poses must use different frame IDs or provide explicit target-atom alignment pairs and an upstream-approved RMSD gate. Failed alignment disables synchronized comparison; display transforms apply to the entire complex, never the ligand alone.
- `coverage.value` must equal the fraction of the declared core in `engaged`; engaged/missed sets must partition that core. Failed or absent measurements use a null value with status/reason, not zero.
- Ranks and their metric come from the pipeline. Failed confidence gates cannot be ranked. Small-molecule ranks require an assessed addressable site. Rejected/unranked records remain inspectable.
- Confidence values retain native name/scale and an explicit gate rule. A percentile is not a probability.
- Decoy arrays, n, percentile, tie rule, and target/signature/modality context are required together. The viewer checks context and displays supplied calibration; it does not run significance tests or infer missing calibration.
- Supply approval jurisdiction/source as human-readable `approval`; do not imply FDA status from a global database flag. Supply citation identifiers and reviewed claim text; absent citations remain absent.
- Provenance is an allowlisted display object, not a raw config dump. Never include credentials, private paths, API responses with authentication headers, or signed URLs.

The current evolving `scripts/interfaces.py` implementation says insertion codes are dropped and contacts are distance approximations. A producer must resolve those mapping ambiguities before a real export and preserve those scientific caveats; this workstream does not modify that code or certify its outputs.

## Setup and use

All commands below run from the repository root. Existing pipeline `env/` is untouched.

```
python3 -m venv visualization/.venv
visualization/.venv/bin/python -m pip install -r visualization/requirements.txt
visualization/.venv/bin/python -m visualization.setup_assets
visualization/.venv/bin/python -m visualization validate --manifest PATH/manifest.json
visualization/.venv/bin/python -m visualization build --manifest PATH/manifest.json --out visualization/out/RUN_ID
visualization/.venv/bin/python -m visualization serve --bundle-dir visualization/out/RUN_ID --port 8000
```

The portal always loads `web/examples.json`: the worked examples from `SUCCESSES.md`, transcribed by hand and rendered as recorded runs (nothing is recomputed in the browser). A built export's `bundle.json` is optional and, when present, is offered alongside them as the saved analysis. To demo the worked examples without a scientific export, serve the source directory: `visualization/.venv/bin/python -m visualization serve --bundle-dir visualization/web --port 8000`. When `SUCCESSES.md` changes, update `examples.json` to match; `tests/test_examples.mjs` checks every recorded number against the document.

The colorectal worked example also carries the six co-folded complexes from that run in `web/poses/`, described by `web/example-structures.json` and displayed through the existing `MolecularStage`, which is reused unchanged. Contact residues were extracted upstream by `demo/contacts.py` at the same 4.5 A heavy-atom cutoff and are only displayed here. Two numbering and frame limits are recorded in that file and surfaced in the UI rather than worked around: residue numbers in the predicted complexes are the model's own (`author_offset` measured +833 but refused it at 83% agreement, because the reference crystal structure is missing 56 loop residues), and the target-only rigid fit failed the 2.0 A gate for every pose (3.29 to 6.40 A), so synchronized cross-molecule comparison stays disabled and each complex keeps its own frame. Compare contact sets, which are frame-independent, not camera positions.

The browser layer uses native ES modules and pinned local 3Dmol 2.4.2. No Node build, external fonts, CDN, account, live scientific service, or GPU is required. `setup_assets` is the sole optional network download for runtime assets; it checks the pinned npm SHA-512 and retains the BSD license. Vendored assets are copied into each export so viewing is offline.

The interactive viewer requires a local HTTP server, not `file://`. `report.html` is a self-contained static snapshot that works via `file://`. It explicitly says figures were not rendered until actual PyMOL PNGs are supplied. Large generated exports, the venv, and browser binaries are ignored only by this subtree's `.gitignore`.

## PyMOL gate

Actual PyMOL is not installed in the verified shell. A licensed standalone application or an approved external open-source installation is needed. Never install Homebrew/conda or purchase a license silently. Node is not needed. Rendering verification must remain blocked until an actual executable runs; unit tests and generated scripts are not proof that PyMOL rendered successfully.

```
visualization/.venv/bin/python -m visualization render-pymol --bundle visualization/out/RUN_ID/bundle.json --executable /APPROVED/PATH/TO/pymol --out visualization/out/RUN_ID/pymol
visualization/.venv/bin/python -m visualization report --bundle visualization/out/RUN_ID/bundle.json --renders visualization/out/RUN_ID/pymol --out visualization/out/RUN_ID/report-with-figures.html
```

PyMOL runs its own interpreter. It consumes the exact normalized display PDBs, residue atom selections, colors and scene presets used by the browser; it does not import science packages. Scenes and a `.pse` session plus 1920×1080 PNGs are generated. Pixel-identical camera projection between engines, movies, and live rendering are not promised. The baseline figure set is target/site/consensus/coverage/missed residues. The optional `design_reference` object (`id`, `name`, `structure_id`, `protocol`, `source`) adds a representative-design comparison and PyMOL scene. Its structure's `ligand_residues` must identify all binder chains/residues. A representative is explicitly not the consensus, and a known ligand must never be relabeled as a designed binder.

## Verification

```
visualization/.venv/bin/python -m pip install -r visualization/requirements-test.txt
PLAYWRIGHT_BROWSERS_PATH="$PWD/visualization/.browser-cache" visualization/.venv/bin/python -m playwright install chromium
visualization/.venv/bin/python -m unittest discover -s visualization/tests -v
VIZ_BROWSER_TEST=1 PLAYWRIGHT_BROWSERS_PATH="$PWD/visualization/.browser-cache" visualization/.venv/bin/python -m unittest discover -s visualization/tests -v
```

Tests use temporary synthetic fixtures and do not call Rowan or mutate scientific outputs. Browser tests cover chapter navigation, residue/candidate selection, filters, empty modalities, static offline reports, WebGL fallback, mobile layout, and absence of external requests. Geometry tests cover rigid transforms and unsafe/ambiguous input rejection. Run the narrow visualization suite, not the other agent's science suite.

## Explicit release limits

This is the first independent visualization implementation, not a completed scientific run. Real-run readiness depends on a valid upstream manifest. Actual PyMOL rendering, real-run frame rate, Safari/WebKit coverage and external scientific validation must not be reported as verified until run. Unsupported convergence/ablation series should be exported as unavailable validation records rather than fabricated plots. No browser code computes scientific metrics, significance, or therapeutic conclusions.
