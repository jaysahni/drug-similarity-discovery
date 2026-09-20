"""AutoRepurpose, lightweight: disease/target -> site -> BoltzGen -> signature -> approved drugs.

A five-stage demo of PROJECT_GOAL.md's pipeline, built on the NovaKit toolkit
(github.com/jaysahni/cheminformatics-kit) instead of this repo's own 12,000 lines
of scripts/. Every scientific step is a toolkit call; this file is orchestration.

  1 research   UniProt + PubMed           -> the target and why
  2 site       Rowan pocket detection     -> which pocket, in author numbering
  3 design     Rowan BoltzGen             -> de novo binders against that pocket
  4 signature  contacts over the designs  -> consensus interface signature
  5 match      Rowan cofold, per drug     -> approved drugs ranked by hotspot coverage

Run:  ./env-kit/bin/python -m demo.pipeline --stage all
Stages are cached by Rowan workflow uuid, so a re-run costs no credits.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from demo import contacts, library, nova

REPO = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# target configuration
# ---------------------------------------------------------------------------

TARGETS = {
    # Same design structure the heavy pipeline used (results/m2_gate_CDK2.json),
    # so the demo's signature is comparable to the one it reports.
    "cdk2": {
        "symbol": "CDK2",
        "uniprot": "P24941",
        "pdb_id": "6Q4G",
        "chain": "A",
        "out": "cdk2",
    },
    # Thrombin, for the direct-matching arms. Chosen because it is the one
    # target where BOTH arms have a real approved positive control: bivalirudin
    # (20 aa, data/approved_biologics.csv) for the peptide arm, and argatroban /
    # ximelagatran / dabigatran etexilate (moa_targets = F2) for the small
    # molecule arm.
    #
    # 1PPB is the canonical Bode structure: chains H (catalytic heavy chain, 259
    # residues = UniProt 364-622) and L (light chain), with PPACK in the active
    # site as a HETATM that stage 2 strips. Structures like 4UD9/4UE7/5AFY were
    # rejected: they carry a hirudin fragment as a POLYMER chain, which
    # remove_heterogens does not remove, so designing against them would repeat
    # the occupied-pocket mistake with a peptide instead of a ligand.
    "thrombin": {
        "symbol": "F2",
        "uniprot": "P00734",
        "pdb_id": "1PPB",
        "chain": "H",
        "out": "thrombin",
        # Anchor the geometric site on the catalytic nucleophile Ser195
        # (UniProt 568) rather than the whole triad: His57 and Asp102 sit BEHIND
        # the serine, so a sphere around all three reaches into the protein core
        # and away from where a ligand binds. Chosen a priori from enzymology,
        # not tuned against the co-crystal ligand.
        "site_anchors": [568],
    },
}


def _select_target() -> dict:
    name = os.environ.get("DEMO_TARGET", "cdk2").lower()
    if name not in TARGETS:
        raise SystemExit(f"DEMO_TARGET={name!r} unknown; choose from {sorted(TARGETS)}")
    return TARGETS[name]


TARGET = _select_target()
OUT = REPO / "results" / "demo" / TARGET["out"]
STRUCTURES = OUT / "structures"

N_DESIGNS = 12          # heavy run used 24 at 55.14 credits; half that here
DESIGN_BUDGET = 12
CORE_THRESHOLD = 0.6    # PROJECT_GOAL.md 4.3 default
COFOLD_MAX_CREDITS = 25  # per drug; observed ~13 on a ~300-residue kinase


def write(name: str, payload: dict) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    path.write_text(json.dumps(payload, indent=1, default=str))
    print(f"  -> {path.relative_to(REPO)}")
    return payload


def read(name: str) -> dict:
    path = OUT / name
    if not path.exists():
        raise SystemExit(f"{path.relative_to(REPO)} missing -- run the earlier stage first.")
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# stage 1 - research
# ---------------------------------------------------------------------------


def stage_research() -> dict:
    """Identify the target and pull its sequence and literature support.

    NovaKit tools: database.fetch_uniprot_entry, literature.search_pubmed.
    """
    print("[1/5] research")
    nova.load_env()

    entry = nova.tool_data(
        "database.fetch_uniprot_entry",
        {"accession": TARGET["uniprot"], "live_network": True},
    )

    query = f'{TARGET["symbol"]}[Title/Abstract] AND (inhibitor[Title/Abstract] OR "drug repurposing"[Title/Abstract])'
    try:
        lit = nova.tool_data(
            "literature.search_pubmed",
            {"query": query, "max_results": 10, "live_network": True,
             "contact_email": "jaysahni70@gmail.com", "timeout_seconds": 30.0},
        )
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        lit = {"error": str(exc), "articles": []}
        print(f"    pubmed unavailable: {exc}")

    # ATP-site annotation, straight from UniProt features. This is the site
    # definition used in stage 2. It is curated from the literature and is
    # independent of the ligand co-crystallised in the test structure, which is
    # why it may drive design where known_ligand_contacts.json may not
    # (that file is validation-only, task B12).
    site_residues, features = _uniprot_binding_site(TARGET["uniprot"])

    return write(
        "01_research.json",
        {
            "target": TARGET,
            "uniprot": {
                "accession": entry["accession"],
                "id": entry["uniprot_id"],
                "name": entry["primary_name"],
                "organism": entry["organism"],
                "sequence_length": len(entry["sequence"]),
                "sequence": entry["sequence"],
            },
            "annotated_site": {
                "source": "UniProt features (ft_binding, ft_act_site)",
                "ligands": sorted({f["ligand"] for f in features if f["ligand"]}),
                "residue_ids_author": site_residues,
                "n": len(site_residues),
                "features": features,
            },
            "literature": {
                "query": query,
                "n_hits": len(lit.get("articles", []) or []),
                "articles": (lit.get("articles") or [])[:10],
                "error": lit.get("error"),
            },
            "toolkit": ["database.fetch_uniprot_entry", "literature.search_pubmed"],
        },
    )


def _uniprot_binding_site(accession: str) -> tuple[list[int], list[dict]]:
    """Curated ATP/Mg binding residues for the target, from UniProt features."""
    import urllib.request

    url = (
        f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        "?fields=ft_binding,ft_act_site"
    )
    with urllib.request.urlopen(url, timeout=60) as handle:
        payload = json.load(handle)

    residues: set[int] = set()
    features = []
    for feature in payload.get("features", []):
        if feature["type"] not in ("Binding site", "Active site"):
            continue
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
        ligand = (feature.get("ligand") or {}).get("name")
        features.append(
            {"type": feature["type"], "start": start, "end": end, "ligand": ligand}
        )
        residues.update(range(start, end + 1))
    return sorted(residues), features


# ---------------------------------------------------------------------------
# stage 2 - site
# ---------------------------------------------------------------------------


def stage_site() -> dict:
    """Detect pockets and pick the one matching the annotated site.

    Toolkit: rowan.detect_pockets (Rowan-hosted), via the rowan SDK NovaKit pins.
    """
    print("[2/5] site")
    research = read("01_research.json")
    annotated_list = research["annotated_site"]["residue_ids_author"]
    annotated = set(annotated_list)
    rowan = nova.rowan_client()

    cached_protein = nova.cached(f'protein_{TARGET["pdb_id"]}')
    if cached_protein:
        raw_uuid = cached_protein["uuid"]
    else:
        protein = rowan.create_protein_from_pdb_id(
            TARGET["pdb_id"], name=f'{TARGET["symbol"]} {TARGET["pdb_id"]} demo'
        )
        raw_uuid = str(protein.uuid)
        nova.cache(f'protein_{TARGET["pdb_id"]}', {"uuid": raw_uuid})

    # 6Q4G is a HOLO structure, and create_protein_from_pdb_id keeps its ligand
    # and waters. Detecting pockets on an occupied pocket, then designing into
    # one, is PROJECT_GOAL.md 3 WS-B's "strip ligands" step skipped -- and it
    # measurably hurt: the first run of this demo designed against a pocket with
    # HJK still in it and every design came back at ipTM 0.13-0.20.
    protein_uuid = _stripped_protein(raw_uuid)
    print(f"  protein {protein_uuid} (apo, stripped from {raw_uuid})")

    chain = resolve_chain(
        _download_structure(protein_uuid, f'target_{TARGET["out"]}'),
        research["uniprot"]["sequence"],
    )
    print(f'  target chain {chain} (config hint was {TARGET["chain"]})')
    index_to_author = _residue_index_map(protein_uuid, chain)

    record = nova.run_workflow(
        f'pockets_{TARGET["out"]}',
        "submit_pocket_detection_workflow",
        label="pocket detection",
        max_credits=5,
        protein=protein_uuid,
        name=f'{TARGET["symbol"]} demo pockets',
    )

    pockets = []
    for i, pocket in enumerate(record["data"]["pockets"]):
        author = sorted(
            index_to_author[r] for r in pocket["residue_numbers"] if r in index_to_author
        )
        pockets.append(
            {
                "rank_by_score": i,
                "score": pocket["score"],
                "volume": pocket["volume"],
                "residue_ids_author": author,
                "residue_indices_0based": sorted(pocket["residue_numbers"]),
                "overlap_with_annotated_site": len(set(author) & annotated),
            }
        )

    chosen = max(pockets, key=lambda p: p["overlap_with_annotated_site"])
    best_overlap = chosen["overlap_with_annotated_site"]
    site_rule = "detected pocket with the most annotated-site residues"
    fallback = None

    if best_overlap < MIN_POCKET_OVERLAP:
        # Pocket detection did not find the functional site. Fall back to the
        # geometry around the annotated catalytic residues, and say so loudly.
        pdb = _download_structure(protein_uuid, f'target_{TARGET["out"]}')
        to_author = contacts.align_to_reference(pdb, chain, research["uniprot"]["sequence"])
        anchors = TARGET.get("site_anchors") or annotated_list
        pairs = geometric_site(pdb, chain, anchors, to_author)
        fallback = {
            "used": True,
            "reason": (
                f"best detected pocket overlapped the annotated site in only "
                f"{best_overlap} of {len(annotated_list)} residues, below the "
                f"threshold of {MIN_POCKET_OVERLAP}"
            ),
            "rule": f"all residues within {SITE_RADIUS_A} A of the annotated residues",
            "anchors_author": anchors,
            "radius_a": SITE_RADIUS_A,
            "n_pockets_containing_each_annotated_residue": {
                str(a): sum(1 for p in pockets if a in p["residue_ids_author"])
                for a in annotated_list
            },
        }
        chosen = {
            "rank_by_score": None,
            "score": None,
            "volume": None,
            "source": "geometric",
            "residue_ids_author": [a for a, _ in pairs],
            "residue_labels_structure": [str(lbl) for _, lbl in pairs],
        }
        site_rule = fallback["rule"]
        print(f"  pocket detection insufficient ({best_overlap}/{len(annotated_list)}); "
              f"using geometry around {anchors} at {SITE_RADIUS_A} A -> {len(pairs)} residues")

    return write(
        "02_site.json",
        {
            "protein_uuid": protein_uuid,
            "pdb_id": TARGET["pdb_id"],
            "target_chain": chain,
            "target_chain_hint": TARGET["chain"],
            "n_residues_in_structure": len(index_to_author),
            "n_pockets": len(pockets),
            "pockets": pockets,
            "chosen": chosen,
            "selection_rule": site_rule,
            "geometric_fallback": fallback,
            "finding": {
                "top_scoring_pocket_rank": 0,
                "top_scoring_pocket_overlap": pockets[0]["overlap_with_annotated_site"],
                "chosen_pocket_rank": chosen["rank_by_score"],
                "chosen_pocket_overlap": chosen.get("overlap_with_annotated_site"),
                "note": (
                    "Rowan's pocket score ranks the ATP site below pockets with no "
                    "annotated-site overlap at all. Taking the top-scoring pocket "
                    "would have designed against the wrong site."
                ),
            },
            "credits_charged": record.get("credits_charged"),
            "toolkit": ["rowan.detect_pockets"],
        },
    )


def _binder_chain(pdb: Path) -> str:
    """The designed chain in a BoltzGen complex.

    Not simply "B": Rowan returns the target as chain A, any retained heterogen
    as its own chain, and the design after that. The design is the polymer chain
    that is not the target.
    """
    protein, _, _ = contacts.parse_pdb(pdb)
    chains = sorted({ch for ch, _ in protein} - {TARGET["chain"]})
    if len(chains) != 1:
        raise ValueError(
            f"{pdb.name}: expected exactly one designed polymer chain beside "
            f'{TARGET["chain"]}, found {chains}'
        )
    return chains[0]


def _stripped_protein(raw_uuid: str) -> str:
    """A copy of the structure with heterogens and waters removed.

    `Protein.prepare` mutates in place, so the raw holo entry is copied first and
    the apo one kept under its own uuid -- the holo structure is still needed by
    the validation path, which reads its ligand.
    """
    hit = nova.cached(f'protein_stripped_{TARGET["out"]}')
    if hit:
        return hit["uuid"]

    rowan = nova.rowan_client()
    raw = rowan.retrieve_protein(raw_uuid)
    apo = rowan.create_protein_from_pdb_id(
        TARGET["pdb_id"], name=f'{TARGET["symbol"]} {TARGET["pdb_id"]} apo'
    )
    apo.prepare(
        remove_heterogens=True,
        keep_waters=False,
        find_missing_residues=False,
        add_missing_atoms=True,
        add_hydrogens=True,
        timeout=900.0,
    )
    del raw
    nova.cache(f'protein_stripped_{TARGET["out"]}', {"uuid": str(apo.uuid), "from": raw_uuid})
    return str(apo.uuid)


# Radius for the geometric site fallback. 6 A around the nucleophile yields ~20
# residues on thrombin, comparable to the 15-residue pocket detection found on
# CDK2 -- i.e. sized like a pocket, not like a domain. Set a priori from that
# size target; the overlap with the co-crystal ligand is reported afterwards as
# validation, never used to pick the radius.
SITE_RADIUS_A = float(os.environ.get("DEMO_SITE_RADIUS", "6.0"))        # around the annotated catalytic residues
MIN_POCKET_OVERLAP = 3      # below this, a detected pocket is not trusted to be the site


def geometric_site(
    pdb: Path, chain: str, anchors_author: list[int], to_author: dict, *, radius: float = SITE_RADIUS_A
) -> list:
    """Residues within `radius` of the annotated catalytic residues.

    The fallback when pocket detection cannot identify the functional site.

    MEASURED on thrombin: Rowan returns 6 pockets for 1PPB chain B and **none**
    contains Ser195 or Asp102; the best annotated-site overlap is 1 of 3, and two
    pockets tie at it. That is not a failure of the detector so much as a fact
    about serine proteases -- the active site is a shallow groove split across
    S1/S2/S3 subsites rather than one enclosed cavity, so a cavity finder
    fragments it. CDK2's ATP site, being a real cavity, overlapped 11 of 20.

    The catalytic triad is a *functional annotation*, not a ligand, so anchoring
    on it does not leak the answer the way using the co-crystal ligand would
    (task B12). This is the same reasoning that let the UniProt site drive CDK2.
    """
    protein, _, _ = contacts.parse_pdb(pdb)
    author_of = {(chain, label): to_author.get(label) for (ch, label) in protein if ch == chain}

    anchor_coords = [
        coords for (ch, label), (_, coords) in protein.items()
        if ch == chain and author_of.get((ch, label)) in set(anchors_author)
    ]
    if not anchor_coords:
        raise SystemExit(
            f"{pdb.name}: none of the annotated residues {anchors_author} could be "
            "located in the structure; refusing to define a site."
        )
    anchor_xyz = np.vstack(anchor_coords)

    site = []
    for (ch, label), (_, coords) in protein.items():
        if ch != chain:
            continue
        author = author_of.get((ch, label))
        if author is None:
            continue
        d = float(np.sqrt(((coords[:, None, :] - anchor_xyz[None, :, :]) ** 2).sum(-1)).min())
        if d <= radius:
            site.append((author, label))
    return sorted(site, key=lambda pair: pair[0])


def resolve_chain(pdb: Path, reference: str) -> str:
    """The chain holding the target domain, found by sequence, not by name.

    Chain letters are not stable across this pipeline. `Protein.prepare` renames
    1PPB's chains alphabetically -- the catalytic heavy chain H becomes B and the
    light chain L becomes A -- so the letter in TARGETS is a hint about the
    original PDB, not a fact about the prepared one. CDK2 is single-chain so this
    returns "A" there either way.

    The chain that aligns to the most of the canonical sequence wins; a tie or a
    failure to align anything is an error rather than a guess.
    """
    protein, _, _ = contacts.parse_pdb(pdb)
    best, best_n = None, 0
    report = {}
    for chain in sorted({ch for ch, _ in protein}):
        try:
            mapped = contacts.align_to_reference(pdb, chain, reference)
        except ValueError:
            report[chain] = "no usable alignment"
            continue
        report[chain] = f"{len(mapped)} residues -> {min(mapped.values())}-{max(mapped.values())}"
        if len(mapped) > best_n:
            best, best_n = chain, len(mapped)

    if best is None:
        raise SystemExit(f"{pdb.name}: no chain aligns to the reference sequence. {report}")
    return best


def _residue_index_map(protein_uuid: str, chain: str | None = None) -> dict[int, int]:
    """0-based index over chain-A residues present in the file -> author number.

    Rowan reports pocket residues by position in the file, and `Protein.prepare`
    renumbers the file contiguously from 1, so neither the file's own numbers nor
    a constant offset gives author numbering. The map is therefore recovered by
    aligning the structure's sequence to the canonical one.
    """
    pdb = _download_structure(protein_uuid, f'target_{TARGET["out"]}')
    sequence = read("01_research.json")["uniprot"]["sequence"]
    chain = chain or resolve_chain(pdb, sequence)
    numbers, _ = contacts.chain_sequence(pdb, chain)
    to_author = contacts.align_to_reference(pdb, chain, sequence)
    return {i: to_author[resseq] for i, resseq in enumerate(numbers) if resseq in to_author}


def _download_structure(uuid: str, name: str) -> Path:
    """Fetch a Rowan structure as PDB, cached on disk."""
    STRUCTURES.mkdir(parents=True, exist_ok=True)
    path = STRUCTURES / f"{name}.pdb"
    if path.exists():
        return path
    rowan = nova.rowan_client()
    rowan.retrieve_protein(uuid).download_pdb_file(STRUCTURES, name)
    return path


# ---------------------------------------------------------------------------
# stage 3 - design
# ---------------------------------------------------------------------------


def stage_design() -> dict:
    """Run BoltzGen against the chosen pocket.

    Toolkit: rowan.design_protein_binder -> Rowan `protein_binder_design`, which
    the account's feature flags identify as BoltzGen (`boltzgen_design_limit_100`).
    """
    print("[3/5] design")
    site = read("02_site.json")

    # BoltzGen indexes residues 1..N over residues present in the file.
    binding = ",".join(str(i + 1) for i in site["chosen"]["residue_indices_0based"])

    binder_design_input = {
        "constraints": [],
        "ligand_entities": [],
        "protein_entities": [{"id": "B", "sequence": "60..90"}],
        "file_entities": [
            {
                "uuid": site["protein_uuid"],
                "design": [],
                "exclude": [],
                "include": [],
                "binding_types": [{"chain_id": TARGET["chain"], "binding": binding}],
                "design_insertions": [],
                "include_proximity": [
                    {
                        "radius": 12,
                        "chain_id": TARGET["chain"],
                        "residue_indices": binding,
                    }
                ],
                "secondary_structure": [],
            }
        ],
    }

    record = nova.run_workflow(
        "designs",
        "submit_protein_binder_design_workflow",
        label=f"BoltzGen x{N_DESIGNS}",
        max_credits=120,
        binder_design_input=binder_design_input,
        protocol="protein-anything",
        num_designs=N_DESIGNS,
        budget=DESIGN_BUDGET,
        name=f'{TARGET["symbol"]} lightweight demo - BoltzGen',
    )

    binders = record["data"]["generated_binders"]
    designs = []
    for i, binder in enumerate(binders):
        scores = binder.get("scores", {})
        designs.append(
            {
                "design_id": i,
                "bound_structure_uuid": binder.get("bound_structure"),
                "sequence": (binder.get("binder_sequences") or [{}])[0].get("sequence"),
                "iptm": scores.get("iptm"),
                "design_to_target_iptm": scores.get("design_to_target_iptm"),
                "quality_score": scores.get("quality_score"),
                "delta_sasa_refolded": scores.get("delta_sasa_refolded"),
                "num_filters_passed": scores.get("num_filters_passed"),
            }
        )

    iptms = [d["iptm"] for d in designs if d["iptm"] is not None]
    return write(
        "03_designs.json",
        {
            "workflow_uuid": record["uuid"],
            "n_designs": len(designs),
            "binding_spec_1based": binding,
            "binding_residues_author": site["chosen"]["residue_ids_author"],
            "designs": designs,
            "iptm": {
                "min": min(iptms) if iptms else None,
                "max": max(iptms) if iptms else None,
                "mean": round(sum(iptms) / len(iptms), 4) if iptms else None,
                "n_above_0.85": sum(1 for v in iptms if v >= 0.85),
                "note": (
                    "PROJECT_GOAL.md 4.3 suggests filtering designs at ipTM 0.85. "
                    "Report the count above it; if it is zero the consensus below "
                    "is UNFILTERED and must be labelled as such."
                ),
            },
            "credits_charged": record.get("credits_charged"),
            "compute_hardware": record.get("compute_hardware"),
            "elapsed_s": record.get("elapsed"),
            "toolkit": ["rowan.design_protein_binder"],
        },
    )


# ---------------------------------------------------------------------------
# stage 4 - signature
# ---------------------------------------------------------------------------


def stage_signature() -> dict:
    """Turn the designs into a consensus target-side interface signature."""
    print("[4/5] signature")
    designs = read("03_designs.json")
    site = read("02_site.json")
    sequence = read("01_research.json")["uniprot"]["sequence"]

    per_design = []
    failures = []
    for design in designs["designs"]:
        uuid = design["bound_structure_uuid"]
        if not uuid:
            failures.append({"design_id": design["design_id"], "reason": "no bound structure"})
            continue
        try:
            pdb = _download_structure(uuid, f'design_{design["design_id"]}')
            binder_chain = _binder_chain(pdb)
            # The design complex renumbers the target by a constant offset, and
            # it is NOT an index into the prepared structure. Measure it against
            # the canonical sequence rather than assuming.
            offset = contacts.author_offset(pdb, TARGET["chain"], sequence)
            found = contacts.contacts_between_chains(
                pdb, target=TARGET["chain"], binder=binder_chain
            )
        except Exception as exc:  # noqa: BLE001
            failures.append({"design_id": design["design_id"], "reason": str(exc)})
            continue

        author = sorted(r + offset for r in found["residues"])
        per_design.append(
            {
                "design_id": design["design_id"],
                "binder_chain": binder_chain,
                "author_offset": offset,
                "residues_author": author,
            }
        )

    if not per_design:
        raise SystemExit(
            f"no design yielded contacts. failures: {failures}. "
            "Refusing to build a signature from nothing."
        )

    signature = contacts.consensus(
        [d["residues_author"] for d in per_design], core_threshold=CORE_THRESHOLD
    )

    # Validation only, never used to build the signature (task B12).
    known = json.loads(
        (REPO / "results/pipeline/cdk2-second-target/target/known_ligand_contacts.json").read_text()
    )
    known_ids = known["known_ligand_residue_ids"]
    annotated = read("01_research.json")["annotated_site"]["residue_ids_author"]

    unfiltered = designs["iptm"]["n_above_0.85"] == 0
    return write(
        "04_signature.json",
        {
            "source": "boltzgen_consensus",
            "n_designs_used": len(per_design),
            "n_designs_failed": len(failures),
            "failures": failures,
            "per_design": per_design,
            **signature,
            "validation": {
                "_note": "computed AFTER the signature, never used to build it",
                "known_ligand_residue_ids": known_ids,
                "jaccard_core_vs_known_ligand": round(
                    contacts.jaccard(signature["core"], known_ids), 4
                ),
                "recall_of_known_ligand_contacts": round(
                    len(set(signature["core"]) & set(known_ids)) / len(known_ids), 4
                )
                if known_ids
                else None,
                "jaccard_chosen_pocket_vs_known_ligand": round(
                    contacts.jaccard(site["chosen"]["residue_ids_author"], known_ids), 4
                ),
                "jaccard_annotated_site_vs_known_ligand": round(
                    contacts.jaccard(annotated, known_ids), 4
                ),
            },
            "caveats": [
                "UNFILTERED consensus: no design reached the ipTM 0.85 gate."
                if unfiltered
                else "filtered at ipTM 0.85",
                f"consensus over {len(per_design)} designs from ONE structure "
                f'({TARGET["pdb_id"]}).',
                "results/m2_gate_CDK2.json measured this arm against pocket geometry "
                "on 31 held-out CDK2 co-crystals and BoltzGen LOST "
                "(Jaccard 0.338 vs 0.494, Holm p=0.0016). Stage 5 therefore scores "
                "drugs against both signatures, not just this one.",
            ],
        },
    )


# ---------------------------------------------------------------------------
# stage 5 - match
# ---------------------------------------------------------------------------


def stage_match() -> dict:
    """Co-fold each shortlisted approved drug and rank by hotspot coverage."""
    print("[5/5] match")
    research = read("01_research.json")
    site = read("02_site.json")
    signature = read("04_signature.json")

    sequence = research["uniprot"]["sequence"]
    shortlist = library.load(TARGET["symbol"])
    print(f"  {len(shortlist)} drugs, ~{COFOLD_MAX_CREDITS * len(shortlist)} credits worst case")

    # Two signatures, scored side by side. The pocket arm is the honest control:
    # results/m2_gate_CDK2.json says it beats the BoltzGen arm on this target.
    boltzgen_core = signature["core"]
    pocket_core = site["chosen"]["residue_ids_author"]

    results = []
    for drug in shortlist:
        key = f'cofold_{drug["name"]}'
        record = nova.run_workflow(
            key,
            "submit_protein_cofolding_workflow",
            label=f'cofold {drug["name"]}',
            max_credits=COFOLD_MAX_CREDITS,
            initial_protein_sequences=[sequence],
            initial_smiles_list=[drug["smiles"]],
            use_msa_server=True,
            name=f'{TARGET["symbol"]} demo cofold - {drug["name"]}',
        )
        data = record["data"]
        scores = data.get("scores") or {}
        structure_uuid = data.get("predicted_structure_uuid")

        pdb = _download_structure(structure_uuid, f'cofold_{drug["name"]}')
        found = contacts.contacts_to_ligand(pdb, chain=TARGET["chain"])
        # Cofold numbers the protein 1..N over the supplied sequence, which is the
        # UniProt canonical sequence -- so residue i is author residue i.
        engaged = found["residues"]

        results.append(
            {
                **drug,
                "workflow_uuid": record["uuid"],
                "structure_uuid": structure_uuid,
                "iptm": scores.get("iptm"),
                "ptm": scores.get("ptm"),
                "confidence_score": scores.get("confidence_score"),
                "n_engaged_residues": len(engaged),
                "engaged_residues": engaged,
                "ligand_atoms": found["n_ligand_atoms"],
                "core_coverage_boltzgen": round(contacts.coverage(engaged, boltzgen_core), 4),
                "core_coverage_pocket": round(contacts.coverage(engaged, pocket_core), 4),
                "engaged_core_boltzgen": sorted(set(engaged) & set(boltzgen_core)),
                "missed_core_boltzgen": sorted(set(boltzgen_core) - set(engaged)),
                "engaged_core_pocket": sorted(set(engaged) & set(pocket_core)),
                "missed_core_pocket": sorted(set(pocket_core) - set(engaged)),
                "credits_charged": record.get("credits_charged"),
            }
        )

    for row in results:
        row["rank_boltzgen"] = 0
        row["rank_pocket"] = 0
    for field, rank_field in (
        ("core_coverage_boltzgen", "rank_boltzgen"),
        ("core_coverage_pocket", "rank_pocket"),
    ):
        for i, row in enumerate(sorted(results, key=lambda r: -r[field]), start=1):
            row[rank_field] = i

    results.sort(key=lambda r: -r["core_coverage_pocket"])

    return write(
        "05_ranked.json",
        {
            "target": TARGET,
            "shortlist": library.summarise(shortlist),
            "signatures": {
                "boltzgen_consensus": {"n_core": len(boltzgen_core), "core": boltzgen_core},
                "pocket_geometry": {"n_core": len(pocket_core), "core": pocket_core},
            },
            "ranked_by": "core_coverage_pocket",
            "results": results,
            "credits_charged_total": round(
                sum(r["credits_charged"] or 0 for r in results), 2
            ),
            "toolkit": ["rowan.cofold"],
        },
    )


# ---------------------------------------------------------------------------


STAGES = {
    "research": stage_research,
    "site": stage_site,
    "design": stage_design,
    "signature": stage_signature,
    "match": stage_match,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", default="all", choices=[*STAGES, "all"], help="stage to run"
    )
    args = parser.parse_args()

    names = list(STAGES) if args.stage == "all" else [args.stage]
    for name in names:
        STAGES[name]()
    print(f"\ncredits charged this run: {nova.spent():.2f}")


if __name__ == "__main__":
    main()
