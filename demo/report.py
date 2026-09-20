"""Render the demo's computed artifacts as a report.

Reads only what the pipeline wrote. Every number here comes from a JSON file in
results/demo/cdk2/; nothing is illustrative and nothing is filled in ahead of a
run (CLAUDE.md working agreement).

Run:  ./env-kit/bin/python -m demo.report
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "demo" / "cdk2"


def load(name: str) -> dict:
    path = OUT / name
    if not path.exists():
        raise SystemExit(f"{path.relative_to(REPO)} missing -- run the pipeline first.")
    return json.loads(path.read_text())


def main() -> None:
    research = load("01_research.json")
    site = load("02_site.json")
    designs = load("03_designs.json")
    signature = load("04_signature.json")
    ranked = load("05_ranked.json")

    target = ranked["target"]
    results = ranked["results"]
    lines: list[str] = []
    add = lines.append

    add(f'# {target["symbol"]} — lightweight AutoRepurpose demo')
    add("")
    add(
        f'Target **{target["symbol"]}** ({target["uniprot"]}, '
        f'{research["uniprot"]["sequence_length"]} aa), design structure '
        f'**{target["pdb_id"]}**. Computed with the `novakit` toolkit; every number '
        "below comes from a file in this directory."
    )
    add("")

    # -- pipeline ---------------------------------------------------------
    add("## What ran")
    add("")
    add("| Stage | Toolkit call | Result | Credits |")
    add("|---|---|---|---|")
    add(
        f'| research | `database.fetch_uniprot_entry`, `literature.search_pubmed` | '
        f'{research["annotated_site"]["n"]} annotated site residues, '
        f'{research["literature"]["n_hits"]} PubMed hits | 0 |'
    )
    add(
        f'| site | `rowan.detect_pockets` | {site["n_pockets"]} pockets, chose rank '
        f'{site["chosen"]["rank_by_score"]} | {site.get("credits_charged")} |'
    )
    add(
        f'| design | `rowan.design_protein_binder` (BoltzGen) | '
        f'{designs["n_designs"]} designs | {designs.get("credits_charged")} |'
    )
    add(
        f'| signature | `demo/contacts.py` | {len(signature["core"])} core residues '
        f'from {signature["n_designs_used"]} designs | 0 |'
    )
    add(
        f'| match | `rowan.cofold` | {len(results)} approved drugs co-folded | '
        f'{ranked["credits_charged_total"]} |'
    )
    add("")

    # -- site finding -----------------------------------------------------
    add("## The ATP site is not the top-scoring pocket")
    add("")
    add("| Pocket | Score | Volume Å³ | Residues | Overlap with annotated site |")
    add("|---|---|---|---|---|")
    for pocket in site["pockets"]:
        mark = " ←chosen" if pocket["rank_by_score"] == site["chosen"]["rank_by_score"] else ""
        add(
            f'| {pocket["rank_by_score"]}{mark} | {pocket["score"]:.2f} | '
            f'{pocket["volume"]:.0f} | {len(pocket["residue_ids_author"])} | '
            f'**{pocket["overlap_with_annotated_site"]}** |'
        )
    add("")
    add(f'{site["finding"]["note"]}')
    add("")
    add(
        f'Chosen site, author numbering: `{site["chosen"]["residue_ids_author"]}`. '
        f'Selection rule: {site["selection_rule"]}'
    )
    add("")

    # -- designs ----------------------------------------------------------
    add("## Designs")
    add("")
    iptm = designs["iptm"]
    add(
        f'{designs["n_designs"]} BoltzGen designs, ipTM '
        f'{iptm["min"]}–{iptm["max"]} (mean {iptm["mean"]}). '
        f'**{iptm["n_above_0.85"]} of {designs["n_designs"]}** reached the ipTM 0.85 '
        "gate PROJECT_GOAL.md §4.3 suggests filtering at."
    )
    add("")
    add(
        f'Consensus over {signature["n_designs_used"]} designs at frequency ≥ '
        f'{signature["core_threshold"]}: **{len(signature["core"])} core residues**. '
        f'Mean pairwise Jaccard between designs: **{signature["mean_pairwise_jaccard"]}** '
        "— how much the designs agree with each other at all."
    )
    add("")
    add(f'Core residues: `{signature["core"]}`')
    add("")

    validation = signature["validation"]
    add("### Scored against the known ligand (validation only)")
    add("")
    add("| Signature | Jaccard vs known-ligand contacts |")
    add("|---|---|")
    add(f'| BoltzGen consensus | {validation["jaccard_core_vs_known_ligand"]} |')
    add(f'| chosen pocket geometry | {validation["jaccard_chosen_pocket_vs_known_ligand"]} |')
    add(f'| UniProt annotated site | {validation["jaccard_annotated_site_vs_known_ligand"]} |')
    add("")
    add(
        f'Known-ligand contacts (n={len(validation["known_ligand_residue_ids"])}) are read '
        "**after** the signature is built and never used to construct it (task B12)."
    )
    add("")

    # -- leaderboard ------------------------------------------------------
    add("## Approved drugs, ranked")
    add("")
    add(
        f'{ranked["shortlist"]["n"]} drugs of '
        f'{ranked["shortlist"]["library_size"]} in `{ranked["shortlist"]["library_file"]}`, '
        f'by tier: {ranked["shortlist"]["by_tier"]}. Co-folded with Boltz-2 via '
        "`rowan.cofold`; score is the fraction of core hotspot residues the docked pose "
        "engages within 4.5 Å."
    )
    add("")
    add("| # | Drug | Tier | CDK2 annotated | Coverage (pocket) | Coverage (BoltzGen) | ipTM | Engaged residues |")
    add("|---|---|---|---|---|---|---|---|")
    for i, row in enumerate(results, start=1):
        add(
            f'| {i} | {row["name"]} | {row["tier"]} | '
            f'{"yes" if row["hits_target"] else "no"} | '
            f'**{row["core_coverage_pocket"]:.2f}** | {row["core_coverage_boltzgen"]:.2f} | '
            f'{row["iptm"]} | {row["n_engaged_residues"]} |'
        )
    add("")

    # -- tier separation --------------------------------------------------
    add("### Does the score separate the tiers?")
    add("")
    add("| Tier | n | Mean coverage (pocket) | Mean coverage (BoltzGen) | Mean ipTM |")
    add("|---|---|---|---|---|")
    for tier in ("positive", "candidate", "decoy"):
        rows = [r for r in results if r["tier"] == tier]
        if not rows:
            continue
        add(
            f"| {tier} | {len(rows)} | "
            f'{sum(r["core_coverage_pocket"] for r in rows) / len(rows):.3f} | '
            f'{sum(r["core_coverage_boltzgen"] for r in rows) / len(rows):.3f} | '
            f'{sum((r["iptm"] or 0) for r in rows) / len(rows):.3f} |'
        )
    add("")

    # -- caveats ----------------------------------------------------------
    add("## Caveats")
    add("")
    for caveat in signature["caveats"]:
        add(f"- {caveat}")
    add(
        "- **Leakage.** Boltz-2 co-folding an approved drug against a target it is "
        "already annotated to bind is very likely reproducing a complex in its "
        "training set. The positive tier's scores are therefore an upper bound, not "
        "a blind prediction."
    )
    add(
        f'- **n = {len(results)} drugs, one target, one structure.** No significance '
        "test is reported because none is meaningful at this n; the tier means above "
        "are descriptive."
    )
    add(
        "- Output is a computationally-ranked, evidence-linked repurposing "
        "hypothesis. It requires experimental validation."
    )
    add("")

    path = OUT / "report.md"
    path.write_text("\n".join(lines) + "\n")
    print(f"-> {path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
