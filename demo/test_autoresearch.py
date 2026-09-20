"""Self-test for demo/autoresearch.py.

The three things that can silently be wrong here, in the order they would hurt:

1. THE POSITIVE-CONTROL GATE. If it lets a target through that has no approved
   drug binding it, the downstream matching arms produce a ranking of 2,153
   approved drugs with no known-correct answer in it, and every percentile they
   report is uninterpretable. The test case is CDK2, because CDK2 *looks* like
   it passes: seven approved drugs carry a CDK2 annotation. None of them has
   CDK2 as its mechanism-of-action target - they are CDK4/6 and multi-kinase
   drugs picking CDK2 up off-target. A gate counting annotations would wave it
   through; the gate here must not.

2. THE PEPTIDE-POLYMER STRUCTURE GUARD. A co-crystallised peptide is ATOM
   records in a polymer chain, so `remove_heterogens` leaves it in the site and
   the design stage designs against an occupied pocket. The test asserts that
   4UD9 - thrombin's 1.12 A entry, the one a resolution-only rule would pick,
   and one this repo already rejected by hand - is rejected for its 12-residue
   hirudin chain, and that the entry actually chosen has no such chain.

3. THE CATALYTIC ANCHOR. demo/pipeline.py anchors thrombin's site on Ser195
   (UniProt 568) alone, by hand, because the His/Asp of the triad sit behind it.
   The rule here must reproduce that from UniProt without being told, and must
   return nothing for CDK2, whose only ACT_SITE is an aspartate proton acceptor.

Plus the boring but load-bearing ones: the `approved` column really is a filter
(4,099 rows, 2,153 approved), a known disease really does return its known
targets, and the emitted config really is the shape demo/pipeline.py consumes.

Run:  ./env-kit/bin/python -m demo.test_autoresearch
Free: no credentials, no credits. Network is used only where data/raw/ has no
cache yet; a second run is fully offline.
"""

from __future__ import annotations

import ast
from pathlib import Path

from demo import autoresearch as A

ROOT = Path(__file__).resolve().parent.parent
PIPELINE = ROOT / "demo" / "pipeline.py"


def pipeline_targets():
    """demo/pipeline.py's TARGETS dict, read without importing the module."""
    tree = ast.parse(PIPELINE.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "TARGETS" for t in node.targets
        ):
            return ast.literal_eval(node.value)
    raise AssertionError("no TARGETS assignment in demo/pipeline.py")


def main() -> int:
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), str(detail)))

    # ------------------------------------------------------------------
    # 1. the corpora, and the `approved` column as a real filter
    # ------------------------------------------------------------------
    c = A.load_corpora()
    check(
        "approved_drugs.csv filtered on approved==1",
        c["n_drug_rows"] == 4099 and c["n_approved"] == 2153,
        f"{c['n_approved']}/{c['n_drug_rows']} rows",
    )
    check(
        "approved_biologics.csv loaded",
        c["n_biologics"] == 262,
        f"n={c['n_biologics']}",
    )
    check(
        "target_annotations.csv restricted to approved + potent",
        0 < c["n_annotation_rows_approved_and_potent"] < c["n_annotation_rows"],
        f"{c['n_annotation_rows_approved_and_potent']}/{c['n_annotation_rows']} rows kept",
    )

    # ------------------------------------------------------------------
    # 2. the positive-control gate: CDK2 must fail, thrombin must pass
    # ------------------------------------------------------------------
    cdk2 = A.approved_binders("CDK2", "P24941", ["CDKN2"], c)
    f2 = A.approved_binders("F2", "P00734", ["PT", "THPH1"], c)

    check(
        "CDK2 has no approved mechanism-of-action drug",
        cdk2["n_approved_binders_moa"] == 0,
        f"moa={cdk2['n_approved_binders_moa']}",
    )
    check(
        "CDK2 would have passed a naive annotation count",
        cdk2["n_approved_binders_annotated"] > 0
        and cdk2["n_approved_binders_bioactivity"] > 0,
        f"annotated={cdk2['n_approved_binders_annotated']} "
        f"({', '.join(cdk2['approved_binders_offtarget_only'][:4])}...), "
        f"pchembl={cdk2['n_approved_binders_bioactivity']}",
    )
    check(
        "F2 has approved mechanism-of-action drugs",
        f2["n_approved_binders_moa"] >= 4
        and {"argatroban", "dabigatran etexilate"} <= set(f2["known_approved_drugs"]),
        f"moa={f2['n_approved_binders_moa']} {f2['known_approved_drugs']}",
    )
    check(
        "F2 has an approved peptide binder for the peptide arm",
        f2["n_approved_peptide_binders"] >= 1
        and any(b["name"] == "bivalirudin" for b in f2["approved_biologic_binders"]),
        f"peptides={f2['n_approved_peptide_binders']} "
        f"{[b['name'] for b in f2['approved_biologic_binders']]}",
    )
    check(
        "CDK2 has no approved biologic binder",
        cdk2["n_approved_biologic_binders"] == 0,
        f"n={cdk2['n_approved_biologic_binders']}",
    )

    # ------------------------------------------------------------------
    # 3. the peptide-polymer structure guard
    # ------------------------------------------------------------------
    ids, n_total, n_holo = A.search_entities_by_resolution("P00734")
    ents = A.fetch_entity_details("P00734", ids)
    entries = {e["pdb_id"]: e for e in A.entry_records(ents, "P00734")}
    chosen, usable, rejected = A.choose_structure(list(entries.values()))
    rej = {r["pdb_id"]: r for r in rejected}

    check(
        "4UD9 (1.12 A, hirudin chain) is rejected",
        "4UD9" in rej
        and rej["4UD9"]["reject_reason"] == "co_crystallised_peptide_polymer_chain",
        rej.get("4UD9", {}).get("reject_reason", "4UD9 not among the entries examined"),
    )
    check(
        "the hirudin chain is what triggered it",
        "4UD9" in rej
        and any(
            "Hirudin" in (d or "") for d in rej["4UD9"].get("reject_detail", [])
        ),
        rej.get("4UD9", {}).get("reject_detail"),
    )
    check(
        "thrombin's own 28-aa light chain does NOT trigger it",
        "4UD9" in entries
        and len(entries["4UD9"]["foreign_peptide_chains"]) == 1
        and any(
            (t["length_aa"] or 0) < A.PEPTIDE_POLYMER_MAX_AA
            for t in entries["4UD9"]["target_entities"]
        ),
        f"{len(entries.get('4UD9', {}).get('foreign_peptide_chains', []))} foreign peptide, "
        f"{len(entries.get('4UD9', {}).get('target_entities', []))} target entities",
    )
    check(
        "a structure was chosen for P00734",
        chosen is not None,
        chosen and f"{chosen['pdb_id']} chain {chosen['chain']} "
        f"{chosen['resolution_a']} A {chosen['apo_or_holo']}",
    )
    check(
        "the chosen entry carries no foreign peptide chain",
        chosen is not None and not chosen["foreign_peptide_chains"],
        chosen and str(chosen["foreign_peptide_chains"]),
    )
    check(
        "the chosen chain is the 259-aa heavy chain, not the light chain",
        chosen is not None and chosen["chain_entity_length_aa"] == 259,
        chosen and f"chain {chosen['chain']} = {chosen['chain_entity_length_aa']} aa",
    )
    check(
        "the guard rejected entries that resolution alone would have taken",
        chosen is not None
        and any(
            r["resolution_a"] is not None and r["resolution_a"] < chosen["resolution_a"]
            for r in rejected
        ),
        f"{len(rejected)}/{len(entries)} entries rejected, "
        f"{sum(1 for r in rejected if (r['resolution_a'] or 9) < (chosen or {}).get('resolution_a', 0))} "
        "of them at better resolution than the pick",
    )

    # CDK2's hand-picked structure must survive the same rule.
    cdk2_struct = A.structure_for("P24941")
    check(
        "CDK2 selection reproduces the hand-picked 6Q4G chain A",
        cdk2_struct["structure"]
        and cdk2_struct["structure"]["pdb_id"] == "6Q4G"
        and cdk2_struct["structure"]["chain"] == "A",
        cdk2_struct["structure"]
        and f"{cdk2_struct['structure']['pdb_id']} {cdk2_struct['structure']['chain']} "
        f"{cdk2_struct['structure']['resolution_a']} A",
    )

    # ------------------------------------------------------------------
    # 4. catalytic anchors
    # ------------------------------------------------------------------
    s_f2 = A.uniprot_sites("P00734")
    s_cdk2 = A.uniprot_sites("P24941")
    act_f2 = {a["position"]: a["residue"] for a in s_f2["active_sites"]}
    check(
        "F2 anchors reproduce pipeline.py's hand-chosen Ser195 (UniProt 568)",
        s_f2["site_anchors"] == [568],
        f"anchors={s_f2['site_anchors']} from ACT_SITE {act_f2}",
    )
    check(
        "the His/Asp of the triad are excluded",
        406 in act_f2 and 462 in act_f2 and 406 not in s_f2["site_anchors"],
        f"406={act_f2.get(406)} 462={act_f2.get(462)} 568={act_f2.get(568)}",
    )
    check(
        "CDK2 gets no anchor (its ACT_SITE is an Asp proton acceptor)",
        s_cdk2["site_anchors"] == []
        and [a["residue"] for a in s_cdk2["active_sites"]] == ["D"],
        f"anchors={s_cdk2['site_anchors']} act_sites="
        f"{[(a['position'], a['residue']) for a in s_cdk2['active_sites']]}",
    )
    check(
        "F2 mature chains recovered (heavy chain 364-622)",
        any(
            ch["start"] == 364 and ch["end"] == 622
            for ch in s_f2.get("mature_chains", [])
        ),
        str(s_f2.get("mature_chains")),
    )

    # ------------------------------------------------------------------
    # 5. circularity arithmetic
    # ------------------------------------------------------------------
    circ = A.circularity({"clinical": 0.9, "literature": 0.3, "genetic_association": 0.3}, "clinical")
    check(
        "known-drug share computed over the datatype scores",
        abs(circ["known_drug_share_of_datatype_sum"] - 0.6) < 1e-6
        and circ["known_drug_dominated"],
        f"share={circ['known_drug_share_of_datatype_sum']} dominated={circ['known_drug_dominated']}",
    )
    circ2 = A.circularity({"literature": 0.9, "clinical": 0.05}, "clinical")
    check(
        "a small known-drug contribution is not flagged",
        not circ2["known_drug_dominated"] and circ2["largest_datatype"] == "literature",
        f"share={circ2['known_drug_share_of_datatype_sum']}",
    )
    circ3 = A.circularity({"genetic_association": 0.8}, None)
    check(
        "a target with no known-drug datatype is not flagged",
        not circ3["known_drug_dominated"] and circ3["known_drug_score"] is None,
        circ3["circularity_note"],
    )

    # ------------------------------------------------------------------
    # 6. a known disease returns its known targets
    # ------------------------------------------------------------------
    # write=False: the self-test must not clobber the run the results directory
    # holds (this one skips literature to stay quick).
    out, _ = A.run("venous thromboembolism", top=25, lit_per_target=0, refresh=False,
                   write=False, verbose=False)
    syms = [t["approved_symbol"] for t in out["targets"]]
    check(
        "VTE resolves to MONDO_0005399 without being told",
        out["meta"]["efo_id"] == "MONDO_0005399",
        f"{out['meta']['efo_id']} {out['meta']['disease_name']!r} "
        f"via {out['meta']['disease_resolution_rule']}",
    )
    coag = {"F2", "F10", "F11", "F5", "PROC", "PROS1", "SERPINC1"}
    check(
        "VTE top-25 is the coagulation cascade",
        len(coag & set(syms)) >= 5,
        f"{sorted(coag & set(syms))} of 25",
    )
    f2h = next(t for t in out["targets"] if t["approved_symbol"] == "F2")
    check(
        "F2 passes both gates",
        f2h["usable_for_pipeline"],
        f"control={f2h['gate_has_approved_positive_control']} "
        f"structure={f2h['gate_has_usable_structure']}",
    )
    unusable = [t for t in out["targets"] if not t["usable_for_pipeline"]]
    check(
        "every unusable target says why, and gets no config",
        unusable
        and all(t["unusable_reasons"] and t["pipeline_config"] is None for t in unusable),
        f"{len(unusable)}/{len(syms)} unusable, "
        f"e.g. {unusable[0]['approved_symbol']}: {unusable[0]['unusable_reasons'][0][:60]}...",
    )
    top_assoc = out["targets"][0]
    check(
        "the top-ranked association is not silently promoted when it fails a gate",
        top_assoc["usable_for_pipeline"]
        or top_assoc["approved_symbol"] in out["meta"]["targets_without_positive_control"]
        + out["meta"]["targets_without_usable_structure"],
        f"rank 1 = {top_assoc['approved_symbol']}, usable={top_assoc['usable_for_pipeline']}",
    )

    # ------------------------------------------------------------------
    # 7. the emitted config is the shape demo/pipeline.py consumes
    # ------------------------------------------------------------------
    pt = pipeline_targets()
    required = {"symbol", "uniprot", "pdb_id", "chain", "out", "protocol", "binder_length"}
    known_keys = set().union(*(set(v) for v in pt.values()))
    cfg = f2h["pipeline_config"]
    check(
        "config has every key pipeline.py's entries share",
        cfg and required <= set(cfg),
        f"missing {sorted(required - set(cfg or {}))}",
    )
    check(
        "config invents no key pipeline.py does not know",
        cfg and set(cfg) <= known_keys,
        f"unknown {sorted(set(cfg or {}) - known_keys)}",
    )
    check(
        "F2 config matches the hand-written thrombin entry on everything but the PDB entry",
        cfg
        and cfg["symbol"] == pt["thrombin"]["symbol"]
        and cfg["uniprot"] == pt["thrombin"]["uniprot"]
        and cfg["protocol"] == pt["thrombin"]["protocol"]
        and cfg["binder_length"] == pt["thrombin"]["binder_length"]
        and cfg.get("site_anchors") == pt["thrombin"]["site_anchors"],
        f"auto={cfg['pdb_id']}/{cfg['chain']} vs hand={pt['thrombin']['pdb_id']}"
        f"/{pt['thrombin']['chain']}; protocol {cfg['protocol']} {cfg['binder_length']}, "
        f"anchors {cfg.get('site_anchors')}",
    )
    check(
        "the peptide protocol follows from an approved peptide binder, not from the symbol",
        cfg["protocol"] == "peptide-anything"
        and f2h["n_approved_peptide_binders"] >= 1,
        f"{f2h['n_approved_peptide_binders']} approved peptide binders",
    )

    # ------------------------------------------------------------------
    failed = 0
    for name, ok, detail in checks:
        print(f'  {"PASS" if ok else "FAIL"}  {name:62s} {detail}')
        failed += not ok
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
