"""Tests for the pipeline's stage interfaces, not its science.

Each test here exists because the thing it checks was actually broken at some
point, and the breakage was silent. They run offline against committed artifacts;
none of them calls Rowan, downloads anything, or runs a model. That distinction
matters: these prove the plumbing, NOT that any prediction is correct.

    ./env/bin/python scripts/test_pipeline.py
"""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# The project venv when it exists, otherwise whatever interpreter is running
# this file. CI has no env/ - it installs into the runner's own Python - so a
# hardcoded path made two of these ten tests fail on every push regardless of
# what changed.
_VENV = ROOT / "env" / "bin" / "python"
PY = str(_VENV) if _VENV.exists() else sys.executable
PIPE = ROOT / "results" / "pipeline" / "colorectal-cancer"
BOARD = ROOT / "results" / "repurpose_colorectal-cancer.json"


def board():
    return json.loads(BOARD.read_text())


# -- stage interface: score refuses a signature it does not have ------------
def test_score_refuses_unavailable_signature():
    """It used to load boltzgen_signature_v4.json (KDR's designs) unconditionally,
    so scoring ANY other target silently ranked it against KDR."""
    r = subprocess.run([PY, str(ROOT / "scripts/repurpose.py"), "score",
                        "--rank-by", "boltzgen_consensus"],
                       capture_output=True, text=True)
    assert r.returncode == 1, f"expected exit 1, got {r.returncode}"
    assert "not available" in (r.stdout + r.stderr), r.stdout + r.stderr


def test_score_runs_without_any_design_run():
    """The documented default path is the pocket signature; it must not need a
    BoltzGen artifact to exist at all."""
    # --out: a test must never overwrite the committed board. Without it this
    # test silently dropped boltzgen_consensus from the shipped artifact.
    with tempfile.TemporaryDirectory() as d:
        r = subprocess.run([PY, str(ROOT / "scripts/repurpose.py"), "score",
                            "--rank-by", "p2rank_geometry",
                            "--out", str(Path(d) / "board.json")],
                           capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-800:] + r.stderr[-800:]
    assert "p2rank_geometry" in r.stdout


# -- resume: re-running must not destroy work ------------------------------
def test_shortlist_append_is_idempotent():
    """Re-running the entry point used to REPLACE the shortlist, dropping the 12
    hard-decoy rows the reported AUC 0.717 rests on (42 -> 30)."""
    state = json.loads((PIPE / "repurpose_state.json").read_text())
    roles = {}
    for row in state["shortlist"]:
        roles[row["role"]] = roles.get(row["role"], 0) + 1
    assert roles.get("hard_decoy", 0) > 0, "fixture lost its hard-decoy arm"
    ids = [r["struct_id"] for r in state["shortlist"]]
    assert len(ids) == len(set(ids)), "shortlist already contains duplicates"


def test_entry_point_passes_append():
    """The orchestrator must never call shortlist destructively."""
    src = (ROOT / "scripts/autorepurpose.py").read_text()
    assert '"--append"' in src, "autorepurpose.py builds a shortlist without --append"


# -- failed candidates stay visible and unscored ---------------------------
def test_failed_candidate_is_visible_and_unscored():
    """A drug whose pose has no ligand must remain in the board with a reason and
    must NOT receive a plausible-looking score."""
    b = board()
    failed = [r for r in b["results"] if r.get("status") != "scored"]
    assert failed, "fixture has no failed candidate to check"
    for r in failed:
        assert r.get("status"), f"{r['name']} failed with no status"
        assert "rank" not in r, f"{r['name']} failed but carries a rank"
        for key in ("core_coverage", "weighted_jaccard", "precision_in_core"):
            assert key not in r, f"{r['name']} failed but carries {key}"
    assert b["n_scored"] == len([r for r in b["results"] if r["status"] == "scored"])
    assert b["n_scored"] < len(b["results"]), "n_scored should exclude failures"


# -- residue mapping -------------------------------------------------------
def test_sequence_offset_maps_to_the_right_residues():
    """cmd_score converts pose residue ids back to author numbering with a start
    offset. If that offset is wrong every contact is silently mislabelled, so it
    is checked against residues whose identity is known: KDR's hinge Cys919 and
    the DFG motif Asp1046/Phe1047."""
    fa = PIPE / "target" / "sequence_kinase_domain.fasta"
    lines = fa.read_text().splitlines()
    header, seq = lines[0], "".join(l.strip() for l in lines[1:])
    start = None
    for tok in header.replace("|", " ").split():
        if "-" in tok and tok.split("-")[0].isdigit():
            start = int(tok.split("-")[0])
            break
    assert start == 834, f"expected the kinase-domain slice to start at 834, got {start}"
    for auth, expect in ((919, "C"), (1046, "D"), (1047, "F"), (885, "E")):
        got = seq[auth - start]
        assert got == expect, f"auth {auth} maps to {got}, expected {expect}"


def test_engaged_residues_lie_inside_the_construct():
    b = board()
    lo, hi = json.loads((PIPE / "target" / "pocket.json").read_text())["construct_range"]
    for r in b["results"]:
        for rid in r.get("engaged_residues", []):
            assert lo <= rid <= hi, f"{r['name']} engages {rid}, outside {lo}-{hi}"


# -- score calculation -----------------------------------------------------
def test_precision_in_core_matches_its_definition():
    """precision_in_core = |engaged AND core| / |engaged|, recomputed from the
    stored residue lists rather than trusting the stored number."""
    b = board()
    core = set(b["signatures"]["p2rank_geometry"]["core_residues_auth"])
    checked = 0
    for r in b["results"]:
        if r.get("status") != "scored":
            continue
        eng = set(r["engaged_residues"])
        want = len(eng & core) / len(eng) if eng else 0.0
        got = r["by_signature"]["p2rank_geometry"]["precision_in_core"]
        assert abs(want - got) < 1e-9, f"{r['name']}: {got} != {want}"
        checked += 1
    assert checked > 20, f"only checked {checked} rows"


def test_ranking_is_consistent_with_the_ranked_by_signature():
    b = board()
    sig, metric = b["ranked_by"], b["ranking_metric"]
    scored = sorted((r for r in b["results"] if r["status"] == "scored"),
                    key=lambda r: r["rank"])
    vals = [r["by_signature"][sig][metric] for r in scored]
    assert vals == sorted(vals, reverse=True), "board is not ordered by its own metric"


def test_provenance_is_recorded():
    """A board must say what produced it."""
    b = board()
    assert b.get("versions", {}).get("rowan-python"), "no rowan version recorded"
    assert "pocket_constrained" in b.get("cofolding", {}), "no constraint flag recorded"
    assert b["cofolding"]["pocket_constrained"] is False, (
        "the shipped run was unconstrained; if this changes, the prose describing "
        "post-hoc scoring must change with it")


TESTS = [
    ("score refuses a signature it does not have", test_score_refuses_unavailable_signature),
    ("score runs with no design run at all", test_score_runs_without_any_design_run),
    ("shortlist has no duplicates and keeps its hard-decoy arm", test_shortlist_append_is_idempotent),
    ("the entry point calls shortlist with --append", test_entry_point_passes_append),
    ("a failed candidate stays visible and unscored", test_failed_candidate_is_visible_and_unscored),
    ("the sequence offset maps to the right residues", test_sequence_offset_maps_to_the_right_residues),
    ("engaged residues lie inside the construct", test_engaged_residues_lie_inside_the_construct),
    ("precision_in_core matches its definition", test_precision_in_core_matches_its_definition),
    ("the board is ordered by its own metric", test_ranking_is_consistent_with_the_ranked_by_signature),
    ("provenance is recorded on the board", test_provenance_is_recorded),
]


def main():
    print(f"scripts/test_pipeline.py - {len(TESTS)} tests "
          f"(offline, against committed artifacts; no model inference)\n")
    failed = []
    for name, fn in TESTS:
        try:
            fn()
            print(f"  ok   {name}")
        except Exception as exc:                              # noqa: BLE001
            failed.append((name, exc))
            print(f"  FAIL {name}\n       {exc}")
    print(f"\n{len(TESTS) - len(failed)} passed, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
