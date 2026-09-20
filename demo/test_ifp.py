"""Check that the interaction fingerprints mean something, not just that they run.

A fingerprint that compiles is worthless; a fingerprint that quietly fires the
wrong bits is worse than worthless, because everything downstream still looks
like a result. So the load-bearing tests here are against facts that were
established before this code existed:

  * 1PPB is thrombin with PPACK (D-Phe-Pro-Arg chloromethylketone) bound in the
    active site, solved by Bode et al. (EMBO J 1989, 8:3467). That paper says
    what PPACK touches: the P1 arginine buries in the S1 pocket and salt-bridges
    Asp189; the ligand backbone pairs with Gly216 in an antiparallel beta-sheet;
    the carbonyl sits in the Gly193/Ser195 oxyanion hole; the D-Phe occupies the
    aryl-binding site against Trp215, Ile174, Leu99 and the 60-loop insertion
    (Tyr60A, Trp60D). If the fingerprint does not reproduce that, it is wrong.

  * A fingerprint is supposed to track where the ligand actually is. Slide the
    real pose out of the pocket and the similarity to its own starting point has
    to fall away. That is asserted as a rank correlation with its p-value over a
    displacement sweep, not as a threshold someone chose after seeing the number.

  * The three moieties of PPACK sit in three named subsites. S2 (proline, under
    the 60-loop) and S3 (D-Phe, in the aryl site) share the shallow upper cleft;
    S1 is the deep specificity pocket below them. The prediction made before
    running it was that S2 and S3 score more alike than either does against S1,
    and it is asserted for all seven similarity measures rather than the one
    that happened to work.

The head-to-head against ECFP4 that this module exists for CANNOT be run: the
batch-docking workflow was submitted with num_poses_to_save = 0, so the 203 Vina
poses were scored and thrown away. test_poses_are_missing_upstream verifies that
rather than assuming it. See docs/14-INTERACTION-FINGERPRINTS.md.

Run:  ./env-kit/bin/python -m demo.test_ifp
No credits. Network only for the optional live Rowan check, which is read-only.
"""

from __future__ import annotations

import json
import os

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr

from demo import ifp
from demo.contacts import THREE_TO_ONE, residue_order

PDB = ifp.REPO / "results/demo/_cache/1PPB.pdb"
DOCKED = ifp.REPO / "results/demo/thrombin/08_docked.json"

# 1PPB ligand 0G6 is PPACK, covalently trapped at the catalytic serine. The PDB
# gives coordinates and element symbols but no bond orders, so the chemistry has
# to be supplied; these 30 heavy atoms are D-Phe-Pro-Arg plus the ketone methyl.
PPACK_SMILES = "NC(Cc1ccccc1)C(=O)N1CCCC1C(=O)NC(CCCNC(N)=[NH2+])C(=O)C"
PPACK_RESNAME = "0G6"
THROMBIN_CHAIN = "H"  # the catalytic heavy chain, chymotrypsin numbering

# PPACK atoms by PDB name, split into the three subsite-occupying moieties.
MOIETIES = {
    "S1_P1_arginine": ["CA2", "CB2", "CG2", "CD3", "NE", "CZ1", "NH1", "NH2"],
    "S2_P2_proline": ["N1", "CA1", "CB1", "CG1", "CD", "C1", "O1"],
    "S3_P3_D_phenylalanine": ["N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
}


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


def load_reference():
    """The real PPACK pose, thrombin, and a site axis defined from the structure."""
    if not PDB.exists():
        raise SystemExit(
            f"{PDB} is missing. It is fetched by demo/selftest.py; this test will "
            "not invent a structure to stand in for it."
        )
    lines = [
        line for line in PDB.read_text().splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == PPACK_RESNAME
    ]
    if not lines:
        raise SystemExit(f"{PDB.name}: no {PPACK_RESNAME} ligand records found")
    pose = ifp.pose_from_pdb_block("\n".join(lines) + "\nEND\n", PPACK_SMILES, "PPACK")
    receptor = ifp.parse_receptor(PDB)

    # The bit space is every chain-H residue with an atom within 8 A of the
    # ligand. Defining it from the structure rather than a hand-typed list keeps
    # the insertion-coded 60-loop in, which is the part that matters most.
    residues = sorted(
        (
            label for (chain, label), entry in receptor.residues.items()
            if chain == THROMBIN_CHAIN
            and ifp._cdist(entry["xyz"], pose.xyz).min() <= 8.0
        ),
        key=residue_order,
    )
    axis = ifp.SiteAxis(chain=THROMBIN_CHAIN, residues=residues)
    return pose, receptor, axis


def moiety_fingerprints(pose, receptor, axis) -> dict[str, ifp.Fingerprint]:
    names = [a.GetPDBResidueInfo().GetName().strip() for a in pose.mol.GetAtoms()]
    out = {}
    for label, wanted in MOIETIES.items():
        ids = [i for i, n in enumerate(names) if n in wanted]
        if len(ids) != len(wanted):
            raise SystemExit(
                f"{label}: expected {len(wanted)} atoms, matched {len(ids)}. The "
                "ligand in the file is not the one these atom names describe."
            )
        out[label] = ifp.fingerprint(ifp.subpose(pose, ids, label), receptor, axis)
    return out


def translate(pose: ifp.Pose, vector: np.ndarray) -> ifp.Pose:
    moved = Chem.Mol(pose.mol)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + vector[0], p.y + vector[1], p.z + vector[2]))
    return ifp.Pose(mol=moved, name=f"{pose.name}+{np.linalg.norm(vector):.1f}A")


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------


def test_measure_formulae() -> bool:
    """The seven similarity formulae, against values worked out by hand.

    A transcription slip in one of these would not crash anything; it would just
    reorder the results table. So they are pinned to an example whose a, b, c, d
    can be counted by eye.
    """
    import math

    x = np.array([1, 1, 1, 0, 0, 0, 0, 0], dtype=bool)
    y = np.array([1, 1, 0, 1, 0, 0, 0, 0], dtype=bool)
    a, b, c, d = ifp.counts(x, y)
    ok = (a, b, c, d) == (2, 1, 1, 4)
    print(f"  counts a={a} b={b} c={c} d={d} (expected 2,1,1,4) {'OK' if ok else 'FAIL'}")

    expected = {
        "tanimoto": 2 / 4,
        "sokal_michener": 6 / 8,
        "rogers_tanimoto": 6 / 10,
        "sokal_sneath_2": 12 / 14,
        "consonni_todeschini_1": math.log(7) / math.log(9),
        "consonni_todeschini_2": (math.log(9) - math.log(3)) / math.log(9),
        "austin_colwell": (2 / math.pi) * math.asin(math.sqrt(6 / 8)),
    }
    for name, want in expected.items():
        got = ifp.MEASURES[name](x, y)
        good = abs(got - want) < 1e-12
        ok &= good
        print(f"  {name:24s} {got:.6f} (expected {want:.6f}) {'OK' if good else 'FAIL'}")

    # Symmetry, and the identity of indiscernibles, on the same vectors.
    for name, fn in ifp.MEASURES.items():
        good = abs(fn(x, y) - fn(y, x)) < 1e-12 and abs(fn(x, x) - 1.0) < 1e-12
        ok &= good
        if not good:
            print(f"  {name}: not symmetric or self-similarity != 1  FAIL")
    print(f"  symmetry and self-similarity over all 7 measures {'OK' if ok else 'FAIL'}")
    return bool(ok)


def test_identity(pose, receptor, axis) -> bool:
    """A real pose against itself is 1.0 on every measure."""
    fp = ifp.fingerprint(pose, receptor, axis)
    sims = ifp.all_similarities(fp, fp)
    ok = all(abs(v - 1.0) < 1e-12 for v in sims.values())
    print(f"  PPACK vs itself: {min(sims.values()):.6f} .. {max(sims.values()):.6f} "
          f"over {len(sims)} measures  {'OK' if ok else 'FAIL'}")
    print(f"  axis: {len(axis.residues)} residues x {len(axis.types)} types = "
          f"{len(axis)} bits; {int(fp.bits.sum())} on")
    return ok


def test_reproduces_published_1ppb_contacts(pose, receptor, axis) -> bool:
    """Every typed channel is checked against Bode et al. 1989, EMBO J 8:3467."""
    fp = ifp.fingerprint(pose, receptor, axis)
    found = fp.by_type()

    # (channel, residue, what the literature says it is)
    claims = [
        ("ionic_cation", 189, "P1 Arg guanidinium salt bridge to Asp189, S1 pocket"),
        ("hbond_donor", 216, "antiparallel beta-sheet pairing with Gly216"),
        ("hbond_acceptor", 193, "carbonyl in the Gly193 oxyanion hole"),
        ("aromatic", 215, "D-Phe ring against Trp215 in the aryl-binding site"),
        ("hydrophobic", "60A", "Tyr60A, the 60-loop lid over S2"),
        ("hydrophobic", "60D", "Trp60D, the 60-loop lid over S2"),
        ("hydrophobic", 99, "Leu99, floor of the aryl-binding site"),
        ("hydrophobic", 174, "Ile174, wall of the aryl-binding site"),
    ]
    ok = True
    for channel, residue, why in claims:
        hit = residue in found[channel]
        ok &= hit
        print(f"  {channel:15s} {str(residue):>4s}  {'OK  ' if hit else 'MISS'}  {why}")

    # Asp189 should be the only ionic contact: thrombin's S1 has exactly one
    # acidic residue at its base, and a second ionic bit would mean the charged
    # group was being matched to something it does not reach.
    only = found["ionic_cation"] == [189] and found["ionic_anion"] == []
    ok &= only
    print(f"  ionic channel is Asp189 and nothing else: "
          f"{found['ionic_cation']} / {found['ionic_anion']}  {'OK' if only else 'FAIL'}")

    for channel, residues in found.items():
        print(f"    {channel:15s} {residues}")
    return bool(ok)


def test_subsites_separate(pose, receptor, axis) -> bool:
    """Moieties in adjacent subsites resemble each other more than distant ones.

    Prediction registered before the numbers were seen: the P2 proline and the
    P3 D-Phe both sit in the shallow upper cleft, while the P1 arginine is buried
    in the deep S1 pocket below, so sim(S2, S3) should exceed both sim(S1, S2)
    and sim(S1, S3). Asserted for all seven measures.
    """
    fps = moiety_fingerprints(pose, receptor, axis)
    s1, s2, s3 = (fps[k] for k in MOIETIES)
    for label, fp in fps.items():
        print(f"  {label:24s} {int(fp.bits.sum())} bits  "
              f"{ {k: v for k, v in fp.by_type().items() if v} }")

    ok = True
    for name, fn in ifp.MEASURES.items():
        near = fn(s2.bits, s3.bits)
        far_a = fn(s1.bits, s2.bits)
        far_b = fn(s1.bits, s3.bits)
        good = near > far_a and near > far_b
        ok &= good
        print(f"  {name:24s} S2-S3={near:.4f}  S1-S2={far_a:.4f}  S1-S3={far_b:.4f}"
              f"  {'OK' if good else 'FAIL'}")

    # And the deep pocket really is a different place: the S1 moiety reaches
    # Asp189, the upper-cleft moieties do not.
    disjoint = 189 in s1.by_type()["ionic_cation"] and not any(
        189 in fp.by_type()["ionic_cation"] for fp in (s2, s3)
    )
    ok &= disjoint
    print(f"  only the P1 arginine reaches Asp189  {'OK' if disjoint else 'FAIL'}")
    return bool(ok)


def test_similarity_decays_with_displacement(pose, receptor, axis) -> bool:
    """Slide the pose out of the pocket; similarity to its own origin must fall.

    Six axis-aligned directions x 0-10 A in 1 A steps, n = 66 fingerprints.
    Reported as Spearman rho with its p-value: the assertion is the direction of
    the correlation, which was fixed in advance, not a cut-off on its size.
    """
    reference = ifp.fingerprint(pose, receptor, axis)
    directions = list(np.eye(3)) + list(-np.eye(3))
    distances = np.arange(0.0, 10.5, 1.0)

    d_col, s_col = [], []
    for direction in directions:
        for d in distances:
            fp = ifp.fingerprint(translate(pose, direction * d), receptor, axis)
            d_col.append(d)
            s_col.append(ifp.tanimoto(reference.bits, fp.bits))

    rho, p = spearmanr(d_col, s_col)
    ok = rho < 0 and p < 0.01
    print(f"  Spearman rho = {rho:.3f}, p = {p:.2e}, n = {len(d_col)} fingerprints"
          f"  {'OK' if ok else 'FAIL'}")

    at_zero = [s for d, s in zip(d_col, s_col) if d == 0.0]
    exact = all(abs(s - 1.0) < 1e-12 for s in at_zero)
    ok &= exact
    print(f"  zero displacement reproduces the pose exactly ({len(at_zero)} runs)"
          f"  {'OK' if exact else 'FAIL'}")

    far = [s for d, s in zip(d_col, s_col) if d == 10.0]
    fell = max(far) < min(at_zero)
    ok &= fell
    print(f"  at 10 A displacement similarity is {min(far):.3f}-{max(far):.3f}"
          f"  {'OK' if fell else 'FAIL'}")
    return bool(ok)


def test_axis_is_enforced(receptor) -> bool:
    """A site residue the structure does not contain must raise, not zero-fill."""
    ok = True
    try:
        ifp.SiteAxis(chain=THROMBIN_CHAIN, residues=[189, 99999])
        receptor.require(THROMBIN_CHAIN, [189, 99999])
        ok = False
        print("  absent residue accepted silently  FAIL")
    except ValueError as exc:
        print(f"  absent site residue refused: {str(exc)[:70]}...  OK")
    try:
        ifp.counts(np.zeros(10, bool), np.zeros(12, bool))
        ok = False
        print("  mismatched fingerprint lengths compared  FAIL")
    except ValueError:
        print("  fingerprints on different axes refused  OK")
    return ok


def test_head_to_head_machinery() -> bool:
    """Code check on synthetic scores. These are NOT results about thrombin.

    The real comparison needs docked poses, which do not exist yet. What can be
    checked now is that the statistics behave: a perfect ranker scores AUC 1.0
    with a permutation p at the floor, a coin flip lands near 0.5 and does not
    reach significance, and the paired bootstrap on the difference separates the
    two while a ranker compared with itself gives exactly zero difference.
    """
    rng = np.random.default_rng(11)
    n, n_actives = 203, 3
    actives = np.zeros(n, dtype=bool)
    actives[rng.choice(n, n_actives, replace=False)] = True

    perfect = np.where(actives, 1.0, 0.0) + rng.normal(0, 0.01, n)
    coin = rng.normal(0, 1, n)

    ok = True
    a_perfect = ifp.auc(perfect, actives)
    a_coin = ifp.auc(coin, actives)
    p_perfect = ifp.permutation_p(perfect, actives, n_perm=2000, seed=1)
    p_coin = ifp.permutation_p(coin, actives, n_perm=2000, seed=1)
    for label, value, cond in [
        (f"perfect ranker AUC = {a_perfect:.3f}", a_perfect, a_perfect == 1.0),
        (f"perfect ranker permutation p = {p_perfect:.2e}", p_perfect, p_perfect < 0.01),
        (f"coin-flip ranker AUC = {a_coin:.3f}", a_coin, 0.0 <= a_coin <= 1.0),
        (f"coin-flip ranker p = {p_coin:.3f}", p_coin, p_coin > 0.01),
    ]:
        ok &= cond
        print(f"  {label}  {'OK' if cond else 'FAIL'}")

    delta = ifp.paired_bootstrap_delta_auc(perfect, coin, actives, n_boot=2000, seed=2)
    good = delta["delta_auc"] > 0 and delta["ci95"][0] > 0
    ok &= good
    print(f"  delta AUC = {delta['delta_auc']:.3f}, 95% CI "
          f"[{delta['ci95'][0]:.3f}, {delta['ci95'][1]:.3f}], "
          f"n_actives = {delta['n_actives']}  {'OK' if good else 'FAIL'}")

    same = ifp.paired_bootstrap_delta_auc(coin, coin, actives, n_boot=500, seed=3)
    good = same["delta_auc"] == 0.0 and same["ci95"] == (0.0, 0.0)
    ok &= good
    print(f"  a ranker against itself: delta = {same['delta_auc']:.3f}, CI "
          f"{same['ci95']}  {'OK' if good else 'FAIL'}")
    return bool(ok)


def test_poses_are_missing_upstream() -> bool:
    """Verify, rather than assume, that there are no poses to fingerprint.

    The offline half reads the docking artefact. The online half asks Rowan and
    is skipped without a key; it is a GET and costs nothing.
    """
    ok = True
    docked = json.loads(DOCKED.read_text())
    pose_like = [k for k in docked if any(
        token in k.lower() for token in ("pose", "structure", "sdf", "pdb", "coord")
    )]
    good = not pose_like
    ok &= good
    print(f"  08_docked.json keys: {sorted(docked)}")
    print(f"  no pose-bearing key among them: {pose_like or 'none'}  "
          f"{'OK' if good else 'FAIL'}")
    print(f"  it holds {len(docked['smiles'])} SMILES and {len(docked['scores'])} "
          "scores, and nothing three-dimensional")

    ifp.load_env_if_available()
    if not os.environ.get("ROWAN_API_KEY"):
        print("  [skipped] live Rowan check: ROWAN_API_KEY not set")
        return ok

    status = ifp.check_saved_poses()
    good = status["num_poses_to_save"] == 0 and status["n_descendant_structures"] == 0
    ok &= good
    print(f"  Rowan workflow {status['workflow_uuid']}: "
          f"num_poses_to_save={status['num_poses_to_save']}, "
          f"max_poses={status['max_poses_per_ligand']}, "
          f"descendant structures={status['n_descendant_structures']}, "
          f"ligands={status['n_ligands']}  {'OK' if good else 'FAIL'}")
    print(f"  poses_available = {status['poses_available']}")
    return bool(ok)


def test_ligand_set() -> bool:
    """The enlarged docking set is repo-sourced, parses, and is labelled by rule."""
    import csv

    ligands = ifp.thrombin_ligand_set()
    ok = True

    n, n_act = len(ligands["smiles"]), sum(ligands["actives"])
    good = n == 216 and n_act == 14 and ligands["n_generated"] == 200
    ok &= good
    print(f"  {n} ligands = {ligands['n_generated']} generated decoys + "
          f"{len(ligands['drugs'])} measured F2 binders; {n_act} actives"
          f"  {'OK' if good else 'FAIL'}")

    # Every SMILES parses. A SMILES that does not is a molecule we do not have.
    bad = [nm for nm, s in zip(ligands["names"], ligands["smiles"])
           if Chem.MolFromSmiles(s) is None]
    ok &= not bad
    print(f"  all {n} SMILES parse: {bad or 'yes'}  {'OK' if not bad else 'FAIL'}")

    # The actives are not hand-picked: the label is the repo's own flag.
    flags = {
        row["name"]: row["passes_threshold"] == "1"
        for row in csv.DictReader((ifp.REPO / "data/target_annotations.csv").open())
        if row["uniprot"] == "P00734" and row["evidence"] == "chembl_pchembl"
    }
    mismatch = [r["name"] for r in ligands["drugs"] if flags[r["name"]] != r["active"]]
    ok &= not mismatch
    print(f"  labels follow data/target_annotations.csv passes_threshold: "
          f"{mismatch or 'all agree'}  {'OK' if not mismatch else 'FAIL'}")

    # The three original controls must carry the exact SMILES they were docked
    # with, or the new run is not against the same background.
    docked = json.loads(DOCKED.read_text())
    controls = dict(zip(docked["controls"], docked["smiles"][200:203]))
    same = all(
        ligands["smiles"][ligands["names"].index(name)] == smiles
        for name, smiles in controls.items()
    )
    ok &= same
    print(f"  original controls keep their exact docked SMILES  "
          f"{'OK' if same else 'FAIL'}")

    # And no generated decoy is secretly one of the drugs.
    canonical = {Chem.MolToSmiles(Chem.MolFromSmiles(s))
                 for s in ligands["smiles"][:ligands["n_generated"]]}
    dupes = [r["name"] for r in ligands["drugs"]
             if Chem.MolToSmiles(Chem.MolFromSmiles(r["smiles"])) in canonical]
    ok &= not dupes
    print(f"  no drug duplicates a generated decoy: {dupes or 'none'}  "
          f"{'OK' if not dupes else 'FAIL'}")

    query_is_active = ligands["actives"][ligands["names"].index(ligands["query"])]
    ok &= query_is_active
    print(f"  query {ligands['query']!r} is in the set and is an active  "
          f"{'OK' if query_is_active else 'FAIL'}")
    return bool(ok)


def test_axes_are_declared_and_differ() -> bool:
    """Both fingerprint axes validate, and the difference between them is real.

    02_site.json's site was built as a 6 A shell around the catalytic Ser195.
    Asserting that it omits Asp189 is not a complaint -- it is the reason the
    box-lining axis is the primary one, and it is checked rather than claimed.
    """
    receptor = ifp.parse_receptor(
        ifp.REPO / "results/demo/thrombin/structures/target_thrombin.pdb"
    )
    site = json.loads((ifp.REPO / "results/demo/thrombin/02_site.json").read_text())
    box = json.loads(DOCKED.read_text())["box"]
    chain = site["target_chain"]

    # The GDSGGP motif fixes the offset from this structure's numbering to
    # chymotrypsin numbering, so Asp189 can be located rather than assumed.
    motif = "".join(
        THREE_TO_ONE[receptor.residues[(chain, i)]["resname"]]
        for i in range(203, 209)
    )
    ok = motif == "GDSGGP"
    print(f"  structure 203-208 spells {motif} (G193-D194-S195-G196-G197-P198) "
          f"-> Asp189 is position 199  {'OK' if ok else 'FAIL'}")
    is_asp = receptor.residues[(chain, 199)]["resname"] == "ASP"
    ok &= is_asp
    print(f"  position 199 is {receptor.residues[(chain, 199)]['resname']}  "
          f"{'OK' if is_asp else 'FAIL'}")

    box_axis = ifp.SiteAxis(chain=chain,
                            residues=ifp.box_lining_residues(receptor, chain, box))
    site_axis = ifp.SiteAxis(
        chain=chain,
        residues=[int(x) for x in site["chosen"]["residue_positions_1based"]],
    )
    receptor.require(chain, box_axis.residues)
    receptor.require(chain, site_axis.residues)

    good = 199 in box_axis.residues and 199 not in site_axis.residues
    ok &= good
    print(f"  box_lining: {len(box_axis.residues)} residues, {len(box_axis)} bits, "
          f"Asp189 present = {199 in box_axis.residues}")
    print(f"  site_02   : {len(site_axis.residues)} residues, {len(site_axis)} bits, "
          f"Asp189 present = {199 in site_axis.residues}  {'OK' if good else 'FAIL'}")

    nested = set(site_axis.residues) <= set(box_axis.residues)
    ok &= nested
    print(f"  the pipeline site is a subset of the box-lining set  "
          f"{'OK' if nested else 'FAIL'}")
    return bool(ok)


def main() -> int:
    pose, receptor, axis = load_reference()

    checks = [
        ("similarity formulae", lambda: test_measure_formulae()),
        ("identity", lambda: test_identity(pose, receptor, axis)),
        ("published 1PPB contacts", lambda: test_reproduces_published_1ppb_contacts(pose, receptor, axis)),
        ("subsite separation", lambda: test_subsites_separate(pose, receptor, axis)),
        ("decay with displacement", lambda: test_similarity_decays_with_displacement(pose, receptor, axis)),
        ("axis is enforced", lambda: test_axis_is_enforced(receptor)),
        ("head-to-head machinery (synthetic)", lambda: test_head_to_head_machinery()),
        ("poses are missing upstream", lambda: test_poses_are_missing_upstream()),
        ("enlarged ligand set", lambda: test_ligand_set()),
        ("fingerprint axes", lambda: test_axes_are_declared_and_differ()),
    ]

    results = []
    for name, fn in checks:
        print(f"\n== {name} ==")
        results.append((name, fn()))

    print("\n== summary ==")
    for name, passed in results:
        print(f"  {'PASS' if passed else 'FAIL'}  {name}")
    failed = [n for n, p in results if not p]
    if failed:
        print(f"\n{len(failed)} check(s) failed: {failed}")
        return 1
    print("\nAll checks passed. The fingerprints are meaningful; there is still "
          "nothing to compare them on.\nSee docs/14-INTERACTION-FINGERPRINTS.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
