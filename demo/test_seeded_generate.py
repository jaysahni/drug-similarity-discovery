"""Check that seeded generation does what it claims -- not merely that it runs.

The failure mode this guards against is a green pipeline that produces nothing
thrombin-shaped: BRICS recombination will happily emit molecules assembled
entirely from the bland end of the fragment pool, and a test that only asserts
"the code returned a list" would pass on those. So the load-bearing checks here
are about chemistry, not control flow:

  * every kept molecule parses, obeys the gate demo/generate.py imposes, and is
    NOT an approved drug -- re-derived here from InChIKeys rather than trusted
    from the generator's own bookkeeping;
  * the arginine mimetic that thrombin's S1 pocket requires is actually present,
    and present far more often than in the unseeded model's output (Fisher, with
    its p);
  * the circularity meter is a real measurement, which means it must read ~1.0
    on a BRICS product AND well below 1.0 on a molecule the seeds had no hand in.
    A meter that says 1.0 for everything measures nothing.
  * the seeded-vs-unseeded comparison is calibrated: handed two samples from the
    same distribution it must NOT report a difference.

The real seeded-vs-unseeded numbers are printed but deliberately not asserted in
a direction. If seeding turns out not to help, that is a finding and the test
suite is not the place to suppress it.

Run:  ./env-kit/bin/python -m demo.test_seeded_generate
No network, no credentials, no credits. Takes ~1-2 minutes (it generates).
"""

from __future__ import annotations

import numpy as np

from demo import generate as unseeded_arm
from demo import seeded_generate as sg

N_TEST = 60          # small: this is a self-test, not the run
BUDGET_S = 45.0


def main() -> int:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED

    RDLogger.DisableLog("rdApp.*")
    checks: list[tuple[str, bool, str]] = []

    approved = unseeded_arm.load_approved()
    seeds = sg.seed_drugs(approved, sg.SYMBOL)

    # 1. The seeds are derived, not asserted from memory.
    names = sorted(s["name"] for s in seeds)
    checks.append(("seeds are the F2-mechanism drugs", len(seeds) == 4, ", ".join(names)))
    checks.append(("every seed parses", all(Chem.MolFromSmiles(s["smiles"]) for s in seeds),
                   f"{len(seeds)} seeds"))

    # 2. The fragment pool carries the pharmacophore the S1 pocket needs.
    pool = sg.fragment_pool(seeds)
    motif_frags = [f for f in pool if sg.motif_hits(f.replace("*", "H"))] or [
        f for f in pool
        if any(Chem.MolFromSmiles(f).HasSubstructMatch(Chem.MolFromSmarts(p))
               for p in sg.ARG_MIMETIC.values())]
    checks.append(("fragment pool non-empty", len(pool) > 0, f"{len(pool)} fragments"))
    checks.append(("pool contains an arginine mimetic", len(motif_frags) > 0,
                   f"{len(motif_frags)}: {motif_frags[:3]}"))

    # 3. Generate.
    gen = sg.generate_seeded(seeds, approved, n=N_TEST, seed=sg.SEED, budget_s=BUDGET_S)
    mols = gen["molecules"]
    checks.append(("generation produced molecules", len(mols) > 0,
                   f'{len(mols)} kept of {gen["build"]["n_raw"]} raw '
                   f'({gen["build"]["stop_reason"]})'))
    if not mols:
        return _report(checks)

    # 4. Valid, and inside the gate -- recomputed, not taken on trust.
    parsed = [Chem.MolFromSmiles(m["smiles"]) for m in mols]
    checks.append(("all valid SMILES", all(p is not None for p in parsed),
                   f"{sum(p is not None for p in parsed)}/{len(parsed)}"))
    in_gate = [
        unseeded_arm.MW_RANGE[0] <= Descriptors.MolWt(p) <= unseeded_arm.MW_RANGE[1]
        and QED.qed(p) >= unseeded_arm.MIN_QED
        for p in parsed if p is not None]
    checks.append(("all inside demo/generate.py's gate", all(in_gate),
                   f"{sum(in_gate)}/{len(in_gate)} MW {unseeded_arm.MW_RANGE} QED>="
                   f"{unseeded_arm.MIN_QED}"))
    canonical = [Chem.MolToSmiles(p) for p in parsed if p is not None]
    checks.append(("deduplicated by canonical SMILES", len(set(canonical)) == len(canonical),
                   f"{len(set(canonical))} unique of {len(canonical)}"))

    # 5. Novel: none of them IS an approved drug. Re-derived from InChIKeys here.
    approved_keys = {r["inchikey"] for r in approved if r["inchikey"]}
    checks.append(("approved InChIKeys were actually loaded", len(approved_keys) > 2000,
                   f"{len(approved_keys)} keys"))
    collisions = [Chem.MolToSmiles(p) for p in parsed
                  if p is not None and Chem.MolToInchiKey(p) in approved_keys]
    checks.append(("none is an approved drug", not collisions,
                   f"{len(collisions)} collisions"))
    #    ... and the check has teeth: an approved drug fed in must be caught.
    sentinel = Chem.MolToInchiKey(Chem.MolFromSmiles(seeds[0]["smiles"])) in approved_keys
    checks.append(("novelty check has teeth (seed drug is caught)", sentinel,
                   f'{seeds[0]["name"]} recognised as approved'))

    # 6. The seeded motif is present, and present far more than in unseeded output.
    from scipy.stats import fisher_exact

    seeded_motif = sum(bool(m["motifs"]) for m in mols)
    unseeded_smiles = sg.load_unseeded()
    unseeded_motif = sum(bool(sg.motif_hits(s)) for s in unseeded_smiles)
    odds, p_motif = fisher_exact(
        [[seeded_motif, len(mols) - seeded_motif],
         [unseeded_motif, len(unseeded_smiles) - unseeded_motif]], alternative="greater")
    #    A quarter, not a majority: BRICSBuild walks depth-first, so the amidine-rich
    #    products come out early and the motif fraction FALLS as the run gets longer
    #    (measured: 78.6% at n=60, 38.6% at n=600). A "majority" threshold would pass
    #    on the small self-test and fail on the real run, which is worse than useless.
    checks.append(("a substantial fraction carry an arginine mimetic",
                   seeded_motif > len(mols) / 4, f"{seeded_motif}/{len(mols)}"))
    checks.append(("motif enriched over unseeded (Fisher)", p_motif < 0.001,
                   f"{seeded_motif}/{len(mols)} vs {unseeded_motif}/"
                   f"{len(unseeded_smiles)}, p={p_motif:.3g}"))

    # 6b. The LibInvent scaffolds are derived from the seeds, not invented. Each one
    #     must be a substructure of a seed drug; this needs no REINVENT4 install.
    provenance = sg.verify_scaffolds(sg.LIBINVENT_SCAFFOLDS, seeds)
    orphans = [k for k, v in provenance.items() if not v]
    checks.append(("every LibInvent scaffold comes from a seed drug", not orphans,
                   "; ".join(f"{k}<-{','.join(v)}" for k, v in provenance.items())))
    #     ... and the provenance check is not vacuous: a scaffold from nowhere is caught.
    bogus = sg.verify_scaffolds({"bogus": "[*:0]c1ccncc1[*:1]S(=O)(=O)C"}, seeds)
    checks.append(("scaffold provenance check has teeth", not bogus["bogus"],
                   f'bogus -> {bogus["bogus"]}'))
    #     ... and every scaffold carries the pharmacophore it is supposed to pin, so
    #     that a LibInvent product cannot come back without it. Checked on the
    #     methyl-capped scaffold, because SMARTS-against-SMARTS matching is not
    #     reliable in RDKit and would quietly pass on anything.
    import re as _re
    scaffold_motif = [
        k for k, s in sg.LIBINVENT_SCAFFOLDS.items()
        if not sg.motif_hits(_re.sub(r"\[\*:\d+\]", "C", s))]
    checks.append(("every LibInvent scaffold pins an arginine mimetic", not scaffold_motif,
                   f"without: {scaffold_motif}" if scaffold_motif
                   else f"all {len(sg.LIBINVENT_SCAFFOLDS)}"))

    # 7. The circularity meter measures something rather than everything. It must
    #    read high on BRICS products, clearly lower on molecules the seeds had no
    #    hand in -- and its size floor must be doing real work, because without one
    #    it saturates at 1.000 for every molecule ever tested.
    cores = sg.fragment_cores(list(pool))
    frac = np.array([m["seed_atom_fraction"] for m in mols], dtype=float)
    background = np.array([sg.seed_atom_fraction(s, cores) for s in unseeded_smiles[:200]],
                          dtype=float)
    checks.append(("seed_atom_fraction separates seeded from unseeded",
                   float(np.median(frac)) > float(np.median(background)) + 0.2,
                   f"seeded={np.median(frac):.3f} unseeded={np.median(background):.3f}"))
    foreign = "CC(=O)Oc1ccccc1C(=O)O"          # aspirin: not made of thrombin fragments
    checks.append(("seed_atom_fraction < 1 on a foreign molecule",
                   sg.seed_atom_fraction(foreign, cores) < 0.9,
                   f"aspirin={sg.seed_atom_fraction(foreign, cores):.3f}"))
    no_floor = sg.fragment_cores(list(pool), min_atoms=1)
    saturated = [sg.seed_atom_fraction(s, no_floor)
                 for s in [foreign, mols[0]["smiles"], unseeded_smiles[0]]]
    checks.append(("the size floor is load-bearing (no floor -> saturates at 1.0)",
                   all(v > 0.999 for v in saturated),
                   f"aspirin/seeded/unseeded = {[round(v, 3) for v in saturated]}"))

    # 8. The comparison is calibrated: same distribution in, no difference out.
    rng = np.random.default_rng(0)
    fake_a = [{"tanimoto": float(v), "nn_is_target_annotated": False}
              for v in rng.random(300)]
    fake_b = [{"tanimoto": float(v), "nn_is_target_annotated": False}
              for v in rng.random(300)]
    null_test = sg.compare(fake_a, fake_b, label="calibration")
    checks.append(("compare() finds nothing between identical distributions",
                   null_test["tanimoto"]["p_seeded_greater"] > 0.05,
                   f'p={null_test["tanimoto"]["p_seeded_greater"]} '
                   f'rb={null_test["tanimoto"]["rank_biserial"]:+.3f}'))
    shifted = [{"tanimoto": r["tanimoto"] + 0.3, "nn_is_target_annotated": False}
               for r in fake_a]
    shift_test = sg.compare(shifted, fake_b, label="calibration, shifted")
    checks.append(("compare() finds a real shift",
                   shift_test["tanimoto"]["p_seeded_greater"] < 1e-6,
                   f'p={shift_test["tanimoto"]["p_seeded_greater"]}'))

    # 9. The real comparison. Reported, not asserted in a direction.
    seeded_rows = sg.nearest_neighbours([m["smiles"] for m in mols], approved, sg.SYMBOL)
    unseeded_rows = sg.nearest_neighbours(unseeded_smiles, approved, sg.SYMBOL)
    real = sg.compare(seeded_rows, unseeded_rows, label="seeded vs unseeded")
    a = sg.summarise(seeded_rows, "seeded")
    b = sg.summarise(unseeded_rows, "unseeded")
    checks.append(("both sides matched against the same corpus",
                   real["n_seeded"] == len(mols) and real["n_reference"] == len(unseeded_smiles),
                   f'n={real["n_seeded"]} vs {real["n_reference"]}'))
    print(f'\n  seeded   n={a["n"]:4d} best={a["best_nn_tanimoto"]} median={a["median_nn_tanimoto"]}'
          f' annotated-NN={a["n_nn_target_annotated"]} top10={a["n_target_annotated_in_top10"]}')
    print(f'  unseeded n={b["n"]:4d} best={b["best_nn_tanimoto"]} median={b["median_nn_tanimoto"]}'
          f' annotated-NN={b["n_nn_target_annotated"]} top10={b["n_target_annotated_in_top10"]}')
    print(f'  Mann-Whitney p={real["tanimoto"]["p_seeded_greater"]} '
          f'rank-biserial={real["tanimoto"]["rank_biserial"]:+.3f}; '
          f'annotated-NN Fisher p={real["nn_target_annotated"]["p_seeded_greater"]}')
    print("  (direction deliberately not asserted -- a null result here is a finding)\n")

    return _report(checks)


def _report(checks: list[tuple[str, bool, str]]) -> int:
    failed = 0
    for name, ok, detail in checks:
        print(f'  {"PASS" if ok else "FAIL"}  {name:46s} {detail}')
        failed += not ok
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
