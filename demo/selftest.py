"""Check demo/contacts.py against the heavy pipeline's ground truth.

contacts.py is the one piece of science the demo implements itself rather than
calling the toolkit for (NovaKit's protein_interactions.analyze_interface needs
Modal engines this demo doesn't have). So it has to be shown to agree with
scripts/interfaces.py rather than merely asserted to.

The check: extract the known-ligand contacts of 6Q4G with this module and compare
against results/pipeline/cdk2-second-target/target/known_ligand_contacts.json,
which the heavy pipeline computed with the same 4.5 A heavy-atom definition.

Run:  ./env-kit/bin/python -m demo.selftest
Needs network (fetches 6Q4G from RCSB), no credentials, no credits.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

from demo import contacts

REPO = Path(__file__).resolve().parent.parent
TRUTH = REPO / "results/pipeline/cdk2-second-target/target/known_ligand_contacts.json"
CACHE = REPO / "results" / "demo" / "_cache"


def fetch_6q4g() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / "6Q4G.pdb"
    if not path.exists():
        url = "https://files.rcsb.org/download/6Q4G.pdb"
        with urllib.request.urlopen(url, timeout=120) as handle:
            path.write_bytes(handle.read())
    return path


def fetch_1ppb() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / "1PPB.pdb"
    if not path.exists():
        url = "https://files.rcsb.org/download/1PPB.pdb"
        with urllib.request.urlopen(url, timeout=120) as handle:
            path.write_bytes(handle.read())
    return path


def insertion_code_checks() -> list[tuple[str, bool, str]]:
    """Thrombin: residues must survive chymotrypsin insertion codes.

    1PPB chain H has 259 residues under 231 distinct resseq values. Keying on
    resseq alone drops 28 of them, including the entire 60-loop that lines the
    active site -- silently, and with plausible-looking output.
    """
    pdb = fetch_1ppb()
    numbers, _ = contacts.chain_sequence(pdb, "H")
    loop = [n for n in numbers if str(n).startswith("60")]
    found = contacts.contacts_to_ligand(pdb, chain="H")["residues"]

    # PPACK is a D-Phe-Pro-Arg chloromethylketone bound in the active site; its
    # contacts must include the catalytic Ser195, the S1 aspartate Asp189, and
    # the 60-loop cage. If insertion codes were dropped, 60A/60D vanish.
    return [
        ("thrombin residues", len(numbers) == 259, f"{len(numbers)} (expect 259)"),
        ("60-loop intact", len(loop) == 10, f"{loop}"),
        ("catalytic Ser195", 195 in found, str(195 in found)),
        ("S1 Asp189", 189 in found, str(189 in found)),
        ("60-loop contacts", "60A" in found and "60D" in found, f'60A={"60A" in found} 60D={"60D" in found}'),
        ("residue ordering", contacts.residue_order(60) < contacts.residue_order("60A")
         < contacts.residue_order("60B") < contacts.residue_order(61), "60 < 60A < 60B < 61"),
    ]


def main() -> int:
    truth = json.loads(TRUTH.read_text())
    expected = truth["known_ligand_residue_ids"]

    got = contacts.contacts_to_ligand(fetch_6q4g(), chain=truth["chain_id"])
    residues = got["residues"]

    checks = [
        ("ligand identified", got["ligand_names"] == [truth["ligand_comp_id"]],
         f'{got["ligand_names"]} vs {truth["ligand_comp_id"]}'),
        ("heavy atom count", got["n_ligand_atoms"] == truth["n_ligand_heavy_atoms"],
         f'{got["n_ligand_atoms"]} vs {truth["n_ligand_heavy_atoms"]}'),
        ("contact residues", residues == expected,
         f"jaccard={contacts.jaccard(residues, expected)}"),
    ]

    # consensus() and coverage() on a known input
    signature = contacts.consensus([[1, 2, 3], [2, 3, 4], [3, 4, 5]], core_threshold=0.6)
    checks.append(("consensus core", signature["core"] == [2, 3, 4], str(signature["core"])))
    checks.append(
        ("coverage", contacts.coverage([2, 3], [2, 3, 4]) == 2 / 3,
         str(contacts.coverage([2, 3], [2, 3, 4]))),
    )

    checks.extend(insertion_code_checks())

    failed = 0
    for name, ok, detail in checks:
        print(f'  {"PASS" if ok else "FAIL"}  {name:20s} {detail}')
        failed += not ok

    print(f'\n{len(checks) - failed}/{len(checks)} checks passed')
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
