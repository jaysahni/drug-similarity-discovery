import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from visualization.bundle import BundleError, atoms_for, build, fit, load_json, safe_asset, validate
from visualization.fixture import make_fixture
from visualization.report import render_report


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.path = make_fixture(self.root / "fixture")
        self.data = load_json(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def valid(self):
        return validate(self.data, self.path.parent)

    def test_fixture_valid(self):
        self.valid()

    def test_reject_fake_coverage(self):
        self.data["candidates"][0]["coverage"]["value"] = 0.99
        with self.assertRaisesRegex(BundleError, "Coverage disagrees"):
            self.valid()

    def test_partition(self):
        self.data["candidates"][0]["missed"] = self.data["candidates"][0]["engaged"]
        with self.assertRaisesRegex(BundleError, "partition"):
            self.valid()

    def test_low_confidence_cannot_rank(self):
        self.data["candidates"][0]["confidence"]["gate"] = "failed"
        with self.assertRaisesRegex(BundleError, "pass confidence"):
            self.valid()

    def test_missing_decoys_not_zero(self):
        self.data["candidates"][0]["decoys"]["percentile"] = 0
        with self.assertRaisesRegex(BundleError, "Unavailable calibration"):
            self.valid()

    def test_decoy_context(self):
        d = self.data["candidates"][0]["decoys"]
        d.update(status="available", percentile=90, n=1, scores=[0.1], context="wrong")
        with self.assertRaisesRegex(BundleError, "context"):
            self.valid()

    def test_no_fabricated_computed_label(self):
        self.data["data_origin"] = "computed"
        with self.assertRaisesRegex(BundleError, "Synthetic"):
            self.valid()

    def test_empty_core_is_unavailable(self):
        for r in self.data["signature"]["residues"]:
            r["core"] = False
        with self.assertRaisesRegex(BundleError, "nonempty core"):
            self.valid()

    def test_unknown_schema(self):
        self.data["schema_version"] = "2.0.0"
        with self.assertRaises(BundleError):
            self.valid()

    def test_suppressed_board(self):
        self.data["signature"]["addressability"]["status"] = "not_addressable"
        with self.assertRaisesRegex(BundleError, "addressability"):
            self.valid()
        for c in self.data["candidates"]:
            if c["modality"] == "small_molecule":
                c["rank"] = None
        self.valid()

    def test_path_traversal(self):
        for path in ("../fixture/toy.pdb", "/etc/passwd", "https://example.invalid/x", "..\\x"):
            with self.subTest(path=path), self.assertRaises(BundleError):
                safe_asset(self.path.parent, path)

    def test_symlink_escape(self):
        (self.root / "secret").write_text("test")
        (self.path.parent / "link").symlink_to(self.root / "secret")
        with self.assertRaises(BundleError):
            safe_asset(self.path.parent, "link")

    def test_checksum(self):
        self.data["structures"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(BundleError, "Checksum"):
            self.valid()

    def test_rigid_alignment(self):
        mobile = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        expected = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        ref = mobile @ expected.T + [4., 5., 6.]
        rotation, translation, rmsd, passed = fit(mobile, ref, 0.01)
        np.testing.assert_allclose(rotation, expected, atol=1e-10)
        np.testing.assert_allclose(mobile @ rotation.T + translation, ref, atol=1e-10)
        self.assertAlmostEqual(np.linalg.det(rotation), 1)
        self.assertTrue(passed)
        self.assertLess(rmsd, 1e-10)

    def test_degenerate_alignment(self):
        with self.assertRaisesRegex(BundleError, "collinear"):
            fit([[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]], 1)

    def test_build_and_no_overwrite(self):
        out = self.root / "out"
        bundle = build(self.path, out)
        self.assertIn("SYNTHETIC FIXTURE", (out / "report.html").read_text())
        self.assertEqual(len(bundle["display"]["toy"]["ligand_atoms"]), 6)
        self.assertIn("toy:A:5", bundle["display"]["toy"]["residue_atoms"])
        self.assertTrue((out / "checksums.json").exists())
        with self.assertRaisesRegex(BundleError, "already exists"):
            build(self.path, out)

    def test_unknown_residue_fails_before_output(self):
        self.data["structures"][0]["residue_map"]["toy:A:5"]["insertion"] = "B"
        self.path.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(BundleError, "not present"):
            build(self.path, self.root / "bad")
        self.assertFalse((self.root / "bad").exists())

    def test_html_escapes_evidence(self):
        self.data["evidence"][0]["title"] = "<script>alert(1)</script>"
        self.path.write_text(json.dumps(self.data))
        bundle = build(self.path, self.root / "escaped")
        report = render_report(bundle)
        self.assertIn("&lt;script&gt;", report)
        self.assertNotIn("<script>", report)

    def test_mmcif_roundtrip(self):
        from Bio.PDB import MMCIFIO, PDBParser
        from visualization.bundle import digest
        structure = PDBParser(QUIET=True).get_structure("toy", str(self.path.parent / "toy.pdb"))
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(str(self.path.parent / "toy.cif"))
        asset = self.data["structures"][0]
        asset.update(path="toy.cif", format="cif", sha256=digest(self.path.parent / "toy.cif"))
        self.path.write_text(json.dumps(self.data))
        exported = build(self.path, self.root / "cif-export")
        self.assertEqual(len(exported["display"]["toy"]["ligand_atoms"]), 6)
        self.assertEqual(len(exported["display"]["toy"]["residue_atoms"]["toy:A:5"]), 5)

    def test_mapping_must_be_one_to_one(self):
        asset = self.data["structures"][0]
        asset["residue_map"]["duplicate"] = dict(asset["residue_map"]["toy:A:1"])
        with self.assertRaisesRegex(BundleError, "one-to-one"):
            atoms_for(asset, self.path.parent)

    def test_target_ligand_overlap_rejected(self):
        asset = self.data["structures"][0]
        asset["ligand_residues"] = [dict(asset["residue_map"]["toy:A:1"])]
        with self.assertRaisesRegex(BundleError, "overlap"):
            atoms_for(asset, self.path.parent)

    def test_display_alignment_does_not_mutate_source(self):
        mobile = copy.deepcopy(self.data["structures"][0])
        mobile.update(id="mobile", frame_id="mobile-native", alignment={"reference_id": "toy", "atom_pairs": [[1, 1], [6, 6], [11, 11]], "max_rmsd": 0.01})
        self.data["structures"].append(mobile)
        self.path.write_text(json.dumps(self.data))
        before = (self.path.parent / "toy.pdb").read_bytes()
        exported = build(self.path, self.root / "aligned")
        self.assertEqual(exported["display"]["mobile"]["alignment"]["status"], "aligned")
        self.assertEqual(exported["display"]["mobile"]["frame_id"], "toy-frame")
        self.assertEqual(before, (self.path.parent / "toy.pdb").read_bytes())

    def test_ligand_alignment_prohibited(self):
        mobile = copy.deepcopy(self.data["structures"][0])
        mobile.update(id="mobile", alignment={"reference_id": "toy", "atom_pairs": [[61, 61], [62, 62], [63, 63]], "max_rmsd": 1})
        self.data["structures"].append(mobile)
        self.path.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(BundleError, "never ligand"):
            build(self.path, self.root / "bad-alignment")

    def test_duplicate_atom_rejected(self):
        from visualization.bundle import digest
        asset = self.data["structures"][0]
        path = self.path.parent / "toy.pdb"
        text = path.read_text()
        path.write_text(text.splitlines()[0] + "\n" + text)
        asset["sha256"] = digest(path)
        with self.assertRaisesRegex(BundleError, "could not be parsed"):
            atoms_for(asset, self.path.parent)

    def test_insertion_code_is_selected_exactly(self):
        from visualization.bundle import digest
        asset = self.data["structures"][0]
        path = self.path.parent / "toy.pdb"
        lines = path.read_text().splitlines()
        lines[:5] = [line[:26] + "B" + line[27:] for line in lines[:5]]
        path.write_text("\n".join(lines) + "\n")
        asset["sha256"] = digest(path)
        asset["residue_map"]["toy:A:1"]["insertion"] = "B"
        self.path.write_text(json.dumps(self.data))
        exported = build(self.path, self.root / "insertion")
        self.assertEqual(exported["display"]["toy"]["residue_atoms"]["toy:A:1"], [1, 2, 3, 4, 5])

    def test_partial_run_has_no_invented_scores(self):
        self.data["candidates"] = []
        self.data["structures"] = []
        self.data["target"]["structure_id"] = None
        self.data["signature"].update(status="not_evaluated", residues=[])
        self.path.write_text(json.dumps(self.data))
        bundle = build(self.path, self.root / "partial")
        self.assertEqual(bundle["display"], {})
        self.assertIn("not_evaluated", (self.root / "partial" / "report.html").read_text())


if __name__ == "__main__":
    unittest.main()
