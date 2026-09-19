import copy
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from visualization.bundle import build
from visualization.fixture import make_fixture
from visualization.pymol_render import render


class RecordingCommands:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def invoke(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            return {"get_position": [0, 0, 0], "get_view": [0] * 18, "get_version": ["TEST DOUBLE — NOT A PYMOL RENDER"]}.get(name)
        return invoke


class PyMOLProtocolTests(unittest.TestCase):
    def test_command_protocol_is_not_render_verification(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = make_fixture(root / "source")
            manifest = json.loads(source.read_text())
            manifest["design_reference"] = {"id": "design-1", "name": "Synthetic design reference", "structure_id": "toy", "protocol": "synthetic", "source": "UI fixture only"}
            source.write_text(json.dumps(manifest))
            bundle = build(source, root / "bundle")
            commands = RecordingCommands()
            with patch.dict(sys.modules, {"pymol": types.SimpleNamespace(cmd=commands)}):
                render(root / "bundle" / "bundle.json", root / "protocol")
            record = json.loads((root / "protocol" / "render_manifest.json").read_text())
            self.assertIn("TEST DOUBLE", record["pymol_version"])
            self.assertEqual(len(record["scenes"]), 12)
            self.assertTrue(any(s["scene"] == "design_design-1" for s in record["scenes"]))
            self.assertTrue(all("SYNTHETIC FIXTURE" in s["caption"] for s in record["scenes"]))
            scene = next(s for s in record["scenes"] if s["scene"] == "candidate_fixture-0")
            self.assertEqual(scene["engaged"], bundle["manifest"]["candidates"][0]["engaged"])
            self.assertEqual(scene["missed"], bundle["manifest"]["candidates"][0]["missed"])
            png_calls = [c for c in commands.calls if c[0] == "png"]
            self.assertEqual(len(png_calls), len(record["scenes"]))
            self.assertTrue(all(c[2]["width"] == 1920 and c[2]["height"] == 1080 for c in png_calls))
            self.assertFalse(list((root / "protocol").glob("*.png")))


if __name__ == "__main__":
    unittest.main()
