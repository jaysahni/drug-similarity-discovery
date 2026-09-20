import argparse
import functools
import http.server
import json
import os
import subprocess
from pathlib import Path
from urllib.parse import unquote

from .errors import BundleError

ROOT = Path(__file__).resolve().parent
WEB = ROOT / "web"


class Handler(http.server.SimpleHTTPRequestHandler):
    """Serves a bundle directory, with the pinned vendor assets mapped in.

    A built export carries its own copy of vendor/. The source tree does not: the browser
    files sit in visualization/web/ while the pinned 3Dmol build is one level up, so
    serving web/ directly would 404 on it and the viewer would report 3D as unavailable.
    """

    protocol_version = "HTTP/1.1"  # keep-alive; the 1.0 default closes every connection

    def __init__(self, *args, vendor_root=None, **kwargs):
        # Set before super().__init__, which serves the request before it returns.
        self.vendor_root = vendor_root
        super().__init__(*args, **kwargs)

    def translate_path(self, path):
        requested = path.split("?", 1)[0].split("#", 1)[0]
        prefix = "/vendor/"
        if self.vendor_root is not None and requested.startswith(prefix):
            candidate = (self.vendor_root / unquote(requested[len(prefix):])).resolve()
            if candidate.is_relative_to(self.vendor_root) and candidate.is_file():
                return str(candidate)
        return super().translate_path(path)


def main():
    parser = argparse.ArgumentParser(prog="python -m visualization")
    commands = parser.add_subparsers(dest="command", required=True)
    fixture = commands.add_parser("fixture")
    fixture.add_argument("--out", type=Path, required=True)
    for name in ("validate", "build"):
        command = commands.add_parser(name)
        command.add_argument("--manifest", type=Path, required=True)
        if name == "build":
            command.add_argument("--out", type=Path, required=True)
    report = commands.add_parser("report")
    report.add_argument("--bundle", type=Path, required=True)
    report.add_argument("--renders", type=Path)
    report.add_argument("--out", type=Path, required=True)
    serve = commands.add_parser("serve")
    serve.add_argument("--bundle-dir", type=Path, default=WEB)
    serve.add_argument("--port", type=int, default=8000)
    render = commands.add_parser("render-pymol")
    render.add_argument("--bundle", type=Path, required=True)
    render.add_argument("--executable", type=Path, required=True)
    render.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "fixture":
            from .fixture import make_fixture
            print(make_fixture(args.out))
        elif args.command == "validate":
            from .bundle import load_json, validate
            data = validate(load_json(args.manifest), args.manifest.resolve().parent)
            print(json.dumps({"valid": True, "run_id": data["run_id"], "origin": data["data_origin"]}))
        elif args.command == "build":
            from .bundle import build
            bundle = build(args.manifest, args.out)
            print(json.dumps({"run_id": bundle["manifest"]["run_id"], "output": str(args.out), "warnings": bundle["warnings"]}))
        elif args.command == "report":
            from .bundle import load_json
            from .report import render_report
            if args.out.exists():
                raise BundleError("Report output already exists; choose a fresh path")
            content = render_report(load_json(args.bundle), args.renders)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(content)
            print(args.out)
        elif args.command == "serve":
            directory = args.bundle_dir.resolve()
            # A built export carries bundle.json; the source web/ directory carries only the recorded worked examples.
            if not (directory / "index.html").is_file() or not ((directory / "bundle.json").is_file() or (directory / "examples.json").is_file()):
                raise BundleError("Serve only a built visualization directory or visualization/web, never the repository root")
            vendor = ROOT / "vendor"
            fallback = None if (directory / "vendor").is_dir() else (vendor if vendor.is_dir() else None)
            handler = functools.partial(Handler, directory=str(directory), vendor_root=fallback)
            with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler) as server:
                saved = "saved analysis and " if (directory / "bundle.json").is_file() else ""
                print(f"Visualization: http://127.0.0.1:{args.port}", flush=True)
                print(f"Serving the {saved}worked examples from {directory}", flush=True)
                server.serve_forever()
        elif args.command == "render-pymol":
            if args.out.exists():
                raise BundleError("Render output already exists; choose a fresh path")
            if not args.executable.is_file():
                raise BundleError("PyMOL executable is unavailable; provide an approved installation")
            from .bundle import load_json, safe_asset
            env = dict(os.environ, AUTOREPURPOSE_VIZ_BUNDLE=str(args.bundle.resolve()), AUTOREPURPOSE_VIZ_RENDER_OUT=str(args.out.resolve()))
            subprocess.run([str(args.executable.resolve()), "-cq", str(Path(__file__).with_name("pymol_render.py"))], env=env, check=True)
            report = args.out / "render_manifest.json"
            if not report.is_file() or not load_json(report)["scenes"]:
                raise BundleError("PyMOL did not produce a nonempty render manifest")
            for scene in load_json(report)["scenes"]:
                if not safe_asset(args.out, scene["file"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"):
                    raise BundleError("PyMOL failed to produce a valid PNG")
            if safe_asset(args.out, "presentation.pse").stat().st_size == 0:
                raise BundleError("PyMOL session is empty")
            print(report)
    except (BundleError, OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"Visualization error: {exc}\n")


if __name__ == "__main__":
    main()
