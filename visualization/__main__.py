import argparse
import functools
import http.server
import json
import os
import subprocess
from pathlib import Path

from .bundle import BundleError, build, load_json, validate


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
    serve.add_argument("--bundle-dir", type=Path, required=True)
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
            data = validate(load_json(args.manifest), args.manifest.resolve().parent)
            print(json.dumps({"valid": True, "run_id": data["run_id"], "origin": data["data_origin"]}))
        elif args.command == "build":
            bundle = build(args.manifest, args.out)
            print(json.dumps({"run_id": bundle["manifest"]["run_id"], "output": str(args.out), "warnings": bundle["warnings"]}))
        elif args.command == "report":
            from .report import render_report
            if args.out.exists():
                raise BundleError("Report output already exists; choose a fresh path")
            content = render_report(load_json(args.bundle), args.renders)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(content)
            print(args.out)
        elif args.command == "serve":
            directory = args.bundle_dir.resolve()
            if not (directory / "bundle.json").is_file() or not (directory / "index.html").is_file():
                raise BundleError("Serve only a built visualization directory, never the repository root")
            handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
            with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler) as server:
                print(f"Visualization: http://127.0.0.1:{args.port}", flush=True)
                server.serve_forever()
        elif args.command == "render-pymol":
            if args.out.exists():
                raise BundleError("Render output already exists; choose a fresh path")
            if not args.executable.is_file():
                raise BundleError("PyMOL executable is unavailable; provide an approved installation")
            env = dict(os.environ, AUTOREPURPOSE_VIZ_BUNDLE=str(args.bundle.resolve()), AUTOREPURPOSE_VIZ_RENDER_OUT=str(args.out.resolve()))
            subprocess.run([str(args.executable.resolve()), "-cq", str(Path(__file__).with_name("pymol_render.py"))], env=env, check=True)
            report = args.out / "render_manifest.json"
            if not report.is_file() or not load_json(report)["scenes"]:
                raise BundleError("PyMOL did not produce a nonempty render manifest")
            from .bundle import safe_asset
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
