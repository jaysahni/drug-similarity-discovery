"""CLI for the agent harness.

    ./env/bin/python -m agent --list
    ./env/bin/python -m agent "which approved drugs engage the KDR site?"

``--list`` needs no optional dependencies and no API key.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.errors import MissingDependencyError  # noqa: E402
from agent.tools import list_tools  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="agent", description=__doc__)
    parser.add_argument("task", nargs="?", help="natural-language task")
    parser.add_argument(
        "--list", action="store_true", help="list stages and exit (no model call)"
    )
    parser.add_argument("--model", default=None, help="model id")
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument(
        "--no-retriever",
        action="store_true",
        help="show the model every stage instead of pre-selecting",
    )
    args = parser.parse_args(argv)

    if args.list:
        for spec in list_tools():
            backends = ",".join(item.value for item in spec.execution_backends)
            print(f"{spec.name:28} [{backends}] {spec.maintenance_status.value}")
            print(f"  scripts/{spec.script}")
            print(f"  {spec.purpose}")
        return 0

    if not args.task:
        parser.error("a task is required unless --list is given")

    try:
        from agent.orchestrator import Orchestrator
    except MissingDependencyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    orchestrator = Orchestrator(
        llm=args.model,
        max_steps=args.max_steps,
        use_tool_retriever=False if args.no_retriever else None,
    )
    for event in orchestrator.go_stream(args.task):
        print(event["output"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
