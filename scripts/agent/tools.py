"""The tool surface the agent is allowed to drive.

This is the repurposed half of the port. Where NovaKit exposed a package-wide
registry of scientific tools, this repo's tools are its pipeline stages: the
runnable entry points under ``scripts/`` that write their outputs to
``results/``.

Purposes are not written here. Each spec's ``purpose`` is read from the target
script's own module docstring at load time (via ``ast``, so nothing is
imported and no script runs as a side effect). A stage whose file is missing is
reported as ``UNAVAILABLE`` rather than being silently dropped, so the agent
sees the gap instead of a shortened menu.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SCRIPTS_DIR.parent


class Backend(Enum):
    """Where a stage does its work."""

    LOCAL = "local"
    ROWAN = "rowan"
    MODAL = "modal"


class MaintenanceStatus(Enum):
    """Whether a stage is expected to work as described."""

    ACTIVE = "active"
    EXPERIMENTAL = "experimental"
    UNAVAILABLE = "unavailable"


@dataclass
class ToolSpec:
    """One runnable pipeline stage."""

    name: str
    script: str
    purpose: str
    tags: list[str] = field(default_factory=list)
    execution_backends: list[Backend] = field(default_factory=lambda: [Backend.LOCAL])
    maintenance_status: MaintenanceStatus = MaintenanceStatus.ACTIVE
    input_schema: dict[str, str] = field(default_factory=dict)

    @property
    def path(self) -> Path:
        return SCRIPTS_DIR / self.script


# (name, script, tags, backends, input_schema). Purpose comes from the docstring.
_STAGES: list[tuple[str, str, list[str], list[Backend], dict[str, str]]] = [
    (
        "query.ask",
        "ask.py",
        ["read-only", "query", "offline"],
        [Backend.LOCAL],
        {
            "target_or_disease": "positional, e.g. 'colorectal cancer' or --target KDR",
            "--top": "int, rows to print",
            "--json": "flag, machine-readable output",
            "--explain": "drug name, per-signature detail",
            "--list": "flag, show prepared targets",
        },
    ),
    (
        "pipeline.run",
        "autorepurpose.py",
        ["pipeline", "end-to-end"],
        [Backend.LOCAL, Backend.ROWAN],
        {"--target": "gene symbol", "--disease": "disease slug"},
    ),
    (
        "repurpose.shortlist",
        "repurpose.py",
        ["pipeline", "library"],
        [Backend.LOCAL],
        {"subcommand": "shortlist", "--slug": "target slug"},
    ),
    (
        "repurpose.submit",
        "repurpose.py",
        ["pipeline", "cofold", "credits"],
        [Backend.ROWAN],
        {"subcommand": "submit", "--slug": "target slug"},
    ),
    (
        "repurpose.collect",
        "repurpose.py",
        ["pipeline", "cofold"],
        [Backend.ROWAN],
        {"subcommand": "collect", "--slug": "target slug"},
    ),
    (
        "repurpose.score",
        "repurpose.py",
        ["pipeline", "ranking"],
        [Backend.LOCAL],
        {"subcommand": "score", "--slug": "target slug"},
    ),
    (
        "target.prep",
        "prep_target.py",
        ["structure", "site"],
        [Backend.LOCAL],
        {"--target": "gene symbol"},
    ),
    (
        "site.p2rank",
        "run_p2rank.py",
        ["structure", "pocket", "baseline"],
        [Backend.LOCAL],
        {"--target": "gene symbol"},
    ),
    (
        "design.signature",
        "boltzgen_signature.py",
        ["design", "signature"],
        [Backend.LOCAL, Backend.MODAL],
        {"--slug": "target slug"},
    ),
    (
        "baseline.chemical",
        "match_candidates.py",
        ["baseline", "cheminformatics"],
        [Backend.LOCAL],
        {"--target": "gene symbol"},
    ),
    (
        "validate.hotspot_recovery",
        "hotspot_recovery.py",
        ["validation", "benchmark"],
        [Backend.LOCAL],
        {},
    ),
    (
        "validate.role_separation",
        "role_separation.py",
        ["validation", "statistics"],
        [Backend.LOCAL],
        {},
    ),
    (
        "validate.metrics",
        "metrics.py",
        ["validation", "statistics"],
        [Backend.LOCAL],
        {},
    ),
    (
        "report.build",
        "report.py",
        ["reporting"],
        [Backend.LOCAL],
        {"--slug": "target slug"},
    ),
    (
        "report.lint_language",
        "lint_language.py",
        ["reporting", "ci", "vocabulary"],
        [Backend.LOCAL],
        {},
    ),
]


def _docstring_of(path: Path) -> str | None:
    """Return a script's module docstring without importing it."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return None
    doc = ast.get_docstring(tree)
    if not doc:
        return None
    return " ".join(doc.strip().split("\n\n")[0].split())


def list_tools() -> list[ToolSpec]:
    """Return every declared stage, each carrying its real docstring purpose."""
    specs: list[ToolSpec] = []
    for name, script, tags, backends, schema in _STAGES:
        path = SCRIPTS_DIR / script
        doc = _docstring_of(path)
        if doc is None:
            status = MaintenanceStatus.UNAVAILABLE
            purpose = (
                f"UNAVAILABLE: scripts/{script} is missing or unparseable in this "
                "checkout; no purpose could be read from it."
            )
        else:
            status = MaintenanceStatus.ACTIVE
            purpose = doc
        specs.append(
            ToolSpec(
                name=name,
                script=script,
                purpose=purpose,
                tags=tags,
                execution_backends=backends,
                maintenance_status=status,
                input_schema=schema,
            )
        )
    return specs


def inspect_tool(name: str) -> ToolSpec:
    """Return one stage by canonical name."""
    for spec in list_tools():
        if spec.name == name:
            return spec
    known = ", ".join(sorted(spec.name for spec in list_tools()))
    raise KeyError(f"Unknown tool {name!r}. Known tools: {known}")
