"""Optional A1-style agent harness for the drug-repurposing pipeline.

Install the extra dependencies first:  pip install -r requirements-agent.txt

``agent.tools`` has no optional dependencies and can be imported on its own to
inspect the stage registry.
"""

from __future__ import annotations

from typing import Any

from agent.tools import Backend, MaintenanceStatus, ToolSpec, inspect_tool, list_tools

__all__ = [
    "Backend",
    "MaintenanceStatus",
    "Orchestrator",
    "ToolSpec",
    "inspect_tool",
    "list_tools",
]


def __getattr__(name: str) -> Any:
    """Import the harness lazily so the registry works without langgraph."""
    if name == "Orchestrator":
        from agent.orchestrator import Orchestrator

        return Orchestrator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
