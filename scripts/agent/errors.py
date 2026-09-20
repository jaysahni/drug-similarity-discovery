"""Errors raised by the agent harness."""

from __future__ import annotations

from typing import Any


class AgentError(Exception):
    """Base class for agent-harness failures."""

    def __init__(
        self,
        message: str,
        *,
        location: str | None = None,
        context: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.location = location
        self.context = context or {}

    def __str__(self) -> str:
        base = super().__str__()
        if self.location:
            base = f"{base} [{self.location}]"
        if self.context:
            detail = ", ".join(f"{k}={v!r}" for k, v in sorted(self.context.items()))
            base = f"{base} ({detail})"
        return base


class ProviderProtocolError(AgentError):
    """The model provider did not honour the agent's response protocol."""


class MissingDependencyError(AgentError):
    """An optional dependency for the agent harness is not installed."""
