"""Configuration for the agent harness.

Adapted from Biomni's Apache-2.0-licensed configuration approach at
https://github.com/snap-stanford/Biomni/tree/400c1f366b96a35ca253e13c9b06c5076af41d65,
by way of the NovaKit adaptation in jaysahni/cheminformatics-kit. Changes here:
no data lake, no keychain, and the environment prefix is ``AUTOREPURPOSE_AGENT_``.
Provider credentials stay provider-owned (for example ``ANTHROPIC_API_KEY``) and
are never stored on this object by default.
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class AgentConfig:
    """Settings for one local agent-harness instance."""

    model: str = "claude-sonnet-5"
    temperature: float = 0.7
    timeout_seconds: int = 600
    use_tool_retriever: bool = True
    source: str | None = None
    base_url: str | None = None
    api_key: str | None = None
    path: str | None = None
    max_steps: int = 100

    def __post_init__(self) -> None:
        self.model = os.getenv("AUTOREPURPOSE_AGENT_LLM", self.model)
        self.temperature = float(
            os.getenv("AUTOREPURPOSE_AGENT_TEMPERATURE", str(self.temperature))
        )
        self.timeout_seconds = int(
            os.getenv("AUTOREPURPOSE_AGENT_TIMEOUT_SECONDS", str(self.timeout_seconds))
        )
        self.use_tool_retriever = (
            os.getenv(
                "AUTOREPURPOSE_AGENT_USE_TOOL_RETRIEVER", str(self.use_tool_retriever)
            ).lower()
            == "true"
        )
        self.source = os.getenv("AUTOREPURPOSE_AGENT_SOURCE", self.source or "") or None
        self.base_url = (
            os.getenv("AUTOREPURPOSE_AGENT_BASE_URL", self.base_url or "") or None
        )
        self.api_key = (
            os.getenv("AUTOREPURPOSE_AGENT_API_KEY", self.api_key or "") or None
        )
        self.path = os.getenv("AUTOREPURPOSE_AGENT_PATH", self.path or "") or None
        self.max_steps = int(
            os.getenv("AUTOREPURPOSE_AGENT_MAX_STEPS", str(self.max_steps))
        )
