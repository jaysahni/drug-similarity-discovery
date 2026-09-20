"""Named LangChain model-provider factory for the agent harness.

Adapted from ``biomni/llm.py`` at commit 400c1f3 (Apache-2.0), by way of the
NovaKit adaptation in jaysahni/cheminformatics-kit. Changes here: a smaller
provider set, ``AUTOREPURPOSE_AGENT_*`` configuration, and provider
credentials read only from the environment.
"""

from __future__ import annotations

import os
from typing import Any, Literal

from agent.config import AgentConfig
from agent.errors import MissingDependencyError

SourceType = Literal["Anthropic", "OpenAI", "Ollama", "Custom"]
ALLOWED_SOURCES = {"Anthropic", "OpenAI", "Ollama", "Custom"}


def _resolve_source(model: str, source: str | None, base_url: str | None) -> str:
    if source is not None:
        if source not in ALLOWED_SOURCES:
            raise ValueError(
                f"Invalid source: {source}. Valid sources: {sorted(ALLOWED_SOURCES)}"
            )
        return source
    configured = os.getenv("AUTOREPURPOSE_AGENT_SOURCE")
    if configured:
        return _resolve_source(model, configured, base_url)
    if model.startswith("claude-"):
        return "Anthropic"
    if model.startswith(("gpt-", "o1-", "o3-")):
        return "OpenAI"
    if base_url is not None:
        return "Custom"
    return "Ollama"


def get_llm(
    model: str,
    *,
    temperature: float = 0.7,
    stop_sequences: list[str] | None = None,
    source: SourceType | None = None,
    base_url: str | None = None,
    api_key: str | None = None,
    config: AgentConfig | None = None,
) -> Any:
    """Return a LangChain chat model for the requested provider."""
    del config  # accepted for signature parity with the upstream harness
    resolved = _resolve_source(model, source, base_url)
    try:
        if resolved == "Anthropic":
            from langchain_anthropic import ChatAnthropic

            return ChatAnthropic(
                model=model,
                temperature=temperature,
                stop=stop_sequences,
                api_key=api_key or os.getenv("ANTHROPIC_API_KEY"),
                base_url=base_url,
                max_tokens=8192,
            )
        if resolved == "OpenAI":
            from langchain_openai import ChatOpenAI

            return ChatOpenAI(
                model=model,
                temperature=temperature,
                stop=stop_sequences,
                api_key=api_key or os.getenv("OPENAI_API_KEY"),
                base_url=base_url,
            )
        if resolved == "Ollama":
            from langchain_ollama import ChatOllama

            return ChatOllama(model=model, temperature=temperature)
        from langchain_openai import ChatOpenAI

        return ChatOpenAI(
            model=model,
            temperature=temperature,
            stop=stop_sequences,
            api_key=api_key or os.getenv("AUTOREPURPOSE_AGENT_API_KEY") or "unset",
            base_url=base_url,
        )
    except ImportError as exc:  # pragma: no cover - depends on the install
        raise MissingDependencyError(
            "The agent harness needs its optional provider dependencies. "
            "Install them with: pip install -r requirements-agent.txt",
            location="agent.llm.get_llm",
            context={"source": resolved, "model": model},
        ) from exc
