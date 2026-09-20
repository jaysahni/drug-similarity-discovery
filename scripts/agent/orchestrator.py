"""A minimal A1-style planning and code-execution harness for this repo.

Adapted from ``biomni/agent/a1.py`` at commit
400c1f366b96a35ca253e13c9b06c5076af41d65 (Apache-2.0), by way of the NovaKit
adaptation in jaysahni/cheminformatics-kit. The loop is ported faithfully:
a LangGraph ``generate`` -> ``execute`` cycle over the ``<execute>`` /
``<solution>`` protocol, with LLM stage selection. What is repurposed is the
tool surface and the operating rules, which are this repo's pipeline stages and
this repo's evidence standard rather than NovaKit's registry.

Warning: this harness deliberately executes model-authored Python, R, and Bash
with the invoking user's permissions. Do not point it at an untrusted model.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from pathlib import Path
from typing import Any, Literal, TypedDict, cast

from agent.config import AgentConfig
from agent.errors import MissingDependencyError, ProviderProtocolError
from agent.execution import LocalExecutor
from agent.llm import SourceType, get_llm
from agent.protocol import validated_response_block
from agent.retriever import ToolRetriever
from agent.tools import REPO_ROOT, ToolSpec, list_tools

try:
    from langchain_core.messages import (
        AIMessage,
        BaseMessage,
        HumanMessage,
        SystemMessage,
    )
    from langchain_core.runnables import RunnableConfig
    from langgraph.graph import END, START, StateGraph
except ImportError as _exc:  # pragma: no cover - depends on the install
    raise MissingDependencyError(
        "The agent harness needs langgraph and langchain-core. "
        "Install them with: pip install -r requirements-agent.txt",
        location="agent.orchestrator",
    ) from _exc


class AgentState(TypedDict):
    """State maintained by the LangGraph execution loop."""

    messages: list[BaseMessage]
    next_step: str | None


class Orchestrator:
    """Plan, execute code, and drive this repo's pipeline from a task string."""

    def __init__(
        self,
        path: str | Path | None = None,
        llm: str | None = None,
        source: SourceType | None = None,
        use_tool_retriever: bool | None = None,
        timeout_seconds: int | None = None,
        base_url: str | None = None,
        api_key: str | None = None,
        max_steps: int | None = None,
    ) -> None:
        """Create an orchestrator for the supplied model and working directory."""
        self.config = AgentConfig()
        if llm is not None:
            self.config.model = llm
        if source is not None:
            self.config.source = source
        if use_tool_retriever is not None:
            self.config.use_tool_retriever = use_tool_retriever
        if timeout_seconds is not None:
            self.config.timeout_seconds = timeout_seconds
        if base_url is not None:
            self.config.base_url = base_url
        if api_key is not None:
            self.config.api_key = api_key
        if path is not None:
            self.config.path = str(path)
        if max_steps is not None:
            self.config.max_steps = max_steps
        self.path = (
            Path(self.config.path).resolve() if self.config.path else REPO_ROOT
        )
        self.timeout_seconds = self.config.timeout_seconds
        self.max_steps = self.config.max_steps
        self.use_tool_retriever = self.config.use_tool_retriever
        self.llm = get_llm(
            self.config.model,
            temperature=self.config.temperature,
            stop_sequences=["</execute>", "</solution>"],
            source=self.config.source,  # type: ignore[arg-type]
            base_url=self.config.base_url,
            api_key=self.config.api_key,
            config=self.config,
        )
        self.executor = LocalExecutor(self.path)
        self.retriever = ToolRetriever()
        self.log: list[str] = []
        self._active_specs: list[ToolSpec] = []
        self.configure()

    @staticmethod
    def _text_content(content: Any) -> str:
        if isinstance(content, str):
            return content
        if isinstance(content, list):
            parts: list[str] = []
            for block in content:
                if isinstance(block, dict) and block.get("type") in {
                    "text",
                    "output_text",
                }:
                    parts.append(str(block.get("text") or block.get("content") or ""))
                elif isinstance(block, str):
                    parts.append(block)
            return "".join(parts)
        return str(content)

    def _tool_prompt(self, specs: list[ToolSpec]) -> str:
        tools = "\n\n".join(
            "\n".join(
                (
                    f"### {spec.name}  (scripts/{spec.script})",
                    spec.purpose,
                    "Backends: "
                    f"{', '.join(item.value for item in spec.execution_backends)}",
                    f"Maintenance: {spec.maintenance_status.value}",
                    f"Arguments: {spec.input_schema}",
                )
            )
            for spec in specs
        )
        return f"""You are a drug-repurposing benchmark agent working inside the
drug-similarity-discovery repository.

Make a numbered checklist plan, update it as steps complete, and ground every
claim in an execution observation. In each response output exactly one of:

1. An <execute>...</execute> block containing Python, R, or Bash code.
   Python is the default. Use #!R or #!BASH as the first line for R or Bash.
   Always save results to a variable and print them.
2. A final <solution>...</solution> block.

The code executes locally with the invoking user's permissions, from the repo
root. For files that must survive across blocks, use ``shared_directory`` in
Python or ``$AGENT_SHARED_DIR`` in Bash.

Stages are the runnable entry points under scripts/. Run one as a subprocess:

    #!BASH
    ./env/bin/python scripts/ask.py --target KDR --top 10 --json

Or inspect the registry from Python, where ``list_tools`` and ``inspect_tool``
are already in the namespace:

    spec = inspect_tool("query.ask")
    print(spec.script, spec.maintenance_status, spec.input_schema)

THE EVIDENCE STANDARD OF THIS REPO IS BINDING ON YOU:

- Never report a number that was not computed in an observation you can see.
  No placeholder metrics, no illustrative values, no "expected" results.
- Keep negative results. A representation that underperforms is a finding and
  stays in the table with its number. Never quietly drop a disappointing arm.
- Every headline claim carries n, and every comparison carries a significance
  test.
- Verify absence programmatically before asserting it (empty column, missing
  file, no overlap). Do not assume.
- If something cannot be evaluated, say so and say why. Do not skip it in
  silence.

VOCABULARY BAN (PROJECT_GOAL.md G8, enforced by scripts/lint_language.py in
CI): no score in this repo is a binding-strength measurement. Never write
"affinity", "Kd", "IC50", "potency", "efficacy", or "binding strength" about a
result. Ranking is interface overlap, not affinity.

COST: stages with the ``rowan`` backend spend real Rowan credits, and
``design.signature`` on ``modal`` spends compute. Do not submit a Rowan job
unless the user explicitly asked to spend credits. To propose one, inspect the
stage and describe the command instead of running it.

Prefer ``query.ask``: it is read-only, offline, and answers from committed
files in about a second. Reach for the expensive stages only when the question
genuinely needs a new computation.

Available stages:
{tools}
"""

    def configure(self, specs: list[ToolSpec] | None = None) -> None:
        """Configure the LangGraph loop for the supplied stage surface."""
        self._active_specs = specs or list_tools()
        self.system_prompt = self._tool_prompt(self._active_specs)

        def generate(state: AgentState) -> AgentState:
            prompt_messages = [
                SystemMessage(content=self.system_prompt),
                *state["messages"],
            ]
            content = ""
            failure_kind = "empty"
            for attempt in range(2):
                messages = prompt_messages
                if attempt:
                    messages = [
                        *prompt_messages,
                        HumanMessage(
                            content=(
                                "Your previous response was empty or malformed. Reply "
                                "with exactly one complete, non-empty "
                                "<execute>...</execute> block or one "
                                "<solution>...</solution> block and no other text."
                            )
                        ),
                    ]
                response = self.llm.invoke(messages)
                candidate = self._text_content(
                    getattr(response, "content", response)
                ).strip()
                validated = validated_response_block(candidate)
                if validated is not None:
                    content = validated
                    break
                failure_kind = "empty" if not candidate else "malformed"
            else:
                raise ProviderProtocolError(
                    "The model provider returned two invalid protocol responses.",
                    location="agent.orchestrator.generate",
                    context={"attempts": 2, "last_response_state": failure_kind},
                )
            next_step: Literal["execute", "end"] = (
                "end" if "<solution>" in content.lower() else "execute"
            )
            return {
                "messages": [*state["messages"], AIMessage(content=content)],
                "next_step": next_step,
            }

        def execute(state: AgentState) -> AgentState:
            content = self._text_content(state["messages"][-1].content)
            match = re.search(
                r"<execute>(.*?)</execute>", content, re.DOTALL | re.IGNORECASE
            )
            code = match.group(1) if match else ""
            output = self.executor.execute(code, self.timeout_seconds)
            if len(output) > 10_000:
                output = output[:10_000] + "\n[output truncated]"
            # Anthropic requires alternating user and assistant messages. An
            # execution result follows the model's assistant response, so it
            # must be represented as the next user turn.
            observation = HumanMessage(content=f"<observation>{output}</observation>")
            return {
                "messages": [*state["messages"], observation],
                "next_step": "generate",
            }

        def route(state: AgentState) -> Literal["execute", "end"]:
            return "execute" if state.get("next_step") == "execute" else "end"

        workflow = StateGraph(AgentState)
        workflow.add_node("generate", generate)
        workflow.add_node("execute", execute)
        workflow.add_edge(START, "generate")
        workflow.add_conditional_edges(
            "generate", route, {"execute": "execute", "end": END}
        )
        workflow.add_edge("execute", "generate")
        self.app = workflow.compile()

    def _select_specs(self, prompt: str) -> list[ToolSpec]:
        specs = list_tools()
        if not self.use_tool_retriever:
            return specs
        return self.retriever.select(prompt, specs, self.llm)

    def _history_messages(
        self, history: list[dict[str, str]] | None
    ) -> list[BaseMessage]:
        messages: list[BaseMessage] = []
        for item in history or []:
            role = item.get("role")
            content = item.get("content")
            if role == "user" and isinstance(content, str):
                messages.append(HumanMessage(content=content))
            elif role == "assistant" and isinstance(content, str):
                messages.append(AIMessage(content=content))
        return messages

    def _run(
        self, prompt: str, history: list[dict[str, str]] | None = None
    ) -> Generator[dict[str, Any], None, None]:
        self.configure(self._select_specs(prompt))
        self.log = []
        inputs: AgentState = {
            "messages": [
                *self._history_messages(history),
                HumanMessage(content=prompt),
            ],
            "next_step": None,
        }
        config: RunnableConfig = {"recursion_limit": self.max_steps * 2 + 2}
        for raw_state in self.app.stream(inputs, stream_mode="values", config=config):
            state = cast(AgentState, raw_state)
            output = self._text_content(state["messages"][-1].content)
            self.log.append(output)
            yield {"output": output}

    def go(
        self, prompt: str, *, history: list[dict[str, str]] | None = None
    ) -> tuple[list[str], str]:
        """Run a task and return the execution log plus final message content."""
        final = ""
        for item in self._run(prompt, history):
            final = item["output"]
        return self.log, final

    def go_stream(
        self, prompt: str, *, history: list[dict[str, str]] | None = None
    ) -> Generator[dict[str, Any], None, None]:
        """Run a task and yield each generated message or execution observation."""
        yield from self._run(prompt, history)
