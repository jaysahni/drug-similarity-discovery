"""LLM-based stage selection for the agent harness.

Adapted from Biomni's Apache-2.0-licensed prompt-based tool retriever at commit
400c1f3, by way of the NovaKit adaptation in jaysahni/cheminformatics-kit.
Only registered ``ToolSpec`` records can be selected.
"""

from __future__ import annotations

import re
from typing import Any

from agent.tools import ToolSpec


class ToolRetriever:
    """Ask an LLM to select relevant pipeline stages by canonical name."""

    def select(self, query: str, specs: list[ToolSpec], llm: Any) -> list[ToolSpec]:
        """Return the stages the model selected, or all of them on any doubt."""
        catalog = "\n".join(
            f"- {spec.name}: {spec.purpose}; tags={','.join(spec.tags)}; "
            f"backends={','.join(item.value for item in spec.execution_backends)}; "
            f"status={spec.maintenance_status.value}"
            for spec in specs
        )
        response = llm.invoke(
            "Select the pipeline stages relevant to this task. Return exactly one "
            "line: TOOLS: [canonical.stage.name, ...]. Do not select unavailable "
            "stages unless the task explicitly asks about them.\n\n"
            f"TASK: {query}\n\nAVAILABLE STAGES:\n{catalog}"
        )
        content = getattr(response, "content", response)
        if isinstance(content, list):
            content = "\n".join(
                str(block.get("text", "")) if isinstance(block, dict) else str(block)
                for block in content
            )
        match = re.search(
            r"TOOLS:\s*\[(.*?)\]", str(content), re.IGNORECASE | re.DOTALL
        )
        if not match:
            return specs
        names = {
            item.strip().strip("'\"")
            for item in match.group(1).split(",")
            if item.strip()
        }
        selected = [spec for spec in specs if spec.name in names]
        return selected or specs
