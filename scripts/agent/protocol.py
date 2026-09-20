"""The <execute> / <solution> response protocol.

Split out of the harness so it carries no optional dependencies and can be
tested without langgraph installed.

Adapted from ``biomni/agent/a1.py`` at commit 400c1f3 (Apache-2.0).
"""

from __future__ import annotations

import re

RESPONSE_BLOCK = re.compile(
    r"\A\s*<(execute|solution)>(.*?)</\1>\s*\Z", re.DOTALL | re.IGNORECASE
)
OPEN_BLOCK = re.compile(r"\A\s*<(execute|solution)>(.*)\Z", re.DOTALL | re.IGNORECASE)


def validated_response_block(content: str) -> str | None:
    """Return one complete, non-empty protocol block or ``None``.

    Some providers omit the closing tag when ``stop_sequences`` include
    ``</execute>`` / ``</solution>``. Reattach that tag when the response is
    otherwise a single opened block with a non-empty body.
    """
    match = RESPONSE_BLOCK.fullmatch(content)
    if match is None:
        opened = OPEN_BLOCK.fullmatch(content)
        if opened is None:
            return None
        kind = opened.group(1).lower()
        inner = opened.group(2)
        if re.search(rf"</{kind}\s*>", inner, re.IGNORECASE):
            return None
        inner = inner.strip()
        if not inner:
            return None
        return f"<{kind}>{inner}</{kind}>"
    if not match.group(2).strip():
        return None
    kind = match.group(1).lower()
    return f"<{kind}>{match.group(2)}</{kind}>"
