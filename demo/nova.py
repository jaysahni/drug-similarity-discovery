"""Thin wrapper around the NovaKit toolkit and Rowan, for the lightweight demo.

Everything the demo needs from the outside world goes through here: credential
loading, tool invocation, workflow polling, a disk cache so a re-run costs zero
credits, and a hard credit guard.

NovaKit (github.com/jaysahni/cheminformatics-kit) requires Python >=3.12,<3.13,
so this module runs under `env-kit/`, NOT the project's `env/` (3.14). See
docs/06-LIGHTWEIGHT-DEMO.md.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "results" / "demo" / "_cache"

# Credit guard. Every billable Rowan call passes max_credits; the running total
# is checked against this ceiling before each submission. Raise it deliberately.
BUDGET_CREDITS = float(os.environ.get("DEMO_BUDGET_CREDITS", "400"))
_spent = 0.0


def load_env() -> None:
    """Read .env into os.environ. Never touches the OS keychain.

    NovaKit falls back to `keyring` when an env var is absent, which blocks on a
    GUI prompt in a non-interactive shell. Setting the var here avoids that path.
    """
    env_file = REPO / ".env"
    if not env_file.exists():
        return
    for line in env_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


def require_key() -> str:
    load_env()
    key = os.environ.get("ROWAN_API_KEY")
    if not key:
        raise SystemExit(
            "ROWAN_API_KEY is not set. Copy .env.example to .env and fill it in.\n"
            "Rowan stages (site, design, cofold) cannot run without it."
        )
    return key


def rowan_client():
    """Return the `rowan` module, authenticated."""
    import rowan

    rowan.api_key = require_key()
    return rowan


# --------------------------------------------------------------------------
# cache
# --------------------------------------------------------------------------


def cached(name: str):
    """Load a cached JSON payload, or None."""
    path = CACHE / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    return None


def cache(name: str, payload: Any) -> Any:
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / f"{name}.json").write_text(json.dumps(payload, indent=1, default=str))
    return payload


# --------------------------------------------------------------------------
# NovaKit tools
# --------------------------------------------------------------------------


def tool(name: str, request: dict[str, Any]) -> Any:
    """Run a local NovaKit tool, returning its ToolResult.

    Preflight first, so a missing dependency is reported as such rather than as
    a stack trace from three layers down.
    """
    import novakit

    pre = novakit.preflight(name, request)
    if not pre.ok:
        raise RuntimeError(
            f"{name} preflight failed: {'; '.join(pre.messages)}"
            + (f" missing deps: {pre.missing_dependencies}" if pre.missing_dependencies else "")
            + (f" missing creds: {pre.missing_credentials}" if pre.missing_credentials else "")
        )
    return novakit.run(name, request)


def tool_data(name: str, request: dict[str, Any]) -> dict[str, Any]:
    return tool(name, request).data


# --------------------------------------------------------------------------
# Rowan workflows
# --------------------------------------------------------------------------


def _charge(credits: float) -> None:
    global _spent
    _spent += credits


def spent() -> float:
    return _spent


def guard(estimate: float) -> None:
    if _spent + estimate > BUDGET_CREDITS:
        raise SystemExit(
            f"credit guard: this call (~{estimate:.0f}) would take the run to "
            f"{_spent + estimate:.0f}, over DEMO_BUDGET_CREDITS={BUDGET_CREDITS:.0f}. "
            "Raise the ceiling deliberately or cut the workload."
        )


def submit(submitter: str, *, max_credits: int, **kwargs) -> str:
    """Submit a Rowan workflow, returning its uuid. Enforces the credit guard."""
    guard(max_credits)
    rowan = rowan_client()
    workflow = getattr(rowan, submitter)(max_credits=max_credits, **kwargs)
    return str(workflow.uuid)


def wait(uuid: str, *, label: str = "", poll: int = 20, timeout: int = 3600) -> dict:
    """Poll a Rowan workflow to completion and return its dumped record.

    Rowan status codes: 2 == completed, 3 == failed, 4 == stopped.
    """
    rowan = rowan_client()
    started = time.time()
    while True:
        workflow = rowan.retrieve_workflow(uuid)
        status = workflow.status
        if status in (2, 3, 4):
            record = workflow.model_dump()
            _charge(float(record.get("credits_charged") or 0.0))
            if status != 2:
                raise RuntimeError(
                    f"Rowan workflow {uuid} ({label}) ended with status {status}. "
                    f"Log tail: {str(record.get('logfile'))[-500:]}"
                )
            return record
        if time.time() - started > timeout:
            raise TimeoutError(f"{label or uuid} still running after {timeout}s")
        elapsed = int(time.time() - started)
        print(f"    [{elapsed:5d}s] {label or uuid}: status={status}", flush=True)
        time.sleep(poll)


def run_workflow(
    cache_key: str, submitter: str, *, label: str, max_credits: int, **kwargs
) -> dict:
    """Submit-and-wait, memoised on disk by cache_key.

    The cache stores the workflow uuid, so a re-run re-reads the finished Rowan
    workflow instead of paying for it twice.
    """
    hit = cached(cache_key)
    if hit and hit.get("uuid"):
        print(f"  cached {label}: {hit['uuid']}", flush=True)
        return wait(hit["uuid"], label=label)
    uuid = submit(submitter, max_credits=max_credits, **kwargs)
    cache(cache_key, {"uuid": uuid, "label": label, "submitter": submitter})
    print(f"  submitted {label}: {uuid}", flush=True)
    return wait(uuid, label=label)
