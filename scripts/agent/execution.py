"""Local Python, R, and Bash execution for the agent harness.

Adapted from Biomni's Apache-2.0-licensed execution helpers at commit 400c1f3,
by way of the NovaKit adaptation in jaysahni/cheminformatics-kit. Changes here:
no viewer bridge, and the persistent namespace is seeded with this repo's
``scripts/`` directory on ``sys.path`` plus its tool registry.

Warning: this executes model-authored code with the invoking user's
permissions. It is not a sandbox.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import threading
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

from agent.tools import REPO_ROOT, SCRIPTS_DIR, inspect_tool, list_tools


class LocalExecutor:
    """Execute agent-authored code with a persistent Python namespace.

    Python cannot safely terminate a running thread. A timed-out Python block is
    therefore reported as timed out but its daemon thread may continue until the
    process exits. R and Bash subprocesses are terminated by ``subprocess``.
    """

    def __init__(self, path: str | Path | None = None) -> None:
        self.path = Path(path).resolve() if path is not None else REPO_ROOT
        self.shared_directory = self.path / ".agent_tmp"
        self.shared_directory.mkdir(parents=True, exist_ok=True)
        if str(SCRIPTS_DIR) not in sys.path:
            sys.path.insert(0, str(SCRIPTS_DIR))
        self.namespace: dict[str, Any] = {
            "__name__": "__autorepurpose_agent__",
            "working_directory": self.path,
            "shared_directory": self.shared_directory,
            "repo_root": REPO_ROOT,
            "scripts_dir": SCRIPTS_DIR,
            "list_tools": list_tools,
            "inspect_tool": inspect_tool,
        }

    def execute(self, code: str, timeout_seconds: int) -> str:
        stripped = code.strip()
        if stripped.startswith("#!R"):
            return self.run_r(stripped.removeprefix("#!R").lstrip(), timeout_seconds)
        if stripped.startswith("#!BASH") or stripped.startswith("#!CLI"):
            marker = "#!BASH" if stripped.startswith("#!BASH") else "#!CLI"
            return self.run_bash(
                stripped.removeprefix(marker).lstrip(), timeout_seconds
            )
        return self.run_python(code, timeout_seconds)

    def run_python(self, code: str, timeout_seconds: int) -> str:
        output = StringIO()
        result: dict[str, str] = {}

        def run() -> None:
            try:
                with redirect_stdout(output):
                    exec(code.strip("`").strip(), self.namespace)  # noqa: S102
                result["value"] = output.getvalue()
            except Exception as exc:  # noqa: BLE001 - surface failures to the agent
                partial = output.getvalue()
                message = f"Error running Python code: {exc!r}"
                result["value"] = f"{partial}\n{message}" if partial else message

        thread = threading.Thread(target=run, daemon=True)
        thread.start()
        thread.join(timeout_seconds)
        if thread.is_alive():
            return f"ERROR: Python execution timed out after {timeout_seconds} seconds."
        return result.get("value", "Error: Python execution completed with no result.")

    def run_bash(self, script: str, timeout_seconds: int) -> str:
        if not script.strip():
            return "Error running Bash code: Empty script"
        self.shared_directory.mkdir(parents=True, exist_ok=True)
        wrapped_script = "\n".join(
            [
                'AGENT_WORKING_DIRECTORY="$PWD"',
                'AGENT_SHARED_DIR="$PWD/.agent_tmp"',
                'mkdir -p "$AGENT_SHARED_DIR"',
                "export AGENT_WORKING_DIRECTORY AGENT_SHARED_DIR",
                'export TMPDIR="$AGENT_SHARED_DIR"',
                script,
            ]
        )
        handle = tempfile.NamedTemporaryFile(
            suffix=".sh",
            prefix="autorepurpose-agent-",
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=self.shared_directory,
            delete=False,
        )
        script_path = Path(handle.name)
        try:
            handle.write(wrapped_script)
            handle.close()
            script_argument = "./" + script_path.relative_to(self.path).as_posix()
            result = subprocess.run(
                ["bash", script_argument],
                cwd=self.path,
                env=os.environ.copy(),
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired:
            return f"ERROR: Bash execution timed out after {timeout_seconds} seconds."
        finally:
            handle.close()
            if script_path.exists():
                script_path.unlink()
        if result.returncode:
            return (
                f"Error running Bash code (exit {result.returncode}):\n{result.stderr}"
            )
        return str(result.stdout)

    def run_r(self, code: str, timeout_seconds: int) -> str:
        handle = tempfile.NamedTemporaryFile(suffix=".R", mode="w", delete=False)
        try:
            handle.write(code)
            handle.close()
            result = subprocess.run(
                ["Rscript", handle.name],
                cwd=self.path,
                env=os.environ.copy(),
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
        except FileNotFoundError:
            return "Error running R code: Rscript is not installed."
        except subprocess.TimeoutExpired:
            return f"ERROR: R execution timed out after {timeout_seconds} seconds."
        finally:
            handle.close()
            if os.path.exists(handle.name):
                os.unlink(handle.name)
        if result.returncode:
            return f"Error running R code (exit {result.returncode}):\n{result.stderr}"
        return str(result.stdout)
