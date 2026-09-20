"""CI entry point for the PROJECT_GOAL.md G8 language ban.

G8 says a CI lint must keep affinity vocabulary out of the results: nothing in
this repo predicts affinity, potency, Kd or IC50, so naming them next to a score
would describe a measurement that was never made. Until now the only enforcement
was inside scripts/report.py at generation time, which does not run in CI and
does not cover the README - so the "CI lints it" claim was untrue. This makes it
true.

The pattern list and the exemption mechanism are imported from report.py rather
than restated, so there is exactly one definition of what is banned.

Exemptions: an HTML report marks disclaimer regions with <!--LANG-EXEMPT-->.
Markdown files use the same markers in an HTML comment. Those regions are the
one place the repo is allowed to name what it does NOT compute.

Usage:
    ./env/bin/python scripts/lint_language.py                # default targets
    ./env/bin/python scripts/lint_language.py README.md ...  # explicit targets
Exit 1 on any violation, so CI fails.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report import BANNED, language_guard  # noqa: E402  single source of truth

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ["README.md", "results/report_colorectal-cancer.html"]


def main():
    targets = sys.argv[1:] or DEFAULT
    total, checked = 0, 0
    for rel in targets:
        p = ROOT / rel
        if not p.exists():
            print(f"  skip   {rel} (not present)")
            continue
        checked += 1
        hits = language_guard(p.read_text())
        if hits:
            total += len(hits)
            print(f"  FAIL   {rel}: {len(hits)} violation(s)")
            for h in hits:
                print(f"           {h}")
        else:
            print(f"  ok     {rel}")
    print(f"\n{checked} file(s) checked against {len(BANNED)} banned patterns; "
          f"{total} violation(s)")
    if total:
        print("\nPROJECT_GOAL.md G8: no affinity vocabulary may appear outside an "
              "explicit disclaimer block. Reword, or wrap a genuine disclaimer in "
              "<!--LANG-EXEMPT--> ... <!--/LANG-EXEMPT-->.")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
