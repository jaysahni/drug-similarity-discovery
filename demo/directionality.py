"""Would inhibiting this target help, or harm?

The target-selection stage gates on "do approved drugs exist for this protein",
which is necessary and not sufficient. Running the pipeline on MSH2 showed why.

MSH2 is rank 1 for colorectal cancer by Open Targets association score, and it
is a **tumour suppressor**: losing it is what causes Lynch syndrome. Association
scores measure how strongly a gene is linked to a disease; they do not
distinguish a gene that drives the disease from a gene whose *loss* causes it.
So the pipeline happily selected a target where the therapeutic goal is to
restore function, and then searched for molecules to plug its active site --
which, if they worked, would make the disease worse.

Small-molecule repurposing almost always means inhibition. This module asks
whether inhibition is the right direction before anything downstream spends
credits on it.

Signals, strongest first:

  approved inhibitors exist   The strongest possible evidence that inhibiting a
                              target is both achievable and tolerated: regulators
                              have already approved doing it in humans.
  UniProt keyword             "Tumor suppressor" (KW-0043) / "Proto-oncogene"
                              (KW-0656). Curated, free, and directly on point.

The verdict is advisory and is reported, never used to silently drop a target --
a tumour suppressor can be a legitimate target for ACTIVATION or for a degrader,
and this module only speaks to inhibition.

Run:  ./env-kit/bin/python -m demo.directionality
"""

from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "results" / "demo" / "_cache" / "uniprot_keywords"

SUPPRESSOR = "Tumor suppressor"
ONCOGENE = {"Proto-oncogene", "Oncogene"}


def keywords(accession: str) -> list[str]:
    """UniProt keywords for an accession, cached on disk."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{accession}.json"
    if path.exists():
        return json.loads(path.read_text())
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json?fields=keyword"
    with urllib.request.urlopen(url, timeout=60) as handle:
        payload = json.load(handle)
    kws = [k["name"] for k in payload.get("keywords", [])]
    path.write_text(json.dumps(kws))
    return kws


def approved_inhibitor_count(symbol: str) -> int:
    """Approved drugs with this target as their mechanism.

    Used as evidence that inhibition is achievable and tolerated. Note it cannot
    distinguish an inhibitor from an agonist -- the corpus carries no action
    type -- so it is treated as weaker evidence than it looks, and a target with
    zero of them is 'no evidence', not 'harmful'.
    """
    path = REPO / "data" / "approved_drugs.csv"
    n = 0
    for row in csv.DictReader(path.open()):
        if row.get("approved") == "1" and symbol in (row.get("moa_targets") or "").split(";"):
            n += 1
    return n


def classify(symbol: str, accession: str) -> dict:
    """inhibit_ok / inhibit_harmful / unknown, with the evidence that decided it."""
    kws = keywords(accession)
    n_moa = approved_inhibitor_count(symbol)
    is_tsg = SUPPRESSOR in kws
    is_onc = bool(ONCOGENE & set(kws))

    if n_moa > 0 and not is_tsg:
        verdict, why = "inhibit_ok", f"{n_moa} approved drug(s) already act on {symbol}"
    elif is_tsg:
        verdict, why = (
            "inhibit_harmful",
            f"UniProt flags {symbol} as a tumour suppressor; loss of function is what "
            "causes disease, so inhibiting it is the wrong direction",
        )
    elif is_onc:
        verdict, why = "inhibit_ok", f"UniProt flags {symbol} as a proto-oncogene"
    else:
        verdict, why = (
            "unknown",
            f"no approved drug acts on {symbol} and UniProt carries no "
            "tumour-suppressor or proto-oncogene keyword",
        )

    # A target can be BOTH: approved drugs exist and it is a suppressor. Say so
    # rather than letting the first rule win silently.
    conflict = is_tsg and n_moa > 0
    return {
        "symbol": symbol,
        "uniprot": accession,
        "verdict": verdict,
        "why": why,
        "n_approved_moa_drugs": n_moa,
        "tumour_suppressor": is_tsg,
        "proto_oncogene": is_onc,
        "conflicting_evidence": conflict,
        "caveat": "speaks only to INHIBITION; activation or targeted degradation "
                  "of a tumour suppressor may still be valid",
    }


def main() -> None:
    import sys

    out_dir = REPO / "results" / "demo" / "autoresearch"
    files = sorted(out_dir.glob("*.json")) if out_dir.exists() else []
    if len(sys.argv) > 1:
        files = [Path(a) for a in sys.argv[1:]]
    if not files:
        raise SystemExit("no autoresearch output found")

    for path in files:
        d = json.loads(path.read_text())
        targets = d.get("targets") or []
        print(f"\n=== {path.stem} ===")
        print(f"{'rank':>4s} {'symbol':10s} {'score':>6s} {'moa':>4s} {'verdict':17s} why")
        revised = []
        for t in targets[:15]:
            sym, acc = t["approved_symbol"], t.get("uniprot")
            if not acc:
                continue
            c = classify(sym, acc)
            revised.append(dict(t_rank=t["rank"], **c))
            print(f"{t['rank']:4d} {sym:10s} {t['association_score']:6.3f} "
                  f"{t.get('n_approved_binders_moa', 0):4d} {c['verdict']:17s} {c['why'][:56]}")
        blocked = [r for r in revised if r["verdict"] == "inhibit_harmful"]
        print(f"\n  would now be blocked as wrong-direction: "
              f"{[r['symbol'] for r in blocked] or 'none'}")
        dest = out_dir / f"{path.stem}_directionality.json"
        dest.write_text(json.dumps({"source": path.name, "targets": revised}, indent=1))
        print(f"  -> {dest.relative_to(REPO)}")


if __name__ == "__main__":
    main()
