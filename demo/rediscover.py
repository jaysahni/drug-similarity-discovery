"""Can the pipeline recover a documented repurposing case with the answer hidden?

This is task T10 in docs/02-PIPELINE-PLAN.md, listed as open and never built.

The logic: a pipeline that cannot rediscover a KNOWN repurposing case will not
find an unknown one, and any "novel" hit it produces is noise. So before
trusting any ranking, hide a case whose answer is documented and see whether it
comes back.

  target   KDR / VEGFR2 (P35968), a soluble kinase domain with good structures
  control  12 approved drugs whose MECHANISM is KDR (sorafenib, sunitinib,
           axitinib, lenvatinib, ...) -- the method must rank these highly or
           nothing else it says is believable
  hidden   mebendazole and niclosamide: approved anthelmintics, both annotated
           to KDR, neither with KDR as its mechanism. Mebendazole-for-cancer is
           a real, documented repurposing story. These are the "approved for
           something else" hits the pipeline is supposed to surface.

KDR was chosen because thrombin structurally could not produce this result: all
four F2-mechanism drugs are anticoagulants, so there is no pool of drugs
approved for something else to find. KDR has 12 mechanism drugs AND 19
off-mechanism binders.

**The ranking never sees an annotation.** It is ECFP4 Tanimoto to a query drug's
structure. Annotations are used only to score the ranking afterwards, so the
test is blind by construction rather than by promise.

Run:  ./env-kit/bin/python -m demo.rediscover
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from demo import match_direct as md  # noqa: E402

OUT = REPO / "results" / "demo" / "rediscover"
SYMBOL = "KDR"
UNIPROT = "P35968"
# Documented repurposing candidates: approved for parasitic infection, annotated
# to KDR, KDR is not their mechanism. Named here so the choice is auditable and
# fixed before the ranking runs -- not picked from the results.
HIDDEN = ["mebendazole", "niclosamide"]


def load_corpus() -> list[dict]:
    """Approved small molecules, base corpus plus the extended additions."""
    rows = []
    for path in (REPO / "data" / "approved_drugs.csv",
                 REPO / "data" / "approved_drugs_extended.csv"):
        if not path.exists():
            continue
        for r in csv.DictReader(path.open()):
            if r.get("approved") != "1" or not r.get("smiles"):
                continue
            rows.append({
                "name": r["name"],
                "smiles": r["smiles"],
                "targets": r.get("targets") or "",
                "moa_targets": r.get("moa_targets") or "",
                "source": path.name,
            })
    seen, out = set(), []
    for r in rows:                       # first file wins on a name collision
        if r["name"].lower() in seen:
            continue
        seen.add(r["name"].lower())
        out.append(r)
    return out


def main() -> None:
    import representations

    corpus = load_corpus()
    names = [r["name"].lower() for r in corpus]
    moa = [i for i, r in enumerate(corpus) if SYMBOL in r["moa_targets"].split(";")]
    ann = [i for i, r in enumerate(corpus) if SYMBOL in r["targets"].split(";")]
    hidden = [names.index(h) for h in HIDDEN if h in names]

    print(f"corpus {len(corpus)} approved drugs")
    print(f"  {SYMBOL} mechanism drugs (positive control): {len(moa)}")
    print(f"  {SYMBOL} annotated (any):                    {len(ann)}")
    print(f"  hidden repurposing cases:                    {[corpus[i]['name'] for i in hidden]}")
    if len(hidden) != len(HIDDEN):
        raise SystemExit(f"not all hidden cases are in the corpus: {HIDDEN}")

    fp = representations.REPRESENTATIONS["morgan"]([r["smiles"] for r in corpus])

    # Two query modes. Single-query is the honest baseline; group-query uses the
    # mechanism drugs EXCLUDING the hidden cases, which is what a real
    # repurposing search would have available.
    queries = {}
    axitinib = names.index("axitinib")
    queries["single (axitinib)"] = md.tanimoto(fp[axitinib], fp)

    group = [i for i in moa if i not in hidden]
    block = fp[group]
    sims = np.vstack([md.tanimoto(block[k], fp) for k in range(len(group))])
    for k, i in enumerate(group):
        sims[k, i] = -np.inf                       # never score a query by itself
    queries[f"group (max over {len(group)} mechanism drugs)"] = sims.max(axis=0)

    # The null: nearest-neighbour similarity for every corpus drug that is NOT
    # annotated to this target. Tells us what rank a random approved drug gets.
    non_target = [i for i in range(len(corpus)) if i not in set(ann)]

    report = {}
    for label, score in queries.items():
        s = np.array(score, dtype=float)
        for i in group if "group" in label else [axitinib]:
            s[i] = -np.inf                          # query itself is not a hit
        order = list(np.argsort(-s))
        rank = {i: order.index(i) + 1 for i in range(len(corpus))}
        n = len(corpus)

        null = np.array([s[i] for i in non_target if np.isfinite(s[i])])
        entry = {
            "n_corpus": n,
            "control_ranks": sorted(rank[i] for i in moa if i not in (group if "group" in label else [axitinib])),
            "hidden": {
                corpus[i]["name"]: {
                    "rank": rank[i],
                    "of": n,
                    "tanimoto": round(float(s[i]), 4),
                    "percentile_vs_null": round(md.percentile_of(float(s[i]), null), 4),
                    "top_percent": round(100.0 * rank[i] / n, 2),
                } for i in hidden
            },
            "top20": [
                {"name": corpus[i]["name"], "tanimoto": round(float(s[i]), 4),
                 "is_moa": i in moa, "is_annotated": i in ann}
                for i in order[:20]
            ],
        }
        hits = sum(1 for i in order[:50] if i in ann)
        base = len(ann) / n
        entry["enrichment_at_50"] = round((hits / 50) / base, 2) if base else None
        entry["n_annotated_in_top50"] = hits
        report[label] = entry

        print(f"\n--- {label} ---")
        print(f"  control ({SYMBOL} mechanism) ranks: {entry['control_ranks'][:10]}")
        print(f"  enrichment of {SYMBOL}-annotated drugs in top 50: {entry['enrichment_at_50']}x "
              f"({hits}/50 vs base rate {base:.3%})")
        for name, h in entry["hidden"].items():
            print(f"  HIDDEN {name:14s} rank {h['rank']:5d}/{n}  top {h['top_percent']:5.2f}%  "
                  f"T={h['tanimoto']:.3f}  pctile {h['percentile_vs_null']:.3f}")

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{SYMBOL}_rediscovery.json"
    path.write_text(json.dumps({
        "task": "T10 rediscovery -- can a documented repurposing case be recovered blind?",
        "target": {"symbol": SYMBOL, "uniprot": UNIPROT},
        "hidden_cases": HIDDEN,
        "hidden_case_rationale": "approved anthelmintics annotated to KDR whose mechanism is not KDR; "
                                 "mebendazole-for-cancer is a documented repurposing story",
        "representation": "ECFP4-2048 (scripts/representations.py)",
        "blindness": "the ranking is structural similarity only and never reads an annotation; "
                     "annotations are used solely to score the ranking afterwards",
        "results": report,
        "caveats": [
            "Rediscovery is not discovery. Recovering a known case is the minimum "
            "bar for trusting a novel hit, not evidence of one.",
            "KDR annotations come from the same curated source used to score, so a "
            "hit means 'the annotation agrees', not 'the drug binds'.",
        ],
    }, indent=1, default=str))
    print(f"\n-> {path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
