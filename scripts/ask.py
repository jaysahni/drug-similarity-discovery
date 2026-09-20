#!/usr/bin/env python3
"""Ask the repo a question, get the ranked drugs back.

    ./env/bin/python scripts/ask.py --list
    ./env/bin/python scripts/ask.py "colorectal cancer"
    ./env/bin/python scripts/ask.py --target KDR --top 10
    ./env/bin/python scripts/ask.py --target KDR --json
    ./env/bin/python scripts/ask.py --target KDR --explain lapatinib
    ./env/bin/python scripts/ask.py --target CDK2

This is the read-only front door. It answers ONLY from results that are already
committed in this repo: no network, no Rowan credits, no GPU, no re-computation.
Everything it prints is read out of a file, and every screenful names the files it
came from so a number can be traced back.

What an answer means (PROJECT_GOAL.md 4.4, G8):
  Each drug was co-folded into the target's binding site and scored on how much
  of that site's core it engages. That is a structural hypothesis about shared
  site engagement. It is not a prediction of clinical benefit, and the rows
  labelled "known binder" are the positive controls of the benchmark - they are
  there to check the ranking, not to be reported as discoveries.

No score anywhere in this repo is a binding-strength measurement; none is read,
predicted or implied (scripts/lint_language.py enforces the vocabulary in CI).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE_DIR = ROOT / "results" / "pipeline"
ROLE_SEP = ROOT / "results" / "role_separation.json"

WIDTH = 82
RULE = "-" * WIDTH

ROLE_LABEL = {
    "known_binder": "known binder (control)",
    "hard_decoy": "HARD DECOY",
    "decoy": "decoy",
}


# --------------------------------------------------------------- small helpers
class AskError(Exception):
    """Something is missing or unreadable. Carries a human fix, never a traceback."""

    def __init__(self, message: str, fix: str = ""):
        super().__init__(message)
        self.message = message
        self.fix = fix


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def load_json(path: Path, what: str, fix: str = "") -> dict:
    """Read a JSON file or raise AskError naming the file and how to make it."""
    p = Path(path)
    if not p.exists():
        raise AskError(f"{what} is missing: {rel(p)}", fix)
    try:
        return json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        raise AskError(f"{what} is not readable JSON: {rel(p)} ({exc})",
                       fix or "the file is corrupt; regenerate it") from None
    except OSError as exc:
        raise AskError(f"{what} could not be read: {rel(p)} ({exc})", fix) from None


def maybe_json(path: Path):
    """Optional file: return its contents, or None. Never raises."""
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def norm(text: str) -> str:
    return " ".join(re.sub(r"[^0-9a-z]+", " ", str(text).lower()).split())


def wrap_list(values, indent: str) -> str:
    """Comma-joined residue list, labelled once and hanging-indented after."""
    body = ", ".join(str(v) for v in values) if values else "(none)"
    return textwrap.fill(body, width=WIDTH, initial_indent=indent,
                         subsequent_indent=" " * len(indent))


def trunc(text: str, n: int) -> str:
    text = str(text)
    return text if len(text) <= n else text[: n - 1] + "…"


# ------------------------------------------------------------------- registry
def build_registry() -> list[dict]:
    """Every prepared target found on disk, with its aliases and its state.

    A target is discovered, not hard-coded: each directory under
    results/pipeline/ that has a target/pocket.json is a target, and its board
    is results/repurpose_<slug>.json if that file exists.
    """
    if not PIPE_DIR.is_dir():
        raise AskError(
            f"no pipeline directory: {rel(PIPE_DIR)}",
            "nothing has been prepared yet - run scripts/autorepurpose.py")

    entries = []
    for d in sorted(p for p in PIPE_DIR.iterdir() if p.is_dir()):
        pocket = maybe_json(d / "target" / "pocket.json")
        structures = maybe_json(d / "target" / "structures.json")
        board_path = ROOT / "results" / f"repurpose_{d.name}.json"
        board_exists = board_path.exists()
        board_head = maybe_json(board_path) if board_exists else None

        uni = ((structures or {}).get("uniprot", {})
               if isinstance(structures, dict) else {})
        gene = uni.get("gene") or (board_head or {}).get("target", {}).get("symbol")
        entry = {
            "slug": d.name,
            "dir": d,
            "gene": gene,
            "protein_name": uni.get("protein_name"),
            "entry_name": uni.get("entry_name"),
            "uniprot": (pocket or {}).get("uniprot_id") or uni.get("accession"),
            "structure": (pocket or {}).get("structure_id")
                         or (structures or {}).get("primary", {}).get("pdb_id"),
            "chain": (pocket or {}).get("chain_id"),
            "pocket": pocket,
            "pocket_path": d / "target" / "pocket.json",
            "board_path": board_path,
            "has_board": bool(board_head and board_head.get("results")),
            "board_exists": board_exists,
            "board_head": board_head,
            "candidates_path": d / "candidates.json",
            "n_designs": (len(list((d / "designs").rglob("design_*")))
                          if (d / "designs").is_dir() else 0),
            "n_cofold_poses": len(list((d / "cofold_poses").iterdir()))
                              if (d / "cofold_poses").is_dir() else 0,
        }
        if entry["has_board"]:
            entry["state"] = "board"
        elif board_exists:
            # the file is there but unreadable or has no results: that is a file
            # problem, and must never be reported as "this was never screened"
            entry["state"] = "board_unusable"
            entry["board_problem"] = ("not readable JSON" if board_head is None
                                      else "readable JSON with no 'results' array")
        elif pocket:
            entry["state"] = "prepared"
        else:
            entry["state"] = "incomplete"

        # aliases: everything a person might reasonably type
        aliases = {entry["slug"], entry["slug"].replace("-", " ")}
        for key in ("gene", "protein_name", "entry_name", "uniprot", "structure"):
            if entry[key]:
                aliases.add(str(entry[key]))
        if entry["entry_name"]:
            aliases.add(str(entry["entry_name"]).split("_")[0])

        targets = maybe_json(d / "targets.json")
        if isinstance(targets, dict):
            meta = targets.get("meta", {})
            for key in ("disease_name", "disease_query"):
                if meta.get(key):
                    aliases.add(meta[key])
                    entry.setdefault("disease", meta[key])
            for row in targets.get("targets", []) or []:
                if row.get("uniprot") and row["uniprot"] == entry["uniprot"]:
                    aliases.update(row.get("symbol_synonyms") or [])
                    if row.get("approved_symbol"):
                        aliases.add(row["approved_symbol"])
                    if row.get("approved_name"):
                        aliases.add(row["approved_name"])
        entry["aliases"] = {norm(a) for a in aliases if norm(a)}
        entries.append(entry)

    if not entries:
        raise AskError(
            f"no prepared targets under {rel(PIPE_DIR)}",
            "run scripts/autorepurpose.py to prepare one")
    return entries


def resolve(query: str, registry: list[dict]) -> tuple[dict, dict]:
    """Map a free-text question or a gene symbol onto one prepared target.

    Deliberately strict: a loose word overlap is NOT a match. "pancreatic cancer"
    shares the token "cancer" with the colorectal run and must still come back as
    "nothing prepared answers that" rather than quietly handing over the KDR
    board. A question matches an alias only if one is a prefix of the other, or
    if every token of one appears in the other. Anything less exact than a
    verbatim alias is reported back to the reader as "read as ...".
    """
    q = norm(query)
    if not q:
        raise AskError("empty question", "try: ask.py --list")

    exact = [e for e in registry if q in e["aliases"]]
    if len(exact) == 1:
        return exact[0], {"exact": True, "matched_alias": q, "query": query}
    if len(exact) > 1:
        raise AskError(
            f'"{query}" matches {len(exact)} targets: '
            + ", ".join(e["gene"] or e["slug"] for e in exact),
            "name one of them with --target")

    qt = set(q.split())
    scored = []
    for e in registry:
        best = (0, None)
        for a in e["aliases"]:
            at = set(a.split())
            if a.startswith(q) or q.startswith(a):
                # credit the overlap, not the alias's length: a one-letter query
                # must stay ambiguous instead of picking the longest name
                score = 3 * min(len(a), len(q))
            elif qt <= at or at <= qt:
                score = 2 * len(" ".join(sorted(qt & at)))
            else:
                continue
            if score > best[0]:
                best = (score, a)
        if best[1]:
            scored.append((best[0], best[1], e))

    if scored:
        scored.sort(key=lambda t: -t[0])
        topscore = scored[0][0]
        top = [t for t in scored if t[0] == topscore]
        if len({t[2]["slug"] for t in top}) == 1:
            _, alias, e = top[0]
            return e, {"exact": False, "matched_alias": alias, "query": query}
        raise AskError(
            f'"{query}" is ambiguous: '
            + ", ".join(sorted({t[2]["gene"] or t[2]["slug"] for t in top})),
            "name one of them with --target")

    known = ", ".join(sorted({e["gene"] or e["slug"] for e in registry}))
    raise AskError(f'nothing prepared here answers to "{query}"',
                   f"prepared targets: {known}. Run ask.py --list to see their "
                   f"state, and what it would take to add one")


# ------------------------------------------------------------------ board maths
def rank_under(board: dict, signature: str) -> dict[str, int]:
    """Rank every scored drug under one signature, by the board's own metric.

    Same rule the board used: sort descending on ranking_metric, stable in file
    order. Verified against the stored `rank` field for the ranked_by signature;
    a disagreement is surfaced rather than hidden.
    """
    metric = board.get("ranking_metric", "precision_in_core")
    scored = [r for r in board.get("results", []) if r.get("status") == "scored"]
    have = [r for r in scored if signature in (r.get("by_signature") or {})]
    have.sort(key=lambda r: -float(r["by_signature"][signature].get(metric) or 0.0))
    return {r["name"]: i + 1 for i, r in enumerate(have)}


def board_rows(board: dict, signature: str) -> tuple[list[dict], list[dict]]:
    """(scored rows in rank order, unscored rows) for one signature."""
    metric = board.get("ranking_metric", "precision_in_core")
    sig_meta = (board.get("signatures") or {}).get(signature, {})
    n_core = sig_meta.get("n_core")
    ranks = rank_under(board, signature)

    scored, unscored = [], []
    for r in board.get("results", []):
        if (r.get("status") != "scored"
                or signature not in (r.get("by_signature") or {})):
            unscored.append({
                "name": r.get("name"),
                "role": r.get("role"),
                "status": r.get("status", "not scored"),
                "struct_id": r.get("struct_id"),
            })
            continue
        s = r["by_signature"][signature]
        engaged_core = s.get("engaged_core") or []
        scored.append({
            "rank": ranks.get(r["name"]),
            "rank_in_file": r.get("rank"),
            "name": r.get("name"),
            "struct_id": r.get("struct_id"),
            "role": r.get("role"),
            "is_known_binder": r.get("role") == "known_binder",
            "is_hard_decoy": r.get("role") == "hard_decoy",
            "score": s.get(metric),
            "score_metric": metric,
            "n_engaged": r.get("n_engaged"),
            "n_core_engaged": len(engaged_core),
            "n_core": n_core,
            "core_coverage": s.get("core_coverage"),
            "engaged_core": engaged_core,
            "missed_core": s.get("missed_core") or [],
            "engaged_residues": r.get("engaged_residues") or [],
        })
    scored.sort(key=lambda r: (r["rank"] is None, r["rank"]))
    return scored, unscored


def limitation(signature: str) -> dict | None:
    """The one measured caveat that matters, read from role_separation.json."""
    data = maybe_json(ROLE_SEP)
    if not isinstance(data, dict):
        return None
    verdicts = data.get("verdict_by_signature") or {}
    if not verdicts:
        return None
    arms = []
    for name, v in verdicts.items():
        auc = v.get("separates_hard_decoy_auc")
        if auc is None:
            continue
        arms.append({
            "signature": name,
            "metric": v.get("metric"),
            "auc": auc,
            "p_two_sided": v.get("separates_hard_decoy_p_two_sided"),
            "separates": bool(v.get("separates_hard_decoy")),
            "n_positive": v.get("separates_hard_decoy_n_positive"),
            "n_negative": v.get("separates_hard_decoy_n_negative"),
        })
    if not arms:
        return None
    arms.sort(key=lambda a: -a["auc"])
    best = arms[0]
    design = next((a for a in arms if a["signature"] == "boltzgen_consensus"), None)
    here = verdicts.get(signature) or {}
    return {
        "source": rel(ROLE_SEP),
        "test": (data.get("test") or {}).get("name", "Mann-Whitney U"),
        "any_signature_separates_hard_decoys": any(a["separates"] for a in arms),
        "best_hard_decoy_arm": best,
        "design_arm": design,
        "this_signature_vs_plain_decoys": {
            "signature": signature,
            "auc": here.get("separates_decoy_auc"),
            "p_holm": here.get("separates_decoy_p_holm"),
            "separates": here.get("separates_decoy"),
            "n_positive": here.get("separates_decoy_n_positive"),
            "n_negative": here.get("separates_decoy_n_negative"),
        } if here else None,
        "arms": arms,
    }


def chem_rank(entry: dict, drug: str):
    """This drug's rank in the chemical-similarity list, if it appears there."""
    data = maybe_json(entry["candidates_path"])
    if not isinstance(data, dict):
        return None
    key = norm(drug)
    for where in ("top", "known_binder_ranks"):
        for row in data.get(where) or []:
            if norm(row.get("name", "")) == key:
                return {
                    "rank": row.get("rank"),
                    "similarity": row.get("similarity"),
                    "representation": data.get("representation"),
                    "n_corpus": data.get("n_corpus"),
                    "reference_ligand": (data.get("ideal_binder")
                                         or {}).get("pref_name"),
                    "listed_in": where,
                    "source": rel(entry["candidates_path"]),
                }
    return {
        "rank": None,
        "not_listed": True,
        "representation": data.get("representation"),
        "n_corpus": data.get("n_corpus"),
        "top_k": data.get("top_k"),
        "reference_ligand": (data.get("ideal_binder") or {}).get("pref_name"),
        "source": rel(entry["candidates_path"]),
    }


# --------------------------------------------------------------------- answers
def answer_board(entry: dict, question: str, top: int, match: dict) -> dict:
    board = load_json(entry["board_path"], "the drug board for this target",
                      f"run: ./env/bin/python scripts/autorepurpose.py run "
                      f"--target {entry['gene'] or entry['slug']} "
                      f"--uniprot {entry['uniprot']}")
    signature = board.get("ranked_by") or "p2rank_geometry"
    sig_meta = (board.get("signatures") or {}).get(signature, {})
    scored, unscored = board_rows(board, signature)
    if not scored:
        raise AskError(
            f"the board {rel(entry['board_path'])} has no scored drugs under "
            f"signature '{signature}'",
            "re-run the score stage of scripts/autorepurpose.py")

    lim = limitation(signature)
    disagree = [r["name"] for r in scored
                if r["rank_in_file"] is not None and r["rank_in_file"] != r["rank"]]
    validation = (board.get("validation_by_signature") or {}).get(signature, {})

    return {
        "ok": True,
        "question": question,
        "answer": "ranked_drugs",
        "resolved_how": match,
        "target": {
            "symbol": entry["gene"], "uniprot": entry["uniprot"],
            "protein_name": entry["protein_name"], "slug": entry["slug"],
            "structure": entry["structure"], "chain": entry["chain"],
            "disease_context": entry.get("disease"),
        },
        "site": {
            "signature": signature,
            "how_defined": sig_meta.get("provenance")
                           or (entry["pocket"] or {}).get("detector"),
            "n_core_residues": sig_meta.get("n_core"),
            "core_residues_auth": sig_meta.get("core_residues_auth"),
            "site_id": (entry["pocket"] or {}).get("site_id"),
        },
        "ranking": {
            "metric": board.get("ranking_metric"),
            "metric_rationale": board.get("ranking_metric_choice"),
            "affinity_used": board.get("affinity_used"),
            "n_rows": len(board.get("results", [])),
            "n_scored": board.get("n_scored", len(scored)),
            "n_known_binders": board.get("n_known_binders"),
            "enrichment_at_top_quartile": validation.get("enrichment_at_top_quartile"),
            "known_binder_ranks": validation.get("known_binder_ranks"),
            "median_known_binder_rank": validation.get("median_rank"),
            "rank_field_disagreements": disagree,
        },
        "matches": scored[:top] if top else scored,
        "n_matches_shown": len(scored[:top] if top else scored),
        "n_matches_total": len(scored),
        "not_scored": unscored,
        "limitation": lim,
        "disclaimer": (
            "Computational hypotheses about shared binding-site engagement, from "
            "co-folded poses. Not a prediction of clinical benefit. Rows labelled "
            "'known binder' are the benchmark's positive controls, not discoveries."
        ),
        "sources": ([rel(entry["board_path"]), rel(entry["pocket_path"])]
                    + ([rel(ROLE_SEP)] if lim else [])),
        "computed_now": False,
    }


def answer_unscreened(entry: dict, question: str, registry: list[dict],
                      match: dict) -> dict:
    """A prepared target with no drug board: say so, say why, give the command."""
    pocket = entry["pocket"] or {}
    ref = next((e for e in registry if e["state"] == "board" and e["board_head"]), None)
    shape = {"n_positives": 10, "n_decoys": 20, "hard_decoys": True,
             "shape_source": "scripts/autorepurpose.py defaults"}
    if ref:
        roles = {}
        for r in ref["board_head"].get("results", []):
            roles[r.get("role")] = roles.get(r.get("role"), 0) + 1
        if roles:
            shape = {
                "n_positives": roles.get("known_binder", 10),
                "n_decoys": roles.get("decoy", 20),
                "hard_decoys": roles.get("hard_decoy", 0) > 0,
                "shape_source": f"same shape as {rel(ref['board_path'])}",
            }
    cmd = (f"./env/bin/python scripts/autorepurpose.py run "
           f"--target {entry['gene'] or entry['slug']} --uniprot {entry['uniprot']} "
           f"--slug {entry['slug']} --n-positives {shape['n_positives']} "
           f"--n-decoys {shape['n_decoys']}"
           + (" --hard-decoys" if shape["hard_decoys"] else ""))

    return {
        "ok": True,
        "question": question,
        "answer": "no_drug_board",
        "resolved_how": match,
        "target": {
            "symbol": entry["gene"], "uniprot": entry["uniprot"],
            "protein_name": entry["protein_name"], "slug": entry["slug"],
            "structure": entry["structure"], "chain": entry["chain"],
        },
        "state": ("prepared, never screened" if not entry["n_cofold_poses"]
                  else "prepared and co-folded, with no scored board on disk"),
        "why": (("the co-folding stage was never run for this target, so no drug "
                 "was ever placed in this pocket and there is no board to rank. "
                 "Nothing is being withheld and nothing is estimated in its place.")
                if not entry["n_cofold_poses"] else
                (f"{entry['n_cofold_poses']} co-folded poses are on disk for this "
                 f"target, but the scored board {rel(entry['board_path'])} is not, "
                 f"so there is nothing to rank. Whether the scoring stage was never "
                 f"run or its output was removed cannot be told from disk, so this "
                 f"script says only what it can see. Nothing is estimated in its "
                 f"place.")),
        "what_exists": {
            "pocket": {
                "path": rel(entry["pocket_path"]),
                "site_id": pocket.get("site_id"),
                "detector": pocket.get("detector"),
                "n_residues": len(pocket.get("residue_ids") or []),
                "residue_ids": pocket.get("residue_ids"),
            },
            "structures": rel(entry["dir"] / "target" / "structures.json"),
            "n_boltzgen_designs": entry["n_designs"],
            "n_cofold_poses": entry["n_cofold_poses"],
        },
        "missing": [rel(entry["board_path"])],
        "to_produce_it": {
            "command": cmd,
            "shape_source": shape["shape_source"],
            "requires": ("network and Rowan credits - the co-folding stage submits "
                         "every drug to Rowan. autorepurpose.py prints its own credit "
                         "estimate (measured 5.4 credits per co-folded drug) and "
                         "refuses to start a stage that exceeds --budget-credits."),
        },
        "sources": [rel(entry["pocket_path"])],
        "computed_now": False,
    }


def answer_explain(entry: dict, drug: str, question: str, match: dict) -> dict:
    board = load_json(entry["board_path"], "the drug board for this target")
    ranked_by = board.get("ranked_by") or "p2rank_geometry"
    metric = board.get("ranking_metric", "precision_in_core")

    key = norm(drug)
    rows = board.get("results", [])
    hit = next((r for r in rows if norm(r.get("name", "")) == key), None)
    if hit is None:
        near = [r["name"] for r in rows if key and key in norm(r.get("name", ""))]
        if len(near) == 1:
            hit = next(r for r in rows if r["name"] == near[0])
        elif near:
            raise AskError(f'"{drug}" matches {len(near)} drugs on this board: '
                           + ", ".join(sorted(near)), "name one exactly")
        else:
            raise AskError(
                f'"{drug}" is not on the {entry["gene"] or entry["slug"]} board',
                f"{len(rows)} drugs were screened; "
                f"run ask.py --target {entry['gene'] or entry['slug']} --top 100 "
                f"to see them all")

    per_sig = {}
    for sig in (board.get("signatures") or {}):
        ranks = rank_under(board, sig)
        s = (hit.get("by_signature") or {}).get(sig)
        if not s:
            per_sig[sig] = {"scored": False,
                            "n_core": (board["signatures"][sig] or {}).get("n_core")}
            continue
        per_sig[sig] = {
            "scored": True,
            "rank": ranks.get(hit["name"]),
            "n_ranked": len(ranks),
            "n_core": (board["signatures"][sig] or {}).get("n_core"),
            "provenance": (board["signatures"][sig] or {}).get("provenance"),
            "score_metric": metric,
            "score": s.get(metric),
            "core_coverage": s.get("core_coverage"),
            "weighted_jaccard": s.get("weighted_jaccard"),
            "f1_engaged_core": s.get("f1_engaged_core"),
            "jaccard_engaged_core": s.get("jaccard_engaged_core"),
            "precision_in_core": s.get("precision_in_core"),
            "engaged_core": s.get("engaged_core") or [],
            "missed_core": s.get("missed_core") or [],
            "is_ranking_signature": sig == ranked_by,
        }

    return {
        "ok": True,
        "question": question,
        "answer": "drug_detail",
        "resolved_how": match,
        "target": {"symbol": entry["gene"], "uniprot": entry["uniprot"],
                   "structure": entry["structure"], "slug": entry["slug"]},
        "drug": {
            "name": hit.get("name"),
            "struct_id": hit.get("struct_id"),
            "role": hit.get("role"),
            "role_label": ROLE_LABEL.get(hit.get("role"), hit.get("role")),
            "status": hit.get("status"),
            "n_engaged": hit.get("n_engaged"),
            "engaged_residues": hit.get("engaged_residues") or [],
            "rank_in_board": hit.get("rank"),
        },
        "ranked_by": ranked_by,
        "by_signature": per_sig,
        "chemical_similarity": chem_rank(entry, hit.get("name", "")),
        "limitation": limitation(ranked_by),
        "limitation_source_unreadable": limitation(ranked_by) is None,
        "disclaimer": (
            "Interface overlap between a co-folded pose and the site core. Not a "
            "prediction of clinical benefit."
        ),
        "sources": [rel(entry["board_path"]), rel(entry["candidates_path"])],
        "computed_now": False,
    }


def answer_list(registry: list[dict]) -> dict:
    out = []
    for e in registry:
        row = {
            "symbol": e["gene"], "uniprot": e["uniprot"], "slug": e["slug"],
            "protein_name": e["protein_name"], "structure": e["structure"],
            "disease_context": e.get("disease"), "state": e["state"],
            "ask": f'ask.py --target {e["gene"] or e["slug"]}',
        }
        if e["state"] == "board" and e["board_head"]:
            b = e["board_head"]
            row.update({
                "board": rel(e["board_path"]),
                "n_drugs": len(b.get("results", [])),
                "n_scored": b.get("n_scored"),
                "ranked_by": b.get("ranked_by"),
                "ranking_metric": b.get("ranking_metric"),
            })
        elif e["state"] == "board_unusable":
            row.update({"board": None,
                        "board_file_present_but_unusable": rel(e["board_path"]),
                        "problem": e.get("board_problem"),
                        "n_boltzgen_designs": e["n_designs"]})
        else:
            row.update({"board": None,
                        "missing": rel(e["board_path"]),
                        "n_boltzgen_designs": e["n_designs"]})
        out.append(row)
    return {"ok": True, "answer": "targets", "n_targets": len(out),
            "targets": out, "computed_now": False,
            "note": ("everything here is read from committed files; no target is "
                     "screened on demand by this script")}


# --------------------------------------------------------------------- printing
def print_resolution(ans: dict) -> None:
    """Say out loud how a question was turned into a target, when it was not exact."""
    m = ans.get("resolved_how") or {}
    name = (ans.get("target") or {}).get("symbol") or ""
    if m.get("query") and not m.get("exact"):
        print(f'read as: "{m["query"]}" -> matched "{m["matched_alias"]}" -> {name}')
    elif m.get("query") and norm(m["query"]) != norm(name):
        print(f'asked as: "{m["query"]}"')


def print_list(ans: dict) -> None:
    print()
    print("WHAT YOU CAN ASK ABOUT")
    print(RULE)
    for t in ans["targets"]:
        head = f"{t['symbol'] or t['slug']}"
        if t["uniprot"]:
            head += f" ({t['uniprot']})"
        if t["protein_name"]:
            head += f" - {trunc(t['protein_name'], 44)}"
        print(head)
        if t["state"] == "board":
            print(f"   READY       {t['n_drugs']} drugs screened, "
                  f"{t['n_scored']} scored, ranked by {t['ranked_by']}")
            ask = f'ask.py "{t["disease_context"]}"' if t.get("disease_context") \
                else t["ask"]
            print(f"   ask         ./env/bin/python scripts/{ask}")
            print(f"   from        {t['board']}")
        elif t["state"] == "prepared":
            print(f"   PREPARED    pocket + structures ready, NOT screened "
                  f"(no co-folding run, so no drug board)")
            print(f"   ask         ./env/bin/python scripts/{t['ask']}   "
                  f"(prints how to produce it)")
            print(f"   missing     {t['missing']}")
        elif t["state"] == "board_unusable":
            print(f"   BROKEN FILE {t['board_file_present_but_unusable']} is "
                  f"{t['problem']}")
            print("   fix         restore it from git; no ranking is shown from a "
                  "file this script could not read")
        else:
            print("   INCOMPLETE  no target/pocket.json; not askable")
        if t.get("structure"):
            print(f"   structure   {t['structure']}")
        print()
    print(RULE)
    print(ans["note"])
    print()


def print_limitation(lim: dict | None) -> None:
    if not lim:
        print("MEASURED LIMITATION - NOT AVAILABLE")
        print(textwrap.fill(
            f"{rel(ROLE_SEP)} could not be read, so the measured caveat that "
            f"belongs here (whether any signature separates known binders from "
            f"hard decoys) is missing. Do not read the ranking above as validated "
            f"until that file is restored and this block prints its numbers.",
            width=WIDTH, initial_indent="  ", subsequent_indent="  "))
        return
    best, design = lim["best_hard_decoy_arm"], lim["design_arm"]
    print("MEASURED LIMITATION - read this before believing a rank")
    if not lim["any_signature_separates_hard_decoys"]:
        body = (f"No signature separates known binders from the HARD DECOYS (drugs "
                f"that engage some other kinase). Best arm: {best['signature']} "
                f"AUC {best['auc']:.3f}, p = {best['p_two_sided']:.3f} "
                f"(n = {best['n_positive']} known binders vs {best['n_negative']} "
                f"hard decoys, two-sided {lim['test']}). The design arm "
                + (f"({design['signature']}) is AUC {design['auc']:.3f}. "
                   if design else "")
                + "So a high rank is a lead to check, not a hit.")
    else:
        body = (f"Best hard-decoy arm: {best['signature']} AUC {best['auc']:.3f}, "
                f"p = {best['p_two_sided']:.3f} (n = {best['n_positive']} vs "
                f"{best['n_negative']}).")
    print(textwrap.fill(body, width=WIDTH, initial_indent="  ",
                        subsequent_indent="  "))
    plain = lim.get("this_signature_vs_plain_decoys")
    if plain and plain.get("auc") is not None:
        verb = ("does separate them" if plain.get("separates")
                else "does not separate them either")
        print(textwrap.fill(
            f"Against ordinary decoys the same score {verb}: AUC "
            f"{plain['auc']:.3f}, Holm p = {plain['p_holm']:.4f} "
            f"(n = {plain['n_positive']} vs {plain['n_negative']}).",
            width=WIDTH, initial_indent="  ", subsequent_indent="  "))
    print(f"  source      {lim['source']}")


def print_board(ans: dict) -> None:
    t, site, rk = ans["target"], ans["site"], ans["ranking"]
    name = t["symbol"] or t["slug"]
    print()
    line = f"{name} ({t['uniprot']})"
    if t["protein_name"]:
        line += f" - {t['protein_name']}"
    print(line)
    print(textwrap.fill(
        f"structure {t['structure']} chain {t['chain']}; site = {site['how_defined']}"
        f"; {site['n_core_residues']} core residues", width=WIDTH))
    print_resolution(ans)
    print(textwrap.fill(
        f"{rk['n_scored']} of {rk['n_rows']} drugs scored, ranked by "
        f"{rk['metric']} under the {site['signature']} signature. Known binders "
        f"land at ranks {rk['known_binder_ranks']} "
        f"(median {rk['median_known_binder_rank']}), enrichment in the top "
        f"quartile {rk['enrichment_at_top_quartile']:.2f}x.",
        width=WIDTH))
    print()
    print(f"{'rank':>4}  {'drug':<26} {'score':>6}  {'core':>7}  role")
    print(RULE)
    for r in ans["matches"]:
        core = f"{r['n_core_engaged']}/{r['n_core']}"
        score = "n/a" if r["score"] is None else f"{r['score']:.3f}"
        print(f"{r['rank']:>4}  {trunc(r['name'], 26):<26} {score:>6}  {core:>7}  "
              f"{ROLE_LABEL.get(r['role'], r['role'])}")
    print(RULE)
    print(f"showing {ans['n_matches_shown']} of {ans['n_matches_total']} scored "
          f"(--top N for more)")
    print(textwrap.fill(
        f"score = {rk['metric']}: of the residues this drug engages, the fraction "
        f"lying in the site core.  core = how many of the site's core residues it "
        f"engages.", width=WIDTH))
    if ans["not_scored"]:
        print()
        print("NOT SCORED (kept, not dropped)")
        for r in ans["not_scored"]:
            label = ROLE_LABEL.get(r["role"], r["role"])
            print(f"  {r['name']} ({label}): {r['status']}")

    print()
    print("WHY THE TOP ROWS RANK WHERE THEY DO - the residues, not just the number")
    for r in ans["matches"][:3]:
        print(f"  #{r['rank']} {r['name']} [{ROLE_LABEL.get(r['role'], r['role'])}] "
              f"- engages {r['n_engaged']} residues, "
              f"{r['n_core_engaged']} of the {r['n_core']} core "
              f"(coverage {r['core_coverage']:.2f})")
        print(wrap_list(r["engaged_core"], "      engaged core: "))
        print(wrap_list(r["missed_core"], "      missed  core: "))
    print()
    print(textwrap.fill("WHAT THIS IS: " + ans["disclaimer"], width=WIDTH))
    print()
    print_limitation(ans["limitation"])
    print(f"  read from   {', '.join(ans['sources'])}")
    print("  no network, no credits, no GPU: every number above was read from those "
          "files.")
    print()


def print_unscreened(ans: dict) -> None:
    t, ex = ans["target"], ans["what_exists"]
    print()
    print(f"{t['symbol'] or t['slug']} ({t['uniprot']})"
          + (f" - {t['protein_name']}" if t["protein_name"] else ""))
    print_resolution(ans)
    print(f"NO DRUG BOARD. This target is {ans['state']}.")
    print()
    print(textwrap.fill("Why: " + ans["why"], width=WIDTH))
    print()
    print("WHAT DOES EXIST")
    p = ex["pocket"]
    print(f"  binding site   {p['site_id']}, {p['n_residues']} residues, "
          f"structure {t['structure']} chain {t['chain']}")
    print(f"                 {p['detector']}")
    print(wrap_list(p["residue_ids"], "                 "))
    print(f"  structures     {ex['structures']}")
    print(f"  designs        {ex['n_boltzgen_designs']} BoltzGen designs on disk")
    gap = "  <- this is the gap"
    zero_cofold = not ex["n_cofold_poses"]
    print(f"  co-folded      {ex['n_cofold_poses']} poses"
          + (gap if zero_cofold else ""))
    print(f"  missing        {', '.join(ans['missing'])}"
          + ("" if zero_cofold else gap))
    print()
    print("TO PRODUCE THE BOARD, run exactly this:")
    print()
    print(textwrap.fill(ans["to_produce_it"]["command"], width=WIDTH - 6,
                        initial_indent="    ", subsequent_indent="        ",
                        break_long_words=False, break_on_hyphens=False)
          .replace("\n", " \\\n"))
    print()
    print(f"  ({ans['to_produce_it']['shape_source']})")
    print(textwrap.fill(ans["to_produce_it"]["requires"], width=WIDTH,
                        initial_indent="  ", subsequent_indent="  "))
    print()
    print(textwrap.fill(
        "Until that has run there is no ranking for this target, and this script "
        "will not invent one.", width=WIDTH))
    print()


def print_explain(ans: dict) -> None:
    d, t = ans["drug"], ans["target"]
    print()
    print(f"{d['name']}  [{d['role_label']}]   on the {t['symbol']} board "
          f"({t['structure']})")
    if d["status"] != "scored":
        print(f"  status: {d['status']} - this drug has no pose to score, so it has "
              f"no rank.")
        print()
        return
    print(f"  engages {d['n_engaged']} residues of {t['symbol']} in its co-folded pose")
    print(wrap_list(d["engaged_residues"], "    "))
    print()
    print("RANK UNDER EACH SIGNATURE  (same metric, three different definitions of "
          "'the site')")
    print(f"  {'signature':<20} {'rank':>8}  {'score':>6}  {'core':>7}  defined by")
    print(RULE)
    for sig, s in ans["by_signature"].items():
        if not s["scored"]:
            print(f"  {sig:<20} {'not scored':>8}")
            continue
        core = f"{len(s['engaged_core'])}/{s['n_core']}"
        mark = " *" if s["is_ranking_signature"] else "  "
        print((f"  {sig:<20} {str(s['rank']) + '/' + str(s['n_ranked']):>8}  "
               f"{s['score']:.3f}  {core:>7}  "
               f"{trunc(s['provenance'], 24)}{mark}").rstrip())
    print(RULE)
    print(f"  * = the signature this board ranks by ({ans['ranked_by']}); "
          f"score = {next(iter(ans['by_signature'].values())).get('score_metric')}")
    print()
    for sig, s in ans["by_signature"].items():
        if not s["scored"]:
            continue
        print(f"  {sig}:")
        print(wrap_list(s["engaged_core"], "      engaged core: "))
        print(wrap_list(s["missed_core"], "      missed  core: "))
    print()
    chem = ans["chemical_similarity"]
    print("CHEMICAL SIMILARITY, for contrast (a different question entirely)")
    if not chem:
        print("  no candidates.json for this target, so no chemical ranking to show.")
    elif chem.get("not_listed"):
        print(textwrap.fill(
            f"this drug is not in the {chem['representation']} top-{chem['top_k']} "
            f"list nor among the known binders ranked in {chem['source']} "
            f"(corpus n = {chem['n_corpus']}), so it has no stored chemical rank.",
            width=WIDTH, initial_indent="  ", subsequent_indent="  "))
    else:
        print(f"  rank {chem['rank']} by {chem['representation']} similarity to "
              f"{chem['reference_ligand']} (corpus n = {chem['n_corpus']}), "
              f"Tanimoto {chem['similarity']:.3f}")
        print(f"  source      {chem['source']}")
    print()
    print(textwrap.fill("WHAT THIS IS: " + ans["disclaimer"], width=WIDTH))
    print()
    print_limitation(ans["limitation"])
    print(f"  read from   {', '.join(ans['sources'])}")
    print()


# -------------------------------------------------------------------------- cli
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="ask.py",
        description="Ask about a disease or a target; get the ranked drugs that "
                    "engage the same binding site. Reads committed results only.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              ask.py --list
              ask.py "colorectal cancer"
              ask.py --target KDR --top 10
              ask.py --target VEGFR2 --json
              ask.py --target KDR --explain lapatinib
              ask.py --target CDK2
            """))
    ap.add_argument("question", nargs="*",
                    help='what to ask about, e.g. "colorectal cancer" or KDR')
    ap.add_argument("--target", help="gene symbol or alias (KDR, VEGFR2, CDK2)")
    ap.add_argument("--list", action="store_true",
                    help="what can be asked about, and the state of each target")
    ap.add_argument("--top", type=int, default=15, metavar="N",
                    help="how many ranked rows to print (default 15; 0 = all)")
    ap.add_argument("--explain", metavar="DRUG",
                    help="one drug's detail: residues engaged and missed, its rank "
                         "under each signature, its chemical-similarity rank")
    ap.add_argument("--json", action="store_true", help="machine-readable answer")
    ap.add_argument("--debug", action="store_true",
                    help="show the traceback on an unexpected error")
    args = ap.parse_args(argv)

    try:
        registry = build_registry()

        if args.list:
            ans = answer_list(registry)
            print(json.dumps(ans, indent=2)) if args.json else print_list(ans)
            return 0

        question = args.target or " ".join(args.question)
        if not question:
            if args.explain:
                boards = [e for e in registry if e["state"] == "board"]
                if len(boards) == 1:
                    entry = boards[0]
                    question = entry["gene"] or entry["slug"]
                else:
                    raise AskError("--explain needs a target",
                                   "add --target KDR, or run ask.py --list")
            else:
                ap.print_help()
                print()
                raise AskError("ask me something",
                               'try: ask.py --list, or ask.py "colorectal cancer"')

        entry, match = resolve(question, registry)

        if entry["state"] == "board_unusable":
            raise AskError(
                f"the drug board for {entry['gene'] or entry['slug']} exists but "
                f"cannot be used: {rel(entry['board_path'])} is "
                f"{entry['board_problem']}",
                "this is a file problem, not a result: restore that file from git "
                "(git checkout -- " + rel(entry["board_path"]) + ") rather than "
                "reading this as 'never screened'")

        if args.explain:
            if entry["state"] != "board":
                raise AskError(
                    f"{entry['gene'] or entry['slug']} has no drug board, so no drug "
                    f"on it can be explained",
                    f"run: ask.py --target {entry['gene'] or entry['slug']}  "
                    f"(it prints how to produce the board)")
            ans = answer_explain(entry, args.explain, question, match)
            print(json.dumps(ans, indent=2)) if args.json else print_explain(ans)
            return 0

        if entry["state"] == "board":
            ans = answer_board(entry, question, max(0, args.top), match)
            print(json.dumps(ans, indent=2)) if args.json else print_board(ans)
            return 0

        if entry["state"] == "prepared":
            ans = answer_unscreened(entry, question, registry, match)
            print(json.dumps(ans, indent=2)) if args.json else print_unscreened(ans)
            return 0

        raise AskError(
            f"{entry['slug']} has no target/pocket.json, so it is not askable yet",
            "run the prep stage of scripts/autorepurpose.py for it")

    except AskError as exc:
        if args.json:
            print(json.dumps({"ok": False, "error": exc.message, "fix": exc.fix},
                             indent=2))
        else:
            print(f"\nask.py: {exc.message}", file=sys.stderr)
            if exc.fix:
                print(f"        {exc.fix}", file=sys.stderr)
            print(file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # never hand a user a traceback
        if args.debug:
            raise
        print(f"\nask.py: unexpected error: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        print("        re-run with --debug for the traceback\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
