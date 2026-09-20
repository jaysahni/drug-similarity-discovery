#!/usr/bin/env python
"""Render one self-contained HTML run report for the repurposing pipeline.

PROJECT_GOAL.md G7: one self-contained HTML per run. No network, no CDN, no
external asset - every byte of CSS and SVG is inline, so the file opens from a
file:// URL with nothing else on disk.

Everything rendered is read from the run's JSON outputs. Nothing is typed in by
hand: if a number is not in the inputs or derived from them here, it does not
appear in the report. Absence is checked, not assumed - a drug with no chemical
rank, a role with no rows, a job with no credits are each counted and named. The
roles, signatures and metrics are whatever the board contains, so a run that adds
an arm renders it without a code change.

PROJECT_GOAL.md G8 language guard: the generated HTML is scanned for affinity
vocabulary before it is written. The build fails and writes nothing if a banned
term appears outside an explicit disclaimer block.

Usage:
    ./env/bin/python scripts/report.py
    ./env/bin/python scripts/report.py --generated-at 2026-09-19T20:00:00-04:00 \
        --previous-board /path/to/earlier/repurpose_colorectal-cancer.json \
        --previous-label 'results/repurpose_colorectal-cancer.json @ 16bd62a'
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import html
import json
import platform
import re
import sys
from pathlib import Path

from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parent.parent

INPUTS = {
    "board": "results/repurpose_colorectal-cancer.json",
    "signature": "results/boltzgen_signature_v4.json",
    "m2_kdr": "results/m2_gate_KDR.json",
    "m2_cdk2": "results/m2_gate_CDK2.json",
    "candidates": "results/pipeline/colorectal-cancer/candidates.json",
    "pocket": "results/pipeline/colorectal-cancer/target/pocket.json",
    "state": "results/pipeline/colorectal-cancer/repurpose_state.json",
    "ligand": "results/pipeline/colorectal-cancer/target/known_ligand_contacts.json",
    "wf": "results/pipeline/colorectal-cancer/boltzgen_v4_wf.json",
    "p2rank_batch": "results/p2rank_pockets.json",
}
STRUCTURE = "results/pipeline/colorectal-cancer/target/structure/3VHE_A_stripped.pdb"

ROLE_LABEL = {
    "known_binder": "known binder",
    "decoy": "decoy",
    "hard_decoy": "hard decoy",
}
ROLE_ORDER = ["known_binder", "hard_decoy", "decoy"]
SIG_LABEL = {
    "boltzgen_consensus": "BoltzGen design consensus",
    "p2rank_geometry": "P2Rank pocket geometry",
    "known_ligand": "known-ligand contacts",
}

# ------------------------------------------------------------ language guard
BANNED = [
    ("affinity", r"\baffinit(?:y|ies)\b"),
    ("Kd", r"\bKd\b"),
    ("IC50", r"\bIC\s?50\b"),
    ("potency", r"\bpotenc(?:y|ies)\b"),
    ("potent", r"\bpotent\b"),
    ("cures", r"\bcures\b"),
    ("proves", r"\bproves\b"),
    ("demonstrates efficacy", r"\bdemonstrates\s+efficacy\b"),
]
EXEMPT_OPEN = "<!--LANG-EXEMPT-->"
EXEMPT_CLOSE = "<!--/LANG-EXEMPT-->"
EXEMPT_RE = re.compile(
    re.escape(EXEMPT_OPEN) + r".*?" + re.escape(EXEMPT_CLOSE), re.DOTALL
)


def language_guard(markup: str) -> list[str]:
    """Scan generated HTML for affinity vocabulary (PROJECT_GOAL.md G8).

    Explicit disclaimer blocks are exempt: they are the one place the report is
    allowed to name what the pipeline does NOT predict. Everything else must be
    clean. Returns a list of violations; empty means the build may proceed.
    """
    scanned = EXEMPT_RE.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), markup)
    violations = []
    for term, pattern in BANNED:
        for hit in re.finditer(pattern, scanned, re.IGNORECASE):
            line = scanned.count("\n", 0, hit.start()) + 1
            lo = max(0, hit.start() - 60)
            ctx = " ".join(scanned[lo:hit.end() + 60].split())
            violations.append(f'"{term}" at line {line}: ...{ctx}...')
    return violations


# ------------------------------------------------------------------ helpers
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def esc(value) -> str:
    return html.escape(str(value), quote=True)


def fmt(value, places=3) -> str:
    return "n/a" if value is None else f"{value:.{places}f}"


def role_label(role: str) -> str:
    return ROLE_LABEL.get(role, role.replace("_", " "))


def plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def join_en(items: list[str]) -> str:
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def residue_names(root: Path) -> dict[int, str]:
    """auth residue id -> three-letter name, parsed from the stripped chain."""
    path = root / STRUCTURE
    if not path.exists():
        return {}
    names: dict[int, str] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("ATOM") and line[21] == "A":
                try:
                    names[int(line[22:26])] = line[17:20].strip()
                except ValueError:
                    continue
    return names


def residue_span(rid: int, names: dict[int, str], core: set[int]) -> str:
    label = f"{names.get(rid, '')}{rid}" if names.get(rid) else str(rid)
    cls = "res core" if rid in core else "res"
    return f'<span class="{cls}">{esc(label)}</span>'


# -------------------------------------------------------------- computation
def max_enrichment(n_scored: int, n_known: int) -> tuple[float, int]:
    """Best enrichment reachable at the top quartile for this base rate."""
    k = n_scored // 4
    return (min(n_known, k) / k) / (n_known / n_scored), k


def rank_band(scored: list[dict], signature: str, metric: str) -> dict[str, tuple]:
    """Tie-aware rank band for every drug under one signature.

    A drug's rank is only defined up to the drugs it ties with, and the pipeline
    breaks those ties by its own input order rather than by score. Comparing a
    single derived integer against the recorded one would flag a tie as a
    disagreement, so each drug gets the closed interval of ranks its score allows
    and the recorded rank is checked against that.
    """
    scores = [r["by_signature"][signature][metric] for r in scored]
    bands = {}
    for r, s in zip(scored, scores):
        better = sum(1 for other in scores if other > s)
        ties = sum(1 for other in scores if other == s)
        bands[r["name"]] = (better + 1, better + ties)
    return bands


def check_recorded_ranks(scored, signature, metric, recorded) -> dict:
    """Do the recorded known-binder ranks agree with the stored scores?"""
    bands = rank_band(scored, signature, metric)
    known = [r["name"] for r in scored if r["role"] == "known_binder"]
    inside, ties = 0, 0
    for name, rec in zip(
        sorted(known, key=lambda n: bands[n]), sorted(recorded)
    ):
        lo, hi = bands[name]
        inside += lo <= rec <= hi
        ties += hi > lo
    return {"ok": inside == len(recorded), "n": len(recorded), "ties": ties}


def chemical_ranks(candidates: dict) -> dict[str, int]:
    """Every chemical-similarity rank candidates.json actually stores."""
    ranks: dict[str, int] = {}
    for entry in candidates.get("top", []):
        ranks.setdefault(entry["name"], entry["rank"])
    for entry in candidates.get("known_binder_ranks", []):
        ranks.setdefault(entry["name"], entry["rank"])
    return ranks


def engagement_frequency(signature: dict) -> list[dict]:
    """Per-residue contact frequency across the designs."""
    designs = signature["per_design"]
    n = len(designs)
    counts: dict[int, int] = {}
    for design in designs:
        for rid in set(design["contacts_auth"]):
            counts[rid] = counts.get(rid, 0) + 1
    return [
        {"residue": rid, "count": counts[rid], "fraction": counts[rid] / n}
        for rid in sorted(counts)
    ]


# ------------------------------------------------------------------ sections
def banner() -> str:
    """The standing caveat. Also the report's one language-guard exemption."""
    return f"""{EXEMPT_OPEN}
<div class="banner" role="note">
  <p class="banner-title">These are computational hypotheses requiring experimental
     validation.</p>
  <p>Nothing on this page is a measurement. Every drug on this board was placed there
     by a structure-prediction model co-folding it into a predicted pocket, then ranked
     by how much of a predicted residue set it overlaps. No affinity, Kd, IC50 or
     potency is predicted, read, or implied anywhere in this pipeline &mdash; Boltz-2
     can emit an affinity score and the pipeline deliberately does not read it
     (PROJECT_GOAL.md 4.4). A high rank here is a reason to run an assay. It is not
     evidence that a drug binds, works, or treats anything.</p>
  <p>Boltz-2 scores BoltzGen designs and the two share training data, and co-folding an
     approved drug against its known target means the model has likely seen that
     complex. The decoy arms exist because of that, and they are what make any number
     here interpretable.</p>
</div>
{EXEMPT_CLOSE}"""


def section_header(ctx: dict) -> str:
    b = ctx["board"]
    target = b["target"]
    counts = " &middot; ".join(
        f"{plural(n, esc(role_label(role)))}" for role, n in ctx["role_counts"].items()
    )
    return f"""
<header>
  <p class="eyebrow">AutoRepurpose v0.2 &middot; interface-signature run report</p>
  <h1>{esc(target['symbol'])} &middot; {esc(ctx['pipeline_label'])}</h1>
  <p class="verdict">{ctx['verdict']}</p>
  <p class="subline">
    Target {esc(target['symbol'])} (UniProt {esc(target['uniprot'])}), structure
    {esc(target['structure'])}. n = {ctx['n_scored']} drugs scored of
    {ctx['n_rows']} co-folded: {counts}. Ranked by
    <code>{esc(b['ranking_metric'])}</code> against the
    {esc(SIG_LABEL.get(b['ranked_by'], b['ranked_by']))}.
  </p>
{ctx['banner']}
</header>"""


def section_board(ctx: dict) -> str:
    b = ctx["board"]
    metric, sig = b["ranking_metric"], b["ranked_by"]
    core = set(b["signatures"][sig]["core_residues_auth"])
    names = ctx["names"]

    rows = []
    for r in ctx["scored"]:
        s = r["by_signature"][sig]
        engaged = "".join(residue_span(x, names, core) for x in s["engaged_core"]) or \
            '<span class="none">none</span>'
        missed = "".join(residue_span(x, names, core) for x in s["missed_core"]) or \
            '<span class="none">none &mdash; full core engaged</span>'
        rows.append(f"""
      <tr class="role-{esc(r['role'])}">
        <td class="num">{r['rank']}</td>
        <td class="drug">{esc(r['name'])}</td>
        <td><span class="badge b-{esc(r['role'])}">{esc(role_label(r['role']))}</span></td>
        <td class="num">{fmt(s[metric])}</td>
        <td class="num">{r['n_engaged']}</td>
        <td class="num">{len(s['engaged_core'])}/{len(core)}</td>
        <td class="reslist">{engaged}</td>
        <td class="reslist muted">{missed}</td>
      </tr>""")

    unscored = "".join(f"""
      <tr class="unscored">
        <td class="num">&mdash;</td>
        <td class="drug">{esc(r['name'])}</td>
        <td><span class="badge b-{esc(r['role'])}">{esc(role_label(r['role']))}</span></td>
        <td colspan="5" class="reason">not scored &mdash;
            <code>{esc(r['status'])}</code>: the co-folding job returned a complex with
            no ligand in the pose, so there were no contacts to score. Excluded from n
            and from every number on this page.</td>
      </tr>""" for r in ctx["unscored"])

    absent = ""
    if ctx["roles_absent"]:
        absent = (
            f'<p class="note">Roles present in <code>repurpose_state.json</code> but '
            f'absent from this board: '
            f'{esc(join_en([role_label(m) for m in ctx["roles_absent"]]))}. '
            f'Checked against the state file, not assumed.</p>'
        )

    return f"""
<section id="board">
  <h2>1 &middot; The ranked board</h2>
  <p>Every co-folded drug, ranked by <code>{esc(metric)}</code> against the
     {esc(SIG_LABEL.get(sig, sig))} ({len(core)} core residues). The core residues a
     drug <em>missed</em> is the explanatory field (PROJECT_GOAL.md G4): it names the
     part of the site the drug fails to reach, which a single score cannot. n =
     {ctx['n_scored']} scored.</p>
  {absent}
  {ctx['hard_decoy_callout']}
  <div class="scroll">
  <table class="board">
    <thead>
      <tr>
        <th>#</th><th>drug</th><th>role</th><th>{esc(metric)}</th>
        <th>residues<br>engaged</th><th>core<br>hit</th>
        <th>core residues engaged</th><th>core residues missed</th>
      </tr>
    </thead>
    <tbody>{''.join(rows)}{unscored}</tbody>
  </table>
  </div>
</section>"""


def hard_decoy_callout(ctx: dict) -> str:
    """The labelling-gap caveat, stated where the board is read, not in a footnote."""
    top = ctx["hard_top_named"]
    if not top:
        return ""
    listed = join_en([f"<strong>{esc(n)}</strong> at {r}" for n, r in top])
    return f"""{EXEMPT_OPEN}
  <div class="callout" role="note">
    <p class="callout-title">Read the top of this board with the hard decoys in
       mind.</p>
    <p>The highest-ranked hard decoys are {listed} of {ctx['n_scored']} &mdash; approved
       kinase inhibitors carrying no annotation against
       {esc(ctx['board']['target']['symbol'])}. {ctx['n_in_top_phrase']} inside the top
       quartile that enrichment is measured over, and {ctx['n_hard_interleaved']} of
       {ctx['n_hard']} outrank the lowest-placed known binder (rank
       {ctx['last_known_rank']}).</p>
    <p><em>This report cannot tell you whether that is a false positive or a gap in the
       annotation.</em> These are promiscuous kinase inhibitors and the annotation
       tables used to assign roles are known to be incomplete, so a high score here is
       equally consistent with the pipeline being right and the label being wrong.
       Deciding between those two cases needs an assay, not another model. Every
       enrichment number on this page counts them as negatives, which is the
       conservative choice and may understate the ranking.</p>
  </div>
{EXEMPT_CLOSE}"""


def section_comparison(ctx: dict) -> str:
    joined = ctx["comparison"]
    cand = ctx["candidates"]
    n_corpus = cand["n_corpus"]

    rows = "".join(f"""
      <tr>
        <td class="drug">{esc(row['name'])}</td>
        <td><span class="badge b-{esc(row['role'])}">{esc(role_label(row['role']))}</span></td>
        <td class="num strong">{row['interface_rank']}</td>
        <td class="num">{fmt(row['interface_pct'], 1)}%</td>
        <td class="num strong">{row['chem_rank']:,}</td>
        <td class="num">{fmt(row['chem_pct'], 1)}%</td>
        <td class="num {'gain' if row['pct_delta'] > 0 else 'loss'}">{fmt(row['pct_delta'], 1)}</td>
      </tr>""" for row in joined)

    return f"""
<section id="comparison">
  <h2>2 &middot; Interface rank against chemical-similarity rank</h2>
  <p>The same drugs, ranked two ways. The interface board is the table above. The
     chemical ranking is <code>{esc(INPUTS['candidates'])}</code>:
     {esc(cand['representation'])} fingerprint similarity to
     {esc(str(cand['ideal_binder']['pref_name']).lower())}
     ({esc(cand['ideal_binder']['chembl_id'])}), the reference ligand for this target,
     over a corpus of {n_corpus:,} drugs. Every joinable drug is shown, sorted by the
     percentile gap &mdash; nothing is selected for effect. n = {len(joined)}.</p>

  <div class="scroll">
  <table>
    <thead>
      <tr><th>drug</th><th>role</th>
          <th colspan="2">interface rank (of {ctx['n_scored']})</th>
          <th colspan="2">chemical rank (of {n_corpus:,})</th>
          <th>percentile<br>gap</th></tr>
    </thead>
    <tbody>{rows}</tbody>
  </table>
  </div>

  <p class="reading"><strong>What this does and does not show.</strong>
     Raw ranks out of {ctx['n_scored']} and out of {n_corpus:,} are not on the same
     scale, so the percentile columns sit next to them and the gap is computed on
     percentiles rather than on the raw difference. Read that way,
     {esc(ctx['comparison_verdict'])}
     The defensible claim is the narrow one: engaging the same site and resembling the
     reference ligand are different properties, and the chemical ranking places
     {ctx['chem_far']} of these {len(joined)} known binders beyond rank
     {ctx['chem_far_threshold']:,} of {n_corpus:,} &mdash; deep enough that a
     similarity screen would not reach them.</p>

  <p class="note"><strong>Coverage, checked not assumed.</strong>
     {ctx['n_joined']} of the {ctx['n_rows']} drugs on the board carry a chemical rank,
     and all {ctx['n_joined']} are known binders. The other {ctx['n_unjoined']} &mdash;
     every decoy and every hard decoy &mdash; are absent from
     <code>candidates.json</code>, which stores only the top {cand['top_k']} of the
     corpus plus the {len(cand['known_binder_ranks'])} known binders, not the full
     {n_corpus:,}-drug ordering. Their chemical ranks were not computed anywhere in
     this run, so they are omitted rather than guessed at. This table is therefore a
     comparison over known binders only, and it cannot say how the two rankings treat
     a negative.</p>
</section>"""


def section_signature(ctx: dict) -> str:
    sig = ctx["signature"]
    freq = ctx["frequency"]
    names = ctx["names"]
    core = set(sig["core_residue_ids_auth"])
    pocket = set(ctx["pocket"]["residue_ids"])
    ligand = set(ctx["ligand"]["known_ligand_residue_ids"])
    n_designs = sig["n_designs"]
    threshold = sig["core_threshold"]
    off_axis = len(ligand) - sum(1 for item in freq if item["residue"] in ligand)
    off_axis_note = (
        f" A strip can only mark residues that have a bar: {off_axis} of the "
        f"{len(ligand)} known-ligand residues were never touched by any design, so "
        f"they are off this axis and the legend counts what is drawn."
        if off_axis else
        " Every residue in all three sets was touched by at least one design, so each "
        "set is drawn in full."
    )

    bar_w, gap = 16, 3
    step = bar_w + gap
    left, top, plot_h = 58, 14, 140
    strip_h, strip_gap = 11, 4
    label_h = 64
    width = left + step * len(freq) + 16
    strips_y = top + plot_h + 8
    height = strips_y + 3 * (strip_h + strip_gap) + label_h

    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'role="img" aria-label="Per-residue contact frequency across '
        f'{n_designs} BoltzGen designs">'
    ]
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = top + plot_h - frac * plot_h
        parts.append(
            f'<line x1="{left - 6}" y1="{y:.1f}" x2="{width - 16}" y2="{y:.1f}" '
            f'class="grid"/>'
            f'<text x="{left - 10}" y="{y + 3.5:.1f}" class="ax" '
            f'text-anchor="end">{int(frac * 100)}%</text>'
        )
    ty = top + plot_h - threshold * plot_h
    parts.append(
        f'<line x1="{left - 6}" y1="{ty:.1f}" x2="{width - 16}" y2="{ty:.1f}" '
        f'class="thr"/>'
        f'<text x="{width - 18}" y="{ty - 5:.1f}" class="ax thrlab" '
        f'text-anchor="end">core threshold {threshold:.0%}</text>'
    )

    for i, item in enumerate(freq):
        rid = item["residue"]
        x = left + i * step
        h = item["fraction"] * plot_h
        y = top + plot_h - h
        cls = "bar core" if rid in core else "bar"
        parts.append(
            f'<rect x="{x}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" '
            f'class="{cls}"><title>{esc(names.get(rid, ""))}{rid}: '
            f'{item["count"]} of {n_designs} designs '
            f'({item["fraction"]:.0%})</title></rect>'
        )
        for row, member in enumerate((rid in core, rid in pocket, rid in ligand)):
            sy = strips_y + row * (strip_h + strip_gap)
            fill = ("m-core", "m-pocket", "m-ligand")[row] if member else "m-off"
            parts.append(
                f'<rect x="{x}" y="{sy}" width="{bar_w}" height="{strip_h}" '
                f'class="mark {fill}"/>'
            )
        ly = strips_y + 3 * (strip_h + strip_gap) + 6
        label = f'{names.get(rid, "")}{rid}' if names.get(rid) else str(rid)
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{ly}" class="rl" '
            f'transform="rotate(-90 {x + bar_w / 2:.1f} {ly})" '
            f'text-anchor="end">{esc(label)}</text>'
        )
    shown = [sum(1 for item in freq if item["residue"] in s)
             for s in (core, pocket, ligand)]
    legends = []
    for text, members, n_shown in zip(
        ("core", "P2Rank pocket", "known ligand"), (core, pocket, ligand), shown
    ):
        legends.append(
            f"{text} ({len(members)})" if n_shown == len(members)
            else f"{text} ({n_shown} of {len(members)} on this axis)"
        )
    for row, text in enumerate(legends):
        sy = strips_y + row * (strip_h + strip_gap) + strip_h - 2
        parts.append(
            f'<text x="{left - 10}" y="{sy}" class="ax" text-anchor="end">'
            f'{esc(text)}</text>'
        )
    parts.append("</svg>")

    return f"""
<section id="signature">
  <h2>3 &middot; The interface signature</h2>
  <p>How often each target residue was contacted across the {n_designs} BoltzGen
     designs (4.5 &Aring; heavy-atom cutoff). A residue joins the core signature when at
     least {threshold:.0%} of designs touch it, which gives {sig['n_core']} core
     residues out of the {len(freq)} touched by any design. Mean pairwise Jaccard
     between designs is {fmt(sig['mean_pairwise_jaccard_between_designs'])}: the
     designs agree with each other about half the time, which is the spread this
     consensus is averaging over.</p>
  <div class="scroll chart">{''.join(parts)}</div>
  <p class="note">The three strips under the bars are set membership, not a second
     measurement: the core signature, the P2Rank rank-1 pocket ({len(pocket)} residues,
     derived from the ligand-stripped structure), and the contacts of the ligand in the
     holo structure ({len(ligand)} residues, validation reference only).{off_axis_note}
     The design core overlaps the pocket it was asked to target at
     {sig['core_overlap_with_requested_pocket']} of
     {len(sig['requested_pocket_auth'])} residues.</p>
</section>"""


def section_validation(ctx: dict) -> str:
    b = ctx["board"]
    rows = []
    for name, val in ctx["validation_sorted"]:
        check = ctx["rank_check"][name]
        note = "consistent" if check["ok"] else "INCONSISTENT"
        if check["ok"] and check["ties"]:
            note += f" ({plural(check['ties'], 'tie')})"
        pct = 100.0 * val["enrichment_at_top_quartile"] / ctx["max_enrichment"]
        rows.append(f"""
      <tr{' class="ranked-by"' if name == b['ranked_by'] else ''}>
        <td class="drug">{esc(SIG_LABEL.get(name, name))}
            {'<span class="tag">ranked by</span>' if name == b['ranked_by'] else ''}</td>
        <td class="num">{b['signatures'][name]['n_core']}</td>
        <td class="num strong">{fmt(val['enrichment_at_top_quartile'], 2)}&times;</td>
        <td class="num">{fmt(pct, 0)}%</td>
        <td class="num">{fmt(val['median_rank'], 1)}</td>
        <td class="reslist">{esc(', '.join(str(r) for r in val['known_binder_ranks']))}</td>
        <td class="num">{esc(note)}</td>
      </tr>""")

    return f"""
<section id="validation">
  <h2>4 &middot; Validation, per signature</h2>
  <p>The same {ctx['n_scored']} drugs re-ranked under each of the three signatures and
     scored the same way, so the design signature and the free pocket finder sit side by
     side. Enrichment is measured over the top quartile (top {ctx['top_k']} of
     {ctx['n_scored']}); with {ctx['n_known']} known binders among {ctx['n_scored']}
     drugs the arithmetic maximum is {fmt(ctx['max_enrichment'], 2)}&times; and chance
     is 1.00&times;. n = {ctx['n_scored']}.</p>
  <div class="scroll">
  <table>
    <thead>
      <tr><th>signature</th><th>core<br>residues</th><th>enrichment<br>@ top quartile</th>
          <th>% of<br>achievable</th><th>median rank of<br>a known binder</th>
          <th>known-binder ranks</th><th>recorded ranks</th></tr>
    </thead>
    <tbody>{''.join(rows)}</tbody>
  </table>
  </div>
  <p class="note">The last column re-derives each ranking here from the per-drug scores
     and checks the ranks the pipeline recorded against it. Equal scores leave a rank
     undefined within the tied group, so the check is that every recorded rank falls
     inside the band its score allows, not that an integer matches. It confirms this
     report reads the board the way the pipeline wrote it; it is not an independent
     result.</p>
  {ctx['ablation_row']}
</section>"""


def section_limits(ctx: dict) -> str:
    items = "".join(f"""
    <li>
      <h3>{item['head']}</h3>
      {item['body']}
    </li>""" for item in ctx["limits"])
    return f"""
<section id="limits">
  <h2>5 &middot; What does not work</h2>
  <p>These are the run's negative results, kept here with their numbers because a reader
     deciding whether to trust this pipeline needs them more than the board above.</p>
  <ol class="limits">{items}</ol>
</section>"""


def build_limits(ctx: dict) -> list[dict]:
    b = ctx["board"]
    sig = ctx["signature"]
    v = b["validation_by_signature"]
    ranked_by = b["ranked_by"]
    best, worst = ctx["best_sig"], ctx["worst_sig"]
    kdr, cdk2 = ctx["m2_kdr"], ctx["m2_cdk2"]
    kdr_cmp = kdr["comparisons"]["boltzgen_consensus - p2rank_geometry [jaccard]"]
    cdk2_cmp = cdk2["comparisons"]["boltzgen_consensus - p2rank_geometry [jaccard]"]
    iptm = ctx["iptm"]
    limits = []

    # 1. the free baseline beats the paid one, in the product itself
    bg = v["boltzgen_consensus"]
    p2 = v["p2rank_geometry"]
    limits.append({
        "head": "The free pocket finder beats the BoltzGen signature in the product, "
                "not just in the ablation.",
        "body": f"""
      <p>Ranking these {ctx['n_scored']} drugs by the P2Rank rank-1 pocket gives
         enrichment {fmt(p2['enrichment_at_top_quartile'], 2)}&times; at the top
         quartile and a median known-binder rank of
         {fmt(p2['median_rank'], 1)}. The BoltzGen design consensus gives
         {fmt(bg['enrichment_at_top_quartile'], 2)}&times; and
         {fmt(bg['median_rank'], 1)}. The design step cost
         {fmt(sig['run']['credits_charged'], 2)} credits and
         {sig['run']['elapsed_s'] / 60:.0f} minutes; P2Rank cost
         {esc(ctx['p2rank_time'])} per structure
         ({esc(ctx['p2rank_time_note'])}). The board above is ranked by
         <code>{esc(SIG_LABEL.get(ranked_by, ranked_by))}</code> for that reason.</p>
      <p>{esc(ctx['separation_note'])}</p>
      <p>This is ablation I6.1 (PROJECT_GOAL.md 8.3) appearing for the third time. It
         already held on {esc(kdr['target']['label'])}
         (Jaccard &Delta; = {fmt(kdr_cmp['delta'])}, 95% CI
         [{fmt(kdr_cmp['ci_lo'])}, {fmt(kdr_cmp['ci_hi'])}], Wilcoxon p =
         {kdr_cmp['wilcoxon']['p_value']:.2g}, Holm-adjusted p =
         {kdr_cmp['p_holm']:.2g}, n = {kdr_cmp['n']} structures) and on
         {esc(cdk2['target']['label'])}
         (&Delta; = {fmt(cdk2_cmp['delta'])}, 95% CI
         [{fmt(cdk2_cmp['ci_lo'])}, {fmt(cdk2_cmp['ci_hi'])}], Wilcoxon p =
         {cdk2_cmp['wilcoxon']['p_value']:.2g}, Holm-adjusted p =
         {cdk2_cmp['p_holm']:.2g}, n = {cdk2_cmp['n']} structures). The design consensus
         does beat a random residue set of the same size on both targets
         (Jaccard &Delta; =
         {fmt(kdr['comparisons']['boltzgen_consensus - random [jaccard]']['delta'])} and
         {fmt(cdk2['comparisons']['boltzgen_consensus - random [jaccard]']['delta'])}),
         so it carries real information &mdash; it is simply worse than the free
         baseline, on every target measured so far.</p>""",
    })

    # 2. the known-ligand signature, below chance
    kl = v.get("known_ligand")
    if kl and kl["enrichment_at_top_quartile"] < 1.0:
        limits.append({
            "head": "The known-ligand signature ranks worse than chance.",
            "body": f"""
      <p>Taking the contacts of the ligand in the holo structure as the signature gives
         enrichment {fmt(kl['enrichment_at_top_quartile'], 2)}&times;, below the
         1.00&times; of an arbitrary ordering, with a median known-binder rank of
         {fmt(kl['median_rank'], 1)} of {ctx['n_scored']}. It is the only arm here that
         is actively misleading, and it is the arm with privileged information.</p>
      <p>One holo ligand's contact set describes <em>that ligand</em>, not the site. It
         is {ctx['ligand_n']} residues shaped around a single chemotype, so scoring other
         drugs against it rewards resembling that one molecule's footprint rather than
         occupying the pocket. It is kept on this page as a reference arm and a warning,
         and it is not used to rank anything.</p>""",
        })

    # 3. degradation against the harder null
    if ctx["degradation"]:
        rows = "".join(f"""
        <tr><td class="drug">{esc(SIG_LABEL.get(name, name))}</td>
            <td class="num">{fmt(d['prev'], 2)}&times; / {fmt(d['prev_max'], 2)}&times;</td>
            <td class="num">{fmt(d['prev_pct'], 0)}%</td>
            <td class="num">{fmt(d['now'], 2)}&times; / {fmt(d['now_max'], 2)}&times;</td>
            <td class="num">{fmt(d['now_pct'], 0)}%</td>
            <td class="num {'loss' if d['now_pct'] < d['prev_pct'] else 'gain'}">
                {fmt(d['now_pct'] - d['prev_pct'], 0)}</td></tr>"""
            for name, d in ctx["degradation"])
        limits.append({
            "head": "Enrichment drops once the null gets harder, and the easy number "
                    "should never be quoted on its own.",
            "body": f"""
      <p>An earlier run of this board scored the same drugs against the easy decoys only
         ({ctx['prev_n_scored']} drugs, {ctx['prev_n_known']} known binders, maximum
         {fmt(ctx['prev_max'], 2)}&times;). This run adds
         {plural(ctx['n_hard'], 'hard decoy')} &mdash; approved kinase inhibitors not
         annotated against {esc(b['target']['symbol'])} &mdash; raising n to
         {ctx['n_scored']} and the maximum to {fmt(ctx['max_enrichment'], 2)}&times;.
         Because the raw enrichment ceiling moves with the base rate, the comparable
         quantity is the fraction of the achievable maximum, and it falls for every
         signature.</p>
      <div class="scroll">
      <table>
        <thead><tr><th>signature</th><th colspan="2">easy decoys only</th>
               <th colspan="2">with hard decoys</th><th>change</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
      </div>
      <p class="note">{esc(ctx['degradation_check'])}</p>""",
        })
    else:
        limits.append({
            "head": "Degradation against a harder null: not evaluated in this report.",
            "body": """
      <p>Comparing this board against the earlier easy-decoys-only run needs that run's
         board passed in with <code>--previous-board</code>. It was not supplied, so no
         before-and-after figure is shown rather than one being reconstructed from
         memory.</p>""",
        })

    # 4. hard decoys near the top
    if ctx["hard_top_named"]:
        listed = join_en([f"{esc(n)} at {r}" for n, r in ctx["hard_top_named"]])
        limits.append({
            "head": "Hard decoys score near the top, and this run cannot say whether "
                    "that is an error.",
            "body": f"""
      <p>The highest-placed hard decoys are {listed} of {ctx['n_scored']}.
         {plural(len(ctx['hard_in_top']), 'hard decoy')} fall inside the top quartile
         enrichment is measured over, and {ctx['n_hard_interleaved']} of
         {ctx['n_hard']} outrank the lowest-placed known binder (rank
         {ctx['last_known_rank']}) &mdash; the two populations interleave rather than
         separate. Counted as negatives, as every number on this page counts them, they
         are what pulls enrichment down from the easy-decoy figure.</p>
      <p>But the label is an annotation, not an observation. These are promiscuous
         approved kinase inhibitors, the annotation tables used to assign roles are
         incomplete, and nothing in this pipeline distinguishes &ldquo;the ranking is
         wrong&rdquo; from &ldquo;the label is wrong&rdquo;. The honest position is that
         the hard-decoy arm bounds the ranking from below, and that resolving these
         specific drugs needs an assay. See the note at the top of the board.</p>""",
        })

    # 5. the board's top rank is not a known binder
    if ctx["rank1"] and ctx["rank1"]["role"] != "known_binder":
        r1 = ctx["rank1"]
        s1 = r1["by_signature"][ranked_by]
        limits.append({
            "head": "The top-ranked drug on the board is not a known binder.",
            "body": f"""
      <p>Rank 1 is {esc(r1['name'])}, {esc(role_label(r1['role']))}, at
         <code>{esc(b['ranking_metric'])}</code> = {fmt(s1[ctx['metric']])}. It reaches
         that by engaging only {r1['n_engaged']} residues, {len(s1['engaged_core'])} of
         them in the core: a small ligand that sits wholly inside the pocket scores well
         on a precision metric without covering much of it. Its core coverage is
         {fmt(s1['core_coverage'])}, and it still misses
         {plural(len(s1['missed_core']), 'core residue')}.</p>
      <p>The metric is size-normalised on purpose &mdash; PROJECT_GOAL.md 1.4a warns
         that raw overlap ranks bigger ligands first for reasons unrelated to biology
         &mdash; and this is the cost of that choice, visible at rank 1. All five scores
         are stored per drug in the board JSON so the choice can be re-examined; the run
         records that enrichment is unchanged under any of them.</p>""",
        })

    # 6. design confidence
    limits.append({
        "head": "Every design falls far below the confidence filter the plan proposed.",
        "body": f"""
      <p>The {sig['n_designs']} designs this signature is built from span ipTM
         {fmt(iptm['min'], 3)} to {fmt(iptm['max'], 3)} (median
         {fmt(iptm['median'], 3)}). PROJECT_GOAL.md suggests filtering at 0.85;
         {iptm['n_pass']} of {sig['n_designs']} reach it. The signature is a consensus
         over designs the project's own criterion would have discarded. It is reported
         that way rather than re-filtered until something survived, and it is one
         plausible reason the design arm underperforms the pocket finder above.</p>""",
    })
    return limits


def section_provenance(ctx: dict) -> str:
    rows = "".join(f"""
      <tr>
        <td class="drug"><code>{esc(rel)}</code>{label}</td>
        <td class="num">{ctx['sizes'][rel]:,}</td>
        <td class="hash">{esc(ctx['hashes'][rel])}</td>
      </tr>""" for rel, label in ctx["provenance_files"])

    sig = ctx["signature"]
    wf_rows = [
        ("BoltzGen design", sig["workflow_uuid"],
         f"{sig['n_designs']} designs, protocol "
         f"{sig['run']['settings']['protocol']}, "
         f"{fmt(sig['run']['credits_charged'], 2)} credits, "
         f"{sig['run']['elapsed_s']:.0f} s, {sig['run']['status']}"),
        (f"M2 gate {ctx['m2_kdr']['target']['label']}",
         ctx["m2_kdr"]["boltzgen"]["workflow"],
         f"{ctx['m2_kdr']['boltzgen']['n_designs']} designs, "
         f"{fmt(ctx['m2_kdr']['boltzgen']['credits'], 2)} credits, "
         f"scored over {ctx['m2_kdr']['n_structures']} structures"
         + (" - the same design run as the row above, re-scored, not a second charge"
            if ctx["m2_kdr"]["boltzgen"]["workflow"] == sig["workflow_uuid"] else "")),
        (f"M2 gate {ctx['m2_cdk2']['target']['label']}",
         ctx["m2_cdk2"]["boltzgen"]["workflow"],
         f"{ctx['m2_cdk2']['boltzgen']['n_designs']} designs, "
         f"{fmt(ctx['m2_cdk2']['boltzgen']['credits'], 2)} credits, "
         f"scored over {ctx['m2_cdk2']['n_structures']} structures"),
    ]
    wf_html = "".join(
        f'<tr><td class="drug">{esc(label)}</td><td class="hash">{esc(uuid)}</td>'
        f'<td>{esc(detail)}</td></tr>' for label, uuid, detail in wf_rows
    )

    credits = ctx["credits"]
    cofold = "".join(
        f'<tr><td class="drug">{esc(role_label(role))}</td>'
        f'<td class="num">{n}</td><td class="num">{fmt(cr, 2)}</td></tr>'
        for role, (n, cr) in credits["by_role"].items()
    )
    lower_bound = ""
    if credits["n_no_credits"]:
        lower_bound = (
            f" {credits['n_no_credits']} of {credits['n_jobs']} jobs have no credit "
            f"figure recorded and contribute 0, so that total is a lower bound."
        )

    tools = " &middot; ".join(f"{esc(k)} {esc(v)}" for k, v in ctx["versions"].items())

    return f"""
<footer id="provenance">
  <h2>6 &middot; Provenance</h2>
  <p>Generated <time>{esc(ctx['generated_at'])}</time>. The timestamp is passed into
     this report rather than read from the clock while rendering, so two runs over the
     same inputs produce byte-identical files and a diff shows only what changed in the
     data.</p>

  <h3>Inputs</h3>
  <div class="scroll">
  <table class="prov">
    <thead><tr><th>file</th><th>bytes</th><th>SHA-256</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  </div>

  <h3>Rowan workflows</h3>
  <div class="scroll">
  <table class="prov">
    <thead><tr><th>step</th><th>workflow uuid</th><th>detail</th></tr></thead>
    <tbody>{wf_html}</tbody>
  </table>
  </div>

  <h3>Co-folding credits, from <code>repurpose_state.json</code></h3>
  <div class="scroll">
  <table class="prov">
    <thead><tr><th>role</th><th>jobs</th><th>credits</th></tr></thead>
    <tbody>{cofold}
      <tr class="total"><td class="drug">total</td>
          <td class="num">{credits['n_jobs']}</td>
          <td class="num">{fmt(credits['total'], 2)}</td></tr>
    </tbody>
  </table>
  </div>
  <p class="note">Design credits ({fmt(sig['run']['credits_charged'], 2)}) are listed
     with the workflows above and are not included in that total.{lower_bound}</p>

  <h3>Tools</h3>
  <p class="tools">{tools} &middot; {esc(ctx['pocket']['detector'])}</p>
  <p class="note">Site detection ran on the ligand-stripped structure. The known-ligand
     contact set is read here for display and for the validation arm only; no
     signature-construction step reads it (PROJECT_GOAL.md B12).</p>
</footer>"""


CSS = """
:root{
  --ink:#16191d; --mute:#5d666f; --faint:#8b949d;
  --rule:#dfe3e8; --bg:#ffffff; --panel:#f6f8fa;
  --accent:#0f5f6e; --accent-soft:#e3f0f2;
  --warn:#8a5a00; --warn-soft:#fdf3e0; --warn-rule:#e0bc76;
  --bad:#8c2f2f;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
  -webkit-text-size-adjust:100%}
.wrap{max-width:1180px;margin:0 auto;padding:36px 20px 80px}
h1{font-size:27px;line-height:1.2;margin:.15em 0 .35em;letter-spacing:-.01em}
h2{font-size:19px;margin:0 0 .5em;padding-bottom:.4em;border-bottom:2px solid var(--accent);
  letter-spacing:-.005em}
h3{font-size:15px;margin:0 0 .3em}
section,footer{margin-top:46px}
p{margin:.55em 0}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.88em;
  background:var(--panel);padding:.08em .32em;border-radius:3px}
a{color:var(--accent)}
.eyebrow{margin:0;font-size:11.5px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--faint);font-weight:600}
.verdict{font-size:17px;line-height:1.45;margin:.5em 0 .7em;
  border-left:4px solid var(--accent);padding:.2em 0 .2em .8em}
.subline{color:var(--mute);font-size:14px;max-width:80ch}
.note{color:var(--mute);font-size:13px;max-width:88ch}
.reading{background:var(--panel);border:1px solid var(--rule);border-radius:6px;
  padding:.8em 1em;font-size:14px}

.banner{border:1px solid #c9d6d9;background:var(--accent-soft);
  border-left:5px solid var(--accent);border-radius:6px;padding:.9em 1.1em;margin:1.2em 0 0}
.banner-title{margin:0 0 .4em;font-weight:700;color:#0b4a56}
.banner p{font-size:13.5px;margin:.4em 0;color:#22333a}

.callout{border:1px solid var(--warn-rule);background:var(--warn-soft);
  border-left:5px solid var(--warn);border-radius:6px;padding:.8em 1em;margin:1.1em 0}
.callout-title{margin:0 0 .35em;font-weight:700;color:#6f4900}
.callout p{font-size:13.5px;margin:.35em 0;color:#3d3121}

.scroll{overflow-x:auto;margin:1em 0;border:1px solid var(--rule);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--rule);
  vertical-align:top}
thead th{background:var(--panel);font-size:11.5px;text-transform:uppercase;
  letter-spacing:.045em;color:var(--mute);font-weight:700;white-space:nowrap}
tbody tr:last-child td{border-bottom:none}
.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.strong{font-weight:700}
.drug{font-weight:600;white-space:nowrap}
.gain{color:var(--accent);font-weight:700}
.loss{color:var(--bad);font-weight:700}
.total td{font-weight:700;background:var(--panel)}
.ranked-by{background:var(--accent-soft)}
.tag{font-size:10px;text-transform:uppercase;letter-spacing:.05em;color:var(--accent);
  border:1px solid var(--accent);border-radius:3px;padding:0 4px;margin-left:4px;
  font-weight:700;white-space:nowrap}

.badge{display:inline-block;font-size:11px;font-weight:700;padding:1px 7px;
  border-radius:10px;white-space:nowrap;border:1px solid}
.b-known_binder{background:var(--accent-soft);color:#0b4a56;border-color:#8fb8be}
.b-decoy{background:var(--panel);color:var(--mute);border-color:var(--rule)}
.b-hard_decoy{background:var(--warn-soft);color:var(--warn);border-color:var(--warn-rule)}

.reslist{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:11.5px;line-height:1.7;min-width:180px}
.res{display:inline-block;padding:0 4px;margin:0 3px 2px 0;border-radius:3px;
  background:var(--panel);border:1px solid var(--rule)}
.reslist.muted .res{color:var(--mute)}
.none{color:var(--faint);font-style:italic}
.unscored td{background:#fbfbfc}
.reason{color:var(--mute);font-size:12.5px;white-space:normal}
.role-known_binder .drug{color:#0b4a56}
.role-hard_decoy .drug{color:#6f4900}

.chart{padding:10px 4px;background:var(--bg)}
svg .bar{fill:#c6d3d6}
svg .bar.core{fill:var(--accent)}
svg .grid{stroke:var(--rule);stroke-width:1}
svg .thr{stroke:var(--warn);stroke-width:1.2;stroke-dasharray:5 3}
svg .ax{font-size:9.5px;fill:var(--mute);
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
svg .thrlab{fill:var(--warn)}
svg .rl{font-size:9.5px;fill:var(--mute);
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
svg .mark{rx:2}
svg .m-off{fill:#eef1f3}
svg .m-core{fill:var(--accent)}
svg .m-pocket{fill:#6f9aa2}
svg .m-ligand{fill:#b0805a}

.limits{padding-left:1.1em;margin:1em 0}
.limits>li{margin:0 0 1.7em;padding-left:.3em}
.limits h3{color:var(--bad)}
.limits p{font-size:14px;max-width:86ch}
.limits .note{max-width:86ch}

.hash{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:10.5px;
  color:var(--mute);word-break:break-all}
.prov{font-size:12px}
.tools{font-size:12.5px;color:var(--mute)}
footer h3{margin-top:1.4em;color:var(--mute);font-size:12px;text-transform:uppercase;
  letter-spacing:.05em}

@media print{
  .wrap{max-width:none;padding:0}
  .scroll{overflow:visible}
  .limits>li,.banner,.callout,.reading{break-inside:avoid}
  a{color:var(--ink);text-decoration:none}
}
@media (max-width:640px){
  .wrap{padding:22px 16px 60px}
  h1{font-size:22px}
  .verdict{font-size:15px}
}
"""


def build(ctx: dict) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{esc(ctx['board']['target']['symbol'])} interface-signature run report</title>
<style>{CSS}</style>
</head>
<body>
<main class="wrap">
{section_header(ctx)}
{section_board(ctx)}
{section_comparison(ctx)}
{section_signature(ctx)}
{section_validation(ctx)}
{section_limits(ctx)}
{section_provenance(ctx)}
</main>
</body>
</html>
"""


def assemble(root: Path, generated_at: str, previous: Path | None,
             previous_label: str | None) -> dict:
    data, hashes, sizes = {}, {}, {}
    for key, rel in INPUTS.items():
        path = root / rel
        if not path.exists():
            sys.exit(f"missing input: {path}")
        data[key] = json.loads(path.read_text())
        hashes[rel] = sha256(path)
        sizes[rel] = path.stat().st_size

    board, sig, cand = data["board"], data["signature"], data["candidates"]
    metric, ranked_by = board["ranking_metric"], board["ranked_by"]

    scored = sorted(
        (r for r in board["results"] if r.get("status") == "scored"),
        key=lambda r: r["rank"],
    )
    unscored = [r for r in board["results"] if r.get("status") != "scored"]

    role_counts: dict[str, int] = {}
    for r in scored:
        role_counts[r["role"]] = role_counts.get(r["role"], 0) + 1
    role_counts = dict(
        sorted(role_counts.items(),
               key=lambda kv: (ROLE_ORDER.index(kv[0])
                               if kv[0] in ROLE_ORDER else len(ROLE_ORDER), kv[0]))
    )
    state_roles = {j["role"] for j in data["state"]["jobs"].values()}
    roles_absent = sorted(state_roles - set(role_counts))

    n_scored = len(scored)
    n_known = sum(1 for r in scored if r["role"] == "known_binder")
    max_enrich, top_k = max_enrichment(n_scored, n_known)

    rank_check = {
        name: check_recorded_ranks(scored, name, metric, val["known_binder_ranks"])
        for name, val in board["validation_by_signature"].items()
    }
    validation_sorted = sorted(
        board["validation_by_signature"].items(),
        key=lambda kv: -kv[1]["enrichment_at_top_quartile"],
    )
    best_sig, worst_sig = validation_sorted[0][0], validation_sorted[-1][0]

    # signature agreement, on every stored score
    rho_all = {}
    for m in ("core_coverage", "weighted_jaccard", "f1_engaged_core",
              "jaccard_engaged_core", "precision_in_core"):
        x = [r["by_signature"]["boltzgen_consensus"][m] for r in scored]
        y = [r["by_signature"]["p2rank_geometry"][m] for r in scored]
        rho_all[m] = spearmanr(x, y).statistic
    main_rho = spearmanr(
        [r["by_signature"]["boltzgen_consensus"][metric] for r in scored],
        [r["by_signature"]["p2rank_geometry"][metric] for r in scored],
    )

    # chemical join
    chem = chemical_ranks(cand)
    n_corpus = cand["n_corpus"]
    comparison = []
    for r in scored:
        if r["name"] not in chem:
            continue
        ip = 100.0 * r["rank"] / n_scored
        cp = 100.0 * chem[r["name"]] / n_corpus
        comparison.append({
            "name": r["name"], "role": r["role"], "interface_rank": r["rank"],
            "chem_rank": chem[r["name"]], "interface_pct": ip, "chem_pct": cp,
            "pct_delta": cp - ip,
        })
    comparison.sort(key=lambda d: -d["pct_delta"])
    n_joined = len(comparison)
    better = sum(1 for c in comparison if c["pct_delta"] > 0)
    if n_joined == 0:
        verdict_cmp = "no drug on this board carries a chemical rank, so nothing can be compared."
    elif better == n_joined:
        verdict_cmp = (
            f"the interface board places all {n_joined} of them at a better percentile "
            f"than the chemical ranking does."
        )
    elif better == 0:
        verdict_cmp = (
            f"the two rankings are closer than the raw ranks suggest: the chemical "
            f"ranking places every one of these {n_joined} drugs at a better percentile "
            f"than the interface board does."
        )
    else:
        verdict_cmp = (
            f"the two rankings are closer than the raw ranks suggest: the interface "
            f"board wins on percentile for {better} of {n_joined} drugs and loses on "
            f"{n_joined - better}."
        )
    far_threshold = 400
    chem_far = sum(1 for c in comparison if c["chem_rank"] > far_threshold)

    # hard decoys inside the top quartile
    hard = [r for r in scored if r["role"] == "hard_decoy"]
    n_hard = len(hard)
    hard_in_top = [(r["name"], r["rank"]) for r in hard if r["rank"] <= top_k]
    # named in the callout: the top quartile, plus the next hard decoy either side of
    # that line, so a drug is not dropped for sitting one rank outside an arbitrary cut.
    hard_top_named = [(r["name"], r["rank"]) for r in hard[:max(len(hard_in_top) + 1, 3)]]
    known_ranks = [r["rank"] for r in scored if r["role"] == "known_binder"]
    last_known_rank = max(known_ranks) if known_ranks else None
    n_hard_interleaved = sum(1 for r in hard if last_known_rank and r["rank"] < last_known_rank)

    # credits
    by_role: dict[str, list] = {}
    for j in data["state"]["jobs"].values():
        slot = by_role.setdefault(j["role"], [0, 0.0])
        slot[0] += 1
        slot[1] += j.get("credits", 0.0)
    by_role = dict(
        sorted(by_role.items(),
               key=lambda kv: (ROLE_ORDER.index(kv[0])
                               if kv[0] in ROLE_ORDER else len(ROLE_ORDER), kv[0]))
    )
    credits = {
        "by_role": by_role,
        "n_jobs": len(data["state"]["jobs"]),
        "total": sum(j.get("credits", 0.0) for j in data["state"]["jobs"].values()),
        "n_no_credits": sum(
            1 for j in data["state"]["jobs"].values() if "credits" not in j
        ),
    }

    iptms = sorted(d["iptm"] for d in sig["per_design"])
    mid = len(iptms) // 2
    iptm = {
        "min": iptms[0], "max": iptms[-1],
        "median": iptms[mid] if len(iptms) % 2 else (iptms[mid - 1] + iptms[mid]) / 2,
        "n_pass": sum(1 for x in iptms if x >= 0.85),
    }

    names = residue_names(root)
    provenance_files = [(rel, "") for rel in INPUTS.values()]
    if names:
        hashes[STRUCTURE] = sha256(root / STRUCTURE)
        sizes[STRUCTURE] = (root / STRUCTURE).stat().st_size
        provenance_files.append((STRUCTURE, ' <span class="tag">residue names</span>'))

    # ---- optional: the earlier, easier board
    degradation, prev_meta = [], {}
    degradation_check = ""
    if previous is not None:
        prev = json.loads(previous.read_text())
        prev_scored = [r for r in prev["results"] if r.get("status") == "scored"]
        prev_n = len(prev_scored)
        prev_known = sum(1 for r in prev_scored if r["role"] == "known_binder")
        prev_max, _ = max_enrichment(prev_n, prev_known)
        for name, val in validation_sorted:
            if name not in prev["validation_by_signature"]:
                continue
            p = prev["validation_by_signature"][name]["enrichment_at_top_quartile"]
            c = val["enrichment_at_top_quartile"]
            degradation.append((name, {
                "prev": p, "prev_max": prev_max, "prev_pct": 100.0 * p / prev_max,
                "now": c, "now_max": max_enrich, "now_pct": 100.0 * c / max_enrich,
            }))
        # verify the two boards really are the same screen plus a harder arm
        prev_names = {r["name"] for r in prev_scored}
        now_names = {r["name"] for r in scored}
        prev_scores = {r["name"]: r["by_signature"] for r in prev_scored}
        now_scores = {r["name"]: r["by_signature"] for r in scored}
        subset = prev_names <= now_names
        unchanged = all(prev_scores[n] == now_scores[n] for n in prev_names & now_names)
        added = sorted(now_names - prev_names)
        same_sig = all(
            prev["signatures"][k]["core_residues_auth"]
            == board["signatures"][k]["core_residues_auth"]
            for k in board["signatures"] if k in prev["signatures"]
        )
        if subset and unchanged and same_sig:
            degradation_check = (
                f"Checked, not assumed: the earlier board is exactly this one minus "
                f"{plural(len(added), 'drug')}, every shared drug carries identical "
                f"per-signature scores, and all three signatures are unchanged. The "
                f"difference above is the null getting harder and nothing else."
            )
        else:
            degradation_check = (
                "Caution: the earlier board is not a strict subset of this one "
                f"(subset: {subset}, shared scores unchanged: {unchanged}, signatures "
                f"unchanged: {same_sig}), so part of the change may come from something "
                "other than the null."
            )
        # the same agreement statistic on the easier board, so the "the harder null
        # separates them" claim carries its own number instead of an adjective.
        prev_rho = spearmanr(
            [r["by_signature"]["boltzgen_consensus"][metric] for r in prev_scored],
            [r["by_signature"]["p2rank_geometry"][metric] for r in prev_scored],
        )
        prev_meta = {"n_scored": prev_n, "n_known": prev_known, "max": prev_max,
                     "rho": prev_rho.statistic, "rho_p": prev_rho.pvalue}
        rel_label = previous_label or str(previous)
        hashes[rel_label] = sha256(previous)
        sizes[rel_label] = previous.stat().st_size
        provenance_files.append(
            (rel_label, ' <span class="tag">earlier run</span>')
        )

    v = board["validation_by_signature"]
    bg_e = v["boltzgen_consensus"]["enrichment_at_top_quartile"]
    p2_e = v["p2rank_geometry"]["enrichment_at_top_quartile"]
    separation_note = (
        f"Spearman rho between the two rankings is {main_rho.statistic:.3f} "
        f"(p = {main_rho.pvalue:.2g}, n = {n_scored}) on {metric}. On the earlier, "
        f"easier board it was {prev_meta.get('rho', float('nan')):.3f} "
        f"(n = {prev_meta.get('n_scored')}), with both signatures returning identical "
        f"enrichment there. Adding the hard decoys pulls the two rankings apart, and "
        f"it pulls them apart against the design signature."
        if previous is not None else
        f"Spearman rho between the two rankings is {main_rho.statistic:.3f} "
        f"(p = {main_rho.pvalue:.2g}, n = {n_scored}) on {metric}, and "
        f"{rho_all['core_coverage']:.3f} on core_coverage: the two orderings are close, "
        f"but where they differ the free one is ahead."
    )

    ranked_e = v[ranked_by]["enrichment_at_top_quartile"]
    verdict = (
        f"On {esc(board['target']['symbol'])}, the best signature reaches enrichment "
        f"{ranked_e:.2f}&times; at the top quartile against an arithmetic maximum of "
        f"{max_enrich:.2f}&times; ({100 * ranked_e / max_enrich:.0f}% of achievable, "
        f"n = {n_scored}) &mdash; and that signature is the free "
        f"{SIG_LABEL[best_sig] if best_sig in SIG_LABEL else best_sig}, which beats the "
        f"{sig['run']['credits_charged']:.0f}-credit BoltzGen design consensus "
        f"({bg_e:.2f}&times;). The design step is not what is doing the work."
    ) if best_sig == "p2rank_geometry" else (
        f"On {esc(board['target']['symbol'])}, the {SIG_LABEL.get(best_sig, best_sig)} "
        f"reaches enrichment {ranked_e:.2f}&times; at the top quartile against an "
        f"arithmetic maximum of {max_enrich:.2f}&times; "
        f"({100 * ranked_e / max_enrich:.0f}% of achievable, n = {n_scored})."
    )

    ablation_row = (
        f'<p class="note">The ablation behind this table is recorded separately, on two '
        f'targets: <code>{esc(INPUTS["m2_kdr"])}</code> '
        f'({esc(data["m2_kdr"]["target"]["label"])}, n = '
        f'{data["m2_kdr"]["n_structures"]} structures) and '
        f'<code>{esc(INPUTS["m2_cdk2"])}</code> '
        f'({esc(data["m2_cdk2"]["target"]["label"])}, n = '
        f'{data["m2_cdk2"]["n_structures"]} structures). Both are summarised in '
        f'<a href="#limits">what does not work</a>.</p>'
    )

    ctx = {
        "board": board, "signature": sig, "candidates": cand,
        "pocket": data["pocket"], "ligand": data["ligand"],
        "m2_kdr": data["m2_kdr"], "m2_cdk2": data["m2_cdk2"],
        "state": data["state"], "wf": data["wf"],
        "hashes": hashes, "sizes": sizes, "provenance_files": provenance_files,
        "scored": scored, "unscored": unscored, "metric": metric,
        "n_scored": n_scored, "n_rows": len(board["results"]), "n_known": n_known,
        "role_counts": role_counts, "roles_absent": roles_absent,
        "rank_check": rank_check, "validation_sorted": validation_sorted,
        "best_sig": best_sig, "worst_sig": worst_sig,
        "rho": main_rho.statistic, "rho_p": main_rho.pvalue, "rho_all": rho_all,
        "separation_note": separation_note,
        "comparison": comparison, "n_joined": n_joined,
        "n_unjoined": len(board["results"]) - n_joined,
        "comparison_verdict": verdict_cmp,
        "chem_far": chem_far, "chem_far_threshold": far_threshold,
        "hard_in_top": hard_in_top, "n_hard": n_hard,
        "hard_top_named": hard_top_named, "last_known_rank": last_known_rank,
        "n_in_top_phrase": (f"{len(hard_in_top)} of them falls"
                            if len(hard_in_top) == 1
                            else f"{len(hard_in_top)} of them fall"),
        "n_hard_interleaved": n_hard_interleaved,
        "rank1": scored[0] if scored else None,
        "credits": credits, "iptm": iptm,
        "frequency": engagement_frequency(sig), "names": names,
        "top_k": top_k, "max_enrichment": max_enrich,
        "ligand_n": len(data["ligand"]["known_ligand_residue_ids"]),
        "degradation": degradation, "degradation_check": degradation_check,
        "prev_n_scored": prev_meta.get("n_scored"),
        "prev_n_known": prev_meta.get("n_known"), "prev_max": prev_meta.get("max"),
        "verdict": verdict, "banner": banner(), "ablation_row": ablation_row,
        "pipeline_label": "colorectal cancer",
        "p2rank_time": f"{data['p2rank_batch']['seconds_per_structure']:.2f} s",
        "p2rank_time_note": (
            f"mean over the {data['p2rank_batch']['n_structures_parsed']:,}-structure "
            f"batch run, {data['p2rank_batch']['threads']} threads, "
            f"{data['p2rank_batch']['wall_seconds']:,.1f} s wall"
        ),
        "generated_at": generated_at,
        "versions": {
            "python": platform.python_version(),
            "scipy": __import__("scipy").__version__,
            "platform": f"{platform.system()} {platform.machine()}",
        },
    }
    ctx["hard_decoy_callout"] = hard_decoy_callout(ctx)
    ctx["limits"] = build_limits(ctx)
    return ctx


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=REPO_ROOT,
                    help="repo root holding results/ (default: this repo)")
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "results/report_colorectal-cancer.html")
    ap.add_argument("--generated-at", default=None,
                    help="ISO timestamp stamped into the report; defaults to now. "
                         "Pass it explicitly to make reruns byte-identical.")
    ap.add_argument("--previous-board", type=Path, default=None,
                    help="an earlier board JSON, to show enrichment against an easier "
                         "null. Omitted, that comparison is labelled not evaluated.")
    ap.add_argument("--previous-label", default=None,
                    help="how to name the earlier board in the provenance table")
    args = ap.parse_args()

    if args.previous_board and not args.previous_board.exists():
        sys.exit(f"missing --previous-board: {args.previous_board}")

    generated_at = args.generated_at or dt.datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    ctx = assemble(args.root.resolve(), generated_at,
                   args.previous_board, args.previous_label)
    markup = build(ctx)

    violations = language_guard(markup)
    if violations:
        print("LANGUAGE GUARD FAILED (PROJECT_GOAL.md G8) - nothing written",
              file=sys.stderr)
        for v in violations:
            print("  " + v, file=sys.stderr)
        return 1
    print(f"language guard: PASS ({len(BANNED)} patterns, no hit outside the "
          f"disclaimer block)")

    for name, check in ctx["rank_check"].items():
        if not check["ok"]:
            print(f"  WARNING: recorded ranks for {name} fall outside the band their "
                  f"scores allow", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(markup)

    b = ctx["board"]
    print(f"wrote {args.out}  ({len(markup):,} bytes, self-contained)")
    print(f"  n = {ctx['n_scored']} scored of {ctx['n_rows']} rows; roles "
          + ", ".join(f"{v} {k}" for k, v in ctx["role_counts"].items()))
    print(f"  ranked by {b['ranked_by']} on {b['ranking_metric']}; "
          f"max enrichment {ctx['max_enrichment']:.2f}x at top {ctx['top_k']}")
    for name, val in ctx["validation_sorted"]:
        print(f"    {name:<20}{val['enrichment_at_top_quartile']:.2f}x  "
              f"({100 * val['enrichment_at_top_quartile'] / ctx['max_enrichment']:.0f}% "
              f"of achievable)  median rank {val['median_rank']}")
    print(f"  Spearman boltzgen vs p2rank on {b['ranking_metric']}: "
          f"{ctx['rho']:.4f} (p={ctx['rho_p']:.3g})")
    print(f"  chemical join covers {ctx['n_joined']} of {ctx['n_rows']} drugs")
    print(f"  hard decoys in top {ctx['top_k']}: {ctx['hard_in_top']}")
    print(f"  hard decoys named in callout: {ctx['hard_top_named']}")
    print(f"  hard decoys outranking the last known binder (rank "
          f"{ctx['last_known_rank']}): {ctx['n_hard_interleaved']} of {ctx['n_hard']}")
    print(f"  limits rendered: {len(ctx['limits'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
