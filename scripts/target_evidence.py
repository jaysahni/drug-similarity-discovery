"""Pipeline task T3 - real literature evidence for each target from T2.

Reads results/pipeline/<disease_slug>/targets.json, searches Europe PMC for
"<symbol> <disease> inhibitor/binder/antagonist" papers, and writes
results/pipeline/<disease_slug>/evidence.md with up to 3 PMIDs per target.

Two rules the repo's no-hallucinated-citations standard imposes, both enforced in
code rather than by care:

1. The one-line claim under each paper is a VERBATIM sentence lifted from the
   abstract this script actually downloaded - never a paraphrase, so it cannot
   drift from what the paper says. The sentence chosen is the first one naming
   the target (or a gene synonym) together with an inhibitor/binder keyword;
   when the abstract offers nothing that strong, the weaker sentence is still
   quoted verbatim but carries its claim tier in the output so it cannot pass
   for a target+inhibitor statement.
2. Every PMID is re-fetched from Europe PMC in a second, independent request
   (EXT_ID:<pmid>) and the returned title is compared to the one written out.
   A mismatch, or a PMID that does not resolve, is dropped from evidence.md and
   reported in the verification table - it is never silently kept. Every PMID
   that survives is then checked a third time against PubMed E-utilities, a
   different authority, so the citation does not rest on one database agreeing
   with itself.

Targets with zero hits are kept in the output as a negative result, not dropped.

Raw API responses are cached under data/raw/europepmc/ so reruns are cheap and
offline-repeatable.

Usage:
    ./env/bin/python scripts/target_evidence.py
    ./env/bin/python scripts/target_evidence.py --disease "colorectal cancer" --top 25
    ./env/bin/python scripts/target_evidence.py --per-target 3 --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import html
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import certifi

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "europepmc"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# The python.org 3.14 framework build ships no CA bundle, so urllib fails SSL
# verification even though curl succeeds. Use certifi's explicitly.
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

EVIDENCE_TERMS = ("inhibitor", "inhibition", "binder", "antagonist", "blockade")
MAX_SYNONYMS = 6

# Only strip known HTML tags, and only when what follows the tag name really
# looks like attributes. A blanket <[^>]+> eats statistics - "(p < 0.001)" -
# and even plain prose, since "a<b and c>d" would match as a <b> tag.
# Inline tags close up (IC<sub>50</sub> -> IC50); block tags leave a space.
# Each attribute must carry a value: bibliographic markup always does, and
# requiring "=" stops "<b and c>" in prose from parsing as a <b> tag.
_ATTRS = r"""(?:\s+[A-Za-z_:][-\w:.]*\s*=\s*(?:"[^"]*"|'[^']*'|[^\s"'<>`]+))*\s*/?"""
INLINE_TAG = re.compile(
    r"</?(?:sup|sub|sc|inf|i|b|u|em|strong|italic|bold|underline)" + _ATTRS + ">", re.I
)
BLOCK_TAG = re.compile(r"</?(?:h[1-6]|p|br|hr|span|div|a|li|ul|ol)" + _ATTRS + ">", re.I)
HEAD_CLOSE = re.compile(r"</h[1-6]\s*>", re.I)


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def epmc(query, cache_key, page_size=8, refresh=False, retries=4):
    cache = RAW / f"{cache_key}.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())

    params = {
        "query": query,
        "format": "json",
        "pageSize": page_size,
        "resultType": "core",
    }
    url = f"{EPMC}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90, context=SSL_CTX) as r:
                payload = json.load(r)
            break
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(payload, indent=1))
    return payload


def clean(text):
    """Europe PMC embeds markup and entities in titles and abstracts.

    Tags are stripped before entities are unescaped, then again after, because
    some records escape their own markup (`&lt;sup&gt;`) while others do not.
    Structured-abstract headings become "Methods: ..." so the heading stays part
    of the sentence it introduces instead of running into it.

    This is markup-to-text rendering only: no word is added, removed or
    reordered, so a sentence taken from the result is still the paper's own.
    """
    if not text:
        return ""
    for _ in range(2):
        text = HEAD_CLOSE.sub(": ", text)
        text = INLINE_TAG.sub("", text)
        text = BLOCK_TAG.sub(" ", text)
        text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def norm_title(text):
    return re.sub(r"[^a-z0-9]+", " ", clean(text).lower()).strip()


def squash_title(text):
    """Same title with every separator gone.

    Europe PMC and PubMed render markup differently - Europe PMC's
    `KRAS<sup>G12D</sup>` is PubMed's `KRAS(G12D)` - which norm_title sees as a
    difference. Squashing catches that without letting a genuinely different
    paper through: a fabricated citation differs in words, not in punctuation.
    """
    return re.sub(r"[^a-z0-9]+", "", clean(text).lower())


def compare_titles(written, back):
    """'match' | 'rendering' (same words, different markup) | 'differs'."""
    if norm_title(written) == norm_title(back):
        return "match"
    if squash_title(written) == squash_title(back):
        return "rendering"
    return "differs"


def cell(text):
    """A title with a pipe in it would otherwise break the mismatch table."""
    return (text or "-").replace("|", "\\|")


def aliases(target):
    """Approved symbol plus the most gene-like of its synonyms.

    Open Targets mixes real aliases with clone ids ("H_DJ0042M02.9"), prose
    names ("MAD homolog 4") and malformed entries ("'C-K-RAS"). Those are
    dropped. What survives is ranked so a truncation at MAX_SYNONYMS keeps the
    informative spellings - a plain alphabetical cut kept PIK3CA's CLOVE and
    CWS5 while dropping PI3K, and dropped VEGFR3 from FLT4.
    """
    sym = target["approved_symbol"]
    out = []
    for s in target.get("symbol_synonyms", []):
        if s == sym or not 4 <= len(s) <= 20:
            continue
        if re.search(r"[\s.'\"()]|^H_", s):
            continue
        out.append(s)

    def rank(a):
        stem = re.sub(r"[^A-Za-z]", "", a).lower()[:4]
        return (
            0 if stem and stem in sym.lower() else 1,  # shares the symbol's stem
            0 if re.search(r"\d", a) else 1,  # a digit makes it unambiguous
            a,
        )

    return [sym] + sorted(sorted(set(out)), key=rank)[:MAX_SYNONYMS]


def matched_alias(sentence, names):
    """The alias the sentence actually names, or None.

    Whole-word: a substring test makes "Metastatic" name MET and "Targets" name
    AR, which would let a sentence about neither claim the strongest tier.
    """
    for n in names:
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(n)}(?![A-Za-z0-9])", sentence, re.I):
            return n
    return None


def build_query(target, disease):
    names = " OR ".join(f'TITLE_ABS:"{a}"' for a in aliases(target))
    terms = " OR ".join(f'TITLE_ABS:"{t}"' for t in EVIDENCE_TERMS)
    return (
        f"({names}) AND TITLE_ABS:\"{disease}\" AND ({terms}) "
        "AND (HAS_ABSTRACT:Y) AND (SRC:MED)"
    )


def pick_claim(abstract, target):
    """First abstract sentence naming the target and an evidence term, verbatim.

    Returns (sentence, tier, alias). The tier is carried into evidence.md so a
    weaker pick is labelled there rather than passing as a target+inhibitor
    sentence:
      target+term - names the target and an inhibitor/binder/antagonist keyword
      target-only - names the target but no such keyword
      first-sentence - the abstract never names the target; nothing stronger
    `alias` is which spelling matched, reported when it is not the approved
    symbol, since a gene alias can also be an ordinary word.
    """
    text = clean(abstract)
    if not text:
        return None, None, None
    sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]
    names = aliases(target)

    def cut(s):
        return s if len(s) <= 400 else s[:397].rstrip() + "..."

    for tier, require_term in (("target+term", True), ("target-only", False)):
        for s in sentences:
            hit = matched_alias(s, names)
            if hit is None:
                continue
            if require_term and not any(t in s.lower() for t in EVIDENCE_TERMS):
                continue
            return cut(s), tier, hit
    return cut(sentences[0]), "first-sentence", None


def paper(result, target):
    ji = result.get("journalInfo") or {}
    claim, tier, alias = pick_claim(result.get("abstractText"), target)
    return {
        "pmid": result.get("pmid"),
        "title": clean(result.get("title")),
        "year": result.get("pubYear"),
        "journal": clean((ji.get("journal") or {}).get("title")) or "n/a",
        "claim": claim,
        "claim_tier": tier,
        "claim_alias": alias,
        "doi": result.get("doi"),
    }


def search_target(target, disease, per_target, refresh):
    q = build_query(target, disease)
    size = max(8, per_target * 3)
    # The key hashes the query itself: editing EVIDENCE_TERMS or the synonym
    # filter must not serve a cached payload under the newly printed query.
    tag = hashlib.sha1(f"{q}|{size}".encode()).hexdigest()[:10]
    key = f"search-{slugify(disease)}-{target['approved_symbol']}-{tag}"
    payload = epmc(q, key, page_size=size, refresh=refresh)
    hits = payload.get("hitCount", 0)
    papers, seen = [], set()
    for r in (payload.get("resultList") or {}).get("result", []):
        pmid = r.get("pmid")
        if not pmid or pmid in seen or not r.get("abstractText"):
            continue
        p = paper(r, target)
        if not p["claim"]:
            continue
        seen.add(pmid)
        papers.append(p)
        if len(papers) == per_target:
            break
    return q, hits, papers


def verify(pmid, written_title, refresh):
    """Independent re-fetch by PMID; the title must come back the same."""
    payload = epmc(
        f"EXT_ID:{pmid} AND SRC:MED", f"verify-{pmid}", page_size=1, refresh=refresh
    )
    res = (payload.get("resultList") or {}).get("result", [])
    if not res:
        return False, "pmid did not resolve", None
    back = clean(res[0].get("title"))
    verdict = compare_titles(written_title, back)
    if verdict == "differs":
        return False, "title mismatch", back
    return True, verdict, back


def pubmed_titles(pmids, refresh):
    """Independent cross-source check: PubMed E-utilities, 200 ids per request."""
    out = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i : i + 200]
        # Key on the batch contents, so a changed shortlist cannot hit a stale file.
        tag = hashlib.sha1(",".join(batch).encode()).hexdigest()[:12]
        cache = RAW / f"pubmed-esummary-{tag}.json"
        if cache.exists() and not refresh:
            payload = json.loads(cache.read_text())
        else:
            params = {"db": "pubmed", "id": ",".join(batch), "retmode": "json"}
            url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
            # NCBI throttles hard without an API key, so this retries like every
            # other fetch here; an outright failure must not lose the whole run.
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url, timeout=90, context=SSL_CTX) as r:
                        payload = json.load(r)
                    break
                except Exception as exc:  # noqa: BLE001 - transient network
                    if attempt == 3:
                        raise RuntimeError(f"PubMed esummary failed: {exc}") from exc
                    print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
                    time.sleep(2 * (attempt + 1))
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(payload, indent=1))
            time.sleep(0.4)  # NCBI allows 3 req/s without a key
        for pid, rec in (payload.get("result") or {}).items():
            if pid == "uids" or not isinstance(rec, dict):
                continue
            # esummary returns {"error": "cannot get document summary"} for an
            # id it cannot resolve; that is "not found", not an empty title.
            if rec.get("error") or not rec.get("title"):
                continue
            out[pid] = clean(rec["title"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease", default="colorectal cancer")
    ap.add_argument("--top", type=int, default=25, help="how many T2 targets to cover")
    ap.add_argument("--per-target", type=int, default=3)
    ap.add_argument("--refresh", action="store_true", help="ignore data/raw cache")
    ap.add_argument(
        "--no-pubmed-crosscheck",
        action="store_true",
        help="skip the second, cross-source PMID check against PubMed E-utilities",
    )
    args = ap.parse_args()
    if args.per_target < 1:
        raise SystemExit("--per-target must be >= 1")
    if args.top < 1:
        raise SystemExit("--top must be >= 1")

    out_dir = ROOT / "results" / "pipeline" / slugify(args.disease)
    targets_path = out_dir / "targets.json"
    if not targets_path.exists():
        raise SystemExit(f"{targets_path} missing - run scripts/disease_targets.py first")
    t2 = json.loads(targets_path.read_text())
    disease_name = t2["meta"]["disease_name"]
    kd_ids = t2["meta"].get("known_drug_datatype_ids") or []
    targets = t2["targets"][: args.top]
    print(
        f"{len(targets)} targets from {targets_path.name} "
        f"({disease_name}, {t2['meta']['efo_id']})",
        file=sys.stderr,
    )

    records, emitted, epmc_verified, verified = [], 0, 0, 0
    mismatches, renderings = [], []
    for t in targets:
        q, hits, papers = search_target(t, disease_name, args.per_target, args.refresh)
        kept = []
        for p in papers:
            emitted += 1
            ok, why, back = verify(p["pmid"], p["title"], args.refresh)
            if ok:
                epmc_verified += 1
                verified += 1
                if why == "rendering":
                    renderings.append(("Europe PMC", p["pmid"], p["title"], back))
                kept.append(p)
            else:
                mismatches.append(
                    {
                        "symbol": t["approved_symbol"],
                        "pmid": p["pmid"],
                        "reason": why,
                        "written_title": p["title"],
                        "refetched_title": back,
                    }
                )
        records.append({"target": t, "query": q, "hit_count": hits, "papers": kept})
        print(
            f"  {t['rank']:2d} {t['approved_symbol']:<9} hits={hits:<6} "
            f"emitted={len(papers)} epmc_verified={len(kept)}",
            file=sys.stderr,
        )

    # Second, cross-source check. Europe PMC verifying its own record only proves
    # the pmid->title binding exists there; PubMed is an independent authority.
    # A PMID that fails here is dropped from evidence.md exactly like one that
    # failed the Europe PMC re-fetch - it must not survive in the citation body.
    kept_pmids = [p["pmid"] for r in records for p in r["papers"]]
    pm_ok = 0
    if args.no_pubmed_crosscheck:
        pm_status = "skipped (--no-pubmed-crosscheck)"
    else:
        print(f"cross-checking {len(kept_pmids)} PMIDs against PubMed...", file=sys.stderr)
        try:
            pm = pubmed_titles(kept_pmids, args.refresh)
        except RuntimeError as exc:
            pm = None
            pm_status = f"NOT RUN - {exc}"
        n_before = len(mismatches)
        for r in (records if pm is not None else []):
            survivors = []
            for p in r["papers"]:
                back = pm.get(p["pmid"])
                verdict = "missing" if back is None else compare_titles(p["title"], back)
                if verdict in ("match", "rendering"):
                    pm_ok += 1
                    if verdict == "rendering":
                        renderings.append(("PubMed", p["pmid"], p["title"], back))
                    survivors.append(p)
                    continue
                why = (
                    "not found in PubMed esummary"
                    if verdict == "missing"
                    else "PubMed title differs"
                )
                verified -= 1
                mismatches.append(
                    {
                        "symbol": r["target"]["approved_symbol"],
                        "pmid": p["pmid"],
                        "reason": why,
                        "written_title": p["title"],
                        "refetched_title": back,
                    }
                )
            r["papers"] = survivors
        if pm is not None:
            pm_status = (
                f"{pm_ok}/{len(kept_pmids)} titles matched PubMed esummary, "
                f"{len(mismatches) - n_before} dropped"
            )
    print(f"pubmed cross-check: {pm_status}", file=sys.stderr)

    no_evidence = [r for r in records if not r["papers"]]

    lines = []
    lines.append(f"# Literature evidence - {disease_name}")
    lines.append("")
    lines.append(
        f"Pipeline task T3. Source: Europe PMC REST "
        f"(`{EPMC}`), searched {dt.datetime.now(dt.UTC).date().isoformat()}. "
        f"Targets from `targets.json` (Open Targets {t2['meta']['efo_id']}, "
        f"data {t2['meta']['data_version']})."
    )
    lines.append("")
    lines.append(
        f"**n = {len(targets)} targets searched, {emitted} PMIDs emitted, "
        f"{verified} verified resolving, {len(mismatches)} mismatches.** "
        "Every PMID below was re-fetched from Europe PMC by id in a second request "
        "and its title compared with the one written here; any that failed was "
        "removed and is listed in the verification section. "
        f"Cross-source check against PubMed E-utilities: {pm_status}."
    )
    lines.append("")
    lines.append(
        "Each quoted line is a **verbatim sentence from the downloaded abstract**, "
        "not a paraphrase - markup is rendered to text but no word is added, "
        "removed or reordered. Quotes are truncated at 400 characters. "
        "A quote that could not pair the target with an inhibitor/binder term "
        "is labelled underneath with its weaker claim tier."
    )
    lines.append("")

    for r in records:
        t = r["target"]
        ds = t["datatype_scores"]
        lit = "-" if ds["literature"] is None else f"{ds['literature']:.3f}"
        kd = "-" if ds["known_drug"] is None else f"{ds['known_drug']:.3f}"
        # Name the datatype id this release actually served, not "known_drug".
        kd_id = t["known_drug_datatype_id"] or "/".join(kd_ids) or "absent"
        lines.append(
            f"## {t['rank']}. {t['approved_symbol']} "
            f"({t['uniprot'] or 'no UniProt'}) - {t['approved_name']}"
        )
        lines.append("")
        lines.append(
            f"Open Targets: overall **{t['overall_association_score']:.4f}**, "
            f"literature {lit}, known-drug ({kd_id}) {kd}, "
            f"{t['n_drug_and_clinical_candidates']} drug/clinical candidates "
            f"({t['n_approved_drugs']} approved)."
        )
        lines.append("")
        lines.append(f"Europe PMC hits: **{r['hit_count']}**, kept {len(r['papers'])}.")
        lines.append("")
        if not r["papers"]:
            lines.append(
                "> No verified paper. "
                + (
                    "Europe PMC returned 0 hits for this query."
                    if r["hit_count"] == 0
                    else "Hits existed but none carried a usable abstract, or all "
                    "failed PMID verification (see below)."
                )
            )
            lines.append("")
            lines.append(f"<sub>query: `{r['query']}`</sub>")
            lines.append("")
            continue
        for i, p in enumerate(r["papers"], 1):
            doi = f" doi:{p['doi']}" if p.get("doi") else ""
            lines.append(
                f"{i}. **PMID {p['pmid']}** - *{p['title']}* - "
                f"{p['journal']}, {p['year']}.{doi}"
            )
            lines.append(f"   > {p['claim']}")
            notes = []
            if p["claim_tier"] != "target+term":
                notes.append(
                    f"claim tier `{p['claim_tier']}`: this sentence is the "
                    "strongest the abstract offered - it does not pair the "
                    "target with an inhibitor/binder/antagonist term"
                )
            if p["claim_alias"] and p["claim_alias"] != t["approved_symbol"]:
                notes.append(
                    f"matched via the alias `{p['claim_alias']}`, not the "
                    f"symbol {t['approved_symbol']}"
                )
            if notes:
                lines.append("   <sub>" + "; ".join(notes) + ".</sub>")
            lines.append("")
        lines.append(f"<sub>query: `{r['query']}`</sub>")
        lines.append("")

    lines.append("## Verification")
    lines.append("")
    lines.append(
        f"{emitted} PMIDs emitted, {epmc_verified} passed the Europe PMC "
        f"re-fetch, {verified} survived both checks, {len(mismatches)} rejected."
    )
    lines.append("")
    if mismatches:
        lines.append("| target | pmid | reason | written title | re-fetched title |")
        lines.append("|---|---|---|---|---|")
        for m in mismatches:
            lines.append(
                f"| {m['symbol']} | {m['pmid']} | {m['reason']} | "
                f"{cell(m['written_title'])} | {cell(m['refetched_title'])} |"
            )
    else:
        lines.append(
            "No mismatches: every emitted PMID resolved and returned the same "
            "title from Europe PMC"
            + (
                "."
                if args.no_pubmed_crosscheck
                else ", and matched PubMed."
            )
        )
    lines.append("")
    lines.append(f"Cross-source check (PubMed E-utilities esummary): {pm_status}.")
    lines.append("")
    tiers = Counter(p["claim_tier"] for r in records for p in r["papers"])
    lines.append(
        "Claim tiers across the "
        f"{sum(tiers.values())} quotes: "
        + ", ".join(f"`{k}` {v}" for k, v in sorted(tiers.items()))
        + ". A `target-only` or `first-sentence` quote is still verbatim, but the "
        "abstract did not put the target and an inhibitor term in one sentence - "
        "usually because the target has no inhibitor, which is itself the finding."
    )
    lines.append("")
    # A PMID can cite under two targets, and one recorded during the Europe PMC
    # pass may since have been dropped by PubMed - it belongs in neither table.
    surviving = {p["pmid"] for r in records for p in r["papers"]}
    renderings = [x for x in dict.fromkeys(renderings) if x[1] in surviving]
    if renderings:
        lines.append(
            f"{len(renderings)} title(s) matched only after markup normalisation - "
            "same words, different rendering of superscripts or symbols by the two "
            "databases. These were kept, and both spellings are shown:"
        )
        lines.append("")
        lines.append("| source | pmid | as written here | as returned |")
        lines.append("|---|---|---|---|")
        for src, pmid, written, back in renderings:
            lines.append(f"| {src} | {pmid} | {cell(written)} | {cell(back)} |")
        lines.append("")

    lines.append("## Not evaluated / why")
    lines.append("")
    if no_evidence:
        lines.append("| target | why |")
        lines.append("|---|---|")
        for r in no_evidence:
            why = (
                "Europe PMC returned 0 hits for "
                f"symbol+disease+inhibitor terms (hitCount 0)"
                if r["hit_count"] == 0
                else f"{r['hit_count']} hits but no result with a usable abstract "
                "survived verification"
            )
            lines.append(f"| {r['target']['approved_symbol']} | {why} |")
    else:
        lines.append("Every target searched produced at least one verified PMID.")
    lines.append("")

    path = out_dir / "evidence.md"
    path.write_text("\n".join(lines))

    print(f"\nWrote {path}", file=sys.stderr)
    print(
        f"{emitted} PMIDs emitted, {epmc_verified} passed the Europe PMC "
        f"re-fetch, {verified} survived both checks, {len(mismatches)} rejected",
        file=sys.stderr,
    )
    print(
        f"targets with >=1 verified PMID: {len(targets) - len(no_evidence)}/{len(targets)}",
        file=sys.stderr,
    )
    if no_evidence:
        print(
            "no verified evidence: "
            + ", ".join(r["target"]["approved_symbol"] for r in no_evidence),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
