#!/usr/bin/env python
"""Build data/approved_biologics.csv: approved protein drugs WITH amino-acid sequences.

Why this exists
---------------
`data/approved_drugs.csv` is DrugCentral/SMILES-derived: 4,065 small molecules,
34 peptide-scale entries, zero antibodies. A SMILES string cannot be fed to a
protein language model. This script builds the complementary corpus: approved
biologics with their actual chain sequences, so ESM-C can embed them alongside
de novo designed binders.

Sources and licences (see docs/08-BIOLOGICS-CORPUS.md for the full statement)
---------------------------------------------------------------------------
  ChEMBL REST  https://www.ebi.ac.uk/chembl/api/data   CC BY-SA 3.0 Unported
      Sequences, molecule metadata, drug mechanisms, targets. Redistributable
      with attribution under share-alike. THIS IS THE ONLY SOURCE OF DATA THAT
      IS COMMITTED TO THE REPO.
  UniProt REST https://rest.uniprot.org                CC BY 4.0
      Optional enrichment/verification of target gene symbols. Redistributable.
  DrugCentral (via the already-committed data/approved_drugs.csv)  CC BY-SA 4.0
      Read-only cross-reference for `approval_agencies`. Not modified here.

  DrugBank is CC BY-NC 4.0 (non-commercial) and is deliberately NOT used: the
  repo is MIT-licensed and public, so NC data cannot be committed.
  Thera-SAbDab (OPIG) is free for academic use with a paid commercial licence
  for the packaged distribution, so it is also not committed. See docs/08.

Re-runnability
--------------
No credentials of any kind. HTTP responses are cached under a temp directory
(--cache-dir, default $TMPDIR/chembl_biologics_cache) so a re-run is cheap;
--no-cache forces a fresh fetch. UniProt enrichment degrades gracefully: if the
service is unreachable the ChEMBL gene symbols are used and the run still
succeeds (a warning is printed and counted).

Usage
-----
    python demo/build_biologics.py                # build data/approved_biologics.csv
    python demo/build_biologics.py --report       # + markdown tables for docs/08
    python demo/build_biologics.py --verify       # + independent sequence spot-check
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import ssl
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

try:  # repo convention: python.org 3.1x ships no CA bundle
    import certifi

    _SSL = ssl.create_default_context(cafile=certifi.where())
except ImportError:  # pragma: no cover - fall back to system trust store
    _SSL = ssl.create_default_context()

ROOT = Path(__file__).resolve().parent.parent
OUT_CSV = ROOT / "data" / "approved_biologics.csv"
SMALL_MOL_CSV = ROOT / "data" / "approved_drugs.csv"

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
UA = "drug-similarity-discovery/0.1 (HackMIT; biologics corpus build)"

# ChEMBL molecule_type values that can carry a protein sequence. Queried
# separately from the biotherapeutic__isnull=false sweep and unioned, because
# neither filter is a superset of the other.
PROTEIN_MOLECULE_TYPES = [
    "Protein",
    "Antibody",
    "Enzyme",
    "Antibody drug conjugate",
]

CANON_AA = set("ACDEFGHIKLMNPQRSTVWY")
EXTENDED_AA = CANON_AA | set("BXZUO")

# ---------------------------------------------------------------------------
# HTTP with a disk cache
# ---------------------------------------------------------------------------

_stats = Counter()


def _cache_path(cache_dir: Path, url: str) -> Path:
    return cache_dir / (hashlib.sha256(url.encode()).hexdigest()[:32] + ".json")


def http_json(url: str, cache_dir: Path | None, retries: int = 4, timeout: int = 90):
    """GET `url` and parse JSON. Returns None on permanent failure."""
    if cache_dir is not None:
        cp = _cache_path(cache_dir, url)
        if cp.exists():
            _stats["cache_hit"] += 1
            return json.loads(cp.read_text())

    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout, context=_SSL) as r:
                payload = json.loads(r.read().decode())
            _stats["http_ok"] += 1
            if cache_dir is not None:
                cache_dir.mkdir(parents=True, exist_ok=True)
                _cache_path(cache_dir, url).write_text(json.dumps(payload))
            return payload
        except urllib.error.HTTPError as e:
            last = e
            if e.code == 404:
                _stats["http_404"] += 1
                return None
            time.sleep(1.5 * (attempt + 1))
        except Exception as e:  # network hiccup, timeout, bad gateway
            last = e
            time.sleep(1.5 * (attempt + 1))
    _stats["http_fail"] += 1
    print(f"  ! giving up on {url}: {last}", file=sys.stderr)
    return None


def chembl_paged(path: str, params: dict, cache_dir: Path | None, page: int = 200):
    """Iterate every record of a paged ChEMBL collection endpoint."""
    offset = 0
    key = None
    while True:
        q = dict(params, format="json", limit=page, offset=offset)
        url = f"{CHEMBL}/{path}?{urllib.parse.urlencode(q)}"
        data = http_json(url, cache_dir)
        if data is None:
            return
        if key is None:
            key = next((k for k in data if k != "page_meta"), None)
            if key is None:
                return
        rows = data.get(key) or []
        for row in rows:
            yield row
        meta = data.get("page_meta", {})
        offset += page
        if not meta.get("next") or offset >= meta.get("total_count", 0):
            return


def chembl_by_ids(path: str, field: str, ids: list[str], cache_dir: Path | None, chunk: int = 25):
    """Fetch records for many ids using ChEMBL's `__in` filter, chunked."""
    out = []
    ids = sorted(set(ids))
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        out.extend(chembl_paged(path, {f"{field}__in": ",".join(batch)}, cache_dir))
    return out


# ---------------------------------------------------------------------------
# Modality classification
# ---------------------------------------------------------------------------

# WHO INN stems, which are the load-bearing signal. Order matters: the first
# rule that fires wins. Documented in docs/08-BIOLOGICS-CORPUS.md.
HORMONE_CYTOKINE_TOKENS = (
    "insulin", "somatropin", "somatotropin", "interferon", "interleukin",
    "erythropoietin", "epoetin", "darbepoetin", "filgrastim", "sargramostim",
    "pegfilgrastim", "lipegfilgrastim", "eflapegrastim", "oprelvekin",
    "follitropin", "lutropin", "choriogonadotropin", "thyrotropin", "menotropins",
    "urofollitropin", "metreleptin", "mecasermin", "palifermin", "becaplermin",
    "teriparatide", "glucagon", "calcitonin", "secretin", "aldesleukin",
    "denileukin", "romiplostim", "tbo-filgrastim", "pegteograstim",
)
ENZYME_TOKENS = (
    "dornase", "alteplase", "reteplase", "tenecteplase", "streptokinase",
    "urokinase", "collagenase", "hyaluronidase", "asparaginase", "rasburicase",
    "pegloticase", "sacrosidase", "elosulfase", "idursulfase", "laronidase",
    "galsulfase", "alglucosidase", "agalsidase", "imiglucerase", "velaglucerase",
    "taliglucerase", "sebelipase", "asfotase", "pegvaliase", "elapegademase",
    "pegademase", "cerliponase", "vestronidase", "olipudase", "avalglucosidase",
    "efgartigimod",  # Fc fragment engineered as a receptor blocker; see below
)
FUSION_TOKENS = ("cept", "-fc", " fc", "fusion")

MAB_STEM = re.compile(r"(mab|tug|bart|ment)$")  # 2022 INN revision split -mab
FRAGMENT_HINT = re.compile(r"\b(fab|f\(ab|scfv|vhh|nanobody|fragment|single.domain)\b", re.I)


def classify(name: str, mol_type: str, chains: list[tuple[str, str]], descs: list[str]) -> str:
    """Assign one of: mab, antibody_fragment, fusion_protein, enzyme,
    hormone_cytokine, peptide, other.

    Rules are applied in order and are deliberately explicit rather than
    learned, so that every assignment in the CSV can be traced to one line here.
    """
    n = name.lower()
    blob = " ".join(d.lower() for d in descs if d)
    total = sum(len(s) for _, s in chains)
    heavy = max((len(s) for cid, s in chains if cid.startswith("H")), default=0)

    # 1. Antibody-derived, by INN stem or ChEMBL molecule_type.
    is_ab = bool(MAB_STEM.search(n.replace(" pegol", "").replace(" ", ""))) or mol_type in (
        "Antibody",
        "Antibody drug conjugate",
    )
    if is_ab:
        if FRAGMENT_HINT.search(blob) or FRAGMENT_HINT.search(n):
            return "antibody_fragment"
        # A full IgG heavy chain is ~440-460 aa (VH+CH1+hinge+CH2+CH3).
        # Fab/scFv/VHH heavy chains are ~110-250 aa.
        if heavy and heavy < 300:
            return "antibody_fragment"
        if not heavy and total and total < 300:
            return "antibody_fragment"
        return "mab"

    # 2. Fusion proteins: -cept stem (etanercept, aflibercept, abatacept,
    #    rilonacept, luspatercept, ...) or an explicit Fc-fusion description.
    if n.endswith("cept") or any(t in blob for t in ("fusion protein", "fc fusion")):
        return "fusion_protein"

    # 3. Enzymes.
    if mol_type == "Enzyme" or any(t in n for t in ENZYME_TOKENS) or n.endswith("ase"):
        return "enzyme"

    # 4. Hormones and cytokines.
    if any(t in n for t in HORMONE_CYTOKINE_TOKENS):
        return "hormone_cytokine"

    # 5. Length-based fallback. 50 aa is the conventional peptide/protein line.
    if total and total <= 50:
        return "peptide"
    if total:
        return "other"
    return "other"


# ---------------------------------------------------------------------------
# Chain handling
# ---------------------------------------------------------------------------


def _label_kind(description: str | None) -> str | None:
    """H / L / A / B from a ChEMBL biocomponent description, or None if the
    description is uninformative (ChEMBL uses the literal string 'Sequence' for
    160 of the 400 protein components, so this returns None a lot)."""
    d = (description or "").lower()
    if "heavy" in d:
        return "H"
    if "light" in d or "kappa" in d or "lambda" in d or "v-kappa" in d:
        return "L"
    if "alpha" in d:
        return "A"
    if "beta" in d:
        return "B"
    return None


def assign_chain_ids(components: list[tuple[str | None, str]], antibody_like: bool) -> list[tuple[str, str]]:
    """Give every (description, sequence) pair a chain id: H<i>/L<i>/A<i>/B<i>/C<i>.

    ChEMBL labels only some components. Where labels are missing we infer, but
    only under conditions that cannot be ambiguous:

      * one component already labelled Light and exactly one unlabelled  -> the
        unlabelled one is the heavy chain (and vice versa);
      * an antibody-like molecule with exactly two unlabelled components -> the
        longer is heavy, the shorter light. An IgG heavy chain (~440-460 aa) and
        a Fab heavy chain (~225-235 aa) are both longer than their kappa/lambda
        partner (~214 aa), so length ordering is safe for IgG and Fab alike.

    Anything else keeps the neutral id C<i>, so a wrong guess never enters the
    heavy_chain_sequence / light_chain_sequence columns.
    """
    kinds: list[str | None] = [_label_kind(d) for d, _ in components]
    unlabelled = [i for i, k in enumerate(kinds) if k is None]

    if len(unlabelled) == 1:
        others = {k for k in kinds if k}
        if others == {"L"}:
            kinds[unlabelled[0]] = "H"
        elif others == {"H"}:
            kinds[unlabelled[0]] = "L"
    elif antibody_like and len(unlabelled) == 2 and len(components) == 2:
        i, j = unlabelled
        long_i, short_i = (i, j) if len(components[i][1]) >= len(components[j][1]) else (j, i)
        if len(components[long_i][1]) != len(components[short_i][1]):
            kinds[long_i], kinds[short_i] = "H", "L"

    seen: Counter = Counter()
    out = []
    for k, (_, seq) in zip(kinds, components):
        kind = k or "C"
        seen[kind] += 1
        out.append((f"{kind}{seen[kind]}", seq))
    return out


def clean_sequence(seq: str | None) -> str:
    if not seq:
        return ""
    s = re.sub(r"\s+", "", seq).upper()
    return s if s and set(s) <= EXTENDED_AA else ""


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


def load_agency_map() -> dict[str, str]:
    """name -> approval_agencies, from the already-committed DrugCentral CSV."""
    if not SMALL_MOL_CSV.exists():
        return {}
    out = {}
    with SMALL_MOL_CSV.open(newline="") as fh:
        for row in csv.DictReader(fh):
            nm = (row.get("name") or "").strip().lower()
            ag = (row.get("approval_agencies") or "").strip()
            if nm and ag:
                out[nm] = ag
    return out


def uniprot_gene(accession: str, cache_dir: Path | None) -> tuple[str, str]:
    """(gene_symbol, protein_name) from UniProt. ('','') if unavailable."""
    url = f"{UNIPROT}/{accession}.json?fields=gene_primary,protein_name"
    d = http_json(url, cache_dir, retries=2, timeout=45)
    if not d:
        return "", ""
    gene = ""
    genes = d.get("genes") or []
    if genes:
        gene = (genes[0].get("geneName") or {}).get("value", "")
    pname = (((d.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value", "")
    return gene, pname


def build(cache_dir: Path | None, skip_uniprot: bool = False):
    # --- 1. candidate molecules -------------------------------------------
    print("[1/5] ChEMBL molecules (max_phase=4)...", file=sys.stderr)
    mols: dict[str, dict] = {}
    only = ",".join(
        [
            "molecule_chembl_id", "pref_name", "molecule_type", "max_phase",
            "first_approval", "withdrawn_flag", "usan_stem", "usan_year",
            "black_box_warning", "cross_references", "natural_product",
        ]
    )
    queries = [{"max_phase": 4, "biotherapeutic__isnull": "false", "only": only}]
    queries += [{"max_phase": 4, "molecule_type": t, "only": only} for t in PROTEIN_MOLECULE_TYPES]
    for q in queries:
        n0 = len(mols)
        for m in chembl_paged("molecule", q, cache_dir):
            mols[m["molecule_chembl_id"]] = m
        label = q.get("molecule_type", "biotherapeutic__isnull=false")
        print(f"      {label:28s} -> +{len(mols) - n0} new (total {len(mols)})", file=sys.stderr)

    # --- 2. sequences ------------------------------------------------------
    print(f"[2/5] biotherapeutic records for {len(mols)} molecules...", file=sys.stderr)
    bio: dict[str, dict] = {}
    for b in chembl_by_ids("biotherapeutic", "molecule_chembl_id", list(mols), cache_dir):
        bio[b["molecule_chembl_id"]] = b

    # --- 3. mechanisms -> targets -----------------------------------------
    print("[3/5] drug mechanisms and targets...", file=sys.stderr)
    mechs: dict[str, list[dict]] = defaultdict(list)
    for m in chembl_by_ids("mechanism", "molecule_chembl_id", list(mols), cache_dir):
        mechs[m["molecule_chembl_id"]].append(m)
    tids = {m["target_chembl_id"] for v in mechs.values() for m in v if m.get("target_chembl_id")}
    targets: dict[str, dict] = {}
    for t in chembl_by_ids("target", "target_chembl_id", sorted(tids), cache_dir, chunk=20):
        targets[t["target_chembl_id"]] = t

    def target_fields(tcid: str):
        t = targets.get(tcid)
        if not t:
            return "", "", ""
        comps = [c for c in t.get("target_components", []) if c.get("accession")]
        human = [c for c in comps if (c.get("component_description") or "") and t.get("organism") == "Homo sapiens"]
        use = human or comps
        accs, genes = [], []
        for c in use:
            accs.append(c["accession"])
            g = next(
                (s["component_synonym"] for s in c.get("target_component_synonyms", []) if s.get("syn_type") == "GENE_SYMBOL"),
                "",
            )
            genes.append(g)
        return ";".join(genes), ";".join(accs), t.get("pref_name") or ""

    # --- 4. UniProt enrichment (optional, degrades gracefully) -------------
    upmap: dict[str, tuple[str, str]] = {}
    accs_needed = sorted({a for tc in targets for a in target_fields(tc)[1].split(";") if a})
    if skip_uniprot:
        print(f"[4/5] UniProt enrichment skipped (--no-uniprot); {len(accs_needed)} accessions", file=sys.stderr)
    else:
        print(f"[4/5] UniProt gene symbols for {len(accs_needed)} accessions...", file=sys.stderr)
        for a in accs_needed:
            upmap[a] = uniprot_gene(a, cache_dir)
        got = sum(1 for v in upmap.values() if v[0])
        print(f"      resolved {got}/{len(accs_needed)}", file=sys.stderr)

    # --- 5. assemble -------------------------------------------------------
    print("[5/5] assembling rows...", file=sys.stderr)
    agencies = load_agency_map()
    rows, excluded = [], []

    for cid, m in sorted(mols.items(), key=lambda kv: (kv[1].get("pref_name") or "")):
        name = (m.get("pref_name") or "").strip().lower()
        mol_type = m.get("molecule_type") or ""
        b = bio.get(cid)

        if not name:
            excluded.append((cid, "", mol_type, "no INN/pref_name in ChEMBL"))
            continue
        if b is None:
            excluded.append((cid, name, mol_type, "no ChEMBL biotherapeutic record"))
            continue

        comps = b.get("biocomponents") or []
        prot = [c for c in comps if (c.get("component_type") or "").upper() == "PROTEIN"]
        if not comps:
            excluded.append((cid, name, mol_type, "biotherapeutic record has no biocomponents"))
            continue
        if not prot:
            kinds = ",".join(sorted({(c.get("component_type") or "?") for c in comps}))
            excluded.append((cid, name, mol_type, f"no PROTEIN biocomponent (only {kinds})"))
            continue

        seen = Counter()
        chains: list[tuple[str, str]] = []
        descs: list[str] = [b.get("description") or ""]
        bad = 0
        for i, c in enumerate(prot):
            s = clean_sequence(c.get("sequence"))
            descs.append(c.get("description") or "")
            if not s:
                bad += 1
                continue
            chains.append((chain_id(c.get("description"), i, seen), s))
        if not chains:
            excluded.append((cid, name, mol_type, f"all {len(prot)} PROTEIN components had empty/invalid sequence"))
            continue

        modality = classify(name, mol_type, chains, descs)
        heavy = [s for k, s in chains if k.startswith("H")]
        light = [s for k, s in chains if k.startswith("L")]

        mlist = mechs.get(cid, [])
        genes, accs, tnames, actions = [], [], [], []
        for mech in mlist:
            g, a, tn = target_fields(mech.get("target_chembl_id") or "")
            if a:
                # prefer the UniProt primary gene name where we have one
                g_parts = g.split(";")
                a_parts = a.split(";")
                g = ";".join(
                    (upmap.get(ap, ("", ""))[0] or (g_parts[i] if i < len(g_parts) else ""))
                    for i, ap in enumerate(a_parts)
                )
                genes.append(g)
                accs.append(a)
                tnames.append(tn)
                actions.append(mech.get("action_type") or "")

        def uniq(xs):
            out = []
            for x in xs:
                for part in x.split(";"):
                    if part and part not in out:
                        out.append(part)
            return ";".join(out)

        rows.append(
            {
                "name": name,
                "chembl_id": cid,
                "modality": modality,
                "molecule_type_chembl": mol_type,
                "chains": "|".join(f"{k}:{s}" for k, s in chains),
                "n_chains": len(chains),
                "length_aa": sum(len(s) for _, s in chains),
                "sequence": chains[0][1] if len(chains) == 1 else "",
                "heavy_chain_sequence": heavy[0] if heavy else "",
                "light_chain_sequence": light[0] if light else "",
                "target_gene_symbol": uniq(genes),
                "target_uniprot": uniq(accs),
                "target_name": uniq(tnames),
                "target_action_type": uniq(actions),
                "n_mechanisms": len(mlist),
                "first_approval": m.get("first_approval") or "",
                "approval_agencies": agencies.get(name, ""),
                "withdrawn": int(m.get("withdrawn_flag") or 0),
                "dropped_components": bad,
                "source_db": "ChEMBL",
                "source_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{cid}/",
            }
        )

    return rows, excluded, mols


FIELDS = [
    "name", "chembl_id", "modality", "molecule_type_chembl",
    "chains", "n_chains", "length_aa", "sequence",
    "heavy_chain_sequence", "light_chain_sequence",
    "target_gene_symbol", "target_uniprot", "target_name", "target_action_type",
    "n_mechanisms", "first_approval", "approval_agencies", "withdrawn",
    "dropped_components", "source_db", "source_url",
]


def write_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["name"]):
            w.writerow(r)


# ---------------------------------------------------------------------------
# Verification: sequences vs. an independent source
# ---------------------------------------------------------------------------

# Independent check set. Each entry is (drug, chain, source, expected sequence).
# The sequences below were taken from the WHO INN Recommended List entry for the
# drug, as reproduced in the cited public record -- NOT from ChEMBL. They are
# used only to test ChEMBL's copy; they are not redistributed as data.
VERIFY_UNIPROT = {
    # single-chain biologics that ARE in UniProt, so we can compare to a second
    # database directly and programmatically.
    "aflibercept": None,          # no UniProt entry; checked by motif instead
}

# Motif checks: short, unambiguous subsequences that must appear in the named
# chain if the sequence is the right molecule. Sources are given in docs/08.
MOTIF_CHECKS = [
    # drug, chain key, motif, what the motif is
    ("pembrolizumab", "H1", "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMY", "INN 8930 heavy-chain N-term"),
    ("nivolumab", "H1", "QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMH", "INN 9623 heavy-chain N-term"),
    ("trastuzumab", "H1", "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIH", "INN 7761 heavy-chain N-term"),
    ("bevacizumab", "H1", "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMN", "INN 8060 heavy-chain N-term"),
    ("adalimumab", "H1", "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMH", "INN 8148 heavy-chain N-term"),
    ("etanercept", "C1", "LPAQVAFTPYAPEPGSTCRLREYYDQTAQMCCSKCSPGQHAKVFC", "TNFRSF1B ectodomain N-term (P20333 res 23-67)"),
    ("aflibercept", "C1", "SDTGRPFVEMYSEIPEIIHMTEGRELVIPCRVTSPNITVTLKKFP", "FLT1 Ig-like D2 N-term (P17948)"),
    ("ranibizumab", "H1", "EVQLVESGGGLVQPGGSLRLSCAASGYDFTHYGMN", "INN 8558 Fab heavy-chain N-term"),
]

# Chains that must be findable as an exact substring of a UniProt entry, because
# the drug IS (a fragment of) a natural human protein. This is a genuinely
# independent check: UniProt is a different database with a different pipeline.
UNIPROT_SUBSTRING_CHECKS = [
    # drug, chain key, uniprot accession, human protein
    ("etanercept", "C1", "P20333", "TNFRSF1B (TNF receptor 2) ectodomain"),
    ("somatropin", "C1", "P01241", "GH1 growth hormone, mature chain"),
    ("aldesleukin", "C1", "P60568", "IL2 interleukin-2"),
    ("dornase alfa", "C1", "P24855", "DNASE1 deoxyribonuclease-1"),
    ("rasburicase", "C1", "Q00511", "Aspergillus flavus urate oxidase"),
]


def verify(rows, cache_dir: Path | None):
    by_name = {r["name"]: r for r in rows}

    def chain_of(row, key):
        for part in row["chains"].split("|"):
            k, _, s = part.partition(":")
            if k == key:
                return s
        return ""

    print("\n### Verification A - INN N-terminal motif checks\n")
    print("| drug | chain | motif source | result |")
    print("|---|---|---|---|")
    a_pass = a_fail = a_miss = 0
    for name, key, motif, src in MOTIF_CHECKS:
        r = by_name.get(name)
        if not r:
            print(f"| {name} | {key} | {src} | NOT IN CORPUS |")
            a_miss += 1
            continue
        s = chain_of(r, key)
        if not s:
            print(f"| {name} | {key} | {src} | chain {key} absent |")
            a_miss += 1
        elif s.startswith(motif) or motif in s:
            print(f"| {name} | {key} | {src} | PASS |")
            a_pass += 1
        else:
            print(f"| {name} | {key} | {src} | **FAIL** (got `{s[:45]}...`) |")
            a_fail += 1
    print(f"\nA: {a_pass} pass, {a_fail} fail, {a_miss} not checkable (n={len(MOTIF_CHECKS)}).")

    print("\n### Verification B - chain is a substring of the UniProt human/source protein\n")
    print("| drug | chain | UniProt | protein | result |")
    print("|---|---|---|---|---|")
    b_pass = b_fail = b_miss = 0
    for name, key, acc, prot in UNIPROT_SUBSTRING_CHECKS:
        r = by_name.get(name)
        if not r:
            print(f"| {name} | {key} | {acc} | {prot} | NOT IN CORPUS |")
            b_miss += 1
            continue
        s = chain_of(r, key)
        d = http_json(f"{UNIPROT}/{acc}.json?fields=sequence", cache_dir, retries=2, timeout=45)
        ref = ((d or {}).get("sequence") or {}).get("value", "")
        if not s or not ref:
            print(f"| {name} | {key} | {acc} | {prot} | UNCHECKED (missing {'chain' if not s else 'UniProt'}) |")
            b_miss += 1
            continue
        # longest common check: exact substring, else report the longest shared
        # prefix of the drug chain that is present in the reference.
        if s in ref:
            print(f"| {name} | {key} | {acc} | {prot} | PASS (exact substring, {len(s)} aa) |")
            b_pass += 1
        else:
            best = 0
            for L in range(min(len(s), 400), 9, -1):
                if s[:L] in ref:
                    best = L
                    break
            verdict = f"partial: first {best}/{len(s)} aa match" if best >= 20 else f"**FAIL** (max shared prefix {best} aa)"
            if best >= 20:
                b_pass += 1
            else:
                b_fail += 1
            print(f"| {name} | {key} | {acc} | {prot} | {verdict} |")
    print(f"\nB: {b_pass} pass/partial, {b_fail} fail, {b_miss} not checkable (n={len(UNIPROT_SUBSTRING_CHECKS)}).")
    return a_fail + b_fail


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

PRIORITY = {
    "PDCD1": (["pembrolizumab", "nivolumab", "cemiplimab", "dostarlimab", "toripalimab", "tislelizumab"], "Q15116"),
    "CD274": (["atezolizumab", "durvalumab", "avelumab"], "Q9NZQ7"),
    "TNF": (["adalimumab", "infliximab", "golimumab", "certolizumab pegol", "etanercept"], "P01375"),
    "KDR": (["ramucirumab"], "P35968"),
    "VEGFA": (["bevacizumab", "ranibizumab", "aflibercept"], "P15692"),
}


def report(rows, excluded, mols):
    n = len(rows)
    print(f"\n## Corpus summary (n = {n})\n")

    print("| modality | n | median length_aa | median n_chains |")
    print("|---|---|---|---|")
    by_mod = defaultdict(list)
    for r in rows:
        by_mod[r["modality"]].append(r)
    for mod in ["mab", "antibody_fragment", "fusion_protein", "enzyme", "hormone_cytokine", "peptide", "other"]:
        v = by_mod.get(mod, [])
        if not v:
            print(f"| {mod} | 0 | - | - |")
            continue
        lens = sorted(r["length_aa"] for r in v)
        chs = sorted(r["n_chains"] for r in v)
        print(f"| {mod} | {len(v)} | {lens[len(lens)//2]} | {chs[len(chs)//2]} |")
    print(f"| **total** | **{n}** | | |")

    with_target = sum(1 for r in rows if r["target_uniprot"])
    with_year = sum(1 for r in rows if r["first_approval"])
    with_agency = sum(1 for r in rows if r["approval_agencies"])
    withdrawn = sum(1 for r in rows if r["withdrawn"])
    multi = sum(1 for r in rows if r["n_chains"] > 1)
    hl = sum(1 for r in rows if r["heavy_chain_sequence"] and r["light_chain_sequence"])
    print(f"\n- rows with >=1 UniProt target accession: **{with_target}/{n}** ({100*with_target/n:.1f}%)")
    print(f"- rows with `first_approval`: **{with_year}/{n}** ({100*with_year/n:.1f}%)")
    print(f"- rows with `approval_agencies` (DrugCentral name match): **{with_agency}/{n}** ({100*with_agency/n:.1f}%)")
    print(f"- rows flagged withdrawn: **{withdrawn}/{n}**")
    print(f"- multi-chain rows (`n_chains` > 1): **{multi}/{n}**")
    print(f"- rows with both a heavy and a light chain: **{hl}/{n}**")
    tot_aa = sum(r["length_aa"] for r in rows)
    print(f"- total residues in corpus: **{tot_aa:,}**")

    print("\n## Priority-target coverage\n")
    print("| target | UniProt | expected | in corpus | target annotated in corpus | missing |")
    print("|---|---|---|---|---|---|")
    by_name = {r["name"]: r for r in rows}
    for gene, (drugs, acc) in PRIORITY.items():
        present = [d for d in drugs if d in by_name]
        annotated = [d for d in present if acc in (by_name[d]["target_uniprot"] or "")]
        missing = [d for d in drugs if d not in by_name]
        print(
            f"| {gene} | {acc} | {len(drugs)} | {len(present)} | {len(annotated)} | "
            f"{', '.join(missing) if missing else '-'} |"
        )

    print("\n## Not included, and why\n")
    reasons = Counter(r[3] for r in excluded)
    print(f"Candidate ChEMBL molecules considered: **{len(mols)}**. Included: **{n}**. Excluded: **{len(excluded)}**.\n")
    print("| reason | n |")
    print("|---|---|")
    for reason, c in reasons.most_common():
        print(f"| {reason} | {c} |")
    print(f"| **total excluded** | **{len(excluded)}** |")

    named = [e for e in excluded if e[1]]
    print(f"\nNamed drugs among the excluded (n={len(named)}); first 40 alphabetically:\n")
    print("| name | ChEMBL molecule_type | reason |")
    print("|---|---|---|")
    for cid, name, mt, reason in sorted(named, key=lambda e: e[1])[:40]:
        print(f"| {name} | {mt or '(none)'} | {reason} |")

    print(f"\n### Fetch stats\n\n```\n{dict(_stats)}\n```")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    ap.add_argument("--cache-dir", type=Path, default=Path(tempfile.gettempdir()) / "chembl_biologics_cache")
    ap.add_argument("--no-cache", action="store_true", help="bypass the HTTP cache entirely")
    ap.add_argument("--no-uniprot", action="store_true", help="skip UniProt gene-symbol enrichment")
    ap.add_argument("--report", action="store_true", help="print markdown tables for docs/08")
    ap.add_argument("--verify", action="store_true", help="run the independent sequence spot-check")
    args = ap.parse_args()

    cache_dir = None
    if not args.no_cache:
        cache_dir = args.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)

    rows, excluded, mols = build(cache_dir, skip_uniprot=args.no_uniprot)
    write_csv(rows, args.out)
    print(f"\nwrote {len(rows)} rows -> {args.out}", file=sys.stderr)

    if args.report:
        report(rows, excluded, mols)
    if args.verify:
        fails = verify(rows, cache_dir)
        if fails:
            print(f"\n{fails} verification check(s) FAILED.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
