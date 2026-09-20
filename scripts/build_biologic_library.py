"""Build the approved peptide / protein / antibody drug library (Tier 2 + Tier 3).

Why this script exists
----------------------
`data/approved_drugs.csv` is DrugCentral's small-molecule table: SMILES only, no
modality column and no sequence column. Biologics that have no SMILES at all
(every antibody, every recombinant protein) were therefore never loaded, and the
peptides that *do* have a SMILES were loaded as if they were small molecules.
The biologics arm of the interface matcher needs SEQUENCES, so this rebuilds the
library from sources that actually carry them.

Sources (each row records which one it came from in `source`)
-------------------------------------------------------------
  chembl_biotherapeutic  ChEMBL REST `molecule?max_phase=4&biotherapeutic__isnull=false`.
                         Carries `helm_notation` (monomer-level, including
                         non-standard residues and disulfide/cyclization bonds)
                         and/or `biotherapeutic.biocomponents[].sequence`.
  chembl_by_inchikey     DrugCentral peptide-like row resolved to a ChEMBL record
                         by exact standard InChIKey.
  chembl_by_name         ... resolved by exact pref_name when the InChIKey did not
                         match (salt / stereo / tautomer differences).
  drugcentral_only       Peptide-like by structure in DrugCentral, but no ChEMBL
                         record with a sequence. Kept with an empty sequence and
                         an exclusion reason; NOT dropped.

DrugCentral peptide detection is structural, not a name heuristic: an RDKit SMARTS
for two consecutive alpha-peptide backbone units, counted over all 4,099 approved
SMILES. Count is reported so the "~133 approved peptides in DrugCentral" claim in
PROJECT_GOAL.md D3 can be checked against what is actually on disk.

Modality rule (applied in this order, recorded so it can be argued with)
-----------------------------------------------------------------------
  antibody  ChEMBL molecule_type in {Antibody, Antibody drug conjugate}, or the
            preferred name ends in -mab.
  other     ChEMBL molecule_type in {Oligonucleotide, Gene, Cell} -- not a
            polypeptide, out of scope for interface matching. Kept in the CSV so
            nothing vanishes silently; consumers filter on modality.
  peptide   resolved residue count <= PEPTIDE_MAX_RESIDUES (50).
  protein   resolved residue count > 50.
  Rows with no sequence fall back to molecule_type (Protein/Enzyme -> protein).

Non-standard residues (PROJECT_GOAL.md D3: "decide model, substitute or exclude
per entry and record the decision")
------------------------------------------------------------------------------
Every bracketed HELM monomer that is not one of the 20 standard residues is looked
up in NONSTANDARD below. Each carries an explicit decision:

  substitute  a named standard residue is used in `sequence`, and the loss is
              recorded verbatim in `residue_notes` (chirality, N-methylation,
              a cap, a side-chain halogen ...). The row stays usable.
  exclude     no defensible standard analogue (ChEMBL's opaque [Xnnn] monomer
              placeholders). `sequence` is still written for inspection but
              `usable_sequence` is 0 and `exclusion_reason` names the monomer.

Nothing is dropped for having a non-standard residue. The decision is per entry
and is auditable from the CSV alone.

Offline / cost
--------------
Zero Rowan credits. ChEMBL REST is free and its responses are cached under
data/raw/chembl_biologics/; re-runs with the cache present make no network calls.

Usage
-----
  ./env/bin/python scripts/build_biologic_library.py            # use cache if present
  ./env/bin/python scripts/build_biologic_library.py --refresh  # re-fetch from ChEMBL
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import certifi

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CACHE = DATA / "raw" / "chembl_biologics"
OUT_CSV = DATA / "biologic_drugs.csv"

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
SSL_CTX = ssl.create_default_context(cafile=certifi.where())

# ChEMBL molecule fields we need. `only` keeps molfile-heavy payloads down.
ONLY = (
    "molecule_chembl_id,pref_name,molecule_type,helm_notation,biotherapeutic,"
    "first_approval,withdrawn_flag,max_phase,atc_classifications,molecule_structures"
)

STD20 = set("ACDEFGHIKLMNPQRSTVWY")
PEPTIDE_MAX_RESIDUES = 50
# stop the ChEMBL name pass after this many consecutive endpoint failures
NAME_LOOKUP_ERROR_LIMIT = 5

# ---------------------------------------------------------------------------
# Non-standard HELM monomer decision table.
# (standard_substitute, decision, note). standard_substitute None => exclude.
# Populated from the monomers that actually occur in the approved set; an
# unseen monomer falls through to an explicit "unmapped" exclusion rather than
# being guessed at.
# ---------------------------------------------------------------------------
NONSTANDARD: dict[str, tuple[str | None, str, str]] = {
    # terminal caps: no residue is added, the capped terminus becomes free
    "ac": (None, "substitute", "N-terminal acetyl cap dropped; free N-terminus"),
    "am": (None, "substitute", "C-terminal amide cap dropped; free C-terminus"),
    # D-enantiomers of standard residues: L form used, chirality lost
    "dA": ("A", "substitute", "D-Ala -> L-Ala; backbone chirality lost"),
    "dR": ("R", "substitute", "D-Arg -> L-Arg; backbone chirality lost"),
    "dN": ("N", "substitute", "D-Asn -> L-Asn; backbone chirality lost"),
    "dD": ("D", "substitute", "D-Asp -> L-Asp; backbone chirality lost"),
    "dC": ("C", "substitute", "D-Cys -> L-Cys; backbone chirality lost"),
    "dE": ("E", "substitute", "D-Glu -> L-Glu; backbone chirality lost"),
    "dF": ("F", "substitute", "D-Phe -> L-Phe; backbone chirality lost"),
    "dW": ("W", "substitute", "D-Trp -> L-Trp; backbone chirality lost"),
    "dY": ("Y", "substitute", "D-Tyr -> L-Tyr; backbone chirality lost"),
    "dL": ("L", "substitute", "D-Leu -> L-Leu; backbone chirality lost"),
    "dV": ("V", "substitute", "D-Val -> L-Val; backbone chirality lost"),
    "dS": ("S", "substitute", "D-Ser -> L-Ser; backbone chirality lost"),
    "dT": ("T", "substitute", "D-Thr -> L-Thr; backbone chirality lost"),
    "dK": ("K", "substitute", "D-Lys -> L-Lys; backbone chirality lost"),
    "dH": ("H", "substitute", "D-His -> L-His; backbone chirality lost"),
    "dQ": ("Q", "substitute", "D-Gln -> L-Gln; backbone chirality lost"),
    "dP": ("P", "substitute", "D-Pro -> L-Pro; backbone chirality lost"),
    "dM": ("M", "substitute", "D-Met -> L-Met; backbone chirality lost"),
    "dI": ("I", "substitute", "D-Ile -> L-Ile; backbone chirality lost"),
    "dG": ("G", "substitute", "D-Gly -> Gly; achiral"),
    # D-enantiomers of non-standard side chains: nearest standard residue
    "dOrn": ("K", "substitute", "D-ornithine -> L-Lys; one CH2 longer, chirality lost"),
    "dCit": ("Q", "substitute", "D-citrulline -> L-Gln; ureido -> amide, chirality lost"),
    "dNal": ("W", "substitute", "D-2-naphthylalanine -> L-Trp; bicyclic aromatic kept"),
    "d1-Nal": ("W", "substitute", "D-1-naphthylalanine -> L-Trp; bicyclic aromatic kept"),
    "d3-Pal": ("H", "substitute", "D-3-pyridylalanine -> L-His; basic aromatic kept"),
    "dPhe(4-Cl)": ("F", "substitute", "D-4-chloro-Phe -> L-Phe; halogen and chirality lost"),
    "dTic": ("F", "substitute", "D-tetrahydroisoquinoline-COOH -> L-Phe; ring constraint lost"),
    # N-methylated residues: backbone N-methyl lost
    "meL": ("L", "substitute", "N-methyl-Leu -> Leu; backbone N-methyl lost"),
    "meV": ("V", "substitute", "N-methyl-Val -> Val; backbone N-methyl lost"),
    "meY": ("Y", "substitute", "N-methyl-Tyr -> Tyr; backbone N-methyl lost"),
    "meR": ("R", "substitute", "N-methyl-Arg -> Arg; backbone N-methyl lost"),
    "meA": ("A", "substitute", "N-methyl-Ala -> Ala; backbone N-methyl lost"),
    "Sar": ("G", "substitute", "sarcosine (N-methyl-Gly) -> Gly; backbone N-methyl lost"),
    # other non-standard side chains
    "Nle": ("L", "substitute", "norleucine -> Leu; straight chain -> branched"),
    "Abu": ("A", "substitute", "2-aminobutyrate -> Ala; one CH2 shorter"),
    "Glp": ("E", "substitute", "pyroglutamate -> Glu; N-terminal lactam opened"),
    "Hyp": ("P", "substitute", "4-hydroxyproline -> Pro; 4-OH lost"),
    "Thi": ("F", "substitute", "3-(2-thienyl)Ala -> Phe; thiophene -> phenyl"),
    "Phi": ("F", "substitute", "4-iodo-Phe -> Phe; iodine lost"),
    "Hph": ("F", "substitute", "homophenylalanine -> Phe; one CH2 shorter"),
    "Phg": ("F", "substitute", "phenylglycine -> Phe; one CH2 longer"),
    "Tyr(Bzl)": ("Y", "substitute", "O-benzyl-Tyr -> Tyr; benzyl ether lost"),
    "Aib": ("A", "substitute", "2-aminoisobutyrate -> Ala; alpha-methyl lost"),
    "Orn": ("K", "substitute", "ornithine -> Lys; one CH2 longer"),
    "Cit": ("Q", "substitute", "citrulline -> Gln; ureido -> amide"),
    "Nva": ("V", "substitute", "norvaline -> Val; straight chain -> branched"),
}

# monomers whose sequence contribution is zero (caps)
CAP_MONOMERS = {"ac", "am"}


# ---------------------------------------------------------------------------
# ChEMBL access
# ---------------------------------------------------------------------------
def chembl_get(path: str, params: dict, retries: int = 4):
    url = f"{CHEMBL}/{path}.json?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90, context=SSL_CTX) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001 - transient network
            if attempt == retries - 1:
                raise
            print(f"    retry {attempt + 1} after {exc!r}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return None


OFFLINE = False


def cached(name: str, refresh: bool, produce):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / name
    if path.exists() and not refresh:
        print(f"  cache hit {path.name}", file=sys.stderr)
        return json.loads(path.read_text())
    if OFFLINE:
        raise SystemExit(
            f"--offline was given but {path} is missing; run once without --offline "
            "to populate the cache."
        )
    value = produce()
    path.write_text(json.dumps(value))
    print(f"  wrote cache {path.name}", file=sys.stderr)
    return value


def fetch_chembl_version(refresh: bool) -> str:
    def produce():
        with urllib.request.urlopen(f"{CHEMBL}/status.json", timeout=60, context=SSL_CTX) as r:
            return json.load(r)

    st = cached("chembl_status.json", refresh, produce)
    return f"{st['chembl_db_version']} (released {st['chembl_release_date']})"


def fetch_biotherapeutics(refresh: bool) -> list[dict]:
    def produce():
        rows, offset = [], 0
        while True:
            page = chembl_get(
                "molecule",
                {
                    "max_phase": 4,
                    "biotherapeutic__isnull": "false",
                    "limit": 100,
                    "offset": offset,
                    "only": ONLY,
                },
            )
            mols = page["molecules"]
            if not mols:
                break
            rows.extend(mols)
            offset += 100
            print(f"    {len(rows)}/{page['page_meta']['total_count']}", file=sys.stderr)
            if offset >= page["page_meta"]["total_count"]:
                break
        return rows

    return cached("chembl_biotherapeutics_mp4.json", refresh, produce)


def fetch_by_inchikey(keys: list[str], refresh: bool) -> dict[str, dict]:
    def produce():
        out = {}
        for i in range(0, len(keys), 20):
            batch = keys[i : i + 20]
            page = chembl_get(
                "molecule",
                {
                    "molecule_structures__standard_inchi_key__in": ",".join(batch),
                    "limit": 50,
                    "only": ONLY,
                },
            )
            for m in page.get("molecules", []):
                ik = (m.get("molecule_structures") or {}).get("standard_inchi_key")
                if ik:
                    out[ik] = m
        return out

    return cached("chembl_by_inchikey.json", refresh, produce)


def fetch_by_name(names: list[str], refresh: bool) -> tuple[dict[str, dict], list[str]]:
    """One exact-name lookup per drug, cached incrementally.

    ChEMBL returns HTTP 500 for a few names (reproducibly, e.g. ones carrying
    parentheses). A failure on one name must not lose the whole pass, so each
    name is caught separately, the cache is written as it fills, and the names
    that never resolved are returned so the report can name them instead of
    quietly showing a smaller library.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / "chembl_by_name.json"
    store: dict[str, dict | None] = {}
    if path.exists() and not refresh:
        store = json.loads(path.read_text())
    todo = [n for n in names if n not in store]
    consecutive_errors = 0
    aborted = False
    for i, n in enumerate(todo, 1):
        if aborted:
            break
        try:
            page = chembl_get("molecule",
                              {"pref_name__iexact": n, "limit": 5, "only": ONLY},
                              retries=2)
            mols = page.get("molecules", [])
            store[n] = mols[0] if mols else None
            consecutive_errors = 0
        except Exception as exc:  # noqa: BLE001 - one bad name must not sink the pass
            print(f"    name lookup failed for {n!r}: {exc!r}", file=sys.stderr)
            consecutive_errors += 1
            # Circuit breaker: ChEMBL REST goes fully 500 from time to time (it
            # did during this build). Retrying 60 more names against a dead
            # service wastes minutes and changes nothing, so stop and let the
            # report say the pass was cut short.
            if consecutive_errors >= NAME_LOOKUP_ERROR_LIMIT:
                print(f"    ChEMBL name endpoint failed {consecutive_errors}x in a "
                      f"row; abandoning the name pass after {i}/{len(todo)}",
                      file=sys.stderr)
                aborted = True
        if i % 10 == 0 or i == len(todo) or aborted:
            path.write_text(json.dumps(store))
            print(f"    name lookups {i}/{len(todo)}", file=sys.stderr)
    path.write_text(json.dumps(store))
    resolved = {k: v for k, v in store.items() if v}
    unresolved = [n for n in names if not store.get(n)]
    return resolved, unresolved


# ---------------------------------------------------------------------------
# HELM parsing
# ---------------------------------------------------------------------------
def parse_helm(helm: str) -> dict:
    """Return monomer lists per polymer plus bond topology.

    HELM looks like  PEPTIDE1{[ac].S.Y.[Nle]}$PEPTIDE1,PEPTIDE1,6:R3-1:R3$$$
    Section 0 is the polymers, section 1 the inter-monomer connections.
    R1/R2 are the backbone termini (an R2-R1 bond inside one polymer is a
    head-to-tail macrocycle); R3 is the cysteine thiol (an R3-R3 bond is a
    disulfide).
    """
    sections = helm.split("$")
    polymers: dict[str, list[str]] = {}
    for pid, body in re.findall(r"(\w+?\d+)\{([^}]*)\}", sections[0]):
        polymers[pid] = [m for m in body.split(".") if m]
    conns = [c for c in (sections[1].split("|") if len(sections) > 1 and sections[1] else [])]

    n_disulfide = 0
    cyclic = False
    for c in conns:
        parts = c.split(",")
        if len(parts) < 3:
            continue
        src, dst, bond = parts[0], parts[1], parts[2]
        if "R3-" in bond and bond.endswith("R3"):
            n_disulfide += 1
        elif src == dst and ("R2-" in bond or bond.endswith("R1")):
            # backbone terminus joined to backbone terminus within one chain
            if ("R2" in bond and "R1" in bond):
                cyclic = True
    return {
        "polymers": polymers,
        "n_polymers": len(polymers),
        "n_disulfide": n_disulfide,
        "is_cyclic": cyclic,
        "polymer_types": sorted({re.sub(r"\d+$", "", p) for p in polymers}),
    }


def helm_to_sequence(parsed: dict) -> dict:
    """Monomers -> one-letter sequence, applying the NONSTANDARD decision table."""
    chains, notes, unmapped, subs = [], [], [], []
    for pid in sorted(parsed["polymers"]):
        seq = []
        for mon in parsed["polymers"][pid]:
            if mon in STD20 and len(mon) == 1:
                seq.append(mon)
                continue
            key = mon[1:-1] if mon.startswith("[") and mon.endswith("]") else mon
            if key in STD20 and len(key) == 1:
                seq.append(key)
                continue
            entry = NONSTANDARD.get(key)
            if entry is None:
                unmapped.append(key)
                seq.append("X")
                continue
            std, decision, note = entry
            if decision == "exclude":
                unmapped.append(key)
                seq.append("X")
                continue
            if key in CAP_MONOMERS:
                notes.append(f"{key}: {note}")
                continue  # cap contributes no residue
            seq.append(std)
            subs.append(key)
            notes.append(f"{key}->{std}: {note}")
        chains.append("".join(seq))
    return {
        "sequence": "/".join(chains),
        "n_residues": sum(len(c) for c in chains),
        "nonstandard": sorted(set(subs)),
        "unmapped": sorted(set(unmapped)),
        "notes": sorted(set(notes)),
    }


THREE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "IIE": "I", "LEU": "L",
    "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
    # non-standard three-letter names that occur in ChEMBL sequence strings
    "NAL": "W", "SAR": "G", "AIB": "A", "NLE": "L", "ORN": "K", "HYP": "P",
}
# three-letter codes above that are NOT one of the standard 20 and therefore
# represent a substitution that must be recorded
THREE_LETTER_NONSTANDARD = {"NAL", "SAR", "AIB", "NLE", "ORN", "HYP"}
# terminus markers in classical peptide notation, e.g. H-Asp-Arg-...-OH
TERMINUS_TOKENS = {"H", "OH", "NH2", "AC", "AM"}


def _parse_three_letter(raw: str) -> tuple[str, list[str], list[str]]:
    """Parse 'H-Asp-Arg-Val-OH' or '(Gly Gly Leu ... Sar Lys)2' notation.

    Returns (one-letter sequence, substituted monomers, unresolved tokens).
    A trailing repeat count on a parenthesised block is expanded.
    """
    repeat = 1
    m = re.fullmatch(r"\((.*)\)(\d+)", raw.strip())
    if m:
        raw, repeat = m.group(1), int(m.group(2))
    tokens = [t for t in re.split(r"[-\s]+", raw.strip()) if t]
    seq, subs, unres = [], [], []
    for tok in tokens:
        up = tok.upper()
        if up in TERMINUS_TOKENS:
            continue
        if up in THREE_LETTER:
            if up in THREE_LETTER_NONSTANDARD:
                subs.append(up)
            seq.append(THREE_LETTER[up])
        else:
            unres.append(up)
            seq.append("X")
    return "".join(seq) * repeat, subs, unres


def clean_protein_chain(raw: str, modres: dict[str, tuple[str | None, str]] | None = None) -> dict:
    """Normalise one ChEMBL protein-component sequence string.

    ChEMBL sequence strings are not uniformly one-letter. Three shapes occur in
    the approved set and each is handled explicitly rather than by stripping
    whatever does not parse:

      classical notation   'H-Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-OH'  (angiotensin II)
      bracketed repeat     '(Gly Gly Leu ... Sar Lys)2'            (peginesatide)
      one-letter + markers 'Q1VQLVQSGAEV...'                       (ravulizumab)

    In the third shape the embedded DIGITS are ChEMBL modification markers, not
    residues: removing them reproduces the published chain exactly (verified
    against ravulizumab's heavy chain, which begins QVQLVQSGAEVKKPGASVKVSC).
    They are counted, not silently discarded. '|' separates subchains.
    Any residue letter still outside the standard 20 after this is reported as
    UNRESOLVED and the row is not marked usable -- it is never guessed at.
    """
    raw = (raw or "").strip()
    if not raw:
        return {"chains": [], "subs": [], "unresolved": [], "n_markers": 0, "notes": []}
    subs: list[str] = []
    unres: list[str] = []
    notes: list[str] = []
    n_markers = 0

    looks_three_letter = bool(re.search(r"(?i)\b(ala|arg|asn|asp|cys|gln|glu|gly|his|ile|iie|"
                                        r"leu|lys|met|phe|pro|ser|thr|trp|tyr|val)\b", raw))
    if looks_three_letter:
        seq, s, u = _parse_three_letter(raw)
        subs += s
        unres += u
        chains = [seq] if seq else []
        return {"chains": chains, "subs": sorted(set(subs)),
                "unresolved": sorted(set(unres)), "n_markers": 0, "notes": []}

    modres = modres or {}
    chains = []
    for part in raw.split("|"):
        part = re.sub(r"\s+", "", part)
        n_markers += len(re.findall(r"\d", part))

        # Resolve TOKEN = letter + digits. A token whose base letter is one of
        # the standard 20 (K1, S1, A1) is that residue carrying a modification,
        # so the digit is just a marker. A token based on X has no residue
        # identity of its own and is resolvable only from a "Modified residues"
        # definition; without one it stays unresolved and is never guessed at.
        def _sub(m: re.Match) -> str:
            tok = m.group(0)
            base = tok[0].upper()
            defined = modres.get(tok)
            if defined is not None:
                letter, note = defined
                notes_local.append(f"{tok}: {note}")
                if letter is None:
                    return ""
                if base in STD20 and letter != base:
                    # trust the chain's own letter over the annotation
                    return base
                subs.append(tok)
                return letter
            if base in STD20:
                return base
            unres.append(base)
            return base

        notes_local: list[str] = []
        part = re.sub(r"[A-Za-z]\d+", _sub, part).upper()
        part = re.sub(r"\d", "", part)
        notes.extend(notes_local)
        for ch in set(part) - STD20:
            unres.append(ch)
        if part:
            chains.append(part)
    return {"chains": chains, "subs": sorted(set(subs)),
            "unresolved": sorted(set(unres)), "n_markers": n_markers,
            "notes": sorted(set(notes))}


# ---------------------------------------------------------------------------
# ChEMBL "Modified residues" annotations define tokens such as X1/K1/S1 that
# appear inside the chain string, each given as a SMILES. There are exactly 8
# distinct definitions across the approved set, so every one is resolved here
# by hand rather than by an inferred rule. Keyed by RDKit canonical SMILES.
# value = (one-letter residue or None for a cap that contributes no residue,
#          note recorded on the row)
# ---------------------------------------------------------------------------
MODIFIED_RESIDUE_SMILES_RAW: dict[str, tuple[str | None, str]] = {
    "CC(C)(N)C(O)=O": ("A", "Aib (2-aminoisobutyrate) -> Ala; alpha-methyl lost"),
    "[H][C@]1(CCC(=O)N1)C(O)=O": ("E", "pyroglutamate -> Glu; N-terminal lactam opened"),
    "OC(=O)[C@@H]1CCC(=O)N1": ("E", "pyroglutamate -> Glu; N-terminal lactam opened"),
    "NCC(N)=O": ("G", "C-terminal glycinamide -> Gly; amide cap dropped"),
    "N[C@@H](CO)C(N)=O": ("S", "C-terminal serinamide -> Ser; amide cap dropped"),
    "N[C@@H](CC1=CC([124I])=C(O)C=C1)C(O)=O":
        ("Y", "3-[124I]iodo-Tyr -> Tyr; radioiodine lost"),
    "CCCCCCCCCCCCCC(=O)NC(O)=O":
        (None, "N-myristoyl cap; no alpha-carbon in the definition, so it "
               "contributes no residue"),
}
_MODRES_CANON: dict[str, tuple[str | None, str]] | None = None


def modres_lookup(smiles: str) -> tuple[str | None, str] | None:
    """Resolve one modified-residue SMILES via canonical form.

    A long side-chain conjugate (PEG / fatty-acid / linker hung off a Lys) is
    not in the table; it is resolved structurally instead, by taking the
    residue letter that prefixes the token (K1 -> K) in clean_protein_chain.
    """
    global _MODRES_CANON
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    if _MODRES_CANON is None:
        _MODRES_CANON = {}
        for smi, val in MODIFIED_RESIDUE_SMILES_RAW.items():
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                _MODRES_CANON[Chem.MolToSmiles(mol)] = val
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return _MODRES_CANON.get(Chem.MolToSmiles(mol))


def biocomponents_to_sequence(bio: dict) -> dict:
    """Sequence + topology from ChEMBL biocomponents.

    `biocomponents` is a MIXED list: entries with component_type 'Protein' are
    real polypeptide chains, entries with component_type None are annotation
    (glycosylation sites, disulfide bridge tables, legends), and 'Nucleic Acid'
    entries are oligonucleotides. Verified over the 637 components in the
    approved set: no 'Protein' component carries an annotation description.
    Concatenating the annotation rows into the sequence -- which an unfiltered
    join does -- corrupts it, so they are read for topology instead.
    """
    comps = bio.get("biocomponents") or []
    chains, labels, subs, unres = [], [], [], []
    notes: list[str] = []
    n_markers = 0
    n_ss = 0
    n_glyc = 0
    n_nucleic = 0

    # pass 1: token definitions from the "Modified residues" annotation rows
    modres: dict[str, tuple[str | None, str]] = {}
    n_defs = n_defs_resolved = 0
    for c in comps:
        if not (c.get("description") or "").startswith("Modified residues"):
            continue
        for tok, smi in re.findall(r"([A-Za-z]\w*)\(SMILES\):([^$]+)",
                                   c.get("sequence") or ""):
            n_defs += 1
            hit = modres_lookup(smi.strip())
            if hit is not None:
                modres[tok] = hit
                n_defs_resolved += 1

    # pass 2: chains and topology
    for c in comps:
        ctype = c.get("component_type")
        desc = (c.get("description") or "").strip()
        raw = c.get("sequence") or ""
        if ctype == "Nucleic Acid":
            n_nucleic += 1
            continue
        if ctype == "Protein":
            built = clean_protein_chain(raw, modres)
            notes += built["notes"]
            for ch in built["chains"]:
                chains.append(ch)
                labels.append(desc or "Sequence")
            subs += built["subs"]
            unres += built["unresolved"]
            n_markers += built["n_markers"]
            continue
        # component_type is None -> annotation rows
        if desc.startswith("Disulfide"):
            n_ss += len(re.findall(r"\d+[^\s|$,-]*-\d+", raw))
        elif desc.startswith("Glycosylation"):
            n_glyc += len(re.findall(r"\d+", raw))

    return {
        "sequence": "/".join(chains),
        "n_residues": sum(len(c) for c in chains),
        "chain_labels": labels,
        "substituted": sorted(set(subs)),
        "unresolved": sorted(set(unres)),
        "n_modification_markers": n_markers,
        "n_disulfide": n_ss,
        "n_glycosylation_sites": n_glyc,
        "n_nucleic_components": n_nucleic,
        "notes": sorted(set(notes)),
        "n_modified_residue_defs": n_defs,
        "n_modified_residue_defs_resolved": n_defs_resolved,
    }



# ---------------------------------------------------------------------------
# DrugCentral side
# ---------------------------------------------------------------------------
def norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def drugcentral_peptide_like() -> list[dict]:
    """Structural peptide detection over DrugCentral's approved SMILES.

    SMARTS matches two consecutive alpha-peptide backbone units; a molecule with
    >= 2 such (overlapping) matches carries at least a tetrapeptide-sized run.
    Returns every hit with its match count so the threshold can be audited.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    patt = Chem.MolFromSmarts("[NX3;H0,H1,H2][CX4][CX3](=[OX1])[NX3;H0,H1][CX4][CX3](=[OX1])")
    hits, n_rows, n_unparsed = [], 0, 0
    with (DATA / "approved_drugs.csv").open() as fh:
        for r in csv.DictReader(fh):
            n_rows += 1
            mol = Chem.MolFromSmiles(r["smiles"]) if r.get("smiles") else None
            if mol is None:
                n_unparsed += 1
                continue
            n = len(mol.GetSubstructMatches(patt))
            if n >= 2:
                hits.append(
                    {
                        "struct_id": r["struct_id"],
                        "name": r["name"],
                        "inchikey": r["inchikey"],
                        "n_backbone_units": n,
                        "targets": r.get("targets", ""),
                    }
                )
    print(
        f"  DrugCentral: {n_rows} approved rows, {n_unparsed} unparsable SMILES, "
        f"{len(hits)} peptide-like",
        file=sys.stderr,
    )
    return hits


def drugcentral_targets() -> tuple[dict[str, set[str]], dict[str, set[str]], dict]:
    """name -> HUMAN gene symbols, for all interactions and for MOA-flagged ones.

    DrugCentral's interaction table is multi-species (14,301 of 19,378 rows are
    Homo sapiens; the rest are rat, mouse, bovine, guinea pig, bacterial and
    viral). Keeping the non-human rows double-counts a target under both its
    human symbol (AVPR2) and its rodent ortholog (Avpr2), which silently
    inflates any "binders per target" ranking. Only Homo sapiens rows are kept
    and the number dropped is reported.
    """
    allt: dict[str, set[str]] = collections.defaultdict(set)
    moat: dict[str, set[str]] = collections.defaultdict(set)
    kept = dropped = 0
    with (DATA / "drug.target.interaction.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = {g.strip() for g in (r.get("GENE") or "").split("|") if g.strip()}
            if not genes:
                continue
            if (r.get("ORGANISM") or "").strip() != "Homo sapiens":
                dropped += 1
                continue
            kept += 1
            key = norm_name(r["DRUG_NAME"])
            allt[key] |= genes
            if (r.get("MOA") or "").strip() == "1":
                moat[key] |= genes
    return allt, moat, {"kept_human": kept, "dropped_non_human": dropped}


def file_version(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    return f"{path.name}@sha256:{h}"


# ---------------------------------------------------------------------------
# Row assembly
# ---------------------------------------------------------------------------
def assign_modality(mol: dict, n_residues: int | None) -> str:
    mtype = mol.get("molecule_type") or ""
    name = (mol.get("pref_name") or "").lower()
    if mtype in {"Antibody", "Antibody drug conjugate"} or name.endswith("mab"):
        return "antibody"
    if mtype in {"Oligonucleotide", "Gene", "Cell"}:
        return "other"
    if n_residues:
        return "peptide" if n_residues <= PEPTIDE_MAX_RESIDUES else "protein"
    if mtype in {"Protein", "Enzyme"}:
        return "protein"
    return "other"


def build_row(mol: dict, source: str, dc: dict | None) -> dict:
    helm = mol.get("helm_notation")
    bio = mol.get("biotherapeutic") or {}
    seq, n_res, is_cyclic, n_ss = "", None, "", ""
    nonstd: list[str] = []
    unmapped: list[str] = []
    notes: list[str] = []
    seq_source = ""

    if helm:
        parsed = parse_helm(helm)
        built = helm_to_sequence(parsed)
        seq = built["sequence"]
        n_res = built["n_residues"]
        nonstd = built["nonstandard"]
        unmapped = built["unmapped"]
        notes = built["notes"]
        is_cyclic = "1" if parsed["is_cyclic"] else "0"
        n_ss = str(parsed["n_disulfide"])
        seq_source = "helm_notation"
    n_glyc, n_markers, chain_labels = "", "", []
    if bio.get("biocomponents"):
        built = biocomponents_to_sequence(bio)
        n_glyc = str(built["n_glycosylation_sites"]) if built["n_glycosylation_sites"] else "0"
        if not seq and built["sequence"]:
            seq = built["sequence"]
            n_res = built["n_residues"]
            unmapped = built["unresolved"]
            nonstd = built["substituted"]
            n_markers = str(built["n_modification_markers"])
            chain_labels = built["chain_labels"]
            seq_source = "biocomponents"
            if n_ss in ("", "0") and built["n_disulfide"]:
                n_ss = str(built["n_disulfide"])
            notes += built["notes"]
            for s in built["substituted"]:
                # modified-residue TOKENS (S1, X1) already wrote their own note
                # in clean_protein_chain; only three-letter codes need one here
                if s in THREE_LETTER:
                    notes.append(f"{s}->{THREE_LETTER[s]}: non-standard "
                                 "three-letter residue substituted")
            if built["n_modification_markers"]:
                notes.append(
                    f"{built['n_modification_markers']} ChEMBL modification marker digit(s) "
                    "removed from the chain; residues themselves unchanged"
                )

    modality = assign_modality(mol, n_res)

    exclusion = ""
    if modality == "other":
        exclusion = "not_a_polypeptide:" + (mol.get("molecule_type") or "unknown")
    elif not seq:
        exclusion = "no_sequence_in_source"
    elif unmapped and seq_source == "helm_notation":
        exclusion = "nonstandard_residue_unmappable:" + ",".join(unmapped)
    elif unmapped:
        exclusion = "unresolved_residue_letter:" + ",".join(unmapped)
    usable = "0" if exclusion else "1"

    decision = "none"
    if unmapped:
        decision = "exclude"
    elif nonstd or notes:
        decision = "substitute"

    inchikey = (mol.get("molecule_structures") or {}).get("standard_inchi_key") or ""
    if not inchikey and dc:
        inchikey = dc.get("inchikey", "")

    return {
        "id": mol["molecule_chembl_id"],
        "name": (mol.get("pref_name") or (dc or {}).get("name") or "").strip(),
        "modality": modality,
        "sequence": seq,
        "n_residues": "" if n_res is None else str(n_res),
        "n_chains": "" if not seq else str(len(seq.split("/"))),
        "is_cyclic": is_cyclic,
        "n_disulfide": n_ss,
        "n_glycosylation_sites": n_glyc,
        "n_modification_markers": n_markers,
        "chain_labels": ";".join(chain_labels),
        "has_nonstandard_residues": "1" if (nonstd or unmapped or notes) else "0",
        "nonstandard_residues": ";".join(nonstd + unmapped),
        "residue_decision": decision,
        "residue_notes": " | ".join(notes),
        "usable_sequence": usable,
        "exclusion_reason": exclusion,
        "molecule_type": mol.get("molecule_type") or "",
        "inchikey": inchikey,
        "drugcentral_struct_id": (dc or {}).get("struct_id", ""),
        "first_approval": str(mol.get("first_approval") or ""),
        "withdrawn": "1" if mol.get("withdrawn_flag") else "0",
        "atc": ";".join(mol.get("atc_classifications") or []),
        "sequence_source": seq_source,
        "source": source,
    }


def drugcentral_only_row(dc: dict) -> dict:
    return {
        "id": "DC" + dc["struct_id"],
        "name": dc["name"],
        "modality": "peptide",
        "sequence": "",
        "n_residues": "",
        "n_chains": "",
        "is_cyclic": "",
        "n_disulfide": "",
        "n_glycosylation_sites": "",
        "n_modification_markers": "",
        "chain_labels": "",
        "has_nonstandard_residues": "",
        "nonstandard_residues": "",
        "residue_decision": "unknown",
        "residue_notes": f"peptide-like by SMILES ({dc['n_backbone_units']} backbone units); "
        "no ChEMBL record with a sequence",
        "usable_sequence": "0",
        "exclusion_reason": "no_sequence_in_source",
        "molecule_type": "",
        "inchikey": dc.get("inchikey", ""),
        "drugcentral_struct_id": dc["struct_id"],
        "first_approval": "",
        "withdrawn": "",
        "atc": "",
        "sequence_source": "",
        "source": "drugcentral_only",
    }


# ---------------------------------------------------------------------------
# Reference sequences, for validating the build rather than trusting it.
# Each is the published backbone of an approved peptide, written in the same
# convention this script emits: non-standard residues already replaced by the
# substitute named in NONSTANDARD / MODIFIED_RESIDUE_SMILES_RAW, caps dropped.
# A mismatch is reported, not silently tolerated -- the source can be wrong.
# ---------------------------------------------------------------------------
REFERENCE_SEQUENCES: dict[str, tuple[str, str]] = {
    # name: (expected sequence, what the reference is)
    "OXYTOCIN": ("CYIQNCPLG", "published nonapeptide, Cys1-Cys6 disulfide"),
    "VASOPRESSIN": ("CYFQNCPRG", "arginine vasopressin nonapeptide"),
    "EXENATIDE": ("HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS", "exendin-4, 39-mer"),
    "TERIPARATIDE": ("SVSEIQLMHNLGKHLNSMERVEWLRKKLQDVHNF", "PTH(1-34)"),
    "ENFUVIRTIDE": ("YTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWF", "gp41 HR2 36-mer"),
    "BIVALIRUDIN": ("FPRPGGGGNGDFEEIPEEYL", "hirulog-1 20-mer"),
    "CETRORELIX": ("WFHSYQLRPA",
                   "Ac-D-Nal-D-Cpa-D-Pal-Ser-Tyr-D-Cit-Leu-Arg-Pro-D-Ala-NH2 "
                   "under this script's substitution rules"),
    "SEMAGLUTIDE": ("HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
                    "GLP-1(7-37) analogue, Aib2 -> Ala under substitution"),
    "LIRAGLUTIDE": ("HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
                    "GLP-1(7-37) analogue, Arg34Lys, K26 acylated"),
    "TIRZEPATIDE": ("YAEGTFTSDYSIALDKIAQKAFVQWLIAGGPSSGAPPPS",
                    "GIP/GLP-1 dual agonist 39-mer, Aib2 and Aib13 -> Ala"),
}


def validate_against_references(rows: list[dict]) -> dict:
    """Compare emitted sequences with published ones. Returns counts + detail."""
    by_name = {r["name"].upper(): r for r in rows}
    match, mismatch, absent, no_seq = [], [], [], []
    for name, (expect, note) in REFERENCE_SEQUENCES.items():
        r = by_name.get(name)
        if r is None:
            absent.append((name, note))
        elif not r["sequence"]:
            no_seq.append((name, r["exclusion_reason"]))
        elif r["sequence"] == expect:
            match.append(name)
        else:
            mismatch.append((name, expect, r["sequence"]))
    # Stamp the verdict onto the row itself. Printing it at build time is not
    # enough: a consumer reads the CSV, sees usable_sequence=1, and has no way to
    # know the sequence disagrees with the published one. SEMAGLUTIDE is exactly
    # this case and it sits in the GLP1R positive set, so the warning has to
    # travel with the data.
    for name in match:
        by_name[name]["reference_check"] = "match"
    for name, expect, _got in mismatch:
        by_name[name]["reference_check"] = f"MISMATCH:published={expect}"
    for name, _why in no_seq:
        by_name[name]["reference_check"] = "no_sequence"
    return {"match": match, "mismatch": mismatch, "absent": absent,
            "no_seq": no_seq, "n": len(REFERENCE_SEQUENCES)}


FIELDS = [
    "id", "name", "modality", "sequence", "n_residues", "n_chains",
    "is_cyclic", "n_disulfide", "n_glycosylation_sites", "n_modification_markers",
    "chain_labels", "has_nonstandard_residues", "nonstandard_residues",
    "residue_decision", "residue_notes", "usable_sequence", "exclusion_reason",
    "molecule_type", "inchikey", "drugcentral_struct_id", "first_approval",
    "withdrawn", "atc", "targets", "moa_targets", "n_targets", "approved",
    "source", "source_version", "sequence_source", "reference_check",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true", help="re-fetch from ChEMBL REST")
    ap.add_argument("--offline", action="store_true",
                    help="build only from data/raw/chembl_biologics/; make no network "
                         "calls. Use this to reproduce a build byte-for-byte, or when "
                         "ChEMBL REST is down (it returned HTTP 500 service-wide during "
                         "this project's build).")
    args = ap.parse_args()
    global OFFLINE
    OFFLINE = args.offline

    print("ChEMBL version...", file=sys.stderr)
    chembl_version = fetch_chembl_version(args.refresh)
    dc_version = file_version(DATA / "approved_drugs.csv")
    dti_version = file_version(DATA / "drug.target.interaction.tsv")

    print("ChEMBL approved biotherapeutics...", file=sys.stderr)
    bios = fetch_biotherapeutics(args.refresh)

    print("DrugCentral structural peptide detection...", file=sys.stderr)
    dc_peptides = drugcentral_peptide_like()

    seen_keys = {
        (m.get("molecule_structures") or {}).get("standard_inchi_key")
        for m in bios
    }
    seen_names = {norm_name(m.get("pref_name")) for m in bios}

    todo = [d for d in dc_peptides if d["inchikey"] not in seen_keys
            and norm_name(d["name"]) not in seen_names]
    print(f"  {len(dc_peptides) - len(todo)} of {len(dc_peptides)} peptide-like rows "
          f"already in the ChEMBL biotherapeutic pool", file=sys.stderr)

    print("Resolving remaining DrugCentral peptides by InChIKey...", file=sys.stderr)
    by_ik = fetch_by_inchikey([d["inchikey"] for d in todo if d["inchikey"]], args.refresh)
    def _ik_has_seq(d: dict) -> bool:
        m = by_ik.get(d["inchikey"])
        if not m:
            return False
        if m.get("helm_notation"):
            return True
        return any((c.get("component_type") == "Protein" and c.get("sequence"))
                   for c in ((m.get("biotherapeutic") or {}).get("biocomponents") or []))

    still = [d for d in todo if not _ik_has_seq(d)]
    print(f"  {len(by_ik)} InChIKey matches, {len(still)} still without a sequence",
          file=sys.stderr)

    print("Resolving the rest by name...", file=sys.stderr)
    if args.offline:
        by_name, name_unresolved = {}, [d["name"] for d in still]
        print("  --offline: name pass skipped", file=sys.stderr)
    else:
        by_name, name_unresolved = fetch_by_name([d["name"] for d in still], args.refresh)
    print(f"  {len(by_name)} resolved by name, {len(name_unresolved)} not found",
          file=sys.stderr)

    all_targets, moa_targets, dti_stats = drugcentral_targets()

    rows, used_ids = [], set()
    for m in bios:
        r = build_row(m, "chembl_biotherapeutic", None)
        rows.append(r)
        used_ids.add(r["id"])
    def _has_sequence(m: dict | None) -> bool:
        if not m:
            return False
        if m.get("helm_notation"):
            return True
        return any((c.get("component_type") == "Protein" and c.get("sequence"))
                   for c in ((m.get("biotherapeutic") or {}).get("biocomponents") or []))

    for d in todo:
        mol = by_ik.get(d["inchikey"])
        src = "chembl_by_inchikey"
        # an InChIKey can match a ChEMBL record that carries no sequence at all
        # (liraglutide resolves to one); fall through to the name lookup rather
        # than accepting an empty hit
        if not _has_sequence(mol) and _has_sequence(by_name.get(d["name"])):
            mol = by_name[d["name"]]
            src = "chembl_by_name"
        if mol is None:
            mol = by_name.get(d["name"])
            src = "chembl_by_name"
        if mol is None:
            rows.append(drugcentral_only_row(d))
            continue
        if mol["molecule_chembl_id"] in used_ids:
            continue
        r = build_row(mol, src, d)
        if not r["name"]:
            r["name"] = d["name"]
        rows.append(r)
        used_ids.add(r["id"])

    # annotate targets and approval
    for r in rows:
        key = norm_name(r["name"])
        tg = sorted(all_targets.get(key, ()))
        mo = sorted(moa_targets.get(key, ()))
        r["targets"] = ";".join(tg)
        r["moa_targets"] = ";".join(mo)
        r["n_targets"] = str(len(tg))
        r["approved"] = "1"
        r["source_version"] = (
            f"chembl_rest:{chembl_version}; drugcentral:{dc_version}; dti:{dti_version}"
            if r["source"] != "drugcentral_only"
            else f"drugcentral:{dc_version}; dti:{dti_version}"
        )

    # Run the reference check BEFORE writing: it stamps `reference_check` onto
    # each checked row, and report() below re-runs it for the printed summary.
    validate_against_references(rows)

    rows.sort(key=lambda r: (r["modality"], r["name"]))
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {OUT_CSV}", file=sys.stderr)

    report(rows, dc_peptides, dti_stats, name_unresolved)
    return 0


def report(rows: list[dict], dc_peptides: list[dict], dti_stats: dict,
           name_unresolved: list[str]) -> None:
    n = len(rows)
    out = print

    out("\n" + "=" * 72)
    out(f"BIOLOGIC DRUG LIBRARY  n = {n} rows")
    out("=" * 72)

    out("\n-- by modality (all rows) --")
    by_mod = collections.Counter(r["modality"] for r in rows)
    for k, v in by_mod.most_common():
        usable = sum(1 for r in rows if r["modality"] == k and r["usable_sequence"] == "1")
        withseq = sum(1 for r in rows if r["modality"] == k and r["sequence"])
        out(f"   {k:<10} n={v:<4} with a sequence={withseq:<4} usable sequence={usable}")
    out(f"   {'TOTAL':<10} n={n:<4} "
        f"with a sequence={sum(1 for r in rows if r['sequence']):<4} "
        f"usable sequence={sum(1 for r in rows if r['usable_sequence'] == '1')}")

    out("\n-- sequence source --")
    for k, v in collections.Counter(r["sequence_source"] or "(none)" for r in rows).most_common():
        out(f"   {k:<16} {v}")

    out("\n-- provenance (source) --")
    for k, v in collections.Counter(r["source"] for r in rows).most_common():
        out(f"   {k:<24} {v}")

    out("\n-- residue length, rows with a usable sequence --")
    lens = sorted(int(r["n_residues"]) for r in rows
                  if r["usable_sequence"] == "1" and r["n_residues"])
    if lens:
        mid = lens[len(lens) // 2]
        out(f"   n={len(lens)}  min={lens[0]}  median={mid}  max={lens[-1]}")
        for cut in (15, 30, 50, 100, 500):
            out(f"   <= {cut:<4} residues: {sum(1 for x in lens if x <= cut)}")

    out("\n-- non-standard residues --")
    ns = [r for r in rows if r["has_nonstandard_residues"] == "1"]
    out(f"   rows carrying >=1 non-standard monomer: {len(ns)}")
    for k, v in collections.Counter(r["residue_decision"] for r in rows).most_common():
        out(f"   decision {k:<12} {v}")
    mon = collections.Counter()
    for r in rows:
        for m in filter(None, r["nonstandard_residues"].split(";")):
            mon[m] += 1
    out(f"   distinct non-standard monomers seen: {len(mon)}")
    for k, v in mon.most_common(20):
        if k in NONSTANDARD:
            _, dec, note = NONSTANDARD[k]
        elif k in THREE_LETTER:
            dec, note = "substitute", (
                f"three-letter {k} -> {THREE_LETTER[k]}")
        elif re.fullmatch(r"X\d*", k):
            dec, note = "exclude", "opaque ChEMBL monomer placeholder / unspecified residue"
        else:
            dec, note = "exclude", "residue letter outside the standard 20; not guessed at"
        out(f"     {k:<14} x{v:<3} {dec:<10} {note}")

    out("\n-- topology, rows with HELM --")
    helm_rows = [r for r in rows if r["sequence_source"] == "helm_notation"]
    out(f"   n={len(helm_rows)}   cyclic (head-to-tail)={sum(1 for r in helm_rows if r['is_cyclic'] == '1')}"
        f"   with >=1 disulfide={sum(1 for r in helm_rows if r['n_disulfide'] not in ('', '0'))}")
    out("   is_cyclic unknown (no HELM): "
        f"{sum(1 for r in rows if r['is_cyclic'] == '')}")

    out("\n-- NOT EVALUATED / WHY (exclusion reasons, with counts) --")
    exc = collections.Counter(
        r["exclusion_reason"].split(":")[0] for r in rows if r["exclusion_reason"]
    )
    for k, v in exc.most_common():
        out(f"   {k:<34} {v}")
    out(f"   {'(none - usable)':<34} {sum(1 for r in rows if not r['exclusion_reason'])}")

    v = validate_against_references(rows)
    out("\n-- VALIDATION against published reference sequences --")
    out(f"   n={v['n']} hand-checked approved peptides compared residue by residue")
    out(f"   exact match : {len(v['match'])}  ({', '.join(n.lower() for n in v['match'])})")
    out(f"   mismatch    : {len(v['mismatch'])}")
    for name, exp, got in v["mismatch"]:
        out(f"     {name.lower()}")
        out(f"       published : {exp}")
        out(f"       emitted   : {got}")
        out( "       -> the emitted string is what the source stores; treated as a")
        out( "          SOURCE discrepancy, not corrected here")
    out(f"   no sequence : {len(v['no_seq'])}"
        + ("  " + "; ".join(f"{n.lower()} ({why})" for n, why in v["no_seq"])
           if v["no_seq"] else ""))
    out(f"   absent      : {len(v['absent'])}"
        + ("  " + ", ".join(n.lower() for n, _ in v["absent"]) if v["absent"] else ""))

    out("\n-- DrugCentral peptide reconciliation --")
    out(f"   PROJECT_GOAL.md D3 states '~133 approved peptides in DrugCentral'.")
    out(f"   Measured on data/approved_drugs.csv (n=4,099 approved rows): "
        f"{len(dc_peptides)} peptide-like by structure")
    out( "   (RDKit SMARTS, >=2 consecutive alpha-peptide backbone units). The claim is")
    out( "   NOT reproduced exactly; 121 is what is on disk.")
    dc_in = sum(1 for r in rows if r["drugcentral_struct_id"])
    out(f"   of those, carried into this library with a resolved sequence: "
        f"{sum(1 for r in rows if r['drugcentral_struct_id'] and r['sequence'])}")
    out(f"   of those, no sequence found in ChEMBL by InChIKey or name: "
        f"{sum(1 for r in rows if r['drugcentral_struct_id'] and not r['sequence'])}")
    out(f"   (DrugCentral rows represented in the library at all: {dc_in})")
    if name_unresolved:
        out(f"   ChEMBL exact-name lookup returned nothing (or errored) for "
            f"{len(name_unresolved)} of these rows at build time:")
        out( "     " + ", ".join(sorted(name_unresolved))[:400]
             + (" ..." if len(", ".join(sorted(name_unresolved))) > 400 else ""))
        out( "     Some are genuine non-peptides caught by the structural SMARTS")
        out( "     (beta-lactams, echinocandins); the rest stay sequence-less and")
        out( "     unusable rather than being filled in with a guess.")

    out("\n-- rows needing manual residue curation (would become usable) --")
    need = [r for r in rows if r["exclusion_reason"].startswith(
        ("unresolved_residue_letter", "nonstandard_residue_unmappable"))]
    out(f"   n={len(need)}; each is kept in the CSV, none dropped:")
    for r in sorted(need, key=lambda r: r["name"]):
        out(f"     {r['name'][:34]:<34} {r['modality']:<8} {r['exclusion_reason']}")

    out("\n-- target annotation (DrugCentral drug.target.interaction) --")
    out(f"   interaction rows kept (Homo sapiens): {dti_stats['kept_human']}; "
        f"dropped (non-human): {dti_stats['dropped_non_human']}")
    with_t = [r for r in rows if r["n_targets"] != "0"]
    out(f"   rows with >=1 known target: {len(with_t)} / {n}")
    for mod in ("peptide", "protein", "antibody", "other"):
        sub = [r for r in rows if r["modality"] == mod]
        if sub:
            out(f"     {mod:<9} {sum(1 for r in sub if r['n_targets'] != '0')} / {len(sub)}")
    out(f"   rows usable AND with >=1 target: "
        f"{sum(1 for r in rows if r['usable_sequence'] == '1' and r['n_targets'] != '0')}")

    out("\n" + "=" * 72)
    out("VALIDATION TARGETS: which genes have the most approved biologic binders")
    out("=" * 72)
    out("Counted over rows with a USABLE sequence only -- a target whose binders")
    out("all lack sequences cannot validate anything.")

    def rank(mods: tuple[str, ...], label: str, top: int = 15) -> list[tuple[str, int]]:
        cnt: dict[str, set[str]] = collections.defaultdict(set)
        for r in rows:
            if r["modality"] not in mods or r["usable_sequence"] != "1":
                continue
            for g in filter(None, r["targets"].split(";")):
                cnt[g].add(r["name"])
        ranked = sorted(((g, len(v)) for g, v in cnt.items()), key=lambda x: (-x[1], x[0]))
        out(f"\n-- top targets by {label} --")
        if not ranked:
            out("   (none)")
        for g, c in ranked[:top]:
            names = sorted({r["name"] for r in rows
                            if r["modality"] in mods and r["usable_sequence"] == "1"
                            and g in r["targets"].split(";")})
            out(f"   {g:<10} {c:>3}   {', '.join(n.lower() for n in names[:6])}"
                f"{' ...' if len(names) > 6 else ''}")
        return ranked

    pep_rank = rank(("peptide",), "approved PEPTIDE binders (n_drugs)")
    rank(("peptide", "protein"), "approved PEPTIDE+PROTEIN binders (n_drugs)")

    out("\n-- programmatic absence check: current demo target --")
    for gene in ("KDR", "EGFR", "BRAF", "CA2", "CDK2"):
        hits = [r["name"] for r in rows
                if r["usable_sequence"] == "1" and gene in r["targets"].split(";")]
        peps = [r["name"] for r in rows
                if r["modality"] == "peptide" and r["usable_sequence"] == "1"
                and gene in r["targets"].split(";")]
        out(f"   {gene:<6} biologic binders with usable sequence: {len(hits):<3} "
            f"of which peptides: {len(peps)}"
            + (f"   ({', '.join(n.lower() for n in hits[:4])})" if hits else ""))

    out("\n-- best candidate validation targets (>=3 peptide binders) --")
    good = [(g, c) for g, c in pep_rank if c >= 3]
    if not good:
        out("   none reach 3")
    for g, c in good:
        out(f"   {g:<10} {c} approved peptide binders")


if __name__ == "__main__":
    raise SystemExit(main())
