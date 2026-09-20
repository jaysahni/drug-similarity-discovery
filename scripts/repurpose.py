"""THE PIPELINE. Binding site -> designed binders -> interface signature -> approved
drugs that engage the same site.

    disease/target + pocket
        -> BoltzGen designs a few binders against that pocket      (scripts/boltzgen_signature.py)
        -> consensus INTERFACE SIGNATURE: which target residues a good binder engages
        -> co-fold each approved drug WITH the target (Boltz-2 via Rowan),
           unconstrained by default, then score the pose against that site
        -> extract each drug's target-side contacts
        -> rank drugs by how much of the signature they engage
        -> check: do the target's KNOWN binders come out on top?

This is PROJECT_GOAL.md's v0.2 pipeline. The match is computed in TARGET-SIDE
coordinates - which residues are engaged - not in chemical space, because a designed
miniprotein and a small molecule share no chemistry to compare. That is the whole
point: scripts/match_candidates.py ranks the same drugs by chemical similarity and
buries sunitinib at 464 and pazopanib at 764, both real binders of this target.

NO AFFINITY IS USED ANYWHERE IN THE RANKING. Boltz-2 can emit an affinity score and
this pipeline deliberately does not read it: PROJECT_GOAL.md 4.4 removes affinity from
the ranking path and G8 bans the vocabulary. Ranking is interface overlap, gated on
structural confidence.

Contacts for designs and for drugs come from the SAME function in interfaces.py
(PROJECT_GOAL.md F5) - if the two drifted apart, every score would be meaningless.

BOTH MODALITIES GO THROUGH THE SAME PATH. A shortlist row carries a `modality`:
a small molecule is a SMILES in the ligand slot, a peptide or protein binder is an
extra CHAIN in initial_protein_sequences next to the target. An absent modality
means small molecule, so every shortlist written before this existed is unchanged.
Only shortlist and submit know the difference; the scoring is in target-side
residues and does not care what the binder is made of. See docs/07-BIOLOGICS.md
for what is implemented, what it is expected to cost, and what is still unverified.

Phases, each resumable; state in results/pipeline/<slug>/repurpose_state.json:
    shortlist   pick the drugs to co-fold, with positive controls and decoys
    submit      send one co-folding job per drug. UNCONSTRAINED by default:
                the pocket is used to SCORE the pose afterwards, not to
                condition the folding. --constrain-pocket changes that
    collect     poll, download poses
    score       extract contacts, score against the signature, rank, validate

Usage:
    ./env/bin/python scripts/repurpose.py shortlist --n-decoys 20
    ./env/bin/python scripts/repurpose.py submit --max-credits 120
    ./env/bin/python scripts/repurpose.py collect
    ./env/bin/python scripts/repurpose.py score --designs v4

    # biologics, costing nothing:
    ./env/bin/python scripts/repurpose.py shortlist --biologics --target GLP1R \
        --out /tmp/bio_shortlist.json
    ./env/bin/python scripts/repurpose.py submit --dry-run \
        --shortlist /tmp/bio_shortlist.json --max-credits 12
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interfaces  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEED = 0
CUTOFF = 4.5
CONFIDENCE_GATE = 0.5      # documented below; records below it are kept and flagged


def load_env():
    env = ROOT / ".env"
    if env.exists():
        for line in env.read_text().splitlines():
            if "=" in line and not line.strip().startswith("#"):
                k, _, v = line.partition("=")
                os.environ.setdefault(k.strip(), v.strip())


def paths(slug):
    p = ROOT / "results" / "pipeline" / slug
    return {"dir": p, "state": p / "repurpose_state.json",
            "poses": p / "cofold_poses", "target": p / "target"}


def state_load(slug):
    f = paths(slug)["state"]
    return json.loads(f.read_text()) if f.exists() else {}


def state_save(slug, st):
    paths(slug)["state"].write_text(json.dumps(st, indent=1))


def target_info(slug):
    t = paths(slug)["target"]
    pocket = json.loads((t / "pocket.json").read_text())
    fa = next(x for x in (t / "sequence_kinase_domain.fasta", t / "sequence.fasta") if x.exists())
    header = fa.read_text().splitlines()[0]
    seq = "".join(l.strip() for l in fa.read_text().splitlines() if not l.startswith(">"))
    # the kinase-domain fasta header records the UniProt slice it was cut from
    start = 1
    for tok in header.replace("|", " ").split():
        if "-" in tok and tok.split("-")[0].isdigit():
            start = int(tok.split("-")[0])
            break
    return {"pocket": pocket, "sequence": seq, "seq_start": start,
            "symbol": pocket.get("gene") or slug, "uniprot": pocket.get("uniprot_id")}


# --------------------------------------------------------------------------
# MODALITY.  A small molecule and a peptide differ here in exactly one thing:
# which co-folding INPUT SLOT the binder occupies.  A small molecule is a SMILES
# string in `initial_smiles_list`; a biologic is a POLYMER CHAIN and belongs in
# `initial_protein_sequences` next to the target.  Everything downstream of the
# pose - consensus_signature(), weighted_jaccard(), coverage() - already works in
# target-side residues and never looks at the binder's chemistry, which is why
# the same scoring machinery covers both modalities (PROJECT_GOAL.md 4.4 gives
# InterfaceMatch a `modality` field for exactly this reason).
# --------------------------------------------------------------------------
SM = "small_molecule"
# the modalities data/biologic_drugs.csv actually emits are peptide / protein /
# antibody / other. "other" is that corpus's label for oligonucleotide, gene and
# cell products - NOT a polypeptide, so it is refused here rather than folded as
# a protein chain.
BIOLOGIC = ("peptide", "protein", "biologic", "antibody")
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

# measured on this target, from results/pipeline/<slug>/repurpose_state.json and
# results/cofold_constraint_evidence.json - see credit_note() for the numbers
MEASURED_SM_CREDITS = 4.55
MEASURED_CREDITS_PER_SECOND = 0.0500


def row_modality(row):
    """Modality of a shortlist row.

    An ABSENT field means small molecule, so every shortlist written before
    biologics existed keeps its meaning and the small-molecule path is untouched.
    """
    m = (row.get("modality") or SM).strip().lower().replace(" ", "_")
    return SM if m in ("", "sm", "smallmolecule", SM) else m


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y", "t")


def binder_chains(row):
    """The binder's chain sequences, or (None, why not).

    A biologic binder is a polymer: one chain for a peptide, several for a Fab
    (heavy + light).  Chains are separated by '/' - which is what
    data/biologic_drugs.csv writes - or by ';' or '+'.

    A residue outside the 20 standard amino acids is REFUSED rather than silently
    mangled.  Many approved peptides carry non-natural residues (Aib, D-amino
    acids, N-methylation, a lipid on a lysine); stjames can express those as
    ResidueModification(position, ccd) but this path does not build them, and
    dropping the modification would submit a different molecule from the drug.
    """
    raw = (row.get("sequence") or "").strip()
    if not raw:
        return None, "no sequence on the row"
    for sep in ("/", "+"):
        raw = raw.replace(sep, ";")
    chains = [c.strip().upper().replace(" ", "")
              for c in raw.split(";") if c.strip()]
    if not chains:
        return None, "no sequence on the row"
    for c in chains:
        bad = sorted(set(c) - AA20)
        if bad:
            return None, (f"non-standard residue(s) {''.join(bad)}: outside the 20 standard "
                          f"amino acids, so it needs a stjames ResidueModification with a CCD "
                          f"code, which this path does not build")
    return chains, None


def cofold_inputs(info, row):
    """Co-folding input slots for one shortlist row.

    Returns (protein_sequences, smiles_list, binder_slot), where binder_slot is
    the (input_type, input_index) that a pocket constraint names as the binder.

    THE TARGET IS ALWAYS PROTEIN INPUT 0.  Every constraint this file builds
    addresses target residues as Token(input_type="protein", input_index=0,
    token_index=auth_resid - seq_start), so the target must keep slot 0 whatever
    the binder is; the binder is APPENDED, never prepended.  input_index is an
    index into the list for that input_type, not a global input counter - the
    evidence for that reading is in stjames'
    BatchProteinCofoldingWorkflow.validate_model_compatibility, which range-checks
    a ligand-side pocket constraint against `ligand_slot_count`, the number of
    LIGAND slots alone.
    """
    mod = row_modality(row)
    name = row.get("name", row.get("struct_id", "?"))
    if mod == SM:
        if not row.get("smiles"):
            sys.exit(f"{name}: small-molecule row with no SMILES")
        return [info["sequence"]], [row["smiles"]], ("ligand", 0)
    if mod not in BIOLOGIC:
        sys.exit(f"{name}: unknown modality {mod!r} (expected {SM} or one of {BIOLOGIC})")
    if row.get("smiles"):
        sys.exit(f"{name}: modality {mod} but the row also carries SMILES. A binder "
                 f"occupies one slot, not two - a peptide sent as both a chain and a "
                 f"SMILES would be co-folded twice over")
    chains, why = binder_chains(row)
    if chains is None:
        sys.exit(f"{name}: {why}")
    if _truthy(row.get("cyclic")):
        # `cyclic` lives on stjames' ProteinSequence, and initial_protein_sequences is
        # typed list[ProteinSequence] | list[str] - a UNION OF HOMOGENEOUS LISTS, so
        # one cyclic binder forces the target into a ProteinSequence too.
        from stjames.protein import ProteinSequence
        seqs = ([ProteinSequence(sequence=info["sequence"])]
                + [ProteinSequence(sequence=c, cyclic=True) for c in chains])
    else:
        seqs = [info["sequence"], *chains]
    return seqs, [], ("protein", 1)


def pocket_tokens(info):
    """0-based token indices, within protein input 0, of the site residues.

    The pocket is stored in author numbering; a token index is the offset into the
    construct.  If the two disagree the constraint would silently point at the
    wrong residues, so a residue outside the construct is a hard stop.
    """
    seq_len, start = len(info["sequence"]), info["seq_start"]
    tok = [a - start for a in info["pocket"]["residue_ids"]]
    bad = [a for a, t in zip(info["pocket"]["residue_ids"], tok) if not 0 <= t < seq_len]
    if bad:
        sys.exit(f"pocket residues {bad} (author numbering) fall outside the {seq_len}-residue "
                 f"construct starting at {start}; the token indices a constraint would carry "
                 f"are wrong, so nothing was submitted")
    return tok


def check_indices(wf):
    """Range-check every constraint index against the inputs actually present.

    stjames does NOT do this for a single co-fold: ProteinCofoldingWorkflow's
    model_validator checks model/feature compatibility only, and the one range
    check in the library (ligand slots) lives in BatchProteinCofoldingWorkflow and
    does not run here.  An out-of-range index would therefore validate locally,
    submit, and burn credits on a constraint pointing at nothing.
    """
    seqs = [s if isinstance(s, str) else s.sequence for s in wf.initial_protein_sequences]
    n_lig = len(wf.initial_smiles_list)
    for pc in wf.pocket_constraints:
        lim = len(seqs) if pc.input_type == "protein" else n_lig
        if not 0 <= pc.input_index < lim:
            sys.exit(f"pocket constraint names {pc.input_type} input {pc.input_index}, but the "
                     f"payload has {lim} {pc.input_type} input(s)")
        for t in pc.contacts:
            if t.input_type != "protein":
                sys.exit(f"pocket-constraint contacts must be polymer residues, got {t.input_type}")
            if not 0 <= t.input_index < len(seqs):
                sys.exit(f"contact token names protein input {t.input_index}, payload has {len(seqs)}")
            if not 0 <= t.token_index < len(seqs[t.input_index]):
                sys.exit(f"contact token_index {t.token_index} is outside protein input "
                         f"{t.input_index} ({len(seqs[t.input_index])} residues)")
    return True


def cofold_payload(name, seqs, smiles, pcs, max_credits):
    """The exact body rowan.submit_protein_cofolding_workflow would POST.

    Mirrors rowan/workflows/protein_cofolding.py (the ProteinCofoldingWorkflow it
    constructs, then its `data` envelope), so --dry-run validates the real thing:
    the pydantic models validate ON CONSTRUCTION, offline, for zero credits.
    """
    from stjames.workflows.protein_cofolding import ProteinCofoldingWorkflow
    wf = ProteinCofoldingWorkflow(
        use_msa_server=True,
        use_potentials=False,
        contact_constraints=[],
        pocket_constraints=pcs or [],
        bond_constraints=[],
        templates=[],
        num_samples=None,
        model="boltz_2",
        ligand_binding_affinity_index=None,   # affinity is never read - PROJECT_GOAL.md 4.4
        initial_smiles_list=smiles or [],
        initial_protein_sequences=seqs,
        initial_dna_sequences=[],
        initial_rna_sequences=[],
        do_pose_refinement=False,
        compute_strain=False,
    )
    check_indices(wf)
    return wf, {"workflow_type": "protein_cofolding",
                "workflow_data": wf.model_dump(serialize_as_any=True, mode="json"),
                "name": name, "folder_uuid": None, "max_credits": max_credits,
                "webhook_url": None, "is_draft": False}


def estimate_credits(args, mod, seqs):
    """(credits, how it was arrived at) for one job.

    A small molecule gets the measured figure. A polymer binder gets the guard
    scaled by token count, because the one thing this project has measured about
    co-folding cost is that it does NOT depend on ligand size while the receptor is
    held fixed - which says the receptor dominates, and says nothing about what
    happens when the input grows by a whole chain. The exponent is an ASSUMPTION
    (--token-exponent, default 2.0 for attention-dominated cost), used only to set
    a budget cap; it is not a measurement and no number here should be reported as
    a cost.
    """
    if mod == SM:
        return args.per_job_credits, (f"guard {args.per_job_credits}; measured "
                                      f"4.55 cr (sd 1.02, n=42) for this receptor")
    lens = [len(x) if isinstance(x, str) else len(x.sequence) for x in seqs]
    n_target, n_binder = lens[0], sum(lens[1:])
    scale = ((n_target + n_binder) / n_target) ** args.token_exponent
    return (args.per_peptide_credits * scale,
            f"EXTRAPOLATED: {n_target} target + {n_binder} binder tokens, "
            f"x{scale:.2f} on an assumed tokens^{args.token_exponent} scaling")


def credit_note():
    """Where the per-job credit numbers come from, and what is extrapolation."""
    return (
        "MEASURED (small molecule, this target): 42 collected co-folds of the same "
        "329-residue receptor cost 4.55 credits each (sd 1.02, 95% CI 4.23-4.87, range "
        "3.06-6.96, n=42), and the cost is FLAT in ligand size across 2-85 heavy atoms "
        "(Pearson r=-0.090, p=0.573, n=42). Four separately timed probes billed at "
        "0.0500 credits/second each (n=4, no spread), so a credit is wall-clock time: "
        "4.55 credits is ~91 s.\n"
        "EXTRAPOLATED (peptide): a 9-50 residue binder adds 9-50 tokens to a "
        "329-token receptor (+3% to +15%), and ligand-size independence above says the "
        "receptor dominates, so ~4.6-5.5 credits is the expectation. The unquantified "
        "term is the SECOND MSA: use_msa_server=True now has two chains to search, and "
        "nothing here measures what that costs. Confidence: LOW on the exact number, "
        "moderate on the order of magnitude. One real peptide co-fold settles it, and "
        "max_credits caps the downside at 3x the guard."
    )


# --------------------------------------------------------------------------
# biologic corpus.  Column names are resolved through aliases because
# data/biologic_drugs.csv is built elsewhere; a missing column is reported with
# the header that WAS found rather than guessed around.
BIO_COLUMNS = {
    "struct_id": ("struct_id", "id", "drug_id", "drugbank_id", "chembl_id"),
    "name": ("name", "drug_name", "generic_name"),
    "sequence": ("sequence", "seq", "aa_sequence", "sequences", "chains"),
    "targets": ("targets", "target", "target_symbols", "gene_targets"),
}
BIO_OPTIONAL = {
    "modality": ("modality", "type", "drug_type", "class"),
    "cyclic": ("cyclic", "is_cyclic"),
    "n_targets": ("n_targets",),
    # the corpus's own per-row verdict on whether its sequence is co-foldable, and
    # why not. It is used rather than re-derived: that decision (which non-standard
    # residue was substituted, which was excluded) belongs to whoever built the
    # corpus and is auditable there.
    "usable": ("usable_sequence", "usable"),
    "exclusion_reason": ("exclusion_reason",),
    # the corpus's residue-by-residue check against a published reference sequence,
    # where it made one. A row can be perfectly co-foldable (usable_sequence=1) and
    # still carry a sequence that disagrees with the published drug, so this has to
    # be surfaced at shortlist time rather than left in the CSV.
    "reference_check": ("reference_check",),
}


def _resolve_columns(fieldnames):
    low = {(f or "").strip().lower(): f for f in (fieldnames or [])}
    got = {}
    for want, aliases in {**BIO_COLUMNS, **BIO_OPTIONAL}.items():
        for a in aliases:
            if a in low:
                got[want] = low[a]
                break
    return got, [k for k in BIO_COLUMNS if k not in got]


def net_charge(seq):
    """K+R minus D+E: the closest peptide analogue of cLogP for decoy matching.

    Ignores histidine and the termini, so it is a rank statistic for matching,
    not a titration curve.
    """
    return sum(seq.count(c) for c in "KR") - sum(seq.count(c) for c in "DE")


def _split_targets(v):
    return [t.strip() for t in str(v or "").replace(",", ";").split(";") if t.strip()]


def _shortlist_biologic(args, info, sym):
    """A peptide/protein shortlist with the SAME role labels as the small-molecule one.

    Positives are approved biologics annotated against the target.  Decoys are other
    APPROVED PEPTIDES, matched on (length, net charge) - the peptide analogue of the
    (MW, cLogP) match used for small molecules.  Length is the confound that matters
    for an interface score: a longer binder buries more target residues for reasons
    that have nothing to do with binding (PROJECT_GOAL.md 1.4a), so a decoy that is
    the same size as the positive it replaces takes that explanation away.  Net charge
    stands in for cLogP because peptide-epitope recognition is electrostatics-heavy
    and there is no meaningful logP for a 30-mer.

    --hard-decoys narrows the pool to peptides that are annotated against SOME other
    target - i.e. known binders of a different site, the peptide analogue of the
    "-tinib that does not hit this kinase" rule.  A peptide that binds nothing in the
    corpus proves less when it fails to score.

    One confound has no small-molecule counterpart and is handled explicitly: a decoy
    that is a SEQUENCE RELATIVE of a positive may well bind the same site, so anything
    above --max-identity to any positive is excluded, and the identity actually
    achieved is recorded per decoy.  (difflib ratio, a crude similarity proxy, not an
    alignment identity - it is used to EXCLUDE, so a loose measure errs safe.)
    """
    import difflib

    src = Path(args.biologics_csv) if args.biologics_csv else ROOT / "data" / "biologic_drugs.csv"
    if not src.exists():
        sys.exit(f"--biologics needs a biologic corpus and {src} does not exist.\n"
                 f"Expected a CSV with columns (aliases in brackets):\n" +
                 "\n".join(f"    {k:<10} [{', '.join(v)}]" for k, v in BIO_COLUMNS.items()) +
                 f"\n    optional   [{', '.join(sum(BIO_OPTIONAL.values(), ()))}]\n"
                 f"`sequence` holds the one-letter chain sequence; several chains "
                 f"(a Fab) are separated by ';'.")
    with src.open(newline="") as fh:
        rdr = csv.DictReader(fh)
        cols, missing = _resolve_columns(rdr.fieldnames or [])
        if missing:
            sys.exit(f"{src} is missing column(s) {missing}. Header found: {rdr.fieldnames}")
        raw = list(rdr)

    usable, dropped = [], []
    for r in raw:
        name = str(r.get(cols["name"], "")).strip()
        if "usable" in cols and str(r.get(cols["usable"], "")).strip() == "0":
            why = (str(r.get(cols["exclusion_reason"], "")).strip() if "exclusion_reason" in cols
                   else "") or "corpus marked the sequence unusable"
            dropped.append({"name": name, "why": f"corpus verdict: {why}"})
            continue
        mod_raw = (str(r.get(cols["modality"], "")).strip().lower() if "modality" in cols else "")
        if mod_raw == "other":
            dropped.append({"name": name, "why": "corpus modality 'other' (oligonucleotide, "
                                                 "gene or cell product): not a polypeptide, so "
                                                 "it cannot be a protein chain in a co-fold"})
            continue
        row = {"struct_id": str(r.get(cols["struct_id"], "")).strip(),
               "name": str(r.get(cols["name"], "")).strip(),
               "sequence": str(r.get(cols["sequence"], "")).strip(),
               "targets": _split_targets(r.get(cols["targets"])),
               "modality": mod_raw or "peptide", "source_modality": mod_raw,
               "cyclic": _truthy(r.get(cols["cyclic"])) if "cyclic" in cols else False,
               "n_targets": r.get(cols["n_targets"]) if "n_targets" in cols else None,
               "reference_check": (str(r.get(cols["reference_check"], "")).strip()
                                   if "reference_check" in cols else "")}
        chains, why = binder_chains(row)
        if chains is None:
            dropped.append({"name": row["name"] or row["struct_id"], "why": why})
            continue
        row["chains"] = chains
        row["n_chains"] = len(chains)
        row["length"] = sum(len(c) for c in chains)
        row["net_charge"] = net_charge("".join(chains))
        usable.append(row)

    known = [r for r in usable if sym in r["targets"]]
    print(f"{src.name}: {len(raw)} rows, {len(usable)} with a co-foldable sequence, "
          f"{len(dropped)} not representable")
    for d in dropped[:6]:
        print(f"   not evaluated  {d['name'][:32]:<32} {d['why'][:88]}")
    if len(dropped) > 6:
        print(f"   ... and {len(dropped) - 6} more (every one of them recorded in the output)")
    if not known:
        seen = sorted({t for r in usable for t in r["targets"]})
        sys.exit(f"\nno biologic in {src.name} is annotated against {sym} "
                 f"({len(usable)} rows scanned, {len(seen)} distinct target symbols present). "
                 f"That is a result, not a bug: without a positive control there is nothing "
                 f"to validate a ranking against, so this target is a 'not evaluated' row. "
                 f"First few symbols present: {seen[:10]}")

    def _nt(r):
        try:
            return float(r["n_targets"])
        except (TypeError, ValueError):
            return 0.0

    known.sort(key=lambda r: (-_nt(r), r["name"]))
    positives = known[:args.n_positives]
    pos_seq = ["".join(r["chains"]) for r in positives]

    pool = []
    for r in usable:
        if sym in r["targets"]:
            continue
        if args.hard_decoys and not r["targets"]:
            continue
        ident = max((difflib.SequenceMatcher(None, "".join(r["chains"]), p).ratio()
                     for p in pos_seq), default=0.0)
        if ident >= args.max_identity:
            dropped.append({"name": r["name"], "why": f"similarity {ident:.2f} to a known binder "
                                                      f">= --max-identity {args.max_identity}"})
            continue
        pool.append((r, ident))

    chosen, used = [], set()
    for p in positives:
        scored = sorted(pool, key=lambda ri: abs(ri[0]["length"] - p["length"]) / 10
                                             + abs(ri[0]["net_charge"] - p["net_charge"]))
        for r, ident in scored:
            if r["struct_id"] not in used:
                used.add(r["struct_id"])
                chosen.append((r, ident, p["name"]))
                break
        if len(chosen) >= args.n_decoys:
            break
    # Top up by DISTANCE to the nearest positive, not in file order and not at random
    # (the small-molecule path tops up with an RNG draw). File order here put 664- and
    # 1334-residue antibodies against 31-44 residue peptide positives - the exact size
    # confound the matching exists to remove, and 10-40x the co-folding cost each.
    rest = sorted(((r, ident) for r, ident in pool if r["struct_id"] not in used),
                  key=lambda ri: min(abs(ri[0]["length"] - p["length"]) / 10
                                     + abs(ri[0]["net_charge"] - p["net_charge"])
                                     for p in positives))
    for r, ident in rest:
        if len(chosen) >= args.n_decoys:
            break
        used.add(r["struct_id"])
        chosen.append((r, ident, "nearest positive overall"))

    def _row(r, role, ident=None, matched=None):
        d = {"struct_id": r["struct_id"], "name": r["name"], "role": role,
             "modality": r["modality"] if r["modality"] in BIOLOGIC else "protein",
             "source_modality": r["source_modality"],
             "sequence": ";".join(r["chains"]), "n_chains": r["n_chains"],
             "length": r["length"], "net_charge": r["net_charge"],
             "annotated_targets": r["targets"]}
        if r["cyclic"]:
            d["cyclic"] = True
        if ident is not None:
            d["identity_to_nearest_positive"] = round(ident, 3)
        if matched:
            d["matched_to"] = matched
        return d

    role = "hard_decoy" if args.hard_decoys else "decoy"
    shortlist = ([_row(r, "known_binder") for r in positives] +
                 [_row(r, role, ident, matched) for r, ident, matched in chosen])
    print(f"\n{sym}: {len(known)} approved biologics annotated against it in {src.name}")
    print(f"shortlist: {len(positives)} known binders + {len(chosen)} {role}s "
          f"= {len(shortlist)} co-folds")
    if chosen:
        print(f"decoy sequence similarity to the nearest positive: "
              f"max {max(i for _, i, _ in chosen):.2f}, "
              f"median {float(np.median([i for _, i, _ in chosen])):.2f} "
              f"(difflib ratio, n={len(chosen)})")
    # A row can be co-foldable and still hold a sequence that disagrees with the
    # published drug. Co-folding it would answer a question about a molecule that
    # is not the approved one, so say so here rather than let it pass silently.
    flagged = [r for r in ([p for p in positives] + [c for c, _, _ in chosen])
               if str(r.get("reference_check", "")).upper().startswith("MISMATCH")]
    for r in flagged:
        print(f"   WARNING  {r['name'][:32]:<32} the corpus sequence DISAGREES with the "
              f"published one ({r['reference_check']}). Co-folding this row models a "
              f"different molecule from the approved drug; fix the row before trusting "
              f"any result that depends on it.")

    extra = {
        "sequence_reference_mismatches": [
            {"name": r["name"], "reference_check": r["reference_check"]} for r in flagged],
        "n_known_in_corpus": len(known),
        "biologic_source": str(src.relative_to(ROOT)) if src.is_relative_to(ROOT) else str(src),
        "decoy_matching": (f"other approved {'target-annotated ' if args.hard_decoys else ''}"
                           f"peptides, nearest neighbour in (length/10, net charge K+R-D-E) to "
                           f"each positive, excluding anything annotated against the target and "
                           f"anything with difflib similarity >= {args.max_identity} to a positive"),
        "biologic_not_evaluated": dropped,
    }
    return shortlist, extra


def _finish_shortlist(args, sym, shortlist, extra):
    st = state_load(args.pipeline)
    if (st.get("jobs") and not args.append and not args.out
            and not args.replace_shortlist):
        sys.exit(f"the state for {args.pipeline!r} already tracks {len(st['jobs'])} job(s) "
                 f"against {st.get('target_symbol')}, and a plain `shortlist` REPLACES the "
                 f"shortlist, orphaning every one of them (this is what "
                 f"test_shortlist_append_is_idempotent exists to catch).\n"
                 f"  --append            add these rows to the existing shortlist\n"
                 f"  --out <path>        write the shortlist to a file and leave the state alone\n"
                 f"  --replace-shortlist do it anyway")
    if args.append and st.get("shortlist"):
        seen = {d["struct_id"] for d in st["shortlist"]}
        shortlist = st["shortlist"] + [d for d in shortlist if d["struct_id"] not in seen]
    payload = {"target_symbol": sym, "shortlist": shortlist, **extra}
    if args.out:
        Path(args.out).write_text(json.dumps(payload, indent=1))
        print(f"wrote {args.out}; the pipeline state was NOT touched "
              f"(--out exists so a dry run cannot mutate a committed artifact)")
    else:
        st.update(payload)
        state_save(args.pipeline, st)
    for s in shortlist[:8]:
        print(f"   {s['role']:<12} {s['name']}")


# --------------------------------------------------------------------------
def cmd_shortlist(args):
    """Positive controls (known binders) + property-matched decoys.

    Both are needed: the known binders say whether the ranking works at all, and the
    decoys give it something to be enriched against. Decoys are matched on MW and
    cLogP so that "looks like a drug that binds kinases" is not what is being measured.

    --biologics builds the same thing out of peptide/protein drugs instead; see
    _shortlist_biologic for what "property-matched" has to mean when the binder is
    a polymer.
    """
    info = target_info(args.pipeline)
    sym = args.target or info["symbol"]
    if args.biologics:
        shortlist, extra = _shortlist_biologic(args, info, sym)
        _finish_shortlist(args, sym, shortlist, extra)
        return

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors
    RDLogger.DisableLog("rdApp.*")

    rows = [r for r in csv.DictReader((ROOT / "data" / "approved_drugs.csv").open())
            if r["smiles"]]
    known = [r for r in rows if sym in (r.get("targets") or "").split(";")]
    if not known:
        sys.exit(f"no approved drug in the corpus is annotated against {sym}")

    def props(smi):
        m = Chem.MolFromSmiles(smi)
        return (Descriptors.MolWt(m), Crippen.MolLogP(m)) if m else None

    kp = [(r, props(r["smiles"])) for r in known]
    kp = [(r, p) for r, p in kp if p]
    kp.sort(key=lambda x: -float(x[0].get("n_targets") or 0))
    positives = [r for r, _ in kp[:args.n_positives]]
    pos_props = [p for _, p in kp[:args.n_positives]]

    # decoys.
    # --hard-decoys picks APPROVED KINASE INHIBITORS that are not annotated against
    # this target (INN convention: -tinib). This is a much harder null than MW/cLogP
    # matching, which happily returns steroids and perfluorocarbons - chemically
    # matched on bulk properties but not plausible ATP-site binders, so separating
    # them proves little. A -tinib decoy binds SOME kinase ATP pocket by construction.
    # Caveat worth stating in any result: DrugCentral annotation is incomplete, so
    # "not annotated against KDR" is not "does not bind KDR" - several of these
    # (crizotinib, bosutinib) are promiscuous. A hard decoy scoring well may be a
    # labelling gap rather than a false positive.
    pool = [(r, props(r["smiles"])) for r in rows
            if sym not in (r.get("targets") or "").split(";")
            and (not args.hard_decoys or r["name"].strip().endswith("tinib"))]
    pool = [(r, p) for r, p in pool if p]
    rng = np.random.default_rng(SEED)
    chosen, used = [], set()
    for mw, lp in pos_props:
        scored = sorted(pool, key=lambda rp: abs(rp[1][0] - mw) / 100 + abs(rp[1][1] - lp))
        for r, _ in scored:
            if r["struct_id"] not in used:
                used.add(r["struct_id"])
                chosen.append(r)
                break
        if len(chosen) >= args.n_decoys:
            break
    while len(chosen) < args.n_decoys and pool:
        r, _ = pool[int(rng.integers(len(pool)))]
        if r["struct_id"] not in used:
            used.add(r["struct_id"])
            chosen.append(r)

    shortlist = ([{"struct_id": r["struct_id"], "name": r["name"], "smiles": r["smiles"],
                   "role": "known_binder"} for r in positives] +
                 [{"struct_id": r["struct_id"], "name": r["name"], "smiles": r["smiles"],
                   "role": "hard_decoy" if args.hard_decoys else "decoy"}
                  for r in chosen])
    print(f"{sym}: {len(known)} approved drugs annotated against it in the corpus")
    print(f"shortlist: {len(positives)} known binders + {len(chosen)} property-matched decoys "
          f"= {len(shortlist)} co-folds")
    _finish_shortlist(args, sym, shortlist,
                      {"n_known_in_corpus": len(known),
                       "decoy_matching": "nearest neighbour in (MW/100, cLogP) to each positive, "
                                         "excluding anything annotated against the target"})


# --------------------------------------------------------------------------
def cmd_submit(args):
    """One co-folding job per shortlist row, small molecule or biologic.

    The modality decides the input slot and NOTHING else: same model, same MSA
    setting, same absence of an affinity index, same scoring downstream.

        small molecule   initial_protein_sequences=[target]
                         initial_smiles_list=[drug]
        biologic         initial_protein_sequences=[target, binder, ...]
                         initial_smiles_list=[]

    --dry-run builds and VALIDATES the payload locally (the stjames pydantic models
    validate on construction) and prints it, without a network call, an API key, or
    a credit.
    """
    from stjames.workflows.protein_cofolding import PocketConstraint, Token

    info = target_info(args.pipeline)
    st = state_load(args.pipeline)
    # What protein this pipeline directory actually holds. pocket.json carries a gene
    # symbol only when the site came from a gene-annotated target, so the pipeline's
    # own recorded symbol is the fallback; target_info() substitutes the SLUG when
    # neither exists, which is not an identity worth comparing against.
    pipeline_symbol = info["pocket"].get("gene") or st.get("target_symbol")
    if args.shortlist:
        if not args.dry_run:
            sys.exit("--shortlist reads a shortlist from a file instead of the pipeline "
                     "state and is a --dry-run input only; a real submission must come "
                     "from the state it will write its job ids back into")
        loaded = json.loads(Path(args.shortlist).read_text())
        st = dict(st)
        st["shortlist"] = loaded["shortlist"] if isinstance(loaded, dict) else loaded
        if isinstance(loaded, dict) and loaded.get("target_symbol"):
            st["target_symbol"] = loaded["target_symbol"]
    if not st.get("shortlist"):
        sys.exit("run `shortlist` first")
    st.setdefault("target_symbol", info["symbol"])
    # The shortlist names a target symbol; the pipeline directory holds ONE construct.
    # If they disagree the run would co-fold this shortlist's drugs against a different
    # protein and never say so - the failure mode that already cost this project a
    # scoring run against the wrong designs.
    mismatch = bool(pipeline_symbol) and st["target_symbol"] != pipeline_symbol
    if mismatch and not args.allow_target_mismatch:
        sys.exit(f"the shortlist is for {st['target_symbol']} but pipeline "
                 f"{args.pipeline!r} holds the {pipeline_symbol} construct "
                 f"({len(info['sequence'])} aa, {info['uniprot']}). Co-folding these "
                 f"binders against that protein would be meaningless. Build the "
                 f"{st['target_symbol']} pipeline, or pass --allow-target-mismatch if you "
                 f"are exercising the payload builder and know the pairing is nonsense.")
    if not pipeline_symbol:
        print(f"note: pipeline {args.pipeline!r} records no gene symbol for its construct "
              f"({len(info['sequence'])} aa, {info['uniprot']}), so the shortlist's target "
              f"{st['target_symbol']} could not be checked against it.")
    # a dry run from a file is a rehearsal of a FRESH submission, so it is not filtered
    # against jobs the state has already submitted for some other shortlist
    jobs = {} if args.shortlist else st.setdefault("jobs", {})
    tok = pocket_tokens(info)

    if not args.dry_run:
        import rowan
        load_env()
        rowan.api_key = os.environ.get("ROWAN_API_KEY") or sys.exit("no ROWAN_API_KEY in .env")

    print(f"per-job credit budget: {args.per_job_credits} small molecule, "
          f"{args.per_peptide_credits} biologic\n{credit_note()}\n")

    spent, submitted = 0.0, 0
    for d in st["shortlist"]:
        if d["struct_id"] in jobs:
            continue
        mod = row_modality(d)
        seqs, smiles, binder_slot = cofold_inputs(info, d)
        per_job, per_job_why = estimate_credits(args, mod, seqs)
        if spent + per_job > args.max_credits:
            print(f"budget guard: stopping at {submitted} {'payloads' if args.dry_run else 'submissions'} "
                  f"(~{spent:.0f} of {args.max_credits} credits)")
            break
        # UNCONSTRAINED BY DEFAULT, and that is a measured decision that departs from
        # PROJECT_GOAL.md F3. F3 argues for a pocket constraint so a promiscuous drug
        # does not dock somewhere irrelevant. Measured on this target, the constraint
        # does the opposite damage - it forces EVERY ligand into the site, so a
        # non-binder scores like a binder and the ranking says nothing:
        #
        #                      constrained            unconstrained
        #   axitinib (binder)  16/17 pocket, 13.1 cr  16/17 pocket, 5.1 cr
        #   aspirin  (control) 11/17 pocket, 9-11 cr   6/17 pocket, 5.7 cr
        #
        # Unconstrained keeps the true binder in the pocket, halves the false signal
        # from the control, doubles the separation, and costs 2.5x less. Pass
        # --constrain-pocket to restore F3's behaviour and see this for yourself.
        # That measurement is SMALL-MOLECULE ONLY (n=2 drugs, one target); nothing
        # here says which way it goes for a peptide, whose binding surface is the
        # size of the site itself.
        pcs = None
        if args.constrain_pocket:
            if len(seqs) > 2:
                sys.exit(f"{d['name']}: --constrain-pocket with a {len(seqs) - 1}-chain binder is "
                         f"refused. A PocketConstraint names ONE binder input, so a Fab would "
                         f"need a per-chain decision (which chain sits on the epitope) that "
                         f"nothing here has measured.")
            pcs = [PocketConstraint(input_type=binder_slot[0], input_index=binder_slot[1],
                                    contacts=[Token(input_type="protein", input_index=0,
                                                    token_index=t) for t in tok],
                                    max_distance=args.pocket_distance, force=False)]
        name = f"{st['target_symbol']} cofold - {d['name'][:40]}"
        max_credits = int(per_job * 3)
        wf, payload = cofold_payload(name, seqs, smiles, pcs, max_credits)

        if args.dry_run:
            chains = [s if isinstance(s, str) else s.sequence for s in seqs]
            print(f"  DRY RUN  {d['role']:<12} {d['name'][:38]:<38} modality={mod}")
            if mismatch:
                print(f"      !! TARGET MISMATCH: the shortlist is {st['target_symbol']}, "
                      f"protein input 0 is {pipeline_symbol} ({info['uniprot']}). "
                      f"Payload shape only; "
                      f"this pairing is not a biological claim")
            for i, c in enumerate(chains):
                print(f"      protein input {i}: {len(c):>4} aa  "
                      f"{'TARGET ' + str(pipeline_symbol or info['uniprot']) if i == 0 else 'BINDER'}")
            for i, s in enumerate(smiles):
                print(f"      ligand  input {i}: {s[:60]}")
            print(f"      constraints: "
                  f"{'none (unconstrained)' if not pcs else f'pocket on {binder_slot[0]} input {binder_slot[1]}, {len(tok)} target contact tokens, <= {args.pocket_distance} A'}")
            print(f"      est. credits {per_job:.1f} ({per_job_why}); "
                  f"max_credits={max_credits}; payload validated by stjames, nothing sent")
            print(indent_json(payload))
        else:
            wfr = rowan.submit_protein_cofolding_workflow(
                initial_protein_sequences=seqs,
                initial_smiles_list=smiles or None,
                pocket_constraints=pcs, model="boltz_2", use_msa_server=True,
                max_credits=max_credits, name=name)
            jobs[d["struct_id"]] = {"uuid": str(wfr.uuid), "name": d["name"],
                                    "role": d["role"], "modality": mod}
            print(f"  submitted {d['role']:<12} {d['name'][:38]}")
        spent += per_job
        submitted += 1
    if args.dry_run:
        print(f"\n{submitted} payload(s) built and validated locally. "
              f"0 credits spent; ~{spent:.1f} credits is what submitting them would cost "
              f"on the estimates above.")
        return
    state_save(args.pipeline, st)
    print(f"\n{submitted} submitted, {len(jobs)} total tracked")


def indent_json(payload):
    return "\n".join("      " + l for l in json.dumps(payload, indent=1).splitlines())


# --------------------------------------------------------------------------
def cmd_collect(args):
    import rowan
    load_env()
    rowan.api_key = os.environ.get("ROWAN_API_KEY") or sys.exit("no ROWAN_API_KEY in .env")
    st = state_load(args.pipeline)
    jobs = st.get("jobs") or sys.exit("no jobs; run `submit` first")
    out = paths(args.pipeline)["poses"]
    out.mkdir(parents=True, exist_ok=True)

    pending = [k for k, v in jobs.items() if not v.get("pose_path")]
    for attempt in range(args.max_polls):
        still = []
        for sid in pending:
            j = jobs[sid]
            wf = rowan.retrieve_workflow(j["uuid"])
            if not wf.is_finished():
                still.append(sid)
                continue
            j["status"] = str(wf.status)
            j["credits"] = wf.credits_charged
            if "OK" not in str(wf.status).upper():
                j["pose_path"] = None
                print(f"  FAILED {j['name']}: {wf.status}")
                continue
            res = (wf.data or {}).get("result") or wf.data
            uid = None
            for key in ("predicted_structure_uuid", "predicted_refined_structure_uuid"):
                uid = (res or {}).get(key) or uid
            if not uid:
                j["pose_path"] = None
                j["why"] = f"no structure uuid in result (keys {list((res or {}).keys())[:8]})"
                print(f"  NO POSE {j['name']}: {j['why']}")
                continue
            d = out / sid
            if not list(d.glob("*.pdb")):
                rowan.retrieve_protein(uid).download_pdb_file(str(d))
            files = list(d.glob("*.pdb"))
            j["pose_path"] = str(files[0].relative_to(ROOT)) if files else None
            j["scores"] = (res or {}).get("scores")
            j["constrained"] = bool((res or {}).get("pocket_constraints"))
            print(f"  got {j['role']:<12} {j['name'][:36]}  credits={wf.credits_charged}")
        pending = still
        state_save(args.pipeline, st)
        if not pending:
            break
        print(f"  {len(pending)} still running; waiting...", flush=True)
        time.sleep(args.poll_seconds)

    done = sum(1 for v in jobs.values() if v.get("pose_path"))
    total_cr = sum(v.get("credits") or 0 for v in jobs.values())
    print(f"\n{done}/{len(jobs)} poses collected, {total_cr:.1f} credits charged")


# --------------------------------------------------------------------------
def _versions():
    """Backend and library versions, so a board says what produced it.

    The audit found no rowan, stjames, Boltz model id or rdkit version recorded
    anywhere in the outputs, which makes a result impossible to attribute to a
    toolchain after the fact.
    """
    import platform
    out = {"python": platform.python_version(),
           "platform": f"{platform.system()} {platform.machine()}"}
    from importlib.metadata import PackageNotFoundError, version
    # distribution name -> import name, because several differ and rowan/stjames
    # expose no __version__ attribute at all
    for dist in ("rowan-python", "stjames", "rdkit", "biopython", "numpy", "scipy"):
        try:
            out[dist] = version(dist)
        except PackageNotFoundError:
            out[dist] = "not installed"
    out["p2rank"] = "2.5"
    return out


def cmd_score(args):
    """Rank the co-folded drugs by how much of the design signature they engage."""
    st = state_load(args.pipeline)
    jobs = st.get("jobs") or sys.exit("no jobs; run submit/collect first")
    info = target_info(args.pipeline)
    start = info["seq_start"]

    # THREE signatures, not one. m2_gate.py measured that the BoltzGen consensus is
    # WORSE than the P2Rank pocket at locating a real ligand's contacts, so ranking
    # drugs only against the design signature would hide the cheaper, better option.
    # PROJECT_GOAL.md 4.3 makes `source` a discriminator precisely so these are swappable.
    # The BoltzGen signature is OPTIONAL. It used to be loaded unconditionally with
    # a default tag of "v4", which is KDR's design run - so scoring any other target
    # silently ranked it against KDR's designs. A signature is only offered when a
    # design run exists FOR THIS PIPELINE (or is named explicitly), and asking to
    # rank by one that does not exist is an error rather than a silent substitution.
    sig = None
    sig_file = None
    if args.designs:
        sig_file = ROOT / "results" / f"boltzgen_signature_{args.designs}.json"
        if not sig_file.exists():
            sys.exit(f"--designs {args.designs!r} but {sig_file} does not exist.\n"
                     f"Run a design run for this target, or drop --designs to use the "
                     f"pocket signature (which measured at least as good).")
        sig = json.loads(sig_file.read_text())

    signatures = {}
    if sig is not None:
        design_sets = [set(d["contacts_auth"]) for d in sig["per_design"]]
        freq = Counter(r for s in design_sets for r in s)
        n_des = len(design_sets)
        signatures["boltzgen_consensus"] = {
            "core": {r for r, k in freq.items() if k / n_des >= args.core_threshold},
            "weights": {r: k / n_des for r, k in freq.items()},
            "provenance": f"{n_des} BoltzGen designs from {sig_file.name}, ipTM "
                          f"{min(d['iptm'] or 0 for d in sig['per_design']):.3f}-"
                          f"{max(d['iptm'] or 0 for d in sig['per_design']):.3f} "
                          f"(all below the 0.85 the plan suggests filtering at)",
            "workflow_uuid": sig.get("workflow_uuid"),
        }
    signatures.update({
        "p2rank_geometry": {
            "core": set(info["pocket"]["residue_ids"]),
            "weights": {r: 1.0 for r in info["pocket"]["residue_ids"]},
            "provenance": ("P2Rank 2.5 rank-1 pocket on the ligand-stripped "
                           "structure. Timing, measured: 2.1-2.5 s wall for a single "
                           "structure in isolation (P2Rank self-reports 1.87 s); "
                           "0.47 s/structure amortised over a 1,531-structure batch at "
                           "12 threads, 250 per JVM. The batch figure is not a "
                           "single-run cost."),
            "site_id": info["pocket"].get("site_id"),
            "structure_id": info["pocket"].get("structure_id"),
        },
    })
    klc = paths(args.pipeline)["target"] / "known_ligand_contacts.json"
    if klc.exists():
        k = json.loads(klc.read_text())
        ids = k.get("residue_ids") or k.get("known_ligand_residue_ids") or []
        if ids:
            signatures["known_ligand"] = {
                "core": {int(x) for x in ids},
                "weights": {int(x): 1.0 for x in ids},
                "provenance": "contacts of the ligand in the primary holo structure; "
                              "VALIDATION REFERENCE ONLY - it is the answer for one ligand",
            }
    for name, sg in signatures.items():
        print(f"signature {name:<20} {len(sg['core']):3d} residues  ({sg['provenance']})")
    print()
    if args.rank_by not in signatures:
        sys.exit(f"--rank-by {args.rank_by} is not available for this run "
                 f"(have: {', '.join(sorted(signatures))}). "
                 f"A BoltzGen signature needs --designs <tag> pointing at a design run "
                 f"for THIS target.")
    core = signatures[args.rank_by]["core"]
    weights = signatures[args.rank_by]["weights"]

    rows = []
    for sid, j in jobs.items():
        if not j.get("pose_path"):
            rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                         "status": "no_pose", "why": j.get("why") or j.get("status")})
            continue
        path = ROOT / j["pose_path"]
        try:
            struct = interfaces.load_structure(path, sid)
            ligs = interfaces.drug_like_ligands(struct)
            if not ligs:
                rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                             "status": "no_ligand_in_pose"})
                continue
            lg = ligs[0]
            con = interfaces.ligand_contacts(struct, lg["resname"], chain=lg.get("chain_id"),
                                             cutoff=CUTOFF)
            engaged = {int(r) + start - 1 for r in con["residue_ids"]}
        except Exception as exc:                        # noqa: BLE001
            rows.append({"struct_id": sid, "name": j["name"], "role": j["role"],
                         "status": "extract_failed", "why": repr(exc)[:120]})
            continue
        row = {"struct_id": sid, "name": j["name"], "role": j["role"],
               "status": "scored", "n_engaged": len(engaged),
               "engaged_residues": sorted(engaged), "by_signature": {}}
        for sname, sg in signatures.items():
            c, w = sg["core"], sg["weights"]
            hit = engaged & c
            num = sum(w.get(r, 0) for r in engaged & set(w))
            den = sum(w.values()) + len(engaged - set(w))
            # FIVE scores, because a bigger ligand touches more residues and so
            # scores higher for reasons that have nothing to do with binding -
            # PROJECT_GOAL.md 1.4a. E4 asks for at least two normalisations
            # compared and the default justified. Measured on this screen
            # (correlation of the score with the number of residues engaged,
            # and where a 39-residue lipopeptide decoy lands out of 29):
            #   core_coverage        r=+0.641  decoy rank  1
            #   weighted_jaccard     r=+0.520  decoy rank  2
            #   f1_engaged_core      r=+0.413  decoy rank  6
            #   jaccard_engaged_core r=+0.409  decoy rank  6
            #   precision_in_core    r=+0.135  decoy rank 14   <- default
            # Enrichment and median known-binder rank are IDENTICAL under all
            # five, so the separation is not an artefact of this choice; the
            # metric only decides where the oversized decoy lands.
            row["by_signature"][sname] = {
                "core_coverage": len(hit) / len(c) if c else float("nan"),
                "weighted_jaccard": num / den if den else 0.0,
                "f1_engaged_core": (2 * len(hit) / (len(engaged) + len(c))
                                    if (engaged or c) else 0.0),
                "jaccard_engaged_core": (len(hit) / len(engaged | c)
                                         if (engaged | c) else 0.0),
                "precision_in_core": len(hit) / len(engaged) if engaged else 0.0,
                "engaged_core": sorted(hit), "missed_core": sorted(c - engaged)}
        # convenience aliases for the signature this run is ranked by, so a reader
        # of a single row does not have to know which signature was primary
        row["core_coverage"] = row["by_signature"][args.rank_by]["core_coverage"]
        row["weighted_jaccard"] = row["by_signature"][args.rank_by]["weighted_jaccard"]
        rows.append(row)

    scored = [r for r in rows if r["status"] == "scored"]
    scored.sort(key=lambda r: -r["by_signature"][args.rank_by][args.metric])
    for i, r in enumerate(scored, 1):
        r["rank"] = i
    print(f"{len(scored)}/{len(rows)} drugs scored "
          f"({sum(1 for r in rows if r['status'] != 'scored')} without a usable pose)\n")
    print(f"{'rank':>4}  {'drug':<34}{'role':<13}"
          f"{args.metric[:9]:>9}{'core_cov':>10}{'engaged':>8}")
    print("-" * 80)
    for r in scored[:30]:
        b = r["by_signature"][args.rank_by]
        print(f"{r['rank']:>4}  {r['name'][:34]:<34}{r['role']:<13}"
              f"{b[args.metric]:9.3f}{b['core_coverage']:10.3f}{r['n_engaged']:8d}")

    # the same validation for EVERY signature, so the cheap arm is not hidden
    print(f"\n{'signature':<22}{'enrichment@25%':>16}{'median rank of a known binder':>32}")
    print("-" * 72)
    per_sig = {}
    for sname in signatures:
        order = sorted(scored, key=lambda r: -r["by_signature"][sname][args.metric])
        rk = [i for i, r in enumerate(order, 1) if r["role"] == "known_binder"]
        npos = len(rk)
        kk = max(1, round(0.25 * len(order)))
        hh = sum(1 for r in order[:kk] if r["role"] == "known_binder")
        bb = npos / len(order) if order else float("nan")
        e = (hh / kk) / bb if bb else float("nan")
        per_sig[sname] = {"enrichment_at_top_quartile": e, "known_binder_ranks": rk,
                          "median_rank": float(np.median(rk)) if rk else None}
        print(f"{sname:<22}{e:16.2f}{(np.median(rk) if rk else float('nan')):32.1f}")

    # validation: do the known binders come out on top?
    n_pos = sum(1 for r in scored if r["role"] == "known_binder")
    k = max(1, round(0.25 * len(scored)))
    hits = sum(1 for r in scored[:k] if r["role"] == "known_binder")
    base = n_pos / len(scored) if scored else float("nan")
    ef = (hits / k) / base if base else float("nan")
    ranks = [r["rank"] for r in scored if r["role"] == "known_binder"]
    print(f"\nVALIDATION: {n_pos} known binders, {len(scored) - n_pos} decoys")
    print(f"  known-binder ranks: {ranks}")
    print(f"  top {k}: {hits}/{k} are known binders (base rate {base:.3f}) "
          f"-> enrichment {ef:.2f}x")
    if ranks:
        print(f"  median rank of a known binder: {int(np.median(ranks))} of {len(scored)}")

    constrained = sorted({bool(j.get("constrained")) for j in jobs.values()})
    out = {
        "pipeline": "binding site -> interface signature -> approved drugs co-folded "
                    "WITH the target (unconstrained by default) -> pose scored against "
                    "the site -> ranked by interface overlap",
        "cofolding": {
            "pocket_constrained": (constrained[0] if len(constrained) == 1 else "mixed"),
            "what_that_means": (
                "unconstrained: the drug was folded with the target WITHOUT pocket "
                "conditioning, and the site was used afterwards to score the pose. "
                "This is post-hoc scoring against a pocket, not directed docking."
                if constrained == [False] else
                "constrained: the fold was conditioned on the pocket, which forces "
                "every ligand into the site and was measured to destroy discrimination"),
            "backend": "rowan submit_protein_cofolding_workflow",
            "model": "boltz_2",
        },
        "versions": _versions(),
        "ranking_metric": args.metric,
        "ranking_metric_choice": ("size-normalised; PROJECT_GOAL.md 1.4a warns raw "
                                  "overlap ranks bigger ligands higher for non-binding "
                                  "reasons. All five scores are stored per drug so the "
                                  "choice can be checked, and enrichment is identical "
                                  "under all of them on this screen."),
        "affinity_used": False,
        "affinity_note": "Boltz-2 can emit an affinity score; it is deliberately not read. "
                         "PROJECT_GOAL.md 4.4 removes affinity from the ranking path.",
        "ranked_by": args.rank_by,
        "signatures": {k: {"n_core": len(v["core"]), "provenance": v["provenance"],
                           "core_residues_auth": sorted(v["core"])}
                       for k, v in signatures.items()},
        "validation_by_signature": per_sig,
        "target": {"symbol": st.get("target_symbol"), "uniprot": info["uniprot"],
                   "structure": info["pocket"].get("structure_id")},
        "n_scored": len(scored), "n_known_binders": n_pos,
        "enrichment_at_top_quartile": ef, "known_binder_ranks": ranks,
        "results": rows,
    }
    dest = Path(args.out) if args.out else ROOT / "results" / f"repurpose_{args.pipeline}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pipeline", default="colorectal-cancer")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("shortlist"); s.set_defaults(fn=cmd_shortlist)
    s.add_argument("--target"); s.add_argument("--n-positives", type=int, default=10)
    s.add_argument("--n-decoys", type=int, default=20)
    s.add_argument("--hard-decoys", action="store_true",
                   help="use approved kinase inhibitors (-tinib) that do not hit the "
                        "target, instead of MW/cLogP-matched drugs")
    s.add_argument("--append", action="store_true",
                   help="add to the existing shortlist rather than replacing it")
    s.add_argument("--biologics", action="store_true",
                   help="build the shortlist from peptide/protein drugs "
                        "(data/biologic_drugs.csv) instead of small molecules")
    s.add_argument("--biologics-csv", default=None,
                   help="path to the biologic corpus; default data/biologic_drugs.csv")
    s.add_argument("--max-identity", type=float, default=0.5,
                   help="biologics only: drop a candidate decoy whose sequence similarity "
                        "(difflib ratio) to any positive reaches this. Two peptides with "
                        "related sequences may bind the same epitope, so a close relative "
                        "of a known binder is not a decoy")
    s.add_argument("--out", default=None,
                   help="write the shortlist here INSTEAD of into the pipeline state, so a "
                        "dry run cannot mutate a committed artifact. Feed it to "
                        "`submit --dry-run --shortlist <path>`")
    s.add_argument("--replace-shortlist", action="store_true",
                   help="allow replacing a shortlist whose jobs have already been submitted")

    s = sub.add_parser("submit"); s.set_defaults(fn=cmd_submit)
    s.add_argument("--max-credits", type=float, default=120)
    s.add_argument("--per-job-credits", type=float, default=6.0,
                   help="budget guard for a small-molecule co-fold; measured 4.55 credits "
                        "(sd 1.02, n=42) on this target")
    s.add_argument("--per-peptide-credits", type=float, default=6.0,
                   help="budget guard for a biologic co-fold. EXTRAPOLATED, not measured - "
                        "see credit_note(); nothing with a peptide binder has been run")
    s.add_argument("--dry-run", action="store_true",
                   help="build and validate each payload locally against the stjames "
                        "models and print it. No network call, no API key, no credits")
    s.add_argument("--shortlist", default=None,
                   help="read the shortlist from this file instead of the pipeline state "
                        "(requires --dry-run)")
    s.add_argument("--allow-target-mismatch", action="store_true",
                   help="submit a shortlist whose target symbol is not the construct in this "
                        "pipeline directory. Only sane for a --dry-run of the payload shape")
    s.add_argument("--token-exponent", type=float, default=2.0,
                   help="ASSUMED scaling of co-folding cost with total token count, used only "
                        "to size the budget guard for a polymer binder. Not measured here")
    s.add_argument("--pocket-distance", type=float, default=6.0)
    s.add_argument("--constrain-pocket", action="store_true",
                   help="restore PROJECT_GOAL.md F3's pocket constraint. Measured to "
                        "destroy discrimination on this target - see the comment in "
                        "cmd_submit before using it")

    s = sub.add_parser("collect"); s.set_defaults(fn=cmd_collect)
    s.add_argument("--poll-seconds", type=int, default=30)
    s.add_argument("--max-polls", type=int, default=120)

    s = sub.add_parser("score"); s.set_defaults(fn=cmd_score)
    s.add_argument("--designs", default=None,
                   help="tag of a BoltzGen design run for THIS target, e.g. --designs v4. "
                        "Omit to score against the pocket signature only. It used to "
                        "default to v4 (KDR), which silently scored other targets "
                        "against KDR's designs")
    s.add_argument("--core-threshold", type=float, default=0.6)
    s.add_argument("--out", default=None,
                   help="write the board here instead of results/repurpose_<pipeline>.json; "
                        "used by tests so they cannot mutate a committed artifact")
    s.add_argument("--rank-by", default="boltzgen_consensus",
                   choices=["boltzgen_consensus", "p2rank_geometry", "known_ligand"])
    s.add_argument("--metric", default="precision_in_core",
                   choices=["precision_in_core", "core_coverage", "weighted_jaccard",
                            "f1_engaged_core", "jaccard_engaged_core"],
                   help="size-normalised by default; see the comment in cmd_score")

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
