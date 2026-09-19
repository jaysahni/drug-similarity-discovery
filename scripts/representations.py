"""Representation zoo for the drug-similarity benchmark.

Every representation is a function  list[smiles] -> ndarray (n, d)  registered in
REPRESENTATIONS by the @representation decorator. Each carries two attributes the
benchmark reads:

    fn.binary : bool   True  -> compare with Tanimoto (bit vectors)
                       False -> compare with cosine
    fn.blurb  : str    one line describing what it encodes, for the results table

Rows are never dropped. A SMILES that will not parse - or a molecule an individual
representation cannot encode (a conformer that will not embed, say) - gets a zero
row, and its index is recorded in FAILURES[name] = (n_failed, [indices]) so the
count can be reported rather than silently lost. All matrices stay row-aligned to
the input SMILES list.

The expensive ones (usrcat conformers, chemberta forward passes) cache their matrix
under data/raw/repcache/ keyed by a hash of the SMILES list, so a benchmark rerun
is instant. Set REPZOO_NO_CACHE=1 to force recomputation.

Usage:
    ./env/bin/python scripts/representations.py          # self-test over data/drugs.csv
    ./env/bin/python scripts/representations.py morgan maccs
"""

from __future__ import annotations

import csv
import functools
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "drugs.csv"
RESULTS = ROOT / "results"
CACHE = ROOT / "data" / "raw" / "repcache"
SEED = 0
FPSIZE = 2048

# ChemBERTa: first choice, then the fallback used if the first will not download.
CHEMBERTA_MODELS = ("DeepChem/ChemBERTa-77M-MLM", "seyonec/ChemBERTa-zinc-base-v1")
CHEMBERTA_BATCH = 64

# Gobbi 2D pharmacophore is cubic in matched features; see the comment in build().
GOBBI_MAX_HEAVY_ATOMS = 150

# Things deliberately left out of the registry, with the reason, per CLAUDE.md.
NOT_EVALUATED = {
    "gobbi_pharm2d over 150 heavy atoms": (
        "23 of 2114 drugs exceed the cap and get a zero row (counted in n_failed). "
        "Measured: <10 ms at the median 24 heavy atoms, 31-60+ s each at 279-342 atoms"
    ),
    "multi-conformer usrcat": (
        "one ETKDGv3 conformer per molecule only; a conformer ensemble is the correct "
        "way to use USRCAT but costs ~20x the CPU and was not run"
    ),
    "GPU-backed molecular encoders": (
        "no GPU and no hosted-inference credential on this machine (docs/03-SCOPE-AND-CONSTRAINTS.md); "
        "chemberta runs on CPU, anything larger does not"
    ),
}


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------
REPRESENTATIONS = {}
FAILURES = {}  # name -> (n_failed, [failing indices])

_running = []  # stack of (name, [failing indices]) for the representation in flight


def _fail(i):
    """Record that input molecule `i` could not be encoded by the running representation."""
    if _running:
        _running[-1][1].append(int(i))


def representation(name, binary=False, blurb=None):
    """Register a representation. Same call API as the decorator in benchmark.py."""

    def deco(fn):
        @functools.wraps(fn)
        def wrapper(smiles):
            _running.append((name, []))
            try:
                mat = fn(smiles)
            finally:
                _, bad = _running.pop()
                FAILURES[name] = (len(bad), sorted(set(bad)))
            assert mat.shape[0] == len(smiles), f"{name} returned {mat.shape[0]} rows for {len(smiles)} smiles"
            return mat

        wrapper.binary = binary
        wrapper.blurb = blurb or (fn.__doc__ or "").strip().splitlines()[0]
        REPRESENTATIONS[name] = wrapper
        return wrapper

    return deco


def _mols(smiles):
    """Parse SMILES, recording unparseable ones as failures. Returns list[Mol | None]."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    out = []
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if s else None
        if m is None:
            _fail(i)
        out.append(m)
    return out


def _bitvec_rep(smiles, size, per_mol):
    """Dense float32 (n, size) from a per-molecule callable returning a bit array."""
    out = np.zeros((len(smiles), size), dtype=np.float32)
    for i, m in enumerate(_mols(smiles)):
        if m is None:
            continue
        out[i] = np.asarray(per_mol(m), dtype=np.float32)
    return out


def _zscore(mat, bad):
    """Z-score columns over the rows that succeeded; failed rows stay exactly zero."""
    keep = np.ones(mat.shape[0], dtype=bool)
    keep[list(bad)] = False
    ref = mat[keep] if keep.any() else mat
    sd = ref.std(0)
    sd[sd == 0] = 1.0
    out = (mat - ref.mean(0)) / sd
    out[~keep] = 0.0
    return np.nan_to_num(out, posinf=0.0, neginf=0.0).astype(np.float32)


# --------------------------------------------------------------------------
# cache
# --------------------------------------------------------------------------
def _key(smiles, tag):
    h = hashlib.sha256(("\n".join(smiles)).encode()).hexdigest()[:16]
    return CACHE / f"{tag}-{h}.npz"


def _cached(tag, smiles, build):
    """Load (mat, failed) from disk if present, else build() and store it."""
    path = _key(smiles, tag)
    if path.exists() and os.environ.get("REPZOO_NO_CACHE") != "1":
        z = np.load(path)
        for i in z["failed"]:
            _fail(i)
        return z["mat"]
    mat, failed = build()
    for i in failed:
        _fail(i)
    CACHE.mkdir(parents=True, exist_ok=True)
    # write-then-rename: another session running this module must never read a half-written file
    tmp = path.with_suffix(f".{os.getpid()}.tmp.npz")
    np.savez_compressed(tmp, mat=mat, failed=np.asarray(sorted(failed), dtype=np.int64))
    tmp.replace(path)
    return mat


# --------------------------------------------------------------------------
# 2D fingerprints
# --------------------------------------------------------------------------
@representation("morgan", binary=True)
def morgan(smiles):
    """Morgan (ECFP4) 2048-bit fingerprints - the baseline every model must beat."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            _fail(i)
            continue
        out[i] = np.asarray(gen.GetFingerprint(m), dtype=np.float32)
    return out


@representation("fcfp4", binary=True)
def fcfp4(smiles):
    """Feature-Morgan (FCFP4) 2048 bits: same radius as ECFP4 but pharmacophoric atom types."""
    from rdkit.Chem import rdFingerprintGenerator as rfg

    gen = rfg.GetMorganGenerator(
        radius=2, fpSize=FPSIZE, atomInvariantsGenerator=rfg.GetMorganFeatureAtomInvGen()
    )
    return _bitvec_rep(smiles, FPSIZE, gen.GetFingerprintAsNumPy)


@representation("morgan_count")
def morgan_count(smiles):
    """Count-based ECFP4 (2048 hashed counts, cosine) - substructure multiplicity, not just presence."""
    from rdkit.Chem import rdFingerprintGenerator as rfg

    gen = rfg.GetMorganGenerator(radius=2, fpSize=FPSIZE)
    return _bitvec_rep(smiles, FPSIZE, gen.GetCountFingerprintAsNumPy)


@representation("maccs", binary=True)
def maccs(smiles):
    """MACCS keys: 167 hand-written SMARTS substructure questions (bit 0 unused)."""
    from rdkit.Chem import MACCSkeys

    return _bitvec_rep(smiles, 167, MACCSkeys.GenMACCSKeys)


@representation("rdkit_fp", binary=True)
def rdkit_fp(smiles):
    """RDKit path fingerprint, 2048 bits: hashed linear subgraphs up to 7 bonds."""
    from rdkit.Chem import rdFingerprintGenerator as rfg

    gen = rfg.GetRDKitFPGenerator(fpSize=FPSIZE)
    return _bitvec_rep(smiles, FPSIZE, gen.GetFingerprintAsNumPy)


@representation("atom_pair", binary=True)
def atom_pair(smiles):
    """Atom-pair fingerprint, 2048 bits: (atom type, atom type, topological distance) triples."""
    from rdkit.Chem import rdFingerprintGenerator as rfg

    gen = rfg.GetAtomPairGenerator(fpSize=FPSIZE)
    return _bitvec_rep(smiles, FPSIZE, gen.GetFingerprintAsNumPy)


@representation("topological_torsion", binary=True)
def topological_torsion(smiles):
    """Topological-torsion fingerprint, 2048 bits: hashed 4-atom paths - local backbone shape."""
    from rdkit.Chem import rdFingerprintGenerator as rfg

    gen = rfg.GetTopologicalTorsionGenerator(fpSize=FPSIZE)
    return _bitvec_rep(smiles, FPSIZE, gen.GetFingerprintAsNumPy)


@representation("avalon", binary=True)
def avalon(smiles):
    """Avalon fingerprint, 1024 bits: a curated enumeration of paths, rings and atom environments."""
    from rdkit.Avalon import pyAvalonTools

    return _bitvec_rep(smiles, 1024, lambda m: pyAvalonTools.GetAvalonFP(m, nBits=1024))


@representation("pattern_fp", binary=True)
def pattern_fp(smiles):
    """RDKit pattern fingerprint, 2048 bits: substructure-screening bits, deliberately coarse."""
    from rdkit import Chem

    return _bitvec_rep(smiles, FPSIZE, lambda m: Chem.PatternFingerprint(m, fpSize=FPSIZE))


@representation("layered_fp", binary=True)
def layered_fp(smiles):
    """RDKit layered fingerprint, 2048 bits: path bits hashed over several atom/bond abstraction layers."""
    from rdkit import Chem

    return _bitvec_rep(smiles, FPSIZE, lambda m: Chem.LayeredFingerprint(m, fpSize=FPSIZE))


# --------------------------------------------------------------------------
# pharmacophore
# --------------------------------------------------------------------------
@representation("gobbi_pharm2d", binary=True)
def gobbi_pharm2d(smiles):
    """Gobbi 2D pharmacophore fingerprint, feature triangles by topological distance bin; >150 heavy atoms skipped."""

    def build():
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

        RDLogger.DisableLog("rdApp.*")
        factory = Gobbi_Pharm2D.factory
        # The signature is ~40k sparse bits; collect on-bits, then keep only the columns
        # some molecule actually sets. Dropping all-zero columns is lossless for both
        # Tanimoto and cosine and keeps the matrix from being 340 MB of mostly zeros.
        on, bad = [], []
        for i, s in enumerate(smiles):
            m = Chem.MolFromSmiles(s) if s else None
            # Cost is cubic in matched features. Measured on this set: <10 ms at the
            # median 24 heavy atoms, but 31-60+ s each for the 279-342 atom entries,
            # which alone would take ~17 min. Over the cap the molecule gets a zero
            # row and is counted in n_failed rather than silently included.
            if m is None or m.GetNumAtoms() > GOBBI_MAX_HEAVY_ATOMS:
                bad.append(i)
                on.append(())
                continue
            on.append(tuple(Generate.Gen2DFingerprint(m, factory).GetOnBits()))
        used = sorted({b for bits in on for b in bits})
        col = {b: j for j, b in enumerate(used)}
        out = np.zeros((len(smiles), len(used)), dtype=np.float32)
        for i, bits in enumerate(on):
            for b in bits:
                out[i, col[b]] = 1.0
        return out, bad

    return _cached("gobbi", smiles, build)


# --------------------------------------------------------------------------
# 3D shape
# --------------------------------------------------------------------------
def _usrcat_one(arg):
    """One ETKDGv3 conformer + MMFF optimisation + USRCAT. Runs in a worker process."""
    i, s = arg
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(s) if s else None
    if m is None:
        return i, None
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = SEED + i  # fixed per molecule, so the result does not depend on scheduling
    p.useRandomCoords = True
    p.maxIterations = 200
    try:
        if AllChem.EmbedMolecule(m, p) < 0:
            return i, None
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)  # returns 1 if not converged; keep it either way
        return i, list(rdMolDescriptors.GetUSRCAT(m))
    except Exception:  # noqa: BLE001 - a single molecule must not kill the pass
        return i, None


@representation("usrcat")
def usrcat(smiles):
    """USRCAT on one ETKDGv3/MMFF conformer: 60 shape+pharmacophore moments, z-scored, cosine."""
    import multiprocessing as mp

    def build():
        out = np.zeros((len(smiles), 60), dtype=np.float64)
        bad = []
        work = list(enumerate(smiles))
        try:
            with mp.Pool(min(15, os.cpu_count() or 1)) as pool:
                done = list(pool.imap_unordered(_usrcat_one, work, chunksize=16))
        except Exception as e:  # noqa: BLE001 - spawn cannot re-import a stdin-driven caller
            print(f"  usrcat: multiprocessing unavailable ({type(e).__name__}), running serially")
            done = [_usrcat_one(w) for w in work]
        for i, d in done:
            if d is None:
                bad.append(i)
            else:
                out[i] = d
        return _zscore(out, bad), bad

    return _cached("usrcat", smiles, build)


# --------------------------------------------------------------------------
# physicochemical
# --------------------------------------------------------------------------
@representation("rdkit_descriptors")
def rdkit_descriptors(smiles):
    """Physicochemical descriptors, z-scored. Tests whether bulk properties suffice."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")
    names = [n for n, _ in Descriptors.descList]
    calc = dict(Descriptors.descList)
    out = np.zeros((len(smiles), len(names)), dtype=np.float64)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            _fail(i)
            continue
        for j, n in enumerate(names):
            try:
                v = calc[n](m)
            except Exception:  # noqa: BLE001 - individual descriptors can fail
                v = 0.0
            out[i, j] = v if np.isfinite(v) else 0.0
    out = np.nan_to_num(out, posinf=0.0, neginf=0.0)
    sd = out.std(0)
    sd[sd == 0] = 1.0
    return ((out - out.mean(0)) / sd).astype(np.float32)


@representation("mw_logp")
def mw_logp(smiles):
    """Molecular weight + cLogP only, z-scored - a deliberately weak 2-feature control."""
    from rdkit.Chem import Descriptors

    out = np.zeros((len(smiles), 2), dtype=np.float64)
    bad = []
    for i, m in enumerate(_mols(smiles)):
        if m is None:
            bad.append(i)
            continue
        out[i] = (Descriptors.MolWt(m), Descriptors.MolLogP(m))
    return _zscore(out, bad)


# --------------------------------------------------------------------------
# scaffold
# --------------------------------------------------------------------------
@representation("murcko_scaffold", binary=True)
def murcko_scaffold(smiles):
    """One-hot Bemis-Murcko scaffold identity: similarity 1 iff same scaffold, else 0."""
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scafs = []
    for i, m in enumerate(_mols(smiles)):
        if m is None:
            scafs.append(None)
            continue
        try:
            # acyclic molecules give "", which becomes its own (large, blunt) bucket
            scafs.append(MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        except Exception:  # noqa: BLE001
            _fail(i)
            scafs.append(None)
    uniq = sorted({s for s in scafs if s is not None})
    col = {s: j for j, s in enumerate(uniq)}
    out = np.zeros((len(smiles), len(uniq)), dtype=np.float32)
    for i, s in enumerate(scafs):
        if s is not None:
            out[i, col[s]] = 1.0
    return out


# --------------------------------------------------------------------------
# neural
# --------------------------------------------------------------------------
@representation("chemberta")
def chemberta(smiles):
    """ChemBERTa-77M-MLM, mean-pooled last hidden state over non-pad tokens (not the CLS token)."""

    def build():
        import torch
        from transformers import AutoModel, AutoTokenizer

        os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
        torch.manual_seed(SEED)
        torch.set_grad_enabled(False)
        errs = []
        tok = mod = None
        for mid in CHEMBERTA_MODELS:
            try:
                tok = AutoTokenizer.from_pretrained(mid)
                mod = AutoModel.from_pretrained(mid, add_pooling_layer=False)
                break
            except Exception as e:  # noqa: BLE001 - network or hub failure, try the fallback
                errs.append(f"{mid}: {type(e).__name__}: {e}")
        if mod is None:
            raise RuntimeError("no ChemBERTa checkpoint could be loaded -> " + " | ".join(errs))
        mod.eval()

        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        dim = mod.config.hidden_size
        out = np.zeros((len(smiles), dim), dtype=np.float32)
        # The tokenizer will happily embed a string RDKit cannot parse. Gate on RDKit
        # anyway, so every representation is zero on exactly the same rows.
        bad = [i for i, s in enumerate(smiles) if not s or Chem.MolFromSmiles(s) is None]
        for start in range(0, len(smiles), CHEMBERTA_BATCH):
            chunk = [s if s else "" for s in smiles[start : start + CHEMBERTA_BATCH]]
            enc = tok(chunk, padding=True, truncation=True, max_length=512, return_tensors="pt")
            h = mod(**enc).last_hidden_state
            mask = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
            pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
            out[start : start + len(chunk)] = pooled.numpy()
        out[list(bad)] = 0.0
        return out, bad

    return _cached("chemberta", smiles, build)


# --------------------------------------------------------------------------
# floor
# --------------------------------------------------------------------------
@representation("random")
def random_rep(smiles):
    """Random vectors - the floor. Any representation must beat this."""
    rng = np.random.default_rng(SEED)
    return rng.standard_normal((len(smiles), 64)).astype(np.float32)


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------
DETERMINISM_CHECK = ["morgan", "maccs", "murcko_scaffold", "usrcat", "chemberta"]


def load_smiles():
    with DATA.open() as fh:
        return [r["smiles"] for r in csv.DictReader(fh) if r["smiles"]]


def main():
    wanted = sys.argv[1:] or list(REPRESENTATIONS)
    unknown = [w for w in wanted if w not in REPRESENTATIONS]
    if unknown:
        sys.exit(f"unknown representation(s): {unknown}\nknown: {list(REPRESENTATIONS)}")

    smiles = load_smiles()
    n_unparseable = sum(1 for m in _mols(smiles) if m is None)
    print(f"{len(smiles)} SMILES from {DATA}")
    print(f"{n_unparseable} of them do not parse in RDKit\n")

    rows = []
    mats = {}
    for name in wanted:
        fn = REPRESENTATIONS[name]
        t = time.time()
        mat = fn(smiles)
        sec = time.time() - t
        n_failed, idx = FAILURES[name]
        mats[name] = mat
        rows.append(
            {
                "name": name,
                "dim": int(mat.shape[1]),
                "binary": bool(fn.binary),
                "seconds": round(sec, 2),
                "n_failed": n_failed,
                "failed_indices": idx,
                "blurb": fn.blurb,
            }
        )
        print(f"  {name:<20} dim={mat.shape[1]:<7} {sec:7.2f}s  n_failed={n_failed}")

    # determinism: recompute with the cache bypassed and require bit-identical matrices
    print("\ndeterminism check (recomputed with REPZOO_NO_CACHE=1):")
    os.environ["REPZOO_NO_CACHE"] = "1"
    det = {}
    for name in DETERMINISM_CHECK:
        if name not in mats:
            continue
        again = REPRESENTATIONS[name](smiles)
        det[name] = bool(np.array_equal(mats[name], again))
        print(f"  {name:<20} identical: {det[name]}")
    os.environ.pop("REPZOO_NO_CACHE", None)

    width = max(len(r["name"]) for r in rows) + 2
    print(f"\n{'name'.ljust(width)}{'dim':>8}{'binary':>9}{'seconds':>10}{'n_failed':>10}")
    print("-" * (width + 37))
    for r in rows:
        print(
            r["name"].ljust(width)
            + f"{r['dim']:>8}{str(r['binary']):>9}{r['seconds']:>10.2f}{r['n_failed']:>10}"
        )

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "representations_selftest.json"
    out.write_text(
        json.dumps(
            {
                "n_molecules": len(smiles),
                "n_unparseable_smiles": n_unparseable,
                "representations": rows,
                "determinism_identical": det,
                "not_evaluated": NOT_EVALUATED,
            },
            indent=2,
        )
    )
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
