# 13 — Scaffold-seeded generation: what seeding buys, and what it only appears to buy

Run of 2026-09-20, thrombin (gene `F2`). Every number in this note was computed by
`demo/seeded_generate.py` and is recorded under `results/demo/seeded/`. Nothing is
estimated.

    ./env-kit/bin/python -m demo.seeded_generate      # the run
    ./env-kit/bin/python -m demo.test_seeded_generate # the self-test: 22/22 pass

Outputs: `10_seeded_generated.json`, `11_seeded_vs_unseeded.json`,
`12_libinvent_generated.json`.

---

## The problem this addresses

`demo/generate.py` samples `entropy/gpt2_zinc_87m`, a SMILES language model that is
**not conditioned on the target**. It samples all of drug-like space, so almost
nothing it produces lands near thrombin chemistry. From its own output
(`results/demo/thrombin/07_generated.json`, n = 591 kept molecules), matched against
the 2,153 approved drugs by ECFP4-2048 Tanimoto:

- best nearest neighbour was **dithranol** (anthralin, for psoriasis) at **0.579**;
- **zero** of the top ten had a thrombin annotation.

The hypothesis under test: if you generate *in the right neighbourhood*, matching
becomes easy. It does. The whole question is whether that means anything.

## The circularity, stated before the results

The seeds are the four approved drugs whose **mechanism** target is F2 —
**argatroban, bivalirudin, dabigatran etexilate, ximelagatran** (derived from
`data/approved_drugs.csv`, not asserted; `seed_drugs()`). Step 1 decomposes them
with BRICS and recombines the fragments. So the products are rearrangements of
approved-drug fragments, which are then matched *back against the approved
library*. A similarity lift is guaranteed by construction. It is arithmetic, not a
result.

Three controls exist to take that lift away again, and the honest reading of this
note is the part that survives them:

| control | what it removes |
|---|---|
| **decoy-seeded panels** | the same BRICS pipeline seeded from MW-matched approved drugs with no F2 annotation. Whatever also shows up here is generic "recombine approved drugs", not thrombin. |
| **leave-the-seeds-out** | the four seed drugs deleted from the corpus before matching. Removes the generator recognising its own parents. |
| **LibInvent (step 2)** | scaffold fixed, but R-groups drawn from a ChEMBL-trained prior instead of from approved drugs. |

`seed_atom_fraction` puts a number on the circularity: the fraction of a molecule's
heavy atoms covered by a seed fragment of at least 5 heavy atoms.

**The size floor on that metric is load-bearing and was nearly missed.** With no
floor, the F2 fragment pool contains `OC` and `NC(=N)N`, which occur in most
drug-like molecules, and the meter reads a median of exactly **1.000** for BRICS
products, for LibInvent products, for the unseeded GPT-2 molecules *and* for
aspirin. All four measured. A metric that cannot tell aspirin from a thrombin
fragment assembly measures nothing. At a 5-heavy-atom floor the same four read
0.700 / 0.682 / 0.286 / 0.538. The self-test asserts the degeneracy explicitly, so
the floor cannot be removed silently.

## Results

All five sets matched against the same 2,153 approved drugs with ECFP4-2048
Tanimoto (`scripts/representations.py` `REPRESENTATIONS["morgan"]`, the same
representation the rest of the repo uses). Base rate of an F2 annotation in the
corpus: **12 / 2,153 = 0.56 %**.

| set | n | best NN | median NN | NN is F2-annotated | of top 10 |
|---|---:|---:|---:|---:|---:|
| **seeded** (BRICS from the 4 F2 drugs) | 378 | 0.6026 | 0.3066 | **151 (40.0 %)** | **10** |
| seeded, 4 seeds deleted from the corpus | 378 | 0.4182 | 0.2872 | 23 (6.1 %) | 2 |
| **LibInvent** (REINVENT4 scaffold decoration) | 81 | 0.5789 | 0.3725 | 1 (1.2 %) | 0 |
| LibInvent, 4 seeds deleted | 81 | 0.5789 | 0.3725 | 1 (1.2 %) | 0 |
| **decoy-seeded**, 3 panels pooled | 726 | **0.8200** | 0.3623 | 1 (0.14 %) | 0 |
| **unseeded** `gpt2_zinc_87m` (baseline) | 591 | 0.5789 | 0.2639 | 5 (0.85 %) | 0 |

Tests (Mann-Whitney one-sided on the NN-Tanimoto distributions, Fisher exact
one-sided on the annotated-NN counts):

| comparison | Tanimoto p | rank-biserial | annotated-NN p |
|---|---:|---:|---:|
| seeded > unseeded | 1.25 × 10⁻⁴⁵ | +0.537 | 3.88 × 10⁻⁶⁴ |
| seeded > decoy-seeded | 1.00 | **−0.452** | 1.80 × 10⁻⁷⁹ |
| decoy-seeded > unseeded | 7.30 × 10⁻¹¹⁷ | +0.734 | 0.992 |
| LibInvent > unseeded | 1.20 × 10⁻³⁴ | +0.837 | 0.539 |
| LibInvent > decoy-seeded | 0.250 | +0.046 | 0.191 |

### Reading 1 — the Tanimoto improvement is entirely generic

Seeding does raise nearest-neighbour similarity, decisively (median 0.3066 vs
0.2639, p = 1.25 × 10⁻⁴⁵). And **it is not a target-specific effect at all**:
BRICS seeded from *random, F2-free* approved drugs does better still — median
0.3623, best 0.8200 against the thrombin-seeded arm's 0.6026, rank-biserial
−0.452 against it. Recombining fragments of approved drugs produces molecules that
resemble approved drugs, and it makes no difference which approved drugs you start
from. **None of the Tanimoto lift is evidence about thrombin.**

This is also why "best nearest neighbour" is a bad headline metric for a seeded
generator, and the reason the decoy panels exist. Without them, "0.60 to argatroban
beats the unseeded model's 0.58 to dithranol" would have read as a result.

### Reading 2 — the annotation enrichment is real, and mostly circular

The count that does separate the arms is *how often the nearest neighbour carries
an F2 annotation*: 40.0 % seeded vs 0.85 % unseeded vs 0.14 % decoy-seeded
(p = 1.8 × 10⁻⁷⁹ against the decoy control). That is a genuinely target-specific
signal — the decoy arm, which enjoys the *larger* Tanimoto lift, does not have it.

Then look at which drugs those neighbours are:

| set | F2-annotated nearest neighbours |
|---|---|
| seeded | argatroban ×115, ximelagatran ×20, captopril ×16 |
| seeded, seeds deleted | captopril ×23 |
| LibInvent | nafamostat ×1 |
| decoy-seeded | cianidanol ×1 |
| unseeded | captopril ×2, apixaban ×1, betrixaban ×1, bortezomib ×1 |

**135 of the 151 hits are the generator recognising its own parents.** All ten of
the seeded arm's top ten are argatroban, a seed. Deleting the four seeds from the
corpus takes 151 down to 23 — still 7.6× the unseeded rate (p = 2.75 × 10⁻⁶) and
47× the decoy rate (p = 2.06 × 10⁻¹⁰) — but every one of those 23 is **captopril**,
an ACE inhibitor whose F2 entry is an off-target annotation, not a thrombin drug.

And the leave-seeds-out control has a hard ceiling that has to be said out loud:
**all four of the corpus's F2-*mechanism* drugs are the seeds**, so once they are
deleted a mechanism-level hit is impossible by construction. What remains reachable
is the 8 off-target F2 annotations (0.37 % of the remaining corpus). The control
therefore cannot distinguish "found a real thrombin analogue" from "found nothing";
it can only rule out the circular hits, which it does.

### Reading 3 — LibInvent is the only non-circular hit, and it is n = 1

LibInvent's single F2-annotated neighbour is **nafamostat**, a guanidinobenzoate
serine-protease inhibitor — chemically the right answer for an S1 pocket that ends
in Asp189, and it is *not* a seed, so it survives the leave-seeds-out control
unchanged. Its R-groups came from a ChEMBL prior rather than from the approved
library.

That is one molecule out of 81, and **it is not significant**: 1/81 vs the unseeded
5/591 gives p = 0.539, and LibInvent's Tanimoto distribution is statistically
indistinguishable from the decoy-seeded arm's (p = 0.250). Reported because it is
the result, not because it is a good one. The obvious next step — more LibInvent
samples — is cheap and was not run.

The gate is the binding constraint on that n: 508 raw samples became 81, with 304
dropped on QED < 0.4 and 122 on MW. LibInvent decorating a small scaffold with
large R-groups produces a lot of things the drug-likeness filter rejects. The
seeded BRICS arm lost 202 of 600 the same way.

### Circularity, quantified

| generator | median `seed_atom_fraction` |
|---|---:|
| seeded (BRICS) | 0.700 |
| LibInvent | 0.682 |
| unseeded `gpt2_zinc_87m` | 0.286 |

Both seeded arms carry roughly 70 % of their heavy atoms from seed fragments,
against a 28.6 % background from a model that never saw a seed drug. LibInvent is
no better on this measure, because its *scaffold* is itself a seed fragment; what
differs is where the decorations come from, and that is exactly where its one
non-circular hit came from.

Motif presence — the arginine mimetic (amidine / guanidine / benzamidine) that
thrombin's S1 pocket requires:

| generator | carries an arginine mimetic |
|---|---:|
| seeded (BRICS) | 146 / 378 (38.6 %) |
| LibInvent | 81 / 81 (100 %, the scaffold pins it) |
| unseeded `gpt2_zinc_87m` | **1 / 591 (0.17 %)** |

Fisher, seeded > unseeded: OR = 371, **p = 1.66 × 10⁻⁶⁷**. This is the one place
where seeding unambiguously does what it was supposed to do: the unseeded model
essentially never proposes the pharmacophore this target needs. Whether that
translates into a useful molecule is what the rest of this note answers, and the
answer is "not yet".

### The bottom line

Seeding made the matching problem easier and did not make it more informative.
Of the improvement over the unseeded baseline:

- the **Tanimoto** part is **0 % target-specific** — the decoy control beats it;
- the **annotation** part is **~89 % circular** (135 of 151 hits are seed drugs),
  and the residue is one off-target annotation repeated 23 times;
- the only hit that is neither circular nor a seed is nafamostat, **n = 1, p = 0.54**.

The honest summary is that the unseeded generator's failure was real, seeding fixes
the *symptom* (nothing near thrombin chemistry) without producing evidence of a
thrombin-relevant molecule, and the experiment's main product is the control design
that shows this.

### Not evaluated

| thing | why |
|---|---|
| Docking the seeded molecules | costs Rowan credits; out of scope for this note, which is about generation |
| REINVENT4 Mol2Mol (similarity-constrained analogues) | prior downloadable from the same Zenodo record; not run |
| REINVENT4 LinkInvent | needs two warheads per line; no warhead pair defined for this target |
| REINVENT4 staged/reinforcement learning | needs a scoring function; a docking-based one would need credits, and a Tanimoto-based one would be circular by construction |
| More LibInvent samples to power the n = 1 hit | ran out of session, not of method |
| Targets other than thrombin | only F2 has 4 mechanism drugs to seed from in this corpus |


## Step 2 — REINVENT4 on Apple silicon

**It installs.** Evidence, on this machine (macOS arm64, Darwin 26.6.2):

    reinvent 4.8.24 | torch 2.12.0 | torch.backends.mps.is_available() True | arm64

- Source: `github.com/MolecularAI/REINVENT4`, Apache-2.0, commit fetched
  2026-09-20 via `git clone --depth 1`.
- Priors: the README's DOI `10.5281/zenodo.15641296` is a *concept* DOI. Resolved
  through `doi.org` it redirects to **`zenodo.org/records/20701824`**, which is the
  record id any download script needs. It holds seven priors; this note uses
  `libinvent.prior` (94.8 MB).

**Two install traps, both hit:**

1. **`uv pip install ./REINVENT4` fails on macOS.** The repo's `pyproject.toml`
   carries `[tool.uv.sources] torch = [{ index = "pytorch" }]` pointing at
   `download.pytorch.org/whl/cu126`, so uv resolves `torch==2.12.0+cu126`, which has
   no `macosx_*_arm64` wheel, and the whole resolution fails. Plain
   `pip install ./REINVENT4` ignores `[tool.uv.sources]`, picks up `torch==2.12.0`
   from PyPI, and succeeds in **42 s**. Their own `install.py mac` shells out to
   `pip`, which is why it works for them and not for uv.
2. **`scipy` is missing from the dependency list.** `reinvent` imports
   `scipy.stats.gaussian_kde` at module import time via
   `runmodes/utils/plot.py`, but `scipy` appears nowhere in `pyproject.toml`. Every
   run mode dies at import until it is installed by hand.

REINVENT4 pins `torch==2.12.0` and `env-kit` is on 2.14, so it gets its own
interpreter and is driven as a subprocess. To reproduce:

    uv venv --python 3.12 --seed env-reinvent
    env-reinvent/bin/python -m pip install /path/to/REINVENT4   # pip, NOT uv pip
    uv pip install --python env-reinvent/bin/python scipy
    mkdir -p priors && curl -L -o priors/libinvent.prior \
      https://zenodo.org/api/records/20701824/files/libinvent.prior/content

`demo/seeded_generate.py` finds them at `env-reinvent/bin/python` and
`priors/libinvent.prior`, overridable with `--reinvent-python` / `--reinvent-prior`
or the `REINVENT_PYTHON` / `LIBINVENT_PRIOR` environment variables. If either is
absent the arm records itself as `not_evaluated` with the reason rather than being
skipped in silence.

### The scaffolds

LibInvent decorates a fixed scaffold with two attachment points. The three used
here are each **verified at run time to be a substructure of a seed drug**
(`verify_scaffolds()`), so "derived from the known binders" is computed, not
claimed:

| scaffold | SMILES | contained in |
|---|---|---|
| `arginine_guanidine` | `[*:0]C(CCCNC(=N)N)[*:1]` | argatroban, bivalirudin |
| `benzamidine_para_amine` | `[*:0]NC(=N)c1ccc(N[*:1])cc1` | dabigatran etexilate |
| `benzylamine_amidine` | `[*:0]Cc1ccc(C(=N)N[*:1])cc1` | ximelagatran |

Each pins an arginine mimetic, which is what thrombin's S1 pocket needs — Asp189 at
the bottom of the pocket is why every approved thrombin inhibitor carries an
amidine, guanidine or amidoxime. A candidate fourth scaffold
(`[*:0]c1ccc(C(=N)N)cc1[*:1]`) matched **none** of the seeds and was dropped; the
self-test keeps a deliberately bogus scaffold around to prove the check still bites.

**Not reproducible bit-for-bit.** REINVENT4's `sampling` run mode exposes no RNG
seed, so the molecule list changes between runs. The `n` and the statistics are
stable; the individual SMILES are not. Recorded in the output as `determinism`.

## Limitations, stated rather than worked around

- **Nothing here is docked.** The pocket enters only through the choice of seeds
  and scaffolds. Tanimoto similarity is not affinity and implies nothing about
  binding.
- **BRICS is enumeration, not learning.** It cannot propose a substructure that was
  not already in a seed drug. "Novel" here means only "not in the approved library
  by InChIKey".
- **Decoy panels are matched to the three small-molecule seeds only.** Bivalirudin
  is a 20-residue peptide (MW 2,180); BRICS on a peptide yields interchangeable
  amide fragments whose recombination generates duplicates indefinitely. Two early
  runs wedged on exactly that, which is why `brics_enumerate` now runs the build in
  a child process under a hard wall-clock deadline and `DECOY_MAX_MW` caps the
  matching at 900.
- **The BRICS motif fraction falls as the run gets longer** (78.6 % at n = 60,
  38.6 % at n = 600) because `BRICSBuild` walks depth-first and the amidine-rich
  products come out early. Any threshold tuned on a short run is a trap; the
  self-test asserts a quarter, not a majority, for this reason.
