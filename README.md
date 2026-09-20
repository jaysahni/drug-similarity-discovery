# Rebind

**Given a disease, find the approved drugs that fit the target's binding site — including
the ones that look nothing like what you'd expect.**

Built at HackMIT 2026 by Evan Xiang, Jay Sahni, Richard Shan and Rachel Chen. MIT
licensed. Every number here was computed by a script in this repo, and you can reproduce
any of them by running it.

## Inspiration

Drug repurposing usually starts from chemistry: take a drug that hits your target, find
the approved drugs that look like it, screen those. It's a reasonable instinct and it
works better than chance.

It also misses most of the answer. We ranked the approved-drug corpus by chemical
similarity to axitinib, the best-characterised VEGFR2 ligand, and looked for the 32
approved drugs known to hit that target. Similarity genuinely helps — 10× enrichment over
base rate. And it still buried sunitinib at rank 464, pazopanib at 764 and vandetanib at
1,760. Those are real, approved VEGFR2 inhibitors that any chemistry-first shortlist would
throw away.

The reason is simple once you see it. **Two drugs that share a binding site don't have to
share any chemistry at all.** So we stopped asking what a molecule looks like and started
asking which residues it touches. That question has an answer you can compare across
molecules with nothing in common — even between a designed miniprotein and a small
molecule, which is what made the whole generative arm possible.

## What it does

Rebind takes a disease and returns a ranked board of approved drugs that engage the same
site on a disease-relevant target, with the evidence for every row attached.

```
disease  ->  ranked targets (literature evidence, weighed by an LLM that can't invent facts)
         ->  structure + binding site
         ->  design binders against that site  ->  consensus INTERFACE SIGNATURE
         ->  co-fold every approved drug with the target (Boltz-2)
         ->  score each pose on how much of the site it actually engages
         ->  a ranked board you can interrogate
```

Ask it something. This runs offline, instantly, with no credits and no GPU:

```bash
./env/bin/python scripts/ask.py --list             # what can be asked about
./env/bin/python scripts/ask.py "colorectal cancer"
./env/bin/python scripts/ask.py --target KDR --json
./env/bin/python scripts/ask.py --explain sunitinib
```

Run it live against a new target:

```bash
./env/bin/python scripts/autorepurpose.py run --target KDR --uniprot P35968
```

The score is `precision_in_core`: of the residues a drug engages, the fraction lying
inside the site core, at a 4.5 Å heavy-atom cutoff in author numbering. Every row names
the residues engaged and the ones missed, so you can open the structure and argue with it.
That's the part we care about most — the output is a hypothesis with its reasoning
attached, not a number you have to take on faith. A drug whose pose has no drug-like
ligand keeps its row with a status and no score, rather than a plausible-looking number.

Three surfaces sit on top of the results: the offline `ask.py` front door, the **Rebind**
web workspace with chats, folders, saved results and a 3Dmol viewer that highlights
engaged versus missed residues on the real complex, and a self-contained HTML report.

<!--LANG-EXEMPT-->
One deliberate omission: **no affinity, Kd, IC50 or potency is predicted or read anywhere
in this pipeline.** Boltz-2 can emit an affinity score and we never read it. Ranking is
interface overlap only, and a CI lint keeps that vocabulary out of the repo. Everything
here ranks structural hypotheses about shared site engagement — nothing predicts binding
strength or clinical benefit, and rows labelled "known binder" are the benchmark's
positive controls rather than discoveries.
<!--/LANG-EXEMPT-->

## How we built it

Python 3.14, RDKit, BioPython, NumPy and SciPy. P2Rank for pocket detection. Rowan as the
hosted GPU backend for BoltzGen design, Boltz-2 co-folding and Vina docking, since we had
no local GPU. ESM-C for peptide embeddings, ChemBERTa for the learned chemical baseline,
a ZINC-trained GPT-2 for molecule generation. DrugCentral for the approved corpus, ChEMBL
for bioactivity, Open Targets and Europe PMC for disease evidence. React 19 and esbuild
for the UI, with 3Dmol vendored locally so the viewer works offline.

We ran four experiments in order, each one gating the next.

**Which representation actually retrieves drugs that share a target?** 16 representations,
2,069 queries, a random floor and a perfect-ranking ceiling, bootstrap confidence intervals
and Holm correction across the whole comparison family. ECFP4 won and we used it
downstream because it measured best, not because it's conventional.

**Is the target-side idea sound?** Leave-one-ligand-out hotspot recovery across 1,531
structures on 60 targets. Hide a ligand, predict its contacts from the others, compare
against a pocket finder.

**Does a generated design beat a free pocket finder?** The decisive ablation, on 38 KDR
co-crystals with every arm scored on the same 38 so the comparisons are paired.

**The product run.** 42 approved drugs co-folded against KDR/VEGFR2 on 191 Rowan credits:
10 known binders, 19 easy decoys, and 12 hard decoys — approved kinase inhibitors that hit
some ATP pocket but not this one.

Alongside that, `demo/` re-implements the entire pipeline in ~1,500 lines on a
cheminformatics toolkit and runs it on CDK2 and thrombin, asking the inverse question:
given a designed binder, which approved drug resembles it?

## Individual Contributions

**Evan Xiang** built the scientific core: the retrieval benchmark and representation zoo,
the metrics layer, target-side contact extraction, the pocket-finder baseline, hotspot
recovery, the BoltzGen arm and the diagnostics that made it aim, the decisive ablation,
the repurposing pipeline end to end, the hard-decoy test, the `ask.py` interface and the
CI language lint.

**Jay Sahni** wrote the project spec and the entire lightweight demo: the direct-matching
harness and its calibrated null, the peptide arm with ESM-C matching and the hub check
that explained its result, the small-molecule generation and docking arm, the approved-
biologics corpus, the insertion-code fix that recovered thrombin's active site, and the
final findings report.

**Richard Shan** handled integration and evidence discipline: the LLM judge step with the
hallucination guard that rejects any citation or gene symbol outside the verified evidence
set, the visualization handoff bundle and wiring results through to the viewer, decoy
relabelling from measured ChEMBL activity, the pipeline and GPU docs, the retractions, and
the cross-branch merges that kept four workstreams on one trunk.

**Rachel Chen** built the Rebind UI: the React pipeline diagram and the workspace layer of
chats, folders, saved results and compact navigation.

## Challenges we ran into

**The decisive experiment came out against us.** On 38 KDR co-crystals a free pocket
finder reaches Jaccard 0.6355 against the ligand's true contacts, while a 24-design,
85-credit BoltzGen consensus reaches 0.3115 (Holm p=0.0016). It replicated on CDK2. The
pocket finder lands within 0.03 of the known-ligand ceiling, so there is barely any
headroom on this target for anything to occupy. Most of the second half of the project was
establishing that carefully enough to report it as a result rather than a failure — and
the pipeline now ships with the cheap arm, which is the better answer.

**Six bugs that produced perfectly ordinary-looking wrong answers.** None threw an error. A
parser ignoring insertion codes silently dropped 28 of thrombin's 259 residues, including
the entire loop lining the active site. Structure prep renamed chains, so downstream stages
pointed at nothing. A co-crystal ligand left in the pocket sent every design to the wrong
site. And a percentile using a strict `<` charged every exact tie against the candidate,
biasing every headline number low and invisibly, until a calibration test caught it. All
six share one root cause — trusting an identifier an upstream tool is free to change — and
all six are now guarded by tests.

**A hardware and packaging floor we couldn't raise.** No local GPU, and every GPU tier on
our hosted fallback was refused pending a payment method. The free credit tier was 500 plus
20 a week. `boltz` won't install on Python 3.14, Vina publishes no wheel for it, and every
pocket-conditioned generator we tried depends on `torch-scatter`, which has no build for
our torch. So the generative arm ended up unconditioned plus a docking filter, which is a
weaker thing than we set out to build, and we said so rather than papering over it.

**BoltzGen doesn't do what its interface implies.** It indexes residues 1..N over residues
present in the file, ignoring author numbering, so our construct's pocket residues were all
out of range and the first run died immediately. And naming the hotspot residues you want
does almost nothing: 0, 2, 0, 0 of 17 designs reached the requested pocket, with an empty
consensus. Restricting the presented surface with a 12 Å proximity window is what worked,
taking designs on target to 6–15 of 17 and consensus core residues from 0 to 29.

**A headline that was wrong twice over.** An earlier version of this README compared a rank
out of a 41-drug board against a rank out of a 4,099-drug corpus and claimed a 100× gap.
The ratio of those populations is almost exactly 100, so the comparison was measuring
population sizes. It also cherry-picked the four drugs with the largest gap. We retracted
it in place and re-ran the honest version rather than quietly dropping it.

**The hard test still isn't passed.** Our enrichment separates known binders from
chemically unrelated drugs very well, but against 12 approved kinase inhibitors that miss
this target, separation isn't significant (AUC 0.717, p=0.092). Worth knowing that two of
those hard decoys are promiscuous inhibitors whose annotations are probably incomplete, so
that number is a lower bound rather than a measurement.

## Accomplishments that we're proud of

**The pipeline runs end to end, and seven of its top ten are known binders.** 42 drugs
co-folded, 41 scored, ranked purely by how much of the binding site each engages.
Enrichment 2.87× at the 25% cut against an arithmetic maximum of 4.10×, median known-binder
rank 8th of 41.

**It surfaces the drugs chemistry throws away.** This is the result we set out to get:

| drug | by binding site | by chemical similarity |
|---|---|---|
| **pazopanib** | **4** | 764 |
| **sunitinib** | **5** | 464 |
| **nintedanib** | **9** | 755 |
| **vandetanib** | **10** | 1,760 |

Structural ranking finds them because it isn't looking at what they're made of. It's
looking at where they sit.

**And it beats the chemical baseline measured properly.** Re-ranking the same 41 drugs both
ways, so the populations match, structural wins on median rank (8.0 vs 12.5) and on AUC
against property-matched decoys (0.913, p=0.0003, vs 0.795, p=0.011). A real win of roughly
1.5×.

**The target-side idea holds up under attack.** Across 1,531 structures on 60 targets, a
consensus of which residues other ligands engage predicts a held-out ligand's contacts at
0.7322 precision against a pocket finder's 0.4616 — +0.27, Wilcoxon p=7.9e-113, and it
holds when every target counts once. Then one of us wrote a 15-finding adversarial critique
against it before anyone was allowed to believe it, and it survived every blocking check:
alternative set rules, three contact cutoffs, resolution confounds, lipid exclusions.

**We decomposed our own win against a free null.** Predicting from one randomly chosen
other ligand already scores 0.6193. So of that +0.27, about +0.16 is knowing that ligands
bind here at all, and +0.11 is what aggregating many binding events actually buys. That
told us more about what the method does than the headline did.

**A benchmark with a floor and a ceiling.** ECFP4 reaches 66% of the way to a perfect
ranking on p@1 from a floor of 7%, and a 2048-bit hash of substructures beats a
77M-parameter transformer at this task. The interesting part was the disagreement: ECFP4
wins at putting one correct drug first while four other representations beat it on
whole-list quality. Reporting one number would have called it a clean sweep and hidden
something real.

**Right-sizing the binder nearly tripled design confidence**, from ipTM 0.162 to 0.598.
Obvious in hindsight — a 60–90-mer can't fit inside an ATP slot, an 8–16-mer sits
comfortably in a protease groove.

**Every arm runs a positive control on the same data with the same metric**, because a
failed search and a broken instrument look identical from outside. Bivalirudin retrieves
lepirudin at rank 1 of 36. Docking ranks three real thrombin drugs above the generated
median in the same batch. Every control passes.

**Honesty as an engineering constraint rather than a disclaimer.** CI lints the vocabulary.
The gate artifact records "neither passed nor failed" rather than claiming the stronger
ablation we didn't run. `results/` holds computed output only, never hand-edited.

## What we learned

**Run a power check before reporting any null.** Our docking arm was written up as "docking
adds nothing" at n=23. At n=199 the same effect is significant. The result existed the
whole time, the molecules cost about 0.1 credit each, and we simply didn't buy the power.
An underpowered null is not a negative result.

**Controls are what make a negative readable.** Several of our findings are negative, and
each is worth reporting only because a control on the same data with the same metric
passes. Without them we'd have had a pile of shrugs.

**How you read a result can invert it.** Read as a top-5 list, our chemical baseline looked
useless on thrombin: the five nearest drugs to argatroban are all peptidomimetics with no
thrombin annotation. Read as ranks, it works — bivalirudin comes back at rank 10 of 4,098
against ~2,050 expected by chance. The top-five reading nearly produced a false negative
about our own control.

**Decompose a win against a free baseline before believing it.** Most of what any arm
predicted was simply that the pocket is the pocket, which a pocket finder returns for
nothing.

**Re-derive identity from content, never from an identifier upstream can change.** One root
cause, six separate incidents, six separate fixes. We should have made it a stated
invariant with a test the first time.

**Size the tool to the site.** Beyond the binder-length fix, pocket detection can't find a
serine protease active site at all: six pockets on thrombin and not one containing the
catalytic residues, because that site is a shallow groove across subsites rather than an
enclosed cavity.

**Establish the question before building the pipeline.** We built footprint matching, with
the designed binder as a discarded probe, before establishing that the intent was direct
matching, with the binder as the query. One question at the outset would have saved roughly
half our working time.

## What's next for our project

1. **A genuinely pocket-conditioned generator.** Our molecule arm generated unconditionally
   and used the pocket only as a docking filter, which is why finding nothing says little.
   DiffSBDD is reachable — MIT licensed, public checkpoints — and its only blocker needs
   roughly a ten-line shim.
2. **The filtered design ensemble, which is genuinely untested.** Zero of our 24 designs
   passed even ipTM 0.5, so what we ran was an unfiltered consensus. Reaching ten filtered
   survivors needs roughly 1,188 credits against the 500 we had. Whether that closes the
   gap is not answered by our result, and shouldn't be claimed either way.
3. **Compare the two arms.** Two independent rankings of approved drugs for the same target
   now sit on disk and have never been compared. Where they agree is more interesting than
   either alone, and the data is already there.
4. **Replicate on a second target end to end.** CDK2 is prepared and never screened. The
   gate wasn't trusted until it replicated; the newer arms are still one target each.
5. **Make the report target-agnostic.** Its inputs are currently hardcoded to KDR.
6. **Correct the spec.** Our plan document specifies a design-confidence filter that this
   repo's own investigation showed to be unreachable by construction, and that has no
   provenance in the underlying model's paper.
7. **Build the golden fixtures** the spec already mandates, whose absence is why a
   statistics block silently vanished the first time a stage re-ran.
