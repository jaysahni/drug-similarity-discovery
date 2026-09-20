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

**Six bugs that produced perfectly ordinary-looking wrong answers.** Not one of them ever
threw an error, which is exactly what made them expensive. A parser that ignored insertion
codes quietly dropped 28 of thrombin's 259 residues, and among the residues it dropped was
the entire loop lining the active site. Structure prep renamed chains partway through, so
the stages downstream went looking for something that was no longer there. A co-crystal
ligand left sitting in the pocket sent every single design off to the wrong site. And a
percentile calculation written with a strict `<` charged every exact tie against the
candidate, which quietly pushed every headline number down until a calibration test finally
caught it. Looking back, all six come from the same bad habit: trusting an identifier that
an upstream tool is free to change underneath you. All six now have tests guarding them.

**The docking API cannot give you back a pose, and it costs you to find out.** We built
an interaction-fingerprint layer to score drugs on *which* residues a pose contacts rather
than a bare similarity number, with six typed interaction channels validated against the
published contacts of a thrombin crystal structure — it reproduces the Asp189 salt bridge,
the oxyanion hole and the 60-loop cage in the right channels. Then we went to fetch the
poses and there were none. `num_poses_to_save` defaults to zero and Rowan's documented
helper does not expose it, so no call through the normal path can keep one. We hand-built
the payload to ask for a pose explicitly, the workflow accepted the field and completed —
and still returned nothing retrievable. Twenty credits to establish that the capability is
not there.

The useful part was what we did before spending it. A simulated power analysis said that
with the three actives we had, the test we wanted to run had a 23% chance of detecting even
a large effect — so re-docking the same ligands would have bought a result that could not
answer the question. We enlarged the active set first. That habit came from earlier in the
project, where we had reported a null at n=23 that turned out to be significant at n=199;
the effect had been there the whole time and we simply could not see it.

**BoltzGen doesn't do what its interface implies.** It numbers residues 1 through N over
whatever happens to be present in the file and ignores author numbering completely, so
every pocket residue we asked for in our construct fell outside the valid range and the
first run died immediately. The subtler problem showed up once that was fixed. Naming the
hotspot residues you want turns out to do almost nothing: across four attempts, 0, 2, 0
and 0 of 17 designs reached the pocket we had asked for, and the consensus came back empty.
What actually worked was restricting the surface we presented to the model with a 12 Å
proximity window. That took designs on target up to somewhere between 6 and 15 of 17, and
consensus core residues from 0 to 29.

## Accomplishments that we're proud of

**The pipeline runs end to end, and seven of its top ten are known binders.** We co-folded
42 approved drugs against the target and scored 41 of them, ranking the whole list purely
by how much of the binding site each one engages. That works out to 2.87× enrichment at the
25% cut, against an arithmetic maximum of 4.10×, with the median known binder landing 8th
of 41.

**It surfaces the drugs chemistry throws away.** This is the result we set out to get:

| drug | by binding site | by chemical similarity |
|---|---|---|
| **pazopanib** | **4** | 764 |
| **sunitinib** | **5** | 464 |
| **nintedanib** | **9** | 755 |
| **vandetanib** | **10** | 1,760 |

Structural ranking finds all four because it never looks at what they are made of. It only
cares about where they sit.

**And it beats the chemical baseline, measured properly.** We re-ranked the same 41 drugs
both ways so the two populations match exactly. Structural ranking wins on median rank, 8.0
against 12.5, and it wins on AUC against property-matched decoys, 0.913 at p=0.0003 versus
0.795 at p=0.011. That is a real win of roughly 1.5×.

**The front of the pipeline picks its own target now, and it has opinions.** Give it a
disease and it ranks Open Targets associations, then gates each candidate on whether an
approved drug exists to validate against and whether a usable structure exists at all. On
venous thromboembolism it examined 967 associated proteins and returned 5 usable; on
colorectal cancer, 16,299 and 10; on rheumatoid arthritis it picked TYK2, which is where
that field actually is.

Two things it did that we did not ask for. It rejected **CDK2 twice over** — not in the
top 500 for colorectal cancer, and zero approved drugs with it as a mechanism — which is
the target we had spent the first half of the project on. And its structure filter threw
out 121 of 1,092 PDB entries, **every single rejection for the same reason**: a short
peptide already sitting in the binding site, which no ligand-stripping step removes. On
thrombin it rejected 37 of 64 entries, 10 of them at better resolution than the one it
kept — independently reproducing, from a rule, a judgement one of us had made by hand and
written in a code comment.

Where it stops short: the literature is fetched, every PMID re-checked against its title,
attached to the output — and still does not move the ranking. We narrowed that gap from
"fetched and unused" to "fetched, verified, attached and unused", which is honest progress
and not a solution.

**Fusing representations buys global ranking and costs you the top of the list.** We
combined ECFP4 with 3D shape and pharmacophore fingerprints by reciprocal rank fusion
across 2,114 drugs, scored with the same bootstrap CIs and Holm correction as everything
else. It beats ECFP4 on AUROC, 0.6839 against 0.6692, and on enrichment at 5%. It loses on
precision-at-1 and nDCG@10. All twelve comparisons are significant after correction, so
this is not noise — it is the same lesson a third time: pick one metric and you will hide
the disagreement.

**Chemistry does find one class of hit, and we pinned down exactly which.** Running the
inverse question blind — hide the answer, rank 2,382 approved drugs by structure alone
against VEGFR2 — brings the known VEGFR2 cancer drugs back at ranks 3, 6, 17 and 22, and
**mebendazole at rank 24**. Mebendazole is a 1974 deworming pill with no cancer
indication, and mebendazole-for-cancer is a real line of research currently in human
trials. Nothing in the ranking ever saw an annotation.

Co-folding it lands it in the same pocket as the approved drug, touching **19 of the same
19 amino acids**, overlap 0.950. Four unrelated approved drugs run as controls reach
0.41–0.64. Worth knowing which number carries that: the confidence score separates
nothing — paracetamol scores 0.970 against mebendazole's 0.990 — so it is the contact
overlap doing the work, not the model's confidence.

EGFR is the cleanest version. Query with one approved EGFR drug and the rest come back at
ranks **1, 2, 3, 4, 5, 7 and 12** of 2,382, 10.8× enriched in the top 50.

**So the two halves of this project cut in opposite directions, and both results are
real.** Chemistry finds mebendazole because mebendazole happens to look like axitinib. It
buries sunitinib at 464 because sunitinib does not. On EGFR the same thing happens to
chlorpromazine, which touches the target on unrelated chemistry and lands at rank 1,013.
Structural ranking catches precisely those. Neither method subsumes the other, and knowing
which one to reach for is most of the value.

**We also ran it where nothing is known, which is the harder honesty.** PADI4 drives
rheumatoid arthritis and has zero approved drugs, so there is no query molecule and no
positive control. Seeded with the arginine mimetic that published PAD4 programmes are
built on, it returns **pentamidine** and **hydroxystilbamidine** — both bis-amidines,
chemically the right class — alongside benzoic acid and phenol, which resemble a small
query only by being small. With no known answer, nothing in that run separates the
coherent hypothesis from the artifact. That is what the output looks like at the edge of
what the method can check, and we would rather show it than crop it.

**The target-side idea holds up under attack.** Across 1,531 structures covering 60 targets,
taking a consensus of the residues that other ligands engage predicts a held-out ligand's
contacts at 0.7322 precision, where a pocket finder manages only 0.4616. That is a gap of
+0.27 with a Wilcoxon p of 7.9e-113, and it still holds when you let every target count
exactly once so the well-studied ones can't dominate. Then one of us sat down and wrote a
15-finding adversarial critique of the whole thing before anybody was allowed to believe
it, and the result came through every blocking check intact: alternative rules for building
the set, three different contact cutoffs, resolution confounds, lipid exclusions.

**We decomposed our own win against a free null.** Predicting from a single randomly chosen
other ligand already gets you to 0.6193. So of that +0.27, roughly +0.16 is simply knowing
that ligands bind here at all, and the remaining +0.11 is what aggregating many binding
events actually buys you. That split told us far more about what the method is doing than
the headline number ever did.

**A benchmark with both a floor and a ceiling.** ECFP4 gets 66% of the way from a random
floor of 7% to a perfect ranking on p@1, and a 2048-bit hash of substructures beats a
77M-parameter transformer at this task. The interesting part wasn't the winner, though, it
was the disagreement: ECFP4 is the best at putting one correct drug first, while four other
representations beat it on the quality of the whole list. Had we reported a single number,
it would have looked like a clean sweep and we would have hidden something real.

**Right-sizing the binder nearly tripled design confidence**, taking ipTM from 0.162 to
0.598. It is obvious in hindsight. A 60–90-mer has nowhere to go inside an ATP slot,
whereas an 8–16-mer sits comfortably in a protease groove.

**Every arm runs a positive control on the same data with the same metric**, for the simple
reason that a failed search and a broken instrument look identical from the outside.
Bivalirudin retrieves lepirudin at rank 1 of 36. Docking puts three real thrombin drugs
above the generated median in the same batch. Every control passes.

**Honesty as an engineering constraint rather than a disclaimer.** CI lints the vocabulary,
so the claims can't quietly drift. The gate artifact records "neither passed nor failed"
instead of claiming the stronger ablation we never actually ran. And `results/` holds
computed output only — nothing in it has ever been hand-edited.

## What we learned

**Run a power check before you report a null.** Our docking arm got written up as "docking
adds nothing" on the strength of 23 molecules. At 199 the very same effect is significant.
The result had been sitting there the whole time, the molecules cost about 0.1 credit each,
and we simply hadn't bought enough power to see it. An underpowered null is not a negative
result.

**How you read a result can invert it.** Read as a top-5 list, our chemical baseline looked
useless on thrombin, since the five nearest drugs to argatroban are all peptidomimetics with
no thrombin annotation at all. Read as ranks, the same baseline works fine: bivalirudin
comes back at rank 10 out of 4,098, where chance would have put it around 2,050. Reading
only the top five nearly produced a false negative about our own control.

## What's next for our project

1. **A genuinely pocket-conditioned generator.** Our molecule arm generated unconditionally
   and only brought the pocket in afterwards, as a docking filter, which is precisely why
   finding nothing there tells you so little. DiffSBDD is within reach, since it's MIT
   licensed with public checkpoints, and the one thing standing in the way needs about a
   ten-line shim.
2. **The filtered design ensemble, which is genuinely untested.** Not one of our 24 designs
   cleared even ipTM 0.5, so what we actually ran was an unfiltered consensus. Getting to
   ten filtered survivors would take somewhere around 1,188 credits and we had 500. Whether
   that would close the gap is not something our result answers, and we don't think it
   should be claimed either way.
3. **Compare the two arms.** We have two independent rankings of approved drugs for the same
   target sitting on disk, and we have never once put them side by side. Where they agree
   is more interesting than either one alone, and the data is already there waiting.
4. **Replicate on a second target, end to end.** CDK2 is prepared and has never been
   screened. We didn't trust the gate until it replicated, and the newer arms are still
   standing on one target each.
