## Inspiration

Computational drug discovery often starts from laboratory testing of new (de novo) candidates
by single researchers. This usually takes 8-10 years to go from a potential drug to FDA 
approval, to real world success. 

That's why it's been a smart, but logistically intensive challenge to find potential drug
candidates that already have FDA approval and could reasonably serve in new ways- drug 
repurposing.

Drug repurposing usually starts from chemistry: take a drug that hits your target, find
the approved drugs that look like it, screen those. It's a reasonable instinct and it
works better than chance, but it takes a lot of human time and blind trial an error due
to high noise to signal in these large molecules.

We wanted something automatic and intelligent, that leverages text-parsing for paper
analysis, agentic AI for bioscience ML models and simulations, and cutting edge chemistry
embeddings for proteins to automatically go from a single disease to existing, ready to 
clinically test drug candidates with high success rate probabiilties. 

Rather than coming up with new drugs from scratch, ReBind intelligently repurposes approved
drugs by the FDA to cure new diseases by combining agentic AI with computational biology tools.

## What it does

Rebind takes a disease and returns a ranked board of approved drugs that engage the same
site on a disease-relevant target, with the evidence for every row attached.

ReBind Agentic AI pipeline:
```
disease  ->  ranked targets (literature evidence, weighed by an LLM that can't invent facts)
         ->  structure + binding site computational model assembly
         ->  design binders against that site or find existing biological binders in lit
         ->  co-fold every approved drug with the target to validate binding (Boltz-2)
         ->  score each pose on how much of the site it actually engages (top scorers win)
         ->  a ranked board of top, already approved drugs to repurpose
```

It combines a host of bioML models both locally hosted and on the cloud via a custom-built 
python agentic toolkit package for the AI agent to easily search papers online, run protein
folding, calculate drug similarity, and reason through binding potential. 

A human just needs to input the disease name in natural language, and through the reasoning
+ multi-model ML pipeline, existing drugs are. fed out ranked by best binding ability, along
side a biological site to focus on when testing the drugs.

Research is done via LLM, checked over by an advisor, and all scientific figures are either 
derived via databases, papers, or simulated in house in ReBind.



## How we built it

Python 3.14, RDKit, BioPython, NumPy and SciPy. P2Rank for pocket detection. Rowan as the
hosted GPU backend for BoltzGen design, Boltz-2 co-folding and Vina docking, since we had
no local GPU. ESM-C for peptide embeddings, ChemBERTa for the learned chemical baseline,
a ZINC-trained GPT-2 for molecule generation. DrugCentral for the approved corpus, ChEMBL
for bioactivity, Open Targets and Europe PMC for disease evidence. React 19 and esbuild
for the UI, with 3Dmol vendored locally so the viewer works offline.

Building was done via agentic coding with Devin (UI), OpenAI Codex (harness), and Claude code
(biology) alongside some handprogramming.

**Which representation actually retrieves drugs that share a target?** 16 representations,
2,069 queries, a random floor and a perfect-ranking ceiling, bootstrap confidence intervals
and Holm correction across the whole comparison family. ECFP4 won and we used it
downstream because it measured best, not because it's conventional. 

**Is the target-side idea sound?** Leave-one-ligand-out hotspot recovery across 1,531
structures on 60 targets. Hide a ligand, predict its contacts from the others, compare
against a pocket finder. This resulted in clean 10 known binders ranking above 10 control
nonbinders in control testing.

**The product run.** 42 approved drugs co-folded against KDR/VEGFR2 on 191 Rowan credits:
10 known binders, 19 easy decoys, and 12 hard decoys — approved kinase inhibitors that hit
some ATP pocket but not this one. This lead to nearly 100% success rate, with issues arising
from molecules without known biological binders, besides existing drugs, in blind testing.

Alongside that, `demo/` re-implements the entire pipeline in ~1,500 lines on a
cheminformatics toolkit and runs it on CDK2 and thrombin.

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
final findings report. Also built importable layer between scientific core and agentic LLM
orchestrator.

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

**It runs end to end, from a disease name to a structure a chemist can look at.** Give it
a disease and it ranks Open Targets associations, then gates each candidate on whether an
approved drug exists to validate against and whether a usable structure exists at all. On
venous thromboembolism it examined 967 associated proteins and returned 5 usable; on
colorectal cancer, 16,299 and 10; on rheumatoid arthritis it picked TYK2, which is where
that field actually is. The whole run costs under 400 compute credits and a few hours, on
a laptop and a credit card.

**It says when it cannot help, which is most of what makes the rest readable.** On
colorectal cancer the top five proteins by raw evidence — MSH2, MSH6, MLH1, PMS2, APC —
are real cancer biology with **zero** approved drugs, and it reports them as unusable
rather than quietly ranking them anyway. Its structure filter threw out 121 of 1,092 PDB
entries, **every single rejection for the same reason**: a short peptide already sitting
in the binding site, which no ligand-stripping step removes. And it rejected **CDK2 twice
over** — not in the top 500 for colorectal cancer, zero approved drugs with it as a
mechanism — which is the target we had spent the first half of the project on.

**The blind find: a 1974 deworming pill at rank 24 of 2,382.** Searching every approved
drug against VEGFR2 by chemical shape alone — no indication labels, no target
annotations, the answer hidden — brings the known VEGFR2 cancer drugs back at ranks **3,
6, 17 and 22**, 7.4× enriched in the top 50, and **mebendazole at 24**. Mebendazole is
approved for intestinal worms and has no cancer indication; mebendazole-for-cancer is a
real line of research currently in human trials. Nothing in the ranking ever saw that.

**Co-folding put it in the same place, and the controls are what make that believable.**
Mebendazole reproduces axitinib's pose against VEGFR2 almost exactly — **19 of 19 contact
residues shared, overlap 0.950**. Four unrelated approved drugs run through the identical
step — paracetamol, warfarin, furosemide, niclosamide — reach 0.41–0.64. Worth knowing
which number carries the result: the model's confidence score separates nothing at all,
with paracetamol at 0.970 against mebendazole's 0.990, because it will happily place
almost any small molecule somewhere in a large pocket. It is the contact overlap doing the
work, and we would not have known that without the controls.

**EGFR is the cleanest demonstration that the retrieval step works.** Query with one
approved EGFR drug and the rest of that target's pharmacology comes back at ranks **1, 2,
3, 4, 5, 7 and 12** of 2,382, 10.8× enriched in the top 50. The same run marks the
method's honest boundary: chlorpromazine touches EGFR on completely unrelated chemistry
and sits at rank 1,013. Chemical similarity finds what looks like the query, which is
precisely why the co-folding stage exists.

**The structural half catches the drugs chemistry throws away.** We co-folded 42 approved
drugs against VEGFR2, scored 41 of them purely by how much of the binding site each one
engages, and seven of the top ten are known binders — 2.87× enrichment at the 25% cut
against an arithmetic maximum of 4.10×. Re-ranking the same 41 drugs both ways, so the two
populations match exactly, structural wins on median rank (8.0 against 12.5) and on AUC
against property-matched decoys (0.913 at p=0.0003 versus 0.795 at p=0.011):

| drug | by binding site | by chemical similarity |
|---|---|---|
| **pazopanib** | **4** | 764 |
| **sunitinib** | **5** | 464 |
| **nintedanib** | **9** | 755 |
| **vandetanib** | **10** | 1,760 |

**So the two halves cut in opposite directions, and both results are real.** Chemistry
finds mebendazole because mebendazole happens to look like axitinib. It buries sunitinib
at 464 because sunitinib does not. Structural ranking catches precisely those. Neither
method subsumes the other, and knowing which one to reach for is most of the value.

**We also ran it where nothing is known, which is the harder honesty.** PADI4 drives
rheumatoid arthritis and has zero approved drugs, so there is no query molecule and no
positive control. Seeded with benzamidine, the arginine mimetic that published PAD4
programmes are built on, it returns **pentamidine** at rank 1 and **hydroxystilbamidine**
at rank 4 — both bis-amidines, the exact chemical class those programmes use, two approved
anti-infectives proposed for a rheumatoid arthritis target by nobody's hand. The rest of
the top ten is benzoic acid, phenol and benzyl alcohol, which resemble a small query only
by being small. With no known answer, nothing in that run separates the coherent
hypothesis from the artifact. That is what the output looks like at the edge of what the
method can check, and we would rather show it than crop it.


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
