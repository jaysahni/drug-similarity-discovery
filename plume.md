# Rebind — HackMIT 2026 submission writeup

Team: Evan Xiang, Jay Sahni, Richard Shan, Rachel Chen.
Repo: `drug-similarity-discovery`. MIT licensed. Every number below was computed by a
script in the repo and can be reproduced by running it.

---

## Inspiration

Drug repurposing usually starts from chemistry. Take a drug known to hit a target, find
the approved drugs that look most like it, screen those. We wanted to check that premise
before building on it, so we measured it first.

It does not hold up. We ranked the approved-drug corpus by ECFP4 similarity to axitinib,
the best literature ligand for VEGFR2, and asked where the 32 approved drugs annotated
against that target actually land. Similarity is genuinely enriched — 10.2× over base
rate — and it still buries sunitinib at rank 464, pazopanib at 764 and vandetanib at
1,760. Those are real, approved VEGFR2 inhibitors. A chemical-similarity shortlist would
never have surfaced them.

The reason is structural rather than statistical: two drugs that share a binding site
need not share any chemistry. So we moved the question to the other side of the
interaction. Instead of asking what a molecule looks like, ask which residues it engages.
That description is comparable across molecules with no chemistry in common — even
between a designed miniprotein and a small molecule, which is what made the generative
arm thinkable at all.

## What it does

Rebind takes a disease and returns a ranked board of approved drugs predicted to engage
the same site on a disease-relevant target, with the evidence for each row attached.

The pipeline runs in target-side coordinates end to end:

```
disease -> ranked targets (Open Targets + Europe PMC + PubMed, LLM weighs verified evidence)
        -> structure + pocket (P2Rank)
        -> BoltzGen designs binders -> consensus INTERFACE SIGNATURE
        -> co-fold each approved drug WITH the target, unconstrained (Boltz-2)
        -> score the pose on how much of that site's core it engages
        -> ranked board + report
```

The score is `precision_in_core`: of the residues a drug engages, the fraction lying in
the site core, at a 4.5 Å heavy-atom cutoff in author numbering. Every row names the core
residues engaged and missed, so a hit can be inspected and argued with rather than taken
on faith.

<!--LANG-EXEMPT-->
**Nothing in the pipeline predicts binding strength.** Boltz-2 can emit an affinity
score; it is deliberately never read. Ranking is interface overlap only. This is enforced
rather than promised — `scripts/lint_language.py` runs in CI and fails the build on
affinity, Kd, IC50 or potency vocabulary outside explicitly marked exemptions.
<!--/LANG-EXEMPT-->

Three surfaces exist on top of the results:

- `scripts/ask.py` — an offline read-only front door. Ask a disease or a target, get the
  board back instantly, with every screenful naming the file each number came from. No
  network, no credits, no recomputation.
- **Rebind (web UI)** — a React and vanilla-JS workspace with chats, folders, saved
  results, a pipeline diagram and a 3Dmol structure viewer that highlights engaged versus
  missed core residues on the actual complex.
- A self-contained HTML report per run, plus a schema-validated JSON manifest that is the
  only integration boundary between the science and the visualization workstreams.

## How we built it

**Stack.** Python 3.14, RDKit, BioPython, NumPy/SciPy. P2Rank for pocket detection (needs
a JVM). Rowan as the hosted GPU backend for BoltzGen design, Boltz-2 co-folding and Vina
docking — we had no local GPU. ESM-C for peptide embeddings, ChemBERTa for the learned
chemical baseline, `entropy/gpt2_zinc_87m` for small-molecule generation. DrugCentral for
the approved corpus, ChEMBL for bioactivity, Open Targets and Europe PMC for disease
evidence. React 19 and esbuild for the UI, 3Dmol 2.4.2 vendored locally so the viewer
works offline.

**Four experiments, run in order, each gating the next.**

*E1 — does any representation beat a fingerprint?* 16 representations on a 2,069-query
retrieval task over 2,114 DrugCentral drugs, with a random floor, a perfect-ranking
ceiling, bootstrap CIs and Holm correction over the 32-comparison family. ECFP4 wins p@1
at 0.6597 and three others tie it. ChemBERTa, a 77M-parameter transformer, is
significantly *worse* at 0.5607. But the metrics disagree: on AUROC, four representations
beat ECFP4. Reporting one number would have called this a clean win and hidden the
disagreement. Replicated on a stricter mechanism-of-action ground truth.

*E2 — is the target side better?* Leave-one-ligand-out hotspot recovery on 1,531
structures across 60 targets, 1,566 distinct ligands, 10,735 P2Rank pockets. A consensus
of which residues the other ligands engage predicts a held-out ligand's contacts at
precision 0.7322 against P2Rank's 0.4616, Δ +0.2706, Wilcoxon p=7.9e-113. It survived a
15-finding adversarial critique written against it before we believed it — alternative
set rules, three contact cutoffs, resolution confound, lipid exclusions.

*The M2 gate — the decisive ablation.* The whole project turns on one question: does a
BoltzGen design consensus beat a free pocket finder at saying where a ligand binds?
Ground truth is the residues the known ligand actually contacts in 38 KDR co-crystals,
all arms scored on the same 38 so comparisons are paired.

*E3 and the product run.* 42 approved drugs co-folded against KDR/VEGFR2 and ranked by
interface overlap, on 191 Rowan credits: 10 known binders, 19 easy decoys, and 12 **hard
decoys** — approved kinase inhibitors that DrugCentral does not annotate against KDR,
each of which binds *some* ATP pocket by construction.

**A parallel lightweight arm.** `demo/` re-implements the whole pipeline in ~1,500 lines
on the `cheminformatics-kit` toolkit rather than the main repo's 12,305-line `scripts/`,
and runs it on CDK2 and thrombin to ask a different question: given a designed binder,
which approved drug resembles it? Three arms, three negatives, each with a positive
control that says the negative is real.

## Individual Contributions

Attribution from git history across 77 commits on four branches.

**Evan Xiang** — the scientific core. The E1 retrieval benchmark and the 17-representation
zoo; the metrics layer (bootstrap CIs, paired tests, Holm); target-side contact extraction
and the `InterfaceSignature` type; the P2Rank baseline arm and its data build; E2 hotspot
recovery; getting BoltzGen to run at all and the targeting diagnostics that made it aim;
the M2 gate and its convergence study; the repurposing pipeline end to end; the hard-decoy
test that changed the project's verdict; the `ask.py` query interface; the CI language
lint; the role-separation significance tests.

**Jay Sahni** — the project spec (`PROJECT_GOAL.md` v0.2) and the entire lightweight demo.
The direct-matching harness and its calibrated null; arm P (peptide design and ESM-C
matching, with the hub check that explained the result); arm S (GPT-2 generation, Vina
docking, ECFP4 matching); the approved-biologics corpus and its builder; the insertion-code
fix that recovered thrombin's active site; the final findings report and retrospective;
the corpus-approval fix that corrected a headline downward.

**Richard Shan** — integration, evidence discipline and the agent step. The LLM judge
between evidence retrieval and structure prep, with the hallucination guard that rejects
any PMID or gene symbol not in the verified evidence set, negative-tested; the
visualization handoff bundle and wiring the ranked board and structures through to the
viewer; the decoy relabelling from measured ChEMBL activity; the GPU handoff and pipeline
docs; the retraction of the 100× comparison and of a batch timing quoted as a single-run
cost; the cross-branch merges that kept four workstreams on one main.

**Rachel Chen** — the Rebind UI. The React pipeline diagram, and the workspace layer:
chats, folders, saved results and compact navigation.

## Challenges we ran into

**The decisive experiment failed.** The M2 gate says BoltzGen loses. A sub-second-per-
structure pocket prediction reaches Jaccard 0.6355 against held-out ligand contacts; a
24-design, 85-credit, 28-minute-of-A100 BoltzGen consensus reaches 0.3115. Δ −0.324, Holm
p=0.0016, n=38 paired. It loses on precision, recall and F1 too, and it replicated on
CDK2 (0.338 vs 0.494). P2Rank lands within 0.03 Jaccard of the known-ligand ceiling —
there is almost no headroom above pocket geometry on this target for anything to occupy.
Most of the project's second half is the work of establishing that carefully enough to
report it as a result rather than a failure.

**Six bugs that produced ordinary-looking wrong answers.** None threw an error. A PDB
parser that ignored insertion codes silently lost 28 of thrombin's 259 residues, including
the entire 60-loop lining the active site. Structure preparation renamed chains (H→B,
L→A), so every downstream stage pointed at nothing. A co-crystal ligand left in the pocket
during design sent every design to the wrong site — Jaccard 0.000 against known ligand
contacts, 52 credits spent. Design structures cached without a UUID meant a re-run scored
against the previous run's complexes. Residues mapped by a constant offset were correct by
accident on one run. And a percentile used a strict `<`, charging every exact Tanimoto tie
against the candidate and biasing every headline number low, invisibly — caught only by a
calibration test (KS D=0.187, p<0.001). All six share one root cause: trusting an
identifier that an upstream tool is free to change. All six are now guarded by tests.

**A hardware and packaging floor we could not raise.** No local GPU. Modal authenticates
and runs CPU, but all four GPU tiers are refused pending a payment method. The Rowan free
tier is 500 credits plus 20 a week. `boltz` cannot install on Python 3.14; `vina` 1.2.7
publishes no cp314 wheel and its sdist fails to build; every pocket-conditioned generator
we tried (DiffSBDD, Pocket2Mol, TargetDiff, PocketFlow, DecompDiff, LiGAN) depends on
`torch-scatter`/`torch-sparse`/`torch-cluster`, which have no build for torch 2.14. So
the generative arm was an unconditioned generator plus a docking filter, which is a
weaker thing than we set out to build, and we said so.

**BoltzGen does not do what its interface suggests.** It indexes residues 1..N over
residues present in the file, ignoring author numbering — our 303-residue construct is
authored 811–1169, so every pocket residue was out of range and the first run died
immediately. And naming hotspot residues via `binding_types` does not aim the design at
all: 0, 2, 0, 0 of 17 designs hit the requested pocket, with an *empty* consensus.
Restricting the presented surface with a 12 Å `include_proximity` is what worked — 6, 12,
12, 15 of 17, and consensus core residues went from 0 to 29.

**A headline that was wrong twice over.** An earlier README claimed known binders rank
"4, 5, 9, 10 structurally versus 464, 755, 764, 1760 chemically". That compared a rank out
of a 41-drug board against a rank out of a 4,099-drug corpus — the ratio 4099/41 ≈ 100 was
almost exactly the apparent gap, so the comparison was measuring population sizes. It also
cherry-picked the four drugs with the largest gap; over all 10, paired Wilcoxon gives
p=0.85. We retracted it in place rather than quietly dropping it, and re-ran the honest
comparison on the same 41 drugs: structural still wins, by about 1.5×, not 100×.

## Accomplishments that we're proud of

**The pipeline runs end to end and it works.** 42 drugs co-folded, 41 scored, ranked by
interface overlap: enrichment 2.87× at 25% against an arithmetic maximum of 4.10×, median
known-binder rank 8.0 of 41, 7 of the top 10 are known binders. Pazopanib ranks 4,
sunitinib 5, nintedanib 9, vandetanib 10 — the same drugs chemical similarity buried at
764, 464, 755 and 1,760.

**And it beats the chemical baseline like-for-like.** Re-ranking the same 41 drugs by
chemical similarity: AUC 0.913 (p=0.0003) structural versus 0.795 (p=0.011) chemical
against property-matched decoys, median rank 8.0 versus 12.5.

**We built the test that broke our own result, and then reported it.** The hard decoys —
approved kinase inhibitors that miss KDR — are the only null that is actually hard, and
against them neither arm is significant: AUC 0.717 (p=0.092) for the pocket signature and
0.529, a coin flip, for the design signature. Enrichment falls from 2.87× to 1.83×. The
headline number is real but it is measuring the easier task, and that sentence is in the
README above the headline.

**Free baselines were tested against us, including the one that wins.** Every co-folded
pose comes back with Boltz-2 confidence fields we did not ask for. On the easy null, ipTM
alone reaches AUC 0.995 against our score's 0.913. That is reported because it is true.
Two things survive it: on the hard null nothing works, including ipTM, which locates the
ceiling in the co-folding backend rather than in our scorer; and interface overlap names
the residues engaged and missed, while ipTM is one number with no mechanism attached.

**Every negative has a positive control on the same data with the same metric.** A failed
search and a broken instrument look identical from outside. bivalirudin retrieves
lepirudin at rank 1 of 36, so ESM-C works and the designs simply fail it. Vina ranks three
real thrombin drugs above the generated median in the same batch. 53% of randomly shuffled
sequences land on abarelix, which is why every design landing there carries no
information. Six controls, all passing, and one that nearly reported itself as broken.

**Honesty as an engineering constraint, not a disclaimer.** CI lints the vocabulary. The
M2 gate artifact records `m2_gate_status: "neither passed nor failed"` rather than
claiming the stronger ablation we did not run. Drugs whose pose produced no drug-like
ligand keep their row with a status and no score, never a plausible-looking number.
`results/` is computed output only, never hand-edited.

## What we learned

**Run a power check before reporting any null.** Our docking arm was written up as
"docking adds nothing" at n=23 (ρ=−0.084, p=0.703). At n=199 the same effect is
significant (ρ=−0.182, p=0.0099; ρ=−0.146, p=0.040 after correcting the corpus to the
2,153 actually-approved structures). The result existed the whole time. Docking cost
about 0.1 credit per molecule — the power was affordable and we simply did not buy it.
An underpowered null is not a negative result.

**Controls are what make a negative readable.** Four of our findings are negative. Every
one is publishable only because a control on the same data with the same metric passes.
Without them we would have had four shrugs.

**How you read a result can invert it.** Read as a top-5 list, ECFP4 looked useless for
thrombin: argatroban's five nearest approved drugs are all peptidomimetics with no
thrombin annotation. Read as ranks, it works — argatroban retrieves bivalirudin at rank 10
of 4,098 against ~2,050 expected by chance. The top-five reading nearly produced a false
negative about our own positive control.

**Aggregation buys less than it looks like it does.** E2's consensus beats a pocket finder
by +0.27 precision. But a null that costs nothing — predict from one randomly chosen other
ligand — already scores 0.6193. So +0.16 of that gap is just "knowing that ligands bind
here at all", and only +0.11 is what aggregating many binding events actually buys.
Decomposing the win against a free null changed what we believed the method was doing.

**Re-derive identity from content, never from an identifier upstream can change.** One
root cause, six separate incidents, six separate fixes. We should have made it a stated
invariant with a test the first time.

**Size the tool to the site.** Switching BoltzGen from `protein-anything` at 60–90
residues to `peptide-anything` at 8–16 raised interface confidence from ipTM 0.162 to
0.598. A 60–90-mer cannot enter an ATP slot; an 8–16-mer fits a protease groove. And
pocket detection cannot find a serine protease active site at all — six pockets on
thrombin, not one containing Ser195 or Asp102, because a protease active site is a shallow
groove across subsites rather than an enclosed cavity.

**Establish the question before building the pipeline.** We built footprint matching —
the designed binder as a discarded probe — before establishing that the intent was direct
matching, the designed binder as the query. One question at the outset would have saved
roughly half the working time.

## What's next for our project

1. **A genuinely pocket-conditioned generator.** Our small-molecule arm generated
   unconditionally and used the pocket only as a docking filter, which is why finding
   nothing says little. DiffSBDD is reachable — MIT licensed, Zenodo checkpoints — and its
   only blocker, `torch-scatter`, needs roughly a 10-line shim. REINVENT4 with Rowan
   docking as the RL reward is the alternative.
2. **The filtered design ensemble, which is genuinely untested.** Zero of 24 designs
   reached even ipTM 0.5 against the plan's 0.85 gate, so what we ran is an *unfiltered*
   consensus. Reaching 10 filtered survivors at the published rate needs ~333 designs ≈
   1,188 credits against a 500-credit tier. Whether a properly filtered ensemble closes a
   0.32 Jaccard gap is not answered by our result and should not be claimed either way.
3. **Compare the two arms.** Two independent rankings of approved drugs for the same
   target now exist on disk and have never been compared. Where they agree is more
   interesting than either alone, and the data is already there.
4. **Replicate on a second target end to end.** CDK2 is prepared and never screened. The
   M2 gate was not trusted until it replicated; the new arms are still n=1 target.
5. **Make the report target-agnostic.** `report.py` inputs are hardcoded to KDR.
6. **Correct the spec.** `PROJECT_GOAL.md` §4.3's ipTM ≥ 0.85 gate has no BoltzGen
   provenance — the paper publishes no ipTM distribution at all, and published miniprotein
   work filters at > 0.5. The plan document specifies a filter this repo's own
   investigation showed to be unreachable by construction.
7. **Golden fixtures**, which the spec already mandates and whose absence is why a
   statistics block silently vanished the first time a stage re-ran.

---

*Rebind ranks structural hypotheses about shared site engagement. Nothing here predicts
clinical benefit or binding strength, and every output requires experimental validation.
Rows labelled "known binder" are the benchmark's positive controls, not discoveries.*
