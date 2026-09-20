# 15 — Autoresearch: disease → target, the stage the demo never had

`demo/pipeline.py` calls its first stage "research". It does three things:

```python
TARGET = {"symbol": "CDK2", "uniprot": "P24941", "pdb_id": "6Q4G", ...}  # a dict literal
entry = nova.tool_data("database.fetch_uniprot_entry", ...)   # fetch its sequence
lit   = nova.tool_data("literature.search_pubmed", ...)       # 10 abstracts, stored, unused
```

The target is a literal a human typed. The literature is fetched and never read
by anything downstream. There is no disease anywhere in the stage. The original
ask was *autoresearch → drug candidate → folding → match against approved FDA
drugs*, and the first word of that was missing.

This note records the stage that fills it: `demo/autoresearch.py`, which takes a
disease and returns a ranked list of targets with the evidence for each, a PDB
id and chain, and — the part that decides whether the rest of the pipeline is
readable at all — a count of approved drugs that already bind the target.

Builder: `demo/autoresearch.py`. Self-test: `demo/test_autoresearch.py` (32
checks, all passing). Output: `results/demo/autoresearch/<disease-slug>.json`.
Every number below is printed by the script on each run; none is estimated.

Cost: zero. Open Targets, RCSB, Europe PMC and UniProt are open APIs with no
credentials. No Rowan credits are spent at this stage. Every response is cached
under `data/raw/` (gitignored), so the second run is offline.

---

## 1. What was reused, and what had to be written

Prior sessions built most of this and it was never wired into the demo. The
reuse is by import, not by copy, so the caches and the verification rules are
shared rather than forked.

| Source | Reused | How |
|---|---|---|
| `scripts/disease_targets.py` (401 ln) | `resolve_disease`, `associated_targets`, `build_target`, `gql`, `api_meta`, `drug_candidates`, `uniprot_rest` | imported. Free-text → EFO/MONDO id, associated targets, per-datatype scores, UniProt accession, Open Targets' own drug counts |
| `scripts/target_evidence.py` (634 ln) | `search_target`, `verify`, `aliases`, `pick_claim`, `epmc` | imported. Europe PMC search, verbatim claim sentence with its tier, every PMID re-fetched and title-compared |
| `scripts/prep_target.py` (932 ln) | `http`, `cached`, `NOT_LIGANDS`, `MIN_LIGAND_MW`, `RCSB_SEARCH`, `RCSB_GRAPHQL` | imported. Retrying HTTP with an on-disk cache; the buffer-vs-ligand rule, so "drug-like ligand" means the same thing here as in the heavy pipeline |

Written fresh (`demo/autoresearch.py`, 1,009 ln) — all of it selection logic,
none of it re-implementation:

- **`approved_binders()`** — the positive-control count, from this repo's own
  three corpora rather than from Open Targets' drug counts (§2).
- **`search_entities_by_resolution()` / `entry_records()` / `choose_structure()`**
  — structure selection with the co-crystallised-peptide guard (§4). The RCSB
  GraphQL query is extended to the entry's *other* polymer entities, which
  `prep_target.py` never needed.
- **`uniprot_sites()`** — catalytic-nucleophile anchors from UniProt `ACT_SITE`
  (§5).
- **`circularity()`** — the known-drug share of the datatype evidence (§3).
- **`site_and_modalities()`**, **`pipeline_config()`** — the
  `TargetHypothesis` fields PROJECT_GOAL §4.1 adds (`site_type_hint`,
  `preferred_modalities`) and the emitter for `demo/pipeline.py`'s `TARGETS`
  dict.
- **`target_extras()`** — Open Targets `tractability` / `safetyLiabilities` /
  `targetClass`, which `disease_targets.py` does not fetch. The field is
  `tractability { label modality value }`; `{ id ... }` returns HTTP 400,
  established by introspecting the live schema rather than guessing.

---

## 2. The gate that matters: does an approved drug already bind it?

A target with no approved binder produces a ranking of 2,153 approved drugs with
no known-correct answer in it. Every percentile the matching arms report is then
uninterpretable — not wrong, *unreadable*. This is the whole difference between
thrombin and CDK2, and it is now a gate rather than a hunch.

The count comes from the same files the arms are scored against, not from Open
Targets:

| Corpus | n | What is counted |
|---|---:|---|
| `data/approved_drugs.csv` | 4,099 rows → **2,153** with `approved == 1` | `moa_targets` (mechanism-of-action) and `targets` (any annotation), split on `;` |
| `data/target_annotations.csv` | 28,337 rows → 20,744 joined to an approved drug and `passes_threshold == 1` | pChEMBL-backed binding, matched on UniProt accession |
| `data/approved_biologics.csv` | 262 | `target_gene_symbol` / `target_uniprot`, with modality, for the peptide arm |

Three counts are kept separate, because collapsing them is exactly the mistake
that would let CDK2 through:

```
CDK2   n_approved_binders_moa          = 0
       n_approved_binders_annotated    = 7   (ceritinib, lapatinib, nintedanib,
                                              palbociclib, ribociclib, sorafenib,
                                              trilaciclib)
       n_approved_binders_bioactivity  = 9   pChEMBL-backed
       n_approved_biologic_binders     = 0

F2     n_approved_binders_moa          = 4   argatroban, bivalirudin,
                                              dabigatran etexilate, ximelagatran
       n_approved_binders_annotated    = 12
       n_approved_binders_bioactivity  = 13
       n_approved_biologic_binders     = 2   bivalirudin (peptide), lepirudin
```

CDK2 has seven approved drugs carrying a CDK2 annotation and nine with measured
potency against it. **None of them is a CDK2 drug** — they are CDK4/6 and
multi-kinase inhibitors picking CDK2 up off-target. A gate counting annotations
waves CDK2 through; a gate counting mechanism-of-action targets does not. The
self-test asserts both halves of that, because the failure mode is a gate that
silently gets looser.

Gate as implemented: **`n_approved_binders_moa >= 1`**.

Its cost is recorded rather than hidden: it fails 14/25 targets for colorectal
cancer and 20/25 for venous thromboembolism, including the top-ranked
association in both cases (MSH2, PROC). Those are kept in the output with
`unusable_reasons` spelling out the count that failed.

---

## 3. And the gate that undermines it: circularity

PROJECT_GOAL §1.3: Open Targets target–disease evidence already contains ChEMBL
known-drug data. A target ranked partly *because* drugs exist for it, then
"discovered" to have those drugs, is a database lookup in a lab coat. §2's gate
selects **for** targets with approved drugs, which is to say it selects for
exactly the evidence that is circular. The two criteria pull against each other
and the honest move is to show both, not to pick one.

So every target reports its full `score_by_datatype` plus:

- `known_drug_score` — the known-drug datatype score. In data release 26.06 the
  datatype id is **`clinical`**, not `known_drug`; `disease_targets.py` probes
  both spellings and records which the live API returned, so a rename does not
  silently turn the flag off.
- `known_drug_share_of_datatype_sum` — the known-drug score over the sum of all
  datatype scores. **This is a descriptive statistic, not a decomposition.**
  Open Targets aggregates datatypes by a weighted harmonic sum, so these do not
  add to `association_score` and must not be read as if they did. It answers one
  question: among the lines of evidence, how big is the known-drug one.
- `known_drug_dominated` — true when the known-drug datatype is the target's
  single largest datatype **and** its share is ≥ 0.25.

Worked examples from the runs:

| Target | Disease | `clinical` | Largest datatype | share | dominated |
|---|---|---:|---|---:|---|
| MSH2 | colorectal | — (absent) | `genetic_association` 0.975 | — | no |
| BRAF | colorectal | 0.952 | `literature` 0.996 | 0.189 | no |
| KDR | colorectal | 0.971 | `clinical` | 0.372 | **yes** |
| F2 | VTE | 0.901 | `clinical` | 0.441 | **yes** |

Result: 2/10 usable colorectal targets and 3/5 usable VTE targets are
known-drug dominated. It is reported next to the pick, never used to filter —
dropping them would delete the only targets with a positive control.

---

## 4. Structure selection, and the mistake this repo already made

Downstream needs a PDB id and a chain, and `demo/pipeline.py` strips
heteroatoms before design. **A co-crystallised peptide is not a heteroatom.** It
is ATOM records in a polymer chain, `remove_heterogens` leaves it exactly where
it is, and a binder designed against that entry is designed against an occupied
site. `demo/pipeline.py`'s thrombin comment records rejecting 4UD9 / 4UE7 / 5AFY
by hand for precisely this. That judgement is now code.

**Rule.** Candidates are X-ray polymer entities mapping to the target's UniProt
accession.

- **Hard reject** — the entry contains a `polypeptide(L)` entity that does *not*
  map to the target accession and is ≤ 50 aa. Thrombin's own 28-aa light chain
  maps to P00734 and therefore does not trigger it; hirudin (P09945, 12 aa)
  does. The self-test asserts both.
- **Hard reject** — no auth chain for the accession, or no reported resolution.
- **Soft penalty** — a larger foreign polymer (a Fab, a partner protein, a
  nucleic acid): flagged as `ppi_interface` evidence and deprioritised, not
  dropped, because the pocket of interest may be elsewhere.
- **Order among survivors** — holo before apo, no foreign polymer before one,
  then best resolution, then PDB id for determinism.
- **Chain** — the auth chain of the *longest* entity mapping to the accession
  (thrombin's 259-aa heavy chain, not its 36-aa light chain).

**Candidate generation is two searches, not one.** Sorting P00734's 900 entities
by resolution and taking the best 60 gives a slice that is almost entirely
hirudin co-complexes; the guard rejects 29 of 30 and the one survivor is apo. So
the first search is biased toward entries that *can* be holo and cannot be a
large complex (`nonpolymer_entity_count > 0`, at most two protein entities —
a bar of one would be wrong for thrombin, whose heavy and light chains are both
P00734), and the unrestricted search follows as the fallback. The bias only
orders candidates; `choose_structure()` decides.

**What it does on thrombin** (P00734, 900 entities, 120 examined, 64 entries):

```
rejected 37/64 entries, all for co_crystallised_peptide_polymer_chain
10 of the rejects are at BETTER resolution than the pick:
  5AFY 1.12  4UD9 1.124  4UE7 1.129  4UEH 1.16  4UDW 1.16
  5AF9 1.18  3RM2 1.23   5AHG 1.24   6FJT 1.27  5MM6 1.29
chosen: 1SHH chain B, 1.55 A, holo, ligand 0G6 (PPACK), 259 aa
```

1SHH carries the *same* ligand as the hand-picked 1PPB (0G6 = PPACK) at 1.55 Å
instead of 1.92 Å. The automatic pick is the hand pick's construct, one entry
better resolved.

**What it does on CDK2** (P24941, 506 entities): **6Q4G chain A, 0.98 Å, holo** —
byte-identical to the hand-written `TARGETS["cdk2"]` entry.

Across both disease runs: 1,092 entries examined, 121 rejected, the peptide
guard firing on 121 of them across 18 targets.

---

## 5. Catalytic anchors, derived rather than typed

`demo/pipeline.py` anchors thrombin's site on `[568]` — Ser195, the catalytic
nucleophile — with a hand-written comment explaining that His57 and Asp102 sit
*behind* the serine, so a sphere over the whole triad reaches into the core
rather than the ligand site.

`uniprot_sites()` reproduces that without being told: it keeps a UniProt
`ACT_SITE` position only when the residue at that position is S, C or T.

```
P00734  ACT_SITE {406: H, 462: D, 568: S}  ->  site_anchors [568]
P24941  ACT_SITE {127: D}                  ->  site_anchors []
```

CDK2's only active site is the aspartate proton acceptor, so it gets no anchor —
matching `TARGETS["cdk2"]`, which has none. The same function returns UniProt
`CHAIN` features, recovering thrombin heavy chain 364–622.

---

## 6. Results

Two diseases, top 25 associated targets each, `--lit 2`.

### 6.1 Venous thromboembolism — MONDO_0005399, 967 associated targets

```
 # symbol   uniprot    assoc   kd%  moa bio pep  pdb  ch   res site          verdict
 1 PROC     P04070    0.7417   37%    0   0   0  1AUT C   2.80 pocket        no control
 2 F10      P00742    0.7258   48%    4   0   0  2JKH A   1.25 pocket        USABLE (known-drug dominated)
 3 F2       P00734    0.7084   44%    4   2   1  1SHH B   1.55 pocket        USABLE (known-drug dominated)
 4 PROS1    P07225    0.6847     -    0   0   0  -    -      - unknown       no control
 5 SERPINC1 P01008    0.6612   63%    1   0   0  1E04 I   2.60 unknown       USABLE (known-drug dominated)
 6 F11      P03951    0.6442   32%    0   0   0  7MBO A   0.92 pocket        no control
14 CPS1     P31327    0.5188     -    1   0   0  6UEL A   1.90 pocket        USABLE
17 KLKB1    P03952    0.5007     -    2   2   0  6I44 A   1.36 pocket        USABLE
```
*(rows 7–13, 15–16, 18–25 all "no control"; full table in the JSON)*

n = 25 examined, **5 usable**; 20 fail the positive-control gate, 3 fail the
structure gate (PROS1, C4BPB, SH2B3 — all three fail both, so the union is 20:
no X-ray entity maps to any of those accessions). 22/25 got a
structure, 17 of them holo. 19/25 targets returned ≥1 Europe PMC paper;
**29 PMIDs, all title-verified, 0 dropped**.

### 6.2 Colorectal cancer — MONDO_0005575, 16,299 associated targets

```
 # symbol   uniprot    assoc   kd%  moa bio pep  pdb  ch   res site          verdict
 1 MSH2     P43246    0.9239     -    0   0   0  8RB1 A   2.08 pocket        no control
 2 MSH6     P52701    0.9157     -    0   0   0  2O8B B   2.75 pocket        no control
 3 MLH1     P40692    0.9142     -    0   0   0  4P7A A   2.30 pocket        no control
 4 PMS2     P54278    0.8825     -    0   0   0  9S89 A   2.00 pocket        no control
 5 APC      P25054    0.8585     -    0   0   0  3NMW A   1.60 unknown       no control
 6 BRAF     P15056    0.8385   19%    5   0   0  8JNB B   1.62 pocket        USABLE
 8 PIK3CA   P42336    0.8178    5%    2   0   0  9CMK A   1.75 pocket        USABLE
 9 KRAS     P01116    0.8093   15%    1   0   0  9IAY A   0.95 pocket        USABLE
15 EGFR     P00533    0.7732   29%   10   2   0  8A27 A   1.07 pocket        USABLE
16 FGFR2    P21802    0.7719   23%    3   1   0  10OU A   1.77 pocket        USABLE
17 POLE     Q07864    0.7718   11%    2   0   0  5VBN B   2.35 ppi_interface USABLE
20 POLD1    P28340    0.7561   14%    2   0   0  -    -      - unknown       no structure
21 FLT4     P35916    0.7524   40%    7   0   0  4BSJ A   2.50 unknown       USABLE (known-drug dominated)
22 ERBB2    P04626    0.7495   20%    5   3   0  8VB5 A   1.48 pocket        USABLE
23 KDR      P35968    0.7495   37%   12   1   0  2XIR A   1.50 pocket        USABLE (known-drug dominated)
24 FGFR3    P22607    0.7483   29%    4   0   0  9KFU A   1.40 pocket        USABLE
```

n = 25 examined, **10 usable**; 14 fail the positive-control gate, 2 fail the
structure gate (ATM, POLD1). 23/25 got a structure, 17 holo. 25/25 targets
returned ≥1 paper; **49 PMIDs, all title-verified, 0 dropped**.

**The top five associations are all unusable.** MSH2, MSH6, MLH1, PMS2 are the
mismatch-repair genes and APC is a tumour suppressor — real colorectal biology,
zero approved drugs binding any of them, and the arms would have nothing to
score against. The stage says so and keeps them in the output; it does not
quietly slide down to rank 6.

---

## 7. Would the pipeline have picked CDK2 or thrombin on its own?

Asked directly, because the hand-picked targets are the thing this stage is
meant to replace.

**Thrombin: yes, at rank 2 of the usable list.** Run on venous thromboembolism,
F2 is association rank 3 and pipeline rank 2, behind Factor Xa (F10) — the other
approved-anticoagulant target, which the stage found on the same evidence and a
human had not considered. The emitted config:

```python
"f2": {
    "symbol": "F2", "uniprot": "P00734",
    "pdb_id": "1SHH", "chain": "B", "out": "f2",
    "protocol": "peptide-anything", "binder_length": "8..16",
    "site_anchors": [568],
}
```

versus the hand-written `TARGETS["thrombin"]`: identical on symbol, uniprot,
protocol, binder length and site anchors; different only in the PDB entry, where
it chose 1SHH/B (1.55 Å, PPACK) over 1PPB/H (1.92 Å, PPACK). The peptide
protocol is not read off the symbol — it follows from `n_approved_peptide_binders
= 1` (bivalirudin, 20 aa, in `approved_biologics.csv`), which is the same fact
that made thrombin the target where the peptide arm had a control.

**CDK2: no, twice over.**

1. It never surfaces. CDK2 is **not in the top 500** of colorectal cancer's
   16,299 associated targets (verified by fetching 500 and searching; CDK4 is at
   361 and CDK6 at 363, CDK2 nowhere).
2. Even if it did, it fails the positive-control gate — `n_approved_binders_moa
   = 0`. Its seven annotations are off-target hits from CDK4/6 and multi-kinase
   drugs.

**On colorectal cancer the pipeline picks BRAF** (association rank 6, pipeline
rank 1): 5 mechanism-of-action approved drugs (dabrafenib, encorafenib,
regorafenib, sorafenib, vemurafenib), 8JNB chain B at 1.62 Å, not known-drug
dominated (`literature` 0.996 beats `clinical` 0.952, share 0.189).

```python
"braf": {
    "symbol": "BRAF", "uniprot": "P15056",
    "pdb_id": "8JNB", "chain": "B", "out": "braf",
    "protocol": "protein-anything", "binder_length": "60..90",
}
```

So: neither of the two hand-picked targets is what this stage would choose for
the demo's own disease, and the one it does re-derive (thrombin) it re-derives
for the right reason and on a better structure.

---

## 8. Ranking

Stated in full, because a composite nobody validated is worse than no composite:

1. **Gate** on `n_approved_binders_moa >= 1` **and** a usable structure.
2. **Order the survivors by Open Targets association score, descending.**

No weighting of tractability, safety or literature into a single number. There
is nothing here measured well enough to justify one, and inventing one would put
an unvalidated score in front of every downstream result. Tractability, safety
liabilities, target class and the per-datatype breakdown are all carried in the
output for a reader to weigh; they do not move the rank.

---

## 9. Not evaluated, and why

| Thing | Status | Why |
|---|---|---|
| Wild-type vs mutant construct | **not checked** | FGFR2's pick, 10OU, is "FGFR2 mutant D650V with compound 12". The rule has no mutation filter, so a mutant entry can win on resolution. Real gap; it needs the `rcsb_polymer_entity_align` mismatch list, which is one more field on the query |
| Construct coverage of the catalytic domain | **not checked** | `prep_target.py` enforces ≥95% coverage of the catalytic domain and relaxes in stated steps. Not ported — it needs the UniProt domain features and an aligned-region intersection per entry. The chosen chain's length is reported (`chain_entity_length_aa`) but nothing gates on it |
| Cryo-EM and NMR entries | **excluded by construction** | the search is `experimental_method == X-ray`. For targets with no crystal structure (ATM, POLD1, PROS1, C4BPB, SH2B3 here) this reports "no usable structure" when a cryo-EM entry may exist. Stated, not silently absorbed |
| AlphaFold fallback | **not implemented** | PROJECT_GOAL §4.2 allows `structure_source: "afdb"`. 5 of 50 targets across the two runs would need it |
| `mechanism_rationale`, `desired_effect` | **not implemented** | `TargetHypothesis` (§4.1) asks for a direction — agonise or antagonise. Every approved binder found here is recorded with its ChEMBL `target_action_type`, but nothing classifies the *desired* direction for the disease, which needs more than a database field |
| Literature → ranking | **fetched, verified, not used** | 78 verified PMIDs across the two runs, every one re-fetched and title-compared. They are attached to each target as evidence a reader can check; they do not move the rank. This is the same gap the old stage had, narrowed from "fetched and unused" to "fetched, verified, attached and unused" |
| `--exclude-known-drug-evidence` re-ranking | **not implemented** | WS-A A4 asks for a drug-naive ranking mode. §3's flag reports the overlap; it does not offer the alternative ranking. Under the current gate that mode would return the empty set, which is itself the finding |
| Entity cap | **capped at 60 + 60 per target** | 1,092 entries examined across 50 targets, out of tens of thousands available. `n_entities_total` vs `n_entities_examined` is written per target so the cap is visible. A target whose only clean entry is at rank 200 by resolution would be missed |

---

## 10. Running it

```bash
./env-kit/bin/python demo/autoresearch.py --disease "venous thromboembolism"
./env-kit/bin/python demo/autoresearch.py --disease "colorectal cancer" --top 25 --lit 2
./env-kit/bin/python demo/autoresearch.py --disease MONDO_0005575 --top 15 --lit 0
./env-kit/bin/python demo/autoresearch.py --disease "venous thromboembolism" --emit-config F2
./env-kit/bin/python -m demo.test_autoresearch
```

`--emit-config SYMBOL` prints one `demo/pipeline.py` `TARGETS` entry, keyed and
shaped for paste-in. It refuses, with the reason, for a target that fails either
gate — a config is never emitted for a target the pipeline cannot read.

The self-test (32 checks) covers: the `approved == 1` filter on its stated
counts; CDK2 failing the positive-control gate while passing a naive annotation
count; thrombin passing with its four MoA drugs and its one approved peptide;
4UD9 rejected for its hirudin chain while thrombin's own light chain does not
trigger the same rule; the chosen chain being the 259-aa heavy chain; CDK2's
selection reproducing 6Q4G/A; the anchors reproducing `[568]` and `[]`; the
circularity arithmetic on three cases; VTE resolving to MONDO_0005399 and
returning the coagulation cascade; every unusable target carrying a reason and
no config; and the emitted config using exactly the keys `demo/pipeline.py`
already knows, parsed out of that file with `ast` rather than assumed.
