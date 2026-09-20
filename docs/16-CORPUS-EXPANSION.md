# 16 — Corpus expansion: measuring and partly filling the approved-drug hole

**Builder:** `demo/expand_corpus.py` (re-runnable; `--offline` works from cache)
**Output:** `data/approved_drugs_extended.csv` — **240 rows**, a *new* file
**Untouched:** `data/approved_drugs.csv`, `data/drugs.csv`

Every number in this note was computed by `demo/expand_corpus.py`. Re-run
`./env-kit/bin/python demo/expand_corpus.py --offline --selftest` to reproduce it.

---

## 1. Why

`data/approved_drugs.csv` is 4,099 DrugCentral structures, 2,153 with
`approved == 1`. Every enrichment factor, base rate and positive-control list in
this project is computed against it. A drug the corpus is missing is therefore a
ceiling on what any method can be *shown* to do — and it fails silently, because
nothing in the pipeline can tell "this method can't find it" from "it isn't in
the corpus".

The motivating case is thrombin (`F2`, UniProt P00734). Measured on the corpus:

| | n |
|---|---|
| approved rows with `F2` in `targets` | 12 |
| approved rows with `F2` in `moa_targets` | 4 (argatroban, bivalirudin, dabigatran etexilate, ximelagatran) |

and six thrombin-active approved or formerly-approved drugs are absent — verified
by substring search over `name` in both `approved_drugs.csv` and
`approved_biologics.csv`, not assumed.

A separate investigation established that annotation enrichment cannot fix this:
ChEMBL holds 3,996 F2 activities at pChEMBL ≥ 6, but they belong to
investigational compounds. Only adding drugs fixes it.

## 2. Method

1. **Fetch** every ChEMBL molecule at `max_phase = 4` — approved by some agency,
   currently or historically. n = **4,225** (API `total_count`; the builder
   refuses to continue if the paged count disagrees).
2. **Standardise** with `scripts/build_drug_corpus.py:181 standardise()`,
   *imported, not reimplemented* — `LargestFragmentChooser(preferOrganic=True)`,
   `Uncharger`, canonical SMILES, InChIKey, and its flag vocabulary
   (`multi_fragment` / `inorganic` / `element` / `charged_parent`). Consistency
   with the existing corpus matters more here than any preference of mine.
3. **Diff** on the 14-character InChIKey skeleton, so salt and stereo variants of
   a parent already in the corpus do not count as missing.
4. **Rescue** ChEMBL entries typed `Small molecule` that carry no ChEMBL
   structure, by name lookup against PubChem, with an independent check (§5).
5. **Annotate** the additions from ChEMBL `mechanism` (→ `moa_targets`) and
   ChEMBL bioactivity (→ `targets`).

Every API response is cached under `data/raw/corpus_expansion/` (227 files,
17 MB, gitignored). A re-run makes no network calls.

## 3. The hole, measured

Of the 4,225 ChEMBL phase-4 molecules:

| | n |
|---|---|
| no SMILES in ChEMBL at all | 808 |
| unparseable / no InChIKey | 0 |
| standardised OK | 3,417 |
| → skeleton present in `approved_drugs.csv` | 3,162 |
| → **absent from `approved_drugs.csv`** | **255** |

So DrugCentral covers **92.5 %** (3,162 / 3,417) of ChEMBL's structurally-defined
approved drugs. The hole is 255 structures, plus a second and larger hole of a
different kind: the 808 phase-4 entries ChEMBL has no structure for.

### 3a. What kind of drugs are the 255?

| ChEMBL `molecule_type` | n |
|---|---|
| Small molecule | 241 |
| Protein | 10 |
| Oligonucleotide | 2 |
| Unknown | 2 |

| `first_approval` decade | n |
|---|---|
| (not recorded) | 42 |
| 1940s–1960s | 4 |
| 1970s | 15 |
| 1980s | 29 |
| 1990s | 9 |
| 2000s | 13 |
| 2010s | 24 |
| 2020s | 119 |

**The dominant cause is recency.** The latest `first_approval` anywhere in
`approved_drugs.csv` is **2021** (17 rows), which matches its provenance:
DrugCentral dbversion 49, dumped 2021-09-24. Of the 240 rows finally added, **86
have `first_approval` ≥ 2021 and 70 are ≥ 2022** — approvals that post-date the
snapshot and could not possibly be in it. 16 were approved in 2021, 16 in 2022,
28 in 2023 and 26 in 2024.

**But recency is not the whole story.** 117 of the 240 additions have a
`first_approval` *before* 2021 and 37 carry no date. Those are genuine coverage
gaps, not staleness: older non-US approvals, withdrawn compounds
(3,3',4',5-tetrachlorosalicylanilide, fenclozic acid, exifone, urethane), and
metal salts and coordination complexes.

Also measured: 15 of the 255 have a `pref_name` that *is* in the corpus already —
the same drug name on a different parent structure (usually a stereochemical or
tautomeric difference that survives skeleton matching).

### 3b. The 808 structureless phase-4 entries — the biologic hole

| ChEMBL `molecule_type` | n |
|---|---|
| Protein | 178 |
| **Small molecule** | **164** |
| Antibody | 150 |
| Unknown | 133 |
| Vaccine component | 54 |
| Enzyme | 51 |
| Oligonucleotide | 25 |
| Gene | 22 |
| Antibody drug conjugate | 13 |
| Oligosaccharide | 9 |
| Cell | 9 |

The 477 proteins, antibodies, enzymes, vaccines, genes, ADCs and cells are
out of scope for a SMILES corpus; `data/approved_biologics.csv` (262 rows with
sequences) is the right home, and §6 records which of them are already there.

The 164 typed `Small molecule` **are** in scope: ChEMBL simply has no structure
for them. Those were the target of the PubChem rescue.

## 4. What was added

**n = 240 rows** in `data/approved_drugs_extended.csv`.

| | n |
|---|---|
| from ChEMBL (structure in ChEMBL) | 180 |
| from PubChem (ChEMBL had no structure) | 60 |
| withdrawn (`withdrawn = 1`) | 16 |
| with ≥ 1 target annotation | 84 |
| with ≥ 1 mechanism-of-action target | 79 |
| **unflagged, drug-like organic** | **101** |
| flagged `metal_stripped` | 63 |
| flagged `inorganic` | 42 |
| flagged `element` | 20 |
| flagged `multi_fragment` | 95 |
| flagged `dup_parent_inchikey` | 59 |

Concatenated with the existing corpus: 4,339 rows, 2,393 with `approved == 1`
(from 4,099 / 2,153). Median MW of the additions is 392 Da against 326 Da for the
corpus — the additions skew to modern, larger drugs.

The additions bring **190 distinct gene symbols**, of which **57 appear in no
`targets` cell of the existing corpus**.

### Schema

The first twelve columns are exactly `approved_drugs.csv`'s, in order, so the two
files concatenate. Then:

| column | meaning |
|---|---|
| `source` | `ChEMBL 36 REST API` or `PubChem PUG REST (name lookup; ChEMBL has no structure)` |
| `source_url` | the exact record the SMILES came from |
| `withdrawn` | 0/1, ChEMBL `withdrawn_flag` |
| `chembl_id` | ChEMBL molecule id |
| `molecule_type` | ChEMBL's typing |
| `chembl_max_phase` | ChEMBL's raw max_phase |
| `approval_basis` | `chembl_max_phase_4` (n=239) or `active_moiety_of_approved_prodrug` (n=1) |
| `target_evidence` | `chembl_mechanism` and/or `chembl_pchembl>=6,n>=2` |
| `notes` | PubChem CID and Title where used; names collapsed onto this parent |

### `struct_id`

Assigned as **`900000 + row index`**, so 900000–900239. The DrugCentral corpus
runs 4–5462, and DrugCentral's own `struct_id` sequence is nowhere near 900000,
so this cannot collide now or after a DrugCentral refresh. The self-test asserts
the two id sets are disjoint.

### `approved` and `approval_basis`

`approved = 1` on every row, on one of two bases, recorded per row:

- **`chembl_max_phase_4`** (n = 239) — ChEMBL records the molecule as approved by
  some agency, currently or historically. Withdrawn drugs keep `approved = 1` and
  carry `withdrawn = 1`; a withdrawn drug is still a legitimate repurposing
  candidate, and ximelagatran, already in the corpus, is one.
- **`active_moiety_of_approved_prodrug`** (n = 1, dabigatran) — see §6.

Anyone who wants the strict regulatory reading can filter on
`approval_basis == "chembl_max_phase_4"`.

### Target annotations

- `moa_targets` — ChEMBL `mechanism` rows, resolved to **human gene symbols** via
  `target_components` `GENE_SYMBOL` synonyms. Non-human targets are dropped,
  matching `build_drug_corpus.load_targets`, which keeps
  `ORGANISM == "Homo sapiens"` only. 139 / 328 candidate molecules had ≥ 1
  mechanism row.
- `targets` — mechanism targets ∪ bioactivity targets, where a bioactivity target
  needs **≥ 2 human measurements at pChEMBL ≥ 6**. Same evidence rule as
  `data/target_annotations.csv` in this repo (`evidence = chembl_pchembl`),
  tightened with the `n ≥ 2` requirement so one stray measurement is not a target
  claim. 93 / 328 had ≥ 1 such target.

The bioactivity arm is not decoration: **melagatran has no ChEMBL mechanism row
at all**, and would have arrived as dead weight without it. It has 13 F2
activities at pChEMBL ≥ 6.

## 5. Structure provenance — and one thing that nearly went wrong

**No SMILES in this file was typed from memory.** Every one was fetched from
ChEMBL or PubChem, the record URL is on the row, and RDKit canonicalised it
locally. The self-test recomputes the InChIKey from the stored SMILES for all 240
rows and requires an exact match.

### The circular check that had to be thrown away

For the 164 structureless `Small molecule` entries, the plan was: resolve by name
against PubChem, keep it if PubChem's `Title` agrees with the name. That rejected
9 of 111 resolutions — including cisplatin (PubChem titles it
`azane;dichloroplatinum`) and oxaliplatin. Real drugs, thrown away by a string
test.

The obvious repair — "accept it if the name we asked for is in the CID's synonym
list" — is **circular and worthless**: PubChem's `/compound/name/` endpoint
resolves *through* that synonym index, so the test passes by construction. It
duly "rescued" `TECHNETIUM MEDRONATE` → CID 8094, **heptanoic acid**, which is
not technetium medronate by any reading.

What replaced it is independent of the name index: **if the name promises an
element, the structure must contain it**, checked on the raw pre-salt-strip
SMILES. Result:

| | n |
|---|---|
| structureless `Small molecule` entries | 164 |
| PubChem has no record under that name | 51 |
| Title agreed | 104 |
| Title disagreed, element check passed → kept, `notes` records `pubchem_name_index_only` | 8 |
| **rejected by the element check** | **1** (`TECHNETIUM MEDRONATE` → heptanoic acid) |
| resolved but already in the corpus | 40 |
| **new structures rescued** | **72** |

The 8 title-disagreements kept are cisplatin, cupric sulfate, oxaliplatin,
Prussian blue insoluble, mecobalamin, thallous chloride, ferric oxyhydroxide and
polidocanol — all correct, all recoverable only because the title test was not
the final word.

### Spot-check: InChIKey → PubChem → name

PubChem is a different database from ChEMBL, so a round-trip through it is
genuinely independent evidence: the InChIKey of the *stored* structure is handed
to PubChem, and PubChem is asked what compound that is.

Run A — `--spotcheck`'s own sample, n = 10 (the two named drugs plus a
deterministic spread across the file):

| row | InChIKey | PubChem CID | PubChem Title | |
|---|---|---|---|---|
| dabigatran | `YBSJFWOBGCMAKL-UHFFFAOYSA-N` | 216210 | Dabigatran | ✓ |
| melagatran | `DKWNMCUOEDMMIN-PKOBYXMFSA-N` | 183797 | Melagatran | ✓ |
| lutetium chloride | `AEDROEGYZIARPU-UHFFFAOYSA-K` | 24919 | Lutetium chloride | ✓ |
| gentamicin | `CEAZRRDELHUEMR-UHFFFAOYSA-N` | 3467 | Gentamicin | ✓ |
| zuranolone | `HARRKNSQXBRBGZ-GVKWWOCJSA-N` | 86294073 | Zuranolone | ✓ |
| maribavir | `KJFBVJALEQWJBS-XUXIUFHCSA-N` | 471161 | Maribavir | ✓ |
| telotristat | `NCLGDOBQAWBXRA-PGRDOPGGSA-N` | 25025298 | Telotristat | ✓ |
| bortezomib d-mannitol | `QDMRNLRJDHCHLB-DNNBANOASA-N` | 12990536 | Bortezomib D-mannitol | ✓ |
| sodium picosulfate | `UJIDKYTZIQTXPM-UHFFFAOYSA-N` | 5243 | Picosulfuric acid | ✓ (sodium stripped) |
| tedizolid | `XFALPSLJIHVRKE-GFCCVEGCSA-N` | 11234049 | Tedizolid | ✓ |

**10 agreed, 0 disagreed, 0 lookup failures.**

Run B — four rows chosen deliberately because an earlier build's random sample
had flagged them, i.e. the adversarial cases rather than the lucky ones:

| row | InChIKey | PubChem CID | PubChem Title | verdict |
|---|---|---|---|---|
| strontium chloride | `PWYYWQHXAPXYMF-UHFFFAOYSA-N` | 104798 | Strontium cation | counter-ion stripped, expected |
| calcium glubionate | `JYTUSYBCFIZPBE-AMTLMPIISA-N` | 7314 | Lactobionic Acid | counter-ion stripped, expected |
| **carboplatin** | `CCQPAEQGAVNNIA-UHFFFAOYSA-N` | 2568 | **1,1-Cyclobutanedicarboxylic acid** | **metal stripped** |
| **auranofin** | `SFOZKJGZNOBSHF-RGDJUOJXSA-N` | 88293 | **1-thio-β-D-glucopyranose tetraacetate** | **metal stripped** |

None of these four is a fetch error. They are `standardise()` doing exactly what
it is designed to do. The last two matter enough to have their own section.

### `metal_stripped` — 63 rows are a ligand, not a drug

`LargestFragmentChooser(preferOrganic=True)` keeps the largest *organic*
fragment. On a coordination complex that discards the metal centre: the stored
carboplatin structure is cyclobutane-1,1-dicarboxylic acid, with no platinum.
Auranofin loses its gold, cyanocobalamin its cobalt.

The fetch was not at fault: carboplatin's `source_url` is PubChem CID 426756,
`Platinum(II),1-cyclobutanedicarboxylato)diammine-, cis-`, which *is* carboplatin
with its platinum. The platinum is lost afterwards, in standardisation.

This is the **existing corpus's own behaviour** — the same imported function — so
the additions are consistent with it rather than better than it. But it is
silent there, and here it is not: the builder compares the metal content of the
raw and standardised structures and sets a **`metal_stripped`** flag. **63 of 240
rows** carry it. Any analysis that would be misled by a platinum drug with no
platinum should filter on it. The 101 rows with no flags at all are the clean,
drug-like subset.

### Deduplication

Deduped against the existing corpus by 14-char skeleton (the strict test; the
full-InChIKey test is weaker and also passes). Within the additions, 88 collisions
were found on the full InChIKey — salt-form siblings collapsed onto one parent by
standardisation (sodium / calcium / magnesium / lithium carbonate all become
carbonic acid). One row per parent is kept, and **the collapsed names are written
into the survivor's `notes`** with `dup_parent_inchikey` on its `flags`, so
nothing disappears without trace.

`magnesium hydroxide` collapses onto `water`. That is a standardisation artifact
and it is flagged `inorganic;element;dup_parent_inchikey`, not hidden.

## 6. The six named drugs

| drug | status | InChIKey | provenance |
|---|---|---|---|
| **dabigatran** | **ADDED** `struct_id=900046` | `YBSJFWOBGCMAKL-UHFFFAOYSA-N` | ChEMBL `CHEMBL48361`; `approval_basis = active_moiety_of_approved_prodrug` |
| **melagatran** | **ADDED** `struct_id=900124` | `DKWNMCUOEDMMIN-PKOBYXMFSA-N` | ChEMBL `CHEMBL266349`, `max_phase = 4`, `first_approval = 2004` |
| lepirudin | not added | — | 65-mer Protein, `CHEMBL1201666`, no SMILES. **Already in `data/approved_biologics.csv`.** |
| desirudin | not added | — | 65-mer Protein, `CHEMBL1201662`, `max_phase = 4`, 1997. **Not in the biologic corpus either — a real gap in `approved_biologics.csv`.** |
| hirudin | not added | — | No ChEMBL molecule with this `pref_name`. The natural leech protein is not a distinct approved drug entity; its recombinant analogues are (lepirudin, desirudin, bivalirudin), and all three are accounted for. |
| heparin | not added | — | Oligosaccharide, `CHEMBL4650299`, no SMILES. A polydisperse polysaccharide with no single structure — it cannot be represented in a SMILES corpus at all, and belongs in the biologic corpus. Not there either. |

### dabigatran: why it needed a special case

ChEMBL files dabigatran at **`max_phase = 3.0`**, so a plain `max_phase = 4` diff
would have missed it. The reason is regulatory, not chemical: the marketed entity
is the prodrug dabigatran etexilate (Pradaxa, FDA 2010), so the free acid never
had an approval of its own. It is nonetheless the species that actually binds F2,
and the only one worth having as a thrombin positive control.

`ACTIVE_MOIETIES` in the builder is the **only** place a human judgement enters
the file, it has one entry, and `approval_basis` carries that judgement onto the
row.

### Verifying dabigatran and melagatran are what they claim to be

Beyond the PubChem round-trip above, the internal consistency checks:

| | formula | MW |
|---|---|---|
| dabigatran (added) | C₂₅H₂₅N₇O₃ | 471.52 |
| dabigatran etexilate (already in corpus) | C₃₄H₄₁N₇O₅ | 627.75 |
| difference | C₉H₁₆O₂ | — |

The difference is exactly the hexyloxycarbonyl group on the amidine plus the
ethyl on the carboxylate — the two prodrug groups, and nothing else. The added
molecule is the etexilate with its prodrug groups removed, which is what
dabigatran is.

| | formula | MW |
|---|---|---|
| melagatran (added) | C₂₂H₃₁N₅O₄ | 429.52 |
| ximelagatran (already in corpus) | C₂₄H₃₅N₅O₅ | 473.57 |
| difference | C₂H₄O | — |

Again exactly the double-prodrug difference: the ethyl ester plus the
hydroxyamidine oxygen. Melagatran is ximelagatran's active metabolite, and the
structures are consistent with that and with each other.

## 7. F2 coverage, before and after

| | before | after |
|---|---|---|
| approved rows with `F2` in `targets` | 12 | **14** |
| approved rows with `F2` in `moa_targets` | 4 | **5** |

Added: **dabigatran** (F2, mechanism + bioactivity) and **melagatran** (F2,
bioactivity only — ChEMBL has no mechanism row for it).

This is a +17 % / +25 % change on a base of 12 / 4. It is small in absolute terms
and it is the honest number: the F2 hole was mostly biologics (lepirudin,
desirudin, hirudin, heparin), and biologics cannot be fixed by a SMILES corpus.
Three of the six named drugs are structurally unrepresentable here. The
structural additions available were two, and both are now present.

## 8. A defect found in `data/approved_drugs.csv` — not fixed here

The builder checks the existing corpus for InChIKey duplicates, because
`build_drug_corpus.py` dedupes by `struct_id`, not by InChIKey. There is exactly
one full-InChIKey duplicate in 4,099 rows, and it is not a salt form:

```
BODYFEUFKHPRCK-ZCZMVWJSSA-N: struct_id 5461 'sotorasib', struct_id 5462 'ibrexafungerp'
```

Both rows carry the identical SMILES, a triterpenoid. That triterpenoid is
**ibrexafungerp**. Sotorasib is a pyridopyrimidinone (PubChem CID 137278711,
`NXQKSXLFSAEQCZ-SFHVURJKSA-N`) and shares no scaffold with it. **Row 5461 has
sotorasib's name on ibrexafungerp's structure.**

The consequence: the corpus contains no correct structure for sotorasib, an
FDA-2021 KRAS G12C inhibitor. The ChEMBL diff confirms it independently —
sotorasib came back as *absent* from the corpus by skeleton match and is now in
`approved_drugs_extended.csv` with the correct structure.

There are also 76 skeletons shared by more than one `struct_id`; those are
ordinary salt and stereoisomer variants and the existing builder already flags
them `dup_parent_inchikey`.

**`data/approved_drugs.csv` was not modified.** It is load-bearing for published
numbers, and editing it in place would silently invalidate them. The defect is
recorded here and printed by every run of the builder.

## 9. Not added / why

| what | n | why |
|---|---|---|
| Proteins, antibodies, enzymes, vaccine components, genes, ADCs, cells (phase-4, no SMILES) | 477 | Not representable as SMILES. `data/approved_biologics.csv` is the right corpus; desirudin and heparin are absent from it and should be added there. |
| Oligosaccharides (incl. heparin) | 9 | Polydisperse polymers with no single structure. |
| Oligonucleotides (phase-4, no SMILES) | 25 | Out of scope for a small-molecule corpus. |
| `Unknown` type, no SMILES | 133 | ChEMBL has neither a structure nor a type; nothing to fetch on. |
| `Small molecule`, no ChEMBL structure, PubChem has no record under the name | 51 | Mostly mixtures, imaging agents and trade-name entries. |
| `TECHNETIUM MEDRONATE` | 1 | PubChem's name index resolved it to heptanoic acid. Rejected by the element check rather than admitted with a wrong structure. |
| 40 PubChem-resolved structures | 40 | Resolved correctly, but the parent is already in the corpus. |
| Remaining `max_phase < 4` molecules | — | Not approved. Only dabigatran was admitted below phase 4, with a stated basis. |

## 10. Reproducing

```bash
./env-kit/bin/python demo/expand_corpus.py             # build (cached fetches)
./env-kit/bin/python demo/expand_corpus.py --offline   # cache only, no network
./env-kit/bin/python demo/expand_corpus.py --selftest  # build + 15 assertions
./env-kit/bin/python demo/expand_corpus.py --spotcheck # + PubChem round-trip
```

The 15 self-test assertions — all passing at n = 240 — cover: the file re-reads;
the first twelve columns are the corpus schema; `source` / `source_url` /
`withdrawn` are present; every SMILES parses; every InChIKey recomputes from its
stored SMILES; InChIKeys are unique within the file; none is empty; no full
InChIKey and no 14-char skeleton collides with the existing corpus; `struct_id`s
are unique, disjoint from the corpus, and above the offset; every row has a
`source_url`; `withdrawn` is 0/1; and `n_targets` matches `targets`.

## 11. Known limitations

- **`approval_agencies` is empty on every added row.** ChEMBL does not carry the
  approving agency, and `first_approval` is a year, not a date. The
  DrugCentral-sourced rows have both; the additions have only the year, and 37
  have not even that. Any comparison that conditions on agency must exclude the
  additions.
- **The hole was measured against ChEMBL only.** The FDA Orange Book and the WHO
  INN recommended lists would each give a different denominator and were not
  used. "255 absent" is *absent relative to ChEMBL max_phase 4*, not an absolute
  count of approved drugs the world knows about.
- **63 of 240 rows have had a metal stripped**, and are the ligand rather than the
  drug. Flagged, not fixed — fixing it would mean diverging from the
  standardisation the existing corpus uses.
- **`targets` coverage is thin**: 84 of 240 rows have any target annotation. The
  156 without are not target-free, they are unannotated in ChEMBL's mechanism
  table and below the bioactivity threshold.
- **The 51 PubChem name-lookup misses were not pursued further** (DrugBank,
  UniChem and the FDA UNII registry would each recover some).
