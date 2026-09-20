# 12 — Target annotations: what measured bioactivity does and does not fix

`data/approved_drugs.csv` carries DrugCentral's curated `targets` / `moa_targets`
columns. A drug with no target can be neither a true positive nor a scored
negative in a retrieval benchmark, so it is dead weight in the denominator. This
note records an attempt to fill that hole from **measured** bioactivity (ChEMBL
activities carrying a pChEMBL value), what it recovered, and — more usefully —
what it proved is not recoverable that way.

Builder: `demo/enrich_targets.py`. Output: `data/target_annotations.csv`.
Every number below is printed by the script on each run; none is estimated.

---

## 0. First, the denominator is not what it looks like

The headline "48% of the library has no target annotation" is true of all 4,099
rows in `approved_drugs.csv`, but 4,099 is not the library. Only **2,153** rows
are `approved == 1` (see commit `145caec`). Splitting the two:

| Population | n | ≥1 curated target | no target |
|---|---:|---:|---:|
| all rows in the file | 4,099 | 2,114 (51.6%) | **1,985 (48.4%)** |
| `approved == 1` only | 2,153 | 1,669 (77.5%) | **484 (22.5%)** |

So 1,501 of the 1,985 unannotated rows are *not* approved drugs. The real
approved library is 77.5% annotated, not 52%. The gap is real but roughly half
the size the aggregate number suggests. Both denominators are reported below
because both are defensible; `approved == 1` is the one the benchmark runs on.

---

## 1. Matching drugs to ChEMBL

Every one of the 4,099 rows has an InChIKey, so InChIKey is the primary key and
name matching is a labelled fallback (`match_method` column), never silently
mixed in.

| Method | n (of 4,099) | share |
|---|---:|---:|
| `inchikey_exact` | 3,857 | 94.1% |
| `name_exact` (fallback) | 139 | 3.4% |
| unmatched | 103 | 2.5% |

Restricted to `approved == 1`: **2,139 / 2,153 matched (99.3%)** — 2,077 by
InChIKey, 62 by name. Only **10** approved drugs fail both. Matching is not the
bottleneck.

---

## 2. The assay-quality filter, and why it is not optional

The first run applied no assay filter and produced visibly wrong annotations:
aspirin → TSHR, fluorouracil → THPO/SMN1/TSHR, nicardipine → CYP2C9/CYP3A4. Two
distinct contaminants, both measured rather than assumed:

* **PubChem/Tox21 qHTS panels.** Sampling ChEMBL target `CHEMBL1963` (TSHR) at
  pChEMBL ≥ 6 returns **1000/1000** rows of `standard_type=Potency`,
  `assay_type=F`, `src_id=7`. `AC50` is the same story: **832/1000** sampled rows
  are `src_id=7`. These are screening readouts on a promiscuity-prone panel, not
  affinity constants.
* **ADMET / CYP inhibition panels** (`assay_type=A`). Real measurements, but of
  metabolism, not target engagement — and near-universal among lipophilic drugs.

For contrast, F2/thrombin (`CHEMBL204`) at pChEMBL ≥ 6 is **779 Ki + 221 IC50**,
all `src_id=1` (literature). The filter below leaves genuine pharmacology alone.

Rule: keep `standard_type ∈ {Ki, Kd, IC50, EC50}`, `assay_type ∈ {B, F}`,
relation `=`. Effect on the 137,958 fetched activity rows:

| dropped | reason |
|---:|---|
| 67,511 | target not a single protein (cell line / complex / no UniProt) |
| 8,342 | `standard_type=AC50` |
| 4,605 | `standard_type=Potency` |
| 1,744 | `assay_type=A` (ADMET) |
| 105 | `assay_type=T` (toxicity) |
| 26 | `standard_type=ED50` |
| 2 | `assay_type=U` |
| **82,335** | **total (59.7%)** |

55,623 rows survive → **13,498 distinct (drug, target) pairs** (10,170 human,
3,328 non-human).

The filter's effect on agreement with curated data, which is the real test:

| | before filter | after filter |
|---|---:|---:|
| drugs sharing ≥1 target with DrugCentral | 83.1% | **96.1%** |
| mean fraction of our targets that are curated | 65.0% | **89.6%** |
| mean MOA-target recall | 73.7% | **80.1%** |

Precision rose sharply and recall rose too, so nothing of value was cut.

---

## 3. Threshold

pChEMBL ≥ **6.0** (1 µM) is the default, the conventional cut and the one
`passes_threshold` flags. Rows are *fetched and written* down to **5.0**, so a
consumer can re-threshold without refetching. Measured sensitivity:

| pChEMBL ≥ | drugs w/ target (human) | rescued (all 4,099) | rescued (approved) | pairs | distinct targets |
|---:|---:|---:|---:|---:|---:|
| 5.0 | 1,437 | 70 | 38 | 10,170 | 1,314 |
| 5.5 | 1,313 | 53 | 24 | 8,454 | 1,174 |
| **6.0** | **1,200** | **33** | **14** | **6,440** | **1,043** |
| 6.5 | 1,071 | 23 | 9 | 4,662 | 878 |
| 7.0 | 968 | 16 | 4 | 3,480 | 740 |
| 8.0 | 739 | 6 | 2 | 1,840 | 464 |

Including non-human targets (§4):

| pChEMBL ≥ | drugs w/ target | rescued (all) | rescued (approved) | pairs | distinct targets |
|---:|---:|---:|---:|---:|---:|
| 5.0 | 1,642 | 177 | 101 | 13,498 | 2,230 |
| **6.0** | **1,373** | **102** | **63** | **8,831** | **1,756** |
| 7.0 | 1,119 | 59 | 39 | 5,052 | 1,274 |
| 8.0 | 861 | 35 | 27 | 2,738 | 805 |

Moving the cut from 6.0 to 5.0 buys 39 more rescued approved drugs at the cost
of admitting sub-micromolar-to-micromolar binding. Caffeine (§6) is the clearest
case of what sits in that band.

---

## 4. The unannotated tail is anti-infectives, so human-only is the wrong filter

After the human-only pass, 460 approved drugs matched ChEMBL but still had no
target. Their names are the explanation: abacavir, alatrofloxacin, albendazole,
amikacin, amphotericin B, artemether, artesunate, azlocillin, aztreonam,
bacampicillin, balofloxacin, benznidazole, biapenem. **The unannotated tail is
anti-infectives**, and their targets are pathogen proteins — excluded by
construction by `target_organism = Homo sapiens`.

Measured, not assumed: of a random sample of **120** of those drugs, **32 (27%)**
have a pChEMBL ≥ 6 affinity hit once any organism is allowed. Organisms seen:
HIV-1, *E. coli*, hepatitis C virus, *P. falciparum*, *L. donovani*, SARS-CoV-2,
*K. pneumoniae*, *T. cruzi*.

For a *drug-similarity* benchmark this matters: two fluoroquinolones sharing
bacterial DNA gyrase is exactly the signal the benchmark is meant to detect.
Non-human targets are therefore kept under a separate `evidence` value,
`chembl_pchembl_nonhuman`, with an `organism` column — opt-in, never silently
mixed into the human tier. This roughly triples the yield: rescued approved drugs
go from 14 to **63**; rescued rows overall from 33 to **102**.

> **Bug found and fixed while doing this.** Keying aggregation on gene symbol
> alone let a rodent ortholog with the same symbol occupy the key and flip a pair
> into the non-human tier, silently costing 157 human pairs (10,170 → 10,013).
> The key is now `(struct_id, organism, gene_symbol)`. Human pair counts are
> byte-identical to the human-only run, which is how the fix was verified.

---

## 5. Coverage, before and after

At pChEMBL ≥ 6.0, counting a drug as annotated if it has *either* a curated
target or a ChEMBL one:

| Population | annotated before | annotated after | rescued | distinct targets | median targets/drug | median drugs/target |
|---|---:|---:|---:|---:|---:|---:|
| all 4,099 — human only | 2,114 (51.6%) | 2,147 (52.4%) | 33 | 1,647 → 1,737 | 3 | 3 |
| all 4,099 — + non-human | 2,114 (51.6%) | 2,216 (54.1%) | **102** | 1,647 → **3,403** | 4 | 3 |
| approved 2,153 — human only | 1,669 (77.5%) | 1,683 (78.2%) | 14 | 1,565 → 1,648 | 4 | 3 |
| approved 2,153 — + non-human | 1,669 (77.5%) | **1,732 (80.4%)** | **63** | 1,565 → **3,206** | 5 | 2 |

**Be clear about the size of this.** Of the 1,985 unannotated rows, 102 were
rescued — **5.1%**. Of the 484 unannotated approved drugs, 63 — **13.0%**.
Coverage of the real approved library moves 77.5% → 80.4%. This is a real but
modest gain. The large number in the brief (1,985) was never 1,985 rescuable
drugs; most of those rows are not approved drugs at all, and most of the rest
have no measured affinity anywhere in ChEMBL.

The one genuinely large change is **distinct targets: 1,565 → 3,206** for the
approved library, because the pathogen proteome is broad. That also drops median
drugs-per-target from 3 to 2 — a wider, sparser target graph. For retrieval
that is a mixed blessing and should be measured, not assumed to help.

---

## 6. Agreement with the curated column, and three disagreements run down

n = 1,167 drugs have both a curated target set and a ChEMBL set at pChEMBL ≥ 6.

| Metric | Value |
|---|---:|
| drugs sharing ≥1 target with DrugCentral | **1,121 / 1,167 (96.1%)** |
| mean recall of curated targets | 57.8% |
| mean fraction of our targets that are curated | 89.6% |
| drugs with ≥1 MOA target recovered (n=851) | 724 (85.1%) |
| mean MOA-target recall | 80.1% |

96.1% concordance says the pipeline is not contradicting curated data. The 57.8%
recall is the interesting half. Three of the worst cases, checked individually:

**fluorouracil → TYMS (curated MOA, we miss it).** ChEMBL holds exactly **4**
5-FU/TYMS activities and **none carries a pChEMBL value** — they are `Activity`
and `Inhibition` percentage readouts. Not a threshold artifact and not a filter
artifact: the concentration-response constant does not exist upstream.
Pharmacologically correct too, since 5-FU is a prodrug and FdUMP is what
inhibits TYMS. Curation captures mechanism; bioactivity captures the assay that
happened to be run. **These are different things and neither subsumes the other.**

**caffeine → ADORA1/ADORA2A (curated, we miss it).** ChEMBL has 22 Ki rows for
caffeine at ADORA2A, spanning 5.0–48 µM, i.e. pChEMBL 4.32–5.30 — every one
below the 6.0 cut. The annotation is right and our threshold excludes it.
Caffeine genuinely is a micromolar adenosine antagonist that matters only because
the dose is large. This is the single best argument for keeping the 5.0–6.0 band
in the file rather than discarding it.

**liothyronine → F2 (we assert it; it is wrong).** Our only net-new approved F2
drug. Its three supporting activities are assays `CHEMBL4629641`, `CHEMBL4629643`
and `CHEMBL4629645`, whose descriptions read *"Agonist activity at GST-labelled
THRbeta"* and *"Agonist activity at human THRbeta"* — thyroid hormone receptor
beta. All three are assigned in ChEMBL to target `CHEMBL204`, which is
**Prothrombin / P00734 / F2**. This is an upstream curation error on a `THR`
abbreviation collision (thrombin vs thyroid hormone receptor), and our pipeline
inherited it faithfully. Flagging it: liothyronine is T3, an EC50 agonist result
against a serine protease is not chemically sensible, and `activity_type=EC50`
is unique among our F2 rows.

---

## 7. F2 (thrombin), specifically

This was the motivating example. It does not improve, and the reason is worth
more than the fix would have been.

| pChEMBL ≥ | curated | ChEMBL | union | net new |
|---:|---:|---:|---:|---:|
| 5.0 | 12 (4 MOA) | 10 | 14 | 2 |
| **6.0** | **12 (4 MOA)** | **8** | **13** | **1** |
| 7.0 | 12 (4 MOA) | 6 | 13 | 1 |

**F2 goes from 12 / 2,153 to 13 / 2,153, and the 13th is liothyronine — the
false positive from §6.** At the honest reading, F2 coverage does not improve at
all. At pChEMBL ≥ 5.0 it reaches 14, adding apixaban (5.51) and edoxaban (5.22),
both real but weak off-target thrombin binders.

The cause is the corpus, not the annotation. Checking which thrombin drugs exist
in `approved_drugs.csv`:

| in the corpus | absent from the corpus |
|---|---|
| dabigatran etexilate, argatroban, bivalirudin, ximelagatran, betrixaban, apixaban, edoxaban, rivaroxaban, warfarin | **dabigatran** (active form), melagatran, lepirudin, desirudin, hirudin, heparin |

ChEMBL holds 3,996 activities at pChEMBL ≥ 6 against F2, but they belong to
investigational compounds, not approved drugs, and the approved thrombin agents
that are missing are missing because they are peptides, biologics, or active
metabolites without a DrugCentral structure record. **No annotation source fixes
this. Only adding drugs to the corpus does.** If thrombin retrieval is a headline
result, the corpus is what needs work — see `data/approved_biologics.csv`.

---

## 8. Not annotated, and why

Verified programmatically, not assumed. Of `approved == 1`, n = 2,153, at
pChEMBL ≥ 6.0 with both tiers:

| n | reason |
|---:|---|
| 373 | matched ChEMBL, but zero single-protein pChEMBL activities (any organism) |
| 38 | ChEMBL activities exist, but all below pChEMBL 6.0 |
| 10 | no ChEMBL molecule matched (InChIKey *and* name both missed) |
| **421** | **total still unannotated (19.6% of approved)** |

Splitting the 411 matched-but-unannotated drugs a second way — by whether ChEMBL
holds *any* pChEMBL activity for them at all — gives the more useful picture:

| n | situation |
|---:|---|
| 212 | **no pChEMBL activity anywhere in ChEMBL**, at any threshold, any organism |
| 199 | pChEMBL activities exist, but every one is removed by the filters |

For those 199, the reason is overwhelmingly one thing (counted at the activity-row
level, n = 3,375 rows):

| rows | reason dropped |
|---:|---|
| 2,867 | **target is not a single protein** — ribosome, cell-wall/PBP complex, cell line |
| 225 | `standard_type=AC50` |
| 209 | `standard_type=Potency` |
| 56 | single protein, but pChEMBL < 6 |
| 18 | `assay_type` A or T |

This is the §4 problem again, one level deeper. Inspecting the 212 "nothing at
all" group by name shows two populations, and only one of them is hopeless:

* **Genuinely target-free**, and nothing will ever annotate them: imaging and
  contrast agents (`AMMONIA N-13`, `Choline C-11`, `Ioflupane I-123`, `UREA C-14`,
  acetrizoic acid, betiatide), sunscreens (amiloxate, bemotrizinol, bisoctrizole),
  nutrients and excipients (betaine, calcium pantothenate, 1-octacosanol),
  antiseptics (benzalkonium), gases (carbon dioxide, ammonia).
* **Anti-infectives whose target is real but not a single protein**: the group is
  thick with cephalosporins — cefamandole, cefathiamidine, cefetamet, ceforanide,
  cefpiramide, cefpirome, cefprozil, ceftaroline, ceftibuten, ceftizoxime,
  ceftolozane — plus azlocillin, carindacillin, carumonam, biapenem, besifloxacin,
  candicidin. β-lactams act on penicillin-binding proteins and the cell wall,
  aminoglycosides on the ribosome, amphotericin on ergosterol: complexes,
  assemblies and a lipid, none of which survive a `SINGLE PROTEIN` filter.

So the residual gap is **not** mostly "drugs with no target". It is roughly half
target-free products and half drugs whose target is not a single protein. The
first group should be **excluded from the retrieval denominator by a flag**
rather than counted as misses; the second needs a target vocabulary that admits
complexes, which is a larger change than this note covers.

Only the 10 unmatched rows are genuine pipeline misses.

---

## 9. What was not attempted, and why

| Source | Status |
|---|---|
| **ChEMBL** REST | Done. 137,958 activity rows fetched, 55,623 kept. |
| **BindingDB** | Not used, and the marginal yield is **not** measured. ChEMBL's source table shows `src_id=37` is *"BindingDB Patent Bioactivity Data"* only — the patent subset, not BindingDB's main literature set. So the common assumption that "BindingDB is already inside ChEMBL" is not supported by what we pulled, and querying it directly could add real rows. This is the largest unexplored source and is recorded as a gap, not as a judgement that it is worthless. |
| **PubChem BioAssay** | Deliberately *excluded*, not skipped. It reaches us as ChEMBL `src_id=7` and §2 shows it is the main source of false annotations. Adding it directly would make the file worse. |
| **DrugCentral bioactivity table** | Not reached. The curated `targets` / `moa_targets` columns already in `approved_drugs.csv` are DrugCentral's, and are carried through as the `drugcentral_moa` / `drugcentral_target` evidence tiers. The separate `drug.target.interaction.tsv` bulk file is gitignored and was not present in the working tree. |
| **UniProt/HGNC normalisation** | Not needed. ChEMBL's `GENE_SYMBOL` component synonyms are already HGNC and joined to the existing `targets` vocabulary without a single manual alias (verified: F2, KCNH2, CYP2D6 etc. match exactly). |

---

## 10. The file

`data/target_annotations.csv`, 28,337 rows, 2,291 distinct `struct_id`.

| `evidence` | rows | meaning |
|---|---:|---|
| `drugcentral_target` | 11,849 | curated, carried through from `approved_drugs.csv` |
| `chembl_pchembl` | 10,170 | measured, human single protein |
| `chembl_pchembl_nonhuman` | 3,328 | measured, pathogen/other organism |
| `drugcentral_moa` | 2,990 | curated mechanism-of-action target |

Columns: `struct_id`, `name`, `inchikey`, `gene_symbol`, `uniprot`, `organism`,
`evidence`, `activity_type`, `activity_value`, `activity_units`, `pchembl`,
`passes_threshold`, `n_activities`, `assay_types`, `molecule_chembl_id`,
`target_chembl_id`, `match_method`, `source_url`.

One row per (drug, organism, target), carrying the **strongest** measurement and
a count of how many activities supported it. Join on `struct_id`.

Integrity, checked: every human row has a gene symbol, a UniProt accession and
`organism == Homo sapiens`; every ChEMBL row has a `source_url`;
`passes_threshold` agrees with `pchembl >= 6.0` on every row. Four rows duplicate
on (`struct_id`, `gene_symbol`, `uniprot`) because ChEMBL carries both a generic
`Bacteria` target and a species-specific one for the same UniProt — dedupe on
(`struct_id`, `uniprot`) if that matters.

### Re-running

```bash
./env-kit/bin/python demo/enrich_targets.py            # rebuild (cache-first)
./env-kit/bin/python demo/enrich_targets.py --offline  # cache only, never network
```

No credentials. Raw API responses are cached gzipped under `data/raw/chembl/`
(already gitignored); a full cold run is ~675 requests and a warm re-run makes
**zero** network calls, which is how the numbers above were reproduced. If ChEMBL
is unreachable the script emits the DrugCentral tier alone and says so in the
provenance block rather than failing or writing a partial file.

---

## 11. Honest summary

* Coverage of the approved library: **77.5% → 80.4%**. 63 of 484 unannotated
  approved drugs rescued (13.0%); 102 of 1,985 rows overall (5.1%).
* Agreement with curated annotations: **96.1%** of drugs share ≥1 target;
  **85.1%** have ≥1 MOA target recovered.
* Distinct targets for the approved library: **1,565 → 3,206**, almost entirely
  from pathogen proteins.
* **F2/thrombin: 12 → 13 of 2,153, and the one addition is a false positive
  traced to a ChEMBL target misassignment. Effectively no change.**
* The biggest remaining limitation is **not** annotation coverage. It is that
  ~212 approved "drugs" (imaging agents, sunscreens, electrolytes, excipients)
  have no protein target at all and are silently counted as retrieval misses;
  that another ~199 are anti-infectives whose target is a ribosome or a cell-wall
  complex rather than a single protein; and that the corpus is missing the
  peptide and biologic drugs carrying the interesting targets.
* Largest unexplored source: **BindingDB's main literature set**. ChEMBL carries
  only its patent subset, so the overlap assumption usually made here is
  unverified.
