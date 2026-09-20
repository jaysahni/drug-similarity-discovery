# 09 — ESM-C protein embeddings: licence, size choice, and what the vectors separate

Covers `demo/esmc.py` and `demo/test_esmc.py`. Every claim is labelled
**MEASURED** (computed here, on this machine, today, with n given), **PUBLISHED**
(cited URL), or **INFERRED** (argued, not measured) — the convention from
`07-BOLTZGEN-BEHAVIOR.md`. No number below is illustrative.

Date of all measurements: **2026-09-19**. Machine: Apple M5 Pro, 24 GB unified
memory, macOS arm64 (`00-ENVIRONMENT.md`). Environment: `env-kit/`
(Python 3.12.13), with `torch 2.14.0` and `transformers 5.17.0` added by
`uv pip install --python env-kit/bin/python torch transformers`.

---

## 0. Verdict in one paragraph

**ESM-C 300M is MIT, it runs on this laptop, and its embeddings separate the
proteins they should.** The licence is clean: the HuggingFace card's second
`other` tag is a pointer to a third-party *dependency* table, not to a second
licence on the weights, and the linked repository's `LICENSE.md` is the
unmodified MIT text, Copyright 2026 Chan Zuckerberg Biohub — no field-of-use
clause, no non-commercial clause. **300M was chosen over 600M on measurement,
not size prejudice:** 600M is 1.6× slower per sequence and 1.7× larger on disk,
and on our 22-pair discrimination panel it *widened* neither the gap nor the
p-value (300M separates the families by 0.043 of cosine, 600M by 0.021; both
p = 1.34 × 10⁻⁵). The one unpleasant surprise is the context limit: **ESM-C does
not raise on a sequence longer than its 2048-token trained context** — it
extrapolates its rotary embeddings and returns finite, ordinary-looking numbers.
`demo/esmc.py` refuses those inputs itself.

---

## 1. Licence — the blocker check

Our repo is MIT and intended to be public (`CLAUDE.md`; Regeneron's challenge
prefers an MIT-licensed public repo). So this was checked before anything else.

**PUBLISHED.** The model card front-matter of
[`biohub/ESMC-300M`](https://huggingface.co/biohub/ESMC-300M/raw/main/README.md)
— identical on `ESMC-600M` and `ESMC-6B` — reads:

```yaml
license:
- mit
- other
license_link: https://github.com/Biohub/esm/blob/main/THIRD_PARTY_NOTICE.md
```

The dual `mit` + `other` tag is what made this worth checking: an `other` entry
can hide an acceptable-use rider. It does not here.

**PUBLISHED.** The linked file,
[`THIRD_PARTY_NOTICE.md`](https://raw.githubusercontent.com/Biohub/esm/main/THIRD_PARTY_NOTICE.md)
(HTTP 200, read in full on 2026-09-19), is *nine table rows listing the licences
of the `esm` package's own dependencies* and nothing else:

| Library | Licence stated |
|---|---|
| flash-attn, PyTorch, xformers, omegaconf, scipy | BSD / BSD-3-Clause |
| jaxtyping, einops, attrs | MIT |
| lightning / torchmetrics | Apache 2.0 |

It grants nothing, restricts nothing, and says nothing about the weights. It is
an attribution notice. **INFERRED, but on a complete reading:** the `other` tag
exists only because HuggingFace requires a tag whenever `license_link` points at
a non-standard file. There is no second licence.

**PUBLISHED.** The licence proper is
[`LICENSE.md`](https://raw.githubusercontent.com/Biohub/esm/main/LICENSE.md) at
the root of `github.com/Biohub/esm` (1095 bytes, read in full): the verbatim MIT
text, headed `**License (MIT)**`, `Copyright 2026 Chan Zuckerberg Biohub, Inc.`,
with the standard "without restriction … use, copy, modify, merge, publish,
distribute, sublicense, and/or sell" grant. No added clause of any kind.

**Finding: MIT confirmed, no blocker.** BSD, MIT and Apache-2.0 dependencies are
all compatible with shipping an MIT repo that calls this model.

Two caveats recorded so nobody re-derives them:

- GitHub's own repo metadata reports `spdx_id: NOASSERTION` for `Biohub/esm`.
  That is only because the file is named `LICENSE.md` with a bold Markdown
  heading, which GitHub's licence detector does not match. The text is MIT.
- The prior session's note that "earlier sources report a tiered/non-commercial
  licence — those are stale" is **corroborated**: nothing reachable from the
  current model cards mentions tiers or non-commercial use. We do not distribute
  the weights, so no notice obligation attaches to this repo beyond citation.

---

## 2. Corrections to the prior session's notes

The starting notes were close but not exact. **MEASURED** against the
HuggingFace API on 2026-09-19:

| Prior note | What is actually true |
|---|---|
| Weights at `biohub/esmc-300m-2024-12` | That repo exists but holds **only** `data/weights/esmc_300m_2024_12_v0.pth` — the legacy `esm`-package format, no tokenizer, no safetensors. Its own card says to use `biohub/ESMC-300M` instead. |
| `EvolutionaryScale/esmc-600m-2024-12` | HTTP **307**, redirecting to `biohub/esmc-600m-2024-12`. The EvolutionaryScale org name now forwards to `biohub`. |
| `pip install esm` needed | **Not needed.** `transformers` 5.17.0 ships native `EsmcModel` / `EsmcTokenizer` / `EsmcConfig`. `trust_remote_code` is not required either. `demo/esmc.py` depends only on `torch` + `transformers`. |
| 300M → 960-d, 600M → 1152-d, 6B → 2560-d | **Confirmed** from each repo's `config.json`: `hidden_size` 960 / 1152 / 2560, `num_hidden_layers` 30 / 36 / 80. |
| 1280-d is ESM-2 650M, a different model | Confirmed by absence: no ESM-C config has `hidden_size: 1280`. |

The canonical repos to use are **`biohub/ESMC-300M`**, `biohub/ESMC-600M`,
`biohub/ESMC-6B`, all `gated: false`.

---

## 3. Size choice — 300M, on measurement

`biohub/ESMC-300M`, 332.0 M parameters (**MEASURED**, counted from the loaded
module), 960-d.

### 3.1 Cost

**MEASURED**, all on this machine today. Download is a one-off over the
university network; load and embed are the numbers that matter per run. Embed
timings are cold-cache (cache redirected to a fresh temp directory to force real
inference), batch of n = 8 real human proteins totalling 2238 residues, one timed
run after a discarded warm-up pass.

| | ESMC-300M | ESMC-600M |
|---|---:|---:|
| parameters (counted) | 332.0 M | 574.0 M |
| `model.safetensors` on the hub | 1.33 GB | 2.30 GB |
| first download, cold | 58.9 s | 99.8 s |
| load from disk to MPS | 3.6 – 4.7 s | 3.6 s |
| embed, **CPU**, n = 8 | 4.32 s → **540 ms/seq** | 7.12 s → **890 ms/seq** |
| embed, **MPS**, n = 8 | 2.18 s → **272 ms/seq** | 3.49 s → **437 ms/seq** |

Load time varied 3.6–4.7 s across runs on both models; it is dominated by file
I/O, not parameter count, so it does not discriminate between the two.

### 3.2 Quality on the one panel we have

**MEASURED**, the discrimination panel of §5.2 (4 human CDKs against 4 unrelated
human proteins; 6 within-family pairs, 16 between-family pairs, n = 22):

| | ESMC-300M | ESMC-600M |
|---|---:|---:|
| mean within-family cosine | 0.9011 | 0.9236 |
| mean between-family cosine | 0.6211 | 0.6665 |
| min within − max between (the gap) | **+0.0426** | +0.0205 |
| every within > every between | yes | yes |
| Mann-Whitney U, one-sided | U = 96, p = 1.34 × 10⁻⁵ | U = 96, p = 1.34 × 10⁻⁵ |

The p-values are identical because the test is on ranks and both models produce
the same complete separation; with n = 22 that is the floor this test can reach.
**600M does not separate this panel better than 300M** — by the only continuous
statistic available here, the gap, it separates it *less* well, though on 22
pairs that difference is not itself significant and we do not claim 300M is the
better model in general. It is enough to say the extra 1.6× of inference time
buys nothing we can measure.

**Decision: 300M.** At the measured 272 ms/sequence a 20-entry biologic
library embeds in 5.4 s of inference, and in 0 s on every later run (§4.2). That is the right trade for a hackathon laptop.

### 3.3 MPS versus CPU — measured, not assumed

**MEASURED.** MPS is **2.0–2.1× faster** than CPU on both sizes (300M:
2.18 s vs 4.32 s; 600M: 3.49 s vs 7.12 s, n = 8 each). Agreement between the two
backends on the full hidden-state tensor is **max |Δ| = 1.32 × 10⁻⁵** in fp32 —
numerically irrelevant at four decimal places of cosine. `demo/esmc.py` therefore
prefers MPS, with `dev=` to override.

---

## 4. Implementation decisions

### 4.1 Pooling

Mean over **residue positions only**: `<cls>`, `<eos>` and padding are masked out
before averaging. **MEASURED**: a padded batch and one-at-a-time inference give
`max |Δ| = 0.0` *exactly* on the hidden states, so batching is free of fidelity
cost and the module groups sequences by length to keep padding small.

### 4.2 Cache

One `.npy` per sequence under `results/demo/_cache/esmc/`, named
`sha256(model_name + "\n" + sequence)[:24]`. The model name is in the key because
300M and 600M vectors are different widths and must never collide.

**MEASURED**: re-embedding the self-test's 11 sequences from a warm cache takes
**1 ms** total and loads no model at all, against 2.61 s (plus a ~4 s model load)
cold. 11 vectors occupy 42 KB.

### 4.3 Multi-chain entities — `combine=`, default `"mean"`

Antibodies are two chains; haemoglobin is two distinct chains. There is no
settled way to make one vector from several, so this is a parameter with four
modes, and the decision is recorded here rather than buried.

**Default is `"mean"`: embed each chain separately, average the vectors, each
chain weighted equally.** Two reasons:

1. **It stays in R⁹⁶⁰.** The point of this module is a library that mixes
   two-chain antibodies with single-chain enzymes and hormones. `"concat"` puts
   a two-chain entity in R¹⁹²⁰, where it is comparable only to other two-chain
   entities in the same order — `embed_entities` raises rather than silently
   stacking ragged widths, and the self-test checks that it does.
2. **Equal weighting, not length weighting.** A full IgG heavy chain is roughly
   twice the light chain's length, but the light chain carries three of the six
   CDRs. Length weighting lets the heavy chain dominate: **MEASURED** on IgG1
   constant regions (IGHG1 399 aa + IGKC 107 aa), the `"length"` combination sits
   at cosine 0.9951 to the heavy chain but only 0.9360 to the light, while
   `"mean"` sits at 0.9725 / 0.9750 — near-equidistant, which is the intent.

`"chainbreak"` joins the chains with `|`, **ESM-C's own chain-break token** —
verified present in the tokenizer vocabulary — and embeds them in a single pass.
It is the only mode where the model attends across both chains at once, so it is
the most principled, and it is offered for that reason; it is not the default
only because the chains must jointly fit the 2046-residue context. **MEASURED**,
it behaves much like `"length"` (IgG1: 0.9921 heavy / 0.9381 light), which is
expected — it *is* a single pooled pass over the concatenated residues.

**Not evaluated:** embedding only the variable region (Fv). It is a real option
and probably the right one for discriminating *between* antibodies, since the
constant regions are near-identical across an IgG1 library and would dominate the
cosine. It is not done here because we have no annotated CDR boundaries in
`data/`, and guessing them would put an unverified number in the pipeline.

### 4.4 Sequence-length limit — ESM-C does not enforce its own

**MEASURED, and the most important behavioural finding here.** ESM-C's
`config.json` gives `max_position_embeddings: 2048`, and the model card states
the context was raised from 512 to 2048 for the final 500k training steps. That
is a *training* context, not a hard cap. Feeding the model inputs of 2046, 2047,
2100 and 4000 residues:

| residues | tokens | result |
|---:|---:|---|
| 2046 | 2048 | OK, finite |
| 2047 | 2049 | **OK, finite** |
| 2100 | 2102 | **OK, finite** |
| 4000 | 4002 | **OK, finite** |

Nothing raises. Nothing warns. The rotary embeddings extrapolate and the model
returns ordinary-looking numbers that were produced outside the regime it was
trained in. This is the failure mode `CLAUDE.md` exists to prevent, so
`demo/esmc.py` refuses anything over **2046 residues** (2048 tokens less `<cls>`
and `<eos>`) with a message saying exactly this, rather than truncating.

**PUBLISHED**, and the prior session's note about ESM-2 is confirmed:
`facebook/esm2_t33_650M_UR50D`'s `config.json` gives `max_position_embeddings:
1026` — 1024 positions less `<cls>`/`<eos>`, so 1022 residues — with
`hidden_size: 1280`, which is the 1280-d figure that must not be confused with
any ESM-C. **ESM-C's usable context is twice ESM-2's, and the practical question
is settled: a full IgG heavy chain is
~450 aa and the longest sequence in the self-test panel is human albumin at
609 aa — both comfortably inside 2046.**

The module also refuses characters outside ESM-C's alphabet. **MEASURED**: the
tokenizer maps lowercase `a` to `<unk>` (id 3) without complaint, so a
lowercased FASTA would embed as a string of unknown residues and produce a
plausible, meaningless vector.

### 4.5 Cosine

`similarity()` is cosine, the conventional similarity for protein language model
embeddings. Mean-pooled PLM vectors vary in norm with length and composition, and
cosine discards exactly that. It raises on a dimension mismatch, because 960-d
and 1152-d vectors are not comparable and silently broadcasting them is the kind
of error that produces a number instead of a crash.

---

## 5. Self-test results

`./env-kit/bin/python -m demo.test_esmc`. Every sequence is a full UniProt
canonical entry with its accession inline in the test file. The orderings
asserted were written before the numbers were seen; none is a fitted threshold.

**MEASURED**, one run, 2026-09-19, ESMC-300M on MPS: loaded in 4.41 s; embedded
**n = 11** sequences (2886 residues) in 2.61 s = **237 ms/sequence**, dim 960.

### 5.1 Identity and cache

- Self-cosine of all 11 vectors with themselves: largest deviation from 1.0 is
  **2.22 × 10⁻¹⁶** (float64 rounding). Passes.
- Re-embedding an already-cached sequence returns a bit-identical vector
  (`max |Δ| = 0.0`). Passes.

### 5.2 CDK family versus unrelated proteins

Four human CDKs — CDK1 (P06493), CDK2 (P24941), CDK4 (P11802), CDK6 (Q00534) —
against four unrelated human proteins — insulin (P01308), lysozyme C (P61626),
albumin (P02768), haemoglobin β (P68871).

| within-family (n = 6) | cosine | | between-family (n = 16) | cosine |
|---|---:|---|---|---:|
| CDK1~CDK2 | 0.9707 | | CDK1~INS | 0.4811 |
| CDK1~CDK4 | **0.8467** ← min | | CDK1~LYZ | 0.7591 |
| CDK1~CDK6 | 0.8872 | | CDK1~ALB | 0.2990 |
| CDK2~CDK4 | 0.8610 | | CDK1~HBB | 0.7896 |
| CDK2~CDK6 | 0.8730 | | CDK2~INS | 0.5466 |
| CDK4~CDK6 | 0.9679 | | CDK2~LYZ | 0.7721 |
| | | | CDK2~ALB | 0.2354 |
| | | | CDK2~HBB | **0.8041** ← max |
| | | | CDK4~INS | 0.6383 |
| | | | CDK4~LYZ | 0.7497 |
| | | | CDK4~ALB | 0.4724 |
| | | | CDK4~HBB | 0.7490 |
| | | | CDK6~INS | 0.6140 |
| | | | CDK6~LYZ | 0.7545 |
| | | | CDK6~ALB | 0.5222 |
| | | | CDK6~HBB | 0.7500 |

- **Complete separation.** Min within-family (0.8467, CDK1~CDK4) exceeds max
  between-family (0.8041, CDK2~HBB). Passes.
- Within mean **0.9011** (n = 6) vs between mean **0.6211** (n = 16),
  **Mann-Whitney U = 96, one-sided p = 1.34 × 10⁻⁵**. Passes.

### 5.3 Shuffled control

CDK2 with its 298 residues shuffled (Python `random.Random(0)`, composition
verified unchanged by a sorted-multiset assertion in the test):

| pair | cosine |
|---|---:|
| CDK2 ~ CDK6 (homolog) | 0.8730 |
| CDK2 ~ insulin (unrelated but real) | 0.5466 |
| CDK2 ~ shuffled CDK2 (same composition, not a protein) | **−0.0115** |

The homolog beats the shuffle by **+0.8845** of cosine. Passes. The shuffle lands
essentially *orthogonal* to the original — strong evidence the embedding encodes
sequence order and motif structure, not amino-acid composition, which is exactly
the thing a bag-of-residues baseline would get wrong.

### 5.4 Multi-chain

Haemoglobin (HBA P69905 + HBB P68871) and the IgG1 constant regions
(IGHG1 P01857 + IGKC P01834), against albumin as the negative control:

| entity | mode | dim | ~chain 1 | ~chain 2 | ~albumin |
|---|---|---:|---:|---:|---:|
| haemoglobin | mean | 960 | 0.9926 | 0.9926 | 0.3137 |
| haemoglobin | length | 960 | 0.9923 | 0.9928 | 0.3138 |
| haemoglobin | chainbreak | 960 | 0.9650 | 0.9678 | 0.4531 |
| haemoglobin | concat | 1920 | — (not in the single-chain space) | | |
| IgG1 constant | mean | 960 | 0.9725 | 0.9750 | 0.5272 |
| IgG1 constant | length | 960 | 0.9951 | 0.9360 | 0.5301 |
| IgG1 constant | chainbreak | 960 | 0.9921 | 0.9381 | 0.5241 |
| IgG1 constant | concat | 1920 | — | | |

Every mode places the combined entity far closer to its own chains than to
albumin. All pass. `combine="concat"` correctly refuses to stack a two-chain
entity with a one-chain one.

### 5.5 Guards

- 2047-residue input refused with the over-context error. Passes.
- Longest real chain in the panel (albumin, 609 aa) is inside the limit. Passes.
- Lowercased insulin refused. Passes.

**All 15 checks pass.** No check was weakened to get there.

---

## 6. Caveats and honest limits

- **The cosines are compressed upward.** CDK2 and haemoglobin β share no
  homology, and they score 0.8041. Mean-pooled PLM embeddings are anisotropic —
  they occupy a narrow cone — so *absolute* cosine is not interpretable as
  "percent related". Only the ordering is. Any downstream threshold must be
  calibrated on this repo's own distribution, never imported from a paper. The
  0.8467 / 0.8041 separation in §5.2 is complete but narrow, and on a bigger,
  more adversarial panel it would not be safe to assume it stays so.
- **n = 22 pairs, one protein family, one species.** This is a sanity check that
  the vectors carry biology, not a benchmark of ESM-C. It licenses the statement
  "the embeddings separate CDKs from non-CDKs"; it does not license "ESM-C ranks
  drug targets well".
- **Albumin is an outlier in both directions** (0.2354–0.5222 to everything). At
  609 aa it is the longest sequence in the panel; mean pooling over more residues
  regresses toward the corpus mean. Worth knowing before a long biologic enters
  the library.

### Not evaluated, and why

| Thing | Why not |
|---|---|
| ESMC-6B (2560-d) | 6 B parameters in fp32 is ~24 GB of weights against 24 GB of unified memory shared with the OS. Not attempted; would need quantisation, and the 300M/600M comparison already shows scale is not the binding constraint here. |
| ESM-C versus the toolkit's ESM-2 | A head-to-head belongs in the benchmark proper, on a labelled task, not in a module self-test. `demo/embed.py` is another agent's file. |
| Fv-only antibody embedding | No annotated CDR/variable-region boundaries in `data/`. See §4.3. |
| `flash_attention_2`, bf16 | The model card recommends them for CUDA. No NVIDIA GPU here (`00-ENVIRONMENT.md`). |
| Real therapeutic antibody sequences (nivolumab, pembrolizumab) | Their variable regions are not in UniProt as such. The self-test uses IGHG1/IGKC constant regions instead, which are verifiable UniProt entries; the multi-chain machinery is what is being tested, and it does not care which antibody. |

---

## 7. Reproducing

```bash
uv pip install --python env-kit/bin/python torch transformers
./env-kit/bin/python -m demo.esmc        # loads, embeds CDK2, prints timings
./env-kit/bin/python -m demo.test_esmc   # the 15 checks above
```

First run downloads 1.33 GB from HuggingFace (58.9 s measured). Everything after
that is local, and any sequence embedded once is free thereafter.
