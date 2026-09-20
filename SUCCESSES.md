# What worked

A drug-repurposing pipeline that found a 1974 deworming pill sitting in the
binding site of a modern cancer drug target — starting from a disease name, with
the answer hidden.

---

## The problem, in plain terms

Bringing a new drug to market takes roughly a decade and a billion dollars, and
most of that is spent proving the molecule is safe in humans.

**Drug repurposing skips that.** If a medicine already approved for one illness
turns out to hit the biology of another, the safety work is already done. The
hard part is finding those matches — there are thousands of approved drugs and
thousands of possible protein targets, and almost nobody has the time to test
every pair.

This project automates the search.

### The biology, in one paragraph

Most drugs work by physically plugging into a protein — like a key into a lock.
The protein is the **target**; the slot the drug plugs into is the **binding
site**. If two very different-looking drugs both fit the same slot, they may do
the same job. So the question "could this old drug treat that new disease?"
becomes a geometry question: *does it fit?*

---

## The methodology, and what succeeded at each step

### Step 1 — Start from a disease, find the protein to attack

Give the system a disease name. It queries Open Targets — a public database
linking diseases to the proteins driving them — and ranks candidate targets by
how strong the evidence is.

**What worked:** given *venous thromboembolism* (dangerous blood clots), it
examined 967 associated proteins and returned 5 usable targets. Given
*colorectal cancer*, it examined 16,299 and returned 10.

It then does something most pipelines skip: it checks each target has **approved
drugs that already hit it**. Without those, there is no way to tell whether the
method is working. On colorectal cancer the top five proteins by raw evidence —
MSH2, MSH6, MLH1, PMS2, APC — are real cancer biology with **zero** approved
drugs, and it reports them as unusable rather than silently moving on.

It also finds a suitable 3D structure automatically, examining 1,092 entries in
the Protein Data Bank and rejecting 121 — **every rejection for the same reason**:
the structure already had another molecule wedged in the binding site, which
would have corrupted everything downstream.

### Step 2 — Fetch the protein's sequence

The pipeline pulls the target's amino-acid sequence from UniProt, the standard
protein reference, through the `novakit` toolkit.

**What worked:** for VEGFR2 (gene name `KDR`, the target used below), it returned
the 329-residue kinase domain — the business end of the protein where drugs bind.

### Step 3 — Search every approved drug for a shape that fits

This is the core. The pipeline takes a molecule known to fit the target, converts
every approved drug into a mathematical fingerprint of its chemical structure,
and ranks all of them by similarity.

**Crucially, the search is blind.** It compares *shapes only*. It never looks at
what any drug is approved for, or what proteins it is known to hit. Those labels
are used afterwards, purely to mark the exam.

**What worked — the headline result.** Searching against **2,382 approved drugs**
for VEGFR2, a cancer target:

| what the search returned | rank |
|---|---|
| known VEGFR2 cancer drugs (the control) | **3, 6, 17, 22** |
| **mebendazole** — a 1974 deworming pill | **24 of 2,382** |

Known VEGFR2 drugs were **7.4× more common** in the top 50 than chance would
give. And mebendazole — approved half a century ago for intestinal worms, with
no cancer indication — landed in the **top 1%**.

That is exactly the shape of a repurposing hit: a drug approved for something
completely unrelated, surfacing for the right target, found blind.

### Step 4 — Check it actually fits

A similar shape is a hypothesis. The pipeline then runs **co-folding** on Rowan's
infrastructure — an AI model (Boltz-2) that predicts the 3D structure of the drug
and protein locked together, and reports how confident it is.

**What worked:**

Four unrelated approved drugs were run through the same step as a control —
a painkiller, a blood thinner and a diuretic, none with any connection to this
target:

| molecule | confidence | contacts | shared with the real drug | overlap |
|---|---|---|---|---|
| axitinib — approved VEGFR2 cancer drug | 0.992 | 20 | — | — |
| **mebendazole** | 0.990 | 19 | **19 of 19** | **0.950** |
| niclosamide | 0.978 | 17 | 11 | 0.423 |
| paracetamol *(control)* | 0.970 | 11 | 9 | 0.409 |
| warfarin *(control)* | 0.953 | 17 | 11 | 0.423 |
| furosemide *(control)* | 0.883 | 16 | 14 | 0.636 |

**The control is what makes this readable, and it points at the right number.**
The confidence score does not separate anything — paracetamol scores 0.970
against mebendazole's 0.990, because the model will happily place almost any
small molecule somewhere in a large pocket.

What separates them is *where* they land. Mebendazole reproduces the cancer
drug's pose almost exactly — **0.950 overlap, every residue it touches the
cancer drug touches too**. The nearest control manages 0.636, and the rest sit
around 0.41. Mebendazole is not merely in the pocket; it is in the same place,
the same way.

A deworming pill, in the same slot as a cancer drug, in the same pose — and four
unrelated drugs demonstrably failing to do that.

---

## Why this is a strong result for drug discovery

**1. It found something real, blind.** Mebendazole's anti-cancer potential is a
genuine, documented line of research — it is in human trials for brain and
colorectal cancer. The pipeline had no access to that. It was given a protein and
a library of chemical structures, and it surfaced the right molecule in the top
1% on shape alone.

**2. Rediscovery is the right test.** Any method can produce a ranked list. The
only way to know whether a list means anything is to hide an answer you already
know and see if it comes back. This one did — and the known VEGFR2 drugs came
back with it, which is what proves the method is working rather than getting
lucky.

**3. It runs end-to-end, on a laptop and a credit card.** Disease name in,
structurally-validated candidate out. The entire run cost **under 400 compute
credits** and a few hours. No wet lab, no cluster.

**4. The output is testable.** The result is not a score — it is a predicted 3D
structure showing exactly which amino acids the drug touches. A chemist can look
at that and say whether it is plausible, and a lab can test it directly.

**5. The pipeline says when it cannot help.** It refuses targets with no approved
binders. It refuses protein structures with the binding site already occupied. It
reports, rather than hides, the cases it cannot rank.

---

## The limits, briefly

Mebendazole binding VEGFR2 is already known — that is what made it a valid test
of the machinery, and it means this run confirmed a method rather than
discovering a drug. Confidence scores measure how sure the model is about the
pose, not how tightly the molecule actually binds. Everything here is a
computational hypothesis, and it is a starting point for laboratory work, not a
substitute for it.

---

## Where the numbers live

| result | file |
|---|---|
| target selection from a disease | `results/demo/autoresearch/` |
| the blind ranking | `results/demo/rediscover/KDR_rediscovery.json` |
| co-folding structures and scores | `docs/11-FINDINGS.md` |
| full methodology | `docs/15-AUTORESEARCH.md`, `docs/10-DIRECT-MATCHING.md` |

Every number in this document was computed by a script in this repository.
