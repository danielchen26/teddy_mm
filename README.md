<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/hero-dark.png">
  <img alt="TEDDY × ANM — keep TEDDY frozen; make its answers editable, honest and auditable with ANM" src="docs/assets/readme/hero-light.png" width="100%">
</picture>

<p>
  <a href="https://danielchen26.github.io/teddy_mm/"><b>Live site</b></a> ·
  <a href="https://danielchen26.github.io/teddy_mm/anm-loop.html"><b>Decision-loop report</b></a> ·
  <a href="#how-anm-and-teddy-connect">How ANM connects to TEDDY</a> ·
  <a href="#quickstart">Quickstart</a> ·
  <a href="#key-terms">Key terms</a> ·
  <a href="docs/reports/">Proof reports</a> ·
  <a href="#claim-boundary">Claim boundary</a>
</p>

<img alt="TEDDY: frozen, never retrained" src="https://img.shields.io/badge/TEDDY-frozen%2C%20never%20retrained-c8641f?style=flat-square">
<img alt="scope: decision layer (Mode B)" src="https://img.shields.io/badge/scope-decision%20layer%20(Mode%20B)-008c7e?style=flat-square">
<img alt="data: CITE-seq, site4, 16,750 cells" src="https://img.shields.io/badge/data-CITE--seq%20%C2%B7%20site4%20%C2%B7%2016%2C750%20cells-6554c9?style=flat-square">
<img alt="numbers: from committed reports" src="https://img.shields.io/badge/numbers-from%20committed%20reports-4c566a?style=flat-square">

</div>

> Linked names (in blue) jump to [Key terms](#key-terms), where each one is explained in plain words. On the [live site](https://danielchen26.github.io/teddy_mm/) they are underlined and explain themselves on hover or tap, and a guide column shows the whole workflow.

## What this is

**[TEDDY](#t-teddy)** ([Merck TEDDY-G 70M](https://huggingface.co/Merck/TEDDY) · [paper](https://arxiv.org/abs/2503.03485)) is a single-cell foundation model. Kept frozen, it reads a cell's RNA and turns it into an [embedding](#t-z512); a small head trained in this repo on those embeddings predicts the cell's surface proteins (134; we use 9) on [CITE-seq](#t-citeseq) data (NeurIPS 2021 BMMC, [held-out site4](#t-holdout), phase-1 [Pearson ≈ 0.61](#t-pearson)).

**[ANM](#t-anm)** (Active Neural Matter) is a decision layer. It reads TEDDY's predictions as [typed evidence](#t-evidence) and answers a written-down question, the [observer](#t-observer): call a lineage, or say an honest [“no call”](#t-abstain). Each answer comes with a [workability](#t-pf) check (can it decide?), an [exactness](#t-qf) check (was it right?) and a per-marker reason.

**This repository** is the bridge between the two, plus six head-to-head tests of where ANM *helps* TEDDY. TEDDY's weights and predictions never change; only the decision layer does.

## How ANM and TEDDY connect

Three roles, one direction. Nothing flows back into TEDDY.

| Role | Reads | Produces | Trained? |
|---|---|---|---|
| **[TEDDY](#t-teddy)** · evidence source | the cell's RNA only | a 512-number [embedding](#t-z512); a small head trained in this repo predicts 134 surface proteins, 9 of which go to ANM | TEDDY-G 70M frozen, never retrained; the head was trained once in this repo |
| **[ANM](#t-anm)** · decision layer | normally TEDDY's outputs only (the 9 predictions as [typed evidence](#t-evidence)); labelled exceptions: the protein-only stand-in when RNA is missing (Block 2, reused in Block 4) and its average with TEDDY in “both”, measured protein in Block 6's complement test. Block 6's look-alike pairs are picked by a separate step outside ANM | a lineage call or [“no call”](#t-abstain), with checks | no learned weights: the [question](#t-observer) and the engine's settings are a few declared numbers |
| **Measured proteins** · answer key | (held out) | the grade for each call, keyed per question (soft/strict: plain mean of each lineage's 3 markers, no key-marker doubling, with that rule's margin; key-marker: weighted key marker), when the top lineage clearly leads | none. The key has 3 classes and no NK class, so most NK cells, which carry CD16, count as “myeloid” (1,342 of 1,690 under the soft rule; [script](scripts/answer_key_vs_annotation.py)) |

**Is this a multimodal model?** TEDDY is not: it reads only RNA, and nothing fuses protein into it. The *data* are multimodal (RNA and protein measured in the same cell), and normally the protein is only the answer key. It also reaches the evidence in labelled places: in Block 2, with RNA missing, this repo's phase-2 model (trained on TEDDY's embedding plus protein, run here with its RNA input off) [stands in](#t-adtonly), and the [both](#t-joint) condition averages TEDDY's predictions with it (Block 4's RNA-missing view reuses the stand-in); Block 6 uses measured protein to pick look-alike pairs and, in its [complement](#t-complement) test, as ANM's evidence.

**What ANM reads in each block**

| Block | ANM's evidence | Graded against |
|---|---|---|
| 1 · Edit the question | TEDDY's 9 predictions, asked three ways (soft · strict · key-marker) | measured proteins, keyed per question |
| 2 · Missing modality | RNA: TEDDY's predictions · RNA missing: the protein-only stand-in, no TEDDY · both: their average | measured proteins |
| 3 · Attribution | TEDDY's 9 predictions, then the same with one marker removed at a time | measured proteins |
| 4 · Scope gate | TEDDY's 9 predictions (RNA-missing view: the stand-in) | measured proteins |
| 5 · Zero-label transfer | TEDDY's 9 predictions, 0 labels. The trained heads it is compared with learn from those 9 plus 32 numbers of TEDDY's embedding (the threshold grid uses the 9 only) | measured proteins |
| 6 · Decision loop | Pairs are picked outside ANM: neighbours in TEDDY's embedding whose measured proteins disagree. ANM then calls each cell from TEDDY's predictions, and from measured protein in the complement test | measured proteins |

**Two separate axes: phases build the evidence, modes say how ANM uses TEDDY**

| | Axis | What it is | Status |
|---|---|---|---|
| [Phase 1](#t-phase1) | building the evidence | Frozen TEDDY embeds RNA once; a small head trained in this repo maps it to all 134 surface proteins (MLP head 0.610 vs latent flow matching 0.595, test Pearson) | done; every block uses its predictions |
| [Phase 2](#t-phase2) | building the evidence | Multimodal fusion model: TEDDY's RNA embedding plus a protein encoder, trained with one modality randomly dropped (RNA in about 15% of cells, protein in 15%, never both). No ANM inside | unfinished: 6-epoch scaffold; its recorded test Pearson (0.25–0.27) scores one random flow-matching sample per cell, and decoded directly it reaches 0.56 from RNA alone (phase 1: 0.61; [check](scripts/phase2_decode_check.py)); used only as the RNA-missing stand-in (Block 2, half of "both", reused in Block 4's RNA-missing view) |
| [Mode B](#t-modeb) | how ANM uses TEDDY | ANM as a decision layer on TEDDY's outputs: the phase-1 predictions, with labelled exceptions that use the phase-2 stand-in or measured protein | done: the six blocks (Block 6 under review) |
| [Mode A](#t-modea) | how ANM uses TEDDY | ANM inside TEDDY: its residual stream as ANM's field, layer Jacobians, in-silico gene perturbations | not done |

*Why two words:* phases are steps in building the model that supplies the evidence (phase 2 was meant to extend phase 1). Modes are not steps but two alternative ways for ANM to engage with TEDDY: on its outputs (B) or inside it (A). The axes are independent. Even unfinished, phase 2 already feeds Mode B as the RNA-missing stand-in (Blocks 2 and 4); a finished fusion model could be used in Mode B as a stronger second evidence channel that ANM weighs against TEDDY's, and that would still be Mode B, not Mode A.

**Two ways to study a foundation model**

| | [Mode B](#t-modeb) · this repo | [Mode A](#t-modea) · not done |
|---|---|---|
| What | Decide on top of TEDDY: TEDDY stays a closed box that supplies evidence | Look inside TEDDY: residual stream, layer Jacobians, in-silico gene perturbations |
| Answers | How should we decide from TEDDY's outputs? (edit the question, decline, explain, trust, label cost) | Why does TEDDY predict what it predicts, e.g. which genes drive its CD16 prediction? |
| Claimed here | Yes | No |

## Why ANM

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/why-dark.png">
  <img alt="Six results: 0 labels to change the question; 0 → 1,664 no-calls when RNA is missing and TEDDY can’t run; 16,433 explained calls; accuracy 0.945 → 0.998 when keeping only confident calls; 0 vs 50–500 labels to adopt a new question; false agreement 0.85 → 0.40" src="docs/assets/readme/why-light.png" width="100%">
</picture>

<details>
<summary>The same results as a table</summary>

| A biologist asks… | TEDDY alone | TEDDY + ANM |
|---|---|---|
| Can I change the question (soft → key-marker rule)? | Keeps its old rule; answers 1,752 cells the new question holds back | Written down with 0 labels; “no calls” follow the question: 317 → 834 → 3,373 |
| What if RNA is missing? | TEDDY can’t run; TEDDY alone’s fixed rule, applied to a weaker protein-only [stand-in](#t-adtonly), answers every cell (0 “no calls”, accuracy 0.63) | ANM on the same stand-in (no TEDDY), told to trust it at a declared 0.35: 1,664 “no calls”; still calls 86.8% of cells |
| Why this call? | No written-down rule to explain | [Leave-one-out](#t-loo) + [flip distance](#t-flip) on 16,433 cells ([permutation test](#t-permutation) p ≈ 0.0099) |
| Which calls can I trust? | Every call looks equally sure | Accuracy of TEDDY's own calls 0.945 → 0.998 keeping the 40% ANM is most [confident](#t-softp) about |
| What does a new question cost? | A [trained head](#t-train) needs 50–500 new labelled cells | 0 labels |
| Can it tell apart cells that look identical to TEDDY? | [False agreement](#t-falseagree) 0.85 on myeloid↔T [must-separate pairs](#t-mustpair) | 0.40 after [veto](#t-veto) → [complement](#t-complement) → [verify](#t-verify) |

Accuracy is not the claim. On Block 1’s held-out cells, compared like for like ([accuracy of calls made](#t-qdec)), ANM is within 0.02 of TEDDY alone (soft 0.950 vs 0.947, strict 0.972 vs 0.973, key-marker 0.960 vs 0.943). Other comparisons differ in either direction; for example, in Block 2 with both inputs under the key-marker rule the fixed rule scores 0.887 and ANM 0.758.

</details>

## How it works

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/how-dark.png">
  <img alt="Step 1: TEDDY stays frozen and hands its predictions on as typed evidence. Step 2: write down the question. Step 3: ANM decides per cell and explains. Compared against TEDDY alone and a trained head." src="docs/assets/readme/how-light.png" width="100%">
</picture>

A question (observer) is a short written declaration, not a trained model. The [strict rule](#t-o1) is the [soft rule](#t-o0) with three numbers changed:

```python
# bridge_anm/lib/lineage_panels.py → CRITERIA  (readable copy: bridge_anm/readouts/cite_lineage_O1.yaml)
"O1": {                          # strict rule
    "score_mode": "equal_panel_mean",
    "readout_threshold": 0.28,   # soft rule (O0): 0.12  → stricter calls
    "key_marker_boost": 2.0,     # soft rule (O0): 1.0   → CD19 / CD3 / CD16 count double
    "expected_margin": 0.12,     # soft rule (O0): 0.05
    ...
},
```

The three questions used throughout: **soft rule** ([O0](#t-o0)), **strict rule** ([O1](#t-o1)) and **key-marker rule** ([O2](#t-o2)).

## How we test "ANM helps TEDDY"

Every block follows the same protocol:

- **Same frozen TEDDY.** `best.pt` is never retrained; TEDDY alone and ANM read the same 9 predicted proteins per cell, and in Blocks 1 and 5 the trained heads also get the first 32 numbers of TEDDY's embedding. Exception: without RNA TEDDY can't run, so a protein-only [stand-in](#t-adtonly) supplies the evidence in Block 2's RNA-missing condition (and Block 4's RNA-missing view, which reuses it); Block 2's [both](#t-joint) condition averages TEDDY's predictions with it.
- **Same cells.** 16,750 [held-out cells](#t-holdout) from site4, plus an [out-of-site check](#t-ood) on 2,000 more.
- **Truth only checks answers.** The measured proteins ([ADT](#t-adt)) score the calls. Labelled exceptions also feed them into the evidence: Block 2's RNA-missing [stand-in](#t-adtonly) reconstructs the panel from them (half of [both](#t-joint) too; reused in Block 4's RNA-missing view), and Block 6's [complement](#t-complement) test uses them directly.
- **Decision metrics first.** Declining when unsure, label cost, per-marker reasons, coverage vs accuracy, false agreement. Accuracy is reported, never the win condition.

| Arm | What it is | Adopting a new question |
|---|---|---|
| TEDDY alone | Predicted proteins → a score per lineage → pick the best one if it clears a threshold (the soft rule) | Re-tune the rule on new labels |
| **TEDDY + ANM** | The same predictions as typed evidence → [ANM's decision engine](#t-finitefield) under a written-down question | **Edit the question: 0 labels, no retrain** |
| [Train + Jev](#t-train) | A classifier trained on labels, standing in for a Jev-class judge: logistic/MLP heads on TEDDY's 9 predictions plus 32 embedding numbers (Blocks 1, 5), heads on the 9 values only (Blocks 2, 4), and a re-tuned threshold grid (Block 5). Not TypeSafe AI's [Jev](#t-jev), which can't be trained on labels and was never called here | Collect labels and retrain |

## Results by block

| # | Block | What we do | Key result | Report |
|---|---|---|---|---|
| 1 | [Edit the question](https://danielchen26.github.io/teddy_mm/#edit) | Write the question three ways (soft → strict → key-marker); see who follows it without new labels | ANM “no calls” 317 → 834 → 3,373; trained head flat at 185; the key-marker rule changes 12.48% of expected answers | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 2 | [Missing modality](https://danielchen26.github.io/teddy_mm/#missing) | Take RNA away, so TEDDY can’t run ([RNA](#t-rnaonly) · [RNA missing](#t-adtonly) · [both](#t-joint)); declare how far to trust the protein-only stand-in (0.35) | RNA missing: “no calls” 0 → 1,664, 86.8% of cells still called; under the key-marker rule the 0.35 trust declines all 12,563 by construction | [MISSING_MODALITY](docs/reports/MISSING_MODALITY_ANM_DEMO.md) |
| 3 | [Attribution](https://danielchen26.github.io/teddy_mm/#attr) | Remove one marker's evidence at a time → the marker that decides → how far from flipping; [bootstrap](#t-bootstrap) 200×, permutation 100× | 16,433 cells · 27.7% of calls flip · p ≈ 0.0099 | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 4 | [Scope gate](https://danielchen26.github.io/teddy_mm/#scope) | Keep TEDDY's call only where ANM is confident ([soft_P](#t-softp) cut-off); compare hard cells with random ones | Accuracy 0.945 → 0.998 at 40% [coverage](#t-coverage); hard cells 0.810 vs 0.946 | [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) |
| 5 | [Zero-label transfer](https://danielchen26.github.io/teddy_mm/#labels) | Train heads on 50 → 10,619 labels until they copy ANM's key-marker answers | Copying needs ≈ 50 / 200 / 500 labels; ANM needs 0 | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 6 | [Decision loop](https://danielchen26.github.io/teddy_mm/anm-loop.html) | Flag look-alike pairs TEDDY's [embedding](#t-z512) can't separate → swap in measured protein as evidence, with the strict rule → check again | False agreement 0.85 → 0.40; accuracy 0.74 → 0.99; the [feature swap z⁺](#t-zplus) makes it worse (−0.27) | [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) |

Charts, data tables and every number's source: **[danielchen26.github.io/teddy_mm](https://danielchen26.github.io/teddy_mm/)**.

## Quickstart

**You need**

- Python 3.11 and the packages in `requirements.txt`
- TEDDY-G 70M weights at `../teddy_mwe/ckpt/teddy_g_70M/` (from [Merck/TEDDY](https://huggingface.co/Merck/TEDDY))
- ANM on `PYTHONPATH` (Active Neural Matter; the research code is currently private, so ask the authors for access)
- About 600 MB for the CITE-seq data (GEO GSE194122)

**1 · Environment**

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export PYTORCH_ENABLE_MPS_FALLBACK=1          # Apple silicon
```

**2 · Build TEDDY's evidence (once)** — TEDDY stays frozen; only a small RNA → protein head is trained

```bash
bash scripts/01_download_cite.sh               # GEO GSE194122, ~587 MB
python scripts/02_prepare_cite.py              # held-out site = site4
python scripts/03_embed_rna.py --device mps --seq-len 1024 --batch-size 16
python scripts/04_train.py --device mps
python scripts/05_eval.py --device mps
```

**3 · Hand the predictions on as typed evidence**

```bash
python bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000
```

**4 · Run the tests**

```bash
export PYTHONPATH=/path/to/ANM:.
python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100   # blocks 1 · 3 · 5
python bridge_anm/run_scope_refine_proof.py                                # block 4
python bridge_anm/export_missing_modality_events.py --n-cells 0            # block 2 (needs phase-2 ckpt)
python bridge_anm/run_missing_modality_demo.py
python bridge_anm/reverse_step1_must_separate.py                           # block 6
python bridge_anm/reverse_loop_small.py
```

Each script writes a Markdown report and a JSON file under `outputs/anm_cite_bridge/<block>/`; the published copies live in [`docs/reports/`](docs/reports/).

**5 · Ask a different question**

Edit a question's entry in `CRITERIA` ([`bridge_anm/lib/lineage_panels.py`](bridge_anm/lib/lineage_panels.py)): its threshold, key-marker boost or lineage weights. Keep its readable YAML copy in [`bridge_anm/readouts/`](bridge_anm/readouts/) in sync, then re-run the proof script. A new question id also needs adding to `CRIT_IDS` at the top of the proof script. TEDDY is not retrained and no labels are used; the held-out proteins only score the result.

## Key terms

<details open>
<summary>Every name used on this page and the site, in plain words (click to collapse)</summary>

| Term | In plain words |
|---|---|
| <a name="t-teddy"></a>**TEDDY** | Merck’s single-cell foundation model (TEDDY-G, 70M parameters). Here it reads a cell’s RNA and, with a small head trained in this repo, predicts the cell’s 134 surface proteins; 9 of them go to ANM. TEDDY itself is never retrained. |
| <a name="t-anm"></a>**ANM · Active Neural Matter** | A decision layer. It reads evidence, applies a written-down question (the observer) and returns a call or an honest “no call”, with checks attached. |
| <a name="t-citeseq"></a>**CITE-seq** | A technology that measures RNA and surface proteins in the same cell. The measured proteins are the answer key we hold out. |
| <a name="t-adt"></a>**ADT · measured surface protein** | CITE-seq’s protein readout (antibody-derived tags). Held out and used only to check answers, except where a test says otherwise. |
| <a name="t-holdout"></a>**Held-out cells (site4)** | 16,750 cells from a site never used in training. Their measured proteins score the answers; labelled tests also use them as input (Block 2’s RNA-missing stand-in, reused in Block 4’s RNA-missing view; Block 6’s complement test). |
| <a name="t-ood"></a>**Out-of-site check** | 2,000 cells from other sites (val_non_site4): a second check that the results hold elsewhere. |
| <a name="t-lineage"></a>**Lineage call** | The decision per cell: B cell, T cell or myeloid, each judged from 3 markers (B: CD19 CD72 CD22 · T: CD3 CD2 CD5 · myeloid: CD16 CD11c CD36). |
| <a name="t-evidence"></a>**Typed evidence (δu)** | A value for each of the 9 panel proteins, scaled 0–1 and handed to ANM as labelled inputs such as “CD19 is high”. Normally these are TEDDY’s predictions; in Block 2’s RNA-missing condition (reused in Block 4’s RNA-missing view) they are the protein-only stand-in’s reconstruction (no TEDDY), and in Block 6’s complement test the measured proteins. |
| <a name="t-observer"></a>**Observer · the question** | The question, written down: which markers count, how strong a signal must be, and when to decline. Changing it needs no labels and no retraining. |
| <a name="t-o0"></a>**Soft rule (O0)** | The default lineage call. All 9 markers weigh equally; call the best lineage if its score reaches 0.12. |
| <a name="t-o1"></a>**Strict rule (O1)** | The same question asked more strictly: key markers CD19, CD3 and CD16 count double and the score must reach 0.28. Same answers where both decide; declines more often. |
| <a name="t-o2"></a>**Key-marker rule (O2)** | A different question: only CD19 (B), CD3 (T) and CD16 (myeloid) count, with B and T weighted up. It changes 12.48% of the expected answers. |
| <a name="t-abstain"></a>**Decline to call (abstain)** | The readout says “no call” because the evidence does not meet the written-down rule: honest silence instead of a guess. |
| <a name="t-pf"></a>**Workability (P_f)** | 1 if the observer can make a call from the evidence, 0 if it declines. Averaged over cells, it is the share of cells called. |
| <a name="t-qf"></a>**Exactness (Q_f)** | 1 if a call matches the held-out protein truth. Averaged over labelled cells, a decline counts as a miss. |
| <a name="t-qdec"></a>**Accuracy of calls made (Q)** | Correct calls divided by calls made, checked against the held-out proteins. Used to compare arms like for like. |
| <a name="t-softp"></a>**Confidence score (soft_P)** | How strong ANM’s best option is for a cell. Keeping only high-score cells trades coverage for accuracy. |
| <a name="t-coverage"></a>**Coverage** | The share of cells that still get a call once a confidence cut is applied. |
| <a name="t-loo"></a>**Leave-one-out (LOO)** | Remove one protein’s evidence and decide again. The marker whose removal moves the score most is the reason for the call. |
| <a name="t-flip"></a>**Flip distance** | How much one marker’s value must change before the call switches lineage. Small means a fragile call. |
| <a name="t-bootstrap"></a>**Bootstrap** | Resample the cells with replacement (200 times) to see how stable a number is; gives a 95% interval. |
| <a name="t-permutation"></a>**Permutation test** | Shuffle the data many times to see how often a pattern this strong appears by chance. Here p ≈ 0.0099 over 100 shuffles. |
| <a name="t-train"></a>**Train + Jev · trained head (Jev-class stand-in)** | The “just train harder” control: a small classifier fitted on labelled cells, playing the role of a Jev-class judge (one probability per option, trained on the task’s own labels, as in the anm-jev repo). It is not TypeSafe AI’s Jev, which is zero-shot, cannot be trained on labels and was never called here. Its inputs depend on the block: logistic and MLP heads on TEDDY’s 9 predictions plus the first 32 embedding numbers (Blocks 1 and 5), heads on the 9 values only (Blocks 2 and 4), and in Block 5 also a re-tuned threshold grid on TEDDY alone’s rule. |
| <a name="t-jev"></a>**Jev (TypeSafe AI)** | TypeSafe AI’s “System One” decision model (jev-1.13), closed and reached only through its API. It answers typed questions (a Choice over named options, a Score, or yes/no) with calibrated probabilities and no rationale text. TypeSafe post-trains it (RLCD); it is used zero-shot and cannot be trained or fine-tuned on your labels. This repo never calls it. |
| <a name="t-z512"></a>**TEDDY embedding (z_512)** | TEDDY’s 512-number summary of a cell: the mean of its last-layer tokens at context length 1024. |
| <a name="t-mustpair"></a>**Must-separate pair** | Two cells TEDDY sees as near-identical (cosine ≥ 0.98 in its embedding) whose measured proteins disagree under the 3-lineage answer key, for example one keyed myeloid and one keyed T. Many of the “myeloid” cells in these pairs are NK cells, which the key has no class for. |
| <a name="t-veto"></a>**Veto** | The loop’s first check, a step outside ANM’s engine: it flags pairs that TEDDY’s embedding puts together (cosine ≥ 0.98) although their measured proteins disagree, so a readout built on the embedding should not be trusted there. |
| <a name="t-complement"></a>**Complement** | Give the decision the evidence it lacks (here, measured protein as typed evidence, in place of TEDDY’s predictions) instead of retraining TEDDY. |
| <a name="t-verify"></a>**Verify** | Score the same pairs again after the change: false agreement, soft separation and accuracy. |
| <a name="t-falseagree"></a>**False agreement** | Among must-separate pairs where both cells get a call, the share given the same call even though they differ. |
| <a name="t-softsep"></a>**Soft separation** | The share of must-separate pairs the readout tells apart: different calls, or one cell declined. |
| <a name="t-corrsep"></a>**Correct separation** | The share of must-separate pairs split into the right, different calls. |
| <a name="t-zplus"></a>**z⁺ · feature swap** | TEDDY’s embedding with protein features glued on: the “just swap features” alternative to typed evidence. |
| <a name="t-rnaonly"></a>**RNA (TEDDY)** | Evidence from TEDDY’s predictions, made from the cell’s RNA. |
| <a name="t-adtonly"></a>**RNA missing (protein-only stand-in)** | Without RNA, TEDDY cannot run on the cell. A second model from this repo stands in: the phase-2 bidirectional model, trained on TEDDY’s RNA embedding plus measured protein, here run with its RNA input switched off. It reads the cell’s 134 measured proteins and reconstructs the 9 panel proteins, poorly (panel Pearson 0.446). No TEDDY prediction is used for these cells. |
| <a name="t-joint"></a>**Both (averaged)** | The average of TEDDY’s prediction and the stand-in’s protein-only reconstruction. |
| <a name="t-finitefield"></a>**Finite field · ANM’s engine** | ANM’s decision engine: it combines the typed evidence under the observer into a score per lineage, then calls or declines. |
| <a name="t-pearson"></a>**Pearson ≈ 0.61** | How well TEDDY’s phase-1 head predicts all measured proteins (correlation on held-out cells). A reference point; we do not try to beat it. |
| <a name="t-phase1"></a>**Phase 1 · TEDDY + small head (done)** | Frozen TEDDY embeds each cell’s RNA once (512 numbers, context 1024). A small head trained in this repo maps that embedding to all 134 surface proteins, fitted with a negative-binomial loss. An MLP head beat a latent flow-matching head (test Pearson 0.610 vs 0.595); every block uses its predictions (labelled exceptions: the RNA-missing stand-in in Blocks 2 and 4, measured protein in Block 6’s complement test). |
| <a name="t-phase2"></a>**Phase 2 · fusion scaffold (unfinished)** | A multimodal model: TEDDY’s RNA embedding plus a protein encoder, trained with one modality randomly dropped (RNA in about 15% of cells, protein in 15%, never both) and a latent flow-matching decoder. Only a 6-epoch scaffold run exists. Its recorded test Pearson (0.25–0.27) scores one random flow-matching sample per cell; decoded directly it reaches 0.56 from RNA alone, still below phase 1’s 0.61 (scripts/phase2_decode_check.py). It contains no ANM; Block 2 uses it as the RNA-missing stand-in (and half of “both”), and Block 4’s RNA-missing view reuses it. |
| <a name="t-modeb"></a>**Mode B** | What this project does: treat TEDDY as an evidence source and study the decision layer on top of it. |
| <a name="t-modea"></a>**Mode A (not done)** | Studying TEDDY’s internals, such as its residual stream, layer Jacobians or in-silico gene perturbations. Not claimed here. |

</details>

## Repository layout

```text
teddy_mm/
├── scripts/                   # phase 1: download → prepare → embed (frozen TEDDY) → train head → eval
├── teddy_mm/                  # TEDDY encoder wrapper, RNA → protein models, phase-2 scaffold
├── bridge_anm/                # TEDDY → ANM bridge
│   ├── export_*.py            #   predictions → typed evidence
│   ├── run_*.py, reverse_*.py #   the six tests
│   ├── lib/lineage_panels.py  #   marker panels + CRITERIA (the questions the code runs)
│   ├── readouts/              #   readable YAML copy of each question
│   └── schemas/               #   ANM finite-field schema
├── configs/                   # phase-1 / phase-2 training configs
├── outputs/anm_cite_bridge/   # generated reports + JSON
└── docs/                      # GitHub Pages site, proof reports, figures
```

## Claim boundary

**Claimed ([Mode B](#t-modeb): TEDDY as an evidence source for ANM decisions)**

- TEDDY's predictions as typed evidence → ANM decisions, with questions that can be rewritten (soft, strict, key-marker) and no TEDDY retrain
- Declining to call when an input is missing; per-marker reasons (leave-one-out, flip distance)
- Confidence cut-offs, and zero-label adoption of a new question vs a trained head
- A decision loop on must-separate pairs: flag → add evidence → check again

**Not claimed**

- [Mode A](#t-modea): studying TEDDY's internals (residual stream, layer Jacobians, in-silico gene perturbations)
- The loop as TEDDY's representation responding to perturbations
- Pearson above phase-1 ≈ 0.61 · fusion audit · ATAC / chromatin · clinical superiority · anything about TypeSafe AI's Jev itself (never called here)

**Honest weak spots.** With RNA missing TEDDY can’t run: Block 2’s evidence then comes from a weak protein-only stand-in (panel Pearson 0.446), and its “no calls” follow a trust factor we declared (0.35). Block 4 trades coverage for accuracy, and under the key-marker rule accuracy dips at very low coverage. Block 6 uses measured protein as the evidence, in place of TEDDY's predictions, in one explicit test; in normal use it stays the answer key.

**TEDDY's embedding** (`z_512`) is the mean of its last-layer tokens at context length 1024 (not pretrain 2048, not a disease token), L2-normalized 512-D (`z_rna_512.npy`). Full scope: [`docs/MODE_B_SCOPE.md`](docs/MODE_B_SCOPE.md).

## Reports

[HARD_PROOF](docs/reports/HARD_PROOF.md) · [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) · [MISSING_MODALITY_ANM_DEMO](docs/reports/MISSING_MODALITY_ANM_DEMO.md) · [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) · [REVERSE_STEP1](docs/reports/REVERSE_STEP1.md) · [ANM_HELPS_TEDDY_LOOP](docs/reports/ANM_HELPS_TEDDY_LOOP.md) · [CHEAP_PROBES_RANKING](docs/reports/CHEAP_PROBES_RANKING.md) · [STAT_PROOF](docs/reports/STAT_PROOF.md) · [BAKEOFF_REPORT](docs/reports/BAKEOFF_REPORT.md)

Regenerate figures: `python scripts/make_story_infographics.py` (site posters) · `node scripts/render_readme_figures.mjs` (README figures, needs Playwright).

## References

- TEDDY: Merck TEDDY-G 70M on [Hugging Face](https://huggingface.co/Merck/TEDDY); paper [arXiv:2503.03485](https://arxiv.org/abs/2503.03485)
- Data: NeurIPS 2021 BMMC CITE-seq (Open Problems multimodal), GEO [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)
- ANM: Active Neural Matter (research code, currently private)

> **中文简介**：TEDDY 冻结不动（RNA→蛋白，Pearson ≈ 0.61）。ANM 在同一份冻结预测上加一层可编辑、会弃权、可审计的决策层。六个模块（改问题 · 缺模态 · 归因 · 范围门控 · 零标签迁移 · 决策闭环）都用同一套对照：TEDDY alone vs TEDDY + ANM vs 训练头，只证明“帮助 TEDDY”，不证明 Pearson 胜利，不重训，非临床。文中蓝色的名词都链接到上方 Key terms 的通俗解释；网站上则是虚线下划线，鼠标移上去即显示解释。

<details>
<summary><b>Phase-1 / phase-2 training notes (中文)</b></summary>

### TEDDY multimodal phase-1

冻结 TEDDY 编 RNA，预测 NeurIPS 2021 BMMC CITE 的表面蛋白。
对照：MLP vs latent flow matching。

为 Mac M4 Max 128GB（MPS）写的。不要在这台机器上预训练第二模态。

#### 一次性安装

```bash
cd teddy_mm
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export PYTORCH_ENABLE_MPS_FALLBACK=1
```

TEDDY-G 70M 权重应在 `../teddy_mwe/ckpt/teddy_g_70M/`（先跑 `teddy_mwe/setup.sh`）。

#### 正式数据（约 587 MB）

```bash
bash scripts/01_download_cite.sh
python scripts/02_prepare_cite.py
python scripts/03_embed_rna.py --device mps --seq-len 1024 --batch-size 16
python scripts/04_train.py --device mps
python scripts/05_eval.py --device mps
```

编码 9 万细胞、seq=1024，在 M4 Max 上大概数小时，只做一次。  
想先冒烟：

```bash
python scripts/02_prepare_cite.py --max-cells 4000
python scripts/03_embed_rna.py --seq-len 256 --batch-size 8 --device mps
python scripts/04_train.py --epochs 8 --device mps
```

没有 GEO 文件时可用合成数据测训练环：

```bash
python scripts/00_smoke_synthetic.py
python scripts/04_train.py --processed data/processed/cite_smoke --out outputs/smoke --epochs 5 --device cpu
```

#### 切分

- test = `site4`（数据集里的held-out site）
- val = 从剩余 donor 里抽 10%
- train = 其余

#### 现在做什么 / 还不做什么

做：RNA → ADT，冻结 TEDDY，MLP 基线 + FM。  
不做：ATAC foundation、解冻 400M、三模态联合生成。

双向 / 模态缺失是下一期：给 ADT 加 encoder，训练时 drop RNA。

#### Phase-2 bidirectional CITE (scaffold)

Independent of full phase-1 90k embed. ADT encoder (shallow MLP on the stored
ADT, transformed exactly once: `--adt-input-transform auto` = none for an
already-normalized pack such as `data/processed/cite`, log1p for integer counts)
+ modality dropout (~18% drop RNA) + latent FM with `cond = z_obs`. Val/test
decoding uses the train-median ADT size factor (`--eval-size-factor`), never the
evaluated cells' own measured ADT depth. `--legacy` restores the pre-fix
behaviour (clr on top of the stored ADT, measured size factor). In `adt_only` /
`joint` the scored cell's measured ADT is the encoder input, so those scores are
reconstructions. No ATAC / perturb-seq / 160M.

**Do not** write into `data/processed/cite/` or `outputs/cite_phase1/` while
phase-1 is running. Use a separate pack:

```bash
# Option A — tiny synthetic (immediate)
python scripts/00_smoke_synthetic.py
python scripts/06_train_bidirectional.py \
  --processed data/processed/cite_smoke --out outputs/cite_phase2 \
  --epochs 5 --device cpu

# Option B — 4k real CITE in a NEW dir (does not touch data/processed/cite)
python scripts/02_prepare_cite.py --max-cells 4000 --out data/processed/cite_phase2_4k
# Prefer CPU while phase-1 owns MPS:
python scripts/03_embed_rna.py --processed data/processed/cite_phase2_4k \
  --device cpu --seq-len 256 --batch-size 4
python scripts/06_train_bidirectional.py \
  --processed data/processed/cite_phase2_4k --out outputs/cite_phase2 \
  --epochs 8 --device mps   # or cpu if MPS still busy

# After full phase-1 z_rna exists, same script on the full pack:
# python scripts/06_train_bidirectional.py \
#   --processed data/processed/cite --out outputs/cite_phase2 --epochs 40 --device mps
```

Config: `configs/phase2_cite.yaml`. Code: `teddy_mm/bidirectional.py`,
`teddy_mm/modality_dropout.py`, `scripts/06_train_bidirectional.py`.

</details>
