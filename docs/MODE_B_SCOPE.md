# Mode B scope (claim boundary)

**Mode B (implemented in this repo):** TEDDY as a **typed evidence source** → ANM **finite_field** decisions on CITE holdout (site4/test). Editable observer YAML (O0/O1/O2); missing-modality abstain demos; LOO/flip attribution on **decisions**.

**Scope in one sentence.** teddy_mm runs ANM's open-loop decision layer: TEDDY is a prescribed evidence source and does not read the field, so ANM's dynamics (state feedback, operator memory, delays) are not tested here.

This is **not** Mode A. Do not read current results as residual-stream / Jacobian / in-silico gene-perturb claims.

> **Corrections in progress (29 Sep 2026).** See the notice at the top of the [README](../README.md#corrections-in-progress-29-sep-2026) and the [site](https://danielchen26.github.io/teddy_mm/#corrections). In short: markers enter the engine one per step in panel order (a hidden weight; rerun with simultaneous markers pending); every prediction was scaled by a per-cell factor computed from measured protein (rerun with a train-median factor pending); the answer key has no NK or out-of-scope class; Experiment 6's headline and the "why this call" headline (CD5 23.6%) are withdrawn pending reruns; TEDDY's own margin gates calls as well as ANM's confidence, so Experiment 4 claims no ANM-specific gate; "0 labels for a new question" holds for any written rule, including TEDDY's fixed rule. Disclosures: donor 15078 is in training at other sites and is 32.6% of the site4 test cells; the "out-of-site check" is held-out donor 18303 at a training site, which was also the validation set; the 9 panel proteins were chosen by test-set Pearson.

---

## z_512 loading factor (required definition)

| Item | Value |
|---|---|
| Symbol | `z_512` |
| Source | Frozen TEDDY-G last-layer hidden states |
| Pool | **Mean-pool over tokens** (attention mask), not a disease token |
| Context length | **512** (CITE embed; see `scripts/03_embed_rna.py --seq-len 512`; the stored embeddings reproduce at 512 with cosine 1.000000, vs median 0.970 at 1024) |
| Not | TEDDY pretrain context **2048** |
| Not | A dedicated disease / CLS token readout |
| Dim | `d_model = 512` after mean-pool |

Bridge sidecar `outputs/anm_cite_bridge/z_rna_512.npy` stores **L2-normalized full** `z_512` (export row order). Compact `z_rna_export.npy` keeps only `z_keep=32` leading dims for typed-event JSONL — **not** the reverse-step1 probe space.

---

## Implemented (Mode B)

- TEDDY RNA→ADT preds as typed δu evidence into ANM finite_field
- Editable observer YAML: `bridge_anm/readouts/cite_lineage_O{0,1,2}.yaml`
- Missing-modality abstain demos (honest silence under weak channels)
- Closed-form LOO / top-1 flip attribution on **declared decisions**
- Scope gate / zero-label observer transfer vs TEDDY + trained classifier (see proof reports)
- Mode B reverse **decision** loop on must-separate pairs (veto→complement→verify)

Verified numbers only — cite [`HARD_PROOF`](reports/HARD_PROOF.md), [`SCOPE_REFINE_PROOF`](reports/SCOPE_REFINE_PROOF.md), [`MISSING_MODALITY_ANM_DEMO`](reports/MISSING_MODALITY_ANM_DEMO.md).

## NOT claimed (Mode A and adjacent)

- Mode A: residual stream as field, layer Jacobian, in-silico gene perturbs
- ANM's dynamics: TEDDY does not read the field, so state feedback, operator memory and delays are not exercised
- Gated ±ε / ±ε/2 G1–G4 protocols
- Mode-A reverse closed loop (representation response to perturbs)
- Fusion audit
- ATAC / chromatin foundation

## What current conclusions *are* about

Current proofs speak to **decision-layer sensitivity to TEDDY evidence** (edit observer → abstain / P_f / Q_f / LOO flips), **not** to TEDDY representation response under gene/protein perturbs.

They **cannot** yet answer whether `z_512` (or the export sidecar) is **sufficient for protein readout** in the strong sense. A reverse step-1 “must-separate” probe (near-identical z, different true lineage/protein labels) is a **sufficiency-for-readout** diagnostic — still Mode B / evidence geometry, **not** Mode A perturb-response.

## Honest weak-model note

**ADT-only** is a weak channel in the missing-modality demo (panel Pearson drops; the fixed rule's abstain can stay 0 while Q is mediocre). The abstain story is intentionally honest: ANM silence under weak evidence beats silent over-answer. Not clinical. Do not claim Pearson > phase-1 full-ADT **~0.61**.

## Cheap next probes (Mode B only)

Four cheap probes ranked in [`docs/reports/CHEAP_PROBES_RANKING.md`](reports/CHEAP_PROBES_RANKING.md). Still **not** Mode A. Gated ±ε G1–G4 here are **decision-layer** nudges on typed evidence, not residual-stream / gene-perturb protocols.


## Reverse decision loop (Mode B — implemented; Experiment 6 headline withdrawn pending redesign)

Script: [`bridge_anm/reverse_loop_small.py`](../bridge_anm/reverse_loop_small.py)

**Veto → complement → verify** on reverse_step1 `z_512` must-pairs. Winner = true-ADT typed evidence + O1 (not `z+` alone). Deep writeup + NotebookLM posters: [`docs/reports/ANM_HELPS_TEDDY_LOOP.md`](reports/ANM_HELPS_TEDDY_LOOP.md). Still **not** Mode A perturb-response.

## Reverse step 1 (sufficiency probe only)

Script: [`bridge_anm/reverse_step1_must_separate.py`](../bridge_anm/reverse_step1_must_separate.py)

Reads full `z_rna_512.npy` (not compact `z_keep=32`) + true ADT panel, finds near-identical z pairs with different lineage/protein labels (**must-separate**). Report: [`docs/reports/REVERSE_STEP1.md`](reports/REVERSE_STEP1.md).

Still Mode B / evidence geometry — **not** Mode A perturb-response.
