# Mode B scope (claim boundary)

**Mode B (implemented in this repo):** TEDDY as a **typed evidence source** → ANM **finite_field** decisions on CITE holdout (site4/test). Editable observer YAML (O0/O1/O2); missing-modality abstain demos; LOO/flip attribution on **decisions**.

This is **not** Mode A. Do not read current results as residual-stream / Jacobian / in-silico gene-perturb claims.

---

## z_512 loading factor (required definition)

| Item | Value |
|---|---|
| Symbol | `z_512` |
| Source | Frozen TEDDY-G last-layer hidden states |
| Pool | **Mean-pool over tokens** (attention mask), not a disease token |
| Context length | **1024** (CITE embed; see `scripts/03_embed_rna.py --seq-len 1024`) |
| Not | TEDDY pretrain context **2048** |
| Not | A dedicated disease / CLS token readout |
| Dim | `d_model = 512` after mean-pool |

Bridge sidecar `outputs/anm_cite_bridge/z_rna_export.npy` stores **L2-normalized** z with `z_keep=32` (first 32 dims) for compact export — a **loading factor / evidence index**, not a claim that 32-D is biologically complete.

---

## Implemented (Mode B)

- TEDDY RNA→ADT preds as typed δu evidence into ANM finite_field
- Editable observer YAML: `bridge_anm/readouts/cite_lineage_O{0,1,2}.yaml`
- Missing-modality abstain demos (honest silence under weak channels)
- Closed-form LOO / top-1 flip attribution on **declared decisions**
- Scope gate / zero-label observer transfer vs Train+Jev (see proof reports)

Verified numbers only — cite [`HARD_PROOF`](reports/HARD_PROOF.md), [`SCOPE_REFINE_PROOF`](reports/SCOPE_REFINE_PROOF.md), [`MISSING_MODALITY_ANM_DEMO`](reports/MISSING_MODALITY_ANM_DEMO.md).

## NOT claimed (Mode A and adjacent)

- Mode A: residual stream as field, layer Jacobian, in-silico gene perturbs
- Gated ±ε / ±ε/2 G1–G4 protocols
- Reverse closed loop (representation response to perturbs)
- Fusion audit
- ATAC / chromatin foundation

## What current conclusions *are* about

Current proofs speak to **decision-layer sensitivity to TEDDY evidence** (edit observer → abstain / P_f / Q_f / LOO flips), **not** to TEDDY representation response under gene/protein perturbs.

They **cannot** yet answer whether `z_512` (or the export sidecar) is **sufficient for protein readout** in the strong sense. A reverse step-1 “must-separate” probe (near-identical z, different true lineage/protein labels) is a **sufficiency-for-readout** diagnostic — still Mode B / evidence geometry, **not** Mode A perturb-response.

## Honest weak-model note

**ADT-only** is a weak channel in the missing-modality demo (panel Pearson drops; TEDDY-alone abstain can stay 0 while Q is mediocre). The abstain story is intentionally honest: ANM silence under weak evidence beats silent over-answer. Not clinical. Do not claim Pearson > phase-1 full-ADT **~0.61**.

## Reverse step 1 (sufficiency probe only)

Script: [`bridge_anm/reverse_step1_must_separate.py`](../bridge_anm/reverse_step1_must_separate.py)

Reads `z_rna_export.npy` + true ADT panel, finds near-identical export-z pairs with different lineage/protein labels (**must-separate**). Report: [`docs/reports/REVERSE_STEP1.md`](reports/REVERSE_STEP1.md).

Still Mode B / evidence geometry — **not** Mode A perturb-response.

