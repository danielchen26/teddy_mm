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

Bridge sidecar `outputs/anm_cite_bridge/z_rna_512.npy` stores **L2-normalized full** `z_512` (export row order). Compact `z_rna_export.npy` keeps only `z_keep=32` leading dims for typed-event JSONL — **not** the reverse-step1 probe space.

---

## Implemented (Mode B)

- TEDDY RNA→ADT preds as typed δu evidence into ANM finite_field
- Editable observer YAML: `bridge_anm/readouts/cite_lineage_O{0,1,2}.yaml`
- Missing-modality abstain demos (honest silence under weak channels)
- Closed-form LOO / top-1 flip attribution on **declared decisions**
- Scope gate / zero-label observer transfer vs Train+Jev (see proof reports)
- Mode B reverse **decision** loop on must-separate pairs (veto→complement→verify)

Verified numbers only — cite [`HARD_PROOF`](reports/HARD_PROOF.md), [`SCOPE_REFINE_PROOF`](reports/SCOPE_REFINE_PROOF.md), [`MISSING_MODALITY_ANM_DEMO`](reports/MISSING_MODALITY_ANM_DEMO.md).

## NOT claimed (Mode A and adjacent)

- Mode A: residual stream as field, layer Jacobian, in-silico gene perturbs
- Gated ±ε / ±ε/2 G1–G4 protocols
- Mode-A reverse closed loop (representation response to perturbs)
- Fusion audit
- ATAC / chromatin foundation

## What current conclusions *are* about

Current proofs speak to **decision-layer sensitivity to TEDDY evidence** (edit observer → abstain / P_f / Q_f / LOO flips), **not** to TEDDY representation response under gene/protein perturbs.

They **cannot** yet answer whether `z_512` (or the export sidecar) is **sufficient for protein readout** in the strong sense. A reverse step-1 “must-separate” probe (near-identical z, different true lineage/protein labels) is a **sufficiency-for-readout** diagnostic — still Mode B / evidence geometry, **not** Mode A perturb-response.

## Honest weak-model note

**ADT-only** is a weak channel in the missing-modality demo (panel Pearson drops; TEDDY-alone abstain can stay 0 while Q is mediocre). The abstain story is intentionally honest: ANM silence under weak evidence beats silent over-answer. Not clinical. Do not claim Pearson > phase-1 full-ADT **~0.61**.

## Cheap next probes (Mode B only)

Four cheap probes ranked in [`docs/reports/CHEAP_PROBES_RANKING.md`](reports/CHEAP_PROBES_RANKING.md). Still **not** Mode A. Gated ±ε G1–G4 here are **decision-layer** nudges on typed evidence, not residual-stream / gene-perturb protocols.


## Reverse decision loop (Mode B — implemented)

Script: [`bridge_anm/reverse_loop_small.py`](../bridge_anm/reverse_loop_small.py)

**Veto → complement → verify** on reverse_step1 `z_512` must-pairs. Winner = true-ADT typed evidence + O1 (not `z+` alone). Deep writeup + NotebookLM posters: [`docs/reports/ANM_HELPS_TEDDY_LOOP.md`](reports/ANM_HELPS_TEDDY_LOOP.md). Still **not** Mode A perturb-response.

## Reverse step 1 (sufficiency probe only)

Script: [`bridge_anm/reverse_step1_must_separate.py`](../bridge_anm/reverse_step1_must_separate.py)

Reads full `z_rna_512.npy` (not compact `z_keep=32`) + true ADT panel, finds near-identical z pairs with different lineage/protein labels (**must-separate**). Report: [`docs/reports/REVERSE_STEP1.md`](reports/REVERSE_STEP1.md).

Still Mode B / evidence geometry — **not** Mode A perturb-response.

