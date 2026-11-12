# Closing the loop: repairing the step the audit names — 12 November 2026

Exploratory. The design, the rule and an audit of the plan are in [`DESIGN.md`](DESIGN.md), committed before any run
(commit `7b1f4c8`). This file reports what came out.

## Outcome by the rule: the loop closed at the located step, on both cell types

The audit of the head as trained names the step from TEDDY's embedding (Z) to the head's latent, that is the head's first
network (the MLP), as the place where the PD-L1 response is lost. Two repairs were fitted on the same cells with the same
loss and budget, and audited on knockouts of genes they had never seen and on control cells no fit had used.

PD-L1, pooled audit. *u* is the share of the measured response that cannot be read at a stage; "added by the head" is *u*
at the head output minus *u* at Z.

**Melanoma, IFN-γ** (1,171 knockout cells of 5 genes, 2,000 audit controls)

| head | Z | head latent | head output | root cause printed by the tool | added by the head | reading |
|---|---|---|---|---|---|---|
| as trained | 0.06 [-0.00, 0.11] unresolved | 0.22 [0.15, 0.27] reject | 0.45 [0.40, 0.49] reject | `relative_loss_at` Z → head latent | 0.39 | located in the head |
| `mlp_only` (prescribed) | same | 0.10 [0.01, 0.18] unresolved | 0.12 [0.05, 0.17] unresolved | `reference_unresolved` | 0.06 | no loss in the head |
| `decoder_only` (ruled out) | same | same as trained | 0.43 [0.37, 0.49] reject | `relative_loss_at` Z → head latent | 0.37 | located in the head |

**THP-1** (5,092 knockout cells of 5 genes, 1,193 audit controls)

| head | Z | head latent | head output | root cause printed by the tool | added by the head | reading |
|---|---|---|---|---|---|---|
| as trained | 0.12 [0.01, 0.29] unresolved | 0.39 [0.31, 0.52] reject | 0.54 [0.47, 0.63] reject | `relative_loss_at` Z → head latent | 0.42 | located in the head |
| `mlp_only` (prescribed) | same | 0.18 [0.11, 0.32] reject | 0.16 [0.05, 0.30] unresolved | `reference_unresolved` | 0.04 | no loss in the head |
| `decoder_only` (ruled out) | same | same as trained | 0.47 [0.41, 0.58] reject | `relative_loss_at` Z → head latent | 0.35 | located in the head |

By the rule written first: head as trained "located in the head", `mlp_only` "no loss in the head", `decoder_only` not
→ **loop closed at the located step**, in both data sets.

What "closed" does and does not say here. The loss the head added (0.39 and 0.42) is gone after the prescribed repair
(0.06 and 0.04), and the tool no longer names a step. The tool does not print "support" for the repaired head: what is
left at its output (0.12 and 0.16) is the size of what is already unresolved at Z on these cells, and in THP-1 the head
latent taken alone is still above the margin (0.18).

## Outside the rule

**Held-out genes one at a time** (30 audits; *u* at the head output and the reading):

| held-out gene | knockout cells | Z | as trained | `mlp_only` | `decoder_only` |
|---|---|---|---|---|---|
| melanoma, IFNGR1 | 198 | 0.14 | 0.46 located in the head | 0.03 no loss in the head | 0.44 located in the head |
| melanoma, IFNGR2 | 278 | -0.02 | 0.52 located in the head | 0.13 not resolved | 0.51 located in the head |
| melanoma, JAK1 | 248 | -0.04 | 0.41 located in the head | -0.04 no loss in the head | 0.40 located in the head |
| melanoma, JAK2 | 244 | 0.18 | 0.50 located in the head | 0.11 no loss in the head | 0.50 located in the head |
| melanoma, STAT1 | 203 | 0.04 | 0.52 located in the head | 0.17 not resolved | 0.46 located in the head |
| THP-1, IFNGR1 | 1,339 | 0.01 | 0.55 located in the head | 0.14 not resolved | 0.43 located in the head |
| THP-1, IFNGR2 | 1,213 | 0.14 | 0.60 located in the head | 0.15 no loss in the head | 0.48 located in the head |
| THP-1, JAK2 | 1,024 | 0.15 | 0.64 located in the head | 0.21 no loss in the head | 0.53 located in the head |
| THP-1, STAT1 | 444 | -0.16 | 0.36 located in the head | 0.02 not resolved | 0.29 located in the head |
| THP-1, IRF1 | 1,072 | 0.33 | 1.12 located in the head | 0.60 not resolved | 1.01 located in the head |

- As trained: the tool names a head step for 10 of 10 genes. After `decoder_only`: still 10 of 10. After `mlp_only`: 0 of
  10; 5 read "no loss in the head" and 5 "not resolved" (the output is more than 0.10 above Z but no step is named). In
  four of those five *u* at the output is 0.02 to 0.17. IRF1 is the one that clearly remains: 1.12 falls to 0.60.
- "No loss in the head" compares with Z. For THP-1 JAK2 the head output taken alone is still rejected with `mlp_only`
  (0.21 [0.16, 0.28]), while Z is at 0.15.
- With about 200 knockout cells the tool does not always name the same step as the pooled audit. For the melanoma head
  as trained it names Z → head latent for one gene, head latent → head output for two, and Z → head output (no single
  step) for two. The prescription was taken from the pooled audit; repairing the MLP removed the loss for all five.

**The average at the endpoint would have passed the wrong repair.** After `decoder_only` the head's own PD-L1 output falls for
the held-out knockouts by about as much as the measured protein (melanoma 1.15 to 1.29 of the measured fall in control
standard deviations, THP-1 0.49 to 1.43; as trained 0.41 to 0.45 and the wrong direction). The audit still rejects that
head at its output (0.43 and 0.47), because the latent it reads from is unchanged. After `mlp_only` the output falls by
1.26 to 1.52 (melanoma) and 0.77 to 1.47 (THP-1) of the measured fall, so both repairs overshoot the average.

**The decoder arm was not cut short.** Every `decoder_only` fit stopped by its own rule after 80 to 200 steps, when the
loss on the set-aside cells no longer fell. Every `mlp_only` fit ran to the limit of 400 steps with that loss still
falling, so the prescribed repair is not converged.

**Readouts the repair was not trained on** (pooled audit):

| readout | measured change | as trained | `mlp_only` | `decoder_only` |
|---|---|---|---|---|
| THP-1, PD-L2 (the head has no PD-L2 output) | -0.31 | latent 0.47 reject, output 0.65 reject; located in the head | latent -0.12, output 0.15, both unresolved; no loss in the head | output 0.53 reject; located in the head |
| melanoma, HLA-E | -0.08 | latent 0.74, output 1.11, both reject; located in the head | latent 0.91, output 0.91, both reject; located in the head | output 1.19 reject; located in the head |
| melanoma, MHC-I | -0.14 | output 0.19 unresolved; not resolved | output -0.00 unresolved; not resolved | output 0.09 unresolved; not resolved |

Repairing the MLP on PD-L1 also removed the loss for PD-L2 in THP-1, which was never in the fit. It did not for HLA-E in
melanoma, whose measured change is small and whose intervals are wide.

**Transfer: none, as expected.** The head repaired on one cell type and audited on all cells of the other is still
"located in the head" for PD-L1 with both repairs: THP-1 cells with the melanoma-fitted `mlp_only` head 0.40 at the output
(0.47 as trained on the same cells in `runs/split_20261111`), melanoma cells with the THP-1-fitted one 0.43 (0.46).

**Cost on the head's original task.** 5,000 test cells of the bone-marrow data the head was trained on; mean over 134
proteins of the correlation between measured counts and the predicted mean (the head's own evaluation; the head as
trained scores 0.602 on these cells and 0.603 on the full test set).

| heads (6 per row) | mean over proteins | PD-L1 alone |
|---|---|---|
| as trained | 0.602 | 0.41 |
| `mlp_only`, fitted on melanoma | 0.560 to 0.567 | 0.20 to 0.21 |
| `mlp_only`, fitted on THP-1 | 0.502 to 0.552 | 0.22 to 0.27 |
| `decoder_only`, fitted on melanoma | 0.579 to 0.595 | 0.36 to 0.39 |
| `decoder_only`, fitted on THP-1 | 0.364 to 0.491 | -0.07 to 0.08 |

Both repairs make the head worse on bone marrow, and the prescribed one halves the correlation for PD-L1 there. The
repaired heads show that the located step is the one to change. They are not heads to use: they were fitted on one protein
in one cell type, without the original training cells in the loss. The other 133 outputs moved on the audited cells by
0.06 (melanoma) and 0.19 (THP-1) root mean square with `mlp_only`, 0.02 to 0.05 with `decoder_only`.

**Floor check.** In one of the 40 audits *u* at the head output is more than 0.10 below *u* at Z: melanoma IFNGR1 with
`mlp_only`, 0.03 against 0.14, intervals [-0.21, 0.09] and [-0.10, 0.23]. It is a single-head audit, so pooling of heads
is not the cause, and the held-out cells and audit controls were in no fit. The likely cause is the reference: with about
2,200 cells the readout models learn less from the 512 numbers of Z than from a latent already fitted to PD-L1 on other
cells, so *u* at Z is too high in small audits (also seen in `runs/exposure_20261111`).

## Follow-up

[`runs/loop_guard_20261112`](../loop_guard_20261112/) repeats the prescribed repair with a term that holds the head's outputs
on bone-marrow training cells. By its rule, written before it ran, the loop then closes within the 0.02 guard on both cell
types.

## Against the step-4 rule sealed on 9 November

The repository holds an earlier and stricter rule for closing the loop
([`declarations/teddy_mm_step4_intervention_rule.md`](../../declarations/teddy_mm_step4_intervention_rule.md)): the head
latent at bounded support in every one of at least three seeds, against a control arm retrained with the same recipe, and
the mean correlation on the bone-marrow test cells no more than 0.02 below. This run is not that step: it fine-tunes one
part of the existing head on one protein, with one seed, and compares with the head as trained. Measured against that
rule the prescribed repair would not be credited. The head latent is unresolved (melanoma) or still rejected (THP-1)
rather than supported, and the bone-marrow score falls by 0.035 to 0.042 (melanoma fits) and 0.050 to 0.100 (THP-1
fits). In that rule's words, the change trades one failure for another.

## Limits

- Five genes per data set are the resampling units in the pooled audits (rule `K1`: read directions, not the intervals);
  the single-gene audits resample 17 to 28 guides in melanoma and 4 in THP-1. One split of the controls. No positive-control node (rule `C5`).
- Two cell types, both IFN-γ-pathway knockouts, one trained readout. The pooled declarations, taken from
  `runs/split_20261111`, pool five genes without declaring it (rule `S2`).
- "No loss in the head" is defined against Z on the same cells. Where Z is itself unresolved (both pooled audits) the
  repaired head is not shown to be right, only to add no loss beyond Z.
- The designer knew an earlier run from another workspace (4b, not in this repository) in which retraining the MLP on the
  melanoma cells removed the loss and retraining the decoder did not. New here: genes held out, a second cell type, the rule
  written first, records in the repository.
- The repair does not transfer across cell types and costs accuracy on bone marrow. A repair meant to be kept would
  retrain the MLP on the original cells together with perturbation cells, and be audited again, bone marrow included.

## Files

40 folders `pooled_…`, `target_…`, `transfer_…`: `declaration.yaml`, `anm_record.json`, `anm_report.md`.
`results/loop.json`: every audit, the readings and the rule. `results/fits.json`: the 24 fits (curves, kept step) and, per
held-out gene, the shift of the head's own PD-L1 output next to the measured one. `results/audit_controls_*.json`: the
control cells no fit used. `results/cost_bone_marrow.json`: the bone-marrow scores of all heads.
`scripts/`: `fit_heads.py` (plain NumPy; every fit asserts that the part not being repaired stays bit for
bit; `--check` runs one short fit per part and writes nothing), `make_declarations.py`, `run_one.py`, `run_all.py`, `analyze.py`,
`cost_bone_marrow.py`. Features and heads are written to `ANM_LOOP_FEATURES` and are not kept in git.

```
cd scripts
python make_declarations.py
python fit_heads.py fit frangieh2021     # 12 fits, about 15 minutes on one core
python fit_heads.py fit papalexi2021     # 12 fits, about 30 minutes
python fit_heads.py features
python run_all.py                        # 40 audits, about an hour on two cores
python analyze.py
ANM_BM_TEST=<bone-marrow test cells> python cost_bone_marrow.py
```
