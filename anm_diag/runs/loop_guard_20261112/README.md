# Closing the loop within the guard — 12 November 2026

> Read again by intervals, with the heads compared pairwise, in [`runs/loop_reading_20261112`](../loop_reading_20261112/)
> (post hoc). The outcome reported in this file is by this run's own rule and stands as recorded.

Exploratory, and a follow-up to [`runs/loop_20261112`](../loop_20261112/). The design and the rule are in
[`DESIGN.md`](DESIGN.md), committed before any guarded head was fitted or audited (commit `e8157e0`). This file reports
what came out.

## Outcome by the rule: closed within the guard, on both cell types

The repair of `runs/loop_20261112` (only the head's MLP, the step the audit names) removed the loss and lowered the head's
score on bone marrow by more than the 0.02 the step-4 rule allows. Here the same repair carries one more term, which
holds the head's outputs on 3,000 bone-marrow training cells where they are. With it the head still adds no loss for PD-L1
on knockouts of genes held out of the fit, and the bone-marrow score stays within the guard.

PD-L1, pooled audit, same cells as in `runs/loop_20261112` (its two rows are repeated for comparison). Bone marrow: the
head's own evaluation on 5,000 test cells, lowest to highest of the heads of a row.

**Melanoma, IFN-γ** (1,171 knockout cells of 5 genes, 2,000 audit controls)

| head | Z | head latent | head output | root cause printed by the tool | added by the head | reading | bone marrow |
|---|---|---|---|---|---|---|---|
| as trained | 0.06 [-0.00, 0.11] unresolved | 0.22 [0.15, 0.27] reject | 0.45 [0.40, 0.49] reject | `relative_loss_at` Z → head latent | 0.39 | located in the head | 0.602 |
| `mlp_only` (no guard) | same | 0.10 [0.01, 0.18] unresolved | 0.12 [0.05, 0.17] unresolved | `reference_unresolved` | 0.06 | no loss in the head | 0.560 to 0.567 |
| `mlp_guarded` | same | 0.10 [0.05, 0.16] unresolved | 0.12 [0.06, 0.17] unresolved | `reference_unresolved` | 0.06 | no loss in the head | 0.603 to 0.604 |

**THP-1** (5,092 knockout cells of 5 genes, 1,193 audit controls)

| head | Z | head latent | head output | root cause printed by the tool | added by the head | reading | bone marrow |
|---|---|---|---|---|---|---|---|
| as trained | 0.12 [0.01, 0.29] unresolved | 0.39 [0.31, 0.52] reject | 0.54 [0.47, 0.63] reject | `relative_loss_at` Z → head latent | 0.42 | located in the head | 0.602 |
| `mlp_only` (no guard) | same | 0.18 [0.11, 0.32] reject | 0.16 [0.05, 0.30] unresolved | `reference_unresolved` | 0.04 | no loss in the head | 0.502 to 0.552 |
| `mlp_guarded` | same | 0.16 [0.07, 0.32] unresolved | 0.18 [0.07, 0.37] unresolved | `reference_unresolved` | 0.06 | no loss in the head | 0.594 to 0.595 |

Mean change of the five heads on bone-marrow test cells against the head as trained (0.602): +0.002 (melanoma
fits) and -0.007 (THP-1 fits); the single head furthest down is at +0.001 and -0.008. The guard allows
−0.020. By the rule written first: pooled PD-L1 "no loss in the head" and a mean fall of at most 0.02 → **closed within
the guard**, in both data sets. The melanoma-fitted heads score 0.001 to 0.002 above the head as trained, a difference too
small to read anything into.

As in the first run, the tool does not print "support" for the repaired head in the pooled audits. What is left at its
output is the size of what is unresolved at Z on these cells.

## Outside the rule

**The pilot that chose the weight** (first target held out; 2,000 bone-marrow validation cells; the head as trained scores
0.7431 on them):

| data set | weight | validation score | bone-marrow outputs moved (root mean square) | set-aside loss, first check to best | chosen |
|---|---|---|---|---|---|
| melanoma | 1 | 0.7426 | 0.031 | 0.0744 to 0.0456 | yes |
| melanoma | 10 | 0.7430 | 0.014 | 0.0787 to 0.0494 |  |
| melanoma | 100 | 0.7430 | 0.007 | 0.0895 to 0.0703 |  |
| THP-1 | 1 | 0.7407 | 0.063 | 0.4436 to 0.1752 | yes |
| THP-1 | 10 | 0.7426 | 0.030 | 0.4802 to 0.1892 |  |
| THP-1 | 100 | 0.7425 | 0.014 | 0.5408 to 0.2249 |  |

The smallest weight of the grid already keeps the validation score within 0.01 on both data sets, so 1 was used
for both. Without the term the same fit scored 0.698 (melanoma) and 0.629 (THP-1) on these cells. The set-aside
loss includes the holding term, so it is not comparable across weights.

**Held-out genes one at a time** (*u* at the head output and the reading; the first three columns are from
`runs/loop_20261112`):

| held-out gene | Z | as trained | `mlp_only` | `mlp_guarded` | bone marrow, this head |
|---|---|---|---|---|---|
| melanoma, IFNGR1 | 0.14 | 0.46 located in the head | 0.03 no loss in the head | 0.02 no loss in the head | 0.603 |
| melanoma, IFNGR2 | -0.02 | 0.52 located in the head | 0.13 not resolved | 0.13 not resolved | 0.604 |
| melanoma, JAK1 | -0.04 | 0.41 located in the head | -0.04 no loss in the head | -0.04 no loss in the head | 0.603 |
| melanoma, JAK2 | 0.18 | 0.50 located in the head | 0.11 no loss in the head | 0.16 no loss in the head | 0.603 |
| melanoma, STAT1 | 0.04 | 0.52 located in the head | 0.17 not resolved | 0.17 not resolved | 0.603 |
| THP-1, IFNGR1 | 0.01 | 0.55 located in the head | 0.14 not resolved | 0.14 not resolved | 0.595 |
| THP-1, IFNGR2 | 0.14 | 0.60 located in the head | 0.15 no loss in the head | 0.14 no loss in the head | 0.594 |
| THP-1, JAK2 | 0.15 | 0.64 located in the head | 0.21 no loss in the head | 0.23 no loss in the head | 0.595 |
| THP-1, STAT1 | -0.16 | 0.36 located in the head | 0.02 not resolved | 0.01 not resolved | 0.594 |
| THP-1, IRF1 | 0.33 | 1.12 located in the head | 0.60 not resolved | 0.60 not resolved | 0.594 |

With the guard the tool names a head step for 0 of the 10 held-out genes; 5 read "no loss in the head" and
5 "not resolved", gene for gene as without the guard. "No loss in the head" compares with Z: for THP-1 JAK2 the head
output taken alone is still rejected (0.23 [0.18, 0.29]), with Z at 0.15. IRF1 keeps a clear remainder (0.60).

**Other readouts** (pooled audit):

| readout | as trained | `mlp_only` | `mlp_guarded` |
|---|---|---|---|
| melanoma, HLA-E | output 1.11 reject; located in the head | output 0.91 reject; located in the head | output 1.04 reject; located in the head |
| melanoma, MHC-I | output 0.19 unresolved; not resolved | output -0.00 unresolved; not resolved | output -0.01 unresolved; not resolved |
| THP-1, PD-L2 | output 0.65 reject; located in the head | output 0.15 unresolved; no loss in the head | output 0.15 unresolved; no loss in the head |

**PD-L1 on bone marrow.** The correlation for PD-L1 alone on the test cells is 0.41 as trained and 0.40 to
0.41 (melanoma fits) and 0.40 to 0.41 (THP-1 fits) with the guard; without it, 0.20 to 0.27.

**The head's own PD-L1 output** falls for the held-out knockouts by 1.26 to 1.51 (melanoma) and 0.72 to 1.45
(THP-1) of the measured fall, in control standard deviations. The other 133 outputs moved on the audited cells by
0.06 (melanoma) and 0.18 to 0.19 (THP-1) root mean square, as without the guard.

**Floor check.** Audits with *u* at the head output more than 0.10 below *u* at Z: `target_frangieh2021_IFNGR1_mlp_guarded/PDL1`.

**Fits.** Kept step: 400 of at most 400. Every fit ran to the step limit with the set-aside loss still falling, so the repair is not converged.

## What this adds to the first run, and what it does not

- The cost seen in `runs/loop_20261112` came from how the repair was made, not from the located step: the MLP can carry the
  PD-L1 response to the output on these cells without moving the head's outputs on bone-marrow cells.
- The holding term keeps the head to its own earlier outputs. It can keep the original task, not improve it.
- Still not the sealed step 4: one seed, a fine-tune of one part on one protein, compared with the head as trained instead
  of a control retrain; the head latent is not at "support". Passing that rule needs the full training recipe, three
  seeds and a control arm.
- Transfer across cell types was not repeated here; in the first run it failed.

## Limits

- Designed after the first run's results were read. The grid of weights was set knowing how far the unguarded heads had
  moved the bone-marrow outputs. The weight was chosen on validation cells from sites the head was trained on; the guard is
  read on test cells from a site it never saw.
- The bone-marrow test cells had been used before, to score the unguarded heads; no choice was made on them.
- Five genes per data set as resampling units in the pooled audits (rule `K1`), 4 guides in the THP-1 single-gene audits,
  one split of the controls, no positive-control node (rule `C5`), two cell types, one trained readout.

## Files

12 folders `pooled_…_mlp_guarded`, `target_…_mlp_guarded`: `declaration.yaml`, `anm_record.json`, `anm_report.md`.
`results/loop_guard.json`: the audits, the readings, the bone-marrow falls and the rule. `results/fits.json`: the pilot
(every weight), the fits and the side measures. `results/bone_marrow.json`: the test scores of the ten heads.
`scripts/`: `fit_guarded.py` (plain NumPy; `--check` compares the gradient with finite differences; every fit asserts that
the decoder stays bit for bit and that the audit controls are those of the first run), `make_declarations.py`,
`run_one.py`, `run_all.py`, `score_bone_marrow.py`, `analyze.py`. Features, heads and the bone-marrow cells are not kept
in git (`ANM_GUARD_FEATURES`, `ANM_BM_FIT`, `ANM_BM_TEST`).

```
cd scripts
python make_declarations.py
python fit_guarded.py pilot frangieh2021 && python fit_guarded.py fit frangieh2021     # about 25 minutes on one core
python fit_guarded.py pilot papalexi2021 && python fit_guarded.py fit papalexi2021     # about 35 minutes
python fit_guarded.py features
python run_all.py                                                                      # 12 audits, about 12 minutes on two cores
python score_bone_marrow.py
python analyze.py
```
