# Closing the loop within the guard: can the located step be repaired without the cost? Design, written before the runs

Status: exploratory, and a follow-up. It was designed after the results of [`runs/loop_20261112`](../loop_20261112/) were
read. This file, the scripts and the 12 declarations were committed before any guarded head was fitted beyond a 20-step
timing check on training cells, and before any audit of this folder was run. The commit history shows the order.

## Why this run

`runs/loop_20261112` repaired the step the audit names (the head's MLP) and, by its rule, closed the loop on both cell
types: the head no longer adds a loss for PD-L1 on knockouts of genes held out of the fit. The same repair lowered the
head's score on bone-marrow test cells, its original task, from 0.602 to 0.560–0.567 (melanoma fits) and 0.502–0.552
(THP-1 fits). The step-4 rule sealed on 9 November
([`declarations/teddy_mm_step4_intervention_rule.md`](../../declarations/teddy_mm_step4_intervention_rule.md)) allows
0.02 and says of a larger fall that "the change trades one failure for another and is not credited". So the problem was
located and removed, and another one was made.

The question here: with the same repair held in place on bone-marrow cells, does the loop still close, and does the
bone-marrow score stay within 0.02?

## The repair

`mlp_guarded`: the `mlp_only` repair of `runs/loop_20261112` with one more term in the loss.

- Same cells, split, audit controls (checked against `loop_20261112/results/audit_controls_*.json`), held-out genes,
  target, optimiser, learning rate, 400 steps, stopping rule. Only the three MLP layers change; the decoder stays bit for
  bit (asserted).
- **Holding term.** 3,000 cells drawn at random (seed 0) from the training split of the bone-marrow data the head was
  trained on. The term is the mean squared change of all 134 outputs on these cells, against the head as trained, times a
  weight. A fifth of the 3,000 is set aside and enters only the stopping check.
- **The weight** is chosen before any audit, without audit cells and without bone-marrow test cells. Grid 1, 10, 100. For
  each data set one fit per weight, with the first target (IFNGR1) held out. Take the smallest weight whose score on 2,000
  bone-marrow validation cells (the head's own validation split) is at most 0.01 below the head as trained on the same
  cells. If none is, take 100 and say so. The four other held-out fits use that weight.

## The audits (12) and the guard

- Pooled, per data set (2): the pooled declaration of `runs/loop_20261112`, same cells; every knockout cell through the
  head fitted with its gene held out, audit controls assigned to the five heads as there. All readouts.
- Single head (10): one held-out gene, its head, all audit controls; PD-L1.
- **Guard:** the head's own evaluation (mean over 134 proteins of the correlation between measured counts and the predicted
  mean) on the 5,000 bone-marrow test cells used in `runs/loop_20261112`, for each of the five heads of a data set.

## Reading and rule

Readings as in `runs/loop_20261112` ("located in the head", "no loss in the head", "not resolved"). Per data set, with
the head as trained "located in the head" there (it was, on these cells):

| pooled PD-L1, `mlp_guarded` | mean fall of the five heads on bone-marrow test cells | outcome |
|---|---|---|
| no loss in the head | at most 0.02 | **closed within the guard** |
| no loss in the head | more than 0.02 | closed, guard broken |
| anything else | at most 0.02 | guard kept, not closed |
| anything else | more than 0.02 | neither |

## Also reported, outside the rule

- The pilot: for each weight, the validation score, how far the bone-marrow outputs moved, the fitting curves.
- The single-head audits; the readouts not trained on (PD-L2, HLA-E, MHC-I); the floor check.
- The largest fall among the five heads, and the correlation for PD-L1 alone on bone-marrow test cells.
- Next to every number, the head as trained and the unguarded `mlp_only` on the same cells, from `runs/loop_20261112`.

## Audit of this plan

| problem | what the design does about it |
|---|---|
| Designed after seeing that the unguarded repair closes the loop and breaks the guard. | Graded exploratory and named a follow-up. The rule, the grid and the selection are fixed here before any guarded audit. |
| The weight could be tuned until the test passes. | Three weights, a selection rule that reads only bone-marrow validation cells, one fit per weight. The test cells and the audit cells are read once, for the chosen weight. |
| The grid was set with knowledge: on the bone-marrow training sample the unguarded heads moved the outputs by 0.19 (melanoma) and 0.42 (THP-1) root mean square and scored 0.698 and 0.629 on the validation cells (0.743 as trained). | Stated. The grid spans two orders of magnitude around the weight at which the holding term is as large as the PD-L1 term. |
| The bone-marrow test cells were already used to score the unguarded heads. | No guarded head has been scored on them. They are not used for any choice. |
| Validation cells come from sites the head was trained on, test cells from a site it never saw; 0.01 on validation may not be 0.02 on test. | The rule reads the test cells; the outcome "closed, guard broken" covers it. |
| Holding the head to itself cannot make it better on bone marrow, only keep it. | The guard is a fall against the head as trained, as is this comparison. |
| A weight large enough to keep the bone-marrow outputs may leave the MLP no room to carry the response. | That is the question. "Guard kept, not closed" would mean that this head cannot have both with this repair. |
| The decoder-only arm is not repeated. | Its test was made in `runs/loop_20261112`; a holding term can only restrict it further. |
| Not the sealed step 4: one seed, a fine-tune, the head as trained as comparator instead of a control retrain; five genes per data set as resampling units (rule `K1`), one split of the controls, no positive-control node (rule `C5`), two cell types, one trained readout. | Stated with the results. |
