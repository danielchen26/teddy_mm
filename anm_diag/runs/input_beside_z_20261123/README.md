# The head reads Z together with a view of the input

23 November 2026. Step 3 of the plan in [`GENERAL.md`](../../GENERAL.md): a repair in front of the head at small
scale. The head is given the 512 components of the model's own input beside the state Z it reads today, and the same
round is run twice, on the same cells, halves and resamples. The design, two hypotheses and eighteen forecasts were
committed before any grid or reading of this run existed ([`DESIGN.md`](DESIGN.md)).

Everything here is exploratory. The cells are development cells that earlier runs have read. No final unit was read,
no state was computed and no model was trained.

## What came out

**The loss before the head comes back when the head sees the input, and the repair gains a step from it. It does not
close.** Both hypotheses are unresolved.

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it.

| target | events | the repaired output, head reads Z | head reads Z and the input | gain | grade with the input |
|---|---|---|---|---|---|
| PD-L1, THP-1 | IRF1, on its own | 0.68 [0.57, 0.79] | 0.39 [0.26, 0.49] | **0.29 [0.23, 0.38]: a step** | not tested (gate 3) |
| MHC-I, melanoma, co-culture | five pathway genes, a pool | 0.36 [0.32, 0.41] | 0.26 [0.21, 0.30] | 0.10 [0.08, 0.12]: unresolved | relative, as before |
| control: PD-L1, THP-1 | four pathway genes | 0.09 [−0.04, 0.16] | −0.02 [−0.11, 0.05] | 0.07 [−0.08, 0.13]: unresolved | not tested (gate 3), as before |
| control: PD-L1, co-culture | five pathway genes | 0.11 [0.01, 0.18] | 0.07 [−0.04, 0.14] | 0.04 [−0.03, 0.07]: no step | not tested (gate 3), as before |
| control: PD-L1, IFN-γ | five pathway genes | 0.04 [−0.03, 0.10] | −0.05 [−0.12, 0.01] | −0.01 [−0.10, 0.09]: no step | **not tested (gate 3); it was relative** |

**IRF1 (H1: unresolved).**

- The place of the loss is confirmed by an intervention. A reader of Z with IRF1's own cells leaves 0.40
  [0.27, 0.50]; a reader of the joint stage leaves 0.21 [0.07, 0.32]. The gain, 0.19 [0.12, 0.26], is a step. The
  joint stage reads what the input view reads on its own (0.18; the difference is 0.03, no added loss). So what the
  encoder dropped for this event is back in front of the head, with no retraining.
- The round's repair gains a step: 0.68 to 0.39. It is still rejected, it adds 0.18 [0.10, 0.27] beyond the reader of
  the joint stage (unresolved), and gate 3 does not pass (on the cells it was fitted on it leaves 0.32). So no grade is reached.
- What is short is the record, not the stage. Fitted on every cell of the fitted half, the same readout leaves 0.06
  [−0.11, 0.18]: no longer rejected, and no worse than the reference. One record of about 540 cells moves 1,024
  weights only part of the way.
- Without IRF1's record the wider stage gives nothing: a readout of control cells leaves 1.00 against 0.90 at Z, and a
  reader that never saw IRF1's cells leaves 0.82 against 0.64. The input view helps an event the repair has a
  record of; it does not carry to one it has not.

**MHC-I in co-culture (H2: unresolved).** The repaired output falls from 0.36 to 0.26. The gain, 0.10 [0.08, 0.12],
lies across the step margin. The reader of the joint stage leaves 0.23, the same as boosting on the input view alone.
The target was relative before and is relative now; it is not supported at either stage.

**The controls.** By the rule written first, none loses a step (3 of 3). One thing was not asked by that rule and is
reported: for PD-L1 under IFN-γ the repair at the joint stage overshoots on the cells it was fitted on (*u* = −0.15
[−0.22, −0.09]), gate 3 no longer passes, and the reference is no longer supported (0.05 [−0.03, 0.11] against 0.03
[−0.04, 0.10]). On the half it was not fitted on the output is where it was. The grade "relative" of that target is
therefore not kept with the wider stage. That target was the one closest to absolute.

**The check.** The readings with the head reading Z repeat those of the first outer step for all five targets
(largest difference of an arm 0.0000).

## What it means for the plan

- **Step 3 is done at small scale, and its answer is partial.** For the one item the loop had placed between the
  input and Z, the lost part is recoverable in front of the head: retraining the encoder is not needed to make it
  readable. It is not enough to close the item, because a single record does not fit a readout of twice the width.
- **It is not a general fix.** Given to every type, the wider stage costs the best target its grade. In the loop's
  own terms it is a change of a view for the types where a loss was placed before the head, not a new head for all.
- **For step 5 (retraining with the TEG) this says where the bar is.** A retrained encoder that keeps what IRF1 needs
  in Z would put the head where the joint stage puts it today: a reader with the event's own cells at about 0.2, and
  a record-based repair at about 0.4. To do better than this run, the retraining has to make the event readable
  without its own record, which is what the arm `controls` measures (1.00 here).
- Nothing here is a test. The cells have been read before, and IRF1 is one event.

## The forecasts

Intervals: 4 of 10 held the value. Outcomes: for 6 of 8 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)). Every miss is on the same side: the gain was larger than forecast.

| missed | forecast | observed |
|---|---|---|
| IRF1: the repaired output at the joint stage | [0.40, 0.70] | 0.39 |
| IRF1: fitted on every cell, at the joint stage | [0.15, 0.40] | 0.06 |
| IRF1: gain of the repair | [0.00, 0.25] | 0.29 |
| IRF1: gain when fitted on every cell | [0.03, 0.27] | 0.36 |
| IRF1: the gain of the repair against the step margin | unresolved 0.55 | a step (given 0.20) |
| MHC-I, co-culture: gain of the repair | [−0.02, 0.10] | 0.10 (0.103) |
| H2 | refuted 0.60 | unresolved (given 0.35) |
| PD-L1, IFN-γ: the repaired output at the joint stage | [−0.02, 0.10] | −0.05 |

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: the joint stage, the targets, the hypotheses and their rules, the forecasts |
| [`results/round.json`](results/round.json) | both rounds for every target, every arm, the gains on paired resamples, the input view alone |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/build_zv.py` | the arrays with the joint stage and the SPEC files; no readout enters |
| `scripts/round.py`, `score.py` | the two rounds and their differences; the score |

To run again (the work folders of `teg_loop_20261120` and `outer_step_20261122` are needed):

    python scripts/build_zv.py --work WORK
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/round.py --work WORK --observe OBSERVE
    python scripts/score.py

## What is next

- More events of IRF1's kind (a factor downstream of the pathway that the encoder drops) would say whether the record
  is what is short. The unread units of the pool hold none that passes admission for PD-L1 (round on batch 1).
- A repair that may use the input view only for types where the loop placed a loss before the head is a change of the
  package (views), to be declared by an outer step before it is read.
