# A trial that trains the last layer of TEDDY-G by the rows of the TEG: declaration

Written on 24 November 2026, **after the scripts were tried on made-up values and after a pilot that was looked at by
its training curve alone (Section 7), and before any state of a trained layer was read**.

**Amended the same day, still before any state of a trained layer was read or kept.** The first training made the
rows swing and was stopped by hand. What was changed is how the gradient of the loss is taken (Section 4, "Steps")
and one arm of part B (`cells5`, Section 5); Section 7 says what happened and what each pilot showed. The loss, the
learning rate, the limit of steps, the stopping rule, the rule of Section 6, the hypotheses and the forecasts are as
first committed.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit is read. One layer of the encoder is trained on the owner's machine; the owner agreed to the trial on
24 November.

## 1. Why this is asked

For IRF1 in THP-1 the model loses a step of the PD-L1 response between its input and the state Z a head reads (0.18
at the input view, 0.40 at Z). The loss was placed: the response does not get into the mean over the gene tokens,
and in the 70M checkpoint the last layer is the one that moves part of it there
([`encoder_stages_20261124`](../encoder_stages_20261124/)). Everything that needs no training was tried and none of
it closed the item: the annotation positions in the input, a head that reads the annotation position or the input
view beside Z, and the two larger checkpoints
([`annotation_input_20261124`](../annotation_input_20261124/),
[`annotation_beside_z_20261124`](../annotation_beside_z_20261124/),
[`input_beside_z_20261123`](../input_beside_z_20261123/), [`model_sizes_20261124`](../model_sizes_20261124/)).

Three numbers say where the item stands, all with the head reading Z:

| | *u* |
|---|---|
| a reader of Z fitted with IRF1's own cells (the reference) | 0.40 |
| the repaired output, made from the control cells and IRF1's one record (`rows`) | 0.68 |
| a readout with no record of IRF1 (`controls`) | 0.90 |

What only training can change is the last: **an event that is readable without a record of its own.** This is also
the measure by which a retraining with the TEG is to be read ([`GENERAL.md`](../../GENERAL.md), Section 7), and it is
the place where the TEG is used as the paper means it, as a constraint in training, and not only as a test.

Two questions:

1. Trained by the rows of the other events of its pathway, is IRF1 readable **without its own record**?
2. Trained by its own row as well, does the encoder's loss close, and does the repair of the head?

## 2. What is trained, and what is not

- **Trained:** every weight of the last encoder layer of TEDDY-G 70M (2.1 million of 71 million).
- **Frozen:** every other layer, the embeddings, the prediction head, and the readout used in training (below).
- The output of the layer before the last is computed once and kept (`scripts/cache.py`, tokens only). The trained
  weights are not kept; what is kept is the state of every cell under them.

## 3. The cells, and the two folds

The 7,478 cells of the THP-1 stratum: 2,386 control cells and the knockout cells of IFNGR1 (1,339), IFNGR2 (1,213),
IRF1 (1,072), JAK2 (1,024) and STAT1 (444). They are split into the two halves of the earlier rounds
(`anm_diag.cycle.halves`, seed 0; 3,739 and 3,739 cells).

Every arm is trained twice: fold 0 on the first half, fold 1 on the second. **A cell takes its state from the fold
that did not hold it** (`scripts/build_tt.py`). So a trained state that is read is the state of a cell the trained
layer never saw, neither its tokens nor its protein. In `scripts/train.py` the measured protein of the other half is
set to NaN before anything is computed.

## 4. The loss

**The readout *r*, fixed before training.** A ridge readout of the state (the mean over the gene tokens, at unit
length), fitted on the control cells of the fitted half at the state as released, by the recipe of the round's arm
`controls` (`anm_diag.cycle.LinearRows`). It reads the protein in units of its spread among those control cells, and
it holds no record of any event. It is not trained.

**The loss** has two terms, with weights 1 and 1; nothing is tuned.

- `rows`. For every event *e* of the arm's record: (the shift of *r* from the control cells to the cells of *e*,
  minus the measured shift)², averaged over the events. This is a row of the TEG written as a term. It uses the
  protein of a knockout cell only through the mean of its event.
- `hold`. The control cells stay where they are: the squared distance of a control cell's state (at unit length)
  from its state as released, in units of the spread of the control cells' states.

A degenerate solution of the rows is closed off by construction: the readout is fixed and was fitted on control
cells, so its slope among them cannot change sign, and the control cells are held.

**Steps** (amended, Section 7). Adam without momentum (betas 0 and 0.999), learning rate 1e-5. A step is 64 cells:
24 control cells and the rest in equal parts from the events of the record. The cells of a step are first read as
the layer stands. That gives the shift of every event in the step: the shift as released, read once on all cells of
the fitted half, plus the mean change of *r* among the step's cells of the event, minus the mean change among the
step's control cells. With these the gradient of the step's loss is exact, and four passes of 16 cells carry it to
the weights. Every 25 steps all cells of the fitted half are read. **Training stops when every row of the record is
read to within 0.03 control SD on the cells of the fitted half, or after 600 steps.** It is not run to a loss of
zero.

(As first committed: Adam with momentum, and the gradient of `rows` taken with a running value of the shifts, made
exact every 25 steps. That made the rows swing.)

## 5. The arms

What the record holds about the knockout cells of the fitted half:

| arm | stage | the record | what it is for |
|---|---|---|---|
| `rows4` | `T4` | the rows of IFNGR1, IFNGR2, JAK2, STAT1. No cell of IRF1 is used | **question 1**: IRF1 is an event without a record |
| `rows5` | `T5` | the same and the row of IRF1 | **question 2** |
| `mixed` | `TM` | as `rows5`, with the event labels permuted among the knockout cells | does the event of a row matter, or only "a knockout shifts" |
| `cells5` | `TC` | the measured protein of every cell of the fitted half, knockout and control, cell by cell, in place of the rows (amended: as first committed, of the knockout cells only) | ordinary fine-tuning with the same cells |
| `prop` | `TP` | `rows4`, and for the cells of IRF1 a row that was **not measured**: the mean of the four measured shifts | the simplest propagation: events of one type shift alike |

**Two parts.** The arms `rows4` and `rows5` are trained and read first (part A, `results/trial.json`): the three
hypotheses are judged there. The other three arms are trained after them and read with all five (part B,
`results/trial_B.json`): it adds the comparisons between arms. Each part is read once.

The state as released (`Z70`) is the arm with no training. The measured shifts, in control SD on all cells:
IFNGR1 −0.93, IFNGR2 −1.05, JAK2 −1.13, STAT1 −0.85, **IRF1 −0.36**. IRF1 shifts a third as far as the other
four, which is how it was found to differ from its type.

## 6. What is read, and the rule

The stages `V` (the input view), `Z70`, `T4`, `T5`, `TM`, `TC`, `TP` and the permuted view, for the same cells. The
grids are made by `anm_diag scope` with the SPEC files of the runs before (IRF1 on its own with its own cells,
linear reader and boosting; the four pathway genes as a pool, linear reader), and the round is read by
`anm_diag.cycle` with the cells, halves and resamples of the earlier rounds (`scripts/read.py`). Margin 0.10, step
margin 0.10, 1,000 resamples; every difference is a difference on paired resamples.

For IRF1, at every stage: the reference reader; its difference from the input view; the three outputs of the round
(`controls`: no record of IRF1 in the readout; `rows`: its one record; `cells`: every cell of the fitted half); the
gain of each over the same output at `Z70` (a step: the interval lies above 0.10; no step: below; unresolved
otherwise); the grade.

- **H1, without its own record.** In the arm `rows4`, where no cell of IRF1 was used in training, the output of a
  readout that has no record of IRF1 (`controls`) gains a step over `Z70`. Confirmed: a step. Refuted: no step.
  Unresolved otherwise. Whether the item is then closed is said beside it (the output supported, inside the margin).
- **H2, the encoder's loss with the record.** In the arm `rows5`, the reference reader of `T5` is no added loss from
  the input view (at `Z70` it is a loss: +0.21 [0.12, 0.32]). Confirmed: no added loss. Refuted: a loss. Unresolved
  otherwise.
- **H3, the repair with the record.** In the arm `rows5`, the repaired output (`rows`) gains a step over `Z70` and
  reaches a grade (relative or absolute). Confirmed: both. Refuted: no step. Unresolved otherwise.
- **The control target.** For the four pathway genes as a pool, at every trained stage: the repaired output does not
  lose a step against `Z70`, and the grade is kept.
- **The guard on the other outputs.** The prediction head as trained is applied to the released state and to every
  trained state. On the control cells its 133 other outputs change, per cell and in units of each output's spread
  among the control cells, by at most 0.10 on average. An arm for which this fails has changed what the head reads
  for unperturbed cells, and is reported as that.
- **Between arms**, on paired resamples, read both ways (brings back a step, costs a step, no change of a step,
  unresolved): `T4` against `T5` (part A); `T5` against `TM`, `T5` against `TC`, `T4` against `TP` (part B).
- **Checks.** With the head reading `Z70` the round repeats the first outer step for both targets (0.005). The state
  from what is kept equals the recorded state (5e-3 per value; 16-bit floats). If a check fails: stop and report.

Described and not judged (`results/trial_manifest.json`): the training curves, how far the states of cells that were
not held moved, and the shift of each fold's readout on those cells (no protein value enters).

What each outcome would mean:

| outcome | what follows |
|---|---|
| H1 confirmed | the records of other events, used as a constraint in training, reach an event that has none. The first thing made after the audit that does. It then needs units that have not been read |
| H1 refuted, H2 or H3 confirmed | training closes what the head could not, and only with the event's own record. The limit "no record, no reading" stands for this event, and what is missing is a type that covers it |
| nothing confirmed | one layer trained on one half of 7,478 cells is not enough to move the response into the mean. Said as that, with the training curves |
| a guard fails | the gain was paid for by other outputs. Reported with the gain, not as a repair |

## 7. What was done before this was written

- **The scripts were tried on made-up values** in a separate folder (a random model of three small layers and made-up
  cells): `cache.py`, `train.py` for every arm and fold, `build_tt.py`, the two grids, `read.py`, `score.py`. One
  rule was changed by it: H2 is not testable when the released state is no loss from the input view.
- **The cache** was made on the owner's machine for the 7,478 cells and checked against the recorded states
  ([`results/cache_manifest.json`](results/cache_manifest.json)).
- **A pilot, looked at by its training curve alone.** One arm and fold (`rows5`, fold 0) was trained with
  `train.py --pilot`, which writes the rows as read on the cells of the fitted half and nothing else: no state of a
  trained layer was computed for the other half, and none was kept. What it showed, and what was set by it
  (the curves: [`results/pilot.json`](results/pilot.json)):

  | learning rate | what the rows did on the fitted cells | control cells moved (in units of their spread) |
  |---|---|---|
  | 1e-4 (the first choice) | after 25 steps the four pathway rows were on the wrong side (+0.16 to +0.34); after 50 they had overshot (−1.04 to −1.38 for rows of −0.83 to −1.17). Stopped | 0.29, then 0.19 |
  | 3e-5 | 25 steps: overshot (−1.10 to −1.48); 50 steps: swung back (−0.52 to −0.62). Stopped | 0.043, then 0.016 |
  | **1e-5** | one overshoot at 25 steps, then settling: after 100 steps the four are at −0.89, −0.89, −1.04, −1.06 for rows of −0.93, −1.02, −1.17, −0.83 | 0.007 falling to 0.002 |

  The learning rate was set to 1e-5. **The limit of steps was raised from 300 to 600**: the row of IRF1 was learned
  slowly (−0.01 at the start, −0.13 after 100 steps, for a row of −0.40), at about 0.001 a step. Nothing else was
  changed. The pilot also says something the forecasts use: the four pathway genes move together, within 25 steps,
  and the cells of IRF1 follow them by about a tenth.
- **Speed** on the owner's machine: 37 cells a second for one training, 60 for two at once. Two folds are trained
  at once.

### What happened after the first commit, and the amendment

All of it before any state of a trained layer was computed for the cells of the other half, written or read. Every
line below is a training curve on the cells of the fitted half
([`results/pilot.json`](results/pilot.json), `after_the_first_commit`).

1. **The first training, as first committed, made the rows swing.** The arm `rows4`, both folds. In fold 0 the row of
   IFNGR1 (−0.93) was read at −0.89, −0.89, −1.02, −0.79, −1.16, −0.57 at steps 100 to 225, and the largest error of
   a row went from 0.15 to 0.50. The pilot of 100 steps had been too short to show it. Both folds were stopped by
   hand before step 250; the script then started the next arm, which was stopped 9 seconds later, before its first
   step. No state was written.
   Why it swung: the gradient of a row was taken with a running value of the shift, which lags, and Adam's momentum
   adds to the lag.
2. **Two replacements of the rows term were tried and dropped.** Both ask every knockout cell, one by one, to move by
   what its row is missing.
   - With the control cells entering only through `hold`, the rows did not move (IFNGR1 from −0.49 to −0.50 in 75
     steps): the layer moves knockout and control cells together, and the row is their difference.
   - Relative to the control cells of the step, and with the exact gradient of a step and no momentum, the curve was
     smooth and stopped short: after 150 steps at 1e-5 the largest error was 0.26 and no longer falling. This term
     asks more than a row does: that the cells of an event move alike.
3. **What is used: the loss as first committed, with its gradient taken exactly.** The cells of a step are read as
   the layer stands, so the shift of every event in the step is known without a running value (Section 4), and Adam
   is run without momentum. Arm `rows5`, fold 0, 200 steps:

   | learning rate | largest error of a row at steps 25, 50, … 200 | control cells moved |
   |---|---|---|
   | **1e-5** | 0.34, 0.31, 0.27, 0.24, 0.22, 0.20, 0.18, 0.17 | 0.001 |
   | 3e-5 | 0.27, 0.22, 0.16, 0.14, 0.14, 0.11, 0.15, 0.21 | 0.002 |

   At 1e-5 the error falls at every reading. At 3e-5 it falls faster and then the four pathway rows drift together
   (0.11 at step 150, 0.21 at step 200). **The learning rate stays 1e-5 and the limit stays 600 steps.** At step 200
   the rows are read at −0.86, −0.87, −0.23, −0.99, −0.91 for IFNGR1, IFNGR2, IRF1, JAK2, STAT1 (the rows: −0.93,
   −1.02, −0.40, −1.17, −0.83). It is not expected that every row is within 0.03 after 600 steps; an arm that is
   not is read as it is, and said (Section 8).
4. **The arm `cells5`** (part B) takes the measured protein of every cell of the fitted half, control cells
   included. As first committed it took the knockout cells only, which is the shape of the first replacement above
   that did not move. No training of this arm was run in either form.
5. **What was not changed:** the readout, the loss and its weights, the cells and folds, the arms of part A, the
   stopping rule, the rule of Section 6, the hypotheses, and the forecasts. The forecasts were written after the
   first pilot and were not revised after the later ones. One thing the later pilots show that the forecasts did not
   have: with the exact gradient the row of IRF1 is learned at about the same pace as before (−0.01 at the start,
   −0.16 after 100 steps, −0.23 after 200).

The training of the trial was started with this script while the amendment was being written, and the amendment was
committed before the first arm had finished, so before any state of a trained layer existed as a file.

Speed with two trainings at once is now 17 cells a second each, the readings of every 25 steps included (the cells
of a step are read once more than before): about 40 minutes for the two folds of an arm.

## 8. What would change the plan

- A check of Section 6 fails: stop and report.
- An arm does not reach the tolerance in 600 steps: it is read as it is, and said.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.

## 9. The forecasts

`results/declared.json`, 29 of them (23 for part A). They were written after the pilot and use what its training
curve showed. Briefly, and why:

- **H1 is more likely refuted than not** (refuted 0.50, unresolved 0.38, confirmed 0.12). In the pilot the cells of
  IRF1 followed the four pathway genes by about a tenth. A layer taught by their rows alone should move IRF1 little.
- **H2: unresolved is the most likely** (0.50). The reference reader has to come within 0.10 of the input view on
  cells the layer did not see, and the row of IRF1 was learned slowly.
- **H3: unresolved** (0.45): a step is likely, a grade is not.
- `mixed` is expected to be worse than `rows5` for IRF1: with the labels mixed no row says that IRF1 is an event.
- `prop` has a target for IRF1 that is three times the measured shift. A layer that learns it slowly may stop near
  the measured value, so its reading depends on how long it trained; its curve is reported with it.
- The control target is not harmed (0.85 for each rows arm) and the guard holds (0.80).

## 10. What this run cannot say

- One event in one kind of cell, on cells read before. A confirmation here is a lead for units that have not been
  read, not evidence.
- One layer of one checkpoint, trained on 3,739 cells. Nothing about a retraining of the encoder on its
  pre-training data.
- The propagation law of the paper's TEG: `prop` is its simplest form only.
