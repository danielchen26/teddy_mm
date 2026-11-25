# A trial that trains the last layer of TEDDY-G by the rows of the TEG

24 November 2026. For IRF1 in THP-1 the 70M checkpoint loses a step of the PD-L1 response between its input and the
state a head reads, and everything that needs no training was tried without closing the item. This run trains for
the first time: the last encoder layer of TEDDY-G 70M (2.1 million of 71 million weights), on the owner's machine,
with the rows of the TEG as the loss. The loss, the arms, the rule, the guards and 29 forecasts were committed before
any state of a trained layer existed ([`DESIGN.md`](DESIGN.md); it was amended once, still before any such state,
and says what was changed and why).

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit was read. One layer was trained, for 600 steps, on half of 7,478 cells at a time. Five arms, in two parts:
part A (`rows4`, `rows5`) by which the three hypotheses are judged, and part B (`mixed`, `cells5`, `prop`).

## What came out

**Trained by the rows, the last layer makes one readout read IRF1 on cells it never saw, and little else.**

- **By the declared reading**, with IRF1's own row in the training that is a step: 0.90 to 0.62 with no record of
  IRF1 at the reading, 0.68 to 0.44 with one. With no cell of IRF1 in the training it is 0.12 more, unresolved.
- **It does not bring back what the encoder loses.** The reference reader stays at 0.40, against 0.18 at the input.
- **The gain belongs to the vector the layer was trained under.** Kept exactly as it was, that readout reads about
  twice the gain. Fitted again on other control cells, a readout reads a third to a half of it.
- **Which event a row belongs to did not matter.** With the event labels shuffled the result for IRF1 is the same.

By the rule: **H1 unresolved, H2 refuted, H3 unresolved.** The guard holds in all five arms.

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. A cell takes its
state from the layer that was trained on the other half of the cells, so every state that is read is the state of a
cell the layer never saw, neither its tokens nor its protein.

**IRF1, on its own.** The input view reads 0.18 [0.04, 0.30].

| the output | as released | `rows4` | `rows5` | `mixed` | `cells5` | `prop` |
|---|---|---|---|---|---|---|
| a readout with no record of IRF1 (`controls`) | 0.90 [0.78, 1.07] | 0.78 [0.66, 0.91] | **0.62 [0.51, 0.73]** | 0.65 [0.55, 0.76] | 0.76 [0.64, 0.89] | 0.32 [0.18, 0.44] |
| the same with IRF1's one record (`rows`) | 0.68 [0.57, 0.79] | 0.58 [0.47, 0.69] | **0.44 [0.31, 0.55]** | 0.47 [0.35, 0.57] | 0.55 [0.44, 0.66] | 0.20 [0.03, 0.33] |
| fitted on every cell of the fitted half (`cells`) | 0.42 [0.30, 0.53] | 0.37 [0.24, 0.48] | 0.25 [0.09, 0.37] | 0.28 [0.13, 0.39] | 0.32 [0.19, 0.44] | 0.12 [−0.07, 0.26] |
| the tool's reference reader, with IRF1's own cells | 0.40 | 0.40 | 0.40 | 0.40 | 0.41 | 0.41 |
| the head as trained, reading the state | 0.95 | 0.95 | 0.95 | 0.95 | 0.95 | 0.95 |

What each arm trains with: `rows4` the rows of the four pathway genes, and no cell of IRF1; `rows5` the same and the
row of IRF1; `mixed` as `rows5` with the event labels shuffled among the knockout cells; `cells5` the measured
protein of every cell, one by one; `prop` as `rows4`, and for IRF1 a row that was not measured, the mean of the four
others.

| against the released state, on paired resamples | `rows4` | `rows5` | `mixed` | `cells5` | `prop` |
|---|---|---|---|---|---|
| no record of IRF1 | +0.12 [0.09, 0.17]: unresolved | **+0.28 [0.22, 0.37]: a step** | +0.25 [0.20, 0.33]: a step | +0.14 [0.10, 0.19]: a step | +0.58 [0.45, 0.76]: a step |
| its one record | +0.10 [0.07, 0.14]: unresolved | **+0.24 [0.19, 0.31]: a step** | +0.21 [0.16, 0.28]: a step | +0.13 [0.09, 0.18]: unresolved | +0.48 [0.37, 0.63]: a step |
| every cell | +0.05 [0.02, 0.09]: no step | +0.17 [0.12, 0.24]: a step | +0.14 [0.10, 0.20]: unresolved | +0.10 [0.05, 0.14]: unresolved | +0.30 [0.20, 0.39]: a step |
| the reference reader | −0.01 [−0.02, 0.01]: no step | −0.01 [−0.03, 0.01]: no step | −0.00 [−0.02, 0.01]: no step | −0.02 [−0.04, 0.01]: no step | −0.02 [−0.09, 0.06]: no step |

| between arms | no record of IRF1 | its one record |
|---|---|---|
| `rows4` against `rows5` | +0.16 [0.12, 0.21]: costs a step | +0.14 [0.11, 0.19]: costs a step |
| `rows5` against `mixed` | −0.03 [−0.06, −0.00]: no change of a step | −0.03 [−0.06, −0.00]: no change of a step |
| `rows5` against `cells5` | −0.14 [−0.20, −0.09]: unresolved | −0.11 [−0.16, −0.06]: unresolved |
| `rows4` against `prop` | +0.46 [0.35, 0.61]: costs a step | +0.38 [0.29, 0.51]: costs a step |

**The hypotheses** (part A):

- **H1, without its own record: unresolved.** In `rows4` no cell of IRF1 was used anywhere in the training, and the
  readout holds no record of it. It gains 0.12 [0.09, 0.17]: a step is not shown and "no step" is not shown either.
  The output is still rejected (0.78), so the item is not closed.
- **H2, the encoder's loss with the record: refuted.** The reference reader of `rows5` is a loss from the input
  view, +0.22 [0.12, 0.33], as it is at release (+0.21). The reader did not move (−0.01 [−0.03, 0.01]).
- **H3, the repair with the record: unresolved.** The repaired output gains a step (0.68 to 0.44) and comes to
  within 0.04 [−0.03, 0.10] of the reference reader. It reaches no grade: gate 3 is not passed, so the grade is
  not tested.

**The three arms of part B:**

- **`mixed`: which event a row belongs to did not matter.** The labels of the knockout cells are shuffled, so every
  "event" is a mix of the five and every row is the mean shift of a knockout cell (−0.83 to −0.96). IRF1 reads 0.65
  against 0.62 with its own row: no change of a step. The layer did not learn five rows. It learned to move a cell
  that looks knocked out along the readout, and it moves the cells of IRF1 less than the others whatever it is
  told.
- **`cells5`: the rows do at least as much as the protein of every cell.** Ordinary fine-tuning with the measured
  protein of the same cells, one by one, gains 0.14; the five rows gain 0.28. The difference, 0.14
  [0.09, 0.20], is not resolved as a step. In 600 steps `cells5` had learned its rows less far (to within 0.24 and
  0.17 on the fitted cells).
- **`prop`: a row that was not measured moves the event, and the value here is not evidence for the rule.** IRF1
  reads 0.32 with no record and no protein of IRF1 anywhere. But the derived row asks for −0.98, which is 2.7 times
  the measured shift (−0.36), because IRF1 differs from the genes of its type: that is why it was flagged in the
  first place. The training stopped at 600 steps on its way there (−0.61 and −0.68 on the fitted cells, and still
  moving), so the reading lands where the training stopped: with the round's readout IRF1 is under-read (0.32), and
  with the readout of the training it is already a little over-read (−0.09, below). Trained to its row, IRF1 would
  be read far beyond what was measured. What the arm shows is the mechanism: a row derived by a rule can be trained
  in. Whether the value is right is the rule's burden, and "the mean of its type" is wrong for IRF1 by a factor of
  2.7.

**The control target** (the four pathway genes as a pool; their rows are in the training of every arm):

| the output | as released | `rows4` | `rows5` | `mixed` | `cells5` | `prop` |
|---|---|---|---|---|---|---|
| no record of these events | 0.46 [0.39, 0.50] | 0.28 [0.22, 0.33] | 0.29 [0.22, 0.34] | 0.26 [0.16, 0.33] | 0.27 [0.20, 0.32] | 0.37 [0.31, 0.42] |
| with their records (`rows`) | 0.09 [−0.04, 0.16]: unresolved | 0.01 [−0.10, 0.09]: supported | 0.02 [−0.10, 0.09]: supported | −0.00 [−0.14, 0.08]: unresolved | −0.02 [−0.14, 0.05]: unresolved | 0.12 [0.01, 0.19]: unresolved |
| the reference reader | 0.08 | 0.09 | 0.09 | 0.09 | 0.07 | 0.10 |

Without a record the pool gains a step in four arms (+0.17 to +0.19) and 0.08, unresolved, in `prop`. With the
records the grade is kept in every arm. "Not harmed" is shown for `rows4`, `rows5` and `prop`; for `mixed` and
`cells5` the interval of the change reaches −0.107, so a loss of a step is not excluded, by 0.007.

**The guard holds in all five arms.** The head as trained, applied to the trained states: its 133 other outputs
change on the control cells by 0.019, 0.024, 0.018, 0.038 and 0.052 of their spread per cell (`rows4`, `rows5`,
`mixed`, `cells5`, `prop`); the limit is 0.10. The head's own PD-L1 output does not change (IRF1: 0.95 before and
after every arm), which is expected: it was not the readout of the training, and it reads none of the response.

**The checks hold.** With the head reading the released state the round repeats the first outer step for both
targets (largest difference 0.0000). The head as trained repeats its recorded output at the released state (7e-6).
Part B repeats part A in every value they share.

**The training**, described and not judged ([`results/trial_manifest.json`](results/trial_manifest.json)). Every
arm and fold ran the 600 steps; none read every row to within 0.03, so each is read as it is. At the stop the rows
were read, on the fitted cells, to within 0.09 and 0.05 (`rows4`), 0.12 and 0.09 (`rows5`), 0.06 and 0.04
(`mixed`), 0.24 and 0.17 (`cells5`), 0.37 and 0.31 (`prop`, where the derived row of IRF1 is the one not reached).
The row of IRF1 was learned in part in `rows5`: −0.30 for a row of −0.40 in one fold, −0.25 for −0.31 in the
other. The states hardly move: in units of the control cells' spread, 0.001 for control cells and 0.003 to 0.004
for the pathway knockouts in `rows4` and `rows5`; at most 0.004 and 0.010 in `prop`.

## A check added after the reading: three readouts

It was not declared and nothing is judged by it ([`scripts/seam_check.py`](scripts/seam_check.py),
[`results/seam_check.json`](results/seam_check.json)). It was added because two readers of the same state disagreed
in part A: the reference reader did not move, and the round's arm fitted on every cell read a step more.

The trial has a seam. The two halves of the cells have their states from two different trained layers, and the round
and the tool both fit a reader on one half and read the other. So there are three readouts to tell apart, all for
the output with no record of the event:

1. **The readout of the training, exactly.** The vector a layer was trained under, kept as it was, applied to the
   cells that layer did not see. Nothing is fitted again.
2. **The round's readout: the declared reading.** The same recipe on the same control cells, fitted again on the
   states those cells have in the trial, which come from the other layer and differ from the released ones by a
   thousandth of their spread.
3. **A fresh readout.** The same recipe on other control cells: inside one half, where every state comes from one
   layer, the half is split again and a readout fitted on one part reads the other. Half the cells, so the released
   values are higher (0.99 and 0.94 for IRF1).

The gain over the released state, IRF1 with no record:

| arm | 1. the readout of the training, exactly | 2. the round's (declared) | 3. a fresh one, first half and second half |
|---|---|---|---|
| `rows4` | +0.23 [0.18, 0.31]: a step (leaves 0.67) | +0.12 [0.09, 0.17]: unresolved (0.78) | +0.05, +0.05: no step |
| `rows5` | +0.53 [0.43, 0.69]: a step (0.37) | +0.28 [0.22, 0.37]: a step (0.62) | +0.08 [0.07, 0.11]: unresolved, +0.14 [0.10, 0.23]: a step |
| `mixed` | +0.50 [0.39, 0.64]: a step (0.41) | +0.25 [0.20, 0.33]: a step (0.65) | +0.08, +0.10: unresolved |
| `cells5` | +0.33 [0.25, 0.43]: a step (0.58) | +0.14 [0.10, 0.19]: a step (0.76) | +0.04: no step, +0.07: unresolved |
| `prop` | +0.81 [0.56, 0.96]: a step (−0.09) | +0.58 [0.45, 0.76]: a step (0.32) | +0.17, +0.23: a step |

The pool of pathway genes, no record: +0.30 to +0.37 with the readout of the training (it leaves 0.09 to 0.16;
`prop` +0.07), +0.17 to +0.19 with the round's, and +0.06 to +0.08 with a fresh one in every arm (`prop` +0.04 to
+0.06).

- **Each time the readout is fitted again, about half of the gain is gone.** The layer was trained under one
  vector, and what it learned is tied to that vector. The loss asks only that the cells of an event move along it.
- **Under the readout it was trained with, the gain is large and it is real on new cells**: 0.90 to 0.37 for IRF1
  with its row in the training, and 0.90 to 0.67, a step, with no cell of IRF1 used. This is the reading that fits
  a system whose head is kept. It was not the declared one, and the hypotheses are not judged by it.
- **With a fresh readout almost nothing is left without IRF1's cells** (+0.05, no step, in both halves), and for
  the pool every arm gains the same 0.06 to 0.08, whatever its loss.
- **The seam does hide something from the reference reader**, and it is small: inside one layer a reader fitted
  with IRF1's own cells gains 0.06 and 0.10 after `rows5`, where across the seam it gains nothing. What is lost
  between the input and the state is 0.21.
- **The first write-up of part A had this wrong.** It called the round's readout "nearly the readout the layer was
  trained under". The two differ by a factor of two in what they read.

## What it means

- **The encoder's loss is not within reach of the last layer.** H2 is refuted, and the check does not change that:
  at most 0.10 of 0.21 comes back to a reader with IRF1's own cells, in one half. This agrees with where the loss
  was placed: the response is already missing from the mean of the gene tokens when it enters the first layer
  ([`encoder_stages_20261124`](../encoder_stages_20261124/), [`model_sizes_20261124`](../model_sizes_20261124/)).
- **The layer learned a vector, not the rows.** Three things say so: the gain halves each time the readout is
  fitted again; shuffling the event labels changes nothing; and the states move by a few thousandths of their
  spread. One layer, 600 steps and a loss on one projection are enough to bend the state toward that projection
  for cells that look knocked out, and not enough to make the response readable in general.
- **For a system whose head is kept, that is still a gain, and this trial did not test it on the system's head.**
  The vector of this trial was a readout fitted for the trial. The head as trained was left alone and reads IRF1 at
  0.95 before and after. The direct follow-up is to fix the head as trained as the readout and train the layer
  under it; then the output that is audited is the system's own. It is a new training, and the owner decides.
- **"No record, no reading" stands for IRF1 by the declared reading**, narrowly: +0.12, unresolved. Under the
  readout of the training it is a step (+0.23) and the output is still rejected (0.67). With a fresh readout it is
  +0.05.
- **A row does in training what the protein of every cell does, or more.** Five numbers against the protein of
  2,546 knockout cells: the rows are not shown to be worse, and read 0.14 more. This is the same thing the repair
  of the head showed with records.
- **A derived row can be trained in, and a wrong rule is trained in just as well.** `prop` moved IRF1 further than
  any other arm, toward a value that is wrong by a factor of 2.7. A propagation rule has to be tested on the record
  before it is given to a training. For IRF1 the rule "the mean of its type" fails that test.
- **Little else was paid.** The guard holds everywhere, and the control target reads better or the same. For
  `mixed` and `cells5` a loss of a step on the control target is not excluded, by 0.007.
- **For the plan** ([`GENERAL.md`](../../GENERAL.md), the step that retrains with the TEG): this is the first time
  "readable without its own record" was measured under a training, on one event. The answer is "a little, not
  resolved". Training one layer on 3,739 cells takes 40 minutes an arm on a laptop, and is stable once the gradient
  is taken exactly.

## The forecasts

Intervals: 8 of 8 held the value. Outcomes: for 16 of 21 the class given the highest probability happened
([`results/outcome_A.json`](results/outcome_A.json), [`results/outcome_B.json`](results/outcome_B.json)).

| missed | forecast | observed |
|---|---|---|
| H1 | refuted 0.50 | unresolved (given 0.38) |
| H2 | unresolved 0.50 | refuted (given 0.32) |
| the gain of the reference reader in `rows5` | unresolved 0.45 | no step (given 0.35) |
| `rows5` against `mixed` | brings back a step 0.50 | no change of a step (given 0.15) |
| `rows5` against `cells5` | no change of a step 0.56 | unresolved (given 0.30) |

Two things I had wrong. I expected the reference reader to move at least a little with IRF1's row in the training;
across the seam it does not move at all. And I expected the event of a row to matter; it does not.

## What was not as written

1. **The declaration was amended once**, before any state of a trained layer existed. The first training made the
   rows swing and was stopped by hand; how the gradient is taken was changed (exact for the cells of a step, no
   momentum), and the arm `cells5` of part B. The loss, the learning rate, the limit of steps, the rule, the
   hypotheses and the forecasts were not changed. [`DESIGN.md`](DESIGN.md), Section 7, reports every attempt in
   between with its training curve.
2. **The reading of part A was made twice.** `scripts/read.py` took the mean of the head's outputs where the
   recorded output is the log mean. The script's own check flagged it (`the_heads_own_output_is_reproduced: false`).
   So in the first reading (kept in the history: commit `d22df6c`) the guard compared a wrong set of 133 outputs and
   the line "the head as trained at this state" read another protein's output. Of its 747 values 34 changed, all in
   those two places. The round, the reference readers, the gains, the comparison between arms and the hypotheses
   are identical in both, and the guard holds in both (0.019 and 0.023 in the first). The script now stops when the
   head does not repeat its recorded output. Part B was read once, with the corrected script.
3. **The check with three readouts was not declared.** It is reported apart and judges nothing. Its first form,
   after part A, had only the fresh readout, and the write-up of part A drew a wrong conclusion from it (above).
4. **No arm read every row to within 0.03 in 600 steps.** Each is read as it is (DESIGN.md, Section 8).
5. **The readings were made in the cloud container from the states and grids of the owner's machine.** The copies
   were checked by hash, and the arrays built again in the container are identical byte for byte.
6. When the grids' status was checked on the owner's machine before the reading of part A, the last line of one
   grid's log was seen ("IRF1: shift: unreadable"). Nothing was decided by it.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration and its amendment: the loss, the arms, the rule, the guards, the forecasts |
| [`results/cache_manifest.json`](results/cache_manifest.json), [`results/pilot.json`](results/pilot.json) | the kept output of the layer before the last; every training curve that was looked at before the trial |
| [`results/trial_manifest.json`](results/trial_manifest.json) | the training of the five arms: the curves, how far the states moved; committed before each reading |
| [`results/trial.json`](results/trial.json), [`results/outcome_A.json`](results/outcome_A.json) | the reading of part A; the forecasts against it |
| [`results/trial_B.json`](results/trial_B.json), [`results/outcome_B.json`](results/outcome_B.json) | the reading of part B (it holds part A again, identical); the forecasts against it |
| [`results/seam_check.json`](results/seam_check.json) | the check added after the reading: three readouts |
| `scripts/cache.py`, `scripts/train.py` | run on the owner's machine: the kept output; the training of one arm and fold |
| `scripts/build_tt.py`, `scripts/read.py`, `scripts/score.py`, `scripts/seam_check.py` | the arrays and SPEC files; the reading; the score; the check |

To run again (the 70M checkpoint of the official release, the tokens and the work folder of
`encoder_stages_20261124`, and the weights of the head as trained are needed):

    python scripts/cache.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded THP1_ENC.npz --out CACHE
    for arm in rows4 rows5 mixed cells5 prop; do for f in 0 1; do
      python scripts/train.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm $arm --fold $f --out STATES
    done; done
    python scripts/build_tt.py --enc WORK_ENC --states STATES --work WORK --arms rows4,rows5,mixed,cells5,prop
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/read.py --work WORK --observe OBSERVE --head HEAD_WEIGHTS.npz --part B
    python scripts/score.py --part B
    python scripts/seam_check.py --work WORK --states STATES

## What this run cannot say

- One event in one kind of cell, on cells read before. A gain here is a lead for units that have not been read,
  not evidence.
- One layer of one checkpoint, trained on 3,739 cells for 600 steps. Nothing about a retraining of the encoder on
  its pre-training data, and nothing about training more layers.
- What the system's own head would read after such a training: it was not the readout of this trial.
- Whether the response can be put into the state for any readout: the loss here asked for one readout only.
- The propagation law of the paper's TEG. The arm `prop` is its simplest form, and for IRF1 that form is wrong.

## What is next

- The same trial with the head as trained as the fixed readout, read at the system's own output. A new training;
  the owner decides.
- If the response is wanted in the state for any readout: train under a readout that is fitted again while the
  layer trains, or under several, and read with a readout fitted on cells the layer has not used.
- A propagation rule is to be tested on the record first, as the proposed types were
  ([`type_growth_20261124`](../type_growth_20261124/)).
