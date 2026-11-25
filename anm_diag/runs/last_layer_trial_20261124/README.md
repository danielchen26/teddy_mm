# A trial that trains the last layer of TEDDY-G by the rows of the TEG

24 November 2026. For IRF1 in THP-1 the 70M checkpoint loses a step of the PD-L1 response between its input and the
state a head reads, and everything that needs no training was tried without closing the item. This run trains for
the first time: the last encoder layer of TEDDY-G 70M (2.1 million of 71 million weights), on the owner's machine,
with the rows of the TEG as the loss. The loss, the arms, the rule, the guards and 29 forecasts were committed before
any state of a trained layer existed ([`DESIGN.md`](DESIGN.md); it was amended once, still before any such state,
and says what was changed and why).

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit was read. One layer was trained, for 600 steps, on half of 7,478 cells at a time.

**Part A is read: the two arms by which the three hypotheses are judged. Part B (three more arms: `mixed`, `cells5`,
`prop`) is being trained and will be added here.**

## What came out

**Trained by the rows, the last layer makes the readout it was trained under read IRF1 on cells it never saw. With
IRF1's own row in the training that is a step: 0.90 to 0.62 with no record of IRF1 at the reading, 0.68 to 0.44 with
one. With no cell of IRF1 in the training it is 0.12 more, which the rule does not resolve. It does not bring back
what the encoder loses: the reference reader stays at 0.40, against 0.18 at the input. And most of the gain belongs
to that one readout: a readout fitted afresh on other control cells sees a third to a half of it.**

By the rule: **H1 unresolved, H2 refuted, H3 unresolved.** The control target is not harmed and the guard holds.

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. A cell takes its
state from the layer that was trained on the other half of the cells, so every state that is read is the state of a
cell the layer never saw, neither its tokens nor its protein.

**IRF1, on its own.** The input view reads 0.18 [0.04, 0.30].

| the output | as released (`Z70`) | `rows4`: trained by the rows of the four pathway genes, no cell of IRF1 used (`T4`) | `rows5`: the same and the row of IRF1 (`T5`) |
|---|---|---|---|
| a readout with no record of IRF1 (`controls`) | 0.90 [0.78, 1.07] | 0.78 [0.66, 0.91] | **0.62 [0.51, 0.73]** |
| the same with IRF1's one record (`rows`) | 0.68 [0.57, 0.79] | 0.58 [0.47, 0.69] | **0.44 [0.31, 0.55]** |
| fitted on every cell of the fitted half (`cells`) | 0.42 [0.30, 0.53] | 0.37 [0.24, 0.48] | 0.25 [0.09, 0.37] |
| the tool's reference reader, with IRF1's own cells | 0.40 [0.27, 0.50] | 0.40 [0.28, 0.51] | 0.40 [0.28, 0.51] |
| the head as trained, reading the state | 0.95 | 0.95 | 0.95 |

| against the released state, on paired resamples | `rows4` | `rows5` |
|---|---|---|
| no record of IRF1 | +0.12 [0.09, 0.17]: unresolved | **+0.28 [0.22, 0.37]: a step** |
| its one record | +0.10 [0.07, 0.14]: unresolved | **+0.24 [0.19, 0.31]: a step** |
| every cell | +0.05 [0.02, 0.09]: no step | +0.17 [0.12, 0.24]: a step |
| the reference reader | −0.01 [−0.02, 0.01]: no step | −0.01 [−0.03, 0.01]: no step |
| the reference reader, from the input view | +0.22 [0.12, 0.33]: a loss | +0.22 [0.12, 0.33]: a loss |

`rows4` against `rows5`: without IRF1's row in the training the readout with no record leaves 0.16 more
[0.12, 0.21], and the one with its record 0.14 more [0.11, 0.19]: the row of IRF1 in the training is worth a step.

- **H1, without its own record: unresolved.** In `rows4` no cell of IRF1 was used anywhere in the training, and the
  readout holds no record of it. It gains 0.12 [0.09, 0.17]: a step is not shown and "no step" is not shown either.
  The output is still rejected (0.78), so the item is not closed.
- **H2, the encoder's loss with the record: refuted.** The reference reader of `T5` is a loss from the input view,
  +0.22 [0.12, 0.33], as it is at release (+0.21). The reader did not move (−0.01 [−0.03, 0.01]).
- **H3, the repair with the record: unresolved.** The repaired output gains a step (0.68 to 0.44) and comes to
  within 0.04 [−0.03, 0.10] of the reference reader. It reaches no grade: gate 3 is not passed, so the grade is
  not tested.

**The control target** (the four pathway genes as a pool; their rows were in the training of both arms):

| the output | as released | `rows4` | `rows5` |
|---|---|---|---|
| no record of these events | 0.46 [0.39, 0.50] | 0.28 [0.22, 0.33] | 0.29 [0.22, 0.34] |
| with their records (`rows`) | 0.09 [−0.04, 0.16]: unresolved | 0.01 [−0.10, 0.09]: supported | 0.02 [−0.10, 0.09]: supported |
| the reference reader | 0.08 | 0.09 | 0.09 |

Without a record the pool gains a step in both arms (+0.18 [0.15, 0.20] and +0.17 [0.15, 0.19]). With the records
nothing is lost, the grade is kept, and the output moves from "unresolved" to "supported" (the gain itself, +0.08,
is not a step).

**The guard holds.** The head as trained, applied to the trained states: its 133 other outputs change on the control
cells by 0.019 (`rows4`) and 0.024 (`rows5`) of their spread per cell; the limit is 0.10. The head's own PD-L1
output does not change (IRF1: 0.95 before and after), which is expected: it reads none of the response either way.

**The checks hold.** With the head reading the released state the round repeats the first outer step for both
targets (largest difference 0.0000). The head as trained repeats its recorded output at the released state (7e-6).

**The training**, described and not judged ([`results/trial_manifest.json`](results/trial_manifest.json)). Every
arm and fold ran the 600 steps; none read every row to within 0.03, so each is read as it is. At the stop the rows
were read to within 0.09 and 0.05 (`rows4`) and 0.12 and 0.09 (`rows5`) on the fitted cells. The row of IRF1 was
learned in part: −0.30 for a row of −0.40 in one fold, −0.25 for −0.31 in the other. The states hardly move: in
units of the control cells' spread, 0.001 for control cells and 0.003 to 0.004 for the pathway knockouts. On the cells a
layer did not see, its own readout reads the pathway genes as it was taught (−0.86 to −1.12, from −0.48 to −0.62
at release), and IRF1 at −0.12 and −0.15 after `rows4` and −0.28 and −0.23 after `rows5` (−0.03 and −0.04 at
release). No protein value enters these.

## A check added after the reading

It was not declared and nothing is judged by it ([`scripts/seam_check.py`](scripts/seam_check.py),
[`results/seam_check.json`](results/seam_check.json)). It was added because two readers of the same state disagreed:
the reference reader did not move, and the round's arm fitted on every cell read a step more.

The trial has a seam. The two halves of the cells have their states from two different trained layers, and the round
and the tool both fit a reader on one half and read the other. For the round's readout this has a consequence that
the declaration did not spell out: **the readout that reads one half is fitted on the control cells of the other
half, which are the cells the layer was trained with. It is, nearly, the readout the layer was trained under** (the
same recipe on the same control cells, which the loss holds in place). So the tables above say: under the readout
it was trained with, on new cells.

The check reads inside one half, where every state comes from one layer and no cell was seen by it. The half is
split again, and the readers are fitted on cells of that half: a readout the layer has not met. IRF1, the gain over
the released state, first half and second half (half the cells, so the released values are higher: 0.99 and 0.94
for the readout with no record):

| the output | `rows4` | `rows5` |
|---|---|---|
| no record of IRF1 | +0.05 [0.03, 0.06], +0.05 [0.03, 0.08]: no step | +0.08 [0.07, 0.11]: unresolved, +0.14 [0.10, 0.23]: a step |
| its one record | +0.04, +0.04: no step | +0.07 [0.06, 0.10]: no step, +0.11 [0.08, 0.18]: unresolved |
| every cell | +0.03, +0.03: no step | +0.06: no step, +0.08 [0.06, 0.13]: unresolved |
| the tool's reader, with IRF1's own cells | +0.03, +0.02: no step | +0.06 [0.04, 0.08]: no step, +0.10 [0.06, 0.16]: unresolved |

For the pool of pathway genes the readout with no record gains 0.07 and 0.06 (no step), where the trial's reading
has 0.18.

- **A third to a half of a gain is left with a readout fitted afresh.** The rest belongs to the one readout the
  layer was trained under. The loss asks only that the cells of an event move along that readout, and that is what
  the layer learned.
- **Without a cell of IRF1 in the training, a fresh readout gains 0.05: no step, in both halves.**
- **The seam does hide something from the reference reader**, and it is small: inside one layer a reader fitted
  with IRF1's own cells gains 0.06 and 0.10 after `rows5`, where across the seam it gains nothing. What is lost
  between the input and the state is 0.21.

## What it means

- **The encoder's loss is not within reach of the last layer.** H2 is refuted, and the check does not change that:
  at most 0.10 of 0.21 comes back to a reader with IRF1's own cells, in one half. This agrees with where the loss
  was placed: the response is already missing from the mean of the gene tokens when it enters the first layer
  ([`encoder_stages_20261124`](../encoder_stages_20261124/), [`model_sizes_20261124`](../model_sizes_20261124/)).
- **A row used as a constraint in training does what a row does in the head's repair, and somewhat more.** With
  IRF1's row in the training, a readout that holds no record of IRF1 reads 0.62 on new cells; at release the
  repaired head needed the record to read 0.68. That is the record moved from the head into the encoder. It is not
  an event read without a record.
- **"No record, no reading" stands for IRF1.** With no cell of IRF1 in the training the gain is 0.12, unresolved,
  under the readout of the training, and 0.05 with a fresh one. The rows of the four other genes of its pathway
  reach IRF1 a little. What is missing is still a type, or a rule, that covers it; the arm `prop` of part B asks the
  simplest form of that.
- **Most of what was gained is the state bent toward one readout.** For a head that is kept, that is a real gain:
  the encoder was trained under it, and it reads more on new cells. For "the response is in the state", it is
  not. A trial that wants the second should train under more than one readout, or fit the readout again while it
  trains, and should read with a readout fitted on cells the layer has not used. The reading declared here did not
  do that, and the check does.
- **Nothing else was paid for it.** The control target reads better, the other 133 outputs of the head stay where
  they were, and the control cells' states move by a thousandth of their spread.
- **For the plan** ([`GENERAL.md`](../../GENERAL.md), the step that retrains with the TEG): this is the first time
  "readable without its own record" was measured under a training, on one event. The answer is "a little, not
  resolved". Training one layer on 3,739 cells takes 40 minutes an arm on a laptop, and is stable once the gradient
  is taken exactly.

## The forecasts (part A)

Intervals: 7 of 7 held the value. Outcomes: for 13 of 16 the class given the highest probability happened
([`results/outcome_A.json`](results/outcome_A.json)).

| missed | forecast | observed |
|---|---|---|
| H1 | refuted 0.50 | unresolved (given 0.38) |
| H2 | unresolved 0.50 | refuted (given 0.32) |
| the gain of the reference reader in `rows5` | unresolved 0.45 | no step (given 0.35) |

I expected the reference reader to move at least a little with IRF1's row in the training. Across the seam it does
not move at all.

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
   head does not repeat its recorded output.
3. **The check of the seam was not declared.** It is reported apart and judges nothing.
4. **No arm read every row to within 0.03 in 600 steps.** Each is read as it is (DESIGN.md, Section 8).
5. **The reading was made in the cloud container from the states and grids of the owner's machine.** The copies were
   checked by hash, and the arrays built again in the container are identical byte for byte.
6. When the grids' status was checked on the owner's machine, the last line of one grid's log was seen before the
   reading ("IRF1: shift: unreadable"). Nothing was decided by it.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration and its amendment: the loss, the arms, the rule, the guards, the forecasts |
| [`results/cache_manifest.json`](results/cache_manifest.json), [`results/pilot.json`](results/pilot.json) | the kept output of the layer before the last; every training curve that was looked at before the trial |
| [`results/trial_manifest.json`](results/trial_manifest.json) | the training of part A: the curves, how far the states moved; committed before the reading |
| [`results/trial.json`](results/trial.json), [`results/outcome_A.json`](results/outcome_A.json) | the reading of part A; the forecasts against it |
| [`results/seam_check.json`](results/seam_check.json) | the check added after the reading |
| `scripts/cache.py`, `scripts/train.py` | run on the owner's machine: the kept output; the training of one arm and fold |
| `scripts/build_tt.py`, `scripts/read.py`, `scripts/score.py`, `scripts/seam_check.py` | the arrays and SPEC files; the reading; the score; the check |

To run again (the 70M checkpoint of the official release, the tokens and the work folder of
`encoder_stages_20261124`, and the weights of the head as trained are needed):

    python scripts/cache.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded THP1_ENC.npz --out CACHE
    for arm in rows4 rows5; do for f in 0 1; do
      python scripts/train.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm $arm --fold $f --out STATES
    done; done
    python scripts/build_tt.py --enc WORK_ENC --states STATES --work WORK --arms rows4,rows5
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/read.py --work WORK --observe OBSERVE --head HEAD_WEIGHTS.npz --part A
    python scripts/score.py --part A
    python scripts/seam_check.py --work WORK

## What this run cannot say

- One event in one kind of cell, on cells read before. A gain here is a lead for units that have not been read,
  not evidence.
- One layer of one checkpoint, trained on 3,739 cells for 600 steps. Nothing about a retraining of the encoder on
  its pre-training data, and nothing about training more layers.
- Whether the response can be put into the state for any readout: the loss here asked for one readout only.
- The propagation law of the paper's TEG. The arm `prop` of part B is its simplest form only.

## What is next

- Part B: the arms `mixed` (does the event of a row matter), `cells5` (ordinary fine-tuning with the same cells) and
  `prop` (a row for IRF1 that was not measured).
- If this is followed: the same trial under a readout that is fitted again while the layer trains, read with a
  readout fitted on cells the layer has not used. It is a new training, and the owner decides.
