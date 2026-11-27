# The learned pooling: what came out

27 November 2026. Run 50, on the owner's word. The first training experiment the loop asked for, as designed in
[`../learned_state_design_20261127`](../learned_state_design_20261127/DESIGN.md). The judge
([`scripts/read.py`](scripts/read.py)), the forecasts and the rule ([`results/declared.json`](results/declared.json))
and the pilot's training curve were committed before any trained state existed (`8f06da8`). The reading was made
once.

Exploratory: the cells are development cells of the THP-1 stratum that earlier runs have read. One event, one
readout. No final unit, no unread unit.

## What was trained

Four vectors of the model's width: 2,048 numbers. They set how the last layer's output is averaged over a cell's
gene positions, into four blocks in place of one mean. The encoder was not touched. The loss was one constraint:
at every step a readout is fitted on control cells alone, and it has to read the measured shift of every recorded
event. 600 steps, a few minutes each on the owner's machine; four arms, each on both halves of the cells.

## By the rule written first: it does not close

*carried*: the share of IRF1's change of PD-L1 that a state carries among cells whose state is alike. Every cell is
read at the state trained on the other half of the cells. A step is 0.10.

| state | carried |
|---|---|
| the input view | 0.89 |
| the mean, as released | 0.52 |
| the same four blocks, nothing trained | 0.50 |
| trained by the rows of the four pathway events, no cell of IRF1 (`pool4`) | 0.59 |
| **trained by the rows of all five (`pool5`)** | **0.66** |
| trained with the event labels shuffled (`mixed`) | 0.64 |
| typed by gene sets declared by hand, nothing trained (run 40) | 0.73 |

| | asked | measured | |
|---|---|---|---|
| H1, the trained state holds more | `pool5` over the mean: at least 0.10, interval above zero | +0.14 [−0.01, 0.31] | **unresolved** |
| H2, the gain belongs to the rows | `pool5` over `mixed`: the same | +0.02 [−0.02, 0.05] | **refuted** |
| H3, an event with no record gains | `pool4` over the mean: the same | +0.07 [−0.08, 0.22] | **unresolved** |

The guards hold: the four pathway events are read no worse (−0.04 to −0.01), a control readout is as good among
control cells (error ratio 1.00), and the arm with no row moved nothing. The author had put 0.35, 0.25 and 0.2 on
the three hypotheses.

## What it does show

- **Training the pooling moves the state the right way, on cells it never saw, and harms nothing that was
  measured.** Against the same four blocks untrained the gain is +0.16 [0.08, 0.25]. What a readout taught by
  control cells alone reads of IRF1's change goes from 0.10 at the mean to 0.32: that is what the loss asked for, and
  it holds on the other half of the cells.
- **It is not enough.** The trained state carries less than the state typed by hand (0.66 against 0.73; −0.07
  [−0.23, 0.08]) and much less than the input (0.89; the input holds 0.24 [0.03, 0.46] more).
- **The rows did not matter as rows.** With the event labels shuffled the result is the same. The last-layer trial
  of 24 November found the same thing. In this stratum that was nearly bound to happen, and the design should have
  said so: all five events lower PD-L1 (rows from −0.40 to −1.17), so shuffled rows (all near −0.85) ask for nearly
  the same thing as the true ones. What was learned is how a knocked-out cell differs from a control cell along the
  readout. Whether the rows of a TEG teach more than that can only be asked where rows differ in sign and in
  readout.

## What follows

The model is still the one part the loop has not corrected. Two ways on, both the owner's decision:

1. **The same experiment where the rows differ**: the melanoma strata, with 20 readouts, about a hundred recorded
   units in each whose rows go both ways, and four of the pairs that lose a part before the state. It needs a cache of the last layer's
   input for those cells (tens of GB on the owner's machine) and can end in a test on unread units (188 are left).
2. **Pooling and last layer together** under this loss, on the THP-1 cells again. Cheap, and limited in what it
   can show for the reason above.

Training the whole encoder with the constraint needs GPUs and the pre-training data.

## Files

| | |
|---|---|
| [`results/declared.json`](results/declared.json) | the forecasts and the rule, committed first |
| [`results/pilot.json`](results/pilot.json) | the pilot's training curve at the declared learning rate |
| [`results/reading.json`](results/reading.json) | the account of IRF1 at every state; the guards; the log of each training |
| [`results/outcome.json`](results/outcome.json) | the hypotheses |
| [`scripts/read.py`](scripts/read.py) | the judge |
| [`../learned_state_design_20261127/scripts/pool.py`](../learned_state_design_20261127/scripts/pool.py) | the training |
