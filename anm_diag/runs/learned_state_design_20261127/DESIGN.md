# A state that is not one mean, trained by the rows of the TEG

27 November 2026. **A design. Not run: nothing here was trained and no cell was read for it.** The script passes a
self-test on made-up cells. It is run on the owner's word, on the owner's machine.

## 1. What the loop says about the model

The loop corrects two of the three places where a change can be lost, and it does so without training:

| where the account places a change | what corrects it | state |
|---|---|---|
| not in the input | an edge on the event graph | done by rule (runs 43, 45, 47) |
| held by the state, not read | a readout taught by records of the type | done by rule (runs 40, 45) |
| **lost before the state** | **the model itself** | **never corrected** |

Where the third one stands:

- **It is measured.** IRF1 on PD-L1 in THP-1: a reader of the input leaves 0.18 of the change unexplained, a reader of
  the state 0.40. Batch 3 found four more pairs that lose a step before the state; batch 4 one (TXNDC17 on CD49f in
  the co-culture, 40% of the change).
- **It is located** ([`encoder_stages_20261124`](../encoder_stages_20261124/)). The tokens the model receives still
  carry IRF1's response (0.21). The first mean over the gene positions does not (0.52), and no layer puts it back
  into its mean except the last, in part (0.40). The loss is the mean: one vector of 512 numbers for a cell.
- **Three repairs were tried and none closes it.** Training the head ([`teg_training_20261119`](../teg_training_20261119/)):
  no arm trains it well. Training the last layer by the rows ([`last_layer_trial_20261124`](../last_layer_trial_20261124/)):
  one readout gains a step, a reader fitted independently gains nothing, and the gain is the same with the event
  labels shuffled. A state kept in blocks by declared gene sets ([`typed_state_20261125`](../typed_state_20261125/)):
  it keeps a part.

Two things follow. The fault is in how the state is formed, so that is what has to be trained. And a loss that
uses a fixed readout can be met by moving cells along that readout, which is what the last-layer trial did.

## 2. The experiment

**What is trained: the pooling, and nothing else.** A cell gets K = 4 vectors in place of one. Block *k* is the
last layer's output averaged over the cell's gene positions with weights softmax(q_k · h_i / √d). The four q are the
only trained numbers (2,048). With q = 0 every block is the mean, so training starts exactly at the released state.
The encoder is frozen, all 71 million weights.

**The loss is one row of the TEG written as a constraint on the state:**

> a readout taught by control cells alone has to read the recorded events.

At every step a ridge readout is fitted on 256 control cells at the state as it stands, in closed form, so the
gradient passes through the fit. For every recorded event: (the shift this readout reads from other control cells to
the event's cells, minus the measured shift)², averaged over the events. Nothing else is in the loss.

Why this loss. It is the quantity the loop already reads as its first system (S0). The readout never holds a record
of an event, at any step, so the rows can be met in one way only: by a state in which a knocked-out cell lies where a
control cell with that much less protein lies. Moving cells along a fixed vector does not help, because no vector
is fixed.

## 3. Cells

The 7,478 development cells of the THP-1 stratum, all read by earlier runs. What enters the last layer is already
kept on the owner's machine (the cache of the last-layer trial, 15.7 GB). The halves are those of the earlier rounds.
The protein of the other half is set to NaN before anything is computed. No final unit, no unread unit.

## 4. Arms

| arm | the record | what it is for |
|---|---|---|
| `none` | no row | nothing may move |
| `pool4` | IFNGR1, IFNGR2, JAK2, STAT1; no cell of IRF1 is used | does an event with no record gain? |
| `pool5` | the same and IRF1 | the main arm |
| `mixed` | as `pool5`, the event labels permuted among knockout cells | is the gain specific to the rows? In the last-layer trial it was not |

Each arm on both halves (two folds); every state that is read belongs to a cell the pooling never saw.

## 5. The judge is the loop, not the loss

On the cells of the other half, `account.account` reads IRF1 on PD-L1 at four states: the input view V, the mean Z,
the state typed by declared gene sets T, and the trained state P. A step is 0.10, as everywhere.

- **H1. The trained state holds more.** In `pool5` the share of the change the state carries is a step higher at P
  than at Z (the interval of the gain above zero and the gain at least 0.10). *Refuted* if the interval lies below
  0.10.
- **H2. The gain belongs to the rows.** `pool5` over `mixed`: a step.
- **H3. An event with no record gains.** `pool4`, for IRF1: a step.

Guards, read the same way: the four pathway events are read no worse by a step at P; among control cells the
cross-validated error of a control readout does not rise by more than 5%; `none` moves nothing.

What the author expects: H1 0.35, H2 0.25, H3 0.2; the guards all hold 0.7. Two priors pull apart: a pooling can
lift only what the last layer's outputs still hold at single positions, and the mean of that layer already recovers
a part; but the typed state, with sets chosen by hand and nothing trained, did keep more.

## 6. Order

1. On the owner's word: a pilot for the learning rate (`--pilot` writes the training curve only; no trained state is
   kept). The declared rate is confirmed or changed, and the change is written here.
2. The script that reads (the account at V, Z, T, P) and the forecasts above, committed before any trained state
   exists.
3. Eight trainings (four arms, two folds). Then the reading, once.

## 7. Cost

One frozen layer forward for 64 cells and the pooling, per step; no gradient passes through the encoder. 600 steps
at most per training. On the owner's machine this is less than one arm of the last-layer trial took. The cache has
to stay where it is until then.

## 8. What comes after, either way

- **If H1 holds:** the same on the co-culture of the melanoma screen (IFNGR1 and JAK1 on CD9, TXNDC17 on CD49f; a
  cache of about 20 GB for its cells), then one batch never read, forecast first, as the test. A state trained
  this way would enter the loop as a new version of the state, under the same gate as any rule.
- **If it does not:** pooling and last layer together under this loss (the loss is new; the trial's was not). If that
  gains nothing either, the outputs of the released encoder do not hold the response in a form that can be lifted,
  and what is left is training the encoder itself with the constraint. That needs GPUs and the pre-training data,
  and cannot be done from here.

## 9. What it cannot show

One event of one readout in one stratum, on cells read before. A gain here is a reason to spend unread units on
it, and not a result about the model.

## Commands

```
python scripts/pool.py --selftest
python scripts/pool.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm pool5 --fold 0 --out OUT --pilot
python scripts/pool.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm pool5 --fold 0 --out OUT
```

The self-test: made-up cells in which three positions hold a quantity the protein follows. Trained by four rows,
a control readout reads an event that had no row, on cells of the other half, at −0.57 where the mean state gives
−0.26 and the measured shift is −0.87.
