# A typed state in place of the mean: does it keep what the encoder loses?

Written on 25 November 2026, **after the scripts were tried on made-up values and before any typed state of a real
cell was computed**. Nothing is trained in this run.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit is read.

## 1. The question

Two lines of work exist side by side.

- In this repository the TEG stands outside the model. It is the test, the source of the records for a repair of the
  head, and once a loss for one layer ([`last_layer_trial_20261124`](../last_layer_trial_20261124/)). For IRF1 in
  THP-1 the item is not closed: the state a head reads leaves 0.40 of the PD-L1 response where the input leaves
  0.18, and no repair after the encoder has brought that back.
- The owner's design of a model built on the TEG has the state itself typed: the input is written as typed events,
  an event changes only the blocks of the state it touches, and whatever reads, reads the state only. It has been
  run on small symbolic tasks, not on cells.

This run brings the second to the first at the one place where it needs no training. TEDDY gives a head **one mean
over all gene positions**. Here the output of the same encoder is **kept by type**: a state of a few blocks that
add up to that mean. The question is whether a reader of the typed state reads what a reader of the mean does not.

What is known already and bears on it ([`encoder_stages_20261124`](../encoder_stages_20261124/)): which genes are
tokens, and in what order, carries the response (0.21 and 0.24); the mean of what enters the first layer does not
(0.52); and eight means over windows of rank positions at the last layer read no more than the single mean (0.42).
So a state that is merely finer is not expected to keep it. What is new here is typing by gene: by the model's own
gene embedding, by the events of the TEG, and by the gene of the measured protein.

## 2. The cells

The THP-1 stratum of the screen, 7,478 cells with their tokens as in the earlier runs: 2,386 control cells, the four
pathway knockouts (IFNGR1, IFNGR2, JAK2, STAT1) and IRF1. The encoder is TEDDY-G 70M as released. The output of the
layer before the last is the one kept for the trial of 24 November; the last layer is applied to it as released.

## 3. The stages

| stage | what | width |
|---|---|---|
| `V` | the input view of the earlier runs | 512 |
| `Z70` | the mean of the last layer's output over a cell's gene positions: what a head reads today | 512 |
| `TK8` | 8 types of genes from the model's own gene embedding (k-means, the genes of the screen in the vocabulary, seed 20261125). Block *k* is the sum of the last layer's output over the positions whose gene is of type *k*, divided by the number of the cell's gene positions | 4,096 |
| **`TA`** | **2 types from the TEG: the 128 genes whose presence among a cell's tokens the events of the pathway change most, and the rest.** Block 1 is the sum over the positions of those genes, divided as above; block 2 is the mean minus block 1 | 1,024 |
| `TR` | the control of `TA`: the same with 128 genes drawn at random | 1,024 |
| `G` | the last layer's output at the position of the gene of the measured protein (CD274 for PD-L1), zeros when the gene is not among the cell's tokens, and one column that says whether it is | 513 |
| **`ZG`** | **`Z70` beside `G`** | 1,025 |
| `EA` | the typing of `TA` at the entry: the same two blocks of what enters the first layer (gene embedding plus position embedding) | 1,024 |
| `PERM` | the input view with its rows permuted | 512 |

- The blocks of `TK8`, `TA` and `TR` add up to `Z70`, and those of `EA` to the mean of what enters the first layer.
  A linear reader of a typed state can do what a linear reader of the mean does, and more only if a block holds
  something that the sum cancels.
- **The genes of `TA`.** Among the genes that are a token in at least 5% of the cells of the stratum: the share of
  the knockout cells of the four pathway genes that have the gene among their tokens, minus the same share of the
  control cells. The 128 genes with the largest absolute difference are the type. The choice reads tokens and the
  labels of those events and of the control cells. **It reads no protein value, and no cell of IRF1 enters it.** For
  the control target (the four pathway genes) the choice has seen the tokens of its own cells; this is said where
  that target is read.
- **The gene of `G`** is given by the record: the protein that is measured and the gene that encodes it. No cell is
  read to choose it.
- `TA` and `ZG` are the two states by which the hypotheses are read. `TK8`, `G` and `EA` are described.

## 4. The reading

As in the runs before: the tool's reference reader of every stage (a reader fitted with IRF1's own cells), what it
adds to the input view, and the round with the head reading the stage (the arms `controls`, with no record of the
event; `rows`, with its one record; `cells`). Every output of a typed stage is compared with the same output at
`Z70` on paired resamples. The halves, the cells and the resamples are those of the earlier rounds.

- The readers of the tool are the linear class only. The rule has chosen it at every stage of every run so far, and
  three of the stages here are wide.
- At `Z70` the round reads the state at unit length, as before. A typed state is given to the round divided by the
  length of that same mean state, so that its blocks add up to what the round reads at `Z70`; the column of `G` that
  says whether the gene is among the tokens is left as it is.

A step is 0.10. "A step" means the interval of a gain lies above 0.10; "no step", below it.

## 5. The hypotheses

For IRF1 on its own, and for each of `TA` and `ZG`:

| | confirmed | refuted |
|---|---|---|
| **the state keeps it** | the reference reader gains a step over `Z70` | no step |
| **the encoder's loss is closed** | the reference reader is no added loss from the input view | it is a loss |
| **the repair with the record** | the output with IRF1's one record gains a step over `Z70` and reaches a grade | no step |
| **without its own record** | the output with no record of IRF1 gains a step over `Z70` | no step |

Unresolved otherwise. And for `TA`: **the choice of genes matters** if `TR` costs a step against `TA` at the
reference reader.

The control target, the four pathway genes as a pool: with their records the output loses no step at a typed stage
and keeps its grade.

## 6. Checks

The script that makes the states stops, and writes none, unless:

- the mean it computes equals the recorded `Z70`, and the mean of what enters the first layer the recorded one, to
  2e-3;
- the blocks of `TK8`, `TA`, `TR` and `EA` add up to their mean, to 1e-4.

The reading stops, and writes nothing, unless the round at `Z70` repeats the first outer step and the reference
readers of `Z70` and of the input view repeat the last reading, to 0.005. The permuted view is expected to read
nothing.

## 7. What each outcome would mean

- **A typed state keeps it and closes the loss.** Then what TEDDY loses for IRF1 is lost by the single mean, the
  typed state is a repair that needs no training of the encoder, and the next step is a head that reads the typed
  state (a training of the head only).
- **A typed state gains a step and does not close.** Part of the loss is the mean. The rest is before the last
  layer's output, and `EA` says whether typing at the entry holds it.
- **No typed state gains a step.** Then the loss is not in how the last layer's output is summed. A typed state
  would have to be written at the input and the layers trained on it: a new encoder, which this machine cannot
  train on the pre-training data.
- **`TA` gains and `TR` gains as much.** Then it is the second block, not the type.

## 8. What this run cannot say

- One event in one cell line, on cells read before.
- Nothing about a state that is learned. The types here are fixed by rules; an event head that finds its own types
  is the part of the owner's design that does not exist yet.
- Nothing about an event without a record of any kind: `TA` is typed by the records of the pathway's events.

## 9. The forecasts

[`results/declared.json`](results/declared.json), 22 of them, written before any typed state of a real cell
existed. In short: I give "the state keeps it" a probability of 0.22 for `TA` and 0.15 for `ZG`; "the loss is
closed" 0.20 and 0.15; and a step without a record 0.25 for `TA`. For `TA` the most likely single outcome of every
hypothesis is "unresolved". For `ZG` I lean to "refuted" on closing the loss (0.45) and on reading without a record
(0.50). For `TK8` I expect no step (0.50), as with the eight windows.
