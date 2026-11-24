# Where inside the encoder: the stages between the input and Z

Written on 24 November 2026, **before the stages were computed and before any reading of this run was made**.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit is read and no model is trained.

## 1. Why this is asked

For IRF1 in THP-1 the loop placed a loss before the head: a reader of the input view with IRF1's own cells leaves 0.18
[0.04, 0.30] of the PD-L1 response, a reader of Z leaves 0.40 [0.27, 0.50], and the difference is 0.21 [0.12, 0.32]
([`outer_step_20261122`](../outer_step_20261122/)). On 23 November this was confirmed by an intervention: with the
input view beside Z the reader is back at 0.21 ([`input_beside_z_20261123`](../input_beside_z_20261123/)).

So the encoder loses it. Where in the encoder is not known, because only its two ends have been read. Between them lie
three different things:

1. **the tokeniser**: only genes that the symbol map and the vocabulary know become tokens; of those the 2,048 with the
   largest normalised value are kept; and the model is given their order, not their values;
2. **twelve layers**;
3. **one mean over positions**, which turns the states of up to 2,048 tokens into the 512 numbers of Z.

What would be done about a loss differs for each: another way of giving the input, a retraining of some layers, or
another summary of the last layer with no training at all. The general code reads a loss step by step along recorded
stages (`cycle.place(between=...)`; checked on planted systems, where the step in which a loss was planted is found).
This run records the stages for one stratum and reads them.

## 2. The stages

All for the same 7,478 cells of the THP-1 stratum (IRF1, IFNGR1, IFNGR2, JAK2, STAT1 and the control cells).

The chain, in the order of the model's own computation:

| stage | what it is | made |
|---|---|---|
| `V` | the input view of the earlier rounds: 512 components of RNA (log1p of counts per 10,000; the 13,289 genes detected in at least 1% of the control cells) | taken as it is |
| `SELV` | the same values, kept only for the genes that are among the cell's tokens: what the tokeniser lets through, with the values still there | here |
| `RANK` | the tokens as the model receives them: a gene's value is 1 for the first token and falls by 1/2,048 per rank. The order is there, the values are not | here |
| `E0` | what enters layer 1 (gene embedding plus position embedding), averaged over the cell's tokens | owner's machine |
| `L1` … `L12` | the output of each layer, averaged over the cell's tokens. `L12` is Z | owner's machine |

Side views, read in the same run:

| stage | what it is |
|---|---|
| `SET` | which genes are tokens, without their order |
| `W0`, `W12` | what enters layer 1, and the output of layer 12, as eight means, one for each window of 256 positions (ranks 1 to 256, 257 to 512, …) |
| `Z` | the stored state, as the earlier rounds read it (must read as `L12`) |
| `PERM` | `V` with its rows permuted: a control that must read nothing |

`SELV`, `RANK`, `SET`, `W0` and `W12` are wide or sparse. Each is reduced to 512 components by the recipe of `V`
(TruncatedSVD, the same seed), fitted on all cells of the stratum. No readout and no label of an event enters any stage.

Two things are known before the run and limit what it can say.

- **A stage here is a summary of 512 numbers, not everything the model holds at that point.** A later summary can read
  more than an earlier one. In particular `E0` carries nothing of the order: the mean of the position embeddings is the
  same for every cell of the same length, and the mean of the gene embeddings depends only on which genes are tokens.
  Order can reach a mean only through the layers. So "a loss from the input at stage s" means: this summary, with
  these readers, does not read it. `W0` and `W12` are there to say whether the single mean is what loses it.
- **512 components are a choice.** `V` was made that way and is the anchor, so the other input views are made the
  same way. Components are ordered by variance, not by what a response needs.

## 3. What is read

| | events | read as | readers | why |
|---|---|---|---|---|
| T | IRF1 | on its own, with its own cells (the scope `own`), as the record | ridge and boosting, as package 2 | the one loss placed between the input and Z |
| C | IFNGR1, IFNGR2, JAK2, STAT1 | a pool (the scope `units`), as the rounds read it | ridge only | the control: for these four nothing was lost before the head. Does a linear reader lose anything along the chain? |

Boosting for the pool at twenty stages is out of reach (it took 22 minutes for two stages on 23 November), so the
control is read with the linear reader only and says only that much.

Margin 0.10, step margin 0.10, 1,000 resamples; every difference is a difference on paired resamples.

## 4. What counts as what

The reader of a stage is chosen by the tool's own rule (`cycle._reader`, "union").

- **From the input.** For every stage s after `V`: |u| at s minus |u| at `V`. It is a **loss** if its interval lies
  above 0.10, **no added loss** if it lies below 0.10, unresolved otherwise.
- **Where the loss appears.** The first stage of the chain that is a loss from the input, and the last stage before it
  with no added loss. The run then says one of:
  - *the tokeniser*, if that first stage is `SELV` (what is let through) or `RANK` (order in place of values);
  - *the entry*, if it is `E0`: one mean of embeddings does not hold it;
  - *the layers*, if it is some `Lk`, with `E0` not a loss.

  If no stage is a loss, nothing is located and that is reported (it would contradict the record at `L12`).
- **Step by step.** The same chain by the general code: a single step from one stage to the next that is a loss is
  named. A loss of 0.21 spread over several steps need not show in any single one.
- **The mean over positions.** At the end: |u| at `L12` minus |u| at `W12`. A loss means the single mean drops at
  least a step that eight window means keep; no added loss means it does not; and `W12` from the input says whether
  the last layer still holds what the input view holds. The same at the entry (`E0` against `W0`).
- **What the layers do.** |u| at `E0` minus |u| at `L12`: above 0.10, the layers bring back at least a step of what
  the entry misses.
- **Order.** |u| at `SET` minus |u| at `RANK`: what the order of the tokens adds over knowing which genes they are.
- **Control.** For C no stage is a loss from the input.
- **Checks.** `L12` equals the stored states to 1e-4 (cells), and *u* at `V`, at `L12` and at `Z` repeat the recorded
  0.184 and 0.396 to 0.005. `PERM` reads nothing (|u| at least 0.9). If a check fails: stop and report.

What each outcome would mean:

| outcome | what follows |
|---|---|
| the tokeniser, at `SELV` | genes the model is never shown carry it: a matter of the vocabulary and the cut, not of training |
| the tokeniser, at `RANK` | the values carry it and their order does not: a matter of how the input is given |
| the entry, and the layers bring part back | the model has to rebuild it from order; what `L12` still misses is what the layers did not rebuild |
| the layers | a retraining aimed at those layers is the candidate (step 5 of the plan) |
| the single mean at the end loses it, `W12` holds it | another summary of the last layer would do, with no training |

## 5. The order

1. Commit this file, the forecasts (`results/declared.json`) and the scripts.
2. `scripts/layers.py` on the owner's machine: tokens and stages of the 7,478 cells; the check against the stored
   states. RNA counts only.
3. `scripts/build_stages.py`: the views and the SPEC files. No readout enters.
4. `anm_diag scope` for the two SPEC files.
5. `scripts/read.py`: one reading; a result file that exists is not overwritten.
6. `scripts/score.py`: the forecasts against what was read.

## 6. The forecasts

`results/declared.json`. Briefly, and why:

- **The tokeniser lets it through** (`SELV` close to `V`): the genes an IRF1 knockout moves are induced genes, and the
  normalisation by gene medians ranks such genes high. Most likely no added loss.
- **Order in place of values costs something, perhaps not a step** (`RANK` 0.18 to 0.50).
- **The entry does not hold it** (`E0` 0.45 to 0.95, a loss from the input with probability 0.75), because it has no
  order at all.
- **The layers bring part back** (`L12` reads less than `E0`), and the first stage that is a loss is most likely `E0`.
- **Whether the single mean at the end is what loses the rest is open**: a loss against eight windows 0.30,
  unresolved 0.40, none 0.30.
- **The control** shows no loss along the chain with probability 0.55; the entry may lose something for it too.

## 7. What would change the plan

- A check of Section 4 fails: stop and report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.
