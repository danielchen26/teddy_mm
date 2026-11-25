# A typed state in place of the mean: what came out

25 November 2026. The question, the stages, the hypotheses and 22 forecasts were committed before any typed state of
a real cell existed ([`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json)). **Nothing is
trained in this run.** The encoder is TEDDY-G 70M as released.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit was read.

> **A correction, made later on 25 November** ([`typed_repair_20261125`](../typed_repair_20261125/)). This file said
> that reading a step more of IRF1 with no record of it (0.90 to 0.73) was new. That holds only for a readout
> fitted on control cells alone. The tool's reader that is fitted with the cells of the four pathway genes and never
> sees a cell of IRF1 leaves 0.64 at the mean state and 0.64 at the typed state (`reader_of_the_stage.units` in
> [`results/typed.json`](results/typed.json)). The typing did not improve the best reading without IRF1's record.
> The two passages that said otherwise are marked below. The run that follows also found that the round's repair
> did not apply its rows, which changes what "with its record" means here.

## The question

TEDDY gives a head one vector for a cell: the mean of the last layer's output over the cell's gene positions. For
IRF1 in THP-1 that mean leaves 0.40 of the PD-L1 response where the model's input leaves 0.18, and no repair after
the encoder has brought the difference back. The owner's design of a model built on the TEG keeps the state by type
instead. This run brings that idea to TEDDY at the one place where it needs no training: the output of the same
encoder is kept in blocks that add up to the mean, and the same readers and the same round read it.

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. A step is 0.10.

## What came out

**By the rule written first: the state typed by the genes that the pathway's events move keeps part of what the mean
loses (unresolved), and the round did not use it (refuted). A check added after the reading shows why the round did
not: the block is small, and the round's linear reader does not bring its columns to one scale. With every column at
one scale the same typed state gains a step in every output of the round, with IRF1's record and without it. Nothing
closes.**

### The reading as declared (IRF1 on its own)

| stage | width | a reader fitted with IRF1's cells (the tool's reference reader) | added to the input view | the round: no record of IRF1 | its one record | its cells |
|---|---|---|---|---|---|---|
| the input view | 512 | 0.18 [0.04, 0.30] | | | | |
| `Z70`: the mean, what a head reads today | 512 | 0.40 [0.27, 0.50] | +0.21: a loss | 0.90 | 0.68 | 0.42 |
| **`TA`: the 128 genes the pathway's events move most, and the rest** | 1,024 | **0.27 [0.13, 0.38]** | +0.09 [−0.00, 0.18]: unresolved | 0.87 | 0.62 | 0.40 |
| `TR`: 128 genes at random, and the rest | 1,024 | 0.43 [0.31, 0.53] | +0.25: a loss | 0.89 | 0.64 | 0.42 |
| `TK8`: 8 types from the model's gene embedding | 4,096 | 0.38 [0.25, 0.49] | +0.20: a loss | 0.93 | 0.69 | 0.42 |
| `G`: the position of the protein's gene (CD274) | 513 | 0.90 [0.83, 0.97] | +0.72: a loss | 0.97 | 0.95 | 0.96 |
| `ZG`: the mean beside that position | 1,025 | 0.43 [0.32, 0.54] | +0.25: a loss | 0.89 | 0.70 | 0.52 |
| `EA`: the typing of `TA`, of what enters the first layer | 1,024 | 0.24 [0.08, 0.37] | +0.06 [−0.04, 0.16]: unresolved | 0.95 | 0.72 | 0.48 |

The hypotheses of the design, for the two states they were written for:

| | `TA` | `ZG` |
|---|---|---|
| the state keeps it (the reference reader gains a step over `Z70`) | **unresolved**: +0.12 [0.08, 0.18] | refuted: −0.04 |
| the encoder's loss is closed (no added loss from the input view) | **unresolved**: +0.09 [−0.00, 0.18] | refuted: a loss |
| the repair with the record (the output with IRF1's record gains a step) | **refuted**: +0.06 [0.04, 0.08] | refuted |
| without its own record (the output with no record gains a step) | **refuted**: +0.03 [0.02, 0.05] | refuted |

- **The choice of genes matters.** `TR` costs a step against `TA` at the reference reader (+0.16 [0.10, 0.22]). A
  second block as such gives nothing; the block of the events' genes does. The 128 genes are interferon-induced
  genes (SERPING1, GBP5, GBP1, HLA-DRA, ETV7, SOCS1, GBP4, IL18BP, …); 122 of them are less often among the tokens of
  a knockout cell of the pathway than of a control cell, and a cell has 61 of them among its tokens on average.
- **The typing at the entry reads as much.** `EA` is the same two blocks of what enters the first layer, before any
  layer has acted: 0.24 against 0.27 for `TA` (no difference resolved). So what the typed state keeps is what the
  input holds and the mean erases. The twelve layers do not add to it.
- **A finer state by the model's own gene embedding gives nothing** (`TK8`, +0.01, no step), as the eight windows of
  rank positions gave nothing in [`encoder_stages_20261124`](../encoder_stages_20261124/).
- **The position of the protein's own gene gives nothing.** CD274 is among the tokens of 17% of the cells (29% of
  control cells, 8% of the pathway's knockout cells). The output at its position reads 0.90, and beside the mean it
  makes the round's arm with IRF1's cells worse (0.52 against 0.42).
- **The control target** (the four pathway genes as a pool, with their records: 0.09 at `Z70`) loses no step at
  `TK8`, `TA`, `TR` and `ZG` and keeps its grade everywhere. The rule reads "harmed" for `G`, which is no state for a
  repair, and for `EA` (0.16; the edge of the interval of the loss is −0.12).
- The checks held: the round at `Z70` and the reference readers of `Z70` and of the input view repeat the earlier
  readings (largest difference 2e-7); the permuted view reads nothing (1.00).

### A check added after the reading

Not declared before the reading. Its scripts and seven forecasts were committed before it was run
([`results/after_declared.json`](results/after_declared.json)); its stages are in
[`scripts/after_build.py`](scripts/after_build.py). It has two questions, both about what the reading above means.

**1. The genes.** The 128 genes of `TA` hold IRF1, JAK2 and STAT1: three genes that an event of the stratum knocks
out. A reader fitted with IRF1's cells could read the knockout from the gene's own transcript. `TAx` is the rule of
`TA` with the five knocked-out genes not eligible (WDFY1, NLRC5 and PTPN6 take the three places).

- It reads the same: 0.268 against 0.272 at the reference reader, and no change of a step in any output of the
  round. The gain does not rest on the knocked-out genes.
- IRF1 is among the tokens of 94% of IRF1's knockout cells and of 97% of control cells, so its own transcript does
  not mark the knockout here.

**2. The scale.** The tool's reference reader brings every column to one scale before it fits. The round's linear
reader centres the columns and penalises all of them alike. In what the round read, a column of the block of the
events' genes has a spread among control cells of 0.0009 where a column of the rest has 0.017: nineteen times
smaller. So the round could hardly use the block. The stages "at one scale" are what the round read, every column
divided by its spread among the control cells; `Zs` is the control, the same scaling of the mean alone.

| stage, at one scale | reference reader | the round: no record of IRF1 | its one record | its cells |
|---|---|---|---|---|
| `Z70` as it was read | 0.40 | 0.90 [0.78, 1.07] | 0.68 [0.57, 0.79] | 0.42 [0.30, 0.53] |
| `Zs`: the mean | 0.40 | 0.90 | 0.67 | 0.41 |
| **`TAs`: the events' genes and the rest** | 0.27 | **0.73 [0.63, 0.85]** | **0.49 [0.38, 0.60]** | **0.24 [0.10, 0.36]** |
| `TRs`: genes at random and the rest | 0.44 | 0.93 | 0.68 | 0.43 |
| `TAxs`: without the knocked-out genes | 0.27 | 0.72 | 0.49 | 0.24 |
| `EAs`: the typing at the entry | 0.25 | 0.79 [0.68, 0.91] | 0.41 [0.28, 0.53] | 0.17 [−0.01, 0.32] |

- **At one scale the round uses the typed state.** Over `Z70` as it was read, `TAs` gains +0.18 [0.13, 0.24] with no
  record of IRF1, +0.19 [0.14, 0.25] with its record and +0.18 [0.12, 0.24] with its cells: a step each. Over `Zs`
  the gains are the same, so it is not the scaling: the mean at one scale reads what the mean read.
- **It is the type.** `TRs` gains nothing and costs a step against `TAs` in every output.
- **The typing at the entry does as much with a record**: 0.41 with IRF1's record (+0.26, a step) and 0.17 with its
  cells, which is the input view (0.18). With no record it gains +0.11 [0.02, 0.22], unresolved.
- **Nothing closes.** Every output with a record is still rejected; the arm with IRF1's cells is no longer rejected
  at `TAs` and `EAs`.
- The control target loses nothing at one scale (with its records 0.03 at `TAs` against 0.09; with no record 0.37
  against 0.46).

Beside the run of 23 November, in which the head read the mean beside the 512 components of the input view
([`input_beside_z_20261123`](../input_beside_z_20261123/); not paired with this one): with IRF1's record that stage
gave 0.39, about what the typed states give here (0.49 and 0.41). **With no record of IRF1 it gave 1.00 against 0.90;
the state typed by the events' genes gives 0.73.** The block is chosen by the records of the pathway's events, and a
readout that has seen only control cells then reads part of the response to an event it has no record of.
*(Corrected, see the head of this file: a reader fitted with the cells of the pathway's events reads 0.64 at the mean
state already, and 0.64 at the typed state.)*

## The forecasts

Declared, before any typed state: **6 of 6 intervals held the value; the class I gave the highest probability
happened in 8 of 16** ([`results/outcome.json`](results/outcome.json)). I was right that `TA` would be unresolved on
keeping and on closing, that `ZG` would not close and that `TK8` would gain no step. I was wrong, among others, on
the choice of genes (I gave 0.20 to "costs a step") and on the round at `TA`, where I expected "unresolved" and the
rule read "refuted".

The check added after the reading: 2 of 2 intervals; 2 of 5 outcomes. I expected the knocked-out genes not to matter
(right) and the scale not to make a step (wrong three times: I gave "a step" 0.20, 0.15 and 0.15).

## What was not as written

1. **The scale of the round was a choice of the design, and it decided two hypotheses.** The design gave a typed
   state to the round "divided by the length of the mean state, so that its blocks add up". That keeps a small block
   small. "The repair with the record" and "without its own record" are refuted by the declared reading and stay so
   as the declared reading. The steps at one scale were read after that result was seen, by a rule written after it.
   They are exploratory and confirm nothing.
2. **`TA` holds three knocked-out genes.** The rule was declared without excluding them. I saw IRF1 and JAK2 in the
   list of genes when the states were made, before the reading, and left the rule as declared. The check above was
   added for it.
3. **Where it was run.** The states, the grids of the tool and both readings were made on the owner's machine, from a
   bundle of this branch's code and not from a clone. The result files are copied here unchanged (`typed.json`
   d6eae6c7…, `after.json` 66dfc752…).
4. A gene is twice among a cell's tokens in 2 cases in the stratum; the design says a gene has one position.
5. The readers of the tool are the linear class only, as declared.
6. The scale of the check uses all control cells of the stratum, those of the half that is read included. No protein
   value and no knockout cell enters it. The choice of the 128 genes reads the tokens of all knockout cells of the
   four pathway genes, both halves; no cell of IRF1.

## What it means

- **For the question of the run.** Keeping TEDDY's output by the types of the TEG, with nothing trained, brings back
  part of what the mean loses for IRF1: about half of it for a reader with IRF1's cells (0.40 to 0.27, with the
  input at 0.18), and, once the columns are at one scale, a step in what the round predicts: 0.68 to 0.49 with one
  record, 0.90 to 0.73 with none. It does not close the item.
- **It is a way around the encoder, not a repair of it.** The same typing of what enters the first layer reads as
  much. The block carries what the input holds about the interferon response; the layers neither add to it nor are
  needed for it. This agrees with the localisation of 24 November: the response is in which genes are tokens, and
  the mean over positions erases it.
- **What is new is the reading without a record.** *(Corrected, see the head of this file: this holds for a readout
  of control cells alone, not for the best reading without IRF1's record.)* Giving the head the whole input view did not carry to an event
  without a record (1.00), and training the last layer without a cell of IRF1 gained 0.12, unresolved
  ([`last_layer_trial_20261124`](../last_layer_trial_20261124/)). A block chosen by the records of the events of one
  type reads a step more of an event that has no record of its own (0.90 to 0.73), with nothing trained. Most of that
  response is still unread, and the step rests on a check that was not declared.
- **The round's reader has a flaw that this run found.** It does not bring the columns of a stage to one scale. For
  the mean this did not matter (`Zs` reads what `Z70` read). For any state made of blocks of unequal size it does.
  A change of the round's rule is the owner's decision and is not made here.
- **What is still missing of the owner's design.** The types here are fixed by a rule and by the records of one
  pathway. A part that finds the types from the input (the event head of the design) does not exist, and nothing
  here says how an event that belongs to no recorded type would be typed.

## What is next

1. **A declared repeat, with the scale in the rule from the start, on events that did not enter this run's check.**
   Nothing is trained. For every pathway gene of a stratum in turn: type by the events of the other genes, and read
   that gene with no record of its own; in THP-1 and in the melanoma strata. If the step without a record repeats,
   the typed state is a repair in front of the head that the loop can use; if not, IRF1 was one case.
2. **Only then a training.** The head that reads the typed state; and, beyond that, a part that learns the typing
   instead of the rule. Both need the owner's word.
3. **Not possible on this machine:** a typed state written at the input with the layers trained on it, on the
   pre-training data.

## Files

- [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json): written first.
- [`results/typed_manifest.json`](results/typed_manifest.json): the states (genes of every type, checks), committed
  before the reading. [`results/typed.json`](results/typed.json): the reading. [`results/outcome.json`](results/outcome.json): the forecasts.
- [`results/after_declared.json`](results/after_declared.json), [`results/after_manifest.json`](results/after_manifest.json),
  [`results/after.json`](results/after.json): the check added after the reading.
- `scripts/`: `typed.py` (the states; needs the model), `build_ty.py`, `read.py`, `score.py`; `after_states.py`,
  `after_build.py`, `after_read.py`.
