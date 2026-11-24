# Where inside the encoder: the stages between the input and Z

24 November 2026. For IRF1 in THP-1 the loop had placed a loss between the model's input and the state Z that its
head reads (0.21 [0.12, 0.32]), and only those two ends had been read. This run records what lies between them for the
7,478 cells of the stratum and reads the loss along the way with the general code. The stages, the rule, the control
and eighteen forecasts were committed before the stages were computed ([`DESIGN.md`](DESIGN.md)).

Everything here is exploratory. The cells are development cells that earlier runs have read, and IRF1 is one event.
No final unit was read and no model was trained.

## What came out

**The tokens the model is given still carry IRF1's response. The first average over them does not, and none of the
twelve layers puts it back into its average, except the last, in part.** The loss is located between the tokens and
the first mean over positions, not in the tokeniser and not in any single layer.

*u* is the share of the measured change that a reader of the stage leaves unexplained: 0 is none of it, 1 is all of
it. IRF1 is read on its own, by a reader that had its own cells. "From the input" is |u| at the stage minus |u| at
`V`, on paired resamples.

| stage | what it is | *u* | from the input | control: four pathway genes |
|---|---|---|---|---|
| `V` | the input view: RNA values | 0.18 [0.04, 0.30] | | 0.08 |
| `SELV` | the values of the genes that become tokens | 0.25 [0.12, 0.36] | +0.07 [0.01, 0.13]: unresolved | 0.10 |
| `RANK` | the tokens as the model receives them: order, no values | 0.21 [0.06, 0.32] | +0.02 [−0.05, 0.09]: **no added loss** | 0.09 |
| `E0` | what enters layer 1, averaged over the tokens | 0.52 [0.40, 0.62] | +0.33 [0.23, 0.45]: **loss** | 0.15 |
| `L1` to `L11` | the output of each layer, averaged | 0.51 to 0.60 | +0.33 to +0.42: a loss at every one | 0.12 to 0.19 |
| `L12` (Z) | the output of the last layer, averaged | 0.40 [0.27, 0.50] | +0.21 [0.12, 0.32]: loss | 0.08 |

Side views:

| stage | what it is | *u* | from the input |
|---|---|---|---|
| `SET` | which genes are tokens, without order | 0.24 [0.10, 0.36] | +0.06 [−0.02, 0.14]: unresolved |
| `W0` | what enters layer 1, as eight window means | 0.44 [0.31, 0.55] | +0.26 [0.15, 0.36]: loss |
| `W12` | the output of the last layer, as eight window means | 0.42 [0.30, 0.52] | +0.23 [0.14, 0.34]: loss |
| `PERM` | the input view with rows permuted | 1.00 | (reads nothing, as it must) |

- **The tokeniser is not where it is lost.** The tokens, with their order and without any value, read as much as the
  input view (0.21 against 0.18). Knowing only which genes are tokens reads 0.24; what the order adds is 0.04
  [−0.04, 0.11], unresolved. What the cut lets through with its values (`SELV`) is unresolved (+0.07): the tokeniser
  shows the model 60% of the values the input view holds, and that may cost a little.
- **The loss appears at the first average.** From the tokens to the mean of what enters layer 1 the reader loses 0.31
  [0.22, 0.42]. It is the only single step of the chain that is a loss, by the general code too
  (`lost_between`: `RANK` to `E0`). `E0` is a function of which genes are tokens and of how many there are, and of nothing else; that set
  reads 0.24: so the content is there, and its average through the gene embeddings does not keep it.
- **The layers do not put it into their averages, until the last.** `L1` to `L11` stay at 0.51 to 0.60, no step among
  them being a loss or a gain. From `L11` to `L12` the reader gains 0.16 [0.11, 0.23]. That comparison was not a rule
  of the design and is reported as an observation. The comparison that was (`E0` against `L12`) gives 0.12
  [0.03, 0.22]: unresolved. After the last layer 0.21 is still missing.
- **Eight windows in place of one mean change nothing at the end.** `L12` against `W12`: −0.02 [−0.07, 0.02], no added
  loss, and `W12` is itself a loss from the input. So the single mean is not what loses it relative to this richer
  summary, and a summary by windows of positions would not bring it back.
- **The control loses nothing.** For the four pathway genes, read as a pool with the linear reader, no stage is a loss
  from the input. The shape is the same on a smaller scale: the averages of the middle layers read less (0.12 to
  0.19, each unresolved) than the tokens (0.09) and than Z (0.08), and the last two layers bring it back.
- **The checks hold.** *u* at `V` and at `L12` repeat the recorded values (differences 0.0000), `L12` equals the
  stored states to 2.8e-5, and the permuted view reads nothing.

At every stage the linear reader read more than boosting, so every reading above is a linear reader's.

## What it means

- **The encoder's fault for this event is now placed more narrowly than "between the input and Z".** It is not a loss
  of input: what the model is given holds the response. It is that the response does not get into the average that
  Z is. The last layer is the one place that moves it there, and it moves part of it (0.16 of the 0.37 that the
  average of layer 11 misses).
- **What follows for a repair**, cheapest first:
  1. *Give the head what the tokens hold, beside Z.* This is what the run of 23 November did with the input view, and
     it brought the reader to 0.21, which is where the tokens read here. No further training is needed for
     that much, and it was not enough to close the item (one record).
  2. *Train the last layer.* It is where the average gets its content. A trial that trains one layer (2.1 million
     weights) is a third of the cost of the four layers proposed before this run, and it has a reason.
  3. *Another summary in place of the mean* would have to be learned (a weighting over tokens). Window means do not do
     it.
- **What is not shown.** A stage here is 512 numbers. The states of single tokens inside the layers were not read,
  so nothing says the middle layers lose anything: only that their averages do not show it. And this is one event
  in one kind of cell.

## The forecasts

Intervals: 8 of 9 held the value. Outcomes: for 5 of 9 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)).

| missed | forecast | observed |
|---|---|---|
| *u* at `SET` | [0.30, 0.85] | 0.24 |
| `SELV` from the input | no added loss 0.60 | unresolved (given 0.30) |
| the layers bring back a step (`E0` against `L12`) | yes 0.45 | unresolved (given 0.35) |
| the single mean at the end against eight windows | unresolved 0.40 | no added loss (given 0.30) |
| `W12` from the input | no added loss 0.35 | loss (given 0.30) |

The place was forecast (`E0`, 0.50), and the reason given for it was wrong. The design said the entry cannot hold the
response "because it has no order at all". The order adds almost nothing: the set of genes without order reads as
well as the tokens. The entry fails although it is given all that the set holds.

## What was not as written

1. **The builder was renamed before any stage was read.** A trial of the scripts on made-up layers (random numbers of
   the right shapes, in a separate folder) showed that `read.py` could not import it: its name clashed with the
   builder of an earlier run. It is `build_stages.py`; both scripts take `--results` so that a trial writes nothing
   here.
2. **One comparison is reported that was not a rule**: the step from `L11` to `L12`. It is marked as an observation
   above.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: the stages, the rule, the control, the checks, the forecasts |
| [`results/stages_manifest.json`](results/stages_manifest.json) | the stages: cells, tokens, the check against the stored states |
| [`results/stages.json`](results/stages.json) | the reading: every stage and reader class, from the input, step by step, the summaries, the control |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/layers.py` | run on the owner's machine: RNA counts, tokens, what the encoder holds at every stage |
| `scripts/build_stages.py`, `read.py`, `score.py` | the views and SPEC files; the reading; the score |

To run again (the exported views of the THP-1 cells and the work folder of `input_beside_z_20261123` are needed):

    python scripts/layers.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cells cells.tsv --out LAYERS
    python scripts/build_stages.py --work WORK --layers LAYERS
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/read.py --observe OBSERVE
    python scripts/score.py

## What is next

- The same reading for a second event that is lost before the head, when one exists. MHC-I in co-culture is the only
  other candidate, and its loss is not resolved.
- A trial that trains the last layer, declared first: what is held out, the control arm, the rule. It is read by the
  arms that ask whether an event is readable without its own record (`controls`, `unseen`).
- Reading the states of single tokens (the token of the measured gene, for one) is a further stage the code can take.
