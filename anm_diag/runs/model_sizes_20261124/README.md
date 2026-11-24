# The same reading on the two larger checkpoints of TEDDY-G

24 November 2026. Every reading so far was made on one checkpoint, TEDDY-G 70M, which loses a step of IRF1's PD-L1
response between its input and the state a head reads. This run asks whether that loss belongs to the small
checkpoint or to the way all of them are built: the same 7,478 cells of the THP-1 stratum, with the same tokens, were
run through the two larger checkpoints of the official release (160M and 400M) and read with the same tool. The
stages, the rule, the checks and 24 forecasts were committed before the states of the larger checkpoints were
computed ([`DESIGN.md`](DESIGN.md)).

Everything here is exploratory. The cells are development cells that earlier runs have read. No final unit was read
and no model was trained. The larger checkpoints have no trained head, so the run asks what a reader of a stage can
read, with the event's own cells: the question in which the loss was found.

## What came out

**Neither larger checkpoint is shown to lose less. Their states leave 0.35 where the small one leaves 0.40, a
difference of 0.04 that the rule does not resolve, and the two larger ones are the same as each other. None closes
the gap to the input (0.18). By the rule the statement is "unresolved" for both.**

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. IRF1 on its own,
read with its own cells; the reader of a stage is chosen by the tool's rule (the linear reader, at every stage).

| stage | 70M | 160M | 400M |
|---|---|---|---|
| what enters the first layer, averaged (`E`) | 0.52 [0.40, 0.62] | 0.45 [0.33, 0.56] | 0.49 [0.37, 0.60] |
| the layer before the last, averaged (`P`) | 0.56 [0.45, 0.66] | 0.38 [0.25, 0.49] | 0.37 [0.25, 0.47] |
| **the state a head would read (`Z`)** | **0.40 [0.27, 0.50]** | **0.35 [0.22, 0.46]** | **0.35 [0.23, 0.46]** |
| the same, reduced to 512 components | | 0.34 [0.21, 0.45] | 0.37 [0.25, 0.48] |
| the first annotation position (`C`) | 0.29 [0.15, 0.40] | 0.28 [0.14, 0.39] | 0.37 [0.26, 0.48] |
| the state with the placeholders in the input (`ZP`) | 0.41 | 0.36 | 0.35 |

The input view reads 0.18 [0.04, 0.30]; the view with its rows permuted reads nothing (1.00).

| | 160M | 400M |
|---|---|---|
| the state, from the input | +0.17 [0.08, 0.27]: unresolved | +0.17 [0.08, 0.28]: unresolved |
| the same, 512 components | +0.16 [0.07, 0.27]: unresolved | +0.19 [0.09, 0.30]: unresolved |
| the state against the small model's | −0.04 [−0.11, 0.04]: unresolved | −0.05 [−0.12, 0.04]: unresolved |
| the same, 512 components | −0.05 [−0.13, 0.03]: unresolved | −0.03 [−0.11, 0.07]: unresolved |
| what enters the first layer, from the input | +0.27 [0.17, 0.38]: a loss | +0.31 [0.20, 0.42]: a loss |
| the last layer against the one before it | −0.03 [−0.07, 0.02]: no change of a step | −0.01 [−0.06, 0.03]: no change of a step |
| the annotation position against the state | −0.08 [−0.14, −0.02]: unresolved | +0.02 [−0.05, 0.09]: no change of a step |
| the placeholders in the input against none | +0.01 [−0.00, 0.01]: no change of a step | −0.00 [−0.01, 0.01]: no change of a step |
| **the statement of the rule** | **unresolved** | **unresolved** |

For the small model the same lines read: the state from the input +0.21 [0.12, 0.32], a loss; the last layer
brings back a step (−0.16 [−0.23, −0.11]); the annotation position against the state −0.11 [−0.19, −0.04],
unresolved. The 400M against the 160M: −0.00 [−0.07, 0.06], no change of a step.

- **Size does not close it, and is not shown to help by a step.** Against the input the larger states leave 0.17
  more, where the small one leaves 0.21. The small model's interval lies above the step margin and theirs reach just
  below it (0.08), which is all that separates "a loss" from "unresolved" here. Against the small model the gain is
  0.04 to 0.05, with intervals from −0.12 to +0.04: a step is not shown and no change of a step is not shown either.
- **Width is not what makes the small gain.** Reduced to the 512 components of the small model, the larger states
  read the same (0.34 and 0.37 against 0.35 and 0.35).
- **The entry is a loss in every size.** What enters the first layer, averaged over the gene tokens, leaves 0.45 to
  0.52: the response is not in the mean of the gene embeddings in any of the three.
- **One thing does differ inside.** In the small model the mean stays empty through the layers and the last layer
  brings back a step (0.56 to 0.40). In the larger ones the layer before the last already reads what the last reads
  (0.38 and 0.37), and the last layer changes nothing. The end is the same.
- **The annotation position is not better in a larger model.** In the 160M it reads 0.28, as in the small one
  (0.29). In the 400M it reads 0.37, no more than its state. The lowest readings of any state of any size are the
  annotation positions of the 160M and of the 70M, and neither is inside the margin.
- **The placeholders in the input change nothing**, in either size, as in the small model.
- **The control holds.** For the four pathway genes, as a pool and with a linear reader, no state of any size is a
  loss from the input (0.03 to 0.15 at every stage, the input view 0.08).
- **The checks hold.** The stages of the 70M repeat their records to 5e-8; the permuted view reads nothing.

Boosting reads less than the linear reader at every stage of every size (for the states: 0.57, 0.48, 0.52) and is
chosen nowhere ([`results/sizes.json`](results/sizes.json)).

## What it means

- **Using a larger checkpoint is not the repair.** The outcome "a larger checkpoint does not lose it" did not
  happen, and "it loses less, by a step" did not either. The states of the three sizes read IRF1 at 0.40, 0.35 and
  0.35 against 0.18 at the input.
- **It points at the design more than at the size, and the rule does not let it be said as decided.** The loss at
  the entry is there in every size, and more than five times the weights move the state by 0.04. But
  the statement "loses it as the small model does" needs a resolved loss from the input and a resolved "no change"
  against the small model, and with one event and intervals of about ±0.10 neither is resolved.
- **For a trial that trains the last layer** (not run; the owner decides): nothing here is a reason to train a
  larger checkpoint. The small one is the cheapest, has the head and every recorded state, and its last layer is the
  one that moves the response into the mean. In the larger ones the last layer does not do that, so "train the last
  layer" would be a different trial there.
- **For the plan** ([`GENERAL.md`](../../GENERAL.md), step 2): this is the lightest form of a second encoder, and it
  is not an independent one. Three checkpoints of one family, trained on the same data by the same recipe, agree on
  where IRF1 stands. Whether the placing holds for a model of another family is still open, and the owner names it.

## What the checkpoints say about these cells

Described, not judged ([`results/sizes_manifest.json`](results/sizes_manifest.json)). The cells are THP-1, a human
monocytic leukaemia line. Given the four placeholders, what each checkpoint's output layer ranks first, as a share
of the 7,478 cells:

| | 70M | 160M | 400M |
|---|---|---|---|
| disease | `infectious_disease` 68%, `cancer_disease` 23% | `infectious_disease` 91% | `healthy` 88% |
| tissue | `respiratory_tissue` 95% | `hematopoietic_tissue` 82% | `respiratory_tissue` 58%, `sensory_tissue` 23% |
| cell type | `immune_cell` 91% | `immune_cell` 96% | `immune_cell` 99.7% |
| sex | `male` 95% | `male` 98% | `male` 94% |

## The forecasts

Intervals: 6 of 6 held the value. Outcomes: for 11 of 18 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)). The statement of the rule was forecast for both sizes
("unresolved", given 0.48 and 0.47).

| missed | forecast | observed |
|---|---|---|
| the state of the 160M, from the input | loss 0.50 | unresolved (given 0.38) |
| the state of the 400M, from the input | loss 0.42 | unresolved (given 0.40) |
| the state of the 160M against the small model's | no change of a step 0.45 | unresolved (given 0.37) |
| the last layer of the 160M | brings back a step 0.45 | no change of a step (given 0.18) |
| the last layer of the 400M | brings back a step 0.40 | no change of a step (given 0.23) |
| the annotation position of the 400M against its state | unresolved 0.52 | no change of a step (given 0.28) |
| what the 400M says at the disease placeholder | `infectious_disease` 0.50 | another (`healthy`; given 0.22) |

Two of the seven are one thing I had wrong: I expected the last layer to do in the larger models what it does in the
small one.

## What was not as written

Nothing in the stages, the rule or the reading. Two things about how it was made:

1. **The scripts were tried on made-up values in the cloud container, not on the owner's machine**, in a separate
   folder and before the grids were made, as DESIGN.md, Section 6 says they would be tried. Nothing had to be changed.
2. **The reading was made in the cloud container from the grids of the owner's machine.** `scripts/read.py` needs
   the result files of two earlier runs, which are in the repository. The states and both grids are the owner's
   machine's (the IRF1 grid 2,906 s under heavy load, the control's 635 s); the copies were checked by hash.

The manifest of the stages ([`results/sizes_manifest.json`](results/sizes_manifest.json)) was committed before the
grid of IRF1 was finished and before the reading; no readout enters it.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: the stages, the rule, the checks, the forecasts |
| [`results/check_70M.json`](results/check_70M.json), [`check_160M.json`](results/check_160M.json), [`check_400M.json`](results/check_400M.json) | the code checks: the official class against this code, the owner's machine against the CPU |
| [`results/sizes_manifest.json`](results/sizes_manifest.json) | the stages and their widths, the checkpoints, what each says at the placeholders |
| [`results/sizes.json`](results/sizes.json) | the reading: every stage, every difference, the statements, the control, the checks |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/sizes.py` | run on the owner's machine: the states of the cells for one checkpoint |
| `scripts/build_sizes.py`, `read.py`, `score.py` | the arrays and SPEC files; the reading; the score |

To run again (the checkpoints are those of the official release; the tokens and the work folder of
`encoder_stages_20261124` and the states of `annotation_input_20261124` are needed):

    python scripts/sizes.py --repo TEDDY_MM --ckpt TEDDY_G_160M --tokens TOKENS.npz --out S160
    python scripts/sizes.py --repo TEDDY_MM --ckpt TEDDY_G_400M --tokens TOKENS.npz --out S400
    python scripts/build_sizes.py --enc WORK_ENC --annot70 ANNOT_70M --s160 S160 --s400 S400 --ckpt TEDDY_G_160M --work WORK
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/read.py --observe OBSERVE
    python scripts/score.py

## What this run cannot say

- One event in one kind of cell, on cells read before. With one event the intervals are about ±0.10, which is the
  size of the step being asked about.
- Nothing about a head: none exists for the larger checkpoints.
- Nothing about a model of another family.

## What is next

- The trial that trains the last layer, on the small checkpoint, if the owner agrees to it.
- Before or beside it, what needs no training: how the schema is revised when an event differs from its type. The
  loop flags such events and has not yet proposed a type with several members and tested it on units that were not
  used ([`GENERAL.md`](../../GENERAL.md), Section 6).
