# TEDDY's annotation positions: three ways of running the encoder, side by side

24 November 2026. TEDDY has four annotation positions (disease, tissue, cell type, sex) in front of the genes, and
this repository had never used them. This run asks whether part of what the encoder loses for IRF1 in THP-1 is lost
only because of that. It runs the encoder without them, with the four placeholders, and with the labels of the cell
line, on the same 7,478 cells, and reads all three side by side with the general code. The ways of running, the rule,
the control and 26 forecasts were committed before the states were computed ([`DESIGN.md`](DESIGN.md)).

It is also the first run whose audit was made on the owner's machine and not in the cloud container
([below](#the-audit-on-the-owners-machine)).

Everything here is exploratory. The cells are development cells that earlier runs have read, and IRF1 is one event.
No final unit was read and no model was trained.

## What came out

**Running TEDDY with its annotation positions changes nothing in the state the head reads. What the encoder loses for
IRF1 is not lost because the annotations were missing. The four annotation positions themselves read the response a
little better than the mean over the genes does, and by the rule that difference is not resolved.**

*u* is the share of the measured change that a reader of the stage leaves unexplained: 0 is none of it, 1 is all of
it. IRF1 is read on its own, by a reader that had its own cells. Differences are on paired resamples.

| stage | what it is | *u* | against the run without annotations | from the input | control: four pathway genes |
|---|---|---|---|---|---|
| `V` | the input view: RNA values | 0.18 [0.04, 0.30] | | | 0.08 |
| `ZN` | the state, **no annotation positions** (the record) | 0.40 [0.27, 0.50] | | +0.21 [0.12, 0.32]: loss | 0.08 |
| `ZP` | the state, **placeholders** in front of the genes | 0.41 [0.29, 0.52] | +0.02 [0.00, 0.04]: **no change of a step** | +0.23 [0.13, 0.34]: loss | 0.09 |
| `ZT` | the state, **the labels of the line** in front of the genes | 0.41 [0.29, 0.52] | +0.02 [−0.00, 0.03]: **no change of a step** | +0.23 [0.13, 0.34]: loss | 0.08 |
| `ZD` | the state, the official defaults (`healthy`, `unknown`, …) | 0.41 [0.28, 0.51] | +0.01 [0.00, 0.02]: no change of a step | +0.22: loss | 0.09 |
| `ZW` | the state, wrong labels (`brain_disease`, …, `female`) | 0.41 [0.29, 0.52] | +0.01 [−0.00, 0.03]: no change of a step | +0.23: loss | 0.09 |

The annotation positions, read in place of the mean:

| stage | what it is | *u* | against the run without annotations | from the input | control |
|---|---|---|---|---|---|
| `AP` | the four positions, with the placeholders | 0.30 [0.16, 0.40] | −0.10 [−0.18, −0.02]: unresolved | +0.11 [0.02, 0.21]: unresolved | 0.10 |
| `CP` | the first of them alone (the official `cell_emb`) | 0.29 [0.15, 0.40] | −0.11 [−0.19, −0.04]: unresolved | +0.10 [0.01, 0.20]: unresolved | 0.04 |
| `AT` | the four positions, with the labels of the line | 0.30 [0.16, 0.41] | −0.10 [−0.18, −0.02]: unresolved | +0.11 [0.02, 0.21]: unresolved | 0.11 |
| `ZAP` | `ZP` beside `AP` | 0.31 [0.18, 0.42] | −0.09 [−0.15, −0.02]: unresolved | +0.13 [0.03, 0.22]: unresolved | 0.07 |
| `PERM` | the input view with rows permuted | 1.00 | | | (reads nothing, as it must) |

- **The statement of the rule: annotations in the input change nothing.** With the placeholders the state reads 0.41,
  with the labels of the line 0.41, without either 0.40. Both intervals end below 0.04, less than half a step. The loss
  from the input is still there in every way of running (0.21 to 0.23).
- **What the labels say does not matter.** The labels of the line against wrong labels: +0.002 [−0.003, 0.007].
  Against the official defaults: +0.006. Against the placeholders: −0.002.
- **The annotation positions read more than the mean, and not by a resolved step.** The first position alone (512
  numbers, what the official class returns as `cell_emb`) leaves 0.29 where the mean over up to 2,048 gene tokens
  leaves 0.40. The difference is −0.11 [−0.19, −0.04]: its interval lies below zero and does not lie below −0.10, so
  by the rule it is **unresolved**, neither "brings back a step" nor "no change". From the input these positions are
  no longer a resolved loss (+0.10 [0.01, 0.20]) and not "no added loss" either. Reading them beside the state adds
  nothing to reading them alone.
- **The control holds.** For the four pathway genes no way of running loses anything or costs a step (0.04 to 0.11).
- **The checks hold.** *u* at `V` and at `ZN` repeat the record (differences below 1e-7), the state without
  annotations equals the recorded states to 2.8e-5 on all cells, and the permuted view reads nothing.

At every stage the linear reader read more than boosting, so every reading above is a linear reader's.

## What was learned about the annotation positions themselves

Described, not judged ([`results/modes_manifest.json`](results/modes_manifest.json)).

- **This code runs them as the official code does.** The official class, loaded from the official release, and this
  code agree exactly at every position, with and without annotations (16 cells, [`results/check.json`](results/check.json)).
  The order of the four tokens does not matter: no position is added to them.
- **The placeholders are how the model is asked, and it answers in kind.** Given the placeholders, the model's own
  output layer puts at least 99.9% of its probability on label tokens at each of the four positions, and ranks a label
  of that position's own kind first. For these cells it says: sex `male` for 95% of the cells (right for THP-1), cell type
  `immune_cell` for 91% (`precursor_cell`, where the official mapping puts monocytes, for 8%), tissue
  `respiratory_tissue` for 95% (the line comes from blood), disease `infectious_disease` for 68% and `cancer_disease`
  for 23% (the line is a leukaemia; the cells of this screen were stimulated, with interferon-γ among the stimuli).
- **The label tokens were an input in training.** Their input embeddings are as large as those of gene tokens (median
  norm 0.36 against 0.35), seven times those of ids that no token has (0.05). So the model has seen labels in front of
  the genes, not only placeholders. DESIGN.md read the paper as saying the labels are mostly an output; the checkpoint
  says they were an input as well. The code that put them together in training is not public.
- **The state moves a little.** With the placeholders the state of a cell moves by 0.26 of its distance from the mean
  of the stratum (median); with the shift common to all cells taken out, by 0.14. Between the labels of the line and
  wrong labels it moves by 0.04, counted the same way.

One thing was looked at after the reading and was not declared ([`results/said_by_unit.json`](results/said_by_unit.json);
counts of cells, no interval, no rule). The disease the model names follows the knockout of the pathway:

| cells | named `infectious_disease` | named `healthy` |
|---|---|---|
| control cells (2,386) | 78% | 0.3% |
| IRF1 knocked out (1,072) | 82% | 0.0% |
| IFNGR1 (1,339) | 63% | 4.5% |
| IFNGR2 (1,213) | 60% | 6.3% |
| JAK2 (1,024) | 59% | 6.0% |
| STAT1 (444) | 29% | 30% |

So the disease position responds to the interferon pathway being cut, most of all for STAT1. For IRF1 the name it
gives does not move, although a reader of the position's state reads IRF1's response somewhat better than a reader of
the mean. Why is not known from this run.

## What it means

- **The missing annotations were not the fault.** The question of the run has an answer for this event: giving TEDDY
  its typed positions, in either form, does not bring the response into the state. Nothing has to be recomputed, and
  the placing of the run before stands: the response does not get into the mean over the gene tokens.
- **The repository's way of running the encoder is confirmed for the state.** It equals the official forward pass,
  and the official tutorial runs the model the same way.
- **There is a second place to read, and it is not yet a result.** The annotation positions are a summary the model
  learned in pre-training, and for IRF1 they hold about half of what the mean loses (0.29 against 0.40; the input view
  0.18). That is one event on cells read before, by a margin the rule does not resolve, and the forecast went the
  other way (the positions were expected to read less). It is a lead to test, declared first: let the head read the
  first annotation position beside Z, on all three strata, as the run of 23 November did with the input view. It
  needs no training of the encoder.
- **For the trial that trains the last layer** (not run; the owner decides): the lead above comes first, because it is
  cheaper and because it changes what such a trial should aim at. If the annotation position carries what the mean
  does not, training has a target inside the model (make the mean, or the head, use it) in place of an outside one.

## The audit on the owner's machine

The cloud container that made every earlier grid has two processors. The owner's machine (Apple M4 Max) now has the
same environment and the audit code of this repository, in a folder of its own outside git. Before this run the two
grids of `encoder_stages_20261124` were made again there and compared with the cloud's
([`results/machine_check.json`](results/machine_check.json)):

| grid | every number of the result, largest difference | time in the cloud | time on the owner's machine |
|---|---|---|---|
| IRF1, 21 stages, linear reader and boosting | 1.7e-6 | 1,519 s | 589 s |
| four pathway genes, 21 stages, linear reader | 4.5e-6 | 675 s | 105 s |

One process per grid and one thread each, on a machine that was under other load. The grids of this run took 533 s
and 82 s there, and the states of the five ways of running 29 minutes.

## The forecasts

Intervals: 7 of 9 held the value. Outcomes: for 10 of 17 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)).

| missed | forecast | observed |
|---|---|---|
| *u* at `AP` | [0.40, 0.95] | 0.30 |
| *u* at `CP` | [0.45, 1.00] | 0.29 |
| `AP` from the input | loss 0.70 | unresolved (given 0.22) |
| `AP` against the run without annotations | costs a step 0.35 | unresolved (given 0.33), in the other direction |
| `ZAP` against the run without annotations | no change of a step 0.60 | unresolved (given 0.27) |
| the disease the model names | `healthy` 0.45, `cancer_disease` 0.40 | `infectious_disease` (another: 0.15) |
| the tissue it names | `hematopoietic_tissue` 0.50 | `respiratory_tissue` (another: 0.30) |
| the cell type it names | `precursor_cell` 0.32 | `immune_cell` (given 0.28) |
| the label tokens look untrained as an input | yes 0.55 | no (given 0.45) |

What was forecast right is what the run was for: nothing changes in the state (0.47, the most likely statement), and
the content of the labels does not matter. What was forecast wrong is everything about the annotation positions as a
place to read: they were expected to hold four coarse labels and little else, and they hold more.

## What was not as written

1. **The limit of the code check was restated before the declaration**, as DESIGN.md, Section 5 says: from every
   single position to what is kept.
2. **The official class needs a version of `transformers` below 5.** With version 5 it does not load. Version 4.57.6
   was used, and the official code was not changed.
3. **One table is reported that was not declared**: the disease the model names, by unit. It is marked above.

The scripts were tried on made-up values in a separate folder before the states were computed. Nothing had to be
changed.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: what the official code and the paper say, the ways of running, the stages, the rule, the control, the forecasts |
| [`results/check.json`](results/check.json) | the code check on 16 cells, made before the declaration |
| [`results/modes_manifest.json`](results/modes_manifest.json) | the ways of running and their tokens; how far the state moves; what the model says at the placeholders; the input embeddings |
| [`results/modes.json`](results/modes.json) | the reading: every stage and reader class, against the run without annotations, from the input, the content, the control |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| [`results/said_by_unit.json`](results/said_by_unit.json) | not declared: the disease the model names, by unit |
| [`results/machine_check.json`](results/machine_check.json) | the grids of the run before, on the owner's machine against the cloud |
| `scripts/annot.py` | run on the owner's machine: the states in the five ways of running; `--check` for the code check |
| `scripts/build_modes.py`, `read.py`, `score.py` | the views and SPEC files; the reading; the score |

To run again (the tokens and the work folder of `encoder_stages_20261124` are needed, and the official release for
`--official`):

    python scripts/annot.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded LAYERS_B.npz --out ANNOT --official OFFICIAL --check 16
    python scripts/annot.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded LAYERS_B.npz --out ANNOT
    python scripts/build_modes.py --enc WORK_ENC --annot ANNOT --ckpt TEDDY_G_70M --work WORK
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/read.py --observe OBSERVE
    python scripts/score.py

## What is next

- **The head reads the first annotation position beside Z**, for the three strata of the knockout screens, declared
  first: does the repaired output for IRF1 gain a step, and does any other target lose its grade? The states with
  placeholders of the two other strata take about ten minutes on the owner's machine.
- The trial that trains the last layer, after that and only if the owner decides to run it.
- Between kinds of cell the labels differ. Whether the annotation positions help what is fitted on one kind of cell
  carry to another was not asked here.
