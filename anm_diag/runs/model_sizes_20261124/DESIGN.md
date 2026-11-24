# The same reading on the two larger checkpoints of TEDDY-G: declaration

Written on 24 November 2026, **after the code checks on a few cells and before the states of the stratum were
computed with the two larger checkpoints, and before any reading of this run was made**.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit is read and no model is trained.

## 1. Why this is asked

Every reading so far was made on one checkpoint, TEDDY-G 70M. For IRF1 in THP-1 that model loses a step of the PD-L1
response between its input and the state Z a head reads (0.18 at the input view, 0.40 at Z), and the loss was placed:
the response does not get into the mean over the gene tokens
([`encoder_stages_20261124`](../encoder_stages_20261124/)). Two things that need no training were then tried and
neither closed it: the annotation positions in the input changed nothing
([`annotation_input_20261124`](../annotation_input_20261124/)), and a head reading the first annotation position
beside Z regained a step ([`annotation_beside_z_20261124`](../annotation_beside_z_20261124/)).

What is left is training, and before a trial trains anything one cheaper question decides what it should train: **is
this loss a property of the small checkpoint, or of the way all of them are built?** The official release holds two
larger checkpoints of the same model, TEDDY-G 160M and 400M. They take the same tokens. So the same cells can be run
through them and read with the same tool.

- If a larger checkpoint does not lose the response, the cheapest repair is to use it, and training the small one is
  beside the point.
- If the larger ones lose it as the small one does, the loss belongs to the design (a mean over tokens as the state),
  not to the size, and a trial on the last layer has a reason and a choice of checkpoint.

It is also the lightest form of step 2 of the plan in [`GENERAL.md`](../../GENERAL.md), a second encoder: until now
nothing says whether the placing holds for more than one model. Two checkpoints of the same family, trained on the
same data by the same recipe, are not an independent second system.

## 2. The checkpoints

From the official release, by the file hashes it publishes.

| | layers | width | heads | weights (encoder and output layer) |
|---|---|---|---|---|
| 70M (used so far) | 12 | 512 | 8 | 71 million |
| 160M | 12 | 768 | 12 | 154 million |
| 400M | 24 | 1,024 | 16 | 394 million |

The checkpoint used so far is byte for byte the official 70M file. The three share one vocabulary (the same file), so
the tokens of a cell are the same for all of them, annotation tokens included.

**No head exists for the larger checkpoints.** The protein head of this repository was trained on the 70M's state.
So this run cannot ask what a head predicts. It asks what the audit asks of a stage: what a reader of that stage,
fitted with the event's own cells, can read. That is the question in which the loss was found.

## 3. The stages

All for the same 7,478 cells of the THP-1 stratum and the same tokens as the runs before. For a size m (70, 160, 400):

| stage | what it is |
|---|---|
| `E<m>` | what enters the first layer, averaged over the cell's gene tokens |
| `P<m>` | the output of the layer before the last, averaged |
| `Z<m>` | the output of the last layer, averaged: the state a head would read |
| `C<m>` | the output of the last layer at the first annotation position, with the four placeholders in front of the genes |
| `ZP<m>` | the state with the placeholders in front of the genes (the two larger sizes; the 70M was read on this day) |

`V` is the input view of the earlier rounds and `PERM` the same with its rows permuted. The stages of the 70M are the
recorded ones.

**Width.** A state of 768 or 1,024 numbers gives a linear reader more to work with than one of 512, whatever it
holds. So `Z160`, `Z400`, `C160` and `C400` are also read reduced to 512 components by the recipe of the input view
(`Z160r`, …), and **a gain or a closure counts only if the reduced state shows it too**.

## 4. What is read, and the rule

As in the runs before: IRF1 on its own with its own cells (scope `own`), linear reader and boosting, the reader of a
stage chosen by the tool's own rule; and as a control the four pathway genes as a pool, linear reader. Margin 0.10,
step margin 0.10, 1,000 resamples; every difference is a difference on paired resamples.

- **From the input.** |u| at a stage minus |u| at `V`: a loss, no added loss, or unresolved.
- **Against the small model.** |u| at `Z<m>` minus |u| at `Z70`: brings back a step (the interval lies below −0.10),
  costs a step (above +0.10), no change of a step (inside), unresolved otherwise.
- **The statement for each larger checkpoint**, in this order:
  - *does not lose it*: `Z<m>` and `Z<m>r` are both no added loss from the input;
  - *loses less, by a step*: `Z<m>` and `Z<m>r` both bring back a step against `Z70`;
  - *loses more, by a step*: `Z<m>` costs a step against `Z70`;
  - *loses it as the small model does*: `Z<m>` is a loss from the input and no change of a step against `Z70`;
  - unresolved otherwise.
- **Where, in each.** `E<m>` from the input (in the small model: a loss); `Z<m>` against `P<m>` (in the small model
  the last layer brings back a step).
- **The annotation position.** `C<m>` against `Z<m>` and from the input, and against the small model's `C70`.
- **The placeholders in the input.** `ZP<m>` against `Z<m>` (in the small model: no change).
- **Control.** For the four pathway genes no state of the three sizes is a loss from the input.
- **Checks.** The stages of the 70M repeat their records to 0.005 (`V` 0.184, `E70` 0.516, `P70` 0.556, `Z70` 0.396,
  `C70` 0.288); `PERM` reads nothing. If a check fails: stop and report.

Described and not judged (`results/sizes_manifest.json`): what each checkpoint's own output layer ranks first at the
four placeholders.

What each outcome would mean:

| outcome | what follows |
|---|---|
| a larger checkpoint does not lose it | the repair is to use that checkpoint: its states for every study, a head trained on them, and the earlier readings repeated. No trial on the small model |
| it loses less, by a step | size helps and does not close it; a trial trains the larger checkpoint's last layer, not the small one's |
| both lose it as the small model does | the loss belongs to the design. A trial on the last layer has its reason, and the small checkpoint is as good a place for it as any |
| unresolved | one event does not decide between size and design; said as that |

## 5. The code checks that were made before this was written

`scripts/sizes.py --check` on the owner's machine; it prints agreement numbers only.
[`results/check_70M.json`](results/check_70M.json), [`check_160M.json`](results/check_160M.json),
[`check_400M.json`](results/check_400M.json).

| checkpoint | cells | the official class against this code, at every position | the owner's machine against the CPU (mean over genes) | the 70M against its recorded state |
|---|---|---|---|---|
| 70M | 16 | 0.0 | 1.1e-5 | 4.8e-6 |
| 160M | 16 | 0.0 | 2.4e-6 | |
| 400M | 8 | 0.0 | 1.5e-6 | |

The script is the forward pass of `annotation_input_20261124/scripts/annot.py`, written for any number of layers
(the layer before the last is found by counting, not by its number).

## 6. The order

1. Commit this file, the forecasts (`results/declared.json`), the code checks and the scripts.
2. `scripts/sizes.py` on the owner's machine: the states of the 7,478 cells for the 160M and the 400M. Tokens only.
3. `scripts/build_sizes.py`: the views and the SPEC files. No readout enters.
4. `anm_diag scope` for the two SPEC files, on the owner's machine.
5. `scripts/read.py`: one reading; a result file that exists is not overwritten.
6. `scripts/score.py`: the forecasts against what was read.

The scripts are tried once on made-up values in a separate folder before step 4.

## 7. The forecasts

`results/declared.json`, 24 of them. Briefly, and why:

- **Unresolved is the most likely statement for both** (0.48 and 0.47): one event, and intervals of about ±0.12.
  After it, *loses it as the small model does* (0.32 for the 160M, 0.25 for the 400M). That a larger checkpoint does
  not lose it at all is given 0.08 and 0.12.
- The states are expected at 0.25 to 0.48 (160M) and 0.20 to 0.48 (400M), against 0.40 for the small model.
- **The entry is a loss in every size** (0.80): it is an average of gene embeddings in all of them.
- The placeholders in the input change nothing (0.85 each), as in the small model.

## 8. What would change the plan

- A check of Section 4 fails: stop and report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.

## 9. What this run cannot say

- One event in one kind of cell, on cells read before.
- Nothing about a head: none exists for the larger checkpoints.
- Nothing about a model of another family. That is step 2 of the plan proper, and the owner names the model.
