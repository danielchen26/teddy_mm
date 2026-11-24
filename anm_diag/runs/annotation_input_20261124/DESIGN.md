# TEDDY's annotation positions: three ways of running the encoder, side by side

Written on 24 November 2026, **after the code was checked on 16 cells and before the states of the stratum were
computed in any way of running with annotations, and before any reading of this run was made**.

Everything here is exploratory. The cells are development cells of the THP-1 screen that earlier runs have read. No
final unit is read and no model is trained.

## 1. Why this is asked

For IRF1 in THP-1 the encoder loses a step of the PD-L1 response between its input and the state Z that the head
reads: a reader of the input view leaves 0.18, a reader of Z leaves 0.40
([`encoder_stages_20261124`](../encoder_stages_20261124/)). All of that was read with the encoder run the way this
repository has always run it: the genes of a cell and nothing else.

TEDDY has a second, typed part of its input sequence that the repository never used: four annotation positions
(disease, tissue, cell type, sex) in front of the genes. The vocabulary of the checkpoint holds the tokens for them.
So one question was left open by every earlier run: **is part of what the encoder loses lost only because it was run
without its annotation positions?** If it is, the repair is to run it properly, with no training at all.

## 2. What the official code and the paper say

Read on 24 November 2026 from the official release (the model code, the tokeniser, the tutorial, the mapping files)
and the paper (Chevalier et al., arXiv 2503.03485).

| | what is there |
|---|---|
| the forward pass | the annotation tokens are embedded and placed **in front of** the genes; no position embedding is added to them; the genes keep positions 0, 1, … (`model.py`, `forward`) |
| the paper | four special tokens, one for each kind of annotation, are put in front of the genes and **the model is trained to say the label at each**. The labels are a supervision signal. The paper adds that giving metadata as context for the gene task "brought no measurable improvements" in early experiments (Section 3.2) |
| the vocabulary of the checkpoint | 11 disease labels, 17 tissue labels, 12 cell-type labels, `female`, `male`, `unknown`, and the four placeholders `<disease>`, `<tissue>`, `<cell>`, `<sex>` |
| the official tokeniser | writes the four labels of a cell through mapping files (for example monocyte → `precursor_cell`, acute myeloid leukemia → `cancer_disease`, blood → `hematopoietic_tissue`). A study without these columns gets `healthy`, `unknown`, `unknown`, `unknown` |
| the official tutorial | makes embeddings **without** any annotation position, and says the collate function "can be further customized to include ontologies" |
| what was not released | the code that put the annotation positions together in training: which token stood at a position when the label was to be said, how often, and whether a label token was ever given as an input. The paper says only that the model is "randomly prompt[ed] with classification tokens" with probabilities that balance the labels |

Two things follow before any number.

- **By the paper, the annotation is mostly an output, not an input.** What is put in is a placeholder; what the model
  gives back at that position is the label. So "running TEDDY with its annotations" has two readings, and this run
  makes both: the placeholders in front of the genes (what the paper describes), and the label tokens themselves in
  front of the genes (what "giving the model the annotation" would mean; whether training ever did that is not public).
- **Every cell of this stratum has the same four labels.** The annotation positions cannot tell one THP-1 cell from
  another. They can only change how the layers treat the genes. A large effect on one knockout would be surprising.

## 3. The ways of running, and the stages

All for the same 7,478 cells of the THP-1 stratum and the same tokens as `encoder_stages_20261124`.

| way | the four tokens in front of the genes |
|---|---|
| **N** | none: how the repository and the official tutorial run the model. This is the record |
| **P** | the placeholders `<disease>`, `<tissue>`, `<cell>`, `<sex>` |
| **T** | the labels of the THP-1 line by the official mapping files: `cancer_disease`, `hematopoietic_tissue`, `precursor_cell`, `male` |
| D (side) | what the official tokeniser writes for a study with no annotation columns: `healthy`, `unknown`, `unknown`, `unknown` |
| W (side) | labels that are wrong for these cells: `brain_disease`, `central_nervous_tissue`, `neural_cell`, `female` |

The stages that are read:

| stage | what it is |
|---|---|
| `V` | the input view of the earlier rounds, as it is |
| `ZN`, `ZP`, `ZT` | the state the head reads (the mean over the cell's gene tokens of the output of layer 12), in each way of running. `ZN` is the record |
| `ZD`, `ZW` | the same for the two side ways |
| `AP`, `AT` | the output of layer 12 **at the four annotation positions**, with the placeholders and with the labels: 2,048 numbers, reduced to 512 components by the recipe of `V` |
| `CP` | the first position with the placeholders, as it is: what the official class returns as `cell_emb` |
| `ZAP` | `ZP` beside `AP` (1,024 numbers): what a head would have if it read both |
| `PERM` | `V` with its rows permuted: a control that must read nothing |

`AP` matters for a reason of its own. The run before this one found that the response does not get into the *mean*
over tokens, and said that another summary "would have to be learned". The annotation positions are a summary the
model did learn in pre-training: each gathers from all gene tokens what it needs to say a label.

No readout and no label of an event enters any stage.

## 4. What is read, and the rule

As in `encoder_stages_20261124`: IRF1 on its own with its own cells (scope `own`), linear reader and boosting, the
reader of a stage chosen by the tool's own rule; and as a control the four pathway genes as a pool, linear reader.
Margin 0.10, step margin 0.10, 1,000 resamples; every difference is a difference on paired resamples.

- **Against the run without annotations.** For a stage s: |u| at s minus |u| at `ZN`. With the tool's step reading in
  both directions it is one of: **brings back a step** (the interval lies below −0.10), **costs a step** (above
  +0.10), **no change of a step** (inside −0.10 to +0.10), unresolved otherwise.
- **From the input.** |u| at s minus |u| at `V`: a loss, no added loss, or unresolved, as before.
- **The statement about annotations in the input**, from `ZP` and `ZT`:
  - *brings the response back*: one of them brings back a step and is no added loss from the input;
  - *brings part of it back*: one of them brings back a step and is still a loss from the input (if what is left
    from the input is unresolved, that is said);
  - *costs a step*: neither brings back a step and one of them costs one;
  - *changes nothing*: both are no change of a step;
  - unresolved otherwise.
- **Content.** `ZT` against `ZW` and against `ZD`: if the labels of the line do something that wrong labels also do,
  it is not their content.
- **The annotation positions.** `AP`, `CP`, `AT`, `ZAP` against `ZN` and from the input, by the same two readings.
- **Control.** For the four pathway genes `ZP` and `ZT` are not a loss from the input and do not cost a step.
- **Checks.** The state without annotations equals the recorded states to 1e-4 on all cells; *u* at `V` and at `ZN`
  repeat the recorded 0.184 and 0.396 to 0.005; `PERM` reads nothing. If a check fails: stop and report.

Described and not judged (written to `results/modes_manifest.json`): how far the state of a cell moves between the
ways of running; what the model's own output layer ranks first at the four positions when it is given the
placeholders; how large the input embeddings of the label tokens are beside those of ids that no token has (an input
that was never used in training keeps a small embedding).

What each outcome would mean:

| outcome | what follows |
|---|---|
| changes nothing, in both ways | the missing annotation positions are not where the response is lost. The placing of the run before stands. Nothing has to be recomputed |
| brings part or all of it back | the repository should run the encoder with the annotation positions; the states of every study would be recomputed, and the earlier readings of this stratum repeated |
| costs a step | giving these tokens harms what the head reads; the repository's way of running is the better one for this purpose |
| `AP` or `ZAP` brings back a step | the head should read the annotation positions beside Z: a repair without training the encoder |
| the labels of the line and wrong labels do the same | whatever changed is the presence of four tokens, not what they say |

## 5. The code check that was made before this was written

`scripts/annot.py --check 16` on the owner's machine, on the 8 longest and the 8 shortest cells; it prints agreement
numbers only (no state, no label, no prediction). [`results/check.json`](results/check.json).

| check | largest difference |
|---|---|
| the official class itself (loaded from the official code, on the CPU) against this code, without annotations, with placeholders, with labels | 0.0 at every position, three times |
| the state without annotations against the recorded state | 4.8e-6 |
| the order of the four annotation tokens (it must not matter: they carry no position) | 1.0e-5 in the mean over genes, 2.0e-5 at the annotation positions |
| the owner's machine (MPS) against the CPU | 1.1e-5 in the mean over genes, 3.7e-5 at the annotation positions |

One thing was changed in the check after its first run. It first compared every single gene position against 1e-4,
the limit used for means; between MPS and the CPU single positions differ by up to 3e-4 (values there reach 16). The
check now judges what is kept (the mean over genes to 1e-4, the four annotation positions to 1e-3) and reports the
largest single position beside them.

## 6. The order

1. Commit this file, the forecasts (`results/declared.json`), the code check and the scripts.
2. `scripts/annot.py` on the owner's machine: the states of the 7,478 cells in the five ways of running. Tokens only.
3. `scripts/build_modes.py`: the views and the SPEC files (the events, readers and decision of the run before, with
   the stages replaced). No readout enters.
4. `anm_diag scope` for the two SPEC files, on the owner's machine.
5. `scripts/read.py`: one reading; a result file that exists is not overwritten.
6. `scripts/score.py`: the forecasts against what was read.

The scripts are tried once on made-up values in a separate folder before step 2, as in the run before.

## 7. The forecasts

`results/declared.json`, 26 of them. Briefly, and why:

- **Nothing changes** is the most likely statement (0.47), then unresolved (0.35). The four tokens are the same for
  every cell.
- `ZP` and `ZT` read 0.30 to 0.50, as `ZN` does (0.40).
- **The annotation positions read less than the mean** (`AP` 0.40 to 0.95; a loss from the input with probability
  0.70): they are trained to say four coarse labels, and those are the same for all cells of this line.
- With the placeholders the model answers in kind at each position (0.75).
- Whether the label tokens were ever an input in training is open: their input embeddings look untrained with
  probability 0.55.
- The control holds (0.80).

## 8. What would change the plan

- A check of Section 4 fails: stop and report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.

## 9. What this run cannot say

- It is one event in one kind of cell, on cells read before.
- Between kinds of cell the labels differ. Whether the annotation positions help what is fitted on one kind of cell
  carry to another is a different question and is not asked here.
- It says nothing about the two larger checkpoints (160M, 400M) of the same release.
