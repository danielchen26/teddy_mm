# A repair computed from the rows of the table, read once on cells not read before. Design, written before the reading

> **Note added after the reading (19 November).** The text below is as committed before the reading. Its statements that
> the co-culture cells had not been read, and that CD274 had not been used anywhere, are true of the run folders of this
> repository and too narrow: earlier work of this project in another workspace read the same cells for other heads.
> None of them entered this repair or its development. See the correction in [`README.md`](README.md).

Status: **a test.** The cells read here (melanoma, co-culture condition) have not been read by any run of this repository.
The repair, the reading, the rule and the forecast below are committed before the reading, and the reading is made once
(`scripts/test.py` refuses to run when its result file exists).

## The question

[`runs/teg_repair_gate_20261119`](../teg_repair_gate_20261119/) computed a repair of the PD-L1 output from the rows of the
table, on cells that had been read before. As development it looked workable in one form: the head's PD-L1 output is
replaced, only inside the cell type of the rows, by a linear readout of TEDDY's embedding that the rows have updated. That
run could not be evidence. This run asks, on cells not read before:

> Does a repair computed from the rows of the table make the audit pass on new cells? And how far beyond its own rows
> does it reach: to a gene of the pathway whose row it has not seen, and to a culture condition it has not seen?

## The cells

| | cells | status |
|---|---|---|
| melanoma, IFN-γ condition | 1,171 knockout cells of IFNGR1, IFNGR2, JAK1, JAK2, STAT1 and 4,000 control cells | read by earlier runs. Here: the source of rows for the transfer settings, never the cells of a reading |
| **melanoma, co-culture condition** | **2,000 knockout cells of the same five genes (321, 475, 390, 418, 396; 1,843 of them with a single guide) and 4,000 control cells** | **not read by any run of this repository. The cells of the reading** |
| melanoma, IFN-γ condition, the rest of the library | 287 knockout cells of CD274; 44,771 knockout cells of 200 other genes | CD274: not used anywhere before (the probe). The 200 other genes: their measured changes and mean displacements were tried as extra rows in development (below); no output was read on them |
| bone marrow | 3,000 training cells (for the gate), 5,000 test cells (for the guard) | as in the earlier runs |

The co-culture cells are every co-culture cell of the five genes that the rule of the declarations admits (a knockout cell
has guides for one gene only, plus any non-targeting guides) and 4,000 of the 5,897 co-culture control cells, with
embeddings made the same way as before (`scripts/cells.py` checks that protein values, targets and embeddings of the
cells already used are identical). [`runs/loop_20261112/DESIGN.md`](../loop_20261112/DESIGN.md) mentions that heads
retrained with IFN-γ cells were once audited on co-culture cells in another workspace; those records are not in this
repository and nothing of them was used here.

THP-1 is not part of this test: the rest of that screen has no embeddings yet.

## The repair, frozen

For a set of control cells C and a set of events (an event is one gene with its knockout cells):

- **Gate.** g(x) between 0 and 1: a logistic regression on embeddings alone, the control cells C against the 3,000
  bone-marrow training cells. No protein enters it.
- **Output.** On the log scale of the head: (1 − g(x)) · head(x) + g(x) · c(x), with head(x) the PD-L1 output as trained
  and c(x) a linear readout of the embedding (unit length, centred on C). The other 133 outputs are not touched.
- **`controls`.** c fitted by ridge on the measured PD-L1 of the control cells C, put on the head's scale.
- **`rows`.** `controls`, then one update: w = w0 + H⁻¹D′(D H⁻¹D′ + S)⁻¹(r − D w0). Each row is one event: the mean
  displacement of its knockout cells from the controls (D), the measured change of PD-L1 (r) and the square of its
  standard error (S, on the scale of the control fit). No knockout cell enters one by one. **This is the repair under test.**
- **`cells`.** Ridge on the control cells and on every knockout cell of the events, cell by cell. The reference.
- **Ridge weight.** Five-fold cross-validation on the cells of the fit, over 10^(k/2), k = −6 … 4, times the mean variance
  of a coordinate (for `rows` the weight of `controls`).
- **Scale.** The own output is read on the log scale, the scale of the measured protein (centred log-ratio).

What development changed since the gate, all of it on cells already read (`scripts/develop_rows.py`,
`develop_weight.py`, and `test.py --development`):

- The gate's grid of ridge weights had steps of a factor of ten and chose 0.1; with half steps the same rule chooses 0.32,
  and `rows` on a held-out gene reads 0.10 on melanoma (0.13 before) and 0.17 on THP-1 (0.19).
- Adding the rows of 200 other genes of the library (most of them with no change of PD-L1) changed nothing: 0.097 against
  0.096. They are not part of the repair.
- Choosing the ridge weight from the rows themselves (leave one row out, read against held-out controls) gave the same
  readings as cross-validation on the controls (0.098 and 0.177). The simpler rule is kept. The same choice read against the
  training controls picks a far too small weight and is not used.
- The code of the test, run on the cells already read with two halves (fit on one, read the other):

| development, `rows` | same events, new cells | a new gene | other four splits |
|---|---|---|---|
| melanoma, IFN-γ | 0.04 [−0.04, 0.10] unresolved | 0.04 [−0.03, 0.11] unresolved | 0.036 to 0.058, and 0.042 to 0.065 |
| THP-1 | 0.13 [0.01, 0.26] unresolved | 0.15 [−0.02, 0.31] unresolved | 0.134 to 0.162, and 0.151 to 0.179 |

  (`cells`: 0.00 and 0.01 on melanoma, supported; 0.12 and 0.13 on THP-1. The head as trained: 0.91 and 0.97.)

## Four settings

Two questions, crossed. Every reading is on co-culture cells.

| | the rows come from the co-culture condition | the rows come from the IFN-γ condition |
|---|---|---|
| the gene's own row is in the table | **1. same events, new cells** | **3. a new condition** |
| the gene's own row is left out | **2. a new gene** | **4. a new gene and a new condition** |

- **Settings 1 and 2** split the co-culture cells in two halves (knockout cells by gene, control cells at random; seed 0).
  The gate, the readout and the rows are computed on one half and read on the other, then the halves swap, so that every
  cell is read by a repair that never saw it. In setting 1 the table holds the rows of all five genes. In setting 2, for
  each gene, the table holds the rows of the other four.
- **Settings 3 and 4** compute everything in the IFN-γ condition (its 2,000 training controls of the earlier split, the rows
  of its knockout cells) and read all co-culture cells. No co-culture cell enters the fit. In setting 3 the table holds
  all five rows; in setting 4, for each gene, the other four.

In settings 2 and 4 the control cells of the reading are divided among the five repairs for the pooled reading, as before.

Also read, with the repair of setting 3 and against the 2,000 IFN-γ control cells that no fit has seen:

- **The probe: CD274**, the gene of PD-L1 itself, knocked out in the IFN-γ condition (287 cells). It is outside the pathway
  of the rows, and RNA need not show that the protein is gone.
- **The 200 other genes** of the library, which hardly move PD-L1: the change the repaired output predicts for each
  (slope among the controls times the shift of the output) beside the measured change. A repair that lowered PD-L1 for
  any knocked-out cell would show here. Described; no rule.

## The reading and the rule

*u* of the own PD-L1 output with `anm_diag.scope`, as in the runs before: 1 − a · shift(output) / shift(measured), a the
slope among the control cells of the reading; pooled over the five genes with 1,000 resamples (knockout cells resampled
by gene, as every pooled reading since `runs/teg_training_20261119`; control cells by cell), and by event (cells resampled
one by one). The interval with knockout cells resampled one by one is printed beside the pooled one; the rule uses the
first. Margin 0.10. A reading is valid only if the slope among the controls is positive and the output shifts the way the
measurement does. Guard: with the repaired output, the head's score on the 5,000 bone-marrow test cells (correlation with
measured counts) must not fall by more than 0.02, for PD-L1 and for the mean over proteins.

For each of the four settings, on `rows`:

- **holds**: the pooled reading is supported (the whole interval inside ±0.10), valid, and the guard is kept.
- **does not hold**: the pooled reading is rejected, or invalid, or the guard is broken.
- **undecided** otherwise.

Reported beside it: `rows` against `controls` and against `cells` on the same resamples (do the rows add anything, and how
far are they from the cell-by-cell reference), the head as trained, the readings by event, the gate on every group of
cells, and settings 1 and 2 with four other splits into halves (seeds 1 to 4, pooled values; the rule reads seed 0).

## Forecast

| | forecast |
|---|---|
| the head as trained, co-culture | rejected, *u* above 0.8 |
| 1. same events, new cells | `rows` between 0.00 and 0.12. Undecided more likely than holds (about 40% holds): with five genes as resampling units the interval is about ±0.07 wide |
| 2. a new gene | within 0.02 of setting 1, same outcome |
| 3. a new condition | `rows` between 0.05 and 0.30. Undecided most likely (60%), does not hold 25%, holds 15% |
| 4. a new gene and a new condition | as setting 3 |
| `rows` against `controls` | lower in settings 1 and 2, interval of the difference below zero |
| `rows` against `cells` | at most 0.05 above |
| guard | PD-L1 on bone marrow changes by less than 0.005 in every setting; the gate is below 0.05 on bone-marrow cells |
| probe CD274 | rejected, *u* above 0.5 |
| the 200 other genes | median predicted change below 0.05 in size; no gene beyond 0.5 |

Three forecasts of this line of work were wrong before. This one may be too.

## What this test cannot show

- Anything about another cell type, another pathway or another protein: one cell type, five genes of one pathway, one output.
- That the repair is better than another method. `cells` is a reference inside the same linear readout, not a competitor;
  the fine-tuned heads of the earlier runs are not read here.
- That the repair lives in the head: it replaces one output of the head, inside one cell type, by a readout of the embedding.
- With five genes as resampling units a pooled interval is approximate.

## How it is run

    python scripts/cells.py --protein <protein .h5ad> --states <states .npz> --name coculture --out <cells> --check <dense.npz>
    python scripts/cells.py --protein <protein .h5ad> --states <states .npz> --name library   --out <cells> --check <dense.npz>
    ANM_TEST_CELLS=<cells> python scripts/test.py --development     # cells already read
    ANM_TEST_CELLS=<cells> python scripts/test.py                   # the reading, once

with `ANM_FEATURES`, `ANM_BM_FIT` and `ANM_BM_TEST` as for `runs/teg_training_20261119`.
