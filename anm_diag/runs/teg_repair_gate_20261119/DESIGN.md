# A repair computed from the rows of the TEG: a gate on cells already read. Design, written before the computation

Status: **development, not evidence.** Every knockout cell and every audit control used here was read for the nine arms of
[`runs/teg_training_20261119`](../teg_training_20261119/), whose design says that no further variant is fitted against
them. With the owner's consent (19 November) they are used once more, as development cells only: nothing read here counts
as a test of anything. A test needs cells that have not been read (the rest of the THP-1 screen, or the melanoma co-culture
cells), a design of its own and one reading.

## The question of the gate

`runs/teg_training_20261119` trained the head on the audit's statistic and did not meet its rule. The next idea is to let
the table compute the repair: a rejected row becomes a constraint (this event has to change the output by this much, with
its sign), and the repair is the smallest change that meets the constraints. Before that is designed and tested, one
assumption has to be checked, because the last run failed on it:

> With only the rows of the events of the fit (for each event, the measured change of PD-L1 with its sign, size and
> interval) and the control cells, does a repair carry to a gene it never saw, as far as a repair that sees the measured
> protein of every knockout cell?

If it does not, a repair from rows alone is not worth designing, and the weight moves to using the table to choose data.

## What is computed

The cells, the split of the controls and the held-out genes are those of `runs/teg_training_20261119` (same code and
seed). For each data set and each held-out gene, everything is fitted without that gene's cells and without the audit
controls, and read on that gene's knockout cells and the audit controls.

Two places for the repair:

- **A, at Z.** The PD-L1 output is replaced by a linear readout of TEDDY's embedding (unit length, as the head sees it).
  This is the most a linear repair can do with the information given, with nothing of the head kept.
- **B, in the head.** The head as trained is kept, and a linear correction from the embedding is added to its PD-L1
  output (in the head's own log scale). The other 133 outputs are untouched by construction. A term holds the corrected
  output on the 3,000 bone-marrow training cells where it was (weight 1, as the guard of the earlier runs).

Three kinds of information, in both places:

| name | what the fit is given |
|---|---|
| `controls` | the measured PD-L1 of the training control cells, cell by cell. Nothing about a knockout |
| `rows` | `controls`, and for each event of the fit one row: the measured mean change against the training controls, with its standard error |
| `cells` | `controls`, and the measured PD-L1 of every knockout cell of the fit. The reference: what the audit's readouts use |

`rows` is computed as an update of `controls`, in closed form. With `H` the ridge-regularised second moment of the
control cells, `w0` the `controls` solution, `D` the displacements of the events (mean of the knockout cells minus mean of
the controls, one row per event), `r` their measured changes and `S` the diagonal of their squared standard errors (on
the scale of the control term):

    w = w0 + H^-1 D' (D H^-1 D' + S)^-1 (r - D w0)

This is the smallest change of the readout, in the metric of the control cells, that brings each row inside its own
uncertainty; it is also the update a Bayesian linear model makes when each row is one observation. In B the same
formula is applied to the correction, with the targets on the head's scale (the measured value put on the scale of the
head's output among the training controls, as in the earlier repairs) and the bone-marrow term added to `H`.

Fixed here, not tuned on the audit cells: the ridge weight is chosen by five-fold cross-validation of `controls` on the
training control cells (grid 10^-3 to 10^2 times the mean variance of a coordinate) and used for all three kinds; for
`cells` it is chosen in the same way on its own training cells. An event is one gene (as in the run before); rows by
guide are computed as a variant (guides with at least 20 cells).

## The reading

As in `runs/teg_training_20261119`: *u* of the own PD-L1 output with `anm_diag.scope`, pooled (each knockout cell through
the fit that never saw its gene, audit controls split among the five fits; 1,000 resamples, knockout cells by gene) and
by event (one gene, all audit controls). A reading is valid only if the slope among the audit controls is positive and
the output shifts the way the measurement does. Without the statistic: correlation with the measured protein, and the
shift of the output over the measured shift. For B: the head's own score on the 5,000 bone-marrow test cells (only the
PD-L1 output can change).

Also, as description: the reading of `rows` with 1, 2, 3 and 4 events of the fit (every subset, point values), and the
angle between the displacement of the held-out event and the span of the events of the fit.

## The decision, written first

The comparison that decides is `rows` against `cells`, in place A, pooled, point values:

- **go**: in both data sets `rows` is valid and its |*u*| is at most 0.05 above that of `cells`. Rows carry about as far
  as cells do; the repair from rows is worth a design and a test on unread cells.
- **no-go**: in either data set `rows` is invalid or its |*u*| is 0.10 or more above that of `cells`. Rows alone are not
  enough; the weight moves to choosing data with the table.
- otherwise **undecided**, reported as such.

Reported beside it, not part of the decision: whether `cells` itself is inside the margin of 0.10 (whether the requirement
can be met by any linear repair on these cells), `rows` in place B against A, the guard, and `controls`.

## Forecast

| | melanoma | THP-1 |
|---|---|---|
| A `cells` | 0.00 to 0.10 | 0.15 to 0.25 pooled, because IRF1 stays near 0.6; the other four near 0.1 |
| A `controls` | 0.3 to 0.6 | 0.4 to 0.8 |
| A `rows` | 0.05 to 0.20 | within 0.10 of `cells` |
| B `rows` | A `rows` plus at most 0.05; guard kept | the same |

Forecast of the decision: go on melanoma, undecided on THP-1. The forecast of the run before, for the arm that mattered,
was wrong; this one is written with that in mind and may be wrong too.

## What this gate cannot show

- That a repair from rows is better than another method: nothing here is a test.
- Anything about genes of another kind than the ones in the table: five genes of one pathway per data set.
- Anything about a repair inside the head's MLP: B adds a linear path and leaves the MLP as it is.
