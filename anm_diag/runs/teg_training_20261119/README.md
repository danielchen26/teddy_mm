# Training the head by the TEG, against four other kinds of knowledge — 19 November 2026

Exploratory. The design, a forecast for every arm and the rule are in [`DESIGN.md`](DESIGN.md), committed before any head
of this run was fitted and before any audit cell was read for it (commit `5001ecb`). This file reports what came out.

## Outcome by the rule: no arm trains the head well, and the rows of the TEG do not either

The example did not show what it set out to show.

- **S1 does not hold.** `rows` does not train the head well. (Its second half holds: none of `class`, `edge1`, `graph`,
  `judge` does.)
- **S2 does not hold.** Against the repair already made, `rows` is lower on melanoma by 0.05 [0.01, 0.10]; on THP-1 the
  difference is 0.18 [−0.05, 0.38], and the reading of `rows` there is not a valid one (below).

*u* of the head's own PD-L1 output on knockouts of genes held out of the fit, margin 0.10. An arm trains the head well if,
in both data sets, the reading is supported, the mean fall on bone marrow is at most 0.02, the slope among the audit
controls is positive and the output shifts the way the measurement does.

**Melanoma, IFN-γ** (1,171 knockout cells of 5 genes, each through the head that never saw its gene; 2,000 audit controls)

| arm | stands for | *u*, own PD-L1 output | reading | slope among audit controls | output shifts the measured way | \|*u*\| minus that of `rows` | bone marrow, mean fall | trains the head well |
|---|---|---|---|---|---|---|---|---|
| head as trained |  | 0.92 [0.90, 0.93] | rejected | 0.85 | yes | 0.67 [0.60, 0.76] | 0.000 | no |
| `data` | the repair already made | 0.30 [0.23, 0.35] | rejected | 3.85 | yes | 0.05 [0.01, 0.10] | −0.002 | no |
| `class` | an ontology | 0.48 [0.42, 0.54] | rejected | 3.86 | yes | 0.24 [0.19, 0.29] | 0.003 | no |
| `edge05` | a knowledge graph, 0.5 SD | 0.42 [0.38, 0.46] | rejected | 3.73 | yes | 0.18 [0.13, 0.24] | −0.003 | no |
| `edge1` | a knowledge graph, 1 SD | 0.33 [0.27, 0.39] | rejected | 3.64 | yes | 0.09 [0.05, 0.15] | −0.000 | no |
| `edge2` | a knowledge graph, 2 SD | 0.32 [0.23, 0.40] | rejected | 2.14 | yes | 0.08 [0.02, 0.14] | 0.006 | no |
| `graph` | a KG embedding | 0.32 [0.26, 0.38] | rejected | 3.65 | yes | 0.08 [0.04, 0.14] | −0.001 | no |
| `judge` | a decision engine | 0.37 [0.31, 0.43] | rejected | 5.04 | yes | 0.13 [0.09, 0.17] | −0.000 | no |
| **`rows`** | the TEG | 0.24 [0.15, 0.32] | rejected | 5.12 | yes |  | −0.002 | no |
| `data_rows` | the repair + the TEG | 0.22 [0.14, 0.29] | rejected | 4.62 | yes | −0.02 [−0.04, −0.01] | −0.002 | no |

**THP-1** (5,092 knockout cells of 5 genes, each through the head that never saw its gene; 1,193 audit controls)

| arm | stands for | *u*, own PD-L1 output | reading | slope among audit controls | output shifts the measured way | \|*u*\| minus that of `rows` | bone marrow, mean fall | trains the head well |
|---|---|---|---|---|---|---|---|---|
| head as trained |  | 0.97 [0.96, 0.99] | rejected | −0.55 | **no** | 0.71 [0.44, 0.94] | 0.000 | no |
| `data` | the repair already made | 0.44 [0.37, 0.53] | rejected | 4.35 | yes | 0.18 [−0.05, 0.38] | 0.007 | no |
| `class` | an ontology | 0.97 [0.92, 1.08] | rejected | 4.93 | yes | 0.71 [0.44, 1.00] | 0.002 | no |
| `edge05` | a knowledge graph, 0.5 SD | 0.67 [0.60, 0.79] | rejected | 4.30 | yes | 0.40 [0.19, 0.66] | 0.001 | no |
| `edge1` | a knowledge graph, 1 SD | 0.57 [0.50, 0.67] | rejected | 3.36 | yes | 0.31 [0.10, 0.53] | 0.002 | no |
| `edge2` | a knowledge graph, 2 SD | 0.61 [0.53, 0.70] | rejected | 1.97 | yes | 0.35 [0.12, 0.57] | 0.004 | no |
| `graph` | a KG embedding | 0.59 [0.51, 0.69] | rejected | 3.26 | yes | 0.32 [0.11, 0.56] | 0.004 | no |
| `judge` | a decision engine | 0.32 [0.20, 0.46] | rejected | 6.23 | yes | 0.06 [−0.12, 0.24] | 0.008 | no |
| **`rows`** | the TEG | 0.27 [−0.12, 0.53] | unresolved | −7.93 | **no** |  | 0.001 | no |
| `data_rows` | the repair + the TEG | 0.64 [0.43, 0.96] | rejected | −11.92 | **no** | 0.38 [0.21, 0.60] | 0.016 | no |

**Melanoma.** Every arm is rejected. `rows` leaves 0.24 [0.15, 0.32]: less than every other single kind of knowledge (by
0.05 to 0.24, every interval above zero), and still more than twice the margin. With the repair added (`data_rows`) it
is 0.22.

**THP-1.** `rows` reads 0.27 [−0.12, 0.53], unresolved, and the reading is an artefact: among the audit controls the
slope is −7.93, and on the knockout cells the output goes up where the protein goes down. The statistic does not change
when the slope and the shift change sign together, and the fit found that solution. `data_rows` did the same. The third
condition of the rule is there for this case. Of the arms with a positive slope, `judge` is lowest at 0.32 [0.20, 0.46]
and the repair already made is at 0.44 [0.37, 0.53]; both are rejected.

The guard on bone marrow is kept by every arm (largest mean fall 0.016, `data_rows` on THP-1; the head as trained scores
0.602).

## The forecast, row by row

| arm | forecast in `DESIGN.md` | came out | did the forecast hold |
|---|---|---|---|
| head as trained | 0.92 and 0.97 | 0.92 and 0.97 | yes (on record) |
| `data` | about 0.3 and 0.4, rejected | 0.30 and 0.44, rejected | yes |
| `class` | far from 0 (above 0.5 in size), rejected | 0.48 and 0.97, rejected | rejected: yes. Melanoma is just under 0.5 |
| `edge05`, `edge1`, `edge2` | not supported; `edge05` above 0, `edge2` below 0 | 0.42, 0.33, 0.32 and 0.67, 0.57, 0.61: all rejected, all above 0 | not supported: yes. `edge2` below 0: no |
| `graph` | as `edge1`; no better on THP-1 | 0.32 against 0.33; 0.59 against 0.57 | yes |
| `judge` | between `data` and the head as trained: 0.3 to 0.7 | 0.37 and 0.32 | inside 0.3 to 0.7: yes. On THP-1 it is below `data`, not between |
| `rows` | supported on melanoma; on THP-1 within ±0.2, possibly unresolved, IRF1 the worst event | 0.24, rejected; 0.27, unresolved, with slope and shift of the wrong sign; IRF1 the worst event (0.81) | no. Only "IRF1 the worst event" held |
| `data_rows` | close to `rows` | 0.22 against 0.24; 0.64 against 0.27 | melanoma: yes. THP-1: no |
| guard | kept by every arm | kept by every arm | yes |
| S2 | holds | does not hold | no |

The forecast for `rows`, the arm the run was made for, was wrong.

## Outside the rule

**The ten held-out events** (*u* of the own output; one gene, the head that never saw it, all audit controls)

| arm | mel IFNGR1 | mel IFNGR2 | mel JAK1 | mel JAK2 | mel STAT1 | THP-1 IFNGR1 | THP-1 IFNGR2 | THP-1 JAK2 | THP-1 STAT1 | THP-1 IRF1 | inside ±0.10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| head as trained | 0.92 | 0.91 | 0.92 | 0.92 | 0.91 | 0.97 | 0.98 | 0.98 | 0.96 | 0.97 | 0 |
| `data` | 0.30 | 0.32 | 0.21 | 0.33 | 0.31 | 0.42 | 0.42 | 0.43 | 0.33 | 0.62 | 0 |
| `class` | 0.49 | 0.49 | 0.44 | 0.48 | 0.50 | 0.94 | 0.95 | 0.96 | 1.22 | 1.00 | 0 |
| `edge05` | 0.42 | 0.43 | 0.39 | 0.43 | 0.43 | 0.63 | 0.63 | 0.65 | 0.73 | 0.99 | 0 |
| `edge1` | 0.34 | 0.32 | 0.31 | 0.35 | 0.35 | 0.54 | 0.56 | 0.58 | 0.54 | 0.82 | 0 |
| `edge2` | 0.32 | 0.30 | 0.30 | 0.34 | 0.33 | 0.60 | 0.64 | 0.66 | 0.55 | 0.78 | 0 |
| `graph` | 0.33 | 0.31 | 0.29 | 0.34 | 0.34 | 0.55 | 0.58 | 0.60 | 0.59 | 0.82 | 0 |
| `judge` | 0.36 | 0.36 | 0.32 | 0.41 | 0.39 | 0.27 | 0.31 | 0.32 | 0.04 | 0.67 | 1 |
| **`rows`** | 0.27 | 0.25 | 0.14 | 0.28 | 0.26 | 0.26 | 0.33 | 0.34 | −0.43 | 0.81 | 0 |
| `data_rows` | 0.23 | 0.23 | 0.12 | 0.26 | 0.24 | 0.61 | 0.60 | 0.60 | 0.19 | 0.62 | 0 |

IRF1 in THP-1 is the worst event of every arm with a positive slope except `class` (0.62 to 0.99). One event of one arm
has its point value inside the margin: `judge`, STAT1 in THP-1, 0.04 [−0.09, 0.15], unresolved.

**The events of the fit against the held-out event** (mean over the five heads of a data set)

| arm | melanoma: *u* on the events of the fit | on the held-out event | THP-1: on the events of the fit | on the held-out event |
|---|---|---|---|---|
| `data` | 0.09 | 0.29 | 0.35 | 0.44 |
| `class` | 0.33 | 0.48 | 1.05 | 1.01 |
| `edge05` | 0.22 | 0.42 | 0.40 | 0.73 |
| `edge1` | 0.13 | 0.33 | 0.20 | 0.61 |
| `edge2` | 0.13 | 0.32 | 0.38 | 0.64 |
| `graph` | 0.12 | 0.32 | 0.24 | 0.63 |
| `judge` | 0.19 | 0.37 | 0.18 | 0.32 |
| **`rows`** | 0.00 | 0.24 | 0.02 | 0.26 |
| `data_rows` | 0.00 | 0.22 | 0.21 | 0.53 |

`rows` fits the events it is given: *u* is 0.00 on them in melanoma and 0.02 in THP-1 (there by the wrong sign). On the
gene it never saw, 0.24 is left. On melanoma the distance between the events of the fit and the held-out event is
between 0.15 and 0.24 in every arm.

**Without the statistic** (mean over the five held-out events of a data set)

| arm | melanoma: correlation with the measured protein | shift of the own output / measured shift | THP-1: correlation | shift / measured shift |
|---|---|---|---|---|
| head as trained | 0.24 | 0.43 | −0.13 | −0.22 |
| `data` | 0.65 | 1.33 | 0.72 | 1.22 |
| `class` | 0.50 | 1.34 | 0.37 | 0.11 |
| `edge05` | 0.62 | 1.11 | 0.56 | 0.67 |
| `edge1` | 0.56 | 1.51 | 0.56 | 1.00 |
| `edge2` | 0.50 | 1.92 | 0.48 | 1.55 |
| `graph` | 0.57 | 1.51 | 0.55 | 1.11 |
| `judge` | 0.59 | 1.30 | 0.67 | 1.67 |
| **`rows`** | 0.60 | 1.52 | −0.52 | −1.77 |
| `data_rows` | 0.61 | 1.56 | −0.42 | −1.28 |

By these numbers, which do not use the statistic, the repair already made is as good as `rows` or better on melanoma, as
`DESIGN.md` said it would be. On THP-1 the output of `rows` is anticorrelated with the measured protein (−0.52): by any
plain measure that head is worse than the repair already made.

**The fits** (five per cell; kept step of 400; bone-marrow outputs moved, root mean square)

| arm | melanoma: kept step | bone-marrow outputs moved | THP-1: kept step | bone-marrow outputs moved |
|---|---|---|---|---|
| `data` | 400 | 0.031 to 0.032 | 400 | 0.062 to 0.064 |
| `class` | 180 to 200 | 0.054 to 0.057 | 300 to 400 | 0.055 to 0.061 |
| `edge05` | 400 | 0.030 to 0.030 | 400 | 0.040 to 0.045 |
| `edge1` | 400 | 0.038 to 0.040 | 380 to 400 | 0.048 to 0.054 |
| `edge2` | 400 | 0.066 to 0.072 | 380 to 400 | 0.070 to 0.081 |
| `graph` | 400 | 0.041 to 0.042 | 400 | 0.050 to 0.059 |
| `judge` | 400 | 0.041 to 0.043 | 400 | 0.055 to 0.057 |
| **`rows`** | 380 to 400 | 0.041 to 0.042 | 80 to 180 | 0.092 to 0.131 |
| `data_rows` | 400 | 0.045 to 0.046 | 40 to 60 | 0.168 to 0.196 |

On THP-1 the fits of `rows` and `data_rows` stopped early and moved the outputs on bone-marrow cells further than any
other arm (0.09 to 0.13 and 0.17 to 0.20, against 0.04 to 0.08), though the score on the bone-marrow test cells stayed
within the guard.

The refit of arm `data` reproduces the heads of [`runs/loop_guard_20261112`](../loop_guard_20261112/): the largest
difference of the PD-L1 output over all cells is 0.002 (melanoma) and 0.003 (THP-1).

The step `DESIGN.md` names for the case that S1 holds (the pooled audits of `runs/loop_guard_20261112` on the heads of
`rows`) was not run, since S1 does not hold.

## What this shows, and what it does not

It shows, for this loss, this budget and four events per fit:

- Writing the audit's statistic as the loss makes the head fit its events and does not carry to an event it never saw.
- The statistic alone can be met the wrong way round. A loss built on it needs the direction fixed by something else.
- No kind of knowledge about four genes, in any of these forms, brings the head's own output within the margin on a fifth.
- The geometry of the graph adds nothing to the signed edge here; the class of a cell alone leaves THP-1 where it was
  (0.97 against 0.97).

A defect of `rows` as written, said after the fact: a row of the TEG carries the sign and the size of the measured change,
and the loss used neither directly, only the ratio *u*. `DESIGN.md` named the risk of a wrong-sign solution and relied on
the cell-by-cell fit of the control cells to prevent it. On THP-1, where the head as trained starts with a negative slope,
that was not enough.

It does not show:

- that training by the TEG cannot work. One loss was tried, with one number per event and without its sign;
- anything about pretraining a foundation model by the TEG. Only the three layers of the head's MLP were trained; TEDDY
  is frozen;
- anything about an ontology, a knowledge graph, a KG embedding or a decision engine as products. The four arms are
  stand-ins: each carries the least that its kind of system carries;
- anything confirmatory. The target was chosen knowing the earlier repair leaves it rejected, these audit cells had been
  read before for other heads, there is one split and one fit per gene, and five genes per data set are few resampling
  units.

These audit cells have now been read for all nine arms. As `DESIGN.md` says, no further variant of an arm is fitted against
them. The 73 reserve donors of the second question were not touched.

## Files

- [`DESIGN.md`](DESIGN.md): the design, the forecast and the rule, committed before the fits.
- `scripts/_common.py`, `fit.py`, `read.py`, `run_all.py`: the run. `scripts/tables.py` prints the tables above from the
  results.
- `results/reading.json`: every number of the reading. `results/fits.json`: the curve, the kept step and the timing of
  each of the 90 fits. The heads are not in git.

The fits were run arm by arm in two parallel queues on two CPU cores (about two minutes a fit, about 85 minutes in all).
Every fit is full-batch from the same start, so the order does not change a head.
