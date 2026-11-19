# A repair computed from the rows of the TEG: the gate — 19 November 2026

**Development, not evidence.** Every cell used here had been read for
[`runs/teg_training_20261119`](../teg_training_20261119/); with the owner's consent they are used once more, as development
cells. The design, the decision rule and a forecast are in [`DESIGN.md`](DESIGN.md), committed before the computation
(commit `f85100e`). Nothing below is a test of anything.

## The question

With only the rows of the events of the fit (one row per event: the measured change of PD-L1, its sign, size and
standard error) and the control cells, does a repair carry to a gene it never saw, as far as a repair that sees the
measured protein of every knockout cell?

## Outcome by the rule written first: undecided

**The gate, place A** (the PD-L1 output as a linear readout of the embedding; pooled)

| given | melanoma | THP-1 | melanoma: \|*u*\| minus `cells` | THP-1: \|*u*\| minus `cells` |
|---|---|---|---|---|
| `controls` | 0.15 [0.08, 0.21] unresolved | 0.45 [0.34, 0.59] rejected | 0.08 [0.05, 0.11] | 0.30 [0.24, 0.35] |
| `rows` | 0.14 [0.07, 0.20] unresolved | 0.20 [0.02, 0.37] unresolved | 0.06 [0.04, 0.09] | 0.05 [−0.02, 0.10] |
| `rows_hard` | 0.16 [0.09, 0.22] unresolved | 0.19 [0.02, 0.38] unresolved | 0.08 [0.06, 0.11] | 0.04 [−0.02, 0.10] |
| `rows_guides` | 0.13 [0.06, 0.20] unresolved | 0.18 [0.02, 0.36] unresolved | 0.06 [0.04, 0.08] | 0.03 [−0.02, 0.08] |
| `cells` | 0.07 [−0.01, 0.15] unresolved | 0.15 [0.01, 0.33] unresolved |  |  |

`rows` is `controls` updated in closed form by one row per event of the fit. By the rule (the gap between `rows` and
`cells`, point values): THP-1 is a **go** (0.05), melanoma is **undecided** (0.06, between the go line of 0.05 and the
no-go line of 0.10). Overall: **undecided**.

- **THP-1.** Four rows take the readout from 0.45 to 0.20; every knockout cell takes it to 0.15. The rows carry.
- **Melanoma.** The control cells alone already give 0.15. The rows add 0.02; the cells add 0.08.
- **The requirement.** No reading is supported. `cells`, the most a linear repair can have, is unresolved in both data
  sets (its point value is inside the margin on melanoma).

The forecast: `cells` and `rows` came out inside their forecast ranges in both data sets, and `controls` on THP-1 did.
`controls` on melanoma was forecast at 0.3 to 0.6 and is 0.15. The decision was forecast as go on melanoma and undecided
on THP-1; it is the other way round. Place B was forecast at most 0.05 above place A and is far above it (below).

**By event, place A** (*u*; one gene, the fit that never saw it, all audit controls)

| given | mel IFNGR1 | mel IFNGR2 | mel JAK1 | mel JAK2 | mel STAT1 | THP-1 IFNGR1 | THP-1 IFNGR2 | THP-1 JAK2 | THP-1 STAT1 | THP-1 IRF1 |
|---|---|---|---|---|---|---|---|---|---|---|
| `controls` | 0.13 | 0.15 | 0.11 | 0.21 | 0.15 | 0.41 | 0.42 | 0.46 | 0.21 | 0.88 |
| `rows` | 0.11 | 0.13 | 0.08 | 0.20 | 0.14 | 0.14 | 0.19 | 0.24 | −0.25 | 0.73 |
| `cells` | 0.06 | 0.09 | −0.03 | 0.15 | 0.08 | 0.09 | 0.13 | 0.15 | −0.13 | 0.64 |

IRF1 is carried by neither: its displacement lies mostly outside the span of the other four (cosine 0.47), where every
other held-out gene lies inside it (0.97 to 1.00).

**`rows` by the number of events of the fit** (place A; mean \|*u*\| on the held-out event, every subset)

| events of the fit | 0 (`controls`) | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| melanoma | 0.15 | 0.14 | 0.13 | 0.13 | 0.13 |
| THP-1 | 0.47 | 0.34 | 0.32 | 0.31 | 0.31 |

**The held-out displacement against the span of the events of the fit** (cosine)

- melanoma: IFNGR1 0.97, IFNGR2 0.98, JAK1 0.97, JAK2 0.98, STAT1 0.98
- THP-1: IFNGR1 1.00, IFNGR2 1.00, JAK2 1.00, STAT1 0.97, IRF1 0.47

One row gives most of what four give. A further row of the same kind adds little.

**The gate, place B** (the head as trained plus a linear correction, bone-marrow term at weight 1; read as the predicted mean, as designed)

| given | melanoma | THP-1 | bone marrow, PD-L1 score (0.408 as trained) |
|---|---|---|---|
| `controls` | 0.68 [0.63, 0.71] rejected | 0.77 [0.73, 0.85] rejected | 0.401 and 0.376 |
| `rows` | 0.62 [0.58, 0.67] rejected | 0.56 [0.49, 0.65] rejected | 0.400 and 0.377 |
| `cells` | 0.62 [0.57, 0.66] rejected | 0.54 [0.48, 0.62] rejected | 0.399 and 0.397 |

Place B, as designed, does not work: the head as trained with a linear correction added to its PD-L1 output stays above
0.5. The two sections below look at why.

## After the gate was read (post hoc)

### The scale on which the own output is read

The measured protein is on a log scale (melanoma: controls 1.21 ± 0.74; THP-1: 5.08 ± 0.79). The reading of the head's own
output introduced in [`runs/scope_20261116`](../scope_20261116/) (`frozen_arms.py`) and used by the rule of
`runs/teg_training_20261119` takes the head's predicted mean, the exponential of its output, and compares it with the
measurement through one slope fitted among the control cells. The statistic does not change under a linear change of the
prediction; it does change under an exponential. `scripts/scale.py` reads the heads already made on both scales, on the
same cells and resamples. No head is fitted.

**The heads of `runs/teg_training_20261119` on two scales** (post hoc; pooled *u* of the own PD-L1 output)

| arm | melanoma, as the predicted mean (as recorded) | melanoma, on the log scale | THP-1, as the predicted mean | THP-1, on the log scale |
|---|---|---|---|---|
| `as_trained` | 0.92 [0.90, 0.93] rejected | 0.91 [0.89, 0.93] rejected | 0.97 [0.96, 0.99] rejected, wrong sign | 0.97 [0.95, 0.98] rejected, wrong sign |
| `data` | 0.30 [0.23, 0.35] rejected | 0.20 [0.13, 0.26] rejected | 0.44 [0.37, 0.53] rejected | 0.24 [0.13, 0.36] rejected |
| `class` | 0.48 [0.42, 0.54] rejected | 0.42 [0.36, 0.49] rejected | 0.97 [0.92, 1.08] rejected | 0.91 [0.86, 1.01] rejected |
| `edge1` | 0.33 [0.27, 0.39] rejected | 0.26 [0.19, 0.33] rejected | 0.57 [0.50, 0.67] rejected | 0.56 [0.40, 0.70] rejected |
| `edge05` | 0.42 [0.38, 0.46] rejected | 0.35 [0.30, 0.40] rejected | 0.67 [0.60, 0.79] rejected | 0.62 [0.53, 0.76] rejected |
| `edge2` | 0.32 [0.23, 0.40] rejected | 0.23 [0.14, 0.32] rejected | 0.61 [0.53, 0.70] rejected | 0.56 [0.37, 0.70] rejected |
| `graph` | 0.32 [0.26, 0.38] rejected | 0.25 [0.18, 0.32] rejected | 0.59 [0.51, 0.69] rejected | 0.48 [0.38, 0.64] rejected |
| `judge` | 0.37 [0.31, 0.43] rejected | 0.32 [0.25, 0.38] rejected | 0.32 [0.20, 0.46] rejected | 0.11 [−0.08, 0.30] unresolved |
| `rows` | 0.24 [0.15, 0.32] rejected | 0.21 [0.10, 0.29] rejected | 0.27 [−0.12, 0.53] unresolved, wrong sign | 0.34 [0.05, 0.57] unresolved, wrong sign |
| `data_rows` | 0.22 [0.14, 0.29] rejected | 0.16 [0.07, 0.23] unresolved | 0.64 [0.43, 0.96] rejected, wrong sign | 0.67 [0.47, 1.01] rejected, wrong sign |

- Read on the scale of the measurement, the repair already made leaves 0.20 and 0.24, not 0.30 and 0.44. The numbers on
  record for the own output of the earlier repairs (0.29 and 0.43) are on mismatched scales in the same way.
- The outcome of `runs/teg_training_20261119` by its rule does not change: on neither scale is any arm supported in
  both data sets.
- One statement of that run does change. "On melanoma `rows` is lower than every other single kind of knowledge" holds
  on the scale of its rule (0.24 against 0.30 for `data`) and not on the log scale (0.21 against 0.20). `rows` was trained
  on the statistic as its rule read it.
- In this gate, place A is on the scale of the measurement and is not affected. Place B was read as the predicted mean.

### Why place B fails, and a correction scoped to the cell type of its rows

**A correction scoped to the cell type of its rows** (post hoc; pooled *u* on the log scale; bone marrow: PD-L1 score, 0.408 as trained)

| correction | given | melanoma | THP-1 | bone marrow |
|---|---|---|---|---|
| added to the head's output, no bone-marrow term | `controls` | 0.52 [0.45, 0.57] rejected | 0.51 [0.42, 0.63] rejected | 0.354 and 0.037 |
| added to the head's output, no bone-marrow term | `rows` | 0.51 [0.44, 0.57] rejected | 0.34 [0.20, 0.48] rejected | 0.354 and 0.037 |
| added to the head's output, no bone-marrow term | `cells` | 0.51 [0.43, 0.57] rejected | 0.31 [0.19, 0.46] rejected | 0.345 and 0.077 |
| added to the head's output, bone-marrow term at weight 1 (place B) | `controls` | 0.64 [0.59, 0.69] rejected | 0.69 [0.64, 0.79] rejected | 0.401 and 0.376 |
| added to the head's output, bone-marrow term at weight 1 (place B) | `rows` | 0.58 [0.51, 0.63] rejected | 0.37 [0.26, 0.50] rejected | 0.400 and 0.377 |
| added to the head's output, bone-marrow term at weight 1 (place B) | `cells` | 0.58 [0.51, 0.63] rejected | 0.38 [0.28, 0.49] rejected | 0.399 and 0.397 |
| added to the head's output, only in this cell type | `controls` | 0.52 [0.45, 0.57] rejected | 0.52 [0.43, 0.64] rejected | 0.408 and 0.408 |
| added to the head's output, only in this cell type | `rows` | 0.51 [0.44, 0.57] rejected | 0.33 [0.20, 0.47] rejected | 0.408 and 0.408 |
| added to the head's output, only in this cell type | `cells` | 0.51 [0.43, 0.57] rejected | 0.32 [0.20, 0.46] rejected | 0.408 and 0.409 |
| replaces the head's output, only in this cell type | `controls` | 0.15 [0.08, 0.21] unresolved | 0.46 [0.36, 0.60] rejected | 0.408 and 0.410 |
| replaces the head's output, only in this cell type | `rows` | 0.13 [0.06, 0.20] unresolved | 0.19 [0.03, 0.37] unresolved | 0.408 and 0.411 |
| replaces the head's output, only in this cell type | `cells` | 0.07 [−0.01, 0.14] unresolved | 0.17 [0.04, 0.33] unresolved | 0.408 and 0.411 |

- **A readout fitted on one cell type cannot replace the output everywhere.** Applied to the bone-marrow test cells, the
  readout of place A scores 0.02 (melanoma fit) and 0.06 (THP-1 fit, `rows`) for PD-L1, against 0.408 for the head as
  trained.
- **Adding a linear correction to the head's output does not work**, with or without the bone-marrow term: the head's
  output on these cells is not a linear function of the embedding, and what the correction cannot cancel stays in.
- **Replacing the output only in the cell type of the rows does.** A gate fitted on embeddings alone (the training control
  cells of the data set against the bone-marrow training cells; no protein and no audit cell) is 0.99 and 0.96 on the
  knockout cells and 0.01 on the bone-marrow test cells. With it, the readings are those of place A, and the head's PD-L1
  score on bone marrow stays at 0.408 to 0.411.

On the same cells and the same scale, the scoped repair computed from rows alone reads 0.13 and 0.19; the head repaired
by fine-tuning its MLP on the measured protein of every cell (`data`) reads 0.20 and 0.24. Both are development numbers.
Most of that difference is not owed to the rows: with the control cells alone the scoped readout is at 0.15 on melanoma.

### The ridge weight

**The ridge weight** (post hoc; place A, pooled point values; the gate used 0.1, chosen by cross-validation on training cells)

| ridge weight | melanoma: `controls` | `rows` | `cells` | THP-1: `controls` | `rows` | `cells` |
|---|---|---|---|---|---|---|
| 0.001 | 0.32 | 0.33 | 0.24 | 0.44 | 0.36 | 0.21 |
| 0.01 | 0.24 | 0.25 | 0.18 | 0.42 | 0.28 | 0.19 |
| 0.1 | 0.15 | 0.14 | 0.07 | 0.45 | 0.21 | 0.15 |
| 1.0 | 0.21 | 0.09 | 0.02 | 0.52 | 0.18 | 0.17 |
| 10.0 | 0.47 | 0.15 | 0.11 | 0.65 | 0.28 | 0.31 |
| 100.0 | 0.80 | 0.40 | 0.44 | 0.71 | 0.41 | 0.47 |

The gate chose 0.1 in every fit, by cross-validation on training cells; the cross-validated error at 1.0 is nearly the
same. Between 0.1 and 10 the picture holds: `rows` within about 0.06 of `cells`, `controls` further away. Below 0.1 every
reading is worse and on melanoma the rows no longer help. The numbers of this gate depend on this choice more than a test
should.

## What this shows, and what it does not

For these two data sets, five genes of one pathway each, as development:

- A repair can be computed from the table. The control cells give the baseline and the metric, each row is one update,
  the cell type of the rows is the scope of the repair, and the audit reads the result. On THP-1 it carries to an unseen
  gene almost as far as a repair that sees every knockout cell.
- It carries inside a class of events and not beyond: IRF1 is reached by neither the rows nor the cells.
- What the rows add depends on the data set: a lot on THP-1, little on melanoma, where the control cells of the cell type
  already carry most of the response.
- The own output has to be read on the scale of the measurement.

It does not show:

- that a repair from rows is better than another method. Every cell here had been read before, the ridge weight matters,
  and the comparison with the fine-tuned head is not like for like (a readout of the embedding against a head);
- that the requirement can be met: no reading is supported, with five genes as resampling units;
- anything about events of another kind, or about a repair inside the head's MLP.

## What a test would need

Cells not yet read (the rest of the THP-1 screen, which also has genes that do not move PD-L1 and genes that raise it;
or the melanoma co-culture cells), the own output read on the scale of the measurement, the ridge weight fixed by a rule
that the design states, the scoped replacement as the repair, and one reading.

## Files

- [`DESIGN.md`](DESIGN.md): the design, the decision rule and the forecast, committed before the computation.
- `scripts/gate.py`: the gate (`results/gate.json`). `scripts/scale.py`, `posthoc.py`, `sensitivity.py`: the three post hoc
  checks (`results/scale.json`, `posthoc.json`, `sensitivity.json`). `scripts/tables.py` prints the tables above.
- `scale.py` reads the heads of `runs/teg_training_20261119`, which are not in git.
