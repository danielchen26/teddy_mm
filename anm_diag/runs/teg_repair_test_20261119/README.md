# A repair computed from the rows of the table, read once on cells that were not used to make it

Status: **a test, read once.** The design, the rule and the forecast are in [`DESIGN.md`](DESIGN.md), committed before the
reading (commit `17038be`). The cells of the reading are 2,000 knockout cells and 4,000 control cells of the melanoma
screen in the co-culture condition. No cell of the reading was used to make or to develop the repair, and the repaired
outputs were read on them once. One cell type, five genes of one pathway, one output (PD-L1). THP-1 is not part of it.

## Correction, added after the reading: these cells had been read before, for other heads

The design says that no run of this repository had read the co-culture cells and that CD274 had not been used anywhere.
That is true of the run folders of this repository and too narrow. Earlier work of this project, made in another
workspace and recorded on the results page and in the ledger of the app (runs 5 to 9 there), read the same cells:

- the same 6,000 co-culture cells for the head as trained (a sealed audit, run 5), for retrained heads (run 6), for the
  head split into stages (run 7) and for another model (run 8);
- the IFN-γ library of 206 targets, CD274 among them, for the head as trained (run 9).

What stays true: none of these cells entered the repair or its development, and this run is the first and only reading
of the repaired outputs on them. What changes:

- "Cells not read before" has to be read as "cells not used to make the repair". The title of this page is changed; the
  design is left as committed, with a note at its top.
- Two results here were already on record for other outputs, and the forecast should have used them. Run 5 found
  that on co-culture cells the embedding itself does not carry all of the change of PD-L1 (*u* 0.19 [0.13, 0.26] at Z,
  rejected; unresolved under the later interval rule). A readout of the embedding was therefore unlikely to come inside
  the margin there, and the 40% given to "holds" in setting 1 was too high. Run 9 found that the knockout of a
  protein's own gene cannot be read from RNA (mean unreadable share 0.92 over 8 such pairs); the probe on CD274 repeats
  that for the repaired output and adds nothing new.

## The question

Does a repair computed from the rows of the table make the audit pass on new cells? And how far beyond its own rows does it
reach: to a gene of the pathway whose row it has not seen, and to a culture condition it has not seen?

The repair: a gate on embeddings decides whether a cell is of the cell type of the rows; inside the gate the head's PD-L1
output is replaced by a linear readout of TEDDY's embedding, fitted on control cells and then updated once by the rows
(for each event: the mean displacement of its knockout cells, the measured change of PD-L1 and its standard error). No
knockout cell enters one by one. `controls` is the same without the rows; `cells` is the reference that sees the measured
protein of every knockout cell.

## Outcome by the rule written first

**Inside the co-culture condition: undecided. From the IFN-γ condition to the co-culture condition: does not hold.**
No setting holds.

**The four settings** (pooled over the five genes; knockout cells resampled by gene)

| setting | the head as trained | `controls` | `rows` | `cells` | outcome of `rows` |
|---|---|---|---|---|---|
| 1. same events, new cells | 1.00 [0.99, 1.00] rejected | 0.27 [0.19, 0.33] rejected | **0.11 [0.01, 0.17] unresolved** | 0.11 [0.01, 0.17] unresolved | **undecided** |
| 2. a new gene | 1.00 [0.99, 1.00] rejected | 0.27 [0.19, 0.33] rejected | **0.11 [0.00, 0.19] unresolved** | 0.11 [-0.01, 0.19] unresolved | **undecided** |
| 3. a new condition | 1.00 [0.99, 1.00] rejected | 0.46 [0.38, 0.51] rejected | **0.43 [0.34, 0.48] rejected** | 0.34 [0.24, 0.40] rejected | **does not hold** |
| 4. a new gene and a new condition | 1.00 [0.99, 1.00] rejected | 0.46 [0.38, 0.51] rejected | **0.43 [0.34, 0.49] rejected** | 0.35 [0.25, 0.42] rejected | **does not hold** |

Settings 1 and 2 fit on one half of the co-culture cells and read the other half. Settings 3 and 4 fit everything in the
IFN-γ condition and read the co-culture cells. *u* is the share of the measured change of PD-L1 that the output does not
account for; the margin is 0.10; a reading is supported when its whole interval is inside the margin.

What the reading says, in order:

1. **The head as trained accounts for nothing of the change in this condition** (1.00).
2. **Five rows do what a thousand knockout cells do.** With the rows the output accounts for all but 0.11 of the change;
   with the measured protein of every knockout cell of the fit, also 0.11. The difference is −0.001 [−0.013, 0.010].
   The rows add 0.16 to what the control cells alone give (0.27).
3. **It is not enough to pass.** The interval of `rows` is [0.01, 0.17], not inside ±0.10. `cells` does not pass either,
   so what is missing is not in the rows.
4. **A gene whose row is not in the table is repaired as well as one whose row is** (0.115 against 0.106), for these five
   genes of one pathway.
5. **The repair does not carry to another condition.** Computed in the IFN-γ condition and read on co-culture cells it
   leaves 0.43, rejected, although the gate is open there (0.99). The rows hardly help (0.46 without them).

**`rows` against the two others** (difference of |u| on the same resamples)

| setting | `rows` − `controls` | `rows` − `cells` | `rows`, knockout cells resampled one by one |
|---|---|---|---|
| 1. same events, new cells | -0.165 [-0.182, -0.150] | -0.001 [-0.013, 0.010] | [0.05, 0.16] unresolved |
| 2. a new gene | -0.156 [-0.187, -0.132] | 0.004 [-0.006, 0.014] | [0.06, 0.17] unresolved |
| 3. a new condition | -0.032 [-0.037, -0.029] | 0.092 [0.082, 0.103] | [0.38, 0.47] rejected |
| 4. a new gene and a new condition | -0.031 [-0.042, -0.024] | 0.079 [0.060, 0.100] | [0.38, 0.47] rejected |

**`rows` by event** (cells resampled one by one)

| setting | IFNGR1 | IFNGR2 | JAK1 | JAK2 | STAT1 |
|---|---|---|---|---|---|
| 1. same events, new cells | 0.11 [0.02, 0.18] unresolved | 0.12 [0.04, 0.20] unresolved | -0.08 [-0.22, 0.03] unresolved | 0.14 [0.06, 0.21] unresolved | 0.19 [0.11, 0.25] rejected |
| 2. a new gene | 0.12 [0.03, 0.19] unresolved | 0.13 [0.05, 0.21] unresolved | -0.10 [-0.24, 0.02] unresolved | 0.15 [0.07, 0.22] unresolved | 0.21 [0.13, 0.27] rejected |
| 3. a new condition | 0.46 [0.39, 0.52] rejected | 0.45 [0.38, 0.51] rejected | 0.26 [0.15, 0.35] rejected | 0.45 [0.39, 0.51] rejected | 0.48 [0.42, 0.53] rejected |
| 4. a new gene and a new condition | 0.46 [0.39, 0.52] rejected | 0.45 [0.38, 0.51] rejected | 0.26 [0.15, 0.35] rejected | 0.45 [0.39, 0.51] rejected | 0.48 [0.42, 0.53] rejected |
| measured change of PD-L1, co-culture | -0.88 | -0.84 | -0.69 | -0.89 | -0.96 |

STAT1 is rejected in every setting; JAK1 is the only gene whose change the repaired output overshoots.

**Other splits into halves** (`rows`, pooled point values; the rule reads seed 0)

| setting | seed 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 1. same events, new cells | 0.106 | 0.114 | 0.073 | 0.106 | 0.101 |
| 2. a new gene | 0.115 | 0.123 | 0.084 | 0.116 | 0.110 |

### The scope of the repair

**Guard and gate**

| setting | PD-L1 on bone-marrow test cells, as trained → with `rows` | fall of the mean over proteins | gate on bone-marrow test cells (mean, highest) |
|---|---|---|---|
| 1. same events, new cells | 0.4077 → 0.4075 | 0.000001 | 0.010, 0.27 |
| 2. a new gene | 0.4077 → 0.4075 | 0.000001 | 0.010, 0.27 |
| 3. a new condition | 0.4077 → 0.4079 | -0.000002 | 0.009, 0.21 |
| 4. a new gene and a new condition | 0.4077 → 0.4079 | -0.000002 | 0.009, 0.21 |

Gate on the co-culture cells of the reading: 0.988 (fitted in the co-culture condition), 0.987 (fitted in the IFN-γ condition); lowest cell 0.63 and 0.73.

The gate is fitted on embeddings alone. It does not tell the two conditions apart: a gate fitted on IFN-γ control cells
is open on co-culture cells, where the repair it guards does not hold.

**The probe: CD274** (287 knockout cells, IFN-γ condition; measured change -1.51)

| the head as trained | `controls` | `rows` | `cells` |
|---|---|---|---|
| 1.00 [0.99, 1.01] rejected | 0.92 [0.88, 0.96] rejected | 0.92 [0.88, 0.96] rejected | 0.92 [0.88, 0.96] rejected |

CD274 is the gene of PD-L1 itself. Its knockout removes more protein than any pathway gene, and no output built on the
embedding sees it. This event is outside what a readout of the embedding can repair. (The same had been found for the
head as trained in earlier work; see the correction above.)

**The 200 other genes** (44,771 knockout cells; measured change: median size 0.044, largest 0.35)

| output | median size of the predicted change | largest | genes beyond 0.2 | correlation with the measured changes |
|---|---|---|---|---|
| as_trained | 0.009 | 0.07 (CCND1) | 0 | 0.13 |
| controls | 0.026 | 0.16 (PFN1) | 0 | 0.66 |
| rows | 0.026 | 0.16 (PFN1) | 0 | 0.66 |
| cells | 0.026 | 0.16 (PFN1) | 0 | 0.67 |

The repaired output does not invent a change where none was measured, and it follows the small measured changes better
than the head as trained. The rows make no difference here. (The measured changes and displacements of these genes had been
tried as rows in development; no output had been read on them.)

## Against the forecast

| | forecast | read |
|---|---|---|
| the head as trained | rejected, above 0.8 | 1.00, rejected |
| 1. same events, new cells | `rows` 0.00 to 0.12; undecided more likely than holds | 0.106, undecided |
| 2. a new gene | within 0.02 of setting 1, same outcome | 0.115, undecided |
| 3. a new condition | 0.05 to 0.30; undecided most likely (60%), does not hold 25% | **0.43, does not hold. The range was wrong** |
| 4. a new gene and a new condition | as setting 3 | **0.43, does not hold. Wrong in the same way** |
| `rows` against `controls` | lower in settings 1 and 2 | −0.165 and −0.156, intervals below zero |
| `rows` against `cells` | at most 0.05 above | −0.001 and 0.004 in settings 1 and 2; **0.09 and 0.08 above in settings 3 and 4** |
| guard | PD-L1 on bone marrow changes by less than 0.005; gate below 0.05 on bone-marrow cells | 0.0002; the gate is 0.01 on average, **0.27 on its highest cell** |
| probe CD274 | rejected, above 0.5 | 0.92, rejected |
| the 200 other genes | median predicted change below 0.05; none beyond 0.5 | 0.026; largest 0.16 |

The forecast inside the condition was right in its range. The forecast across conditions was too hopeful. The 40% for
"holds" in setting 1 did not use what an earlier audit of these cells had found at the embedding (the correction above).

## After the reading (post hoc)

The co-culture cells had been read when these were computed. Nothing below is evidence; it says where to look.

**Post hoc 1: the cells of the fit against the cells of the reading** (setting 1, pooled point values)

| | `controls` | `rows` | `cells` |
|---|---|---|---|
| on the half that was fitted | 0.17 | -0.03 | 0.01 |
| on the other half (the reading) | 0.27 | 0.11 | 0.11 |

On the cells it was fitted on, the repair by rows leaves nothing. The 0.11 appears when the readout is carried to cells it
has not seen, and it appears for `cells` in the same size. It is not something the rows lack. Part of it is a cost of
fitting a readout of 512 coordinates on a few thousand cells; part of it may be the embedding itself, which an earlier
audit found short on these cells (the correction above). This computation cannot separate the two.

**Post hoc 2: by ridge weight** (setting 1, pooled point values; the rule chose 0.32 and 0.1 in the two halves)

| ridge weight | `controls` | `rows` | `cells` |
|---|---|---|---|
| 0.001 | 0.42 | 0.29 | 0.26 |
| 0.00316 | 0.40 | 0.27 | 0.25 |
| 0.01 | 0.36 | 0.23 | 0.22 |
| 0.0316 | 0.32 | 0.18 | 0.18 |
| 0.1 | 0.28 | 0.13 | 0.13 |
| 0.316 | 0.26 | 0.08 | 0.08 |
| 1 | 0.31 | 0.05 | 0.05 |
| 3.16 | 0.44 | 0.08 | 0.07 |
| 10 | 0.63 | 0.18 | 0.18 |
| 31.6 | 0.81 | 0.34 | 0.37 |
| 100 | 0.90 | 0.49 | 0.56 |

The rule (cross-validation of the control fit) chose 0.32 and 0.1. At a weight of 1 the reading would have been 0.05, for
`rows` and for `cells` alike; the same was true of the IFN-γ cells in development (0.08 at 1 against 0.09 at 0.32). A
weight chosen for predicting PD-L1 among control cells is smaller than the weight that carries a change best. A rule that
chooses the weight for the change would have to be stated before another test and read on cells not read before.

Post hoc 3, five folds: `controls` 0.25 [0.17, 0.32] rejected; `rows` 0.09 [-0.01, 0.17] unresolved; `cells` 0.09 [0.00, 0.16] unresolved


**Post hoc 4: the two conditions**

| gene | measured change, IFN-γ | co-culture | cosine of the two displacements | length, IFN-γ | co-culture |
|---|---|---|---|---|---|
| IFNGR1 | -0.90 | -0.88 | 0.82 | 0.119 | 0.113 |
| IFNGR2 | -0.82 | -0.84 | 0.85 | 0.116 | 0.108 |
| JAK1 | -0.69 | -0.69 | 0.84 | 0.116 | 0.116 |
| JAK2 | -0.91 | -0.89 | 0.80 | 0.115 | 0.112 |
| STAT1 | -0.96 | -0.96 | 0.83 | 0.126 | 0.115 |

Distance between the control means of the two conditions: 0.25; mean length of a displacement: 0.11. Cosine of the readouts fitted in the two conditions: `controls` 0.36, `rows` 0.38, `cells` 0.43.
Mixed repair (control cells of the co-culture condition, rows of the IFN-γ condition): 0.16 [0.07, 0.23] unresolved; by event IFNGR1 0.17, IFNGR2 0.18, JAK1 -0.02, JAK2 0.19, STAT1 0.24.

The rows of the two conditions look alike: each gene changes PD-L1 by the same amount, and its knockout displaces the
embedding in nearly the same direction with the same length. What differs is the baseline: the control cells of the two
conditions lie twice as far apart as a knockout moves a cell, and the readouts fitted on them point in different
directions. With the control cells of the new condition and the rows of the old one, the reading is 0.16: between
`controls` alone (0.27) and the condition's own rows (0.11). A row carries part of the way across conditions once the
baseline is fitted where it is used.

## What this shows, and what it does not

For one cell type, five genes of one pathway and one output, on cells not read before:

- A repair can be computed from the table alone. Five rows and the control cells bring the unexplained share from 1.00 to
  0.11, the same as the measured protein of a thousand knockout cells, with bone marrow unchanged and no invented change
  for 200 other genes.
- By the rule it does not pass: the interval is not inside the margin. What is missing is the same with rows and with
  cells.
- A row belongs to its condition. Computed in one condition, the repair is rejected in the other; the gate, which sees
  only the cell type, does not notice.
- An event of another kind (the knockout of the protein's own gene) is not reached at all.

It does not show:

- that the repair is better than another method: `cells` is a reference inside the same readout, and the fine-tuned heads
  of the earlier runs were not read here;
- that the repair lives in the head: one output of the head is replaced, inside one cell type, by a readout of the
  embedding;
- anything about another cell type, pathway or protein. With five genes as resampling units the pooled intervals are
  approximate.

## What a next test would need

- For passing: a rule for the ridge weight aimed at the change (stated first), and cells not read before. The melanoma
  cells with embeddings are now all read; the rest of the THP-1 screen needs embeddings first. If the embedding is what
  falls short in this condition, no readout of it will pass, and the repair has to move below the head.
- For the condition: the condition as part of the type of a row and of the gate, and the control cells of the new
  condition as the baseline.
- For CD274: another input than the embedding. This row goes to the part of the loop that asks what to measure.

## Files

- [`DESIGN.md`](DESIGN.md): the design, the rule and the forecast, committed before the reading.
- `scripts/cells.py`: the features of the cells (no model output). `scripts/test.py`: the repair, the development run on
  cells already read (`results/development.json`) and the one reading (`results/reading.json`).
- `scripts/develop_rows.py`, `develop_weight.py`: development before the design (`results/develop_rows.json`,
  `develop_weight.json`). `scripts/posthoc.py`: after the reading (`results/posthoc.json`). `scripts/tables.py` prints the
  tables above.
- The embeddings and the protein file are not in git.
