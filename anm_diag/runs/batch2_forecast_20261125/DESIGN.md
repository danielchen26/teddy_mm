# Batch 2: forecasts written first, then 33 units never read are opened once

25 November 2026. The owner approved three things after run 41 (`runs/account_20261125`): the readings of that run
enter the method; a small batch of units never read is opened, with the forecasts committed first; nothing is trained
now, and whether to train is decided at the end, on the evidence. This document is the design of the second and the
third. It is committed, with the forecasts, before any protein value of the batch is opened.

Nothing is trained. No final unit and no reserve donor is touched.

## 1. The question

Run 41 rebuilt the reading of a change (three readings, four parts) and the law of the rows of a library (a readout of
control cells, updated by the recorded events, with edges between an event and a readout). Both were built on cells that
earlier runs had read. Two things have never been tested on units nobody had read:

1. Does the law forecast the change of a unit **from the RNA of its cells alone**, with intervals that hold?
2. Does the account place a change where the design says it must be: for a knockout of the gene behind the measured
   protein, *not in the input, and not forecast by other events*, so that only an edge can carry it?

A third question is the owner's, for the decision on training: **how often is a change lost inside the model**, between
the input and the state, and how much of it does a state typed by the TEG bring back? Run 41 measured that on one event
(IRF1) and on fifteen melanoma accounts. This batch adds the units of batch 1 and of batch 2.

## 2. The units

`results/roles.json`, written by `scripts/roles.py` from the census of 22 November
(`runs/outer_step_20261122/results/census.json`). 33 units, every one *unread and not final* in the census. The script
stops if a unit is final, was read, or is unknown. The 133 final units stay unopened; 225 unread units that are not
final are left after this batch.

| stratum | targeted, and why | calibration |
|---|---|---|
| melanoma, co-culture | CD274 (own gene of PD-L1); HLA-B, HLA-C (own genes of MHC-I); HLA-E (own gene of HLA-E); STAT3 (flagged on MHC-I by the law of the library under interferon); PFN1 (named in the plan; no edge, no flag) | CTSA, CTSL, STX7, TOP1MT |
| melanoma, no signal given | CD274; HLA-A, HLA-B; HLA-E (own genes); SSR2, CD44 (flagged on MHC-I under interferon); IFNGR1, IFNGR2, JAK2 (the pathway, with no interferon given) | S100A6, SCARB2, SLC19A1, TMEM173 |
| THP-1 | none (the screen holds no knockout of CD274; its own-gene units are final) | ATF2, CAV1, ETV7, IRF7, MARCH8, MYC, SMAD4, STAT5A, TNFRSF14, UBE2L6: every unit left |

Calibration units in melanoma are drawn by a seeded rule (the docstring of `roles.py`): four per condition from the
units that are unread, not final, not of batch 1 and not targeted. Readouts: PD-L1, MHC-I and HLA-E in melanoma, PD-L1 in
THP-1. That is 23 × 3 + 10 = **79 pairs** of an event and a readout.

No new control cell is opened. The control cells are those already read: 5,276 in co-culture (4,000 of development and
1,276 of batch 1), 2,645 in the condition with no signal (batch 1), 2,386 in THP-1.

## 3. What was known before, said plainly

- **A screen of 13 November read the mean shift of every measured protein for every unit, from the protein files alone**
  (`runs/loop_reading_20261112/results/data_census.json`, in this repository). It lists "clear responses". So the *size*
  of the change of the own-gene units is not unknown: for instance CD274 in co-culture is listed with −2.1 spreads on
  its antibody, and HLA-B, HLA-C and HLA-E are listed as clear there. The author of this design looked at the
  co-culture list while writing it, and not at the list of the condition with no signal. The census of 22 November says
  that none of the ten THP-1 units is listed as clear.
- The forecasts of Section 5 are the output of an algorithm that takes nothing from that screen. Its parts were committed
  in run 41 (the law, the table of edges); what is added here is written in Section 5 before it is run.
- What nobody has read: the protein of these cells **cell by cell against their states**. The account (Section 6) and the
  calibration of the intervals are new. The size forecasts of the eight own-gene pairs are a test of the algorithm, not a
  blind prediction, and are reported as such.
- Batch 1 (23 November) read eight units and found three changes "open before the model": CMTM6 in THP-1, TAPBP in
  co-culture, B2M in the condition with no signal. Those units are recorded rows here.

## 4. The cells and their states (RNA only)

`scripts/cells.py` lists the cells of the batch from the record of the screens (names, guides, condition; no protein
value) and writes through the cells read before: 31,298 cells, 11,061 of them of the batch. `scripts/states.py` reads
the RNA counts of all of them and runs the encoder as released: the state (the mean of the last layer over a cell's gene
positions) and the blocks of the typed state for the sets of genes of run 40 (`A` and `R` of each half, chosen there on
the fitted half of the development cells with tokens and labels only). States of cells read before are computed again,
so that every cell of a stratum has its state from one run; the script compares them with the stored states and writes
nothing if the largest difference is above 1e-4.

## 5. The forecasts

`scripts/forecast.py`, run once, writes `results/forecasts.json`. In spreads of a readout among the control cells of
the stratum:

    forecast(event, readout) = shift of the state . control readout          (read)
                             + what the recorded events add through the state (recorded)
                             + the part by its edge                           (edge)

- **Per stratum, the law of its recorded rows** (`anm_diag.account.Law`, at the state at unit length): the control
  cells, and the recorded events with their measured changes: in co-culture the five pathway events and CD44, TAPBP,
  CDH19; in the condition with no signal B2M and KLF4; in THP-1 the four pathway events, IRF1, BRD4, CMTM6 and POU2F2.
- **A unit of the batch** is forecast by `Law.new` (added to the module with this design, with a test: on the cells of an
  event it gives the forecast of that event left out). The error of a row not yet measured is taken at the spread among
  control cells.
- **The edges known beforehand** are the table of run 41: the event's site is a gene of the measured protein (CD274 for
  PD-L1; HLA-A, HLA-B, HLA-C, B2M for MHC-I; HLA-E, B2M for HLA-E).
- **What a stratum cannot tell from its own rows it takes from the library under interferon** (206 events, all
  development). The spread of a direct part without an edge: from the stratum's own rows when it has at least eight
  recorded pairs without an edge (co-culture, THP-1), else from the library (the condition with no signal). The class
  of edges: the direct parts of the recorded pairs of the class, from the library and from the strata of the same
  study, pooled; a new pair of the class is forecast at their mean, with their spread.

### The laws as fitted (recorded rows only; no value of the batch)

| stratum | control cells | recorded events | spread of a direct part without an edge | widening of the control prior | half-width of a 95% interval, pairs without an edge | recorded pairs, each left out, inside |
|---|---|---|---|---|---|---|
| co-culture | 5,276 | 8 | 0.041 (own rows) | PD-L1 63, MHC-I 2,000, HLA-E 5 | PD-L1 0.28, MHC-I **1.2**, HLA-E 0.18 | 23 of 24 (TAPBP on MHC-I is outside) |
| no signal | 2,645 | 2 | 0.029 (the library) | PD-L1 1, MHC-I 50, HLA-E 250 | PD-L1 0.19, MHC-I 0.33, HLA-E 0.51 | 4 of 4 without an edge |
| THP-1 | 2,386 | 8 | 0.19 (own rows) | PD-L1 25 | 0.41 | 7 of 8 (CMTM6 is outside) |

The class of edges: nine recorded direct parts (seven of the library, two of B2M with no signal), between −0.07 and
−2.8, mean **−0.99**; a new pair of the class is forecast there with a spread of 1.03, so its interval is about −3 to +1.

Two things to keep in mind when the batch is read. In co-culture the recorded pathway events tell the law that the control
readout of MHC-I is far from enough, and it widens that prior two-thousand-fold: an interval of ±1.2 spreads says little,
and HB1 is easy there. With no signal the spread without an edge is the library's 0.029: if the intervals fail, they
fail there.

Before a value is opened the forecasts are committed and pushed, and the owner is told what they are.

## 6. After the commit: the values are opened once, and read

`scripts/pack.py` exports the readouts of the cells of the batch (as batch 1: in melanoma the centred log ratio over
the 24 antibody columns; in THP-1 the log of one plus the PD-L1 counts). It refuses to run without the committed file
of forecasts, and adds the values of 300 cells read before per study, which `read.py` requires to be identical to the
values the earlier runs read. Then `scripts/read.py`, once (it refuses if its output exists):

1. **Every pair**: the measured change in the law's units and its error; whether it is inside the forecast interval; the
   error of the forecast against two forecasts that use less: no change, and the control readout alone.
2. **The account** (`anm_diag.account.reading` and `chain`, as run 41) of every pair of the batch, and of every recorded
   pair of batch 1, whose change is at least four errors: at the input view `V` (512 components of the RNA of the
   stratum's cells), the state `Z`, and, in THP-1 and co-culture, the typed state `T` and the control `R` (genes at
   random). `others`: the recorded events of the stratum. The halves of the development cells are those of run 40, so
   that a typed set was never chosen on the half it is read on; the other cells are split by `cycle.halves`.

`scripts/score.py` compares with `results/declared.json`.

## 7. Hypotheses

A step is 0.10. "Inside" means inside the 95% interval of the forecast.

**HB1. The intervals hold.** Of the 79 pairs, at least 72 (90%) are inside. *Refuted* at 67 or fewer; between, unresolved.

**HB2. The forecast uses what it adds.** Over the 79 pairs the root mean square error of the forecast is smaller than
that of "no change" and than that of the control readout alone. *Refuted* if smaller than neither.

**HB3. An edge known in one condition forecasts in another.** Of the eight own-gene pairs (four per condition), at
least six are inside, and their error with the edge is smaller than without it. The intervals of this class are wide
(its spread is about one spread of the readout), and the sign was known from the screen: this is a weak test, and the
table of the eight is what to read.

**HB4. The account places an own-gene change outside the input.** For the own-gene pairs whose change is at least four
errors (at least four such pairs are needed): the control readout at the state reads less than half (upper end of its
interval below 0.50), the largest part of the chain is "not in the input", and the recorded events do not forecast it
(`others` at the state: upper end below 0.50). *Confirmed* if that holds for at least three quarters of them, *refuted*
below half.

**HB5. A flag made in one condition holds in another.** STAT3 on MHC-I (here in co-culture), SSR2 and CD44 on MHC-I (here
with no signal) were flagged by the law of the library under interferon. With no edge given, at least two of the three
fall outside their interval. *Refuted* if none does.

**HB6. The typed state where the mean state loses.** Over the pairs of point 2 of Section 6 in THP-1 and co-culture that
are not own-gene pairs: (a) for no pair is the interval of `carried(T) − carried(Z)` entirely below −0.10; (b) in every
pair where the state carries a step less than the input (interval of `carried(V) − carried(Z)` above 0.10), the typed
state brings back at least half of it. If no pair loses a step, (b) is *not tested*, and that count is itself what the
owner's decision on training needs.

Forecasts of these, with probabilities, are in `results/declared.json`.

## 8. What this can and cannot show

- 79 pairs, most of them expected to show no change: HB1 is mostly about intervals around zero. The pairs with a change
  are few, and are listed one by one.
- The condition with no signal has two recorded events. Its spread without an edge is the library's (a different
  condition). If the intervals fail there, that is the finding: the spread does not carry across conditions.
- The class of edges has nine recorded pairs. Its forecast says "a drop of about one spread, anywhere between none and
  three". An edge with a size of its own (how much of the protein a knockout removes) is not built.
- THP-1 has no own-gene unit in the batch: the edge across cell types is not tested.
- A change is read for a pair when its measured change is at least four errors, on the same cells; the shares of a pair
  just above that line lean high.

## 9. Order

1. This document, `roles.json`, the scripts, `forecasts.json`, `declared.json`: committed and pushed.
2. The owner is told the forecasts.
3. `pack.py` on the owner's machine; `read.py`; `score.py`.
4. Write-up; the method's text (`LOOP_PROTOCOL.md`) with the readings of run 41; the app.
