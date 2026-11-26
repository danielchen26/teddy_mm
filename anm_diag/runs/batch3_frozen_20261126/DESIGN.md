# Batch 3: one turn of the loop with the rules frozen

26 November 2026. The owner approved the route of the evening of 25 November: fix the rules, read one batch under
them, then the final audit. This document is the design of the batch. It is committed before any value of the batch
is opened; the forecasts are committed before the values they forecast.

Nothing is trained. No final unit and no reserve donor is touched.

## 1. The question

A loop is closed when an error produces a correction by rule, the correction makes the error smaller on data it never
saw, and the rules do not move meanwhile. Batch 2 showed the second once: with the event graph the forecasts were
better than the model alone (0.231 against 0.292). It also showed faults, and the corrections of the afternoon were
made on the batch itself. So:

1. Does the system of version 5, **with no rule changed**, forecast values never read better than the system as it
   stood for batch 2, and better than the model alone, with intervals that hold?
2. Does each lever do what the method says of it, on values never read: an edge with a size (lever G), the records of
   a type (lever E), the state kept by type (lever M)?

The rules are those of commit `d777070` ([`runs/freeze_v5_20261126`](../freeze_v5_20261126/)). If the batch forces a
rule to change, it counts as development, and the final audit is not ready.

## 2. What is opened

**Seven units never read** ([`results/roles.json`](results/roles.json), 1,396 cells): in the co-culture CD47, CD58,
CD59 (the own gene of a readout never read), ITGA3 (no declared edge) and NONO (an edge to PD-L1 that the library's
law proposed in run 43 and that had no row elsewhere); with no signal given CD47 and NONO. `roles.py` stops if a unit
is final or was read. After the batch 218 units that are unread and not final are left. The 133 final units and the
73 reserve donors stay unopened.

**On cells read before, the antibody columns nobody read cell by cell**: 17 of the 20 antibodies of the melanoma
panel (PD-L1, MHC-I and HLA-E were the readouts so far; the four isotype antibodies give the background), and 3 of
the 4 of the THP-1 panel (PD-L2, CD86, CD366).

**What was known before, said plainly.** The screen of 13 November read the mean shift of every antibody for every
unit, from the protein files alone (`runs/loop_reading_20261112/results/data_census.json`). It lists a clear change
for 56 pairs on these columns among the cells read before, and for 5 pairs of the seven units (CD47, CD58 and CD59
on their own antibody; ITGA3 on CD49f). So the size of a change is not unknown, as in batch 2. The forecasts are the
output of code that takes nothing from that screen. What nobody has read: these columns cell by cell against the
states.

## 3. Two stages

`scripts/pack.py`, on the owner's machine.

**Stage A, what the forecasts stand on** (opened first, no forecast is made of it): the control cells of every
stratum; the library under interferon (206 events, development); and the *teachers*: so that the records of a type
can be tested on the new readouts, the recorded pathway events of a stratum are split by a seeded rule (the
docstring of `roles.py`). Teachers: IFNGR2, JAK2 and STAT1 in the co-culture; JAK2 and STAT1 in THP-1.

**Stage B, what is forecast**: every other cell. The script refuses to run without the committed forecasts.

| stratum | units forecast | readouts | pairs |
|---|---|---|---|
| co-culture | IFNGR1, JAK1 (pathway); CD44, TAPBP, CDH19 (batch 1); the ten of batch 2; the five new | 17 new; for the five new units also the 3 read before | 20 × 17 + 5 × 3 = 355 |
| no signal | B2M, KLF4 (batch 1); the thirteen of batch 2; the two new | 17 new; for the two new units also the 3 read before | 17 × 17 + 2 × 3 = 295 |
| THP-1 | IFNGR1, IFNGR2 (pathway); IRF1; BRD4, CMTM6, POU2F2 (batch 1); the ten of batch 2 | 3 new | 16 × 3 = 48 |

698 pairs of an event and a readout.

## 4. Readouts, edges, types (fixed in `runs/freeze_v5_20261126/scripts/_v5.py`)

- A readout of the melanoma screen is the centred log-ratio of one antibody over the 24 columns; of THP-1, the log of
  one plus the counts. Changes are in spreads of the readout among the control cells of the stratum.
- Declared edges: the site is a gene of the protein the antibody measures (the table in `_v5.py`), and B2M for MHC-I
  and HLA-E. Nothing else is declared. The edges admitted in run 43 are in the graph for the readouts read before.
- The scale of an edge (rule 4k) and whether it measures anything (rule 5c) are decided from control cells in stage A.
- One declared type: the five events of the interferon pathway.

## 5. Three systems, all forecast before stage B is opened

`scripts/forecast.py` writes `results/forecasts.json`, once.

| | what it is | how |
|---|---|---|
| S0 | the model alone | the state of the unit's cells, read by a readout of the control cells of its stratum |
| S1 | the system as it forecast batch 2 | the code of commit `903532a` and the steps of batch 2's `forecast.py`: one type for all recorded rows, one size per class of edge, the spread with no edge fitted over all pairs (`scripts/s1.py`) |
| S2 | version 5 | the code of commit `d777070`: `_v5.py` |

For a readout read before, the recorded rows of a stratum are all its units read so far. For a readout never read
they are the teachers alone. The library under interferon gives each system what a stratum cannot tell from its own
rows, as before.

## 6. After the forecasts are committed

`pack.py --stage B`, then `scripts/read.py`, once (it refuses if its output exists):

1. Every pair: the measured change and its error; inside or outside the interval of S1 and of S2; the errors of the
   three systems.
2. The flags of the library's law on every readout are proposed as edges (`events.propose_edges`) and asked of the
   rows of stage B, forecast without them (`events.test_edges`, rule G2).
3. The account (`account.account`) of every pair whose change is at least four errors, at the input view and the
   state, and in the co-culture and THP-1 also at the typed state and its control; `others`: the teachers of the
   stratum. Then `account.plan` and `outer.open_parts`.

Then, by rule and by nothing else: the rows of the batch are records; the admitted proposals are edges with the share
they showed. That system is S3, and it is what a final audit would be forecast with.

## 7. Hypotheses

A pair is *inside* when its measured change lies in the 95% interval of the forecast. Error: the root mean square of
measured change minus forecast. A pair with no forecast of its size (an edge of a class with too few recorded pairs)
is counted apart.

**H1. The intervals of version 5 hold.** Of the pairs with a forecast, at least 90% are inside, and the mean of their
squared *z* is between 0.5 and 2. *Refuted* below 85% inside or with the mean outside 0.33 to 3; between, unresolved.

**H2. The loop corrects.** The error of S2 is smaller than that of S1 and than that of S0, over all pairs, and over
the pairs on which the systems can differ (a pair with an edge, or an event of the declared type). *Refuted* if S2
is not the smallest on the second set.

**H3. An edge has the size of its pair (lever G).** The pairs with a declared edge on a scale that measures, for
which the library holds the pair's own share (at least four such pairs are needed): the error of S2 is smaller than
that of S1 and of S0, and at least three quarters of them are inside the interval of S2.

**H4. The records of a type forecast its other events (lever E).** The forecast pathway events on the new readouts,
pairs with no edge: the error of S2 is smaller than that of S0, and at least 80% are inside.

**H5. They move no event of another type.** Over the pairs of events with no type and no edge on the new readouts,
the error of S2 is no larger than that of S1. (S2 equals S0 there by construction.)

**H6. A proposed edge holds elsewhere.** Of the proposals that have a row in stage B, at least half are admitted. At
least four tests are needed.

**H7. The state kept by type (lever M).** In THP-1 and the co-culture, over the pairs of point 3 that have no edge:
(a) for no pair is the interval of `carried(T) − carried(Z)` entirely below −0.10; (b) in every pair where the state
carries a step less than the input, the typed state brings back at least half of it. If no pair loses a step, (b) is
*not tested*.

**The gate.** The rules stand if H1 and H2 hold and no error of the method is found that changes a verdict. If H1 or
H2 fails, or a rule has to change, the batch is development and another batch is read before a final audit.

What the author expects, written before stage A: H1 holds 0.55, H2 0.75, H3 0.6, H4 0.5, H5 0.8, H6 0.4; H7 (a)
holds 0.8, (b) not tested 0.6. The gate is passed: 0.45.

## 8. What this can and cannot show

- Most of the 698 pairs are expected to show no change: H1 is mostly about intervals around zero. The pairs with a
  change are listed one by one.
- The sizes were seen in the screen. This tests the algorithm, as batch 2 did.
- Splitting a type into teachers and forecast events halves its rows: H4 is asked with three teachers and with two.
- The spread with no edge on the new readouts is the library's, from another condition; with no signal given that
  was a little tight in batch 2.
- An edge at the background of its panel has no forecast from a scale; whether the new readouts stand above their
  background is not known yet.
- Seven units and one study. Nothing here says anything about donors, or about a second model.

## 9. Order

1. This document, `roles.json`, `pack.py`: committed and pushed.
2. Stage A on the owner's machine. `forecast.py`. The forecasts committed and pushed.
3. Stage B. `read.py`, `score.py`. The write-up, the protocol's text, the app.
