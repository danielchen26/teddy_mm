# Batch 4: one turn of the loop under version 6

27 November 2026. The owner approved on 27 November: freeze version 6, open the 30 units, forecasts first, one
reading. This document is committed before any value of the batch is opened; the forecasts are committed before the
values they forecast. Nothing is trained. No final unit and no reserve donor is touched; the final audit is not part
of this run and waits for the owner's word.

## 1. The question

Batch 3 was one turn of the loop under frozen rules, and it did not pass: the intervals of version 5 did not hold
(82% of the pairs inside where 90% was asked). Does the system of version 6, **with no rule changed**, forecast
values never read with intervals that hold?

**Version 6 is the code of commit `b4f6196`**: the package `anm_diag` and
[`runs/v6_candidate_20261126/scripts/_v6.py`](../v6_candidate_20261126/scripts/_v6.py), with the declarations of
`runs/freeze_v5_20261126/scripts/_v5.py`. Before it was frozen it passed a planted twin with the real measurement
(94% and 95% of the pairs inside where version 5 has 85% and 91%). If the batch forces a rule to change, it counts
as development, as batch 3 did.

## 2. What is opened

Thirty units never read, 15 in the co-culture and 15 with no signal given, on all 20 readouts of the melanoma panel:
600 pairs. They were chosen by a written rule before the owner approved
([`../v6_candidate_20261126/scripts/roles4.py`](../v6_candidate_20261126/scripts/roles4.py),
[`results/roles.json`](results/roles.json)): six because an item on record asks for them (CXCR4 in both conditions
and KDR in the co-culture, genes of proteins of the panel; DAG1 in the co-culture and TXNDC17 in both, sites of
edges proposed with no row to test them yet), 24 by lot. After the batch 188 units are unread and not final.

One stage: every stratum holds recorded units of no type on every readout (23 units in the co-culture, 17 with no
signal, from batches 1 to 3), so the forecasts stand on records alone (rule 6b). The states of the new cells are
computed from their RNA before the forecasts; their antibody columns are exported once, after the forecasts are
committed (`scripts/pack.py` refuses without them).

**What was known before.** The screen of 13 November read the mean shift of every antibody for every unit from the
protein files. The forecasts are the output of code that takes nothing from that screen.

## 3. Three systems, forecast before the values are opened

`scripts/forecast.py` writes `results/forecasts.json`, once.

| | what it is |
|---|---|
| S0 | the model alone: the state of the unit's cells, read by a readout of the control cells of its stratum |
| S2 | version 5: `_v5.py` |
| S3 | version 6: `_v6.py` |

Each stratum has one law over its 20 readouts; its recorded rows are all its units read so far. The graph: the
declared edges, and the edges admitted in run 43 and in batch 3. A proposed edge enters no forecast.

## 4. The gate, written first

A pair is *inside* when its measured change lies in the 95% interval of S3.

**H1. The intervals of version 6 hold.** At least 90% of the 600 pairs are inside, and the median of squared *z* is
between 0.23 and 0.91 (the spread of the bulk within a factor of two in variance). *Refuted* below 85% inside or
with the median outside 0.15 to 1.4; between, unresolved. Reported apart: the 480 pairs of the units chosen by lot,
and each stratum.

**The correction is not asked of this batch.** Its units are not chosen for what the graph holds: for nearly every
pair the graph has no edge and no type, so the three systems give the same forecast and differ only in their
intervals. That the loop corrects was shown on batches 2 and 3 and is asked again in the final audit. The errors of
the three systems are reported.

**The gate.** H1 holds, the turn is completed by rule (Section 5), and no fault of the method is found that changes
a verdict. If not, the batch is development.

What the author expects, written before any state of the batch is computed: H1 holds 0.7, unresolved 0.2, refuted
0.1. The three systems' errors agree to two digits: 0.9. Of the three tests of proposed edges at least one admits
the edge: 0.5.

## 5. After the forecasts are committed

`scripts/pack.py` on the owner's machine, then `scripts/read.py`, each part once:

1. Every pair: the measured change and its error; inside or outside the intervals of S2 and S3; the errors of the
   three systems.
2. The edges proposed earlier with no row yet (DAG1 on PD-L1; TXNDC17 on CD49f, in both conditions), asked of the
   rows of the batch, which were forecast without them (`events.test_edges`, rule G2).
3. The account of every pair whose change is at least four errors, as in batch 3.

Then, by rule and by nothing else: the rows are records; each stratum's law is fitted again with them, and what it
flags is proposed as an edge and asked of the library's row of the same pair, which was read under interferon and
took no part in the flag.

## 6. What this cannot show

- The units that remove a large protein are used up among the unread units that are not final. The batch tests the
  width of the intervals on mostly plain units; the arithmetic of the readout meets data never read in the final
  audit.
- Recorded units were chosen for what they might show, the units of the batch mostly by lot. If the recorded pairs
  stand further out than pairs drawn by lot, the intervals come out wide; the median of squared *z* in H1 is there
  to see that.
- One study, one model, knockouts.

## 7. Order

1. This document, `results/roles.json`, `scripts/pack.py`, the protocol's line that freezes version 6: committed and pushed.
2. The list of cells and their states, on the owner's machine (RNA only). `forecast.py`. The forecasts committed and pushed.
3. `pack.py`. `read.py`, `score.py`, written and committed before the reading. The write-up, the app.
