# Batch 4: a design, not yet approved

27 November 2026. Nothing of this batch is opened. Opening it is the owner's decision; this document is what would
be committed, with the frozen version 6, before any of its values is read.

## The question

Batch 3 was one turn of the loop under frozen rules, and it did not pass: the intervals did not hold. Does the
system of version 6, **with no rule changed**, forecast values never read with intervals that hold?

## What would be opened

Thirty units never read, 15 in the co-culture and 15 with no signal given, on all 20 readouts: 600 pairs. They are
chosen by a written rule ([`scripts/roles4.py`](scripts/roles4.py), [`results/batch4_roles.json`](results/batch4_roles.json)):
six because an item on record asks for them (CXCR4 in both conditions and KDR in the co-culture, genes of proteins
of the panel; DAG1 and TXNDC17, sites of edges proposed with no row to test them yet), 24 by lot. After the batch
188 units are unread and not final. The 133 final units and the 73 reserve donors stay unopened.

One stage is enough: every stratum now holds recorded units of no type on every readout (23 units in the co-culture,
17 with no signal), so nothing has to be opened first for the forecasts to stand on (rule 6b).

## Three systems, forecast before the values are opened

| | |
|---|---|
| S0 | the model alone |
| S2 | version 5, the code of commit `d777070` |
| S3 | version 6, frozen at the commit of that day |

## The gate, written first

**H1. The intervals of version 6 hold.** At least 90% of the pairs are inside, and the median of squared *z* is
between 0.23 and 0.91 (the spread of the bulk within a factor of two in variance). *Refuted* below 85% inside or
with the median outside 0.15 to 1.4; between, unresolved. Reported apart: the 480 pairs of the units chosen by lot.

**The correction is not asked of this batch.** Its units are not chosen for what the graph holds: for nearly every
pair the graph has no edge and no type, so the three systems give the same forecast and differ only in their
intervals. That the loop corrects was shown on batches 2 and 3, and it is asked again in the final audit, which
holds own-gene and pathway units. The errors of the three systems are reported.

**The gate.** H1 holds, the turn is completed by rule (the rows become records, the pairs far outside become
proposals of edges), and no fault of the method is found that changes a verdict. Then the final audit is forecast
with the same frozen version. If not, the batch is development, as batch 3 was.

Also read, outside the gate: the edges on record asked of the new rows (rule G2), and the account of every pair
with a change of at least four errors.

## Before the forecasts are made

- Version 6 is frozen at one commit. The candidate has passed the planted twin (94% and 95% inside where version 5
  has 85% and 91%) and the replay of batch 3 (385 of 402 inside with recorded units of no type).
- The owner's machine computes the states of the 30 units and exports their antibody columns once, after the
  forecasts are committed.

## What it cannot show

- The units that remove a large protein are used up among the unread units that are not final. This batch tests the
  width of the intervals on mostly plain units; the arithmetic of the readout meets data never read in the final audit.
- One study, one model, knockouts.
