# The final audit

27 November 2026. Run 49. **The owner gave the word on 27 November.** This document is committed before any cell of a
final unit is listed, before any state of one is computed and before any value is read. The forecasts are committed
before the values they forecast, and the scripts that read and score before the reading. Nothing is trained.

The final units can be read once. A verdict is a verdict: nothing is fitted again, no rule is changed, and the
result is written up as it comes.

## 1. What is opened

| stratum | final units | read against |
|---|---|---|
| melanoma, co-culture | 64 | the 621 final control cells of the stratum |
| melanoma, no signal given | 63 | its 1,333 final control cells |
| THP-1 | 6 | its control cells (every one was read before the census; it has no final control cell) |

133 units, chosen by lot on 22 November (the first third of each stratum's unread units, seed 20261122;
[`results/roles.json`](results/roles.json)). The melanoma units on 20 readouts and the THP-1 units on 4: 2,564 pairs.
A final control cell is an unread control cell whose name hashes to a multiple of three; no run has read one.

Not opened: the 73 reserve donors. The method for donors has had one development round, on eight donors whose donor
and batch cannot be told apart, read in retrospect. Final data is spent on a claim whose method passed a
development test, so they stay closed.

## 2. Under which rules

Version 6 as frozen at `b4f6196`, with the records and the graph as batch 4 left them by rule: every unit of batches
1 to 4 is a record; the graph holds the declared edges and, per stratum, the edges admitted so far. No rule changed
between batch 4 and this audit. Types by rule (run 48) is not in.

Three systems forecast every pair before any value is opened (`scripts/forecast.py`):

| | |
|---|---|
| S0 | the model alone: the state of the unit's cells, read by a readout of the control cells of its stratum |
| S2 | version 5 |
| S3 | version 6 |

THP-1 is forecast as the frozen system has always forecast it: one law over its four readouts, the one declared
type, and no edge. Its panel has no isotype antibodies, so no scale is defined there and no edge has ever been
forecast in it. Two of its final units remove a protein of its own panel (CD86, PDCD1LG2); those two pairs are
forecast with no edge and are expected to stand outside. They count in F1 as any pair does.

**The reference.** A final unit is read against final control cells, which the laws never saw. The forecast is made
against the same reference: the law reads how far the unit's cells stand from the reference cells, and the spread of
the forecast carries the error of a reference of that size (`_fa.relative`, `_fa.against`). The change of a pair is
counted in spreads of the readout among the control cells the law was fitted on, fixed before the reading.

## 3. Before anything is opened

`scripts/preflight.py`, on cells read before ([`results/preflight.json`](results/preflight.json)):

1. the code that forecasts the audit gives the committed forecasts of batch 4 again, number for number;
2. batch 4 again with one control cell in three taken out of everything and used as the reference;
3. six THP-1 units read in batch 2, forecast again from the rest by the law over four readouts.

It finds faults of the code. It is no evidence.

## 4. What is asked

| | hypothesis | holds if |
|---|---|---|
| **F1** | **The intervals hold.** | at least 90% of all pairs inside the 95% interval of S3, and the median of squared *z* between 0.23 and 0.91. *Refuted* below 85% inside or with the median outside 0.15 to 1.4. |
| **F2** | **The loop corrects.** | on the pairs where the systems can differ, the squared error of S3 is below that of S0, and the fall stands beyond chance (below); and S3 is not worse than S2. *Refuted* if S3's error is not below S0's, or if S3 is worse than S2 beyond chance. |
| F3 | An edge has the size of its pair. | the pairs with a declared edge whose scale measures: error of S3 below S0's, and three in four inside S3 |
| F4 | The arithmetic of the readout holds on data never read (rule 6a). | the pairs with no edge of the melanoma units that have such an edge: at least 90% inside S3 |
| F5 | An edge proposed elsewhere holds here. | of the edges the library flagged that have a row here in a stratum where they were not asked before, at least half are admitted (by test and by proposal). Not tested with fewer than four tests. |

*The pairs where the systems can differ*: every pair of a unit that has an edge on any readout of its stratum, or a
declared type. *Beyond chance*: the units of that set are resampled 2,000 times (seed 20261127); the 95% interval
of the difference of the summed squared errors lies above zero (S0 minus S3), and does not lie below zero (S2 minus
S3). Between holds and refuted a hypothesis is unresolved. The sets are written into `forecasts.json` before the
values are opened.

**The audit is passed if F1 and F2 hold. It is not passed if either is refuted.** F3 to F5 are reported with
their verdicts and change nothing in that.

**Changed from the draft of the same day** ([`../final_audit_design_20261127`](../final_audit_design_20261127/)),
before anything was opened, and for a reason known beforehand. The draft asked in F2 that the error of S3 be below
that of S2 as well. Versions 5 and 6 give the same value for most pairs, and in batch 3 a comparison of two systems
that are nearly the same failed by 0.0004 with nothing behind it. So S3 is asked to be not worse than S2 beyond
chance, and better than S0 beyond chance; where versions 5 and 6 differ in substance is asked in F4. The draft also
did not say how THP-1's own-gene pairs are forecast, and did not say what a final unit is read against; both are
said above.

What the author expects, written before any cell is listed: F1 holds 0.8, unresolved 0.12, refuted 0.08. F2 holds
0.7. The audit is passed 0.6. F3 holds 0.5; F4 holds 0.65; F5 holds 0.65.

## 5. What can be claimed

If it passes: in one study and for one model, for knockouts, the event graph forecasts protein changes never read
with intervals that hold, and reduces the model's error where it says it can. If it does not: what failed, by how
much, and on which pairs.

Either way, not: another study, another model, donors, disease, time, or that the model itself was corrected.

## 6. Order

1. This document, `results/roles.json`, the scripts of the forecasts and `results/preflight.json`: committed and pushed.
2. On the owner's machine: the list of cells (`scripts/cells.py`; names only) and their states from RNA.
   `scripts/forecast.py`. The forecasts committed and pushed.
3. `scripts/read.py` and `scripts/score.py`, rehearsed on made-up values, committed and pushed. Then
   `scripts/pack.py` exports the antibody columns of the final cells, once.
4. The reading, each part once: every pair and the edge tests; then the account of every pair whose change is at
   least four errors.
5. The write-up.
