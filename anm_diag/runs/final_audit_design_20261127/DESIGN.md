# The final audit: a design, waiting for the owner's word

27 November 2026. **A draft. Nothing was opened: no cell of a final unit was listed, no state computed, no value
read.** Names only, from the census of 22 November. It becomes the design when the owner says go; it is then
committed as it stands, before anything is opened, and the forecasts before the values.

## 1. Why now, and why it needs a word

By the protocol the final audit may be run once a batch read with no rule changed passes its gate. Batch 4 did
(27 November: 581 of 600 pairs inside, median of squared *z* 0.48). The final units can be read once. After that
there is no data in this study on which nothing was fitted. That cannot be undone, so it is the owner's decision.

## 2. What would be opened

| stratum | final units | among them |
|---|---|---|
| melanoma, co-culture | 64 | B2M, HLA-A (genes of large proteins of the panel); SSR2, NPC1, ACSL3 |
| melanoma, no signal given | 63 | CD58, CD59, HLA-C, KDR (genes of proteins of the panel); JAK1, STAT1; ITGA3, TAPBP, STAT3, DAG1 (sites of edges on record) |
| THP-1 | 6 | CD86, PDCD1LG2 (genes of two of its four proteins) |
| **all** | **133** | chosen by lot on 22 November (the first third of each stratum's unread units, seed 20261122) |

The melanoma units on 20 readouts and the THP-1 units on 4: about 2,560 pairs. The final control cells (one unread
control cell in three, by a hash of its name) are the reference.

**Not part of this audit: the 73 reserve donors.** The method for donors has had one development round, on eight
donors whose donor and batch cannot be told apart, read in retrospect. It has not passed a development gate. Final
data is spent only on a claim whose method has passed one, so the reserve donors stay closed.

## 3. Under which rules

Version 6 as frozen at `b4f6196`, with the records and the graph as batch 4 left them by rule. No rule changes
between batch 4 and the audit. Types by rule (run 48) is not in: it proposes the same one type that is declared, so
the forecasts would be the same.

Three systems forecast every pair before any value is opened: S0, the model alone; S2, version 5; S3, version 6.

## 4. What is asked, written first

This is where the final audit differs from batch 4: it holds the units batch 4 could not have.

| | hypothesis | holds if |
|---|---|---|
| **F1** | **The intervals hold.** | at least 90% of all pairs inside the interval of S3, and the median of squared *z* between 0.23 and 0.91 |
| **F2** | **The loop corrects.** | on the pairs where the systems can differ (an edge, a type, or a part known by the readout's arithmetic), the error of S3 is below that of S0 and of S2 |
| F3 | An edge has the size of its pair. | own-gene pairs: error of S3 below S0's, and three in four inside |
| F4 | The arithmetic of the readout holds on data never read (rule 6a). | the pairs with no edge of the units that remove a protein of the panel: at least 90% inside S3 |
| F5 | An edge proposed elsewhere holds here. | of the edges on record that have a row here, at least half are admitted |

**The audit is passed if F1 and F2 hold.** F3 to F5 are reported with their verdicts and change nothing in that.
Refuted and unresolved are defined as in batches 3 and 4. A verdict is a verdict: nothing is fitted again, no rule
is changed, and the result is written up as it comes.

What can be claimed if it passes: on one study and one model, for knockouts, the event graph forecasts protein
changes never read with intervals that hold, and reduces the model's error where it says it can. What cannot:
another study, another model, donors, disease, time.

## 5. Cost and order

1. This document and the list of units, committed. (The owner's word comes before this step.)
2. The list of cells and their states from RNA, on the owner's machine: about 20,000 cells, about a quarter of an
   hour. The forecasts, committed.
3. The scripts that read and score, rehearsed on made-up values, committed. The export, once.
4. The reading, once: every pair, the edge tests, the accounts of the pairs with a change (one to two hours).
5. The write-up, whatever it shows.

The scripts are those of batch 4 with a third stratum.

## 6. Two things the owner may want decided first

- **Now, or after the training experiment** ([`../learned_state_design_20261127`](../learned_state_design_20261127/))?
  They do not depend on each other. The audit tests the loop around the released model. A trained state would be
  a new version and would need unread units of its own (188 are left that are not final).
- **The six THP-1 units.** Their stratum has few records, so some of their intervals would be marked as assumed. They
  hold the only own-gene units of that screen. The draft includes them and reports them apart.
