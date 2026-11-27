# Batch 4: one turn of the loop under version 6. The gate is passed

27 November 2026. Run 47. The design, the forecasts and the scripts that read and score were each committed before
what they depend on was opened ([`DESIGN.md`](DESIGN.md); commits `936426f`, `5b8f447`, `208e165`; the batch was
exported after the last). Every part was read once. Nothing was trained. No final unit and no reserve donor was
touched, and the final audit was not run: it waits for the owner's word.

## The result

Thirty units never read, 15 in the co-culture and 15 with no signal given, on the 20 readouts of the melanoma panel:
600 pairs. Version 6 as frozen at `b4f6196`, no rule changed.

| | asked beforehand | measured |
|---|---|---|
| pairs inside the 95% interval of version 6 | at least 90% | **581 of 600 (96.8%)** |
| median of squared *z* | 0.23 to 0.91 (0.455 for an interval of exactly the right width) | **0.48** |
| **H1, the gate** | | **holds** |

The turn was completed by rule (below), and no fault of the method was found. **This is the first batch read with
no rule changed that passes its gate.**

| | pairs | inside (v6) | median *z*² (v6) | inside (v5) | median *z*² (v5) |
|---|---|---|---|---|---|
| all | 600 | 581 (96.8%) | 0.48 | 584 (97.3%) | 0.36 |
| units chosen by lot | 480 | 472 (98.3%) | 0.44 | 474 (98.8%) | 0.32 |
| units chosen by an item on record | 120 | 109 (90.8%) | 0.65 | 110 (91.7%) | 0.47 |
| co-culture | 300 | 292 (97.3%) | 0.42 | 293 (97.7%) | 0.30 |
| no signal given | 300 | 289 (96.3%) | 0.52 | 291 (97.0%) | 0.43 |

Root mean square error of the forecast, in spreads of the readout: 0.097 for "no change", 0.086 for the model alone,
0.086 for version 5 and for version 6. The three systems give the same value on this batch, as expected: the graph
holds almost nothing about these units.

## How to read it

**What passed is the width of an interval on plain units.** Two of the 600 pairs changed by four errors or more.
The batch asks whether the system knows how far a pair moves when nothing known acts on it, and it does: the bulk is
right (median 0.48 where 0.455 is exact) and the tails are covered.

**Version 5 would have passed this batch too.** What failed in batch 3 had two causes. One was a spread borrowed
from another condition; in this batch every stratum has its own recorded units of no type (rule 6b), and then both
versions take the spread from them. The other was the arithmetic of the readout on units that remove a large
protein (rule 6a); this batch has no such unit, because none is left among the unread units that are not final.
So the batch tests 6b, 6c and 6d, and it does not test 6a. Where the two versions differ, version 6 is the better
one: its bulk is right (0.48), version 5's intervals are wider than they need to be (0.36), with the same share
inside.

**The intervals are a little wide for units drawn by lot** (98.3% inside). The design said beforehand that this
could happen: the recorded units were chosen for what they might show, so they stand further out than units drawn
by lot, and the width is learned from them. The median, which the gate asks for this reason, is right.

**One unit does something, and the system had said where.** TXNDC17 holds 10 of the 19 pairs outside. Its edge to
CD49f was proposed from the library under interferon before any of its rows here was read; a proposed edge enters
no forecast.

## The edges asked of the batch

| edge | proposed with *z* | *z* in the batch | | direct part as a share of the scale |
|---|---|---|---|---|
| TXNDC17 → CD49f, co-culture | −4.5 | −5.5 | **admitted** | −0.12 (under interferon: −0.10) |
| TXNDC17 → CD49f, no signal | −4.5 | −3.9 | **admitted** | −0.11 |
| DAG1 → PD-L1, co-culture | −3.9 | −1.4 | refused | |

An edge flagged in one condition holds in two others on rows never read, with the same share of its scale (−0.10,
−0.12, −0.11). The account of the two pairs: 80% and 59% of the change is not in the input (the RNA does not show
it), which is what an edge on the graph is for.

## The turn, by rule

- **Records.** 30 units, 600 rows.
- **The laws fitted again** with them: co-culture 38 recorded units, spread of a pair with no edge 0.035; no signal
  32 units, 0.046.
- **Flags → proposals → test.** The laws flag two pairs, both of units recorded earlier (TAPBP → CD9, PFN1 → CD9 in
  the co-culture). Both were proposed as edges and asked of the library's rows of the same pairs, which took no part
  in the flag: both refused.
- **Accounts.** Two pairs changed by four errors or more (both TXNDC17 → CD49f). Their six open parts each have a
  place: two on the graph (the edge now admitted), three at the margin, one open (lost before the state, 40% of the
  change in the co-culture: the lever is the model).
- **Pairs measured twice.** Three more (the two TXNDC17 pairs and CXCR4 → CD184 with no signal); their shares moved
  by 0.016 of a scale.

What the author expected, written before any state was computed: H1 holds 0.7 (it holds); the three errors agree to
two digits 0.9 (they do); at least one of three edge tests admits 0.5 (two do).

## What this does not show

- **The arithmetic of the readout on data never read** (rule 6a). It was checked on the planted twin and on batch 3
  again; its first test on unread data is the final audit, which holds the units that remove a large protein.
- **The correction.** The three systems give the same values here. That the loop reduces the error was shown on
  batches 2 and 3 and is asked again in the final audit.
- A typed event: no unit of the batch has a type.
- Anything beyond one study, one model, knockouts.

## What comes next

The gate is passed, so by the protocol the final audit may be run: 133 units and 73 reserve donors, untouched. It
is run only on the owner's word. After this batch 188 units are unread and not final.

## Files

| | |
|---|---|
| [`DESIGN.md`](DESIGN.md) | the question, the gate and the order, committed first |
| [`results/roles.json`](results/roles.json) | the 30 units and the rule that chose them |
| [`results/forecasts.json`](results/forecasts.json) | 600 pairs by three systems, committed before the export |
| [`results/reading.json`](results/reading.json) | every pair as measured; the edge tests; the records |
| [`results/by_rule.json`](results/by_rule.json) | the laws fitted again, their flags, the proposals and their test |
| [`results/accounts.json`](results/accounts.json) | the accounts of the pairs with a change; the open parts |
| [`results/outcome.json`](results/outcome.json) | the gate and everything reported beside it |
| `scripts/` | `forecast.py`, `pack.py` (the export, on the owner's machine), `rehearsal.py`, `read.py`, `score.py` |

```
python scripts/forecast.py --cells batch4_cells.tsv --states B2 B3 B4 --library LIBRARY.npz --packs batch3_pack_A.npz batch3_pack_B.npz
python scripts/pack.py --repo TEDDY_MM --cells batch4_cells.tsv --forecasts results/forecasts.json --out OUT_DIR
python scripts/read.py --part pairs|rule|accounts --cells batch4_cells.tsv --states B2 B3 B4 --packs batch3_pack_A.npz batch3_pack_B.npz batch4_pack.npz [--library LIBRARY.npz] [--counts B2 B3 B4]
python scripts/score.py
```
