# The candidate for version 6, checked before it meets data never read (run 46)

26 and 27 November 2026. Nothing is opened, nothing is trained, and nothing here is a test: the candidate was made
from batch 3. It is not frozen. It becomes version 6 at the commit made before a fourth batch is forecast.

## Why this run

Batch 3 was read under frozen rules ([`runs/batch3_frozen_20261126`](../batch3_frozen_20261126/)). The correction
held; the intervals did not, and the gate was not passed. Three things lie behind that:

1. **The uncertainty was assumed or borrowed, not measured where it is used.** The arithmetic of the readout was not
   on the graph. The spread of a pair with no edge came from another condition. And the law assumes that some pairs
   hold a direct part of their own, but its interval was 1.96 spreads of the others: narrower than its own assumption.
2. **Every batch so far forced a rule to change.** Data never read was being used to find faults of the method,
   because the planted system of run 43 does not measure a readout as the real one is measured.
3. **The author chose the units of every batch.**

This run answers the first two before any new unit is opened, and prepares a fourth batch whose units are chosen by
a written rule.

## What changes against version 5

| | change | where |
|---|---|---|
| 6a | **The arithmetic of a readout is on the graph.** A readout that is one column minus the mean of *K* columns moves when another column moves. For every event with an edge, the direct part of the edge's readout (measured where recorded, else as the edge forecasts it) gives what every other readout of the set moves by, with no constant. A recorded row is freed of that part before it teaches; a forecast carries it. | `events.derived`, `Law(known=)`, `Law.new(known=)` |
| 6b | **A batch has recorded pairs of no type in every stratum and on every readout it forecasts**, and the spread with no edge comes from them. A design rule; rule 5a itself is unchanged. | the design of a batch |
| 6c | **The gate asks for the share inside and the median of squared *z*.** A few missing edges move the mean more than all other pairs together. | the design of a batch |
| 6d | **The width of the interval of a pair with no edge is measured.** Each recorded pair with no edge is left out and forecast as a new event is; how far it then stands, in spreads of its forecast, is its score. The interval is as wide as the score that 5% of them exceed. A new pair that is exchangeable with them is outside with probability 5% at most, whatever the shape of their distribution. A kind (of a type, of no type) needs 40 pairs from 8 events; else the library's pairs, else the stratum's of any kind; else the law's own mixture, and the row says that its width is assumed. | `Law.conformal`, `Law.new(width=)` |
| 6e | **How far a share carries is known from few pairs**, so its spread is widened by Student's *t* for their number. | `scripts/_v6.py` |

Version 5 is untouched: every addition is an argument that is off unless given, and `_v5.py` gives what it gave.

## A planted twin with the real measurement

`scripts/twin.py`, [`results/twin.json`](results/twin.json). The cells are real control cells with their states and
their 24 antibody columns, all opened before; so the readouts are centred log-ratios with the real noise and the
real background. The events are made up and planted in the log counts, so the answer is known: own-gene units that
remove a share of a protein, units of one type that move the state as a real pathway knockout does, a small direct
part for every pair, wider in the stratum than in the library, and one pair in 25 with an edge the graph does not
hold. Eight worlds, two layouts.

| layout | pairs | inside, version 5 | inside, candidate | worlds at 90% or more: version 5, candidate |
|---|---|---|---|---|
| as in batch 3 (the stratum's only recorded units are three of the type) | 5,120 | 84.6% | **94.3%** | 0 of 8, 7 of 8 |
| with recorded units of no type | 3,200 | 90.6% | **95.1%** | 6 of 8, 8 of 8 |

**Version 5 fails on the twin as it failed on batch 3** (84.6% there, 82.4% on the batch), and in the same place:
of the pairs with no edge of an own-gene unit 69% are inside. With the candidate, 98%. The pairs with an edge the
graph does not hold stay outside under both, as they should: they are what the loop proposes as edges.

| kind of pair, layout as in batch 3 | pairs | inside, version 5 | inside, candidate |
|---|---|---|---|
| no edge, of an own-gene unit | 1,459 | 69.0% | 98.4% |
| no edge, of a plain unit | 3,067 | 96.0% | 98.0% |
| no edge, of a unit of the type | 307 | 96.4% | 98.0% |
| with its edge | 80 | 100% | 100% |
| an edge the graph does not hold | 207 | 2.4% | 3.4% |

The median of squared *z* of the candidate is 0.30 to 0.69 over the sixteen worlds (0.45 if exact); the mean is 4.6
and 6.4, moved by the planted missing edges alone. From now on a change of the method is run on the twin before it
meets a unit never read.

## Batches 3 and 2 again

`scripts/replay.py`, `scripts/replay2.py`. Development and retrospect.

| | pairs | inside: version 5, candidate | median squared *z* | error |
|---|---|---|---|---|
| batch 3, from stage A alone, as it was forecast | 698 | 575, **629** (90.1%) | 0.81, 0.66 | 0.138, 0.127 |
| batch 3, with eight recorded units of no type per stratum | 402 | 315, **385** (95.8%) | 0.93, 0.56 | 0.146, 0.134 |
| batch 2 | 79 | 73, 76 | 0.17, 0.16 | 0.099, 0.095 |

With recorded units of no type the spread with no edge is the stratum's own (0.041 and 0.040, where the library
gave 0.015), and the width comes from the stratum's own 135 pairs. Of the 17 pairs then outside, the furthest are
ITGA3 on CD49f (an edge the batch admitted), IFNGR1 on CD119 (the share that did not carry) and SSR2 on CD9 (admitted).

## What the candidate does not settle

- THP-1 has too few recorded pairs of its own kind: its widths are the law's mixture, marked as assumed (43 of 48
  and 21 of 24 inside).
- A pair far outside is an edge the graph lacks, and stays outside until it is proposed and admitted.
- The units that remove a large protein, where the arithmetic matters most, are no longer among the unread units
  that are not final. A fourth batch tests mostly the width; the arithmetic meets data never read in the final
  audit, which holds five own-gene units.
- The twin plants what batch 3 showed. A fault of another kind would not be in it.

## A fourth batch, prepared and not opened

[`BATCH4_DESIGN.md`](BATCH4_DESIGN.md), `scripts/roles4.py`, [`results/batch4_roles.json`](results/batch4_roles.json):
30 units of the 218 that are unread and not final, chosen by a written rule (6 by an item on record, 24 by lot).
It is the owner's decision.

## Files

| | |
|---|---|
| `scripts/_v6.py` | the candidate: what changes against `_v5.py` |
| `scripts/twin.py`, `results/twin.json` | the planted twin |
| `scripts/replay.py`, `results/replay.json` | batch 3 again |
| `scripts/replay2.py`, `results/replay_batch2.json` | batch 2 again |
| `scripts/roles4.py`, `results/batch4_roles.json`, `BATCH4_DESIGN.md` | a fourth batch, not opened |
| `../../tests/test_version6.py` | the arithmetic, the known part, the width |
