# Batch 3: one turn of the loop with the rules frozen (run 45)

26 November 2026. Design: [`DESIGN.md`](DESIGN.md), committed before any value of the batch was opened. Nothing is
trained. No final unit and no reserve donor is touched. After the batch 218 units are unread and not final.

## The result

1. **The correction holds on values never read, with no rule changed.** Over the 698 pairs the error falls from 0.247
   (the model alone) to 0.174 (the system as it forecast batch 2) to 0.138 (version 5). On the 97 pairs where the
   systems can differ: 0.589, 0.354, 0.209. H2 holds.
2. **An edge has the size of its pair, and proposed edges hold elsewhere.** H3 holds: on the eight own-gene pairs the
   error is 0.65, against 1.19 and 1.98; six are inside their interval. H6 holds: of ten tests of edges proposed
   before the batch, eight are admitted (six of seven pairs).
3. **The intervals of version 5 do not hold.** 575 of 698 pairs are inside (82%); the mean of squared *z* is 3.08. H1
   is refuted. The older system's wider intervals would have met the same criterion (631 inside, 1.69).
4. **The gate is not passed.** As written before the batch, it counts as development, and the final audit is not
   ready. Another batch is read first.

| | hypothesis | expected before stage A | outcome |
|---|---|---|---|
| H1 | the intervals of version 5 hold: at least 90% inside, mean squared *z* in 0.5 to 2 | 0.55 | **refuted**: 82.4%, 3.08 |
| H2 | the error of version 5 is the smallest, over all pairs and where the systems can differ | 0.75 | **holds**: 0.138 < 0.174 < 0.247; 0.209 < 0.354 < 0.589 |
| H3 | an edge has the size of its pair (8 pairs) | 0.6 | **holds**: 0.65 < 1.19 < 1.98; 6 of 8 inside |
| H4 | the records of a type forecast its other events (39 pairs) | 0.5 | **does not hold**: the error halves (0.112 against 0.218), 30 of 39 inside where 32 were asked |
| H5 | they move no event of another type (580 pairs) | 0.8 | **does not hold** as written: 0.12234 against 0.12196 for the older system |
| H6 | a proposed edge holds elsewhere | 0.4 | **holds**: 8 of 10 tests, 6 of 7 pairs |
| H7 | the state kept by type (36 pairs with a change and no edge) | 0.8 / 0.6 | (a) **holds**: in no pair the typed state carries less. (b) **does not hold**: in four pairs the state carries a step less than the input, and the typed state brings back 15% to 30% of it |
| gate | H1 and H2 | 0.45 | **not passed** |

H5 fails by 0.0004: the older system let every recorded row teach every event, which moved the events of no type by
0.005 on average and, on this batch, not for the worse. The rule decides, and it says "does not hold".

## The order of events

| time (EST, 26 November) | what | commit |
|---|---|---|
| 01:14 | the design, the units, the export script | `f1916d2` |
| 01:14 | stage A exported on the owner's machine: control cells, the library under interferon, the teachers | |
| 01:25 | the forecasts of three systems for 698 pairs, the sets the hypotheses are scored on, nine edges to test | `9697c3d` |
| 01:25 | stage B exported (it refuses without the committed forecasts) | |
| 01:41 | `read.py` and `score.py`, rehearsed on a made-up stage-B file built from control cells | `7471cec` |
| 01:42 | the reading, once: the pairs, then the accounts (19 minutes) | the commit of the results |

One fault was found in `forecast.py` after its first run and before stage B: it proposed again an edge that run 43
had refused. The list of edges to test was corrected and the script run again; the 698 forecasts are identical.
`read.py` has two parts, the pairs and the accounts, so that the hypotheses on the pairs could be scored while the
accounts were read; each part ran once.

## By stratum

| | pairs | error: no change, model alone, older system, version 5 | inside: older system, version 5 | mean squared *z*: older, version 5 |
|---|---|---|---|---|
| co-culture | 355 | 0.317, 0.287, 0.204, **0.158** | 313, 274 | 2.07, 4.28 |
| no signal | 295 | 0.213, 0.209, 0.144, **0.117** | 274, 258 | 1.36, 1.84 |
| THP-1 | 48 | 0.138, 0.090, 0.073, 0.076 | 44, 43 | 0.92, 1.83 |
| pairs with a change of at least four errors | 48 | 0.938, 0.869, 0.556, **0.380** | 25, 21 | |

The eight own-gene pairs (change; forecast of version 5 and its interval):

| pair | co-culture | no signal |
|---|---|---|
| CD44 on CD44 | −2.38; −1.39 [−2.80, 0.01] | −1.62; −1.16 [−2.31, −0.02] |
| CD47 on CD47 | −3.32; −2.58 [−3.32, −1.83] (outside by 0.002) | −2.31; −2.30 [−2.99, −1.62] |
| CD58 on CD58 | −1.71; −1.55 [−1.81, −1.28] | (final unit, not opened) |
| CD59 on CD59 | −1.08; −0.89 [−1.18, −0.59] | (final unit, not opened) |
| IFNGR1 on CD119 | −1.43; −0.18 [−0.52, 0.15] **outside** | −1.10; −0.89 [−1.25, −0.54] |

IFNGR1 on its own receptor in the co-culture is the one pair whose share did not carry: 33% of the scale under
interferon, 90% here.

## Why the intervals did not hold

Looked at after the reading ([`scripts/explore.py`](scripts/explore.py), [`results/explored.json`](results/explored.json)).
This is no test.

**1. The arithmetic of the readout.** A readout of the melanoma screen is the log count of one antibody minus the mean
over the 24 columns of the cell. A unit that removes a large protein therefore raises every other readout, by the
change of that protein divided by the number of other columns. No system here forecast that, and nothing has to be
fitted to forecast it.

| unit | mean miss of version 5 on its other readouts | what the arithmetic gives | pairs outside: before, after |
|---|---|---|---|
| CD274, co-culture | +0.147 | +0.161 | 11, 2 |
| CD47, co-culture | +0.156 | +0.117 | 13, 1 |
| CD58, co-culture | +0.085 | +0.082 | 8, 2 |
| CD59, co-culture | +0.056 | +0.046 | 4, 0 |
| TAPBP, co-culture | +0.095 | +0.090 | 5, 4 |
| B2M, no signal | +0.122 | +0.101 | 6, 5 |
| CD47, no signal | +0.088 | +0.086 | 4, 0 |

With that term, from what was known before stage B, 619 of the 698 pairs are inside (88.7%) and the mean of squared
*z* is 2.46. `score.py` had declared one other look at this before the reading: the same forecasts against the change
measured with the isotype antibodies as the reference. That gives 559 of the 650 melanoma pairs inside against 532,
with a slightly larger error (0.148 against 0.141): the isotype reference is noisier and moves itself for some units.

**2. A spread borrowed from somewhere else.** The spread of a pair with no edge has to come from pairs like it. For
the readouts never read no stratum had such a pair of its own, because the teachers were all of the pathway type, so
rule 5a fell back on the library under interferon (0.018). The rows of this batch give 0.048 in the co-culture, 0.050
with no signal and 0.030 in THP-1 by the same rule; for the readouts read before, 0.063 where the stratum's own
recorded rows had given 0. With those spreads 93% of the melanoma pairs with no edge and no type are inside (in-sample).

**3. Edges that are not on the graph.** After the first term 25 pairs stand more than three errors out. Of the six
furthest, five are edges the library had proposed before the batch and the batch admitted (ITGA3 on CD49f at 19.6
errors, CD44 on CD47 and on CD29, SSR2 on CD9). The sixth is IFNGR1 on CD119, above. A few pairs this far out move
the mean of squared *z* a great deal: that criterion of H1 mixes how wide the intervals are with how many edges are
missing.

## What the batch leaves, by rule

- **Records**: 58 units on every readout of their panel (872 pairs).
- **Edges**: eight admitted records of six pairs, with the share of the scale they showed here and under interferon:
  ITGA3 on CD49f (+0.33, +0.25), CD44 on CD47 (+0.10 and +0.08; +0.06), CD44 on CD29 (+0.04 and +0.03; +0.03), SSR2 on
  CD9 (+0.10, +0.10), SSR2 on CD59 (−0.06, −0.07), NONO on PD-L1 in the co-culture (+0.04, +0.08). Refused: NONO on
  PD-L1 with no signal, SSR2 on CD44. Still proposed, with no row to test them: TXNDC17 on CD49f, DAG1 on PD-L1.
- **Pairs measured twice**: 16 more. The shares of the eight own-gene pairs here against under interferon: CD44 −0.19
  and −0.16 for −0.11; CD47 −0.49 and −0.38 for −0.38; CD58 −0.80 for −0.73; CD59 −0.43 for −0.35; IFNGR1 −0.90 and
  −0.40 for −0.33. Over the 16 the spread of a carried share is 0.15, where the seven earlier pairs gave 0.06.
- **Open parts** ([`results/accounts.json`](results/accounts.json)): of the 48 changes of at least four errors, the
  largest open part is outside the input in 31 (the eight own-gene pairs among them, 89% to 101%, on the graph), lost
  before the state in 4 and held but not read in 1; 12 have only unresolved parts. `outer.open_parts` lists 99 parts:
  23 open (lever G 8, E 10, M 5), 13 on the graph, 63 at the margin. The four pairs that lose a step before the
  state are IFNGR1 and JAK1 on CD9 in the co-culture (the input carries 0.77 and 0.81 of the change, the state 0.51,
  the typed state 0.59 and 0.56), IRF1 on CD366 and BRD4 on CD86 in THP-1.

## What this does not show

- The sizes of the changes were seen in the screen of 13 November. This tests the algorithm, as batch 2 did.
- Seven new units and one study. Nothing here is about donors or a second model.
- Everything under "Why the intervals did not hold" was found on the batch itself.

## Next

Proposed, and the owner's decision, since it opens units never read:

1. Version 6, two changes: the arithmetic of a readout that is a ratio enters the graph as derived edges (no constant);
   a batch is designed so that every stratum has teachers of no type on every readout it forecasts, and the spread
   comes from them. H1 keeps its share inside; the mean of squared *z* is replaced by the median, which a few missing
   edges do not move.
2. Batch 4 under the frozen version 6: about 30 units of the 218, forecasts first, the same gate.
3. Then the final audit.

Also the owner's: the four pairs in which the state loses a step that the typed state does not bring back are what
the protocol sends to training of what comes before the state. Nothing was trained.

## Files

| | |
|---|---|
| [`DESIGN.md`](DESIGN.md) | the design, with the hypotheses and the gate |
| `scripts/roles.py`, `results/roles.json` | the seven units and the teachers |
| `scripts/pack.py` | the export on the owner's machine, in two stages |
| `scripts/_b3.py`, `scripts/s1.py`, `scripts/forecast.py`, `results/forecasts.json` | the three systems and their forecasts |
| `scripts/rehearsal.py` | a made-up stage-B file for finding faults of the code |
| `scripts/read.py`, `results/reading.json`, `results/accounts.json` | the reading |
| `scripts/score.py`, `results/outcome.json` | the hypotheses against the reading |
| `scripts/explore.py`, `results/explored.json` | explored after the reading |
