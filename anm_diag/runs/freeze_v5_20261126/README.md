# The rules are frozen: version 5 of the protocol

26 November 2026. Nothing is trained. No unit is opened: every value here was read before.

## Why

Since run 41 every turn of the loop also changed a rule; the afternoon of 25 November changed three. So far it could
not be told whether the loop corrects itself or its author corrects the loop. Closing the loop needs one full turn in
which the rules do not move. This run finishes what was still open, fixes the rules, and tags them
(`anm-protocol-v5`). Batch 3 ([`runs/batch3_frozen_20261126`](../batch3_frozen_20261126/)) is read under them.

## What was finished (rules 5a to 5f of [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md))

| | what was open | now |
|---|---|---|
| 5a | The spread of a direct part with no edge was fitted over all pairs, and a few pairs with a large part of their own widened it. | The spread at which the recorded pairs with no edge and no type, each left out, are as wide as they should be. |
| 5b | Once that spread is small, a row that stands far out taught the readout as much as any other. | Such a pair is an edge the graph does not hold yet; by how likely that is, its row teaches less. |
| 5c | A scale that measures nothing (a readout at the background of its panel) was found after the reading. | Said beforehand from control cells: below one spread, the edge has no scale and is of another class. |
| 5d | The size of an edge. | The pair's own share where the library holds one, with the spread by which shares of pairs measured twice moved; else its class's. |
| 5e | A donor could be placed by the account, and not forecast. | `Law(control=None)`: a unit is forecast from its state and the rows of the other units. |
| 5f | What a batch has to show. | Three systems forecast first; the frozen one calibrated and with the smallest error; then records and proposals by rule only. |

The five older proposals for rules of judgment are not part of version 5; they wait for the owner.

## Checks

**Systems with a known answer** (`tests/`, 203 tests pass, 1 skipped):

- a planted library: the pairs with no direct part are calibrated (inside their interval 95% to 98% of the time, mean
  squared *z* 0.7 to 1.0), and the planted pairs with a part of their own stand far outside and are flagged;
- a planted study of units, each unit left out: units that differ by a factor inside units are forecast, ordinary
  units stay inside their interval, units whose difference is not in the state stand far outside and do not bend the
  forecasts of the others.

**Cells read before** ([`scripts/check.py`](scripts/check.py), [`results/check.json`](results/check.json);
retrospective, nothing here is a forecast). The system exactly as batch 3 uses it ([`scripts/_v5.py`](scripts/_v5.py)):

| batch 2 again | pairs | inside 95% | mean squared *z* | error: no change, the model alone, the system |
|---|---|---|---|---|
| as committed before it was read | 79 | 77 | 0.30 | 0.314, 0.292, 0.231 |
| version 5, all | 79 | 73 | 1.16 | 0.314, 0.292, **0.099** |
| version 5, co-culture | 30 | 28 | 0.89 | 0.450, 0.425, 0.108 |
| version 5, no signal | 39 | 35 | 1.58 | 0.204, 0.186, 0.105 |
| version 5, THP-1 | 10 | 10 | 0.35 | 0.081, 0.024, 0.024 |

The six own-gene pairs on a scale that measures are forecast at their own share under interferon: −1.90 for −2.12,
−0.89 for −0.83, −0.94 for −0.63, −0.78 for −0.73, −0.31 for −0.33, −0.62 for −0.59. The two on HLA-E are of the class
"at the background", forecast in spreads from the two recorded pairs of that class: one inside, one outside. The
library under interferon now flags six pairs (the five of run 43 and DAG1 on PD-L1).

This table was reached on batch 2 itself, after it was read and after three rules were changed on it. It shows that
the system is not broken; it does not test it.

**Eight donors, each left out** (`Law(control=None)`; 26 proteins, 208 pairs): 192 inside (92%), mean squared *z* 1.3;
error 0.146 against 0.198 for the readout inside units alone and 0.237 for no change. On the 26 pairs run 32 admitted
for their size, 20 inside, error 0.28 against 0.35 and 0.47: a new donor is forecast better than by nothing, and not
well where the account says its difference is not in the input.

## What stays as it is, known

- With no signal given the intervals are a little tight (mean squared *z* 1.6): the spread with no edge is the
  library's, from another condition.
- A pair never measured has only its class's share, which is weak. A class "at the background" has two recorded pairs.
- Types are declared for one pathway.

## Files

`scripts/_v5.py` (the declarations of the two screens and the forecasting system), `scripts/check.py`, `results/check.json`.
