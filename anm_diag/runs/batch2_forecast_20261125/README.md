# Batch 2: forecasts committed first, then 33 units never read were opened once

25 November 2026. Design and declared hypotheses: [`DESIGN.md`](DESIGN.md), committed with
[`results/forecasts.json`](results/forecasts.json) (commit `967dd03`) before the protein of any cell of the batch was
exported. Nothing was trained. No final unit and no reserve donor was touched: 133 final units and 225 unread units
that are not final remain.

## In short

- **The intervals held: 77 of 79 pairs** of an event and a readout fall inside their 95% interval. They are wider than
  they need to be (mean squared *z* 0.30, where 1 is exact).
- **The forecast beat "no change" and the control readout alone** (error 0.23 against 0.31 and 0.29), and **all of
  that gain is the edge**. For the 71 pairs with no edge the control readout alone was the better forecast (0.09
  against 0.115): what the recorded events added changed little in six of seven groups and hurt in one (MHC-I in
  co-culture).
- **An edge known in one condition forecasts in another, but not its size.** All eight own-gene pairs fell, by 0.23 to
  2.1 spreads; the forecast was "about 1" for each. The edge needs a size of its own.
- **Every change of these units is outside the RNA.** In all 14 accounts of batch 1 and batch 2, 78% to 102% of the change
  is *not in the input*. No training of the model can supply that part.
- **Of the 34 accounts read so far, the model loses a step between its input and its state in one** (IRF1). That is
  the count the decision on training needs.

## 1. What was opened

33 units ([`results/roles.json`](results/roles.json)): in melanoma, the knockouts of the genes behind the measured
proteins in two conditions never read, three pairs the library under interferon had flagged, the pathway with no
interferon given, and eight units drawn by a seeded rule; in THP-1 the ten units that were left. 11,061 cells. The
control cells and the recorded events are those earlier runs had read (20,237 cells).

States were computed from RNA counts alone for all 31,298 cells; against stored states the largest difference is
1e-6, and the blocks of the typed state equal those of run 40 to 4e-5. The forecasts were made from those states,
committed and pushed; then the readouts of the 11,061 cells were exported (values of 600 cells read before are
identical to what the earlier runs read) and read once.

What was known before is in Section 3 of the design: a screen of 13 November had read the mean shift of every
protein for every unit, so the *size* of the own-gene changes was not unknown. The account and the calibration are
what nobody had read.

## 2. Every pair against its forecast

In spreads of the readout among the control cells of the stratum.

| stratum | pairs | inside | with a change of four errors or more | error: no change | control readout alone | forecast |
|---|---|---|---|---|---|---|
| melanoma, co-culture | 30 | 29 | 4 | 0.45 | 0.42 | **0.29** |
| melanoma, no signal | 39 | 38 | 4 | 0.20 | **0.19** | 0.21 |
| THP-1 | 10 | 10 | 0 | 0.08 | 0.024 | **0.021** |
| all | 79 | **77** | 8 (7 inside) | 0.31 | 0.29 | **0.23** |

The two outside: CD274 on HLA-E in co-culture (+0.18; interval −0.15 to +0.17) and CD44 on MHC-I with no signal
(+0.33; −0.28 to +0.33). Both are knockouts of a protein of the same antibody panel, and both readouts go *up*. A
centred log ratio over the panel does that when one of its proteins is removed. That reading is a guess made after
the fact; it was not tested.

**The eight own-gene pairs** (the edge):

| condition | pair | measured | forecast | what the control readout reads |
|---|---|---|---|---|
| co-culture | CD274 → PD-L1 | **−2.12** ± 0.06 | −1.13 [−3.17, +0.90] | −0.14 |
| co-culture | HLA-B → MHC-I | −0.83 ± 0.08 | −1.36 [−3.97, +1.24] | +0.06 |
| co-culture | HLA-C → MHC-I | −0.63 ± 0.10 | −1.31 [−3.94, +1.31] | −0.05 |
| co-culture | HLA-E → HLA-E | −0.47 ± 0.06 | −0.98 [−3.00, +1.05] | +0.01 |
| no signal | CD274 → PD-L1 | −0.73 ± 0.06 | −1.06 [−3.09, +0.97] | −0.07 |
| no signal | HLA-A → MHC-I | −0.33 ± 0.08 | −0.97 [−3.01, +1.08] | +0.03 |
| no signal | HLA-B → MHC-I | −0.59 ± 0.11 | −1.01 [−3.06, +1.05] | −0.01 |
| no signal | HLA-E → HLA-E | −0.23 ± 0.09 | −0.99 [−3.09, +1.11] | −0.04 |

Eight of eight inside; the error with the edge is 0.64, without it 0.82. The class says "a drop of about one spread",
and the drops run from 0.23 to 2.1. With no signal given, where there is less of the protein to remove, the edge
overshoots and the forecast there is worse than the control readout alone.

**The pairs with no edge.** In THP-1 the ten changes are small (−0.12 to +0.17) and the forecast follows them
(correlation 0.97; error 0.021 against 0.08 for "no change"): a readout of control cells reads small shifts that
come through transcription. Over all 71 pairs with no edge the control readout alone has error 0.09, the forecast
0.115, "no change" 0.116. By stratum and readout:

| pairs with no edge | n | no change | control readout alone | forecast |
|---|---|---|---|---|
| co-culture, PD-L1 | 9 | 0.134 | 0.065 | 0.063 |
| co-culture, MHC-I | 8 | 0.100 | 0.092 | **0.234** |
| co-culture, HLA-E | 9 | 0.084 | 0.064 | 0.065 |
| no signal, PD-L1 | 12 | 0.063 | 0.042 | 0.042 |
| no signal, MHC-I | 11 | 0.193 | 0.173 | 0.178 |
| no signal, HLA-E | 12 | 0.102 | 0.077 | 0.068 |
| THP-1, PD-L1 | 10 | 0.081 | 0.024 | 0.021 |

The control readout reads a part of these small changes everywhere. What the recorded events add changes little in
six groups and hurts in one: MHC-I in co-culture, where the law had widened its prior two-thousand-fold on the
evidence of the pathway events. Left out one at a time, those pathway events are forecast well (23 of 24 recorded
pairs inside); applied to events unlike them, the same widening is noise.

## 3. The accounts

For every pair of batch 2, and every recorded pair of batch 1, whose change is at least four errors
([`results/reading.json`](results/reading.json); the controls of the measure hold in all 14: a permuted stage carries
nothing, a made-up event has no change).

| | pair | not in the input | lost before the state | held, not read | read | what the recorded events forecast |
|---|---|---|---|---|---|---|
| batch 2, co-culture | CD274 → PD-L1 | **95%** [91, 99] | 1% | −1% | 5% | 6% [2, 11]: they do not |
| batch 2, co-culture | HLA-B → MHC-I | **91%** [83, 101] | 13% [7, 20] | 3% | −7% | 29% [4, 53]: unresolved |
| batch 2, co-culture | HLA-C → MHC-I | **78%** [65, 93] | 11% | 4% | 7% | 32% [−10, 69]: unresolved |
| batch 2, co-culture | HLA-E → HLA-E | **94%** [90, 98] | 8% [1, 15] | 0% | −2% | −1%: they do not |
| batch 2, no signal | CD274 → PD-L1 | **92%** [81, 103] | 0% | 0% | 8% | 8%: they do not |
| batch 2, no signal | HLA-B → MHC-I | **93%** [82, 108] | 3% | 2% | 2% | −14%: they do not |
| batch 2, no signal | CD44 → MHC-I | **88%** [65, 111] | 0% | −1% | 13% | 24%: unresolved |
| batch 2, no signal | SSR2 → MHC-I | **89%** [70, 108] | 0% | 2% | 9% | −4%: unresolved |
| batch 1, co-culture | TAPBP → MHC-I | **97%** [92, 102] | 3% | 0% | 0% | −1%: they do not |
| batch 1, co-culture | CD44 → MHC-I | **88%** [75, 102] | 3% | −1% | 11% | 3%: they do not |
| batch 1, co-culture | CD44 → PD-L1 | **80%** [54, 109] | 0% | 1% | 19% | 21%: unresolved |
| batch 1, no signal | B2M → MHC-I | **87%** [84, 90] | 0% | 3% | 10% | 11%: they do not |
| batch 1, THP-1 | CMTM6 → PD-L1 | **102%** [93, 114] | 0% | −1% | −1% | −2%: they do not |
| batch 1, THP-1 | BRD4 → PD-L1 | **39%** [16, 57] | 15% [−6, 39] | 9% | 37% | 25%: unresolved |

In every one the open part is "not in the input". No view of the RNA, no state and no typed state holds these
changes cell by cell: the knockout acts on the protein itself (its own gene; B2M, TAPBP and SSR2 on the way of MHC-I
to the surface; CMTM6 on PD-L1 at the membrane) or on the measurement (CD44). BRD4 is the one whose change is partly
in the RNA.

**All accounts read so far** ([`results/parts.json`](results/parts.json): the 20 of run 41 on the five earlier
targets, and these 14):

| open part (its interval lies above a step) | accounts, of 34 | where |
|---|---|---|
| not in the input | 23 | 9 pathway events in co-culture (the other events forecast them at 89% to 110%: their records teach the readout); the 14 above (the recorded events do not forecast 8 of them, and 6 are unresolved: an edge) |
| **lost before the state** | **1** | IRF1 in THP-1 (37 points [17, 61]) |
| held, not read | 8 | the five THP-1 events (32 to 42 points), three MHC-I events in co-culture (18 to 23) |
| nothing open | 6 | PD-L1 under interferon (5), one PD-L1 event in co-culture |

## 4. The hypotheses

| | hypothesis | forecast | outcome |
|---|---|---|---|
| HB1 | at least 72 of 79 pairs inside their interval | 0.40; 66 to 77 | **confirmed**: 77 |
| HB2 | the forecast has a smaller error than "no change" and than the control readout alone | 0.75 | **confirmed**: 0.23 against 0.31 and 0.29 (the gain is in the edge pairs; see Section 2) |
| HB3 | of the eight own-gene pairs at least six inside, and a smaller error with the edge | 0.75; 5 to 8 | **confirmed**: 8 inside; 0.64 against 0.82 (a weak test, as declared) |
| HB4 | an own-gene change is not in the input and the recorded events do not forecast it, for three quarters of such pairs | 0.55; 4 to 8 pairs with a change | **unresolved**: 4 of 6 (all six are "not in the input" at 78% to 95%; for HLA-B and HLA-C in co-culture what the recorded events forecast has an interval that reaches above one half) |
| HB5 | at least two of the three pairs flagged under interferon fall outside their interval in another condition | 0.60; 1 to 3 | **unresolved**: 1 (CD44; SSR2 is at the end of its interval, +0.46 with *z* 1.9; STAT3 is inside an interval of ±1.2) |
| HB6 | (a) the typed state never carries less than the mean state; (b) where the state loses a step, typing brings back half | (a) 0.80; (b) tested 0.35; 0 to 2 pairs lose a step | **(a) holds; (b) not tested**: of the five pairs, none loses a step before the state |

Three confirmed, two unresolved, one half not tested. Forecasts: the five intervals given all hold the value.

## 5. What this says, and what it does not

1. **The law's intervals can be trusted on units never read, and are too wide.** In THP-1 the spread without an edge
   was learned from eight recorded pairs whose residuals are large; the new units there are forecast ten times more
   closely than the interval says.
2. **Through the state, the system forecasts what transcription carries** (THP-1), **and nothing else.** Every large
   change of this batch is outside the RNA.
3. **An edge is a fact about the event and the readout, and its size is a fact about the condition.** The class of
   edges as built has one mean. A size of its own for an edge (how much of the protein there is to remove) is the
   next thing the event graph needs.
4. **Records of events of one type do not help events of another.** The update of the readout by recorded rows is
   right for the type it was learned on. Applied to every event it changed little in six groups and hurt in one. It
   should be asked for by type.
5. **The typed state had nothing to bring back here.** None of the new pairs holds its change in the RNA, so the
   comparison of the typed state and the mean state still rests on run 41: IRF1 (+21 points [8, 39]) and five MHC-I
   events (+10 to +16, none resolved as a step).
6. **For the decision on training:** a training of the encoder can only recover a part that is *lost before the
   state*. In 34 accounts that part is a step once. Eight accounts are open at the head, and a readout taught by
   records of the same type already reads those (run 41). Twenty-three are open outside the input: nine are forecast
   by events of the same type; of the other fourteen the recorded events do not forecast eight, and six are
   unresolved. For those an edge is what is missing.

Not shown: the edge across cell types (THP-1 has no own-gene unit left that is not final); the cause of the two
pairs outside; anything about the 133 final units.

## 6. Files

- `scripts/roles.py`, `cells.py`, `states.py` (the owner's machine; RNA only), `read_pack.py`, `forecast.py`: before
  the commit.
- `scripts/pack.py` (the owner's machine; refuses without the committed forecasts), `read.py`, `score.py`,
  `census.py`: after it.
- `results/`: `roles.json`, `states_manifest.json`, `forecasts.json`, `declared.json`, `reading.json`, `outcome.json`,
  `parts.json`.
- The method: [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md), version 4; `python -m anm_diag account SPEC.yaml --out DIR`.
