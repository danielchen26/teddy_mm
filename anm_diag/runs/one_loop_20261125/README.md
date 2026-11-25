# One loop for every shape of event: the account, the edges and the step after a rejection in the forms the method had

25 November 2026, afternoon. Design: [`DESIGN.md`](DESIGN.md), committed (`7c09912`) before the runs on real cells of
parts B and C. Nothing was trained. No unit was opened: every value here was read before. The 133 final units, the 73
reserve donors and the 225 unread units that are not final are untouched.

## In short

The owner asked why three things were listed as missing when the method already had a form for each. He was right:
the forms were there and version 4 had not been carried into them. They are connected now, and the connection was
checked on systems with a known answer and run on cells read before.

- **One reading for three shapes of event** (rule 4i). An event with control cells, a unit against the other units
  (a donor, a batch, a study), a group of units against another (a disease): the same three readings and four parts.
  On a planted study of units, 216 placements: **174 right, 42 undecided, none wrong**; both controls of the measure
  hold in 144 of 144.
- **An edge is a record of the event graph** (rule 4j): declared from the run's declarations, or proposed from
  measured rows and admitted or refused by rows it did not see. The library under interferon proposes five edges.
  Four of them can be asked of another condition: **three are admitted** (one of them in both conditions), one is
  refused. The fifth has no row elsewhere.
- **The law by type** (rule 4k). The recorded pathway events still forecast each other (error 0.11 against 0.64 for
  the control readout in co-culture) and no longer move an event of another type. With the error of a row counted
  once, the intervals of batch 2 are about right: mean squared *z* **1.02** (it was 0.30), 74 of 79 inside.
- **The size of an edge: the class is the wrong carrier, the pair is the right one.** The declared hypothesis (the
  class's share of a scale) **does not hold**: error 0.84, worse than one size per class (0.68). Explored afterwards:
  the *pair's own* share, measured under interferon and multiplied by the scale of the other condition, forecasts
  the ten own-gene changes with an error of **0.23**. It works for PD-L1 and MHC-I and fails for HLA-E, which stands
  at the background of the antibody panel.
- **The step after a rejection** (rule 4l): every open part of an account names its lever, its proposal and what
  judges it. For the 14 accounts of batch 1 and batch 2, 11 of the 14 parts outside the input are now *on the graph*;
  three are open (one to propose an edge, two that need more records of their type).
- Three faults of the method were found on the way and corrected (Section 6). One of them changes three readings of
  run 41 from "forecast" to "unresolved".

What is still not done is in Section 7. Nothing here is a confirmation on units never read.

## 1. What was connected

| what the method had | what version 4 did with it | now |
|---|---|---|
| events with no control cells (`events.contrast_events`, `cycle.one_units`, run 32 on eight donors) | `account.reading` asked for control cells and resampled cells | `account.reading(control=None, cluster=...)`: a unit is read against the cells of every other unit; an arm made of units is resampled by unit; the readout "with no record" is one of the variation inside a unit |
| the event graph as one table of typed records, rule S3 (`events.py`) | the edges were a table in a run's script; S3 sent an event that acts on the readout itself to the claim as a limit | `events.edge`, `own_edges`, `propose_edges`, `test_edges`, `edge_table`, `check_edges` (rules G1 to G5); an edge carries a class, a scale and the share its pair showed |
| the step after a rejection (`outer.ROUTES`, `outer.open_items`) | the routes read the verdicts of version 3; "not in the input" had no lever | `account.plan`, `outer.open_parts`, `outer.ROUTES4`, lever **G** (the graph) |

The rules are 4i to 4l of [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md).

## 2. Systems with a known answer

`tests/test_one_loop.py` (nine tests) and [`scripts/planted_units.py`](scripts/planted_units.py): the planted study of
units of `toy.py` (ordinary units and nine that each differ in a known place), every special unit read against the
other units, at two sizes and twelve seeds. The scoring rule is in the script's docstring and was written before the
run ([`results/planted_units.json`](results/planted_units.json)).

| size | right | undecided | wrong | controls of the measure |
|---|---|---|---|---|
| 40 ordinary units of 300 cells | 104 | 4 | 0 | 72 of 72 |
| 16 ordinary units of 150 cells | 70 | 38 | 0 | 72 of 72 |

At the small size the reading rarely decides for a unit whose change is read (33 of the 38 undecided): with two dozen
units resampled by unit the intervals are wide. That is the power of the reading, and it says nothing false.

The first run ([`results/planted_units_first.json`](results/planted_units_first.json)) had **two wrong** of 216:
a unit told apart completely by a marker, where a readout of the other units forecast 14 and 26 times its change and
the rule read "the other events forecast it". The rule was one-sided. It is two-sided now (Section 6).

## 3. Real cells: eight donors, each against the other seven (explored)

[`scripts/donors.py`](scripts/donors.py), [`results/donors.json`](results/donors.json). The cells of run 32 and the 26
pairs of a donor and a protein it admitted. *This was seen before the design was written; it has no hypothesis and is
counted nowhere.* A donor is also a batch of the measurement here: its change cannot be told from its batch's.

| what the account says of the 26 pairs | pairs |
|---|---|
| not in the input | 11 |
| held by the state, not read | 6 |
| held, not read, and a part not in the input | 4 |
| held, not read, and a part lost before the state (CD20, donor 8: 18%) | 1 |
| nothing open (CD19, CD11c, CD21, IgD) | 4 |

The step after a rejection: 11 parts go to the head (a readout taught by the rows of the other donors reads 23% to 90%
of those changes on cells it was not fitted on); 7 parts outside the input that the other donors do not forecast
become **proposed edges** between a donor and a protein, and stay proposed: this study has no second measurement of a
donor to test them on; 8 need more donors to say whether they are forecast; one goes to the model before the state.

Version 3 read these donors by one number per pair. The same cells now say where each change is.

## 4. Batch 2 again, with types and a scale (retrospective)

[`scripts/retro.py`](scripts/retro.py), [`results/retro.json`](results/retro.json). The batch was read this morning.
Nothing here is a forecast.

**The scale** ([`scripts/scale.py`](scripts/scale.py); control cells read before, their isotype columns are new and are
no readout of any unit): how far a readout stands above the background of the panel among control cells, in spreads
of the readout.

| | interferon | co-culture | no signal |
|---|---|---|---|
| PD-L1 | 3.23 | 3.13 | 1.26 |
| MHC-I | 4.95 | 6.77 | 4.29 |
| HLA-E | 0.63 | 0.46 | **−0.25** |

HLA-E stands at the background of the four isotype antibodies, and below it with no signal.

**R1 holds.** With types, the pathway events forecast each other: co-culture 15 of 15 pairs inside, error 0.11
against 0.64 for the control readout; THP-1 4 of 4, 0.18 against 0.38. No event of another type is moved: the 71 pairs
with no edge have the error of the control readout alone, 0.090 (as committed, what the pathway rows added made it
0.115).

**The calibration of the 79 pairs**, under the law as it now is:

| | inside 95% | mean squared *z* (1 is exact) |
|---|---|---|
| as committed before the batch was read | 77 | 0.30 |
| now | 74 | **1.02** |
| now, the 71 pairs with no edge | 68 | 0.85 |

By stratum it is uneven: co-culture 29 of 30; with no signal 35 of 39 (mean squared *z* 1.5: a little tight; two of
the four outside are the pairs part C admits as edges); THP-1 10 of 10 with intervals far too wide (0.01).

**R2 does not hold.** The class's share from the seven pairs of the library, times the scale: error **0.84** on the
ten own-gene pairs of the other two conditions, against 0.68 for one size per class; 8 of 10 inside against 10 of 10.
The shares inside the class run from −0.08 to −0.58: one gene of several behind an antibody (HLA-A, HLA-B, HLA-C)
removes a small part of what it measures, the only gene (CD274) or a subunit nothing is made without (B2M) most of it.

**R3 does not hold as written**: the change divided by the scale varies less across conditions than the change for
5 of 7 pairs, not for HLA-C (two conditions, both nearly equal) and not for HLA-E on its own readout.

**Explored after seeing this: the pair's own share carries.**

| pair | share under interferon | co-culture | no signal |
|---|---|---|---|
| CD274 → PD-L1 | −0.56 | −0.63 | −0.52 |
| B2M → MHC-I | −0.58 | | −0.47 |
| HLA-B → MHC-I | −0.14 | −0.13 | −0.13 |
| HLA-A → MHC-I | −0.08 | | −0.08 |
| HLA-C → MHC-I | −0.13 | −0.08 | |
| HLA-E → HLA-E | −0.37 | −1.04 | +0.73 |
| B2M → HLA-E | −0.30 | | +0.20 |

Forecasting each of the ten pairs at *its own share under interferon × the scale of its condition*:

| error on the ten pairs | all | PD-L1 (2) | MHC-I (5) | HLA-E (3) |
|---|---|---|---|---|
| no edge | 1.01 | 1.47 | 1.06 | 0.30 |
| one size per class (batch 2) | 0.68 | 0.75 | 0.60 | 0.75 |
| the class's share × the scale (R2) | 0.84 | 0.74 | 1.07 | 0.26 |
| the pair's own share × the scale (explored) | **0.23** | 0.16 | 0.25 | 0.25 |

This was looked for after R2 failed, on the same ten pairs. It is a hypothesis for units never read, not a result.
`Law.new(shares=)` and the `share` of an edge record are there so that it can be forecast before a reading.

## 5. Proposed edges, and the step after a rejection (retrospective)

**C.** The law of the library under interferon, given the declared edges, flags five pairs. Each is proposed as an
edge and asked of the rows of the other two conditions, forecast without it (rule G2):

| proposed | asked of | change there | forecast without it | | share of the scale: interferon, there (explored) |
|---|---|---|---|---|---|
| TAPBP → MHC-I | co-culture (batch 1) | −1.72 | +0.02 ± 0.13 | **admitted** | −0.33, −0.26 |
| CD44 → MHC-I | co-culture (development) | +0.41 | +0.04 ± 0.12 | **admitted** | +0.08, +0.05 |
| CD44 → MHC-I | no signal (batch 2) | +0.33 | +0.03 ± 0.08 | **admitted** | +0.08, +0.07 |
| SSR2 → MHC-I | no signal (batch 2) | +0.46 | +0.05 ± 0.10 | **admitted** | +0.09, +0.10 |
| STAT3 → MHC-I | co-culture (batch 2) | +0.19 | +0.01 ± 0.13 | refused | +0.11, +0.03 |
| NONO → PD-L1 | no row elsewhere | | | stays proposed | |

Five proposed pairs, five tests on four of them: three pairs admitted (CD44 in both conditions), one refused, one
with no row to be asked of. With the intervals as committed for batch 2, two of the tests stood outside (TAPBP, and CD44 with no signal); hypothesis
HB5 of batch 2 stays as it was scored, unresolved. The first version of the rule for the widening refused all five
(Section 6): what is in this table was reached after seeing these rows, and is not a test of the rule.

**D.** The 14 accounts of batch 1 and batch 2, through `account.plan` and `outer.open_parts` with the edges of the graph:

| the part outside the input | accounts | |
|---|---|---|
| on the graph: a declared edge (the knockouts of the genes behind the readouts) | 7 | the edge carries its share |
| on the graph: an edge admitted in part C (TAPBP; CD44 in two conditions; SSR2) | 4 | |
| open, not forecast by the recorded events: CMTM6 → PD-L1 in THP-1 | 1 | lever G: propose an edge |
| open, undecided: CD44 → PD-L1 in co-culture, BRD4 → PD-L1 in THP-1 | 2 | lever E: more records of the type |

The other open parts of these accounts (held and not read, lost before the state) are all at the margin: none above a step.

## 6. What was corrected on the way

| found by | the fault | the correction |
|---|---|---|
| the planted study of units (2 wrong of 216) | "The other events forecast it" was read from the lower end of `others` alone: a forecast of 14 times the change counted. | Two-sided: forecast if the interval of `others` lies within a quarter of the change on either side; not forecast if it lies beyond a half on one side. Three readings of run 41 move from "forecast" to "unresolved" (JAK1 on PD-L1 in co-culture, 1.25 [1.10, 1.42]; JAK1 under interferon, 1.14 [0.98, 1.33]; STAT1 in THP-1, 1.54 [1.36, 1.79]); none of batch 2 moves. |
| a planted library (pairs with no change had a mean squared *z* of 0.2) | The error of a row beside its readout added the error of the change and that of the shift. Both are means over the same cells: the shift's error is in both and cancels. The intervals were too wide by that much. | The error is that of what the readout leaves in the cells (`Law`, and the rows of other events in `reading`). This is why the intervals of batch 2 were too wide. |
| part C on real rows (all five proposals refused) | With types, the events of no type chose one widening of the control prior together. Three such events in co-culture, one of them TAPBP, widened it five-thousand-fold, and the interval hid the change. | A widening is chosen by the rows of a declared type. For an event of no type the control prior stands as the control cells give it; what the event shows beyond it is its direct part and is flagged. Tested on a planted library. Changed after seeing part C. |

Also: the flags of the law now carry what `propose_edges` reads (they did not: the two had never been run together),
and the parts of an account have their own table of routes (`ROUTES4`).

## 7. What this does not show, and what is not built

- **Nothing here is confirmed on units never read.** Parts 4 and 5 are retrospective, part 3 is explored, and three
  rules were changed after seeing data of this run.
- The pair's own share needs the pair to have been measured once. A pair never measured has only its class, and the
  class is weak (shares from −0.08 to −0.58). What sets a pair's share (how many genes stand behind an antibody) is
  not in the graph.
- The scale fails where a readout stands at the background of the panel (HLA-E). It is built for one class of edge
  in one study.
- The spread of a direct part with no edge is still one number per stratum: too large in THP-1, a little small with
  no signal.
- Types are declared for one pathway. The law is not run on units with no control cells (the account is).
- A proposed edge between a donor and a protein has nothing to be tested on in a study that measures each donor once.

## 8. Files

`DESIGN.md`; `scripts/scale.py` (the owner's machine), `planted_units.py`, `retro.py`, `donors.py`; `results/planted_units.json`,
`planted_units_first.json`, `retro.json`, `donors.json`.
