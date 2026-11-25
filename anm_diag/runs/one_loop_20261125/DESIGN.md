# One loop for every shape of event: the account carried into the forms that were already there

25 November 2026, afternoon. Nothing is trained. No unit that was not read before is opened.

## 1. What was wrong

The owner asked why the account of a change (run 41, accepted into the method as version 4 of the protocol) was still
listed as missing three things that the method already had a form for. He was right. The forms were there; the account
had not been carried into them.

| what the method already had | where | what version 4 did with it |
|---|---|---|
| events with no control cells: a donor, a batch or a study read against the other units; a disease read as a group of units against another | `events.contrast_events` and rules C1 to C6; `cycle.one_units`, `UnitRows`; planted systems in `toy.py`; run 32 on eight donors | `account.reading` asked for control cells, and resampled cells. Donors were still read by the reading of version 3. |
| the event graph as one table of typed records, with a rule (S3) for an event that acts on the readout itself | `events.py` | The edges of run 41 and of batch 2 were a table in a run's script. Rule S3 still called such an event leakage and sent it to the claim, as a limit. |
| the step after a rejection: every open reading sent to one of the levers M, E, V, S, C, P, with the admissible actions | `outer.ROUTES`, `outer.open_items` | The routes read the verdicts of version 3 ("unreadable", "upstream", "self event"). The four parts of an account were not routed, and "not in the input" had no lever to go to. |

Batch 2 added two faults of the law itself: what the recorded events teach a readout was applied to every event, of
whatever type, and made one readout worse; and an edge had one size per class, where the measured sizes ran from 0.23 to
2.1 spreads.

## 2. What is changed (rules 4i to 4l of the protocol)

**4i. One reading for the three shapes of event** (`account.reading`, `account.account`, `python -m anm_diag account`).
An event with control cells is read as before (on the IRF1 target the numbers are those of run 41 to 1e-13). Where no cell
is a control cell (`control=None`), a unit is read against the cells of every other unit. Where an arm is made of units
(`cluster`), it is resampled by unit, and the readout "with no record of the event" is one of the variation inside a
unit: with no record, nothing says how one unit differs from another. The two controls of the measure are made the same
way. What the rows of other events teach a readout now takes their errors into account (`_teacher`): a row inside its
error teaches nothing, so units that do not differ can sit among the others; with eight rows or more a row may hold a
part of its own, as in the law. With exact rows this is the closing of run 41.

**4j. An edge is a record of the event graph** (`events.edge`, `own_edges`, `propose_edges`, `test_edges`, `edge_table`,
`check_edges`, rules G1 to G5). An edge that follows from the declarations (the site is a gene of the readout's own
measurement, or a declared partner of it: the fields rule S3 already reads) is *declared*. An edge that rests on measured
rows (a flag of the law, or an account that places a change outside the input with no other event forecasting it) is
*proposed*, enters no forecast, and is *admitted* or *refused* by rows it did not see. The law takes its edges from this
table.

**4k. The law by type, and an edge with a size** (`account.Law(types=, scale=)`, `Law.new(types=, scale=)`). The error of a
readout is a property of a type of event: the rows of an event update the readout for the events of its type and for no
other; an event with no type is taught by no row. The direct part of an edge is its class's share of a scale that is
measured in each condition. For the class "the site is a gene of the measured protein" the scale is how far the readout
stands above the background of the antibody panel among control cells: what the knockout can remove.

    scale(readout, condition) = mean over control cells of [log(1 + counts of the readout's antibody) − background] × (1 − 1/24) / spread of the readout among control cells
    background = the mean of log(1 + counts) over the four isotype antibodies of the panel, in the same cell

**4l. The step after a rejection** (`account.plan`, `outer.open_parts`, lever G). Every open part of an account names its
lever and its proposal, and is credited with what the readings already made allow:

| open part | lever | proposal | judged by |
|---|---|---|---|
| held, not read | M, the head | a readout taught by the records of the events of the same type | what it reads on cells it was not fitted on, up to what the state carries |
| held only as the identity of the event | E, records | its own record answers it: the graph's memory | whether other events forecast it |
| lost before the state | M, before the state | a state kept by type; then training (the owner's decision) | what the candidate carries against the mean state and a matched control |
| outside the input, forecast by its type | E, records | the records of the type, asked for by type | the event left out |
| outside the input, not forecast | **G, the graph** | an edge: declared if it follows from the declarations, else proposed | rows the proposal did not see; its size from its class and its scale |
| outside the input, undecided | E, records | more recorded events of the type | the interval of `others` |

An event with a declared or admitted edge to the readout is "on the graph": its part outside the input is no open item.

## 3. Checks on systems with a known answer

`tests/test_one_loop.py` (seven tests): a unit against the other units is placed where it was planted (a factor that
varies inside units: read; dropped by the model: lost before the state; the readout itself: outside the input, not
forecast, lever G; told apart completely by a marker: held as its identity); a disease is resampled by donor; the account
with control cells is unchanged; edges are declared, proposed and tested; records teach within a type and no other; an
edge with a scale is forecast at its share; every open part goes to its lever.

`scripts/planted_units.py`: the planted study of units of `toy.py` at two sizes and twelve seeds each, scored by the rule
in its docstring (written before it is run).

## 4. Real cells, all read before

Everything here is exploratory or retrospective, and labelled so. No unit is opened.

**A. Eight donors** (the cells and the 26 pairs of a donor and a protein that run 32 admitted on a selection third; the
two thirds that were read). The account of every pair at the input view and the state, each donor against the other
seven, the other donors' rows as `others`. *Seen before this design was written*: the first run of the code on real
donors was this. It is reported as explored, with no hypothesis.

**B. Batch 2 again, with types and a scale: retrospective.** The batch was read this morning; nothing about it is a
forecast any more. What is asked is whether the two changes of 4k would have done what the batch said was missing.
Types: the events of the interferon pathway are one type per stratum (IFNGR1, IFNGR2, JAK1, JAK2, STAT1, as declared in
run 40; no other event has a type). Scale: as in 4k, from control cells of each condition (`scripts/scale.py`: control
cells already read; the isotype columns are new, and are no readout of any unit).

- R1. With types, the pathway events of a stratum still forecast each other (left out one at a time, in co-culture and in
  THP-1), and no event without a type is moved by them: the error of the 71 pairs with no edge is that of the control
  readout alone (0.090), by construction.
- R2. The share of the class is taken from the seven pairs of the library under interferon alone. For the ten own-gene
  pairs measured in the other two conditions (eight of batch 2; B2M on MHC-I and HLA-E of batch 1), the forecast "share ×
  scale" has a smaller error than the forecast "one size per class" (that was 0.64 for the eight), and at least eight of
  the ten fall inside their interval.
- R3. The measured change divided by the scale varies less across conditions, for the same gene, than the measured
  change itself.

If R2 or R3 fails, the scale as written is not what sets the size, and that is the finding.

**C. Proposed edges, tested on rows they did not see: retrospective.** The law of the library under interferon, given
only the declared edges, flags pairs. Each flagged pair is proposed as an edge and asked of the rows of the other two
conditions (batch 1 and batch 2), forecast without it: admitted, refused or undecided by rule G2.

## 5. What this does not do

No repair is confirmed on cells never read: this library has none left that could carry one (README of batch 2). The
scale is built for one class of edge in one study. Types are declared for one pathway. The law is not yet run on units
with no control cells (the account is). Nothing is trained.
