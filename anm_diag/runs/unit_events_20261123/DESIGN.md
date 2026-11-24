# Events without control cells in the round: design

Written on 23 November 2026. It is step 1 of the plan in [`GENERAL.md`](../../GENERAL.md) (Section 7): "events without
control cells (donor, disease) in the round: design, planted donors, forecasts written first". Parts 1 to 4 say what
the round does for such events. Part 5 is the rule by which the planted systems are scored, written before the scored
runs. Part 6 declares one round on real cells, with its forecasts, **before that round reads anything**.

## 1. What was missing

A knockout comes with control cells: cells of the same experiment that did not get the event. The round uses them
three times. Admission compares an event's cells with them. The repair starts from a readout fitted on them. And every
reading of *u* takes its slope from them. A donor has no such cells. The audit tool reads a donor against the other
donors and resamples both sides by unit; the round could not (GENERAL.md, Section 6, first item).

There are two shapes, and they need different things.

| | a knockout | a unit against the other units (a donor, a batch, a study) | a group of units against another (a disease) |
|---|---|---|---|
| the other side of an event | control cells | the cells of every other unit | the cells of the healthy donors: they are the control cells |
| how the other side is resampled | by cell | by unit | by unit |
| how the event's side is resampled | by cell, or by event in a pool | by cell | by donor (the patient donors are the events of a pool) |
| admission (A1) | the change against its interval, cells as units | the same, the other side counted by unit | the same |
| what the repair starts from | a readout fitted on control cells | a readout of the variation among the cells of a unit | a readout fitted on the healthy donors' cells |
| what a row is | an event against the control cells | a unit against the other units of the fit | a patient donor against the healthy donors' cells |
| code | `cycle.one_type` | `cycle.one_units`, `cycle.UnitRows`, `select(control=None)` | `cycle.one_type(control_cluster=...)` |

The second column is new. The third is the knockout round with one change: the side an event is read against counts
by unit, in admission and in every reading.

## 2. The repair for a unit against the other units (`cycle.UnitRows`)

A plug-in of the same shape as `LinearRows`: a linear readout of the stage the model reads, on the scale of the
system's own output. Every unit is in the fit; a fit of one half holds half the cells of every unit.

| arm | what it is | what it answers |
|---|---|---|
| as trained | the system's own output | what the model does today |
| `controls` | a readout of the variation among the cells of a unit (the stage and the readout centred by unit) | what a readout reads with no record of how any unit differs from the others. It takes the place of the readout of control cells |
| `rows` | the same, updated once by the row of every unit, the event's own among them | the round's repair, as for a knockout: the event's own record is in it |
| `cells` | fitted on every cell of the fitted half, not centred | the ceiling of this kind of repair |
| `unseen` | the arm `rows` of a fit that holds no cell of the event's unit: its readout, its rows and its centre are those of the other units | does the repair carry to a unit it has not seen |

`rows`, `controls` and `cells` are read against a reader of the stage that had the unit's own cells (the scope `own` of
the audit tool). `unseen` is read against a reader that never saw the unit (the scope `units`). Each event is read on
its own, its cells by cell and the other units by unit, on the resamples of the audit tool.

## 3. What counts as what

Margin 0.10, step margin 0.10, 1,000 resamples, as for every round. The grades and statuses are those of the protocol
(version 3); nothing is amended.

- **The model's part is closed for a unit**: A3 holds, gate 3 passes, and `rows` adds no loss beyond its reference
  (relative), or is supported (absolute).
- **The repair carries to a unit it has not seen**: `unseen` adds no loss beyond its own reference.
- **The unit's own row is needed**: `rows` reads at least a step less than `unseen`, on paired resamples. The outer
  step then raises an item for that unit, "pooled, and differs" (lever S): the type table treats the units as one
  type, and the rows of the other units do not reach this one.
- **Where a loss before the model is**: placed as for a knockout, from the grid of the event.

A unit event is always read on its own. Nothing is said about the family of units from fewer than ten of them.

## 4. What changed in the general code

- `scope.py`: with several readouts a run can also return the grid of every readout on its own (`each_readout`), on
  the same resamples, so that one run of the tool serves a panel of proteins; a ridge reader can choose its weight per
  readout (`alpha_per_target`). Control cells that come in units of their own keep them (`keys.control_block`). All
  three are off unless a SPEC asks.
- `cycle.py`: `screen` and `select` count the other side by unit where it is made of units; `output_u` and `read`
  resample it by unit (`control_cluster`); `UnitRows`; `one_units`; `one_type(control_cluster=...)`; a reading of what
  the own row adds over `unseen`.
- `outer.py`: the item named above.
- `toy.py`: a planted study of units and a planted disease.

No name of a system enters any of them. Every earlier result is unchanged: the options default to what the code did
(the tests of the first system and of the planted knockouts are run again).

## 5. The planted systems, and the rule they are scored by

`anm_diag/toy.py`. The model, its two steps and its head are those of the planted knockouts. One more block of
modality a carries a marker that varies only between units.

**A study of units** (`run_units`). Most units are ordinary. Nine differ from the others, each in a known place:

| unit | differs in | the loop should say |
|---|---|---|
| h1, h2, h3 | factor 1, which the head does not read | A3 holds; closed by the repair; `unseen` closes it too (no item that the unit differs) |
| early | factor 5, dropped in the first step | upstream, lost between the input and mid |
| enc | factor 2, dropped in the second step | upstream, lost between mid and z |
| modb | factor 3, which only modality b carries | unreadable |
| post | the readout itself | unreadable |
| solo | the readout itself, and the marker with it | a loss in the head (A3 holds) and nothing placed before the model; the rows of the other units do not carry: "pooled, and differs" |
| weak | factor 1, barely | refused at selection |

**A disease** (`run_groups`): twelve healthy donors and eight patients. Every donor has a level of its own; the disease
moves factor 1. The loop should read the patients as a pool, by donor on both sides, and close it by the repair.

A second version of each has nothing wrong (the head reads the right factor, and only units it can read differ).

The rule, as for the planted knockouts ([`framework_generalization_20261122`](../framework_generalization_20261122/)):

- **right**: the loop says what was planted;
- **undecided**: the loop says less (no item, an interval that does not decide, a weaker status, a unit not admitted);
- **wrong**: the loop says something else that is definite.

For h1, h2 and h3 a grade other than relative or absolute is undecided: the loop then claims no closure. "Upstream"
or "unreadable" said of a unit whose loss is in the head (h1, h2, h3, solo) is wrong. For `solo` no closure is asked:
one row updates the readout only part of the way (on seed 0, while the code was written, its own row removed about
two thirds; the same was seen for IRF1 on real cells). What is scored for it is that the head has a step to remove
and nothing is placed before the model, and, on its own, the item that it differs: present is right, absent is
undecided; the same item for another unit is wrong.

Facts scored for the study of units: one per row of the table above (h1, h2 and h3 separately; for `early` and `enc`
also between which two stages; for `solo` the two named above), "no other definite item", and
for the version with nothing wrong "nothing is open and no repair is read". For the disease: "closed by the repair"
and, for the version with nothing wrong, "nothing is open". A unit that is not admitted is undecided, except `weak`,
for which admitted is wrong.

Three sizes, twelve seeds (1 to 12; seed 0 was used while the code was written and is not scored):

| size | ordinary units | cells per unit (every unit the same) | healthy and patient donors | cells per donor |
|---|---|---|---|---|
| large | 72 | 300 | 12 and 8 | 400 |
| medium | 40 | 200 | 10 and 6 | 200 |
| small | 20 | 100 | 8 and 5 | 80 |

Why so many ordinary units: the other side of a unit event counts by unit, and a unit that the stage cannot read moves
the mean of the other side when it is drawn twice or not at all. With few units no reading reaches the margin, and the
loop then says nothing definite (which is the rule working, and is scored "undecided").

The property asked for is the one asked of the planted knockouts: `wrong` stays at zero as the data shrink.

The readers of the planted units are ridge only. A unit event is fitted once for every pair of units, and boosting
would take hours. So a reading that needs a second reader class (naming another measured input) is not scored here;
it was scored on the planted knockouts.

## 6. One round on real cells: eight donors

**What is read.** The PBMC of eight healthy donors of Hao et al. 2021 (GSE164378), 5,000 cells each, 112 surface
proteins measured by CITE-seq. The system has never been trained on a cell of this study (its head was fitted on bone
marrow). These are the development donors of round 1 ([`round1_coverage_20261116`](../round1_coverage_20261116/)) and
of the scope run of 16 November: they have been read before, and **everything here is exploratory**. No final unit and
no reserve donor is read. No state is computed and no model is trained.

**What is not read, and why.** The bone-marrow donors: the head was trained on them, and donor and site vary
together there. The COVID-19 study of round 1's test (Stephenson et al. 2021): it is the final study of that round;
40 of its donors were read once as a test, 73 are a reserve that is the owner's to open. So the disease shape is
checked on planted cells only, and no disease event is read on real cells in this run.

**The record.** A donor is a unit; no cell is a control cell. Donors P1 to P4 were processed in batch 1 and P5 to P8
in batch 2: **donor and batch vary together**. A donor read against the other seven differs from them by its own
level and by about half the difference between the batches. Every target of this round carries that as an open item
("confounded", lever E), and nothing here is attributed to the donor alone.

**Roles.** All eight donors are development units. A cell is in the selection third if SHA-256 of
`20261123-units|<barcode>` is divisible by 3; those cells enter no reading.

**The order.**

1. Commit this file and the forecasts (`results/declared.json`).
2. `scripts/build.py`: the arrays (the state, 512 components of RNA, the same with rows permuted, the 112 proteins, the
   head's own outputs) without the selection third; admission of every donor for every protein on the selection third
   (`cycle.select(control=None)`); the SPEC.
3. `anm_diag scope`, one run with `each_readout`: the stages are the 512 components (the model's input view), the state
   Z (what the head reads) and the permuted view (control); the scopes are the system's own output, a readout of its
   training cells, of the other donors, and with the donor's own cells. Ridge only, with the ridge weight chosen for
   every protein on its own, so that a protein's reader is the one a run with that protein alone would fit. Boosting
   is fitted per readout, and 112 proteins times 36 fits is out of reach. So a reading of "unreadable" here rests on a
   linear reader, and the outer step says so in the item.
4. `scripts/round.py`: for every protein with an admitted donor, `cycle.one_units`; the targets and the open items.
   One reading; a result file that exists is not overwritten.
5. `scripts/score.py`: the forecasts against what was read.

**The quantities that are forecast** (each over the admitted pairs of a donor and a protein):

| | quantity |
|---|---|
| N | the number of admitted pairs, and of proteins with at least one |
| T | the share whose output as trained is rejected |
| Z | the share whose reference (a reader of Z with the donor's own cells) is supported; the share rejected |
| A3 | the share for which the head has a step to remove |
| R | among those: the share whose grade is relative or absolute; the share absolute; the median fall of \|u\| from the output as trained to `rows` |
| U | among those: the share for which the donor's own row is needed (a step over `unseen`); the share for which it is not (no step); the median of \|u\| of `unseen` minus \|u\| of `rows` |
| C | among those: the share for which the rows add a step over `controls` |
| P | among the pairs whose reference is rejected: the share read as lost at the input view ("unreadable") and the share read as lost between the input view and Z ("upstream") |
| K | the control view reads nothing for every pair (\|u\| at least 0.9); some reading of the family is supported |

**Why the forecasts are what they are.** On 16 November, over all 112 proteins at once, the head left 0.78 to 0.93 of
a donor's shift, a reader of Z fitted on the other donors 0.30 to 0.52, and one with the donor's own cells 0.20 to
0.34. So the head loses most of what Z carries (A3 should hold often), Z itself rarely reaches the margin (few
absolute grades), and a donor's own cells add something over the other donors' (the own row should matter for a
good part of the pairs). With eight donors the other side has seven units: intervals are wide and every one is
approximate by the protocol's own number.

## 7. What would change the plan

- The tool's reading of the output as trained and the round's do not pair for a protein: stop, report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.
