# The loop: protocol

Written 12 November 2026, after three exploratory rounds on the teddy_mm head
([`runs/loop_20261112`](runs/loop_20261112/), [`runs/loop_guard_20261112`](runs/loop_guard_20261112/),
[`runs/repair_20261112`](runs/repair_20261112/)) and a re-reading of their records
([`runs/loop_reading_20261112`](runs/loop_reading_20261112/)). Those rounds decided pass or fail by comparing a point
value with a threshold. Their audits can show that a system loses a response. They can rarely show that a changed system
no longer loses it. This protocol separates the two questions and puts a cheap check in front of every expensive step.

It is general: a system with a chain of stages, a source, a readout, and a change to the system. The teddy_mm head is the
first instance. The readings are made by [`anm_diag/loop.py`](anm_diag/loop.py) from finished runs; it makes no fits.

## Two loops

An audit locates a step: the response is readable before it and not after it. Two different things can be done next, and
they need different evidence.

| | diagnostic loop | repair loop |
|---|---|---|
| question | Is the located step the one that matters for this response? | Is the changed system free of the loss, and still good at its own task? |
| what is compared | two changes against one baseline: one made at the located step, one made elsewhere | the changed system against a control that differs only by the change |
| quantity | the share of the baseline's loss over a span that each change removes, and the difference, each with an interval | the loss the changed system still adds over the span, and the verdict of the span's end stage |
| kind of claim | one change removes more than the other | the loss is absent |
| what it needs | a baseline loss that is large against its interval | intervals narrow against the margin; data no design decision has seen |
| command | `anm_diag loop diagnose` | `anm_diag loop closure` |

The diagnostic loop tests the audit's statement. The repair loop tests the changed system. Passing the first says where to
work. It does not say the system is repaired.

## Quantities

- *u* at a stage: the unexplained share of the response (as in every run). Stage verdicts in the loop are two-sided:
  support needs the whole interval of some declared readout class inside (−margin, margin).
- **Loss added over a span** a → b: |*u*<sub>b</sub>| − |*u*<sub>a</sub>|, computed on every bootstrap draw. Read three ways
  against the step margin: **loss** (lower bound above it), **no added loss** (upper bound below it), **unresolved**.
- **Removed by a change**: the baseline's added loss minus the changed system's, draw by draw. Runs of different systems on
  the same cells, folds and seed share their resamples; the tool checks this and says when a contrast is not paired.
- **Share removed**: one minus the ratio of the two, with its interval.
- The span of a loop runs from the stage before the located step to the system's output. A repair counts by what reaches
  the output. Single steps are reported beside it.

## Readings

`diagnose`, per readout:

- **located step confirmed**: the baseline adds a resolved loss over the span; the change at the located step removes a
  resolved part of it, not confined below the step margin; and it removes more than the control change (interval of the
  difference above zero). *Selective* if whatever the control change removes lies below the step margin.
- **located step not confirmed**: what the located change removes lies below the step margin, or the control removes more.
- **unresolved**, or **not testable** (the baseline adds no resolved loss over the span).

`closure`, per readout, strongest first:

- **absolute**: the end stage of the changed system is supported and the control's is rejected. This is the paper's
  closing of a loop (a rejected stage becomes supported, and the control's does not).
- **relative**: the changed system adds no loss over the span, by its interval, and the control adds a loss.
- **loss no longer resolved**: the changed system's span is unresolved. The loss is not shown. Its absence is not shown
  either. This is not a closure.
- **not closed**: the changed system still adds a resolved loss. **Not testable**: the control shows no loss to close.

## Gates

Each gate is passed before the next step is paid for.

| gate | what is checked | by | cost |
|---|---|---|---|
| 0 admission | The response is resolved and large against its interval. The reference stage is not rejected. The span adds a resolved loss as it stands. How narrow the intervals are: whether "absolute" or "relative" closure could be reached at all. The number of clusters. | `anm_diag loop admit` on the baseline run | none: one existing run |
| 1 location | Sealed declaration, run once. | `anm_diag seal`, `run` | one audit |
| 2 prescription | One change per round. A control: the same procedure without the change, or the same change made elsewhere. The level of generalization claimed (cells, sources, readouts, cell types never used by the change). The span, the readouts and the rule, written before any run. | a design file, committed first | none |
| 3 admission of the change | It converged (a stated criterion, not a step limit). It reached its own objective on the data it was fitted on. No audit cell entered the fit (checked by identifiers). The guard holds on validation data. | the round's fitting script | low; no audit |
| 4 development audit | `diagnose` and `closure` on development data. Rounds may repeat here. Everything read here is exploratory. | `anm_diag run`, `loop` | audits |
| 5 final audit | The recipe is frozen. At least three training seeds per arm. Sealed declaration, run once, on data no earlier gate has used. The outcome must hold in every seed. | `seal`, `run`, `loop closure` | training and audits |
| 6 use | The guard on the system's own task: the mean, and the worst single item. A comparison with a baseline that does not use the changed part. | the system's own evaluation | low |

A readout that fails gate 0 is not a target and is not counted in any outcome. A change that fails gate 3 was not tested.
It did not fail an audit.

## Rules

1. **Decide by intervals.** No outcome rests on a point value compared with a threshold.
2. **One change per round.** Two changes need a factorial design with every arm audited.
3. **Keep the audit's data for the audit.** Prefer a change that does not use the audited sources. If it must use them,
   cross-fit: every audited cell passes through a system fitted without it, and every control cell is used.
4. **Three sets of data.** Fit, development and final. Final data are touched once, under a sealed declaration. Without
   final data the loop stays exploratory, and says so.
5. **Replication is part of the decision.** The rule names the seeds or splits it will be read on.
6. **Stop.** A round states how many further rounds the development data will carry.
7. **One ledger.** Every change, audit and decision is a typed event (`loop.intervention_event`, `audit_event`,
   `decision_event`), so a run is not made twice and a decision can be traced to its runs.

## Evidence grades

Confirmatory: gate 5, sealed, enough clusters. Everything else is exploratory, including every re-reading of a finished
run under this protocol. A re-reading changes no outcome recorded under the rule a run was made with.

## Chosen numbers

| number | value | where |
|---|---|---|
| margin on *u*, step margin | 0.10 | each declaration |
| clusters for intervals that are not marked approximate | 10 | `loop.DEFAULTS` |
| response against the half-width of its interval | 3 | `loop.DEFAULTS` |
| span half-width for "relative" closure to be attainable | half the step margin | `loop.admit` |
| guard on the system's own task | 0.02 on the mean | [`declarations/teddy_mm_step4_intervention_rule.md`](declarations/teddy_mm_step4_intervention_rule.md) |

None of these is derived. Each is reported with every reading.

## What the first three rounds lacked

| lacking | consequence | gate or rule |
|---|---|---|
| a check that a repaired head could be recognized | on THP-1 the reference stage is unresolved, so no round could close there; on melanoma the intervals are about as wide as half the margin | gate 0 |
| readouts admitted by their response | HLA-E and MHC-I were made targets with a response too small to follow | gate 0 |
| a rule on intervals | round 3 "failed" by 0.009, inside every interval | rule 1 |
| a check of the fit | the fits stopped at a step limit; HLA-E was not learned even where it was trained on | gate 3 |
| one change per round | round 3 added a cell type and three outputs together | rule 2 |
| audit data kept for the audit | half the controls went to the fits, and the reference stage fell from supported to unresolved | rule 3 |
| final data | the same audit cells served every round | rule 4 |
| the sealed design of 9 November | three seeds, a matched control retrain and an interval rule were already written there | gate 5 |

---

# Version 2: scope, hypotheses and tests

Added 16 November 2026. Version 1 above stays as it is: it locates the stage where a response is lost, for a
perturbation, and reads a finished loop after the fact. Three things were missing for a loop whose aim is a system that
works on units it has not seen (donors, batches, studies, diseases):

| missing in version 1 | added | where |
|---|---|---|
| events other than a perturbation | contrasts between groups nobody assigned, with their own rules | `events.contrast_events`, `events.check_contrast` |
| what is missing, not only where | a second axis at every stage: the scope of the readout | `scope.py`, `anm_diag scope` |
| a forecast written before a change, and a rule to read it | hypothesis and test events | `loop.hypothesis_event`, `loop.test`, `loop.ledger_check` |

Everything is added to the same package and the same ledger. No model is trained by any of it.

## One quantity

*u*(event, readout, target, stage, scope): the share of a response that a prediction leaves unexplained.

- **event**: a treated unit against control cells. A knockout against unassigned cells; a donor against the other donors
  of a study; a batch; a study; patients against healthy donors.
- **target**: `shift`, the difference between the two sides (what version 1 audits), or `variation`, the differences
  among the cells of the treated unit (what a system's ordinary score measures).
- **stage**: a view of the cell along the chain, as in version 1.
- **scope**: the cells the readout was learned from.

| scope | learned from | a response read only from here on means |
|---|---|---|
| `frozen` | nothing: the system's own output | — |
| `train` | the system's training cells | the map is at fault, not its data: **function** |
| `other` | another study of the same kind of event | (reported beside the others: does one study teach another) |
| `units` | the other units of the audited study, never the unit that is scored | cells of this kind are needed: **coverage** |
| `own` | the audited unit's own cells too, cross-fitted by cell | it is specific to the unit: **unit** |

If the unit's own cells do not read it at the stage the system reads, and they do at the input view, the loss is
**upstream**. If they do not there either, it is **unreadable** by the declared views and classes. Version 1's
cross-fit by cluster is the `units` scope; its verdicts carry that scope in the ledger.

One statistic serves every scope. With *P* the scope's out-of-sample prediction of the readout *O*:
θ = Δ(*O*) − *a*·Δ(*P*), with *a* the slope of *O* on *P* among the control cells, and *u* = θ / Δ(*O*). The slope
puts a prediction made on another scale (a head trained on another study's protein scale) on the scale of the
measurement. For several readouts, *u* is the norm of the unexplained differences over the norm of the differences,
each in units of the control spread; for one readout that is |*u*|. The same margin (0.10) applies to every scope.

## From a reading to an action

| reading | what to do | what not to do |
|---|---|---|
| no deficit | nothing | — |
| function | change the map (architecture, objective), same data | add data |
| coverage | put cells of this kind in the training | change the map |
| unit | calibrate at use, with a few labelled cells of the unit | train for it |
| upstream | change what comes before the stage (the encoder, the pooling) | change the head |
| unreadable | change the measurement, or drop the readout from the score | train |

The gaps between consecutive scopes are reported with paired intervals whether or not a verdict is reached. They are
what a forecast is made of: `scope.forecast` gives, per kind of change, how much of |*u*| at the system's output it
could remove. A forecast is made with readouts. What a changed system delivers has to be audited.

## Events that nobody assigned

A knockout is assigned by the experiment. A donor, a batch, a study or a disease is not, so a contrast must say what
else differs between its two sides (`events.check_contrast`):

- C1 the kind is one of `donor`, `batch`, `study`, `disease`, `condition`
- C2 the labels come from the design, not from the measured response
- C3 `confounded_with` is declared, even if empty
- C4 both sides come from one study, unless the kind is `study`
- C5 pooled units are declared as pooled
- C6 fewer than 10 units on a side: an interval over units is marked approximate

## Hypotheses and tests

A diagnosis is written as a hypothesis before any change is made:

- what is measured on every arm, as a gain;
- the claim (function, coverage, unit, upstream, unreadable, or a location);
- a forecast interval for every arm;
- a rule: the own arm's gain lies above a stated value and every control arm's gain lies below its bound;
- the runs the forecast was read from.

A test reads the arms against the rule: **confirmed**, **refuted**, **unresolved**, or **not testable** (an arm is
missing, or a change did not pass gate 3). Every arm is also scored against its forecast, and `loop.forecast_record`
keeps the record of how the forecasts held up: the auditor is audited by its own ledger.

At least two arms. The usual three for a coverage claim: cells of the missing kind added (own); the same number of
cells of a kind already covered (amount control); another map on the same cells (the rival, function). This is the
shape of the knockout loop, where a repair at the located step was compared with a repair elsewhere.

Ledger rules (`loop.ledger_check`): H1 a test names a hypothesis in the ledger; H2 a test uses no run its hypothesis
was read from, and a hypothesis written after its outcome was known is marked retrospective; H3 a test has the own arm
and a control arm; H4 an arm that was run has a forecast; H5 a hypothesis without a test is listed as open.

## Intervals and units

- An event is resampled by its clusters where it has several, otherwise by cell; controls by block where they come in
  blocks (donors), otherwise by cell. The same resamples serve every stage and scope of an event, so gaps are paired.
- A statement about donors or studies rests on the events of a family (`scope.family`): the mean over events with an
  interval from resampling the events, and how each event was read. Below 10 events it is marked approximate.
- A forecast and its test use different units. One primary hypothesis per round.
- Units set aside as final are named under `roles.final`. A run that touches them is refused unless it is a final audit
  that names the hypotheses it tests.

## The loop in five actions

1. **Observe.** The grid over stages and scopes, for every event of a kind. No training.
2. **Hypothesize.** Rival readings, each with a forecast for every arm.
3. **Intervene.** One arm per reading, and an amount control. Gate 3 for every arm.
4. **Test.** Which reading survives; how far each forecast was off.
5. **Update.** The grid again on the changed system; the next round, or the final audit of version 1 (gate 5).

Version 1's two loops are still two claims. The diagnostic loop is steps 1 to 4: a test of what the audit said. The
repair loop is a claim about the changed system and needs final data.

## The tool is tested on planted problems

`tests/test_scope.py` builds six situations in which the answer is known (no deficit, function, coverage, unit,
upstream, unreadable) and requires the reading to name each. `tests/test_hypothesis.py` writes one study twice, as a
perturbation audit and as donors of a study, and requires the same file format and the same code to read both.

## What version 2 does not do

- It does not make a system general. It can also find that every study needs its own cells.
- With few studies, a statement about studies in general is descriptive. Ten is the number chosen for an interval.
- A forecast comes from ridge or boosting readouts. A trained system can do better or worse than its forecast; that is
  what the forecast record is for.
- A contrast nobody assigned is confounded. The tool asks for the confounders to be named. It cannot remove them.

## After the first use (runs/scope_20261116)

Two finished pieces of work were read again with this version: the knockout loop and the donors. Four things changed in
how it is to be used.

1. **A forecast is a ceiling.** It is what a readout of a scope reads at the stage the system reads. The repairs that
   had been made delivered about two thirds of it. A rule asks for less than the ceiling, and the forecast record
   keeps the share that was delivered.
2. **A forecast and its test are about one quantity.** For the system's own answer that is the `frozen` scope before
   and after the change. The loss at a stage in the `units` scope (version 1) is another quantity and has its own
   forecast.
3. **A margin is declared per target.** 0.10 is defined for *u* of one response. For the variation among a unit's
   cells, and for a shift over many readouts, no scope came within it, and every donor event was read `unreadable`.
   Where no scope is supported the reading is in the gaps, and a round states the margin it will use on that scale.
4. **A readout is named in the events.** `readout_name` in a SPEC; without it verdicts of different readouts share a
   key in the ledger.

# Version 3: five amendments (22 November 2026)

The round of 20 November ([`runs/teg_loop_20261120`](runs/teg_loop_20261120/)) met five cases that versions 1 and 2
could not evaluate as written. [`OUTER_LOOP.md`](OUTER_LOOP.md) (Section 8.5) proposed an amendment for each. They
apply from the outer step of 22 November ([`runs/outer_step_20261122`](runs/outer_step_20261122/)) on. They change no
outcome recorded under the rule a run was made with.

| | the case | the rule from now on |
|---|---|---|
| 3a | Gate 3 asked a change to bring the output inside the margin on the cells it was fitted on. Where the reference stage itself leaves a third of the response, no change of the head can do that. | **Gate 3 is stated against the reference stage.** On the cells it was fitted on, the changed output adds no loss beyond the reference stage: the interval of \|*u* of the output\| − \|*u* of the reference\| lies below the step margin. The output as trained must fail the same gate, or there was nothing to change. The absolute value is still reported. |
| 3b | A hypothesis about a repair by rows had one control, the repair of another group. Where the other group had no admitted event there was no control, and the test could not be made. | **Two controls.** The readout fitted on the group's control cells without rows always exists; it shows what the rows add. The repair of another group shows what the type adds, where there is one. The rule of the test is unchanged: the own arm above its bound, every control arm below. The difference between the rows arm and the controls arm, on paired resamples, is reported and read against the step margin. |
| 3c | A group passed admission as a pool while no single event did. | **Not a target.** A repair by rows needs a row. The reason is recorded as an open item (a small response). |
| 3d | The response of an event was screened on the cells that were then audited. | **A selection split.** One third of an event's cells, and of the control cells, drawn by cell before anything is read, serves admission (A1 of the round, the no-leakage rules). Those cells enter no reading of *u*. Events admitted before this version keep their admission, and how they fare on a third is reported. |
| 3e | The type of a cell was assigned by a reader of its state, and four in a hundred knockout cells of one group got the wrong one. | **The type comes from the record** (the design labels). The reader of types is a component of the system only where no record exists at use, and is then audited like one. |

# Version 4: the account of a change (25 November 2026)

Run 41 ([`runs/account_20261125`](runs/account_20261125/)) found that the reading of a repair used since run 40 could
not see what a stage holds, and that the events that need the event graph had never been targets. It rebuilt the
reading ([`anm_diag/account.py`](anm_diag/account.py); [`GENERAL.md`](GENERAL.md), Section 8). The owner accepted it
into the method on 25 November. It applies from batch 2 ([`runs/batch2_forecast_20261125`](runs/batch2_forecast_20261125/))
on, and changes no outcome recorded under the rule a run was made with.

| | the case | the rule from now on |
|---|---|---|
| 4a | A repair was read by letting a readout follow the mean change of an event on one half of its cells and reading the other half. That number is near one whenever the event's shift of a stage can be found again, whatever the stage holds. | **The reading of a round is the account.** For an event, a readout and a stage, on halves of cells: `read` (a readout of control cells), `carried` (what the stage holds, from the evidence inside the arms), `others` (what the recorded events forecast). The readout made to follow the event's own mean is still reported, as `found_again`, and is not evidence that a stage holds a change. |
| 4b | Gate 3 asked whether a changed head adds no loss beyond a reference stage. A head that follows the mean passes it on a stage that holds nothing. | **A repair is credited up to what the stage carries** (`account.credit`). What it closes beyond `carried` at the stage it reads is not credited; a step beyond means it has followed the mean. Gate 3 of version 3 is kept as a report. |
| 4c | A deficit was placed by comparing *u* at stages (head, state, upstream, unreadable). | **A change has four parts, cell by cell** (`account.chain`): not in the input, lost before the state, held and not read, read. A part is *open* when its interval lies above a step, *unresolved* when it holds a step. Each open part names its place: an edge, the model before the state, the head. |
| 4d | "Unreadable" was a dead end, and rule S3 kept an event that acts on the measured protein out as leakage. | **A part that is not in the input is asked of the other events.** If they forecast the change (`others` at the state: the interval within a quarter of the change on either side, 0.75 to 1.25), their records teach the readout. If they do not (the interval beyond a half on one side: below 0.50 or above 1.50), no change of the model can supply it: it needs an edge between the event and the readout in the event graph. S3 stays a rule of admission for a *repair of the model*; such an event is a target *of the graph*. |
| 4e | A library was read one target at a time. | **A library is read by the law of its rows** (`account.Law`): every pair of an event and a readout is forecast with its event left out, with an interval; pairs with no edge whose direct part stands out are flagged (a missing edge when the event does not move the stage, a direction the law does not cover when it does). |
| 4f | Units never read were opened after a screen of their protein. | **A batch never read is opened after its forecasts are committed.** `Law.new` forecasts every pair from the RNA of the cells alone; the commit precedes the export of any measured value; the batch is read once and scored by the share inside the intervals and by the error against "no change" and against the control readout alone. What was known of the batch before is written in the design. |
| 4g | Where every cell of an arm is edited, no cell of one arm is like a cell of the other, and `carried` read from cells says "direct". | **Below an overlap of 0.10 a stage is said to hold the change as the identity of the event**, and whether it is more than the identity is asked of the other events. |
| 4h | The selection third (3d) cost a third of the cells of every event. | **The account needs no admission.** Its shares carry the error of the change in their intervals; a pair is read when its change is at least four errors, and that line is reported. 3d stays for the admission of a target of a repair. |
| 4i | The account asked for control cells and resampled cells. A donor, a batch or a study has no control cells, and a disease is a group of donors: they were still read by version 3 (run 32). | **One reading for three shapes of event** (`account.reading`). With control cells: as before. A unit against the other units (`control=None`): the unit is read against the cells of every other unit. A group of units against another (`cluster`): an arm made of units is resampled by unit, and the readout "with no record of the event" is one of the variation inside a unit. The rows of other events teach a readout by their errors (a row inside its error teaches nothing). |
| 4j | The edges of run 41 and of batch 2 were a table in a script, and rule S3 sent an event that acts on the readout itself to the claim, as a limit. | **An edge is a record of the event graph** (`events.edge`; rules G1 to G5 of `events.check_edges`). *Declared* if it follows from the run's declarations (the site is a gene of the readout's own measurement, or a declared partner). *Proposed* if it rests on measured rows (a flag of the law; an account that places a change outside the input with no other event forecasting it): it enters no forecast until rows it did not see were forecast without it, and is then *admitted* (outside again, on the same side) or *refused* (inside). The law takes its edges from this table. |
| 4k | What the recorded events taught a readout was applied to every event, and made one readout worse in batch 2; the widening of the control prior was one number; an edge had one size per class, where the sizes ran from 0.23 to 2.1 spreads; the error of a row beside its readout counted the error of the shift twice. | **The law is by type** (`account.Law(types=, scale=)`). The rows of an event update the readout for the events of its type and for no other, and choose that type's widening; an event with no type keeps the control prior as the control cells give it, and what it shows beyond is its direct part. **An edge has a size: a share of a scale** measured in each condition (for a knockout of a gene behind the measured protein: how far the readout stands above the background of the panel among control cells). The share is the pair's own where the pair was measured in another condition (`Law.new(shares=)`; explored in run 43, not yet confirmed on units never read), else its class's, which is weak. The error of a row is that of what the readout leaves in its cells. |
| 4l | The routes of the outer step read the verdicts of version 3. The four parts of an account were not routed, and "not in the input" had no lever to go to. | **Every open part names its lever, its proposal and its judge** (`account.plan`, `outer.open_parts`, `outer.ROUTES4`). Held and not read: the head, a readout taught by the records of the type, credited up to what the state carries. Held as the identity: its own record, the graph's memory. Lost before the state: a state kept by type, then training (the owner's decision). Outside the input and forecast by its type: the records of the type. Outside the input and not forecast: **lever G, the graph**, an edge declared or proposed (4j). Undecided: more records of the type. An event with a declared or admitted edge to the readout is *on the graph*: its part outside the input is no open item. |

Rules 4i to 4l were written on 25 November after the owner pointed out that the method already had a form for each
(events with no control cells, the event graph, the routes of the outer step) and that version 4 had not been carried
into them ([`runs/one_loop_20261125`](runs/one_loop_20261125/)). The same run made the rule of 4d two-sided (a
forecast of many times the change is no forecast); three readings of run 41 move from "forecast" to "unresolved" under
it, none of batch 2.

What does not change: *u*, the margin and the step (0.10); the roles of units (final, development, unread); a
declaration committed before a reading; a reading made once; nothing final is read except as a final audit.

**Version 5 (26 November 2026): the rules are frozen** (commit `d777070` on the main branch;
[`runs/freeze_v5_20261126`](runs/freeze_v5_20261126/)). Until now every turn of the loop also changed a rule, so
that it could not be told whether the loop corrects itself or its author corrects the loop. From this version on a
batch is read with the rules as they stand at that commit. A rule changes only between batches and by a new version, and
a batch that forces a change counts as development, not as a test.

| rule | what was open in version 4 | version 5 |
|---|---|---|
| 5a | The spread of a direct part with no edge was fitted over all pairs: a few pairs with a large part of their own, or a stratum whose recorded events were chosen for their effect, widened it (THP-1: mean squared *z* 0.01), and it fell slowly where it should be near zero. | **By calibration** (`account._calibrated_spread`): the spread at which the recorded pairs with no edge and no type, each left out, are as wide as they should be (the median of their squared *z* is that of a standard normal). A stratum with fewer than eight such pairs takes the library's (`Law.fit(fallback=)`). |
| 5b | A row that stands far out still taught the readout once the spread was small. | **A pair with a direct part of its own is an edge the graph does not hold yet**, and by how likely that is its row teaches less: left out, it stands where a part of the order of one spread of the readout would put it (`OWN_PART`). It is flagged; it does not bend the forecasts of the others. |
| 5c | The scale failed where a readout stands at the background of its panel (HLA-E), and that was found after the reading. | **Said beforehand, from control cells** (`events.by_scale`): where the scale is below one spread of the readout (`MIN_SCALE`) the edge is of the class "…, at the background", has no scale, and its size is asked of the recorded pairs of that class, in spreads. |
| 5d | The class's share of a scale was refuted (run 43); the pair's own share was explored. | **The size of an edge**: the pair's own share where the library holds one on a scale that measures, with the spread by which the shares of pairs measured in two conditions moved (`events.carried_spread`; the class's spread while fewer than three pairs were measured twice); else its class's share. Explored when frozen: batch 3 is its first test. |
| 5e | The law asked for control cells: a donor could be placed by the account and not forecast. | **`Law(control=None)`**: every unit is an event, read against the mean of the units; the readout "with no record" is one of the variation inside a unit; the rows of the units teach each other. |
| 5f | What a batch has to show was written batch by batch. | **Three systems are forecast before a batch is opened**: the model alone (a readout of control cells), the system as it stood before the last change, and the frozen one. The frozen one has to be calibrated and to have the smallest error. Then the batch's rows become records, and its flags proposals, by rule and by nothing else. |

Not in version 5: the five older proposals for rules of judgment that wait for the owner's decision.

**Batch 3 under version 5** ([`runs/batch3_frozen_20261126`](runs/batch3_frozen_20261126/), 26 November). Three
systems forecast 698 pairs before the values were opened. With no rule changed the error falls from 0.247 (the model
alone) to 0.174 (the system as it forecast batch 2) to 0.138 (version 5); an edge has the size of its pair (error
0.65 on eight own-gene pairs, against 1.19 and 1.98); eight of ten tests of edges proposed beforehand admit the edge.
The intervals of version 5 do not hold: 82% of the pairs are inside, the mean of squared *z* is 3.08. So the gate is
not passed, the batch counts as development, and the final audit waits for another batch.

**Version 6 (27 November 2026): frozen at commit `b4f6196` on the main branch**, the candidate of
[`runs/v6_candidate_20261126`](runs/v6_candidate_20261126/) as it passed the planted twin. It was made from batch 3;
batch 4 is read under it ([`runs/batch4_frozen_20261127`](runs/batch4_frozen_20261127/)).

| | what batch 3 showed | the candidate |
|---|---|---|
| 6a | A readout that is a ratio over a set of columns couples them: an event with an edge to one column moves every other readout by arithmetic. No system forecast that. | **Derived edges.** The declaration of a readout names the columns it is a ratio over. For every event with an edge, the direct part of the edge's readout, measured where recorded and else as forecast, gives what the other readouts of the set move by, with no constant (`events.derived`). A recorded row is freed of that part before it teaches; a forecast carries it (`Law(known=)`). |
| 6b | For the readouts never read no stratum had a recorded pair with no edge and no type, and rule 5a fell back on a library in another condition: a spread of 0.018 where the batch's own rows give 0.03 to 0.05. | **A batch has recorded pairs of no type** in every stratum and on every readout it forecasts, and the spread comes from them. |
| 6c | A few edges missing from the graph move the mean of squared *z* more than all other pairs together. | **The gate asks for the share inside and the median of squared *z*.** A pair far outside is counted as the proposal of an edge. |
| 6d | The law assumes that some pairs with no edge hold a direct part of their own, and gave them an interval of 1.96 spreads of the others. | **The width is measured.** Each recorded pair with no edge is left out and forecast as a new event; the interval is as wide as the score that 5% of them exceed (`Law.conformal`). It needs 40 pairs from 8 events of the kind; else the library's; else the law's mixture, marked as assumed. |
| 6e | The spread of a carried share came from seven pairs and was taken as exact. | **Student's *t*** for the number of pairs measured twice. |

**A change of the method is run on a planted twin before it meets a unit never read.** The twin
(`runs/v6_candidate_20261126/scripts/twin.py`) uses real control cells, so its readouts are measured as the real
ones are, and plants events whose effect is known. Version 5 fails on it as it failed on batch 3 (85% of the pairs
inside); the candidate has 94%. **The units of a batch are chosen by a written rule**: by an item on record, and
the rest by lot.

**Batch 4 under version 6** ([`runs/batch4_frozen_20261127`](runs/batch4_frozen_20261127/), 27 November). Thirty
units never read, chosen by the written rule, 600 pairs, forecasts committed first, no rule changed. 581 pairs
(96.8%) are inside the interval of version 6 and the median of squared *z* is 0.48, where at least 90% and 0.23 to
0.91 were asked: **the gate is passed, the first time for a batch read with no rule changed.** By rule the batch
leaves 600 rows as records, the edge TXNDC17 → CD49f admitted in both conditions (proposed under interferon; the
same share of its scale in all three), DAG1 → PD-L1 refused, and two pairs flagged by the laws fitted again, proposed
and refused on the library's rows. What it does not test: rule 6a on data never read (no unread unit that is not
final removes a large protein) and the correction (the graph holds almost nothing about these units, so the three
systems give the same values). Both are asked in the final audit, which is run only on the owner's word.

**A candidate for the version after: a type is proposed and tested as an edge is**
([`runs/types_by_rule_20261127`](runs/types_by_rule_20261127/), 27 November; exploratory; version 6 is unchanged).
The one type of the system was declared by the author. Rule 7a, not frozen:

| | until now | the candidate |
|---|---|---|
| 7a | A type is declared. | **Flag, propose, test.** An event flagged as moving the stage where the control readout misreads it seeds a type; every event that moves the stage the same way joins (cosine of the shifts at least 0.5; `events.propose_types`). Rows that did not propose it admit the type if, on the pairs it teaches, the forecast with it is closer more often than chance and the squared error falls; else they refuse it (`account.test_types`). A proposed type enters no forecast. |

On a planted twin (16 worlds run once after the rule was settled) every planted type is proposed with exactly its
events, no wrong type is admitted, and a true one is refused in 5 of 32 tests. Told of no type, the real library
proposes exactly the five events declared by hand, and the co-culture's rows admit them (retrospective). Its test on
data never read needs units that move the state.

**The final audit under version 6** ([`runs/final_audit_20261127`](runs/final_audit_20261127/), 27 November, on
the owner's word). 133 final units, 2,564 pairs, read once against control cells no run had read, no rule changed.
The intervals hold (95.3% inside, median of squared *z* 0.41), and on the 160 pairs where the systems can differ the
error falls from 0.336 for the model alone to 0.140: **the audit is passed.** An edge has the size of its pair (5 of
5 inside), and 5 of 8 edges flagged under interferon are admitted. Rule 6a does not meet its own condition (84 of
94 pairs inside, 85 asked): it fails for B2M in the co-culture and holds for the four other units. Version 5 would
have passed as well. What the audit shows holds for one study, one model and knockouts.

**The first training experiment** ([`runs/learned_state_20261127`](runs/learned_state_20261127/), 27 November, on
the owner's word; exploratory; designed in [`runs/learned_state_design_20261127`](runs/learned_state_design_20261127/)).
The pooling alone was trained (2,048 numbers, the encoder frozen) by one constraint: a readout taught by control
cells has to read the recorded events. The loop judged it on cells the pooling never saw. For IRF1 on PD-L1 the
trained state carries 0.66 of the change where the mean carries 0.52 (+0.14 [−0.01, 0.31]: unresolved), the state
typed by hand 0.73 and the input 0.89; with the event labels shuffled the gain is the same (refuted). The model is
still the part the loop has not corrected.

**A check with a known answer** ([`runs/known_answer_20261127`](runs/known_answer_20261127/), 27 November, on the
owner's word; nothing trained; the design and the pass line committed before any value was read). Does a knockout
compose with IFN-γ by addition? For each of 95 knockouts, the departure from addition on 3,303 genes and 20 proteins,
in two halves of the cells; the end of the interval from the mirror side (products of the halves that disagree in
sign), over the knockouts. Refused: exactly the five of the pathway, each with witnesses; the other 90 carry over;
one class, the five, proposed from the witnesses by `events.propose_types`. Passed (P1 to P3). Of the part that does
not add, TEDDY's state carries 0.93 on PD-L1 and its head's output as trained 0.09. Every cell had been read before;
it checks the method and is not confirmatory.

## Commands

```
python -m anm_diag scope SPEC.yaml --out DIR                        the grid, the readings, the typed events
python -m anm_diag account SPEC.yaml --out DIR                      version 4: read, carried, others at every stage; the four parts
python -m anm_diag loop test --hypothesis H.json --arms ARMS.json   a hypothesis against the arms of an intervention
```
