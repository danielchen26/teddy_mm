# The outer loop: design

Status: **proposed** on 20 November 2026. Nothing in this document has been run on cells. It was checked in two ways
that read no cell: a dry run of its routing on the finished round of 20 November, and a back-test on earlier runs
([`runs/outer_loop_design_20261120`](runs/outer_loop_design_20261120/)). What that check changed is listed there.

[`LOOP_PROTOCOL.md`](LOOP_PROTOCOL.md) describes one round: an audit, one change to the system, an audit again. That
is the inner loop. This document describes what happens between rounds when the inner loop has done what it can and
something is still left.

## 1. What is missing

The round of 20 November ([`runs/teg_loop_20261120`](runs/teg_loop_20261120/)) repaired the head under the protocol.
Wherever an event was admitted, the repaired head adds no loss beyond what the embedding carries. The whole chain did
not close, and nothing in the protocol says what to do with what was left:

- intervals that do not decide (PD-L1 in three groups);
- a loss that lies before the head (MHC-I in co-culture);
- responses too small to follow (HLA-E; MHC-I in the IFN-γ condition);
- one event that reads differently from its pool (IRF1);
- rules that could not be evaluated (gate 3, a missing control arm);
- no cell left that had not been read.

None of these is removed by another change to the model. Each asks for a change to something a round holds fixed. The
outer loop is the procedure that decides which fixed thing changes, on what evidence, and when to stop.

## 2. Terms

- **Event**: a typed record of a source. Its fields are time, site, channel, sign, amplitude, stratum, provenance, the
  pool it is merged into, an admission entry and a no-leakage entry (`events.py`).
- **Row**: an admitted event with what was measured for it: the change of the readout, its interval, the displacement at
  each stage, *u* at each stage and scope, and a reading.
- **Schema**: the fields of an event, their allowed values, and the declared pooling (which values are treated as one).
- **Type**: a class of events that the schema treats as interchangeable for a readout. A type is a class of events; a
  field is what generates the classes. The two are kept apart below.
- **TEG** in this repository: the schema, the rows, and the ledger of verdicts joined to them. The paper's typed-event
  graph (main text, Section 4) is another realization, with a prescribed propagation law. What is taken from the paper
  here is the source chain (SI S1.4), the sufficiency rule (SI S1.4.1) and the way a revision is tested (Sections 5.3
  to 5.5 and 6). The propagation law is not used.
- **View**: a stage on the system's own computational path, with a reader specification (classes, their settings, the
  resampling unit). A compression is a stage and is declared.
- **Observer**: the readouts, the continuations (conditions under which they are read) and the tolerance. What a state
  must retain is relative to the observer.
- **Pool**: the real data that exist, by unit (a knockout target, a donor, a study). No experiment is assumed. Each unit
  has a role: unused, fit, development, or final.
- **Package**: what one iteration seals (Section 4).
- **Open item**: a typed record of something an iteration left undecided or unrepaired (Section 5).

## 3. Principles

Each has a source and a consequence for the design.

**P1. Only the two ends are measured outside the system.** The source and the readout are real. A stage has no truth
of its own: it is tested by whether the readout's response to the source can be recovered from it (paper, Section 2).
*Consequence*: a revision is certified against a measured readout, never against the system's own output. A view added
to the input is never the readout or computed from it. The state view, an added view and the readout are three
measurements.

**P2. A stage cannot carry more than the stage before it.** Each stage is computed from the one before.
*Consequences*: (a) the inner loop can bring the output's loss down to what the input view carries, and not below;
what is left at the input view is outside it. (b) If a declared reader of an earlier stage reads worse than a reader of
a later stage, the earlier reader is the weaker one and its reading says nothing about that stage (a *dominated view*).
This is used as a check and never as a shortcut: every stage is measured, none is inferred, because with finite readers
and finite cells the readings need not be ordered.

**P3. A verdict needs resolution.** Support, rejection and unresolved are three outcomes (paper, Section 2; protocol
rule 1). Unresolved is not a loss and not its absence. An interval is wide because the units are few or because
they differ among themselves. *Consequence*: an unresolved item is moved by more independent units of the same type,
or, where the units of a pool differ, by a split of the pool, which is an item of its own with its own test. No change
to the model or the views addresses it. If the units do not exist, the estimate with its interval is the result.

**P4. A distinction is kept or dropped on new histories.** A field may be dropped only when every declared
continuation leaves it invisible to the readout (SI S1.4.1). A finer description is adopted only when, on new
histories fixed in advance, it agrees with the responses *and* the coarser one is violated; when both agree, no benefit
is established and the coarser one stays (Section 5.5). A rejection does not name the missing variable (Introduction;
Section 5.3), and fields that vary together cannot be told apart. *Consequence*: Section 8.3.

**P5. Proposing and certifying are different acts.** Anything may be proposed from development data, by a person or a
program. Only a test fixed in advance, on units not used for the proposal, certifies. Data that were used to choose are
exploratory however often they are read (protocol rule 4). *Consequence*: units are set aside as final before the
first reading and are read once, for the whole claim (Section 11).

**P6. A diagnosis is a hypothesis.** A resolved violation constrains a description; it does not by itself locate a
faulty component or identify a mechanism (paper, Section 9). Outside the declared class of systems a mismatch calls
for competing mechanisms, not a forced label (SI S1.3). *Consequence*: sending an open item to a lever is a claim with
a forecast, a rival and a control, tested like any other hypothesis (Section 10).

**P7. A claim has a domain and an observer, and narrowing either is a result.** The paper abstains where a rule does
not decide (SI S1.3) and separates continued operation from correct operation (Section 3.3). *Consequence*: dropping a
readout, a type or a grade is recorded as a limit with its evidence. It is never reported as closing. A closure holds
for the observer it was read under; when the observer changes, the affected closures are read again.

## 4. An iteration and what it seals

```
the pool of real data
      │  select units, assign roles
      ▼
package t = (E, V, S, C, P), sealed ──────────────────────────────────────────────┐
      │                                                                           │
      │  inner loop (LOOP_PROTOCOL.md): only the model M changes                  │
      │  observe → admit → hypothesize → one change → test and audit again        │
      ▼                                                                           │
every target has one status; everything left is an open item                      │
      │                                                                           │
      │  outer step: check the readings → route each item to a lever →            │
      │  propose → write the test → seal                                          │
      ▼                                                                           │
package t+1 ◀─────────────────────────────────────────────────────────────────────┘
      │
      ⋮   until every item has a terminal state
      ▼
final audit, once, on the units set aside at the start
```

A package holds the five things a sealed declaration fixes. The inner loop changes none of them.

| | sealed thing | holds | changed by the outer step when |
|---|---|---|---|
| **E** | events | the admitted rows; which units of the pool are in which role | an interval does not decide, a response is too small, a type has no row, two fields vary together, nothing is unread |
| **V** | views | the stages that are read, their reader specification, the input view, the readouts | a reading cannot be used, or a resolved loss sits at the input view |
| **S** | schema | the fields of an event, the pooling, how a cell gets its type | a pooled distinction is visible, a kept distinction looks invisible, cells are mistyped |
| **C** | claim | the targets (type × readout), the grade and tolerance asked, the guard, the domain | an event is outside what the input can carry, a limit has been shown, or the guard cannot be held |
| **P** | protocol | the gates, the rules, the chosen numbers | a rule could not be evaluated as written |

The schema is sealed with a version. It does not change inside an iteration. It changes only in an outer step
(Section 8.3 says when a change may be operated and when it counts as adopted).

How deep the inner loop reaches is independent of this design. Today it changes the head. Training the encoder under
the rows of the table is the same inner loop at another depth.

## 5. Open items

At the end of an iteration every target (a type and a readout) has exactly one status:

| status | meaning |
|---|---|
| closed (absolute) | the system's own output is supported, and the control's is rejected |
| the model's part closed (relative) | the model adds no loss beyond the reference stage, and the control does; what the reference stage leaves is an open item |
| not graded | a gate could not be evaluated; the reason is an open item |
| not a target | the response cannot be followed; the reason is an open item |

Everything else that was left is written as an open item, a typed event in the same ledger:

`id · iteration · readout · type · event (pool or site) · quantity · value and interval · resampling units · reading ·
lever · admissible actions · units needed · units in the pool · evidence (the runs it was read from) · state`

An item ends in one of five terminal states:

| terminal state | reached when |
|---|---|
| closed | the quantity is decided on the side that was asked, on final units |
| limit | a loss before the model is shown by an admissible reading (Section 6), and the claim is narrowed by it; the limit names what would reopen it |
| at the margin | the estimate sits so close to the bound that the pool cannot decide it; the estimate and its interval are the result |
| not asked | the readout or type is taken out of the claim by its owner, before a result or with the result shown |
| on the experiment list | the pool has no unit that could decide it; the entry says what would have to be measured and how many units |

"The whole loop closed" is not one of them. The outcome of a campaign is the table of terminal states (Section 12).

## 6. First: can the reading be used

Before an item is sent anywhere, its reading is checked. An item that fails is read `invalid`, and the only action is
to repair the reading.

| check | fails when | repair |
|---|---|---|
| dominated view (P2) | the gap from an earlier stage to a later one is resolved below zero: the later stage reads better | a richer reader of the earlier stage; until then the earlier stage is bounded by the later one and no `upstream` or `unreadable` reading is drawn from it |
| positive control | no reader reaches support for any unit of the family on this scale | the tolerance is declared for this target (protocol v2, "a margin is declared per target"); no `unreadable` reading before that |
| one reader specification | two runs give opposite verdicts for the same event, readout and stage, with different reader settings | one specification, sealed in the package; a bridge reading (Section 9) |
| resampling | fewer units than the chosen number, or draws that are not paired | the interval is marked approximate; a statement about the family waits for more units |
| reader classes | `unreadable` rests on one class | at least two classes, one of them not linear |
| readout control | a readout that should not respond (an isotype control) responds | the readout is repaired or dropped before anything is read from it |

## 7. From a reading to a lever

The inner loop keeps its table (function: change the map; coverage: add cells of the kind, with rows; unit: calibrate at
use). The outer step adds the rest. This table is `ROUTES` in
[`dryrun.py`](runs/outer_loop_design_20261120/scripts/dryrun.py).

| reading | what it means | lever | actions, cheapest first | what not to do |
|---|---|---|---|---|
| invalid | the reading cannot be used (Section 6) | V | repair the reading before any action, as Section 6 says for the check that failed | act on the reading |
| unresolved | the interval does not decide between support and rejection | E | 1. if the units of the pool differ among themselves, the split comes first; 2. more independent units of the same type from the pool; 3. if the pool cannot supply them, report the estimate as it is: at the margin | change the model or the views |
| small response | the response is too small against its interval to be followed | E | 1. more cells per event, or more events that move this readout; 2. otherwise the readout is not asked for this type | make it a target |
| unreadable | a resolved loss at the input view, with the event's own cells, in a view that is not dominated | V | 1. a richer reader of the same input; 2. another input already measured on the same cells, never the readout itself; 3. a new measurement; 4. a stated limit of the claim | train |
| upstream | the input view reads it and the stage the head reads does not | M | 1. change what comes before the stage: an inner loop of another depth; 2. if that part cannot be changed, the item stays open under its owner | change the head |
| no rows | cells of this kind are needed and the table has no admitted event of the type | E | events of the type from the pool | carry rows over from another type |
| pooled, and differs | a distinction the schema pools is visible in the readings | S | split along a field from the record, with rival fields, tested on new units | repair across the pool |
| pooled, may be the same | a distinction the schema keeps looks invisible to the readout | S | merge after bounded support on new units, twice | merge on development readings |
| mistyped | cells get a type from a reader of their state, and some get the wrong one | S | take the type from the record | tune the reader |
| self event | the event acts on the readout's own gene: an RNA-only system may have nothing to read (rule S3) | C | 1. a stated limit for this type; 2. or the assignment as an input, which is another system: the owner's decision | use it as a row |
| rule | a gate or a test could not be evaluated as written | P | amend the protocol before the next iteration | re-grade the finished run |
| selection | an event was selected by a measurement on the cells that were then audited | P | a selection split | keep selecting on audited cells |
| no final data | no unit of the claim is unread | E | a census of the pool; units set aside as final before anything else is read | call a development reading a confirmation |
| confounded | two fields vary together in the units that were read | E | 1. units of the pool in which they vary apart; 2. otherwise the claim names their product | attribute the difference to one of them |
| guard conflict | every change that removes the loss breaks the guard on the system's own task | C | 1. a change that holds the guard: the inner loop again; 2. if none does, the owner decides between the target and the guard | widen the guard after the result |

A row of the table therefore says which kind of action comes next: its reading, with its type, selects the lever. It
does not say which measurement to add, which units to take or which field to split on. Those are proposals, and
Section 8 says how each is made and tested.

## 8. The five levers

### 8.1 E: which units are read next

The pool is finite and real. Selecting from it is the active part of the loop.

0. **Census.** Before the first outer step, every unit of the pool is listed with its type (from the record), its
   cells, whether the system's states exist for it, its role and how often it has been read. A census reads no
   response.
1. **Final units first.** A share of the units of the claim is drawn at random and named under `roles.final`. The
   draw is stratified by every field of the record, not only by the types the schema has today, so that a later split
   still has final units on both sides. A run that touches them is refused unless it is the final audit (protocol v2).
2. **Candidates for an item** are units of the item's type. For a question about a field, they are units that vary
   that field while the others are held (Section 8.3).
3. **A selection split.** A stated share of each candidate's cells, and of the control cells, is used to screen it
   (the response against its interval, the no-leakage rules). Those cells never enter a reading of *u*. An event is
   never selected by the system's error on it. The split costs the audit that share of the cells; this is the paper's
   preparation step, which keeps the endpoint out of what admits an event (SI S1.4).
4. **How many.** For an unresolved item with *n* units, estimate *q* and bound *b*, the units needed are about
   *n*·(half-width / |*q* − *b*|)², if the estimate holds and the half-width falls as 1/√*n*. This is a planning
   number with its own forecast record. Where it exceeds what the pool holds, or four times *n*, the item is at the
   margin.
5. **A calibration sample.** A stated share of every batch is drawn at random from the claim's domain, whatever the
   open items ask for. Units chosen for an open item answer whether that item can be decided. Units drawn at random
   answer how often the system fails in the domain. A statement about the domain rests on the second kind and on the
   final units.
6. **Roles before reading.** The batch is split by unit into fit and development, and written as a sealed list of
   unit identifiers, before states are computed or anything is read. Control cells are assigned to roles like any
   others: a fit that uses them takes resolution from the audit (protocol rule 3).
7. **What the pool cannot supply** goes on the experiment list with the number of units it would take.

### 8.2 V: what is read, and with what

- A view lies on the system's own path. An approximation of a stage is not that stage.
- The reader specification is part of the package. Readings of two iterations are compared under one specification.
- `unreadable` is drawn only at the input view, with the event's own cells, in a view that is not dominated, by at
  least two reader classes of which one is not linear, and where the positive control holds (Section 6).
- The ladder, cheapest first: (0) a richer reader of the same input, which needs no new data; (1) another input
  already measured on the same cells; (2) a new measurement; (3) a limit of the claim.
- Independence: an added input is not the readout and is not computed from it. Where the added input is the
  assignment of the event itself, the system has become one that is told the perturbation. That is another system, and
  the change is the owner's decision.
- A baseline reader that does not use the model (protocol gate 6) stays in every package. Where it reads as well as
  the model's stage, the model adds nothing for that target, and that is reported.

### 8.3 S: the schema

The paper's own procedure for revising a description (Sections 5.3 to 5.5) is applied to the fields of an event.

**Signals.** A split is proposed when a pool hides a difference: an event's interval lies outside its pool's;
`events.pooling_evidence` reads `differs`; or a repair or a row carried across the values of a field adds a resolved
loss while within a value it does not. A merge is proposed when rows carried across the values of a field add no loss.

**Candidates.** A rejection does not name the missing variable. For every split the proposal lists rival fields from
the record that differ between the events in conflict. A field that varies together with another in the pool cannot be
separated from it: it is named as their product (cell type × study), and the claim says so (`events.check_contrast`,
C3).

**Test of a split**, in the next iteration, on units not used for the proposal, fixed in advance: the finer schema is
adopted if, within each new type, the readings agree with the type's rows (a repair made with the rows of the type adds
no loss), and the coarser schema is violated on the same units (the repair made with the pooled rows adds a resolved
loss). If both agree, the coarser schema stays. If the finer one fails too, the missing distinction is another one, and
the rivals are tested.

**Test of a merge.** Bounded support on new units in which the field varies and the others are held: the interval of
the difference inside the step margin for every declared continuation (readout, stage, arm). Twice, on two sets of
units. A merge asks for more than a split because it adds a constraint the system must then satisfy, and because five
events cannot show that nothing differs.

**Operated and adopted.** A split may be operated at once, marked provisional: reading each side separately never
widens a claim, it only leaves fewer units on each side. It counts as adopted, and the version is bumped, when its test
has passed. A merge is not operated before its test has passed twice. So iteration *t*+1 runs under a sealed schema
that may carry provisional splits, with the coarser schema read beside it as the other arm of the test.

**What may be merged depends on what is carried.** The round of 20 November is an example, on development cells:
across the two culture conditions of one cell type, the rows carried and the baseline did not. A pooling declaration
names the quantity it pools.

**Relations: the edges.** A candidate field may be a relation taken from outside knowledge, such as a pathway graph:
"acts on the readout's own gene", "one step before it", "in the same complex". Such a field is a proposal like any
other and is tested in the same way. A knowledge graph is a source of candidates and never of verdicts. The edges of
the table itself are the transfers that have been tested: rows of one type used for another, with the outcome. They
change only as a merge or a split changes.

**Ties.** With a tolerance, A may be indistinguishable from B and B from C while A and C differ. Then none of the three
is merged.

**Type from the record.** A cell's type is read from the design labels. Where no record exists at use, the reader that
assigns types is a part of the system and is audited like one.

**No claim of a minimum.** What can be said is: the coarsest schema that has passed the tests run so far, for these
readouts and continuations.

### 8.4 C: the claim

Actions: take a type out of the domain; ask a lower grade; declare a tolerance for a target; drop a readout from the
score. Each is a `limit` event with the reading it rests on and what would reopen it.

Two rules keep this lever from closing a loop by shrinking it. A limit needs an admissible reading (Section 6): an
unresolved item cannot become a limit. And every report carries the coverage of the claim: of the targets declared at
the start, how many are in each terminal state.

A change of the observer (a readout, a continuation, a tolerance) voids the closures it touches until they are read
again (P7).

### 8.5 P: the protocol

An amendment names the case that could not be evaluated, applies from the next iteration, and changes no outcome
recorded under the rule a run was made with (protocol, "Evidence grades"). The round of 20 November asks for five:

| case | amendment proposed |
|---|---|
| gate 3 asks a repair to come inside the margin on its own cells, which is impossible where the reference stage itself leaves a third | state gate 3 against the reference stage: no added loss on the cells it was fitted on |
| the hypothesis had no control arm, because the other group had no admitted event | name a control that always exists: the same update with rows of events that do not move the readout |
| a group ready as a pool with no admitted event | not a target (as resolved in the round) |
| the response of an event was screened on the cells that were then audited | the selection split of Section 8.1 |
| the type of a cell was guessed from its state | the type comes from the record (Section 8.3) |

## 9. The order of an outer step

One item, one lever per iteration, so that a changed reading can be attributed. Changes that touch every item come
first and are separated from new data by a bridge reading.

1. Check the readings (Section 6).
2. Amend the protocol (P). Write the schema proposals (S) with their tests; provisional splits and the type from the
   record take effect here, merges do not.
3. **Bridge reading.** The cells and the model of iteration *t* are read again under the new rules, the schema as it
   will be operated and the new reader specification. No new data, no change to the model. It shows what the change of
   rules alone does. It is exploratory.
4. Richer readers (V, step 0), which need no new data.
5. Select units (E) and, where the owner has decided it, add a view (V, steps 1 and 2).
6. Write the hypotheses and seal package *t*+1.
7. The inner loop runs under package *t*+1.

## 10. The next package

Every action of an outer step is written as a hypothesis before it is taken, in the form of `loop.hypothesis_event`:

- the open item it answers;
- the claim: resolution, view, split, merge, rows, or limit;
- a forecast interval for the deciding quantity after the action;
- a rule, stated against a reference read on the same units (a bound carried over from other units can be out of
  reach there: run 22). Where the reference shows nothing to remove, the outcome is "not testable", not "refuted";
- a rival: the next lever that could explain the item;
- a control, by the kind of claim:

| claim | control | if it fails |
|---|---|---|
| resolution (more units of the type) | the new units read inside the interval of the old ones | the type hides a distinction: the rival is a split |
| rows (events of a type that had none) | the three arms of the protocol: the rows; the same number of cells of a kind already covered; another map | the reading was not coverage |
| view (a richer reader, an added input) | a view of the same size with its rows permuted across cells gains nothing; the added view alone, without the model's stage | the gain is not information about the event |
| split | a random split of the same sizes does not pass the same test | the test passes anything |
| merge | a field known to matter fails the same test on the same number of units | the test could not have seen a difference |
| limit | the positive control of Section 6 | the scale, not the system, is at fault |

Ledger rules for the outer step, beside H1 to H5:

- **O1** every target has exactly one status at the end of an iteration;
- **O2** every open item has a reading that passed Section 6, a lever and an owner, or a terminal state with its evidence;
- **O3** every difference between package *t* and package *t*+1 names the open items it answers;
- **O4** every action has a hypothesis written before it, and its test uses no unit the proposal was read from;
- **O5** no unit under `roles.final` is touched, and the number of times each development unit has been read is recorded;
- **O6** a change to the claim is a limit event, and the coverage table is computed again.

## 11. Final data and certification

- The final units are fixed at the start, by unit, at random, stratified by the fields of the record (Section 8.1).
- One final audit for the whole claim: the model and the package frozen, a sealed declaration, run once, at least
  three seeds where anything was trained, the outcome required in every seed (protocol gate 5).
- The comparisons of the final audit are declared as a family, and the error is allocated across them, as the paper
  does for its simultaneous bounds (Sections 5.3 and 6.2).
- The guard on the system's own task and the baseline comparison are read on final units too. A closed loss is a
  statement about fidelity to the measured response. It is not a statement that the system is good at its task.
- Until then every status is exploratory, and says so.

## 12. Stopping

A campaign stops when no open item has an admissible action in the pool; or when the development units have carried
the number of iterations stated for them (protocol rule 6); or, for one item, when two hypotheses in a row about it
were refuted.

It terminates: for every item it acts on, an outer step either uses a unit that was unused, or adds a view that was
not read, or takes an action allowed once per item (an amendment, a richer reader). The pool and the measured inputs
are finite. A schema decision is reversed only on new units.

The outcome is reported as one table: every target declared at the start, its status, and for every open item its
terminal state. The campaign has closed the loop only if every target is closed (absolute) on final units. Otherwise
it has closed some targets and stated, for each of the others, a limit, a value at the margin, or what would have to be
measured.

## 13. Who does what

| act | by | on |
|---|---|---|
| propose a reading, a route, a field, a batch of units | anyone: a person, a program, a model | development data |
| check the ledger, admit, grade, bump the schema version | a deterministic program (`loop.ledger_check`, the sealed declaration) | the rules of this document |
| decide a new measurement, a system that is told the perturbation, a limit of the claim, opening final units | the owner | the open item and its evidence |
| certify | the final audit | final units, once |

Later, the proposing roles can be separate agents that communicate only through the ledger. The certifying role stays
a program.

## 14. What comes from where

| | from the paper | from this repository | new here |
|---|---|---|---|
| typed events, the source chain, admission | SI S1.4 | `events.py` | |
| when a field may be dropped | SI S1.4.1 | `events.closure`, `pooling_evidence` | the tests of a split and of a merge as an operating rule (8.3) |
| three outcomes, tolerance fixed first | Section 2 | `scope.py`, `loop.py` | |
| a revision is tested on new histories, with rivals | Sections 5.3 to 5.5, 6 | | applied to the schema, the views and the units (8, 10) |
| a violation does not name the missing variable or the faulty component | Introduction, Section 9 | the diagnostic loop | a route is a hypothesis with a rival and a control (P6) |
| abstaining; fidelity is not task correctness | SI S1.3, Section 3.3 | the guard (gate 6) | limits and the coverage of the claim (8.4) |
| a change to the system | named as an intervention on the system, with an analytic example and no measured repair (Section 2, SI S1.3) | `LOOP_PROTOCOL.md` | |
| readings to actions | | protocol v2, six readings | fifteen readings for the outer step, five levers (7) |
| which data are read next | | `roles.final` | the census, the selection split, the calibration sample, units needed (8.1) |
| when to stop | | protocol rule 6 | terminal states and the outcome table (5, 12) |

## 15. Chosen numbers

| number | value |
|---|---|
| at the margin | more units needed than the pool holds, or more than four times the present number |
| confirmations for a merge | two sets of new units |
| reader classes for `unreadable` | two, one of them not linear |
| calibration sample | one third of a batch |
| selection split | one third of a unit's cells, and of the control cells |
| refuted hypotheses in a row before an item stops | two |

None of these is derived. Each is reported with every use, as in the protocol.

## 16. What this design does not do

- It does not say which new measurement to make. It says when one is the only action left.
- It does not remove confounding. Units of the pool differ in more than the field in question; the design asks for the
  confounded fields to be named and for crossed units where the pool has them.
- It does not make a schema minimal, or unique.
- It has no test yet for the outer rules O1 to O6, and no code beyond the dry run. `loop.hypothesis_event` does not yet
  accept the outer claims.
- It has one observer per package. A range of observers is several claims, and their joint treatment is not designed.
- Four of its fifteen readings (`upstream`, `no rows`, `confounded`, `guard conflict`) did not occur in the round it was
  dry-run on. The last two are taken from earlier runs (22 and 16).
- It has not been run. Whether it converges in a useful number of iterations on this pool is not known.
