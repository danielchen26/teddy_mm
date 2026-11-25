# The loop for any system

Written on 22 November 2026, after the first outer step. Run and checks:
[`runs/framework_generalization_20261122`](runs/framework_generalization_20261122/).
Updated on 23 November 2026: events without control cells are in the round
([`runs/unit_events_20261123`](runs/unit_events_20261123/)), the head was given a view of the input beside Z
([`runs/input_beside_z_20261123`](runs/input_beside_z_20261123/)), and the round on batch 1 was run
([`runs/batch1_round_20261123`](runs/batch1_round_20261123/)). Sections 2 to 4, 6 and 7 carry the changes.
Updated on 24 November 2026: the stages inside the encoder were read for one stratum
([`runs/encoder_stages_20261124`](runs/encoder_stages_20261124/)), and the encoder was run with its annotation
positions ([`runs/annotation_input_20261124`](runs/annotation_input_20261124/)), and the head was given the first
of them beside Z ([`runs/annotation_beside_z_20261124`](runs/annotation_beside_z_20261124/)); Sections 2 and 7.

The loop was built on one system: TEDDY with a prediction head, RNA as its input, proteins as readouts, knockouts as
events. This document says what in it does not depend on that system, what a system has to expose to be put through
it, how that was checked, and what is still tied to the first system. It ends with the plan.

Three documents describe the loop itself and are not repeated here: [`LOOP_PROTOCOL.md`](LOOP_PROTOCOL.md) (one
round: the inner loop), [`OUTER_LOOP.md`](OUTER_LOOP.md) (what happens between rounds), and the ANM paper
(arXiv 2610.14549) for the test.

## 1. Three parts, and where each is

| part | what it does | code | names a system? |
|---|---|---|---|
| the audit | for an event, a readout and a stage: the share of the response that a reader of the stage leaves unexplained (u), by scope | `anm_diag/scope.py` (`run_spec`) | no |
| one round | admits events on a share of their cells, fits a repair on one half and reads it on the other, compares every arm with the stage the model reads, places what is lost before the model | `anm_diag/cycle.py` | no |
| between rounds | turns the state of the targets into typed open items, checks each reading, sends each item to one of five levers, checks the ledger | `anm_diag/outer.py` | no |

A test (`tests/test_cycle.py`) fails if a name of the first system appears in `cycle.py`, `outer.py` or the planted
system. Everything a run knows about its system sits in the run's own scripts: how the cells are loaded, what the
stages are called, which repair is used.

## 2. What a system has to expose

| | what the loop needs | first system (TEDDY) | a model with several modalities | planted system (`anm_diag/toy.py`) |
|---|---|---|---|---|
| cells with a record | for every cell: the event it belongs to, whether it is a control cell, and its type, all from the record of the experiment and none from a reader of the cell's state | guide assignment; cell line and condition | the same, whatever the assay | `unit`, `treated`, type one or two |
| a readout | measured outside the system, never an input, and not used to prepare an input | surface protein (CITE-seq) | any measurement the model does not receive | `y` |
| the system's own output | what the system as trained says for that readout, on a scale that moves with it | the head's protein output | the output of the model for that readout | `o` |
| the stage the model reads | an array per cell: the last thing before the part of the model that a round may change | the embedding Z | the fused embedding | `z` |
| views of the model's own input | one or more arrays per cell, read by at least two reader classes of which one is not linear | 512 components of RNA; 4,000 genes | each modality the model receives, alone and together | `view_a` |
| stages inside the model | optional: arrays per cell at points between the input and the stage the model reads, in order | for one stratum (24 November): the tokens, what enters layer 1 and the output of each of the twelve layers, each averaged over a cell's tokens; the four annotation positions | the output of each modality's encoder; the fused state | `mid` |
| other measured inputs | arrays the model does not receive, never the readout itself | none was read | a modality that was measured and left out | `view_b` |
| a control view | a view that must read nothing | the RNA view with its rows permuted | the same | `perm` |
| events with control cells | units that change the readout, each with cells that did not get the event | knockouts against non-targeting guides | perturbations of any kind | planted shifts |
| events without control cells | a unit read against the other units (a donor, a batch, a study), or a group of units against another (a disease): the unit of every cell, from the record | eight donors of a study the head was not fitted on | the same | `run_units`, `run_groups` |
| a repair | an object fitted on control cells and on the rows of the admitted events, with one output per arm (§3) | a linear readout of Z, updated once by the rows | a change to whatever part the round may change | the same linear readout |
| a census | which units of the claim exist, which were read, which are set aside as final | 615 units, 133 final | the same | not used |

Two things are the system owner's and stay outside the code: the claim (which readouts, which types, what tolerance)
and which part of the model a round may change.

## 3. The calls

    sel = cycle.select(ids, unit, y, seed)                      # A1 and A2 on a third of the cells; those cells are read no further
    grid = cycle.grid(scope.run_spec(spec, folder))             # the audit, or cycle.load(folder) for a grid on disk
    T = cycle.one_type(readout, group, type, x, y, o, treated, unit, sites, *grid,
                       model_stage, inputs=[...], between=[...], others=[...], control=..., alone=[...])
    T = cycle.one_units(readout, group, type, x, y, o, unit, events, *grid,       # events without control cells
                        model_stage, inputs=[...], control=...)
    items = outer.open_items([T["target"], ...], world)

`one_type` is `cross_fit` (the repair, every cell read by a fit that did not see it), `read` (A3, every arm against the
reference stage, gate 3, the two controls), `place` (where a loss before the model is) and `target`. Where a run gives
stages inside the model, `place` reads the loss step by step along them and an `upstream` item says between which two
stages the response stops being readable. That is a statement about two recorded stages and the declared readers, not
about which weights are at fault. A type with fewer
than three admitted events is read event by event and resampled by cell (`target_by_event`): two events say nothing
about a family. `anm_diag/toy.py` is a complete run in about a hundred lines and is the pattern for a new system.

A repair is a plug-in: any class with the constructor `(x, y, o, control_cells, [(event, cells)], other)`, a dict `th`
with one entry per arm (`controls` and `rows` at least) and a method `out(arm, x)`. `cycle.LinearRows` is the one used
so far. A repair of an encoder is another class of the same shape and has not been written.

**Events without control cells** (since 23 November; [`runs/unit_events_20261123`](runs/unit_events_20261123/)). A
unit read against the other units goes through `one_units`: `select(control=None)` screens every unit against the
others, which count as units; the repair is `cycle.UnitRows`, the same linear readout where no cell is a control
cell (it starts from the variation among the cells of a unit and is updated once by the row of every unit); every
event is read on its own, its cells by cell and the other units by unit. Beside the arms of a knockout it has the
arm `unseen`: the same repair fitted without any cell of the unit, read against a reader that never saw the unit.
Where the unit's own row reads at least a step more than `unseen`, the outer step raises "pooled, and differs" for
that unit. A group of units against another (a disease) is `one_type(control_cluster=unit)`: the patient donors are
the events of a pool, the healthy donors' cells are the control cells, and both sides are resampled by unit.

## 4. How it was checked

**A planted system, where the place of every loss is known.** Two modalities, two steps inside the model, a head, and
two types of cell. Planted in it: a loss in the head; a loss in the first step and one in the second; a response that
only the second modality carries; an event that acts on the readout itself; an event too weak to follow; a second type
whose repair is another one; and one event that sits in a pool it does not belong to. A second version of the same
system has nothing wrong. Every run is scored by a rule written before the runs: right, undecided (the loop says
less), or wrong (the loop says something else that is definite). Twelve seeds at three sizes, 504 scores:

| size | control cells per type | cells per event | right | undecided | wrong |
|---|---|---|---|---|---|
| large | 4,000 | 800 | 167 | 1 | 0 |
| medium | 1,500 | 300 | 158 | 10 | 0 |
| small | 500 | 100 | 140 | 27 | 1 |

Where a loss is (upstream or unreadable, in which step of the model, whether another input reads it) is found at
every size: 287 of 288 right, one undecided. What shrinks with the cells is the absolute statement: "closed" and
"nothing is open" need an interval inside the margin, and with few cells the loop says "undecided" instead. The one
wrong score is a second item for an event whose loss was placed correctly: with 67 cells a reader of the input left
0.16 of the response, and the loop read that as a loss at the input.

**The first system, computed again.** The general code was run on the cells, halves, seeds and grids of the first
outer step and compared with what the step's own scripts recorded:

2,732 values in eleven comparisons (six readings of the bridge, two events of the round under package 2, three
placings), largest difference 0; the same statuses for the seven targets; the same 17 open items.

Both checks, with their scripts and result files, are in the run folder.

**Events without control cells, planted** ([`runs/unit_events_20261123`](runs/unit_events_20261123/), 23 November).
A study of units in which no cell is a control cell, with nine units that differ in a known place (three in the
head, one in each step of the model, one that only the second modality carries, one that acts on the readout, one
that only its own row reaches, one too weak to follow), and a disease, patients against healthy donors. Each also
with nothing wrong. The rule was written before the scored runs. Twelve seeds at three sizes, 576 scores:

| size | ordinary units, cells per unit | right | undecided | wrong |
|---|---|---|---|---|
| large | 72, 300 | 192 | 0 | 0 |
| medium | 40, 200 | 191 | 1 | 0 |
| small | 20, 100 | 112 | 80 | 0 |

With 20 units a loss before the model is placed in 3 of 12 runs; in the others no reading is supported and the loop
says `invalid`. Nothing false is said at any size.

## 5. What the check changed

The planted system exposed six places where the loop said something the plant did not support. None had shown on
the first system, and the first step's recorded items and statuses are the same after the fixes (they were
regenerated and compared).

| | what was wrong | what the rule is now |
|---|---|---|
| 1 | a loss that sits both at the input and between the input and the stage the model reads was reported as one of the two | both are reported |
| 2 | an input view was checked, and could raise an item, where the stage the model reads already carries the response | where the reference stage is supported the input view is not consulted |
| 3 | an event read on its own was called "pooled, and differs" | only an event inside a pool can be; and an event inside a pool is flagged when it is rejected at the reference stage while others of the pool are supported there, even if no repair was read |
| 4 | a reader of a stage that failed for an event raised an item although the system's own output reads that event | an output inside the margin has lost nothing on the way: no earlier stage is consulted for it |
| 5 | "another measured input reads it" asked for support at that input | it is also met when that input reads at least a step more than the model's own input, on paired resamples |
| 6 | a system with nothing wrong, read on fewer cells, was called "open before the model" while every reading was only undecided | a loss before the model is said only where a stage before it is rejected; otherwise the status is "undecided before the model" (a new, eighth status) |

Four further things were missing rather than wrong. The general code had no A3 (does the output as trained lose a
step more than the reference stage?), so a repair could have been graded where the model had nothing to remove; it is
now part of `read`, and a target without it is described by the output as trained. The reader of a view was chosen by
the smallest |u| while the audit tool uses another rule (supported if any declared reader supports, rejected only if
all reject); `place` now uses the tool's rule, and on every reading the first step recorded the two rules give the
same reader and verdict. `place` read one gap, from the input to the stage the head reads; it now reads the loss step
by step along any stages inside the model that a run records. And a type with fewer than three admitted events had no
path of its own; it is read event by event.

Two things were seen and not changed, because each alters how a reading of real cells would be used: the check for a
dominated view can fire between two stages that both read nothing, and a reader that fails on few cells is read as
an input that lacks the response. Both are in the run folder with a proposal.

## 6. What is not general yet

- **Events without control cells of their own: in the round, and thin on real cells.** Since 23 November the round
  reads a unit against the other units and a group of units against another (Section 3), and both are checked on
  planted systems (Section 4). On real cells there is one round, on eight donors that had been read before, in which
  donor and batch vary together. A disease has not been read on real cells: the one disease study at hand is a final
  study. A unit event is read with linear readers only when a panel of readouts is read at once. Nothing is said
  about a family of units from fewer than ten of them.
- **A repair is made for what it has seen.** Across kinds of cell the rows of the events carried and the readout of
  control cells did not ([`runs/outer_step_20261122`](runs/outer_step_20261122/)). Across donors the row of a unit
  does not carry to another unit: for 10 of 16 pairs of a donor and a protein the repair needs the donor's own row
  ([`runs/unit_events_20261123`](runs/unit_events_20261123/)). No repair so far holds for a type or a unit it has no
  cells of.
- **Every view on the same cells.** Gaps between stages are read on paired resamples, so every stage has to exist for
  every cell that is read. A modality that is missing for part of the cells is not handled.
- **Events are yes or no, and a type has one partner.** There is no dose. What carries across types is tested against
  one other type per run.
- **One quantity.** The shift of the mean of a readout. The audit also reads the variation among the cells of a unit;
  the round does not use it.
- **One repair.** A linear readout of a stage. On 23 November the same readout was given a wider stage, the state the
  head reads beside a view of the input ([`runs/input_beside_z_20261123`](runs/input_beside_z_20261123/)): one item
  gained a step and did not close, and one control target lost its grade. Nothing that changes an encoder has been
  put through the loop, on any system.
- **A limit of u.** The slope that puts an output on the readout's scale is taken among control cells. A response
  along a direction that does not vary among control cells cannot be read as supported by any output. This is why the
  planted system has no case in which the rows of events add a step over a readout of control cells; that case has
  only been met on real cells.
- **Readings never met.** Of the fifteen readings of the routing table, three have been produced only by constructed
  targets in the unit tests: no rows, confounded, guard conflict. One more (pooled, may be the same) was met on real
  cells and is not planted.
- **One real system.** That the code names no system, places planted losses and reproduces the first system is not the
  same as a second real system passing through it.
- **What the paper has and the loop does not use.** The propagation law of the TEG, and the TEG as a constraint in
  training. Here the TEG is the test and the source of rows for a repair of the head.
  A first trial of the TEG as a constraint in training was made on 24 November, on one layer
  ([`runs/last_layer_trial_20261124`](runs/last_layer_trial_20261124/)); what it says is in Section 7. One of its
  arms gave IRF1 a row derived by the simplest propagation rule, the mean of its type. Such a row can be trained
  in, and for IRF1 its value is wrong by a factor of 2.7: a rule has to be tested on the record before it is given
  to a training.
- **The schema is not revised by the loop yet.** The rules of a split and of a merge are written
  ([`OUTER_LOOP.md`](OUTER_LOOP.md), Section 8.3), and the loop flags an event that differs from its type: IRF1 in
  THP-1, and ten items for eight donors read as one type. What follows a flag has not been done. Every new type so
  far was proposed by hand from a field of the record; the type made for IRF1 has one member, which is the same as
  the event having its own record, and it did not help (0.68 with its own row, 0.54 with the pooled rows); and the
  test of a split on units not used for the proposal could not be run. This is why the readings so far look like
  "readable with a record, not without": a type that covers an event without a record of its own needs several
  members, and no such type has been proposed and tested.
  A first test was made on 24 November, for the donors ([`runs/type_growth_20261124`](runs/type_growth_20261124/)):
  a type as a grouping of the eight donors into two groups of four, the candidate from the record (the batch) and
  the one found in the states with no protein (the same grouping), against all 35 groupings. The batch is the best
  of them and does not cover a member: the type-mates' record leaves 0.82 of a donor's change, the repair with it
  gains no step over the repair without (0.77 against 0.75; 0.45 with the donor's own row), and the own row still
  adds a step for 10 of 10 flagged pairs. So for these donors the items do not close by the schema. Not tried: a
  type of another shape, a dimension that is not a grouping, and a type for knockouts.

## 7. Plan

| | step | needs | whose |
|---|---|---|---|
| done | the general cycle, the planted system, the comparison on the first system | CPU | |
| done, 23 November | events without control cells (donor, disease) in the round: [`runs/unit_events_20261123`](runs/unit_events_20261123/). The general code, planted units and a planted disease (576 scores, none wrong), and one round on eight real donors: 7 of 16 pairs relative, none absolute; for 10 of 16 the repair needs the donor's own row. A disease on real cells is not read | CPU | |
| 2 | a second system by swapping the encoder on the same cells: the loop is run unchanged, only the stage arrays differ. The lightest form was read on 24 November, two larger checkpoints of the same family ([`runs/model_sizes_20261124`](runs/model_sizes_20261124/)): not an independent system, and the step stays open | inference | the owner names the encoder |
| done, 23 November | a repair in front of the head at small scale: the head reads Z together with a view of the input ([`runs/input_beside_z_20261123`](runs/input_beside_z_20261123/)). For IRF1 the repaired output goes from 0.68 to 0.39, a step, and does not close; given to every type it costs one control target its grade | CPU | |
| done, 23 November | the states of batch 1 (6,759 cells) and the round on them: [`runs/batch1_round_20261123`](runs/batch1_round_20261123/). Three events admitted, none gives the model anything to repair; the repair itself is not tested by this batch | inference on the owner's machine | |
| 5 | **retraining TEDDY with the TEG as a constraint in pre-training.** The rows of the TEG enter the training loss; the schema is sealed for an iteration; the retrained encoder is one more system for the same loop, judged by the same items | several GPUs, the pre-training data | the owner; after 1 to 3 |
| 6 | the final audit on the 133 units set aside, once and with the hypotheses named | inference | the owner |

Step 5 was placed after 1 to 3 on purpose: they cost no GPU, and 3 says whether the upstream loss is worth a
retraining at all. Steps 1, 3 and 4 are done. What step 3 says: the part that the encoder drops for IRF1 is in the
input and comes back in front of the head, but a repair made from one record does not close it, and without the
event's own record the wider stage gives nothing. What steps 1 and 3 say together: no repair made after training
holds for a type, a unit or an event it has no cells of. That is what a retraining has to change, and it is the
measure step 5 will be read by (the arms `controls` and `unseen`).

Every step that is left needs the owner: the encoder for step 2, GPUs and the pre-training data for step 5, the
agreement to read for step 6.

On 24 November the loss of IRF1 was read along the stages inside the encoder
([`runs/encoder_stages_20261124`](runs/encoder_stages_20261124/)). The tokens the model receives still carry the
response (0.21 against 0.18 at the input view). The average of what enters layer 1 does not (0.52), the averages of
layers 1 to 11 do not, and the last layer brings part back (Z: 0.40). The tokeniser is not where it is lost, and
eight window means read no more than one mean. For step 5 this says what a retraining has to change, and where to
start: the response has to reach the average that Z is, and the last layer is the one that moves it there. A trial
that trains the last layer on the owner's machine is within reach (measured: 23 cells per second) and comes before
any retraining on GPUs.

On the same day the encoder was run with the part of its input that this repository had never given it: four
annotation positions in front of the genes ([`runs/annotation_input_20261124`](runs/annotation_input_20261124/)).
With the placeholders or with the labels of the cell line the state reads IRF1 as it does without them (0.41 against
0.40, no change of a step), so the loss is not a matter of a missing input and nothing has to be recomputed. The
annotation positions themselves, a summary the model learned in pre-training, read more than the mean (0.29 against
0.40); by the rule that difference is unresolved. It is a lead and comes before the trial: let the head read the
first annotation position beside Z on the three strata, which needs no training. From this run on the grids are made
on the owner's machine, where they equal the cloud's to 5e-6 and take less than half the time.

That was done the same day ([`runs/annotation_beside_z_20261124`](runs/annotation_beside_z_20261124/)), with the
five targets, cells, halves and resamples of the run of 23 November. For IRF1 the repaired output goes from 0.68 to
0.51, a step; it does not close (the input view beside Z gave 0.39). MHC-I in co-culture gains no step. No control
loses and every grade is kept; PD-L1 under IFN-γ goes from relative to absolute, where the input view had cost it its
grade. So the annotation position is the first wider stage that could be given to every type: a candidate for an
adjustment between rounds, to be declared and then read on units that have not been read. What it leaves for a
retraining is the part no stage after training has reached: an event read without its own record (0.82 here).

The same reading was then made on the two larger checkpoints of the official release, 160M and 400M
([`runs/model_sizes_20261124`](runs/model_sizes_20261124/)), with the same cells and tokens. Their states leave 0.35
and 0.35 where the small one leaves 0.40; against the small model that is −0.04 [−0.11, 0.04] and −0.05
[−0.12, 0.04], unresolved, and the two larger ones do not differ. From the input they leave 0.17 more
([0.08, 0.27]), unresolved where the small model's 0.21 is a resolved loss. What enters the first layer is a loss in
every size. So a larger checkpoint is not the repair, and nothing is a reason to make the trial on one. Inside, one
thing differs: in the larger models the layer before the last already reads what the last reads, and the last layer
does not bring back a step as it does in the small one.

A trial at the smallest scale was then made, on the owner's laptop
([`runs/last_layer_trial_20261124`](runs/last_layer_trial_20261124/)): only the last layer of the 70M checkpoint is
trained, the loss is the rows of the TEG under a readout fixed beforehand, and every cell is read with a state from
a layer that did not see it. It is not step 5, which trains in pre-training; it asks what step 5 is to be judged
by. What the encoder loses does not come back: the reference reader stays at 0.40 (H2 refuted). With IRF1's row in
the training the readout with no record of it reads a step more (0.90 to 0.62); with no cell of IRF1 in the
training, 0.12 more, unresolved (H1). The gain belongs to the vector the layer was trained under: that readout kept
exactly reads twice the gain, a readout fitted on other control cells a third to a half, and shuffling the event
labels changes nothing. Three things follow for a retraining. The last layer alone is not where the loss can be
undone. A loss on one projection teaches one projection: if the response is wanted in the state, the readout has to
be fitted again while the encoder trains, or be several; if it is wanted at the system's output, the readout has to
be the system's head, which this trial did not use. And the reading has to use a readout fitted on cells the
trained layers have not used.

The owner's design of a model built on the TEG keeps the state by type and lets every reader read the state only.
It was brought to TEDDY at the one place where it needs no training
([`runs/typed_state_20261125`](runs/typed_state_20261125/)): the output of the last layer is kept in two blocks, the
128 genes that the pathway's events move most and the rest, which add up to the mean a head reads today. By the rule
written first, a reader with IRF1's cells leaves 0.27 at the typed state where the mean leaves 0.40 and the input
0.18 (a gain of 0.12 [0.08, 0.18], unresolved), 128 genes at random give nothing, and the round did not use the
typed state (refuted). A check added after that reading found the reason in the round: its linear reader does not
bring the columns of a stage to one scale, and the block is nineteen times smaller than the rest. At one scale the
typed state gains a step in every output of the round: 0.90 to 0.73 with no record of IRF1, 0.68 to 0.49 with its
record. The same typing of what enters the first layer reads as much with a record, so this is a way around the
encoder, not a repair of it; and nothing closes. The round's reader needs a rule for the scale of a stage before
any state made of blocks is judged by it. (That run also said the typed state reads a step more of an event without
a record. That held for a readout of control cells alone: a reader fitted with the cells of the pathway's events
leaves 0.64 at the mean state and at the typed state.)

The two open steps were then designed again and the audit was run again
([`runs/typed_repair_20261125`](runs/typed_repair_20261125/)); nothing was trained. An exploration on the THP-1
cells found that **the repair of the round did not apply its rows**: its gain weighs a row against a control fit that
is sure of cells it never saw, and for IRF1 it closed 27% of the row on the cells it was fitted on, so gate 3 never
passed. With the uncertainty of the control fit widened until the rows are not surprising, IRF1 with its own record
leaves 0.17 [−0.01, 0.29] at the plain mean state where it left 0.68; gate 3 passes and the grade is relative. A row
derived from the RNA of IRF1's cells, by a rule that projects the event's shift of the input view on the shift of
the recorded events of its type, is 92% of the measured one and leaves 0.22. Those THP-1 readings are not counted:
the rule was found there. On the melanoma strata, declared first, the picture is mixed. Every MHC-I event read on
its own gains a step with its record (5 of 5) and with a derived row (5 of 5); the rule gives 75% to 107% of the
measured change for 15 of 15 events with one left out; a permuted state does not close. But the pooled MHC-I target
does not move (0.36: its rows were already applied, and the state is the limit), the guard on control cells is
tripped there, and PD-L1 under IFN-γ loses its grade because gate 3 reads an output that overshoots on the cells it
was fitted on. The typed state, with its genes now chosen on the fitted half, removes about 0.05 where the state is
the limit and never a step; four programs found among control cells with no record remove nothing; and the last
layers trained on 24 November add nothing once the repair applies a record. Three things follow. The item that
"one record cannot close" was a property of the repair, not of the encoder, and the earlier localisation of a loss
at the mean (0.40 against 0.18) is a statement about ridge readers. An event without a protein record can be
reached with no training, by a row derived from its RNA, as far as that rule holds; it has been tested on members
of a type and on one event that is unlike its type. And four rules of the loop are now in question (the scale of a
stage, the weight of a row, gate 3 for rows that are closed, the reader classes of the tool); they are the owner's
decision, as is a batch of events that no run has read, on which all of this would be tested.
