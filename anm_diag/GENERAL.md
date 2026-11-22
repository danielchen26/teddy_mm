# The loop for any system

Written on 22 November 2026, after the first outer step. Run and checks:
[`runs/framework_generalization_20261122`](runs/framework_generalization_20261122/).

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
| stages inside the model | optional: arrays per cell at points between the input and the stage the model reads, in order | none was read | the output of each modality's encoder; the fused state | `mid` |
| other measured inputs | arrays the model does not receive, never the readout itself | none was read | a modality that was measured and left out | `view_b` |
| a control view | a view that must read nothing | the RNA view with its rows permuted | the same | `perm` |
| events with control cells | units that change the readout, each with cells that did not get the event | knockouts against non-targeting guides | perturbations of any kind | planted shifts |
| a repair | an object fitted on control cells and on the rows of the admitted events, with one output per arm (§3) | a linear readout of Z, updated once by the rows | a change to whatever part the round may change | the same linear readout |
| a census | which units of the claim exist, which were read, which are set aside as final | 615 units, 133 final | the same | not used |

Two things are the system owner's and stay outside the code: the claim (which readouts, which types, what tolerance)
and which part of the model a round may change.

## 3. The calls

    sel = cycle.select(ids, unit, y, seed)                      # A1 and A2 on a third of the cells; those cells are read no further
    grid = cycle.grid(scope.run_spec(spec, folder))             # the audit, or cycle.load(folder) for a grid on disk
    T = cycle.one_type(readout, group, type, x, y, o, treated, unit, sites, *grid,
                       model_stage, inputs=[...], between=[...], others=[...], control=..., alone=[...])
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

- **Events without control cells of their own.** A donor, a batch, a study or a disease is read against the other
  units, and both sides are resampled by unit. The audit reads such events. The round does not: its repair is fitted on
  control cells, and its halves are halves of cells. Before donors and diseases enter the loop this has to be designed
  and planted like the rest.
- **Every view on the same cells.** Gaps between stages are read on paired resamples, so every stage has to exist for
  every cell that is read. A modality that is missing for part of the cells is not handled.
- **Events are yes or no, and a type has one partner.** There is no dose. What carries across types is tested against
  one other type per run.
- **One quantity.** The shift of the mean of a readout. The audit also reads the variation among the cells of a unit;
  the round does not use it.
- **One repair.** A linear readout of the stage the model reads. Nothing that changes an encoder has been put through
  the loop, on any system.
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

## 7. Plan

| | step | needs | whose |
|---|---|---|---|
| done | the general cycle, the planted system, the comparison on the first system | CPU | |
| 1 | events without control cells (donor, disease) in the round: design, planted donors, forecasts written first | CPU | |
| 2 | a second system by swapping the encoder on the same cells: the loop is run unchanged, only the stage arrays differ | inference | the owner names the encoder |
| 3 | a repair in front of the head at small scale: the head reads Z together with a view of the input. It tests whether a loss placed "upstream" (IRF1) closes without retraining | CPU | |
| 4 | the states of batch 1 (6,759 cells) and the round on them | inference on the owner's machine | the owner |
| 5 | **retraining TEDDY with the TEG as a constraint in pre-training.** The rows of the TEG enter the training loss; the schema is sealed for an iteration; the retrained encoder is one more system for the same loop, judged by the same items | several GPUs, the pre-training data | the owner; after 1 to 3 |
| 6 | the final audit on the 133 units set aside, once and with the hypotheses named | inference | the owner |

Step 5 is placed after 1 to 3 on purpose: they cost no GPU, and 3 says whether the upstream loss is worth a
retraining at all.
