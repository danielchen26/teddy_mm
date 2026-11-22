# The design of the outer loop, checked three times before it was run

Status: the checks of [`OUTER_LOOP.md`](../../OUTER_LOOP.md), which was proposed on 20 November and adopted on
22 November. Three things were done:

1. **A dry run** (20 November). The rules of the design were applied to the finished round of 20 November
   ([`teg_loop_20261120`](../teg_loop_20261120/)): does every target get one status, does every thing that was left
   get one reading and one lever, and what would the next step need? No cell was read.
2. **A back-test** (20 November). Twelve earlier results of this project were put to the routing table: would it have
   named the action that was later found to be the right one? No cell was read.
3. **A check of the design's own claims** (22 November, Section 5). What the dry run had said about the pool was
   checked against the library's knockouts (cells that earlier runs had read); the rules were written as code
   ([`anm_diag/outer.py`](../../anm_diag/outer.py)) and every reading was produced from a planted case
   ([`tests/test_outer.py`](../../tests/test_outer.py)); the text was read once more.

The first two checks changed the design in eleven places. The third changed it in seven more, and it withdrew one
statement of the dry run. All are listed below, with a correction of the round of 20 November.

## 1. The dry run

`scripts/dryrun.py` turns the three result files of the round into the target format of `anm_diag.outer` and writes
[`results/open_items.json`](results/open_items.json) and [`results/open_items.md`](results/open_items.md). The tables
below are the dry run as it stands after the third check.

**Targets**

| readout | group | part | status | own output after the repair | added beyond Z |
|---|---|---|---|---|---|
| PD-L1 | melanoma, IFN-γ | development | the model's part closed (relative) | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] no added loss |
| PD-L1 | melanoma, co-culture | development | the model's part closed (relative) | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] no added loss |
| PD-L1 | THP-1 | development | the model's part closed (relative) | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] no added loss |
| MHC-I | melanoma, IFN-γ | replication | not a target | | |
| MHC-I | melanoma, co-culture | replication | not graded | 0.39 [0.34, 0.44] reject | 0.07 [0.05, 0.09] no added loss |
| HLA-E | melanoma, IFN-γ | replication | not a target | | |
| HLA-E | melanoma, co-culture | replication | not a target | | |

**Open items: 21.** Observations that bear on one question are one item.

| id | readout | group | reading | lever | question | value | units now → needed | the next step needs |
|---|---|---|---|---|---|---|---|---|
| O01 | PD-L1 | melanoma, IFN-γ | self event | C | the knockout of CD274 itself | 0.91 [0.87, 0.95] at Z | | the owner's decision |
| O02 | PD-L1 | melanoma, IFN-γ | unresolved | E | is the repaired output inside the margin | 0.04 [−0.02, 0.11] | 5 → 6 | no new data first: the bridge reading with the type from the record (O18) changes the cells this is read on. After that one more event would be needed, and the library has none that passes admission |
| O03 | PD-L1 | melanoma, both conditions | pooled, and differs | S | does the whole repair carry across the two conditions | adds 0.24 and 0.23 beyond Z | | the schema already keeps them apart; new units to confirm |
| O04 | PD-L1 | melanoma, both conditions | pooled, may be the same | S | can the rows of the two conditions be merged | adds −0.01 and 0.05 | | units of the pool, twice |
| O05 | PD-L1 | melanoma, co-culture | invalid | V | the RNA view reads worse than Z (dominated view) | −0.09 [−0.13, −0.04] | | no new data: a richer reader of RNA |
| O06 | PD-L1 | melanoma, co-culture | invalid | V | no reader of any stage reaches support for any event (positive control) | | | no new data first (a richer reader may reach support); then the owner's decision: a tolerance for this target |
| O07 | PD-L1 | melanoma, co-culture | unresolved | E | is anything lost before the head | Z: 0.12 [0.01, 0.19] | 5 → 213 | at the margin |
| O08 | PD-L1 | THP-1 | unresolved | E | is anything lost before the head | Z: 0.13 [0.02, 0.31] | 5 → 67 | at the margin; the split of O09 comes first |
| O09 | PD-L1 | THP-1 | pooled, and differs | S | IRF1 reads differently from its pool | 0.54 [0.43, 0.65]; pool 0.16 | | units of the pool: events of IRF1's kind |
| O10 | PD-L1 | THP-1 and melanoma | pooled, and differs | S | does the whole repair carry across the two cell types | adds 0.62 | | the schema already keeps them apart; new units to confirm |
| O11 | PD-L1 | THP-1 and melanoma | unresolved | E | can the rows of the two cell types be merged | adds 0.14 [0.09, 0.18] | 5 → 9 | units of the pool, after a census; the split of O09 comes first |
| O12 | MHC-I | melanoma, IFN-γ | small response | E | no event of the round passes on its own (pooled 3.1) | | | units at hand: B2M and TAPBP pass the screen, three more are near the threshold |
| O13 | MHC-I | melanoma, co-culture | invalid | V | no reader of any stage reaches support for any event (positive control) | every event 0.24 to 0.35, rejected, at Z and at the RNA view | | no new data first (a richer reader); then the owner's decision: a tolerance for this target |
| O14 | MHC-I | melanoma, co-culture | rule | P | gate 3 is absolute | 0.28 on its own cells | | no new data: an amendment |
| O15 | MHC-I | melanoma, co-culture | rule | P | the hypothesis had no control arm | | | no new data: an amendment |
| O16 | HLA-E | melanoma, IFN-γ | small response | E | pooled 2.2, wanted 3 | | 5 → 10 | the owner's decision: the library has no event that moves HLA-E (largest 1.9) |
| O17 | HLA-E | melanoma, co-culture | small response | E | no event passes on its own (pooled 3.1) | | | units of the pool, or the owner's decision |
| O18 | all | melanoma, IFN-γ | mistyped | S | how a cell gets its type | 96% of knockout cells typed right | | no new data: the type from the record |
| O19 | all | | selection | P | events were screened on the cells that were audited | | | no new data: an amendment |
| O20 | all | | rule | P | a group ready as a pool with no admitted event | | | no new data: an amendment |
| O21 | all | | no final data | E | no unit of the knockout claim is unread | 0 | | units of the pool: a census, final units first |

**Counts**

| by lever | | by what the next step needs | |
|---|---|---|---|
| the model (M) | 0 | no new data | 9 |
| events (E) | 8 | units at hand | 1 |
| views (V) | 3 | units of the pool | 7 |
| schema (S) | 5 | at the margin | 2 |
| claim (C) | 1 | the owner's decision alone | 2 |
| protocol (P) | 4 | | |

**What the dry run shows**

- Every target got one status and every item one reading with a lever. Five of the fifteen readings did not occur
  (`unreadable`, `upstream`, `no rows`, `confounded`, `guard conflict`): those branches are exercised only by planted
  cases.
- No item is sent to the model. That agrees with the round's own result: after the repair the head adds nothing, so
  another change to the head has nothing to remove.
- Nine items need no new data: the amendments to the protocol, the type from the record, richer readers of RNA and
  of Z, and the question whether PD-L1 in the IFN-γ condition is inside the margin, which waits for the bridge
  reading.
- **In co-culture nothing can be said yet about what RNA carries.** For neither readout does any reader of any stage
  reach support for any event (O06, O13). By Section 6 of the design a rejection on such a scale is not read as
  `unreadable`: first a richer reader, then a tolerance declared for the target.
- **For PD-L1 in co-culture and in THP-1 an absolute closure is not within reach at a tolerance of 0.10.** The
  embedding's own estimate sits at 0.12 and 0.13, next to the tolerance, and tens to hundreds of events would be
  needed to decide on which side. By the design these two end as "at the margin": the estimate and its interval are
  the result. A wider tolerance would decide them, but the numbers have been seen; a tolerance chosen now can only be
  tested on new units. **For PD-L1 in the IFN-γ condition it is open**: the estimate misses by 0.005, on cells of
  which 4% were given the wrong type. (The first version of this page said "not within reach" for all three. That was
  withdrawn in the third check.)
- One target can be reached with states at hand: MHC-I in the IFN-γ condition, with events the round did not use.
- What new units could decide: whether rows can be merged across conditions and across cell types (O04, O11), what
  kind of event IRF1 is (O09), and whether the small response of HLA-E in co-culture can be followed (O17). These need
  the system's states for the rest of the two screens, which do not exist yet.
- There are no final units for the knockout claim (O21). Under the design nothing certified can be said about it until
  some are set aside, before anything else is read.

## 2. The back-test

"Run *n*" is the number in the project's ledger. Where a run has a folder here it is linked; runs 5 and 9 belong to
the earlier workspace and are cited as in [`teg_repair_test_20261119`](../teg_repair_test_20261119/).

| | what happened | what the routing says | |
|---|---|---|---|
| 1 | HLA-E and MHC-I were made targets of a repair with a response too small to follow ([`loop_20261112`](../loop_20261112/), [`repair_20261112`](../repair_20261112/)) | `small response` → E, or not asked; never a target | agrees with what the protocol added (gate 0) |
| 2 | Round 3 changed a cell type and three outputs together | one item, one lever per iteration | agrees (protocol rule 2) |
| 3 | A repair removed the loss and pushed the bone-marrow score beyond the guard; the next round held the guard with a term in the loss (runs 16, 17) | first draft: no reading for it | **added**: `guard conflict` |
| 4 | Heads trained with IFN-γ cells kept little of their gain on co-culture cells (runs 4, 6); the repair computed in one condition failed in the other ([`teg_repair_test_20261119`](../teg_repair_test_20261119/)) | `pooled, and differs` → S: the condition stays a field | agrees |
| 5 | The same test asked co-culture for an absolute pass where the embedding itself was not supported | the grade that can be reached is fixed at admission; what the reference stage leaves is its own item | agrees (as the round of 20 November did) |
| 6 | The knockout of CD274 was probed as if a repair could read it; run 9 had found that a knockout of a protein's own gene cannot be read from RNA | `self event` → C: a limit for the type | agrees |
| 7 | Run 9 screened pairs on half the cells and tested them on the other half; the round of 20 November screened events (A1) on the cells it then audited | first draft: nothing | **added**: the `selection` reading and the selection split. No reported number of the round depends on it: where the response was near the threshold no event was admitted |
| 8 | Run 9 used isotype controls: none of 824 control pairs responded | first draft: nothing | **added**: the readout control of Section 6 |
| 9 | Every donor event was read `unreadable` at a tolerance of 0.10, because that tolerance is not defined for the variation among a donor's cells ([`scope_20261116`](../scope_20261116/)) | first draft: `unreadable` → V, "change the measurement". Wrong | **added**: the positive control of Section 6. Without a reader that reaches support somewhere on the scale, the reading is `invalid` and the tolerance is declared first |
| 10 | The final test of the donor round asked for a fall of 0.034, fixed on the development study; on the final study the ceiling was 0.006 ([`round1_coverage_20261116`](../round1_coverage_20261116/)) | first draft: a rule with a fixed bound | **added**: a rule is stated against a reference read on the same units, and "not testable" where there is nothing to remove |
| 11 | In that final study laboratory, assay, antibodies and disease changed together | first draft: named under the schema only | **added**: the `confounded` reading |
| 12 | The round of 20 November read "not in RNA" for PD-L1 in co-culture from an RNA view that reads worse than Z; run 5, sealed, had the opposite order on the same cells | first draft: `unreadable` → V | **added**: the dominated view and "one reader specification" in Section 6. See the correction below |

Seven of the twelve cases were handled wrongly or not at all by the first draft of the table. Each of them is a case
in which this project had already made the mistake once.

## 3. Against the paper

Each principle of Section 3 of the design is tied to a passage of the paper (arXiv 2610.14549), and Section 14 lists
what is taken from it and what is not. Three points where the first draft went beyond the paper and was pulled back:

- The draft spoke of "the minimum schema". The paper's coarsest sufficient representation is an abstract quotient that
  needs every equivalence, and its tests "constrain adequacy without uniquely identifying a missing variable". The
  design now claims only the coarsest schema that has passed the tests run so far.
- The draft sent an item to a lever as a finding. In the paper a resolved violation "does not by itself locate a faulty
  component". A route is now a hypothesis with a rival and a control.
- The draft adopted a split on a rejection of the pooled description. The paper's criterion for a refinement asks for
  two things on new histories: the finer description agrees and the coarser one is violated. Both are now required.

## 4. What the check changed in the design

1. A step before any routing: can the reading be used (dominated view, positive control, one reader specification,
   readout control, reader classes).
2. A terminal state "at the margin", computed from the units an item would need.
3. Rules stated against a reference on the same units; "not testable" as an outcome.
4. A selection split.
5. Two readings: `confounded`, `guard conflict`.
6. Observations grouped into questions: the first count had 24 items for 20 questions.
7. Final units stratified by every field of the record, so that a later split still has final units on both sides.
8. A control for every kind of outer hypothesis, in place of one control for all.

A second reading of the written design found three more things, now changed: an interval can be wide because the units
of a pool differ, so a split comes before more units (the first text said no change to the schema can move an
unresolved item); a split may be operated provisionally and a merge may not, which settles which schema a package
seals while a revision is still being tested; and the final units were described in two ways.

## 5. The third check (22 November)

It was asked for before the design was allowed to run. Three things were done.

**The dry run's statements about the pool were checked against data.** `scripts/pool_screen.py` screens every
knockout of the library in the IFN-γ condition for each readout
([`results/pool_screen.json`](results/pool_screen.json)); these are cells earlier runs had read, all cells of an event,
with an analytic interval.

| readout | events that pass A1 | of them not the readout's own gene | of them not used by the round |
|---|---|---|---|
| PD-L1 | 6 | the five pathway genes | none (the next is at 2.1) |
| MHC-I | 7 | B2M (13.9), TAPBP (9.7), SSR2 (3.5), STAT3 (3.2), CD44 (3.1) | all five |
| HLA-E | 0 | | none (the largest is 1.9, HLA-E's own gene) |

**The rules were written as code and every reading was planted.** `anm_diag/outer.py` derives the open items from a
plain description of the targets; `tests/test_outer.py` builds one target for each of the fifteen readings and checks
that it gives that reading and no other, and that each of the rules O1 to O6 catches the fault it is written for. The
dry run above is that code applied to the round.

**The text was read once more.**

What this changed:

1. *A statement of the dry run was not established and was partly wrong.* It said the pool had no further event "of
   the kind" for PD-L1, from numbers of another run and without the design's own screen. Screened, it holds for PD-L1.
   It does not hold as a statement about the pool: for MHC-I in the IFN-γ condition the library has two strong events,
   B2M and TAPBP. "Not a target" there was a consequence of the five events the round chose. A census now lists, for
   units that have been read, which events pass the screen for each readout.
2. *Forecasts were made as if the items were independent.* Whether PD-L1 in the IFN-γ condition is inside the margin
   (0.04 [−0.02, 0.11]) was read on cells of which 4% had been given the wrong type, and the type from the record is
   one of the changes that need no new data. The headline "an absolute closure is not within reach" was withdrawn for
   that condition. The order of an outer step now derives the open items again after the bridge reading.
3. *An amendment named a control that does not always exist.* "Rows of events that do not move the readout" exist
   with states for one group of three. The control that always exists is the readout fitted on the group's control
   cells without rows, which the round already had as an arm.
4. *The statuses were incomplete.* A target whose model still adds a loss, and one for which that is undecided, had no
   status; an item decided on cells that had been read had no state. Both are added.
5. *A loss before the model with no input view read* had no reading. It is `invalid`: the loss is not placed until the
   input view is read.
6. *A census of this pool already existed, and the design and its dry run did not use it.*
   [`runs/loop_reading_20261112`](../loop_reading_20261112/) (Section 6 there, `results/data_census.json`, written on
   13 November) lists, for both screens and every condition, the knockouts that change a measured protein clearly. It
   already names B2M and TAPBP for MHC-I in all three conditions of the melanoma screen, CMTM6 and BRD4 for PD-L1 in
   THP-1, and NPC1 for PD-L1 in co-culture. Point 1 above was found again by a new screen when it was on file. Two
   things follow. The census step of the design starts from that file. And the design's sentence "a census reads no
   response of a unit that has not been read" cannot be kept for this pool: the census of 13 November read the protein
   response of every target on all its cells. What has not been read for the units outside the audits is the system:
   no state was computed for their cells and no *u*. So an event near the admission threshold was chosen by a
   statistic of all its cells, and an audit on those cells would be biased towards a loss. Final units are final in
   the sense that the system was never read on them, and a final audit admits its events on a selection split drawn
   then (Section 8.1 of the design, changed accordingly).
7. *One check of Section 6 was not applied by the dry run.* The positive control ("no reader reaches support for any
   unit of the family on this scale") was in the design and in the code, but the dry run's adapter did not supply
   it, and a check whose input is missing passes. Supplied, it fails for both readouts in co-culture: of the
   readings of the five events and their pool, at Z and at the RNA view, none is supported (PD-L1: unresolved or
   rejected, the nearest is JAK1 at −0.06 [−0.19, 0.05]; MHC-I: every one rejected, 0.24 to 0.35). So the item that
   read `unreadable` for MHC-I in co-culture is `invalid` (O13), and PD-L1 in co-culture has a second invalid reading
   (O06). The code now draws no reading of where a loss is from a view that fails any check, and a planted case
   tests it. This was found on 22 November while the first outer step was being prepared, after the three checks.

What it did not change: the five levers, the checks of Section 6, the routing of the other readings, the tests of a
split and of a merge, the final units.

## 6. A correction to the round of 20 November

The README of [`teg_loop_20261120`](../teg_loop_20261120/) says that in co-culture what remains "is not in RNA either".
By the check of Section 6 that does not stand for PD-L1:

- In that round the RNA view (512 components, a linear reader) leaves 0.19 [0.11, 0.25] with the event's own cells.
  Z, which is computed from the same RNA, leaves 0.10 [0.01, 0.17]. The gap is resolved below zero
  (−0.09 [−0.13, −0.04]): the RNA view is the weaker reader, and its rejection says nothing about RNA.
- Run 5, sealed, on the same cells, had the opposite order: RNA and tokens 0.12, unresolved; Z 0.19.
- What stands for PD-L1 in co-culture: no reader of RNA or of Z that was tried comes inside the margin; the lowest
  estimate is about 0.10, with an interval from 0.01 to 0.17. Whether anything is missing from RNA is not decided.

For MHC-I in co-culture the two views agree (0.31 at both, each rejected; gap −0.00 [−0.04, 0.03]), so the view is not
dominated. Only a linear reader was used at either stage. The statement that stands is: a third of the change is read
by no linear reader of RNA or of Z. A richer reader comes before any conclusion about the input.

A note with this correction is added to that README. Its numbers do not change.

## 7. What was not checked

- Whether the loop reaches terminal states in a useful number of iterations on this pool is not known.
- The dry run exercises ten of fifteen readings, on one round, of one head. The planted cases show that the rules
  are consistent with each other, not that they are the right rules.
- The same person wrote the design and the check. The back-test uses this project's own history, which the design was
  written with in mind.
- The counts of units needed assume that the estimate holds and that new units are like the old ones. Every pooled
  interval of the round rests on five events and is approximate by the protocol's own number (ten).

## Reproduce

    python scripts/pool_screen.py --work WORK --out results        # WORK: the work folder of runs/teg_loop_20261120
    python scripts/dryrun.py --round ../teg_loop_20261120/results --out results
    python scripts/dryrun.py --routes        # the routing table; fails if OUTER_LOOP.md prints another one
