# The design of the outer loop, checked without reading a cell

Status: a check of a **proposed** design, [`OUTER_LOOP.md`](../../OUTER_LOOP.md). No cell was read and nothing was
fitted. Two things were done:

1. **A dry run.** The rules of the design were applied to the finished round of 20 November
   ([`teg_loop_20261120`](../teg_loop_20261120/)): does every target get one status, does every thing that was left
   get one reading and one lever, and what would the next step need?
2. **A back-test.** Twelve earlier results of this project were put to the routing table: would it have named the
   action that was later found to be the right one?

The check changed the design in eleven places, and it corrects one statement of the round of 20 November. Both are
listed below.

## 1. The dry run

`scripts/dryrun.py` reads the three result files of the round and writes [`results/open_items.json`](results/open_items.json)
and [`results/open_items.md`](results/open_items.md).

**Targets**

| readout | group | part | status | own output after the repair | added beyond Z |
|---|---|---|---|---|---|
| PD-L1 | melanoma, IFN-γ | development | the head's part closed (relative) | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] no added loss |
| PD-L1 | melanoma, co-culture | development | the head's part closed (relative) | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] no added loss |
| PD-L1 | THP-1 | development | the head's part closed (relative) | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] no added loss |
| MHC-I | melanoma, IFN-γ | replication | not a target | | |
| MHC-I | melanoma, co-culture | replication | not graded | 0.39 [0.34, 0.44] reject | 0.07 [0.05, 0.09] no added loss |
| HLA-E | melanoma, IFN-γ | replication | not a target | | |
| HLA-E | melanoma, co-culture | replication | not a target | | |

**Open items: 20.** Observations that bear on one question are one item.

| id | readout | group | reading | lever | question | value | units now → needed | the next step needs |
|---|---|---|---|---|---|---|---|---|
| O01 | PD-L1 | melanoma, IFN-γ | unresolved | E | is the repaired output inside the margin | 0.04 [−0.02, 0.11] | 5 → 6 | at the margin: the pool has no further unit of the kind |
| O02 | PD-L1 | melanoma, both conditions | pooled, and differs | S | does the whole repair carry across the two conditions | adds 0.24 and 0.23 beyond Z | | the schema already keeps them apart; new units to confirm |
| O03 | PD-L1 | melanoma, both conditions | pooled, may be the same | S | can the rows of the two conditions be merged | adds −0.01 and 0.05 | | units of the pool, twice |
| O04 | PD-L1 | melanoma, co-culture | unresolved | E | is anything lost before the head | Z: 0.12 [0.01, 0.19] | 5 → 213 | at the margin |
| O05 | PD-L1 | melanoma, co-culture | invalid | V | the RNA view reads worse than Z | −0.09 [−0.13, −0.04] | | no new data: a richer reader of RNA |
| O06 | PD-L1 | THP-1 | unresolved | E | is anything lost before the head | Z: 0.13 [0.02, 0.31] | 5 → 67 | at the margin; the split of O07 comes first |
| O07 | PD-L1 | THP-1 | pooled, and differs | S | IRF1 reads differently from its pool | 0.54 [0.43, 0.65]; pool 0.16 | | units of the pool: events of IRF1's kind |
| O08 | PD-L1 | THP-1 and melanoma | pooled, and differs | S | does the whole repair carry across the two cell types | adds 0.62 | | the schema already keeps them apart; new units to confirm |
| O09 | PD-L1 | THP-1 and melanoma | unresolved | E | can the rows of the two cell types be merged | adds 0.14 [0.09, 0.18] | 5 → 9 | units of the pool, after a census |
| O10 | PD-L1 | melanoma, IFN-γ | self event | C | the knockout of CD274 itself | 0.91 [0.87, 0.95] at Z | | the owner's decision |
| O11 | all | melanoma, IFN-γ | mistyped | S | how a cell gets its type | 96% of knockout cells typed right | | no new data: the type from the record |
| O12 | MHC-I | melanoma, IFN-γ | small response | E | no event passes on its own (pooled 3.1) | | | units of the pool, or the owner's decision |
| O13 | MHC-I | melanoma, co-culture | rule | P | gate 3 is absolute | 0.28 on its own cells | | no new data: an amendment |
| O14 | MHC-I | melanoma, co-culture | rule | P | the hypothesis had no control arm | | | no new data: an amendment |
| O15 | MHC-I | melanoma, co-culture | unreadable | V | a third is read by no linear reader of RNA or of Z | 0.31 [0.26, 0.36] | | no new data first (a richer reader); then the owner's decision |
| O16 | HLA-E | melanoma, IFN-γ | small response | E | pooled 2.2, wanted 3 | | 5 → 10 | units of the pool, or the owner's decision |
| O17 | HLA-E | melanoma, co-culture | small response | E | no event passes on its own (pooled 3.1) | | | units of the pool, or the owner's decision |
| O18 | all | | selection | P | events were screened on the cells that were audited | | | no new data: an amendment |
| O19 | all | | rule | P | a group ready as a pool with no admitted event | | | no new data: an amendment |
| O20 | all | | no final data | E | no unit of the knockout claim is unread | 0 | | units of the pool: a census, final units first |

**Counts**

| by lever | | by what the next step needs | |
|---|---|---|---|
| the model (M) | 0 | no new data | 7 |
| events (E) | 8 | units of the pool | 9 |
| views (V) | 2 | at the margin | 3 |
| schema (S) | 5 | the owner's decision alone | 1 |
| claim (C) | 1 | | |
| protocol (P) | 4 | | |

**What the dry run shows**

- Every target got one status and every item one reading with a lever. Four of the fifteen readings did not occur
  (`upstream`, `no rows`, `confounded`, `guard conflict`): those branches are not exercised by this round.
- No item is sent to the model. That agrees with the round's own result: after the repair the head adds nothing, so
  another change to the head has nothing to remove.
- Seven items need no new data. They are the five amendments to the protocol, the type from the record, and a richer
  reader of RNA. They come first in the order of Section 9, with a bridge reading.
- **For PD-L1 an absolute closure is not within reach of the present pool at a tolerance of 0.10.** In the IFN-γ
  condition one more event of the kind would be needed if the estimate holds, and the library has none: on the
  readout's own scale its other 200 targets change PD-L1 by 0.35 at most (median 0.04), where the five pathway genes
  change it by 0.69 to 0.96. In co-culture and in THP-1 the
  embedding's own estimate sits at 0.12 and 0.13, next to the tolerance, and tens to hundreds of events would be needed
  to decide on which side. By the design these three end as "at the margin": the estimate and its interval are the
  result. A wider tolerance would decide them, but the numbers have been seen; a tolerance chosen now can only be
  tested on new units.
- What new units could decide: whether rows can be merged across conditions and across cell types (O03, O09), what
  kind of event IRF1 is (O07), and whether the small responses of MHC-I and HLA-E can be followed (O12, O16, O17). All
  of these need the system's states for the rest of the two screens, which do not exist yet.
- There are no final units for the knockout claim (O20). Under the design nothing certified can be said about it until
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

## 5. A correction to the round of 20 November

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

## 6. What was not checked

- Nothing was run. Whether the loop reaches terminal states in a useful number of iterations on this pool is not known.
- The dry run exercises eleven of fifteen readings, on one round, of one head.
- The same person wrote the design and the check. The back-test uses this project's own history, which the design was
  written with in mind.
- The counts of units needed assume that the estimate holds and that new units are like the old ones. Every pooled
  interval of the round rests on five events and is approximate by the protocol's own number (ten).

## Reproduce

    python scripts/dryrun.py --round ../teg_loop_20261120/results --out results
    python scripts/dryrun.py --routes        # the routing table; fails if OUTER_LOOP.md prints another one
