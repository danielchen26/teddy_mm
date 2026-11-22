# The first outer step: declaration

Written on 22 November 2026, **before any reading of this step was made**. It applies
[`OUTER_LOOP.md`](../../OUTER_LOOP.md) to what the round of 20 November
([`teg_loop_20261120`](../teg_loop_20261120/)) left: 7 targets, 21 open items
([`outer_loop_design_20261120/results/open_items.md`](../outer_loop_design_20261120/results/open_items.md)).

The step changes the package, not the model. It ends with package 2 sealed, and one round of the inner loop is then
run under it, on the one target that states at hand can reach.

Everything here is exploratory. No unit of the claim is final yet; setting some aside is one of the actions.

## 1. What preparing the step found

Two things, both in the design's own dry run, both recorded in
[`outer_loop_design_20261120`](../outer_loop_design_20261120/) (Section 5, points 6 and 7):

- A census of the pool was already on file ([`loop_reading_20261112`](../loop_reading_20261112/), 13 November). It
  had read the protein response of every knockout of both screens on all cells. So "unread" can only mean: the system
  was never read on the unit, and every admission from now on uses a selection split.
- The positive control of Section 6 had not been applied. Applied, it fails for both readouts in co-culture: no
  reader of any stage reaches support for any event. The item that read `unreadable` is `invalid`.

## 2. Package 2 against package 1

| | package 1 (the round of 20 November) | package 2 | items |
|---|---|---|---|
| **P** protocol | version 2 | version 3: five amendments ([`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md), 3a to 3e) | O14 O15 O19 O20 |
| **S** schema | a type is cell type × condition; a reader of the state gives each cell its type; the five THP-1 events are one pool | the type comes from the record (3e); the THP-1 pool is read as two types, provisionally: four genes that signal to STAT1, and IRF1; "relation to the readout" is named as a candidate field, not operated | O18 O02 O09 O08 O11 |
| **V** views | Z: ridge and boosting (ridge only for MHC-I in co-culture). Input view: 512 components of RNA, ridge | one reader specification for every target: Z by ridge and boosting; the input view by ridge and boosting on 512 components and by ridge on 4,000 genes; a permuted view as control; the positive control evaluated | O05 O06 O13 |
| **E** events | five pathway genes per group, admitted on all their cells | the same, kept; for MHC-I in melanoma under IFN-γ the events are selected from the library on a selection third; a census from the record; a third of the unread units of every stratum set aside as final | O12 O21 |
| **C** claim | readouts PD-L1, MHC-I, HLA-E; margin 0.10; guard 0.02 | unchanged | |
| **M** model | the head with a repaired output per group | unchanged in this step | |

Not acted on: O01 (the owner's decision), O03 O04 O10 (need units the system has not been read on; none has states),
O07 (read again after V1), O16 O17 (the owner's decision, or units without states).

## 3. The order

As Section 9 of the design. Scripts are in `scripts/`; each writes one file in `results/`.

1. **Check the readings.** Done in the dry run, with the positive control added.
2. **Amend the protocol and the schema** (this file, `LOOP_PROTOCOL.md` version 3).
3. **Bridge reading** (`bridge.py`). The cells and the model of the round, read again under the new rules. No new
   data, no change to the model. Four targets: PD-L1 in the three groups and MHC-I in co-culture; and the THP-1 pool
   under the provisional split, beside the unsplit pool.
4. **Richer readers** (`readers.py`), for the two co-culture targets.
5. **Derive the open items again** (`items.py`).
6. **Select units** (`selection.py`, `census.py`): the selection third for MHC-I in the IFN-γ library; the census from
   the record; the final units.
7. **Seal package 2** (`results/package2.json`).
8. **One round under package 2** (`round2.py`): MHC-I in melanoma under IFN-γ with the admitted events.

## 4. How things are done

- **Type from the record.** Every cell of a group is read by the readout of its group. A bone-marrow cell keeps the
  head's own output, so the guard on the head's task cannot move. The reader of types is reported as a component
  (from the round: 95.6% to 100% of held-out cells given their own type).
- **Halves.** As in the round: every cell is read by a readout fitted without it; the same halves, the same seeds, so
  the bridge reading differs from the round by the rules only.
- **Gate 3, amended.** On the cells a repair was fitted on: the interval of |*u* of the output| − |*u* of the
  reference stage| below the step margin. The output as trained is put through the same gate as the control.
- **The test, amended.** `loop.test` with the rows arm as the own arm and two controls, each with the bound of the
  round (the grid's ceiling minus the step margin). The difference between the rows arm and the controls-only arm is
  also read on paired resamples against the step margin.
- **The split.** Type A: the THP-1 control cells and the cells of IFNGR1, IFNGR2, JAK2, STAT1; a new grid, because
  its reference (a reader of the other genes of the type) may not contain IRF1. Type B: IRF1 with the same control
  cells; one event, so its reference is a reader of Z with its own cells (the grid of the round). The unsplit pool is
  read beside both.
- **The input view.** RNA as log1p of counts per 10,000, genes detected in at least 1% of the group's control cells
  (as the round). Two views of it: 512 components (as the round), and the 4,000 genes with the largest variance among
  the group's control cells, unscaled. No label and no readout enters either. Readers: ridge (13 weights from 0.1 to
  100,000, chosen inside the fit) and gradient boosting (150 trees, 15 leaves, rate 0.05), as the round; boosting is
  not run on the 4,000 genes. The permuted view is the 512 components with their rows permuted across the cells of
  the group (seed 20261122), read by ridge.
- **A gap between two stages** is the difference of |*u*| on the tool's paired resamples, own-cells scope, each stage
  under its best reader.
- **The positive control** of a target holds if any reader class, at any stage, supports any event of the family or
  their pool.
- **The selection third.** A cell is in it if the SHA-256 of `20261122|<cell name>` is divisible by 3. For a
  candidate event: its cells in the third against the control cells in the third; admitted if the change is resolved
  and at least 3 times the half-width of its interval (analytic, cells as units), and the gene is not one of the
  readout's own (HLA-A, HLA-B, HLA-C for MHC-I). Candidates: every knockout of the IFN-γ library. The cells of the
  third enter no reading of *u*, for this readout, from now on.
- **The new target.** Cells, all outside the selection third: the control cells, the admitted events' cells, and the
  cells of the five pathway genes as other units of the type (they are no events here: none passed A1 for MHC-I).
  Two events give no statement about a family; every reading is per event, resampled by cell.
- **Final units.** From the record only (`record_census.py`: cell names, guides, design fields). A unit is a target
  gene in a condition of a study with at least 100 cells that carry guides for that gene only. A unit is unread if no
  state of the system exists for its cells. Per stratum (study × condition), the unread units in alphabetical order
  are permuted with seed 20261122 and the first third (rounded) is final. An unread control cell is final if the
  SHA-256 of `final|<cell name>` is divisible by 3.

## 5. The hypotheses

Eleven, one per action, each with a forecast, a rule, a rival and a control:
[`results/declared.md`](results/declared.md) (`scripts/declare.py`, which also checks that every one of the 21 items
is either acted on or listed with the reason why not). The forecasts that matter most:

| | forecast |
|---|---|
| S1 | With the type from the record, PD-L1 in melanoma under IFN-γ reads 0.02 (−0.01 to 0.05) and is supported with probability 0.6. Then that target is closed (absolute), on development cells |
| P1 | Gate 3 as amended passes for all four targets; MHC-I in co-culture becomes graded |
| P2 | With the controls-only arm as a control, the test of "rows" is unresolved for PD-L1 under IFN-γ (0.85): a readout fitted on that group's control cells already removes most of the loss. The rows add a resolved step in co-culture, in THP-1 and for MHC-I |
| S2 | Without IRF1 the THP-1 pool reads about 0.10 and stays unresolved (0.7). IRF1 alone stays rejected (0.9); whether its loss is upstream of Z is even odds |
| V1 | A reader of 4,000 genes brings the input view of PD-L1 in co-culture to 0.08–0.20; the view is no longer dominated with probability 0.6. A positive control is not reached (0.85) |
| V2 | For MHC-I in co-culture every reader at every stage still rejects (0.8); no positive control (0.9). Then the next action is the owner's |
| E1 | B2M and TAPBP are admitted on a third (0.95); the three near the threshold are not (0.85 each) |
| E2 | The head as trained loses nearly all of both (0.85 to 1.0). A reader of Z with the event's own cells still rejects B2M (0.85): these knockouts may act like a knockout of the readout's own gene. Whether the group is a target for a repair by rows is slightly less than even (0.45) |

## 6. What this step does not do

- No final audit, and nothing certified.
- No test of the split and no merge: both need units the system has not been read on, and none has states. Which
  units to compute states for is an output of this step.
- It does not decide what is the owner's: a tolerance per target, a limit for knockouts of a readout's own gene, a
  system that is told the perturbation, a new measurement.
- It does not change the model, except inside the round on the new target.

## 7. How it is graded

Every forecast is scored when its reading exists (`items.py` writes `results/outcome.json`): inside or outside its
interval, and for outcomes the probability that was given to what happened. The ledger rules O1 to O6 are checked on
the declared diff and on the items as they stand after the step. An item read again on the units of its proposal can
become "decided on development units" and no more.
