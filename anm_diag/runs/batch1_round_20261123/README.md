# The round on batch 1

23 November 2026. Step 4 of the plan in [`GENERAL.md`](../../GENERAL.md): the system's states were computed for the
6,759 cells of batch 1, and one round was run on them. It is the first round on units the system had never been read
on, and the first one made entirely by the general code (`anm_diag.cycle`, `anm_diag.outer`). The design, the roles
of the units and 29 forecasts were committed before the states existed ([`DESIGN.md`](DESIGN.md)); the scripts before
any protein value of these cells was read.

Everything here is development. No final unit was touched.

## What came out

**No event of the batch gives the model anything to repair.** Three of the eight units pass admission, and all three
act on the measured protein itself. Neither the model nor any reader of its embedding or of its input reads them.

| stratum, readout | event (change, control SD) | the output as trained | a reader of Z with the event's own cells | the output loses beyond that reader | status |
|---|---|---|---|---|---|
| THP-1, PD-L1 | CMTM6 (−0.59) | 1.00 [0.97, 1.04] | 1.05 [0.93, 1.22] | −0.04 [−0.21, 0.08]: no added loss | open before the model |
| melanoma, co-culture, MHC-I | TAPBP (−1.97) | 1.02 [0.99, 1.04] | 1.01 [0.95, 1.08] | 0.01 [−0.06, 0.07]: no added loss | open before the model |
| melanoma, control condition, MHC-I | B2M (−2.04) | 1.00 [0.99, 1.00] | 0.89 [0.86, 0.93] | 0.11 [0.07, 0.14]: unresolved | undecided in the model |

*u* is the share of the measured change that is left unexplained: 1 is all of it.

- **A3 holds for none of them**: the output as trained does not lose a step more than the reference. So no repair is
  graded. The arms were still computed: every one stays where the reference is (CMTM6 1.02 to 1.06, TAPBP 1.00 to
  1.03, B2M 0.89 to 0.93), and the event's own record adds no step over a readout of control cells.
- **The input reads none of them either.** The best reader of RNA with the event's own cells leaves 0.99 (CMTM6),
  0.99 (TAPBP), 0.95 (B2M). No reading of any stage, scope or reader class is supported (0 of 20 for each target), so
  the positive control fails and the items are `invalid`, not `unreadable`. The permuted view reads nothing (1.00),
  as it must.
- **B2M in the control condition is the one place with a trace.** A reader of Z with its own cells reads about a
  tenth of the change (0.89). The head's own output reads none of it. The difference, 0.11, lies across the step
  margin, so nothing is decided.
- **Five of eight units are not targets.** BRD4 (THP-1, PD-L1) changes PD-L1 by 0.28 control SD on the selection
  third, 1.4 half-widths; the screen of 13 November had 0.67 on all its cells. CD44 (co-culture, MHC-I): 0.39, 1.7
  half-widths. The three calibration units (POU2F2, CDH19, KLF4) pass for no readout. No unit passes for PD-L1 or
  HLA-E in the two melanoma strata.

## What it means

- It **replicates a limit** on unread units, in two more conditions and for one more mechanism. Under IFN-γ the first
  outer step found that B2M and TAPBP (which build and load the MHC-I complex) and CD274 (which encodes PD-L1) are
  invisible to a system that sees RNA only. The same now holds for TAPBP in co-culture, B2M in the control
  condition, and CMTM6, which keeps PD-L1 on the membrane after it is made. Whether such knockouts belong to the
  claim is the owner's decision; it was already open.
- It is **not a confirmation of the repair**. The batch held one unit that could have tested it, BRD4, which acts on
  transcription. It was not admitted. The statement "the model's part is closed (relative)" for PD-L1 and MHC-I
  still rests on development cells that had been read before.
- A unit chosen because an earlier screen found it clear is **not thereby admitted**. Two of the five chosen units
  fell below the bar on a third of their cells. Choosing by a screen and admitting on a split are different steps,
  and the split did its work.
- The **calibration third** says: of three units drawn at random, none is a target. With three units this is a first
  count, not a rate.

## The forecasts

Intervals: 6 of 7 held the value. Outcomes: for 13 of 16 the class given the highest probability happened. Six
forecasts were conditional on BRD4 being admitted and are not scored. ([`results/outcome.json`](results/outcome.json))

| missed | forecast | observed |
|---|---|---|
| BRD4 is admitted | yes 0.5, no 0.5 | no |
| CMTM6: *u* at Z, own cells | [0.45, 1.0] | 1.05 |
| CMTM6: the earlier repair adds no loss beyond the reference | yes 0.45 | yes (nothing is read, so nothing can be added) |
| MHC-I, control condition: status | open before the model 0.55 | undecided in the model |

## Open items of the batch

Nine, derived by `outer.open_items` ([`results/round.json`](results/round.json), `open_items`): three `invalid`
readings sent to the views (no reader of any stage reads CMTM6, TAPBP or B2M), one `unresolved` (B2M in the control
condition: 0.11 at the step margin), and five `small response` (strata and readouts for which no unit of the batch
passes). The items of the first outer step that waited for new THP-1 events (O08, O11, O12) are not moved: CMTM6 is
lost before the model and BRD4 is not a target.

## What was not as written

1. **The first admission did not screen THP-1.** `build.py` wrote the arrays of a stratum under the name of a file
   of the earlier rounds that was linked into the work folder. A trial of the script on made-up protein values had
   therefore replaced the THP-1 file before the real run, which then found no control cell in the selection third
   and admitted nothing for THP-1. The earlier files were regenerated from their sources (an untouched one is
   reproduced byte for byte), the arrays of this run were renamed, and admission was run again. The file of the
   failed attempt is kept ([`results/selection_first_attempt.json`](results/selection_first_attempt.json)); the
   tables of the two melanoma strata are identical in both.
2. **The states were computed in a new environment** (Python 3.9, torch 2.8.0, MPS), because the one the earlier
   states were made with was no longer on the machine. The check passed: on 192 cells that have stored states the
   largest difference is 1.2e-5 ([`results/states_manifest.json`](results/states_manifest.json)).
3. The symbol map of the melanoma screen was read from a stored copy of the same file.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/roles.json`](results/roles.json), [`results/declared.json`](results/declared.json) | the declaration: what is asked, the roles, the forecasts |
| [`results/states_manifest.json`](results/states_manifest.json) | the states: cells, environment, the check against stored states |
| [`results/selection.json`](results/selection.json) | admission on the selection third, every unit and readout |
| [`results/round.json`](results/round.json) | the reading: every event, every arm, the placing, the targets, the open items |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/states.py`, `pack.py` | run on the owner's machine: RNA counts, tokens, states; the record and the measured protein |
| `scripts/build.py`, `round.py`, `score.py` | admission and arrays; the round; the score |

To run again (the work folders of `teg_loop_20261120` and `outer_step_20261122` are needed for the earlier events):

    python scripts/states.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cells batch1_cells.tsv --out OUT
    python scripts/pack.py --repo TEDDY_MM --cells batch1_cells.tsv --out OUT
    python scripts/build.py --states OUT/batch1_states.npz --pack OUT/batch1_pack.npz --counts-frangieh OUT/counts_frangieh2021.npz --counts-papalexi OUT/counts_papalexi2021.npz --work WORK
    for s in WORK/b1_*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/round.py --work WORK --observe OBSERVE
    python scripts/score.py

## What is next

- A test of the repair on unread units needs events that act through transcription and pass admission. The pool has
  258 unread units left that are not final. A screen that is itself made on a split would propose them without the
  shrinkage seen here.
- The owner's decision on knockouts that act on the measured protein is now backed by six events in four conditions.
- Steps 1 to 3 of the plan are unchanged and need nothing from the owner.
