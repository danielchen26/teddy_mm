# The loop without its first system: a general cycle, a planted system, and a re-audit (22 November 2026)

Status: code and checks. No new cell of the first system was read, the model was not changed, and no forecast was at
stake. What a system has to expose, and the plan, are in [`GENERAL.md`](../../GENERAL.md).

The owner's instruction: retraining TEDDY with the TEG at the level of pre-training goes into the plan as a later
step; first the framework is made general, and audited again.

## 1. What was done

1. **A general cycle** ([`anm_diag/cycle.py`](../../anm_diag/cycle.py)). One round of the loop, with no model,
   modality, readout or stage named in it: selection on a share of the cells, a repair as a plug-in, the reading of a
   round (A3, every arm against the reference stage, gate 3, the two controls), the placing of what is lost before the
   model, and the target that the outer step reads. Until now these lived in the scripts of two run folders.
2. **A planted system** ([`anm_diag/toy.py`](../../anm_diag/toy.py)) in which the place of every loss is known, and a
   second version of it with nothing wrong.
3. **A score over seeds and sizes** (`scripts/planted.py`), by a rule written before the runs.
4. **A comparison on the first system** (`scripts/regress.py`): the general code on the cells of the first outer step,
   against every number that step recorded.
5. **Six rules changed** where the planted system showed the loop saying something the plant did not support
   ([`OUTER_LOOP.md`](../../OUTER_LOOP.md), Section 17), and two parts of the round that the general code lacked.

## 2. The planted system

| | |
|---|---|
| modality a | carries factors 1, 2, 4 and 5; it is the model's input |
| modality b | carries factor 3; measured, and not given to the model |
| first step of the model | keeps factors 1, 2 and 4 and drops factor 5 (stage `mid`) |
| second step | keeps factors 1 and 4 and drops factor 2 (stage `z`, which the head reads) |
| head as trained | follows an unrelated direction of `z` |
| readout | type one: factors 1 + 2 + 3 + 5; type two: factors 4 + 2 + 3 + 5; measured outside the system |

| type | event | acts on | what the loop should say |
|---|---|---|---|
| one | p1, p2, p3 | factor 1 | the head loses it; a repair of the head closes it |
| one | early | factor 5 | upstream, lost in the first step |
| one | enc | factor 2 | upstream, lost in the second step |
| one | modb | factor 3 | unreadable from the model's input; another measured input reads it |
| one | post | the readout itself | unreadable; no input reads it |
| one | weak | factor 1, by 0.03 | refused at selection |
| two | q1, q2, q3 | factor 4 | the head loses it along another direction: the repair of type one does not carry over |
| two | odd | factor 2, and pooled with q1 to q3 | one event inside the pool differs; upstream, lost in the second step |

The sound version has a head that reads the right factor and only the events it can read. The loop should leave no
open item and read no repair.

## 3. The score

Every planted fact is scored for every run: **right** (the loop says what was planted), **undecided** (it says less:
no item, an interval that does not decide, a reading it refuses to use, a weaker status) or **wrong** (it says
something else that is definite). The rule is in the docstring of `scripts/planted.py` and was written before the
runs. Twelve seeds at three sizes, each for the planted and for the sound system: 72 runs, 504 scored facts.

| size | control cells per type | cells per event | right | undecided | wrong |
|---|---|---|---|---|---|
| large | 4,000 | 800 | 167 | 1 | 0 |
| medium | 1,500 | 300 | 158 | 10 | 0 |
| small | 500 | 100 | 140 | 27 | 1 |

A third of the cells serves selection and enters no reading, so the readings rest on two thirds of these numbers.

| planted fact | large: right / undecided / wrong | medium | small |
|---|---|---|---|
| head of type one: closed by its repair | 12 / 0 / 0 | 11 / 1 / 0 | 2 / 10 / 0 |
| early: upstream | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| early: lost in the first step | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| enc: upstream | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| enc: lost in the second step | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| modb: unreadable, and the other input reads it | 11 / 1 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| post: unreadable, and no input reads it | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| weak: refused at selection | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| the repair of one type does not carry to the other | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| odd: differs inside its pool | 12 / 0 / 0 | 9 / 3 / 0 | 7 / 5 / 0 |
| odd: upstream | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| odd: lost in the second step | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| no other definite item | 12 / 0 / 0 | 12 / 0 / 0 | 11 / 0 / 1 |
| a system with nothing wrong: nothing is open | 12 / 0 / 0 | 6 / 6 / 0 | 0 / 12 / 0 |

What the table says:

- **Where a loss is, is found at every size.** Upstream or unreadable, in which step of the model, and whether another
  input reads it (eight facts, 288 scores): 287 right, one undecided, none wrong.
- **What shrinks with the cells is the absolute statement.** "Closed" and "nothing is open" need an interval inside a
  margin of 0.10. With 67 cells per event the interval is wider than that, and the loop says "the model's part
  closed" or "undecided" instead. It does not call a sound system open.
- **The one undecided at the large size** (seed 2, modb): both the model's input and `z` read nothing of that event,
  and the difference between the two readings of nothing came out resolved below zero. The check for a dominated view
  then refused the input view. This is the check firing on noise. A change is proposed in Section 6 and not applied.
- **The one wrong** (small, seed 3): beside the right item for `early` (upstream, first step) the loop raised a second
  one, `unreadable`, because the reader of the model's input left 0.16 [0.100, 0.219] of that event with 67 cells.
  The input carries the factor in full; the reader did not learn it from so few cells. A reader that fails cannot be
  told from an input that lacks the response by its reading alone. Section 6 has the proposal.

Two things were set in the planted system so that what is planted holds at every seed, after first runs showed
otherwise: the loadings of a factor have a fixed size (a weak random loading left a stage carrying its factor with a
floor near the margin), and the factors the head does not read have a small spread among control cells (their
variance is noise for the slope that puts an output on the readout's scale). Both are in `toy.py` with a comment.

## 4. The first system, computed again

`scripts/regress.py` runs `anm_diag.cycle` on the cells, halves, seeds, grids and admitted events of
[`outer_step_20261122`](../outer_step_20261122/) and compares with that step's result files. It first checks that
the comparison sees a planted difference.

| what | values compared | largest difference | unequal |
|---|---|---|---|
| bridge: PD-L1, melanoma under IFN-γ | 304 | 0 | 0 |
| bridge: PD-L1, melanoma in co-culture | 304 | 0 | 0 |
| bridge: PD-L1, THP-1 | 304 | 0 | 0 |
| bridge: PD-L1, THP-1 without IRF1 | 263 | 0 | 0 |
| bridge: PD-L1, THP-1, IRF1 alone | 86 | 0 | 0 |
| bridge: MHC-I, melanoma in co-culture | 198 | 0 | 0 |
| round under package 2: B2M | 110 | 0 | 0 |
| round under package 2: TAPBP | 110 | 0 | 0 |
| readers: MHC-I under IFN-γ | 213 | 0 | 0 |
| readers: PD-L1 in co-culture | 420 | 0 | 0 |
| readers: MHC-I in co-culture | 420 | 0 | 0 |
| **all** | **2,732** | **0** | **0** |

The seven targets get the same statuses, and the outer module derives the same 17 open items from the targets the
general code builds. The rule for choosing the reader of a view was changed to the audit tool's own (Section 5); on
all three placings the two rules pick the same reader with the same verdict.

## 5. What the re-audit changed

Six rules said something the plant did not support. Each is changed in `outer.py`, has a test, and is listed in
[`OUTER_LOOP.md`](../../OUTER_LOOP.md), Section 17. The items and statuses of the first step were derived again after
every change and are the same.

| | found by | what was wrong | now |
|---|---|---|---|
| 17a | an event lost partly at the input and mostly after it | only `unreadable` was reported | both places are reported |
| 17b | a pool that the reference stage reads | its input view was checked and could raise an item | where the reference stage is supported the input view is not consulted |
| 17c | events read on their own; and `odd` inside its pool | the first were called "pooled, and differs"; the second went unflagged where no repair was read | only an event inside a pool can differ, and it is flagged when it is rejected at the reference stage while others of the pool are supported |
| 17d | the sound system | a reader of a stage that failed for one event raised an item although the output reads that event | an output inside the margin has lost nothing on the way |
| 17e | `modb` | "another measured input reads it" asked for support at that input | also met when that input reads a step more than the model's own |
| 17f | the sound system at the medium size | "open before the model" although every reading was only undecided | said only where a stage before the model is rejected; otherwise "undecided before the model", an eighth status |

Missing from the general code, and added:

- **A3.** Whether the output as trained loses a step more than the reference stage. The first system's scripts applied
  it before a round; the general round had no such step, so a repair could have been graded where the model had
  nothing to remove. It is part of `read` now, and a target without it is described by the output as trained.
- **One rule for choosing a reader.** `place` took the reader with the smallest |u|; the audit tool takes a view as
  supported if any declared reader supports it and as rejected only if all reject. `place` uses the tool's rule.
- **Stages inside the model.** `place` read one gap, from the input to the stage the head reads. A model with an
  encoder per modality and a fusion step has more than one. `place` now reads the loss step by step along any stages
  a run records, and an `upstream` item says between which two.
- **Too few events to pool.** Fewer than three admitted events are read one by one and resampled by cell
  (`target_by_event`), as the round under package 2 did by hand.
- **Words.** The two places where an item was worded for the first system take their words from the run.

## 6. Found and not changed

Proposals for the next iteration. Each changes how a reading of real cells would be used, so none is applied here.

| | what was seen | proposal |
|---|---|---|
| 1 | gate 3 as amended fails where the output overshoots on its own cells (first outer step) | one-sided: only a shortfall counts |
| 2 | a view that reads no more than its permuted copy still needs a positive control (first outer step) | such a view needs none: it reads nothing, and that is the finding |
| 3 | the reference stage has two reader classes, the head's repair one (first outer step) | a reader class at Z that is not linear, as a repair |
| 4 | the dominated-view check fired between two stages that both read nothing (planted, large, seed 2) | a view is dominated only if the later stage has a better verdict or reads a step more |
| 5 | a reader that failed on 67 cells was read as an input that lacks the response (planted, small, seed 3) | the item carries the number of cells, and the same reader on half the cells is a control before `unreadable` is used |

## 7. What this does not show

- A planted system shows that the rules are consistent with one another and with a known answer. It was written by
  the author of the rules. It is linear, its factors are independent, and its losses are complete (a stage carries a
  factor or does not). Real losses are partial.
- In the planted system a readout fitted on control cells already closes the head's loss; the rows of events add no
  step. The case where rows add a step has only been met on real cells.
- Three readings of the routing table have still been produced only by constructed targets: no rows, confounded,
  guard conflict.
- The comparison on the first system shows that the general code is the same computation as the scripts it replaces.
  It is not a second system.
- Events without control cells of their own (donor, batch, study, disease) are read by the audit tool and are not yet
  handled by the round. The other limits are in [`GENERAL.md`](../../GENERAL.md), Section 6.

## 8. Run it

    python -m anm_diag.toy                                           # one planted system, printed
    python -m pytest -q tests/test_cycle.py tests/test_outer.py      # 13 and 45 tests
    python scripts/planted.py --seeds 12 --part 0 --of 2 --out results/planted_0.json
    python scripts/planted.py --seeds 12 --part 1 --of 2 --out results/planted_1.json
    python scripts/planted.py --merge results/planted_0.json results/planted_1.json      # writes results/planted.json
    python scripts/regress.py --observe OBSERVE1 --observe2 OBSERVE2 --work2 WORK2        # folders and environment as outer_step_20261122

Results: `results/planted.json` (every run with its statuses and items), `results/regression.json`.
