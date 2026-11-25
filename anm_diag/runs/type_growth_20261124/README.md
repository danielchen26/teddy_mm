# Can a finer type cover a donor that has no record of its own?

24 November 2026. The loop flags an event that differs from its type; for the eight donors of the round of
23 November it did so ten times ([`unit_events_20261123`](../unit_events_20261123/)). This run takes the step after
the flag, for the first time: it proposes a type **with several members** and tests whether the records of the other
members cover a member that has none. The candidates, the two readings, the rule, the control and 12 forecasts were
committed before the script was run on the donors ([`DESIGN.md`](DESIGN.md)).

Everything here is exploratory. The donors are development donors that earlier runs have read. No reserve donor was
opened, no final unit was read and no model was trained.

## What came out

**No grouping of these donors is a type that covers a member. The batch is the best of the 35 groupings, and the
states find it without any protein; but a donor's three batch-mates say about a fifth of what the donor measured,
and the repair gains nothing from it. For all ten flagged pairs the donor's own row still adds a step.** H1 and H2
are refuted for both candidates, which are the same grouping.

A type here is a grouping of the eight donors into two groups of four. *u* is the share of the measured change that
is left unexplained: 0 is none of it, 1 is all of it.

- **The candidate found in the data is the candidate from the record.** Of the 35 groupings, the one that best
  separates the donors' mean states is the batch (donors 1 to 4 against 5 to 8), by a wide margin: 8.6, where the
  median of the 35 is 0.5. No protein entered that.
- **Part 1, the record alone.** For the ten flagged pairs the type-mates' record leaves a median of 0.82
  [0.80, 0.84] of the donor's change: rejected for 10 of 10. The batch is first of the 35 groupings (their median is
  0.99, the worst 1.11), so what is shared is shared within a batch, and it is small. Over all 26 admitted pairs no
  pair is below 0.59.
- **Part 2, with the repair.** With the donor's mean set to its type-mates' mean (`typed`), the output leaves a
  median of 0.77. With no record of the donor and no type (`unseen`) it leaves 0.75; with the donor's own row
  (`rows`) 0.45. The gain of the median over `unseen` is −0.02 [−0.07, 0.03]: no step. No pair gains a step (no step
  for 9, unresolved for 1), and the own row adds a step over `typed` for 10 of 10. Among the 35 groupings the batch
  is second (0.767 against 0.757 for the first; the median is 0.87).

| pair | the type's record alone | `typed` | `unseen` | `rows` | the own row still adds |
|---|---|---|---|---|---|
| CD274, donor 8 | 0.78 [0.74, 0.82] | 0.87 | 0.93 | 0.58 | +0.29: a step |
| CD47, donor 6 | 0.73 [0.71, 0.76] | 0.77 | 0.80 | 0.54 | +0.23: a step |
| CD25, donor 6 | 1.00 [0.92, 1.08] | 0.79 | 0.65 | 0.46 | +0.33: a step |
| CD103, donor 6 | 0.95 [0.86, 1.02] | 1.01 | 1.07 | 0.66 | +0.35: a step |
| CD69, donor 8 | 1.15 [1.10, 1.20] | 0.88 | 0.75 | 0.45 | +0.43: a step |
| CD107a, donor 8 | 0.91 [0.86, 0.95] | 0.69 | 0.61 | 0.36 | +0.33: a step |
| integrin β7, donor 7 | 0.72 [0.66, 0.76] | 0.74 | 0.75 | 0.42 | +0.32: a step |
| CD38, donor 5 | 0.86 [0.83, 0.88] | 0.76 | 0.70 | 0.44 | +0.32: a step |
| CD63, donor 7 | 0.65 [0.62, 0.69] | 0.64 | 0.67 | 0.31 | +0.33: a step |
| CD49a, donor 7 | 0.72 [0.69, 0.74] | 0.75 | 0.80 | 0.38 | +0.37: a step |

The other two sets read the same way. For the 16 pairs with a step to remove: the record alone 0.76, rejected for
16 of 16; `typed` 0.71 against `unseen` 0.66 and `rows` 0.40. For all 26: 0.86, rejected for 26 of 26; `typed`
0.53 against 0.49 and 0.35 ([`results/types.json`](results/types.json)).

**The checks hold.** The arms `rows` and `unseen`, computed again here, repeat the round of 23 November for every
pair (largest difference 0.0). The arm `unseen` as computed for `typed` equals the package's.

## What it means

- **The step after a flag can be taken, and here it says no.** A type was proposed by a rule, compared with every
  other grouping of the same shape, and read with and without the repair. That is the test Section 8.3 of
  [`OUTER_LOOP.md`](../../OUTER_LOOP.md) asks of a split, made for the first time on real units.
- **The schema lever does not close these ten items.** What a donor shows on these proteins is the donor's own.
  The type the data points to, the batch, is real (the states separate by it and it is the best grouping) and it
  covers a fifth of the change. By the routing of the loop the items then go elsewhere: a donor needs its own
  record, which is a limit to state in the claim, or a field the record does not hold.
- **"Readable with a record, not without" is, for these donors, what the data say and not only a step that was
  missing.** That sentence was open as long as no proposed type had been tested. One kind of type has now been
  tested on one set of units.
- **What it does not close.** Two groups of four is one shape of type. A type with two members was not tried, and
  neither was a dimension that is not a grouping, such as the share of a kind of cell in a donor's sample. Donor 6
  is admitted for the B-cell markers together, which is what such a dimension would look like. And eight donors of
  one study are few.

## The forecasts

Intervals: 3 of 3 held the value. Outcomes: for 8 of 9 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)). The one miss: with the repair the batch is second of the 35
groupings, and I gave "first or second" 0.35.

## What was not as written

Nothing.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: the candidates, the two parts, the rule, the control, the forecasts |
| [`results/types.json`](results/types.json) | the reading: every pair, every one of the 35 groupings, the statements, the checks |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/type_test.py`, `scripts/score.py` | the reading; the score |

To run again (the work folder of `unit_events_20261123` is needed):

    python scripts/type_test.py --work WORK
    python scripts/score.py

## What is next

- A dimension that is not a grouping: for each donor, what its sample is made of, read from the states with no
  protein, as a candidate field for the same test.
- The same test where a type has many members and a record of its own: the knockouts, once a second event of
  IRF1's kind exists among the units that have not been read.
- Whether a donor needing its own record is stated as a limit of the claim is the owner's decision.
