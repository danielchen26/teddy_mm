# Can a finer type cover a donor that has no record of its own? Declaration

Written on 24 November 2026, **after the script was tried on made-up values and before it was run on the donors**.

Everything here is exploratory. The donors are the eight development donors of the round of 23 November, whose cells
earlier runs have read. No reserve donor is opened, no final unit is read and no model is trained.

## 1. Why this is asked

The loop flags an event that differs from its type. For donors it did so ten times: the type table reads the eight
donors as one type, and for ten pairs of a donor and a protein the rows of the other donors do not reach the donor.
Its own row is needed ([`unit_events_20261123`](../unit_events_20261123/)).

What follows a flag has not been done ([`GENERAL.md`](../../GENERAL.md), Section 6): proposing a type **with several
members**, so that the records of the other members cover a member that has none, and testing it. This run does it
for the donors, where members are many. It is the first time a proposed type is put to its test.

## 2. The data

The cells of the round of 23 November that lie outside its selection third: 26,706 cells of eight donors, with the
state the head reads, 112 measured proteins and the head's own outputs. The pairs are the 26 that round admitted
(one donor per protein). Three sets of them:

| set | pairs | |
|---|---|---|
| `own_row_needed` | 10 | the head has a step to remove, and the donor's own row reads a step more than the rows of the others: the ten flagged items |
| `a_step_to_remove` | 16 | the head has a step to remove |
| `admitted` | 26 | all |

The anchors, from that round, for the ten: with the rows of the other seven donors and no record of the donor
(`unseen`) the output leaves a median of 0.75; with the donor's own row (`rows`) 0.45.

## 3. What a type is here, and the candidates

A type is a grouping of the eight donors into two groups of four. There are 35 such groupings.

- **From the record: the batch.** The only field of the record beyond the donor. Donors 1 to 4 are one batch, 5 to 8
  the other.
- **From the data, with no protein: the state.** Of the 35 groupings, the one that best separates the donors' mean
  states (the squared distance between the two group means over the mean squared distance of a donor from its group
  mean). It is found first and compared with the record after.

The other groupings are the control: **a grouping that is a type should do better than groupings that are not**.
With 35 of them the comparison is exact, and coarse: the best possible rank is 1 of 35.

## 4. Part 1: the record alone

For a pair of a donor U and a protein, U's type-mates are the three other donors of its group.

    u_type = 1 − (the mean of the mates' cells − the mean of the seven other donors' cells)
                 / (the mean of U's cells − the mean of the seven other donors' cells)

No state and no readout enter: this asks only whether the members of a type shift alike. 0: the mates' record says
exactly what U measured. 1: it says nothing. Cells are resampled within a donor (1,000 resamples, the same for every
protein). Supported, unresolved and rejected are read at the margin 0.10, as everywhere.

## 5. Part 2: with the repair

The round's repair for units (`anm_diag.cycle.UnitRows`), with one more arm, `typed`: a fit with **no protein of U**,
in which the readout of the variation among cells is that of the other seven donors (as in `unseen`) and the rows
are those of every donor against the others, **with the mean of U set to the mean of its type-mates**. So U has a
row, and it is its type's, not its own. It is read as the round reads an output (its cells by cell, the other donors
by unit, the halves and resamples of the round), beside `unseen` and `rows`, on paired resamples.

## 6. The rule

For each candidate, on the ten pairs of `own_row_needed`:

- **H1, the type's record covers the member.** Confirmed: `u_type` is supported for at least half of the pairs.
  Refuted: rejected for more than half. Unresolved otherwise.
- **H2, the repair with the type's record gains a step.** Confirmed: `typed` gains a step over `unseen` for at least
  half of the pairs. Refuted: no step for more than half. Unresolved otherwise.
- **Not a random grouping.** The candidate's median |u| is first or second of the 35 groupings (by the record, and
  with the repair).

Also reported: the same for the other two sets; for every pair whether the donor's own row still adds a step over
`typed`; the best of the 35 groupings and the median of the 35.

**Checks.** The arms `rows` and `unseen`, computed again here, repeat the round of 23 November for every pair (0.005).
The arm `unseen` as computed for `typed` equals the package's. If a check fails: stop and report.

What each outcome would mean:

| outcome | what follows |
|---|---|
| H1 and H2 confirmed, not a random grouping | a type with four members covers a member without a record. It is then a proposal for the schema, to be admitted on donors that were not used: the reserve donors, which are the owner's to open |
| H1 confirmed, H2 not | the members shift alike and the repair does not use it: the repair, not the type, is what is short |
| neither, and no grouping does much better | for these proteins a donor differs from every other donor on its own. No grouping of eight donors is a type that covers them; a donor needs its own record, and that is a limit to state, or a field the record does not hold |

## 7. What was done before this was written

`scripts/type_test.py` was tried on two sets of made-up donors in a separate folder. In one, three proteins shift by
batch only and three by one donor only: the record of the batch explains the first three (*u* 0.14 to 0.32) and none
of the others (1.0), and the batch is first of the 35 groupings. In the other every donor also differs on its own,
and no grouping covers a member. Nothing was changed after the trials except the name of the script.

## 8. The forecasts

`results/declared.json`, 12 of them. Briefly, and why:

- **H1 is more likely refuted** (0.60) for the batch. The admitted pairs are not spread evenly: 23 of 26 are donors
  of one batch, and three donors hold 21 of them with proteins that go together (for one donor the B-cell markers).
  That looks like donors that differ on their own, not like a batch.
- **The batch is still likely to be among the best groupings** (first or second: 0.50): what little is shared is
  probably shared within a batch.
- **H2: refuted** (0.65).
- The grouping found in the states is the batch (0.55).

## 9. What this run cannot say

- Eight donors, one study, one field of the record. A type that fails here is this type on these donors.
- A grouping into two groups of four is the only shape tried. A type with two members, or a dimension that is not a
  grouping (a share of a cell type, a time), is not.
- Passing here would not admit a type: that needs donors that were not used.
