# The repair by rows as one round of the loop protocol. Design

Status: **exploratory**, as [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md) grades every round without final data. Two parts:

- **Part A, development** (already run when this was written): PD-L1 on three groups of cells. Every cell and this
  readout had been read before. Its hypotheses are marked retrospective.
- **Part B, replication** (not yet run when this was committed): the same code, unchanged, on two outputs of the head
  that no repair by rows has touched, MHC-I and HLA-E. The rule and my forecast are below. One run.

## Why this run

The three runs of 19 November that computed a repair from "rows of the table"
([`teg_training`](../teg_training_20261119/), [`teg_repair_gate`](../teg_repair_gate_20261119/),
[`teg_repair_test`](../teg_repair_test_20261119/)) did not go through the protocol of this repository or the source
chain of the paper. They imported neither `anm_diag/events.py` nor `anm_diag/loop.py`. A "row" there was three numbers
computed in a script. Nothing was admitted, no event carried its condition or its type, and closure was a rule written
for the occasion. The test then failed in the places those steps cover: it asked for an absolute pass where the
embedding itself was not supported; it carried rows across conditions without the condition being a field; and it
probed an event on the readout's own gene.

This run does the same repair as one round of the protocol.

## The round

A **group** is a cell type in a culture condition: melanoma in the IFN-γ condition (`mel_ifn`), melanoma in co-culture
(`mel_co`), THP-1 (`thp1`). Every event of a group shares that type.

1. **Observe.** `anm_diag scope` on the head as trained, per group and readout (`scripts/build.py` writes the SPEC
   files): the head's own output (`frozen`, on the log scale), and readouts of Z learned from the head's training
   cells (`train`), from another group (`other`), from the other genes of the group (`units`) and with the gene's own
   cells (`own`). Where it exists, the same from RNA (`scripts/build_rna.py`), as the input view.
2. **Admit**, by the source chain of the paper (SI S1.4):
   - *propose*: one typed event per gene, written by `events.source_events` (time, site, channel, sign, amplitude,
     stratum, provenance, pool, no-leakage entry);
   - *prepare*: only the half of the cells that is fitted; no measurement of the half that is read;
   - *admit*: A1 the response is resolved and at least three times the half-width of its interval; A2 the event does
     not act on the readout's own gene (rule S3); A3 the head loses something that cells of this kind read (the gap
     from the head's own output to the `units` scope at Z lies above the step margin); A4 an event enters the repair of
     its own group only, and a pool across conditions is flagged by `events.closure`;
   - *assign*: the displacement in Z, the measured change and its standard error, computed in the fit.
   The grade that can be reached is fixed here: **absolute** if the `units` scope at Z is supported, otherwise
   **relative** (the output can be brought to what Z carries, not further). What Z leaves is written as an open item,
   placed by the input view: upstream if RNA carries it and Z does not, not readable if RNA does not either.
3. **Hypothesize.** Claim: coverage (it takes cells of this kind). Forecast for the repair by rows: the grid's ceiling
   (the head's own output against the `units` scope). Rule: the repair by rows comes within the step margin of the
   ceiling, and the repair of another group, applied here as it is, does not.
4. **Intervene.** One change. A gate over bone marrow and the groups assigns each cell to one class. Behind the gate
   of a group, the head's output is replaced by a linear readout of Z fitted on the group's control cells and updated
   once by its admitted events (the update of `teg_repair_test`; ridge weight by five-fold cross-validation on the
   control cells, as there). Bone-marrow cells and groups without a repair keep the head's output. Two halves: every
   cell is read by a system fitted without it. Gate 3: the repair reaches its objective on the cells it was fitted on,
   no cell of a reading is in its fit, the guard holds on bone-marrow validation cells.
5. **Test and update.** The head's own output again, on the same resamples as the grid. Closure of the own output:
   - **absolute**: supported (two-sided, margin 0.10), and the output as trained is rejected;
   - **relative**: no loss added beyond the `units` scope at Z (interval of the difference below the step margin 0.10),
     and the output as trained adds one;
   - **loss no longer resolved**, **not closed**: otherwise, by that interval.
   Gate 6: the guard on bone-marrow test cells (the readout's own score, the mean over proteins, the worst single
   protein). Every event, admission, change, hypothesis, test and closure is a line of one ledger, checked by
   `loop.ledger_check`.

Changes against `teg_repair_test`, all part of the one change of this round: admission; one system for all groups
with a gate typed by cell type and condition; closure read by the protocol's grades.

## Part A: development on PD-L1

All three groups are admitted. `scripts/round.py --readouts PDL1 --tag development`; numbers in
`results/round_development.json`.

| group | `units` at Z (grid) | attainable | own output as trained | with the repair by rows | added beyond Z | closure | test of the hypothesis |
|---|---|---|---|---|---|---|---|
| melanoma, IFN-γ | 0.03 support | absolute | 0.91 | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] none | **relative** | confirmed |
| melanoma, co-culture | 0.12 unresolved | relative | 1.00 | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] none | **relative** | confirmed |
| THP-1 | 0.13 unresolved | relative | 0.97 | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] none | **relative** | unresolved |

- In every group the repaired output adds no loss beyond what Z carries. Two groups reach the grade that was
  attainable. The IFN-γ group misses "absolute" by 0.005 on the upper bound of its interval.
- Developed here, on cells already read: the gate assigns one class per cell. With the readouts mixed by the gate's
  probabilities (`--soft-gate`) the IFN-γ group reads 0.03 [−0.03, 0.10], supported, and THP-1 reads 0.29 and fails
  gate 3: a few percent of another group's readout, evaluated far from its own cells, is enough to spoil it. One class
  per cell is kept for Part B. The rule of the hypothesis was first written as "two thirds of the ceiling" and then
  changed to the step margin, after the development numbers were seen.
- The input view: in co-culture a readout of RNA with the gene's own cells leaves 0.19 (rejected), more than Z leaves.
  What remains there is not readable from RNA by these readouts. It is not lost by TEDDY.
- The probe, CD274: not admitted (A2). No scope reads it at Z (0.92 with the other genes, 0.91 with its own cells).

## Part B: replication on MHC-I and HLA-E

Groups `mel_ifn` and `mel_co` (THP-1 has no such measurement). The head's outputs HLA-A-B-C and HLA-E. Nothing of the
code changes: `scripts/round.py --readouts MHC1 HLAE --tag replication`, which refuses to run twice. The grid uses the
ridge class only (declared in the SPEC files).

What had been seen of these readouts: earlier rounds of the loop retrained heads on them in the IFN-γ condition and
found the pooled response too small to follow (`LOOP_PROTOCOL.md`); the census of the screen showed a large change of
MHC-I in co-culture. No repair by rows was computed for them and no grid of this kind was read. The cells are those of
Part A.

**Rule.** Targets are the groups that admission calls ready. The loop closes on a readout if every target reaches the
grade admission called attainable, the guard is kept, and the ledger has no error. A group that is not admitted is not
a target and is not counted, as the protocol says.

**Forecast**, written before the grids of these readouts were run:

| | forecast |
|---|---|
| MHC-I, co-culture | admitted; attainable relative; the repair by rows adds no loss beyond Z (80%); hypothesis confirmed (65%) |
| HLA-E, co-culture | admitted with probability one half; if admitted, relative closure (70%) |
| MHC-I and HLA-E, IFN-γ | more likely not admitted (60%): the response is small against its interval. If admitted, relative closure (60%) |
| guard | kept in every case |
| the loop closes on both readouts | 55% |

## What this run cannot show

- A confirmation. There are no final data: every cell here has been read before, and Part B shares its cells with
  Part A. A round that closes here is a round that closes in development.
- Anything beyond one head, three groups of cells, five genes of one pathway per group.
- That the gate would hold on cells of a kind it has not seen.
- With five genes as resampling units the pooled intervals are approximate.
