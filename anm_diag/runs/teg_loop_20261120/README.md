# The repair by rows as one round of the loop protocol

> **Correction, 20 November, after a check of the design that follows this round**
> (proposed the same day, as the design of the outer loop). Below, this page says that in co-culture what
> remains "is not in RNA either" and is "not readable from RNA". For PD-L1 that does not stand. The RNA view used here
> (512 components, a linear reader) leaves 0.19 where Z, which is computed from the same RNA, leaves 0.10 with the
> event's own cells, and the gap is resolved below zero (−0.09 [−0.13, −0.04]): the RNA view is the weaker reader, and
> its rejection says nothing about RNA. Run 5, sealed, on the same cells, had the opposite order (RNA 0.12, Z 0.19).
> What stands for PD-L1 in co-culture: no reader of RNA or of Z that was tried comes inside the margin; the lowest
> estimate is about 0.10 [0.01, 0.17]; whether anything is missing from RNA is not decided. For MHC-I in co-culture
> the two views agree (0.31 at both, each rejected) but only a linear reader was used: a third of the change is read
> by no linear reader of RNA or of Z, and a richer reader comes before any conclusion about the input. Also found by
> that check: the response of an event (A1) was screened on the cells that were then audited. No number below depends
> on it, since no event near the threshold was admitted. The numbers below are unchanged.

Status: **exploratory**, as [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md) grades every round without final data. The design
is [`DESIGN.md`](DESIGN.md). Part A (PD-L1) was development on cells already read. Part B (MHC-I, HLA-E) was designed
first (commit `e68603d`) and run once.

The three repair runs of 19 November computed a repair from "rows of the table" without the protocol of this repository
or the source chain of the paper: nothing was admitted, an event carried neither its condition nor its type, and closure
was a rule written for the occasion. This run does the same repair as one round of the protocol: the grid of
`anm_diag scope`, typed events through propose, prepare, admit and assign, a gate that gives every cell one type (cell
type and condition), closure by the protocol's grades, one ledger.

*u* is the share of a measured change that an output leaves unexplained; the margin and the step margin are 0.10. "Added
beyond Z" is the difference between *u* of the head's own output and *u* of a readout of Z learned from the other genes
of the same group: what the head still loses of what its input carries.

## Outcome

**The head's part of the loop closes wherever an event is admitted. The whole chain does not: what remains lies before
the head, and in the co-culture condition it is not in RNA either.**

- **Part A, PD-L1, three groups.** After the repair the head's own output adds no loss beyond Z in any group (0.01, 0.00,
  0.03, each interval below 0.10). Two groups reach the grade that admission called attainable. The IFN-γ group, the
  only one where "absolute" was attainable, misses it by 0.005 on the upper bound of its interval.
- **Part B, by the rule written first: the loop does not close.** MHC-I has one target (melanoma in co-culture). Its
  repair did not pass gate 3 as written (on the cells it was fitted on it still leaves 0.28), so its closure was not
  graded, and its hypothesis had no control arm to be tested against. HLA-E has no target: no event passes admission.
- Beside the rule, for that one target: the repaired output adds 0.07 [0.05, 0.09] beyond Z, below the step margin, as
  in Part A. Z itself leaves 0.32 of the change of MHC-I, and so does RNA (0.31). A third of that change is not readable
  from RNA by these readouts, so no repair of the head or of the encoder can supply it.

## Part A: PD-L1 (development)

| group | `units` at Z | attainable | own output as trained | with the repair by rows | added beyond Z | closure | test of the hypothesis |
|---|---|---|---|---|---|---|---|
| melanoma, IFN-γ | 0.03 [−0.04, 0.10] support | absolute | 0.91 | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] no added loss | **relative** | confirmed |
| melanoma, co-culture | 0.12 [0.01, 0.19] unresolved | relative | 1.00 | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] no added loss | **relative** | confirmed |
| THP-1 | 0.13 [0.02, 0.31] unresolved | relative | 0.97 | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] no added loss | **relative** | unresolved |

**melanoma, IFN-γ**

| output | own-output *u* | gain over the head as trained | added beyond Z | by gene |
|---|---|---|---|---|
| the head as trained | 0.91 [0.89, 0.92] reject | | 0.87 | |
| the repair of the other group, as it is | 0.28 [0.20, 0.34] reject | 0.63 [0.56, 0.70] | 0.24 [0.16, 0.28] loss | IFNGR1 0.28, IFNGR2 0.25, JAK1 0.16, JAK2 0.36, STAT1 0.32 |
| this group's control cells only | 0.13 [0.07, 0.18] unresolved | 0.78 [0.72, 0.84] | 0.09 [0.03, 0.12] unresolved | IFNGR1 0.10, IFNGR2 0.12, JAK1 0.07, JAK2 0.21, STAT1 0.13 |
| this group's controls, the other group's rows | 0.03 [−0.04, 0.09] support | 0.88 [0.81, 0.91] | −0.01 [−0.03, 0.02] no added loss | IFNGR1 0.00, IFNGR2 0.02, JAK1 −0.06, JAK2 0.11, STAT1 0.03 |
| **this group's controls and admitted rows** | 0.04 [−0.02, 0.11] unresolved | 0.86 [0.80, 0.91] | 0.01 [−0.03, 0.03] no added loss | IFNGR1 0.02, IFNGR2 0.04, JAK1 −0.04, JAK2 0.13, STAT1 0.04 |
| this group's controls and every knockout cell (reference) | 0.01 [−0.06, 0.07] support | 0.90 [0.83, 0.91] | −0.03 [−0.04, 0.03] no added loss | IFNGR1 −0.01, IFNGR2 0.00, JAK1 −0.09, JAK2 0.09, STAT1 0.02 |

**melanoma, co-culture**

| output | own-output *u* | gain over the head as trained | added beyond Z | by gene |
|---|---|---|---|---|
| the head as trained | 1.00 [0.99, 1.00] reject | | 0.88 | |
| the repair of the other group, as it is | 0.35 [0.26, 0.41] reject | 0.65 [0.58, 0.73] | 0.23 [0.19, 0.28] loss | IFNGR1 0.37, IFNGR2 0.38, JAK1 0.18, JAK2 0.36, STAT1 0.41 |
| this group's control cells only | 0.28 [0.20, 0.34] reject | 0.71 [0.66, 0.79] | 0.16 [0.14, 0.20] loss | IFNGR1 0.28, IFNGR2 0.29, JAK1 0.14, JAK2 0.30, STAT1 0.35 |
| this group's controls, the other group's rows | 0.17 [0.07, 0.23] unresolved | 0.83 [0.76, 0.92] | 0.05 [0.03, 0.07] no added loss | IFNGR1 0.17, IFNGR2 0.18, JAK1 −0.01, JAK2 0.19, STAT1 0.25 |
| **this group's controls and admitted rows** | 0.12 [0.02, 0.19] unresolved | 0.88 [0.81, 0.97] | 0.00 [−0.02, 0.02] no added loss | IFNGR1 0.12, IFNGR2 0.13, JAK1 −0.06, JAK2 0.14, STAT1 0.20 |
| this group's controls and every knockout cell (reference) | 0.12 [0.03, 0.18] unresolved | 0.88 [0.81, 0.97] | 0.00 [−0.02, 0.02] no added loss | IFNGR1 0.11, IFNGR2 0.13, JAK1 −0.05, JAK2 0.15, STAT1 0.19 |

**THP-1**

| output | own-output *u* | gain over the head as trained | added beyond Z | by gene |
|---|---|---|---|---|
| the head as trained | 0.97 [0.95, 0.98] reject | | 0.83 | |
| the repair of the other group, as it is | 0.76 [0.68, 0.82] reject | 0.21 [0.14, 0.28] | 0.62 [0.48, 0.71] loss | IFNGR1 0.75, IFNGR2 0.77, JAK2 0.76, STAT1 0.54, IRF1 0.90 |
| this group's control cells only | 0.53 [0.47, 0.66] reject | 0.44 [0.30, 0.50] | 0.40 [0.32, 0.48] loss | IFNGR1 0.47, IFNGR2 0.49, JAK2 0.52, STAT1 0.50, IRF1 0.90 |
| this group's controls, the other group's rows | 0.27 [0.16, 0.41] reject | 0.69 [0.55, 0.80] | 0.14 [0.09, 0.18] unresolved | IFNGR1 0.24, IFNGR2 0.26, JAK2 0.29, STAT1 −0.04, IRF1 0.70 |
| **this group's controls and admitted rows** | 0.16 [0.08, 0.29] unresolved | 0.80 [0.67, 0.88] | 0.03 [−0.03, 0.08] no added loss | IFNGR1 0.10, IFNGR2 0.14, JAK2 0.17, STAT1 0.02, IRF1 0.54 |
| this group's controls and every knockout cell (reference) | 0.14 [0.08, 0.25] unresolved | 0.82 [0.71, 0.89] | 0.01 [−0.07, 0.09] no added loss | IFNGR1 0.10, IFNGR2 0.12, JAK2 0.13, STAT1 0.10, IRF1 0.43 |

Three things in these tables:

- **Rows and cells agree again.** In each group the repair by rows and the reference that sees every knockout cell
  read the same within their intervals.
- **The condition matters for the baseline more than for the rows.** The repair of the other condition, applied as it
  is, leaves a loss beyond Z (0.24 and 0.23). With this group's own control cells and the *other* condition's rows, the
  loss added is below the step margin in both directions (−0.01 and 0.05). For THP-1 the rows of melanoma do less
  (0.14, unresolved). The field "condition" could be dropped for the rows of one cell type here; it cannot be dropped
  for the baseline. This is the statement of SI S1.4.1, read off a table.
- **THP-1 is held back by one event.** IRF1 stays at 0.54 with the rows and 0.43 with its own cells.

**The grid on the head as trained, PD-L1** (pooled events; `frozen` is the head's own output)

| group | response (control SD) | frozen | train | other | units | own | reading of the grid |
|---|---|---|---|---|---|---|---|
| melanoma, IFN-γ | −1.14 | 0.91 [0.89, 0.92] reject | 0.99 [0.99, 1.00] reject | 0.32 [0.24, 0.38] reject | 0.03 [−0.04, 0.10] support | 0.03 [−0.05, 0.09] support | coverage |
| melanoma, co-culture | −1.12 | 1.00 [0.99, 1.00] reject | 0.99 [0.98, 1.00] reject | 0.29 [0.20, 0.35] reject | 0.12 [0.01, 0.19] unresolved | 0.10 [0.01, 0.17] unresolved | unresolved |
| THP-1 | −0.87 | 0.97 [0.95, 0.98] reject | 0.96 [0.94, 0.99] reject | 0.72 [0.64, 0.79] reject | 0.13 [0.02, 0.31] unresolved | 0.11 [0.02, 0.20] unresolved | unresolved |

**The input view** (the same events read from RNA, 512 components)

| group | from RNA, other genes | from RNA, own cells | added from RNA to Z (own cells) | what Z leaves is |
|---|---|---|---|---|
| melanoma, IFN-γ | 0.08 [0.00, 0.14] unresolved | 0.09 [0.02, 0.14] unresolved | −0.06 [−0.09, 0.02] no added loss | nothing: Z is supported |
| melanoma, co-culture | 0.18 [0.08, 0.25] unresolved | 0.19 [0.11, 0.25] reject | −0.09 [−0.13, −0.04] no added loss | the `units` scope at Z is unresolved. From RNA, a readout with the gene's own cells is reject; from RNA to Z: no added loss. Read as not readable from RNA by these readouts. |
| THP-1 | 0.14 [0.03, 0.35] unresolved | 0.09 [0.02, 0.14] unresolved | 0.02 [−0.04, 0.10] unresolved | the `units` scope at Z is unresolved. From RNA, a readout with the gene's own cells is unresolved; from RNA to Z: unresolved. Not placed: neither view reaches a verdict. |

Where Z is not supported, the input view places what is left. In co-culture a readout of RNA with the gene's own cells
leaves 0.19 (rejected), more than Z leaves: what remains is not readable from RNA by these readouts, and it is not lost
by TEDDY. For THP-1 neither view reaches a verdict. The RNA view is 512 principal components; a richer readout of RNA
could read more.

**Gates 3 and 6, and the gate of the system**

| group | on the cells it was fitted on | guard, validation cells (fall) | guard, test cells: readout, mean, worst protein | own class given to control cells | to knockout cells |
|---|---|---|---|---|---|
| melanoma, IFN-γ | −0.05 | 0.0000 | 0.0003, 0.00000, 0.0003 | 0.995 | 0.956 |
| melanoma, co-culture | −0.01 | 0.0000 | 0.0003, 0.00000, 0.0003 | 0.988 | 0.993 |
| THP-1 | 0.06 | 0.0000 | 0.0003, 0.00000, 0.0003 | 0.999 | 0.985 |

Bone-marrow test cells given to a group by the gate: 1 of 5,000.

**With the readouts mixed by the gate's probabilities** (`--soft-gate`, a sensitivity)

| group | with the repair by rows | on the cells it was fitted on | closure |
|---|---|---|---|
| melanoma, IFN-γ | 0.03 [−0.03, 0.10] support | −0.06 | absolute |
| melanoma, co-culture | 0.12 [0.02, 0.19] unresolved | −0.01 | relative |
| THP-1 | 0.29 [0.23, 0.39] reject | 0.21 | not tested (gate 3) |

The gate assigns each cell to one class. Mixing the readouts by the gate's probabilities puts a few percent of
another group's readout, evaluated far from its own cells, into the output: it lifts the IFN-γ group over the line
(0.03, supported) and spoils THP-1 (0.29, gate 3 not passed). The choice of one class per cell was made on these cells,
after seeing both. That the IFN-γ group reads "absolute" with one gate and "relative" with the other is what "missed by
0.005" means: it is on the line.

About 4% of the IFN-γ knockout cells are given to the co-culture class. A knockout of the pathway makes a cell look
less like its condition: a type read from the embedding is not independent of the event. In the table of the paper a
type comes from the record, not from the state.

Probe CD274 (response -1.99 control SD): frozen 1.00 [0.99, 1.01] reject, train 1.00 [1.00, 1.00] reject, units 0.92 [0.87, 0.96] reject, own 0.91 [0.87, 0.95] reject. not admitted: A2: the event acts on the own gene of readout PDL1 (rule S3).

An event on the readout's own gene is stopped by admission before any repair. The grid agrees: no readout of Z carries
it, with the other genes or with its own cells.

Forecasts of Part A (the grid's ceilings; retrospective, the outcomes were known):

| hypothesis | arm | forecast | measured gain | inside |
|---|---|---|---|---|
| coverage:mel_ifn:PDL1:development | other | [0.52, 0.66] | 0.63 | yes |
| coverage:mel_ifn:PDL1:development | rows | [0.80, 0.91] | 0.86 | yes |
| coverage:mel_co:PDL1:development | other | [0.65, 0.79] | 0.65 | no |
| coverage:mel_co:PDL1:development | rows | [0.80, 0.98] | 0.88 | yes |
| coverage:thp1:PDL1:development | other | [0.17, 0.32] | 0.21 | yes |
| coverage:thp1:PDL1:development | rows | [0.65, 0.94] | 0.80 | yes |

## Part B: MHC-I and HLA-E (designed first, run once)

**Admission**

| readout | group | pooled response (control SD) | times its half-width | events that pass A1 on their own | `units` at Z | admission |
|---|---|---|---|---|---|---|
| MHC1 | melanoma, IFN-γ | −0.27 | 3.1 | 0 of 5 | −0.18 [−0.58, 0.08] unresolved | not a target: no event of the group passes A1 on its own: there is no row to repair with |
| MHC1 | melanoma, co-culture | −1.97 | 16.5 | 5 of 5 | 0.32 [0.27, 0.37] reject | **target**, attainable relative |
| HLAE | melanoma, IFN-γ | −0.14 | 2.2 | 0 of 5 | 0.24 [−0.36, 0.52] unresolved | not a target: A1: the response is small against its interval (|change| / half-width = 2.2, wanted 3) |
| HLAE | melanoma, co-culture | −0.16 | 3.1 | 0 of 5 | 0.43 [0.04, 0.61] unresolved | not a target: no event of the group passes A1 on its own: there is no row to repair with |

MHC-I in co-culture, the one target:

**melanoma, co-culture**

| output | own-output *u* | gain over the head as trained | added beyond Z | by gene |
|---|---|---|---|---|
| the head as trained | 0.97 [0.96, 0.98] reject | | 0.65 | |
| this group's control cells only | 0.67 [0.64, 0.69] reject | 0.30 [0.28, 0.33] | 0.35 [0.31, 0.38] loss | IFNGR1 0.65, IFNGR2 0.68, JAK1 0.66, JAK2 0.65, STAT1 0.68 |
| **this group's controls and admitted rows** | 0.39 [0.34, 0.44] reject | 0.58 [0.53, 0.62] | 0.07 [0.05, 0.09] no added loss | IFNGR1 0.33, IFNGR2 0.40, JAK1 0.40, JAK2 0.37, STAT1 0.43 |
| this group's controls and every knockout cell (reference) | 0.37 [0.32, 0.42] reject | 0.60 [0.55, 0.65] | 0.05 [0.03, 0.07] no added loss | IFNGR1 0.31, IFNGR2 0.37, JAK1 0.38, JAK2 0.35, STAT1 0.40 |

- **Gate 3.** On the cells it was fitted on the repair leaves 0.28. The check in the code asks for less than the margin
  there, so the change counts as not tested and its closure is not graded.
- **The hypothesis** had no control arm: the other group has no admitted event, so there was no repair of another group
  to apply. `loop.ledger_check` reports it as an error (H3).
- **After the run**, the same events from RNA: other genes 0.30 [0.25, 0.35], own cells 0.31 [0.27, 0.36], both rejected;
  from RNA to Z nothing is added (0.00 [−0.04, 0.03]).

Against the forecast:

| | forecast | read |
|---|---|---|
| MHC-I, co-culture | admitted; attainable relative; no loss added beyond Z (80%); hypothesis confirmed (65%) | admitted, relative attainable, 0.07 [0.05, 0.09] added: as forecast. **Not foreseen: gate 3 not passed, and no control arm, so no closure and no test** |
| HLA-E, co-culture | admitted with probability one half | not a target: no single event passes A1 |
| MHC-I and HLA-E, IFN-γ | more likely not admitted (60%) | not targets |
| guard | kept | kept |
| the loop closes on both readouts | 55% | **it does not close on either** |

Two things about Part B that were not in the design:

- The first call stopped on an error before any output was computed or shown: a group was ready as a pool while none of
  its events passed A1 on its own. The design did not name this case. It was resolved as "not a target" (commit
  `2a93211`), and the run was then made.
- The check of gate 3 is written for a readout that Z carries in full. Where Z itself leaves a third, a repair cannot
  reach the margin even on its own cells. A check relative to Z would have passed here; changing it now would be
  choosing the rule after the outcome. It is listed below as something to fix before the next round.

## What the round shows, and what it does not

For one head, three groups of cells and five genes of one pathway per group, on cells that had all been read before:

- Run under the protocol, the repair by rows removes the head's loss in every admitted case: four of four outputs add
  nothing beyond Z after it (0.01, 0.00, 0.03, 0.07).
- Admission does what the test of 19 November lacked. It says beforehand which grade can be reached, stops the event on
  the readout's own gene, and keeps readouts with a response too small to follow out of the round.
- What is left after the repair is placed: supported at Z in one group; not readable from RNA in the two co-culture
  cases; not placed for THP-1.
- "The whole loop closed" cannot be said. The one group where an absolute pass was attainable is on the line, and the
  replication did not close by its rule.

It does not show a confirmation (no final data; Part B shares its cells with Part A), anything about another head,
pathway or cell type, or that the gate holds on cells of a kind it has not seen.

## Before a next round

- Gate 3 stated relative to what Z carries, and a rule for a hypothesis whose control arm cannot be run.
- The type of a cell taken from its record where there is one; the gate only where there is none.
- For what RNA does not carry (PD-L1 and MHC-I in co-culture, CD274): another measurement or input. This belongs to the
  part of the loop that asks what to measure, which is not designed.
- Final data. The melanoma cells with embeddings are all read; the rest of the THP-1 screen needs embeddings.
- The repair as a module of the package with its own tests, once a round has closed on final data.

## Files

- [`DESIGN.md`](DESIGN.md): the round, Part A as it stood, the rule and forecast of Part B.
- `scripts/build.py`, `build_rna.py`: arrays and SPEC files for `anm_diag scope` (the grids themselves are not in git).
  `scripts/round.py`: admission, the repair, the reading, the ledger. `scripts/after.py`: MHC-I from RNA, after the run.
  `scripts/tables.py` prints the tables above.
- `results/round_development.json`, `round_development_soft.json`, `ledger_development.jsonl`: Part A.
  `results/round_replication.json`, `ledger_replication.jsonl`, `after_replication.json`: Part B.
