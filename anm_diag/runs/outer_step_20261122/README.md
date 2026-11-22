# The first outer step: results

22 November 2026. The step was declared in [`DESIGN.md`](DESIGN.md), with eleven hypotheses and their forecasts
([`results/declared.md`](results/declared.md)), and that declaration was merged before any reading was made. This
page says what was read. Everything here is exploratory: the cells are development cells, and no final unit has been
read.

The step changed the package, not the model: the rules (protocol version 3), how a cell gets its type, the readers,
and which events are read. Then one round was run under the new package, on the one target that states at hand could
reach.

## What came out

1. **Nothing closed absolutely.** Seven targets: three with the model's part closed (relative), one open before the
   model, one not graded, two not targets. The step started from 21 open items and ends with 17, of which four are at
   the margin.
2. **PD-L1 under IFN-γ misses the margin by 0.001.** With the type of a cell taken from the record the repaired
   output reads 0.037 [−0.034, 0.101] (in the round: 0.043 [−0.021, 0.105]). On a selection third no further
   knockout of the library passes for PD-L1, so the pool cannot decide it. It ends at the margin.
3. **Under IFN-γ the rows were not what was missing.** A readout fitted on that group's control cells, with no
   knockout row, already removes 0.79 of the 0.91 the head loses. The rows add 0.080 [0.019, 0.085], less than the
   step margin. In co-culture, in THP-1 and for MHC-I the rows add a resolved step (0.17, 0.37, 0.28).
4. **One of my five amendments does not work as written.** Gate 3 against the reference stage passes for three
   targets and fails for PD-L1 in co-culture (−0.09 [−0.15, 0.13]): on the cells it was fitted on, the repaired
   output overshoots in some resamples, |*u*| counts an overshoot as a loss, and the interval rests on five genes. By
   the rule that target is not graded, although its audit half reads no added loss (−0.008 [−0.027, 0.010]). The
   forecast for this was wrong (0.95 that it would pass). An amendment for the next iteration is proposed below.
5. **The split of the THP-1 pool is not supported in the form proposed.** With its own row IRF1 leaves 0.68; with the
   pooled rows it leaves 0.54. Its loss is before the head: a linear reader of RNA with IRF1's own cells leaves 0.18,
   a reader of Z leaves 0.40, and the difference is 0.21 [0.12, 0.32]. This is the first item sent to the model, and
   to its encoder, not its head.
6. **In co-culture the richer readers changed what can be said.**
   PD-L1: boosting on the 512 components reads 0.11, the same as Z (difference −0.00 [−0.04, 0.03]). The RNA view is
   no longer dominated, and nothing is lost between RNA and Z. Nothing reaches support either: 0 of 60 readings.
   MHC-I: boosting on the 512 components leaves 0.19; every reader of Z leaves 0.31 or more. Z reads 0.12
   [0.09, 0.15] less than that reader of RNA. This is not drawn as "upstream": the interval does not clear the step
   margin, and no reading of that target is supported (0 of 60), so the positive control fails. It is the next thing
   to test.
7. **The new target is lost wholly before the model.** MHC-I in melanoma under IFN-γ, with B2M and TAPBP admitted on
   a third of their cells. The knockouts change MHC-I by −2.9 and −1.6 control SD. The head as trained reads none of
   it (1.00, 1.00). Neither does any reader of Z or of RNA with the event's own cells (B2M: 0.91 at Z, 0.99 at RNA;
   TAPBP: 1.00 and 1.02); a permuted view reads 1.00. The head has nothing to repair, and the repair by rows does
   nothing (0.99, 1.01). Both genes act on the measured complex itself: B2M is its light chain, TAPBP loads it. This
   is what the knockout of CD274 showed for PD-L1 (0.91).
8. **Final units exist now.** From the record: 615 units, of which the system has been read on 216. Of the 399
   unread ones, 133 are set aside as final. A first batch is listed for which states are to be computed: 8 units and
   control cells of two conditions, 6,759 cells.
9. **The forecasts.** 19 of 29 intervals held the value. For 37 of 45 outcomes the class given the highest
   probability happened. One outcome had been given no probability at all (the test for PD-L1 in co-culture came out
   "not testable", because of point 4).

## Targets, before and after

| target | after the round of 20 November | after this step | why |
|---|---|---|---|
| PD-L1, melanoma, IFN-γ | the model's part closed (relative) | the same; the absolute question is at the margin | 0.037 [−0.034, 0.101] |
| PD-L1, melanoma, co-culture | the model's part closed (relative) | not graded | gate 3 as amended |
| PD-L1, THP-1 | the model's part closed (relative) | the same | IRF1 is read upstream |
| MHC-I, melanoma, co-culture | not graded | the model's part closed (relative) | gate 3 as amended passes; the test with two controls is confirmed |
| MHC-I, melanoma, IFN-γ | not a target | open before the model | B2M and TAPBP admitted; nothing reads them |
| HLA-E, melanoma, IFN-γ | not a target | not a target | on a third no knockout passes (the largest is 1.6) |
| HLA-E, melanoma, co-culture | not a target | not a target | not acted on |

The open items after the step are in [`results/open_items.md`](results/open_items.md), with what became of each of
the 21 the step started from: 6 closed (five amendments in force, final units set aside), 2 at the margin, 3 decided
on development units, 10 still open.

| lever | open items before | after |
|---|---|---|
| M model | 0 | 1 |
| E events | 8 | 6 |
| V views | 3 | 4 |
| S schema | 5 | 4 |
| C claim | 1 | 1 |
| P protocol | 4 | 1 |

## What is yours to decide

- **Knockouts that act on the measured protein or complex** (CD274 for PD-L1; B2M and TAPBP for MHC-I). No reader of
  RNA carries them. Either the claim says so (B2M would be declared one of the readout's own genes, and the rule that
  refuses CD274 would refuse it), or the system is told what was knocked out, which is another system.
- **A tolerance for the two co-culture targets.** No reading there reaches support at 0.10. Without a declared
  tolerance nothing can be said about what RNA carries.
- **HLA-E.** The library has no knockout that moves it enough. Not asked, or a new experiment.
- **States for batch 1.** `batch1_cells.tsv` (written beside the screens' files, not in this repository) lists the
  cells: in co-culture CD44, TAPBP and CDH19 with 1,276 control cells; in the control condition B2M and KLF4 with
  2,645 control cells; in THP-1 BRD4, CMTM6 and POU2F2. A third of them is a calibration sample drawn at random.
  B2M and NPC1 in co-culture, TAPBP and STAT3 in the control condition and CD86 in THP-1 fell to the final units and
  are not to be read before the final audit.

## What the step showed about the design

Running it found five things the design did not have. The first three are now in
[`OUTER_LOOP.md`](../../OUTER_LOOP.md) and in the code; the last two are proposals for the next iteration and were
not applied here.

- A status for a target that is admitted and lost wholly before the model ("open before the model").
- When an item about the rules is closed: when its amendment is in force.
- A loss before the model can be placed for a single event, with its own reference and input view (IRF1).
- *Proposed, gate 3:* read only a shortfall on the fitted cells (the signed *u* of the output against the reference).
  An overshoot is for the audit half to judge.
- *Proposed, the positive control:* a view that reads no more than its permuted copy (*u* at 1) needs no positive
  control to be read as "nothing is read here". As the rule stands, B2M and TAPBP cannot be written as a limit,
  because no reading of their family is supported.

Two smaller corrections made on the way. The round's script first applied the test although the head had nothing to
remove, which gave "refuted" from a bound below zero; the design says "not testable" there, the script was corrected
and run again, and no number changed. And the test of the split could not be run at all: THP-1 has no unit with
states that was not used for the proposal.

## Limits

- Development cells only. The bridge reading and the readers read cells that had been read before; such an item can
  become "decided on development units" and no more.
- Every pooled interval rests on five genes (four without IRF1) and is approximate by the protocol's own number.
- The new target has two events: every reading is per event, by cell, and says nothing about a family.
- For an event the reader of the input view is the best of the readers tried. That favours the input view. It is the
  cautious direction for "dominated" and for "nothing is read", and the incautious one for "upstream".
- THP-1 has no unread control cells: a final audit there will reuse control cells that development has fitted on.
- The IFN-γ condition has no unread unit: nothing final can be said about it.
- The 4,000-gene view was read by ridge only.

## Files

| | |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.md`](results/declared.md) | the declaration, merged before the readings |
| [`results/bridge.json`](results/bridge.json) | the round of 20 November read again under package 2 (`scripts/bridge.py`) |
| [`results/readers.json`](results/readers.json) | every stage and reader for the two co-culture targets and the new one (`scripts/readers.py`) |
| [`results/selection.json`](results/selection.json) | the selection third (`scripts/selection.py`) |
| [`results/round2.json`](results/round2.json) | the round under package 2 on the new target (`scripts/round2.py`) |
| [`results/census.json`](results/census.json), [`results/record_census.json`](results/record_census.json) | the units from the record, the final units, batch 1 (`scripts/record_census.py`, `census.py`, `batch_cells.py`) |
| [`results/open_items.md`](results/open_items.md), [`open_items.json`](results/open_items.json) | targets and open items after the step; what became of the 21 before |
| [`results/outcome.json`](results/outcome.json) | every forecast against what was read |
| [`results/package2.json`](results/package2.json) | package 2 as sealed |
| [`results/tables.md`](results/tables.md) | every table below, generated (`scripts/tables.py`) |

## Reproduce

Environment as [`teg_loop_20261120`](../teg_loop_20261120/). `WORK2` and `OBS2` are new folders; `OBS1` holds the
grids of that round.

    python scripts/declare.py
    python scripts/selection.py --work $ANM_LOOP_WORK --out results
    python scripts/build2.py --work WORK2 --part bridge --coculture-counts COUNTS_COCULTURE.npz
    python scripts/build2.py --work WORK2 --part target --selection results/selection.json --library-counts COUNTS_LIBRARY.npz
    for s in WORK2/*.yaml; do python -m anm_diag scope $s --out OBS2/$(basename $s .yaml); done
    python scripts/bridge.py --observe OBS1 --observe2 OBS2 --work2 WORK2
    python scripts/readers.py --target "PDL1 mel_co" --z OBS1/mel_co_PDL1 --in512 OBS2/mel_co_PDL1_in512 --in4k OBS2/mel_co_PDL1_in4k --perm OBS2/mel_co_PDL1_perm --out results/readers.json
    python scripts/round2.py --observe2 OBS2 --work2 WORK2 --selection results/selection.json
    python scripts/census.py --record results/record_census.json --screen ../loop_reading_20261112/results/data_census.json --out results
    python scripts/items.py && python scripts/tables.py > results/tables.md

The two grids on 4,000 genes for co-culture take about 50 minutes each on one core; everything else takes minutes.

## Tables

### Bridge reading: the repaired output, the round against package 2

| target | reference at Z | output in the round (reader of types) | output with the type from the record | beyond the reference | grade |
|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.032 [−0.044, 0.100] support | 0.043 [−0.021, 0.105] unresolved | 0.037 [−0.034, 0.101] unresolved | 0.005 [−0.021, 0.024] no added loss | relative |
| PD-L1, melanoma, co-culture | 0.117 [0.007, 0.193] unresolved | 0.118 [0.020, 0.187] unresolved | 0.108 [0.009, 0.177] unresolved | −0.008 [−0.027, 0.010] no added loss | not tested (gate 3) |
| PD-L1, THP-1 (five events) | 0.132 [0.016, 0.312] unresolved | 0.163 [0.079, 0.290] unresolved | 0.122 [0.000, 0.262] unresolved | −0.010 [−0.064, 0.024] no added loss | relative |
| MHC-I, melanoma, co-culture | 0.321 [0.265, 0.369] reject | 0.391 [0.342, 0.439] reject | 0.362 [0.315, 0.407] reject | 0.040 [0.025, 0.058] no added loss | relative |

### Gate 3 as amended, and the test with two controls

| target | gate 3: output − reference on fitted cells | passes | the output as trained at the same gate | rows − controls-only (gain) | reading | test |
|---|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.03 [−0.09, 0.09] | yes | 0.87 [0.80, 0.91]: fails | 0.080 [0.019, 0.085] | no step | unresolved |
| PD-L1, melanoma, co-culture | −0.09 [−0.15, 0.13] | **no** | 0.88 [0.80, 0.98]: fails | 0.165 [0.150, 0.182] | a step | not testable |
| PD-L1, THP-1 (five events) | −0.12 [−0.16, 0.09] | yes | 0.83 [0.65, 0.94]: fails | 0.372 [0.338, 0.416] | a step | unresolved |
| MHC-I, melanoma, co-culture | −0.08 [−0.09, −0.06] | yes | 0.65 [0.60, 0.70]: fails | 0.276 [0.250, 0.304] | a step | confirmed |
| PD-L1, THP-1 without IRF1 | −0.06 [−0.12, 0.13] | **no** | 0.89 [0.82, 0.95]: fails | 0.364 [0.322, 0.402] | a step | not testable |

### What each arm removes (gain: fall of \|*u*\| of the output)

| target | as trained | controls only | rows | every cell | rows of the other group | the other group's repair |
|---|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.91 | 0.79 (leaves 0.12) | 0.87 (leaves 0.04) | 0.90 (leaves 0.00) | 0.89 (leaves 0.02) | 0.63 (leaves 0.28) |
| PD-L1, melanoma, co-culture | 1.00 | 0.72 (leaves 0.27) | 0.89 (leaves 0.11) | 0.89 (leaves 0.11) | 0.84 (leaves 0.16) | 0.65 (leaves 0.35) |
| PD-L1, THP-1 (five events) | 0.97 | 0.47 (leaves 0.49) | 0.84 (leaves 0.12) | 0.86 (leaves 0.10) | 0.73 (leaves 0.23) | 0.21 (leaves 0.76) |
| MHC-I, melanoma, co-culture | 0.97 | 0.33 (leaves 0.64) | 0.61 (leaves 0.36) | 0.63 (leaves 0.34) |  |  |

### The THP-1 pool, three ways

| read as | reference | repaired output | beyond the reference |
|---|---|---|---|
| PD-L1, THP-1 (five events) | 0.13 [0.02, 0.31] unresolved | 0.12 [0.00, 0.26] unresolved | −0.01 [−0.06, 0.02] no added loss |
| PD-L1, THP-1 without IRF1 | 0.08 [−0.03, 0.15] unresolved | 0.09 [−0.04, 0.16] unresolved | 0.01 [−0.02, 0.04] no added loss |
| IRF1 with its own row | 0.40 [0.27, 0.50] reject (a reader of Z with its own cells) | 0.68 [0.57, 0.79] reject | 0.28 [0.20, 0.39] loss |
| IRF1 with the pooled rows | the same | 0.54 [0.43, 0.65] reject | 0.15 [0.08, 0.22] unresolved |
| IRF1, every cell of the fitted half | the same | 0.42 [0.30, 0.53] reject | 0.02 [−0.04, 0.09] no added loss |

IRF1 before the head: the RNA view (512 components, ridge, own cells) leaves 0.18 [0.04, 0.30] unresolved; Z leaves 0.40 [0.27, 0.50]; from the RNA view to Z 0.21 [0.12, 0.32]: loss.

### Readers: PD-L1 in co-culture (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| IFNGR1 | 0.10 unr | 0.18 unr | 0.19 rej | 0.09 unr | 0.13 unr | 1.00 | 0.01 [−0.04, 0.05] |
| IFNGR2 | 0.11 unr | 0.17 unr | 0.21 rej | 0.14 unr | 0.16 unr | 1.00 | −0.04 [−0.08, 0.01] |
| JAK1 | −0.06 unr | −0.06 unr | 0.04 unr | −0.07 unr | −0.07 unr | 1.00 | 0.02 [−0.11, 0.15] |
| JAK2 | 0.14 unr | 0.18 rej | 0.23 rej | 0.15 unr | 0.19 rej | 1.00 | −0.01 [−0.05, 0.03] |
| STAT1 | 0.18 rej | 0.21 rej | 0.23 rej | 0.16 unr | 0.19 rej | 1.00 | 0.01 [−0.03, 0.05] |
| pool | 0.10 unr | 0.14 unr | 0.19 rej | 0.11 unr | 0.13 unr | 1.00 | −0.00 [−0.04, 0.03] |

Positive control: of 60 readings (two scopes, five readers, every event and the pool) 0 supported, 40 unresolved, 20 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### Readers: MHC-I in co-culture (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| IFNGR1 | 0.25 rej | 0.29 rej | 0.27 rej | 0.13 unr | 0.27 rej | 1.00 | 0.12 [0.07, 0.16] |
| IFNGR2 | 0.32 rej | 0.37 rej | 0.33 rej | 0.22 rej | 0.34 rej | 1.00 | 0.10 [0.07, 0.14] |
| JAK1 | 0.34 rej | 0.39 rej | 0.32 rej | 0.20 rej | 0.32 rej | 1.00 | 0.13 [0.10, 0.17] |
| JAK2 | 0.29 rej | 0.36 rej | 0.30 rej | 0.18 rej | 0.31 rej | 1.00 | 0.11 [0.07, 0.15] |
| STAT1 | 0.35 rej | 0.39 rej | 0.33 rej | 0.22 rej | 0.35 rej | 1.00 | 0.12 [0.09, 0.16] |
| pool | 0.31 rej | 0.37 rej | 0.31 rej | 0.19 rej | 0.32 rej | 1.00 | 0.12 [0.09, 0.15] |

Positive control: of 60 readings (two scopes, five readers, every event and the pool) 0 supported, 2 unresolved, 58 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### Readers: MHC-I under IFN-γ: B2M and TAPBP (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| B2M | 0.91 rej | 0.97 rej | 1.00 rej | 1.00 rej | 0.99 rej | 1.00 | −0.07 [−0.10, −0.04] dominated |
| TAPBP | 1.00 rej | 1.00 rej | 1.02 rej | 1.03 rej | 1.04 rej | 1.00 | −0.02 [−0.07, 0.02] |
| pool | 0.94 rej | 0.98 rej | 1.00 rej | 1.01 rej | 1.00 rej | 1.00 | −0.06 [−0.09, −0.01] dominated |

Positive control: of 30 readings (two scopes, five readers, every event and the pool) 0 supported, 0 unresolved, 30 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### The selection third

| readout | knockouts that pass A1 on a third | admitted | the next one |
|---|---|---|---|
| PDL1 | CD274, IFNGR1, STAT1, IFNGR2, JAK2, JAK1 | IFNGR1, STAT1, IFNGR2, JAK2, JAK1 | NSG1 (1.5) |
| MHC1 | B2M, TAPBP | B2M, TAPBP | HLA-B (2.3) |
| HLAE | none | none | STAT1 (1.6) |

| events of the round, on a third | |
|---|---|
| PDL1 mel_ifn | IFNGR1 5.6, IFNGR2 4.8, JAK1 3.1, JAK2 4.6, STAT1 5.0 |
| PDL1 mel_co | IFNGR1 5.9, IFNGR2 5.5, JAK1 5.0, JAK2 5.6, STAT1 6.9 |
| MHC1 mel_co | IFNGR1 8.5, IFNGR2 7.6, JAK1 9.4, JAK2 8.2, STAT1 10.9 |
| PDL1 thp1 | IFNGR1 7.6, IFNGR2 8.6, IRF1 2.6, JAK2 7.6, STAT1 4.3 |

### The round under package 2: MHC-I in melanoma under IFN-γ

| event | cells | change (control SD) | head as trained | Z, other cells of the type | Z, own cells | A3 | rows | controls only | every cell | test |
|---|---|---|---|---|---|---|---|---|---|---|
| B2M | 168 | −2.88 | 1.00 [0.99, 1.00] | 1.00 [0.97, 1.02] | 0.91 [0.88, 0.94] | fails (−0.00) | 0.99 [0.97, 1.02] | 1.01 [0.99, 1.03] | 0.92 [0.89, 0.95] | not testable |
| TAPBP | 153 | −1.60 | 1.00 [0.99, 1.02] | 1.01 [0.97, 1.04] | 1.00 [0.96, 1.04] | fails (−0.00) | 1.01 [0.97, 1.06] | 1.01 [0.97, 1.06] | 1.01 [0.97, 1.05] | not testable |

### The census and the final units

| stratum | units (at least 100 cells) | read by the system | unread | final | control cells | with states | unread |
|---|---|---|---|---|---|---|---|
| frangieh2021, Co-culture | 197 | 5 | 192 | 64 | 5897 | 4000 | 1897 |
| frangieh2021, Control | 188 | 0 | 188 | 63 | 3978 | 0 | 3978 |
| frangieh2021, IFNγ | 206 | 206 | 0 | 0 | 6444 | 4000 | 2444 |
| papalexi2021, all | 24 | 5 | 19 | 6 | 2386 | 2386 | 0 |

| stratum | batch 1: chosen for an open item | calibration | screened as clear, but final |
|---|---|---|---|
| frangieh2021, Co-culture | CD44, TAPBP | CDH19 | B2M, NPC1 |
| frangieh2021, Control | B2M | KLF4 | STAT3, TAPBP |
| papalexi2021, all | BRD4, CMTM6 | POU2F2 | CD86 |

### The forecasts

| hypothesis | intervals inside | outcomes: the class given the highest probability happened |
|---|---|---|
| E1 | 0 of 0 | 5 of 5 |
| E2 | 1 of 5 | 6 of 6 |
| E3 | 3 of 3 | 0 of 0 |
| P1 | 2 of 4 | 7 of 8 |
| P2 | 4 of 4 | 6 of 8 |
| P4 | 1 of 1 | 4 of 4 |
| S1 | 2 of 3 | 0 of 2 |
| S2 | 3 of 4 | 4 of 6 |
| V1 | 1 of 2 | 3 of 3 |
| V2 | 2 of 3 | 2 of 3 |
| all | 19 of 29 | 37 of 45 |

| missed | forecast | observed |
|---|---|---|
| P1: gate 3 passes, PDL1 mel_co | yes 0.95, no 0.05 | no |
| P1: gate 3 against the reference, PDL1 thp1 | [-0.11, -0.03] | −0.117 |
| P1: gate 3 against the reference, MHC1 mel_co | [-0.07, -0.01] | −0.079 |
| P2: outcome of the test with two controls, PDL1 mel_co | unresolved 0.65, confirmed 0.35 | not testable |
| P2: outcome of the test with two controls, MHC1 mel_co | unresolved 0.7, confirmed 0.3 | confirmed |
| S1: PDL1 mel_ifn verdict | support 0.6, unresolved 0.4 | unresolved |
| S1: PDL1 thp1 (five events pooled), repaired output | [0.14, 0.17] | 0.122 |
| S1: PDL1 mel_ifn closed (absolute) on development cells | yes 0.55, no 0.45 | no |
| S2: type B (IRF1, its own row) repaired output | [0.35, 0.55] | 0.678 |
| S2: type B adds no loss beyond a reader of Z with its own cells | yes 0.6, no 0.4 | no |
| S2: type B reading upstream | yes 0.45, no 0.55 | yes |
| V1: PDL1 mel_co, 512 components, boosting, own cells | [0.15, 0.3] | 0.107 |
| V2: MHC1 mel_co, 512 components, boosting, own cells | [0.22, 0.42] | 0.194 |
| V2: boosting reads at least 0.10 better than ridge at a stage | yes 0.2, no 0.8 | yes |
| E2: head as trained, TAPBP | [0.85, 1.0] | 1.002 |
| E2: Z, own cells, B2M | [0.3, 0.9] | 0.913 |
| E2: Z, own cells, TAPBP | [0.4, 0.95] | 0.997 |
| E2: input view, own cells, B2M (best reader) | [0.2, 0.85] | 0.986 |
