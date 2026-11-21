## Targets at the end of the round

| readout | group | part | status | own output after the repair | added beyond Z |
|---|---|---|---|---|---|
| PDL1 | mel_ifn | development | the head's part closed (relative) | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] no added loss |
| PDL1 | mel_co | development | the head's part closed (relative) | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] no added loss |
| PDL1 | thp1 | development | the head's part closed (relative) | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] no added loss |
| MHC1 | mel_ifn | replication | not a target |  |  |
| MHC1 | mel_co | replication | not graded | 0.39 [0.34, 0.44] reject | 0.07 [0.05, 0.09] no added loss |
| HLAE | mel_ifn | replication | not a target |  |  |
| HLAE | mel_co | replication | not a target |  |  |

## Open items

| id | readout | group | event | reading | lever | quantity | value | units now → needed | the next step needs |
|---|---|---|---|---|---|---|---|---|---|
| O01 | PDL1 | mel_ifn | pool | unresolved | E | u of the repaired output, against the margin | 0.04 [−0.02, 0.11] | 5 → 6 | at the margin: the pool has no further unit of the kind |
| O02 | PDL1 | mel_ifn | pool | pooled, and differs | S | the whole repair across: melanoma IFN-gamma and melanoma co-culture | 0.24 [0.16, 0.28] (mel_ifn); 0.23 [0.19, 0.28] (mel_co) |  | units of the pool |
| O03 | PDL1 | mel_ifn | pool | pooled, may be the same | S | the rows across: melanoma IFN-gamma and melanoma co-culture | −0.01 [−0.03, 0.02] (mel_ifn); 0.05 [0.03, 0.07] (mel_co) |  | units of the pool |
| O04 | PDL1 | mel_co | pool | unresolved | E | u at Z, read from the other genes of the group: is anything lost before the head | 0.12 [0.01, 0.19] | 5 → 213 | at the margin: more than four times the units would be needed |
| O05 | PDL1 | mel_co | pool | invalid | V | from the input view to Z, own scope: a later stage reads better than an earlier one | −0.09 [−0.13, −0.04] |  | no new data |
| O06 | PDL1 | thp1 | pool | unresolved | E | u at Z, read from the other genes of the group: is anything lost before the head | 0.13 [0.02, 0.31] | 5 → 67 | at the margin: more than four times the units would be needed |
| O07 | PDL1 | thp1 | IRF1 | pooled, and differs | S | one event inside the pool: IRF1, THP-1 stimulated | 0.54 [0.43, 0.65] |  | units of the pool |
| O08 | PDL1 | thp1 | pool | pooled, and differs | S | the whole repair across: THP-1 stimulated and melanoma IFN-gamma | 0.62 [0.48, 0.71] |  | units of the pool |
| O09 | PDL1 | thp1 | pool | unresolved | E | the rows across: THP-1 stimulated and melanoma IFN-gamma | 0.14 [0.09, 0.18] | 5 → 9 | within reach if the rest of the screen holds units of the type: a census first |
| O10 | PDL1 | mel_ifn | CD274 | self event | C | u at Z for the knockout of the readout's own gene | 0.91 [0.87, 0.95] |  | the owner's decision |
| O11 | all | mel_ifn | pool | mistyped | S | how a cell gets its type | 0.96 |  | no new data |
| O12 | MHC1 | mel_ifn | pool | small response | E | the change of the readout against the half-width of its interval, pooled | 3.09 | 5 → see note | units of the pool, or the owner's decision |
| O13 | MHC1 | mel_co | pool | rule | P | gate 3: the change reaches its objective on the cells it was fitted on | 0.28 |  | no new data |
| O14 | MHC1 | mel_co | pool | rule | P | the test of the hypothesis |  |  | no new data |
| O15 | MHC1 | mel_co | pool | unreadable | V | u at Z and at the input view, with the event's own cells | 0.31 [0.26, 0.36] |  | no new data first; then the owner's decision |
| O16 | HLAE | mel_ifn | pool | small response | E | the change of the readout against the half-width of its interval, pooled | 2.21 | 5 → 10 | units of the pool, or the owner's decision |
| O17 | HLAE | mel_co | pool | small response | E | the change of the readout against the half-width of its interval, pooled | 3.12 | 5 → see note | units of the pool, or the owner's decision |
| O18 | all |  | pool | selection | P | A1 (the response against its interval) was read on the cells that were then audited |  |  | no new data |
| O19 | all |  | pool | rule | P | a group ready as a pool with no admitted event |  |  | no new data |
| O20 | all |  | pool | no final data | E | units of the knockout claim that no run has read | 0 |  | units of the pool |

## Count by lever

| lever | what it is | open items |
|---|---|---|
| M | the model (inner loop) | 0 |
| E | events: which units of the pool are read next, and in which role | 8 |
| V | views: the readers of a stage, the input, the readouts | 2 |
| S | schema: the fields of an event and what is pooled | 5 |
| C | claim: the targets, the grade and the tolerance asked for | 1 |
| P | protocol: the gates and the rules | 4 |
