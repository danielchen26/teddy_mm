## Targets at the end of the round

| readout | group | part | status | own output after the repair | added beyond Z |
|---|---|---|---|---|---|
| PDL1 | mel_ifn | development | the model's part closed (relative) | 0.04 [−0.02, 0.11] unresolved | 0.01 [−0.03, 0.03] no added loss |
| PDL1 | mel_co | development | the model's part closed (relative) | 0.12 [0.02, 0.19] unresolved | 0.00 [−0.02, 0.02] no added loss |
| PDL1 | thp1 | development | the model's part closed (relative) | 0.16 [0.08, 0.29] unresolved | 0.03 [−0.03, 0.08] no added loss |
| MHC1 | mel_ifn | replication | not a target |  |  |
| MHC1 | mel_co | replication | not graded | 0.39 [0.34, 0.44] reject | 0.07 [0.05, 0.09] no added loss |
| HLAE | mel_ifn | replication | not a target |  |  |
| HLAE | mel_co | replication | not a target |  |  |

## Open items

| id | readout | group | event | reading | lever | question | value | units now → needed | after | the next step needs |
|---|---|---|---|---|---|---|---|---|---|---|
| O01 | PDL1 | mel_ifn | CD274 | self event | C | the knockout of the readout's own gene | 0.91 [0.87, 0.95] |  |  | the owner's decision |
| O02 | PDL1 | mel_ifn | pool | unresolved | E | is the output of the model inside the margin | 0.04 [−0.02, 0.11] | 5 → 6 | O18 | first the bridge reading with the type from the record, which changes the cells this estimate is read on; then: at the margin: the pool has no further unit that passes admission |
| O03 | PDL1 | mel_ifn | pool | pooled, and differs | S | the whole repair across: melanoma IFN-gamma and melanoma co-culture | 0.24 [0.16, 0.28] (mel_ifn); 0.23 [0.19, 0.28] (mel_co) |  |  | units of the pool |
| O04 | PDL1 | mel_ifn | pool | pooled, may be the same | S | the rows across: melanoma IFN-gamma and melanoma co-culture | −0.01 [−0.03, 0.02] (mel_ifn); 0.05 [0.03, 0.07] (mel_co) |  |  | units of the pool |
| O05 | PDL1 | mel_co | pool | invalid | V | the reading of the input view fails the check: dominated view | −0.09 [−0.13, −0.04] |  |  | no new data |
| O06 | PDL1 | mel_co | pool | invalid | V | the reading of the input view fails the check: positive control |  |  |  | no new data first; then the owner's decision |
| O07 | PDL1 | mel_co | pool | unresolved | E | is anything lost before the model | 0.12 [0.01, 0.19] | 5 → 213 |  | at the margin: more than 4 times the units would be needed |
| O08 | PDL1 | thp1 | pool | unresolved | E | is anything lost before the model | 0.13 [0.02, 0.31] | 5 → 67 | O09 | at the margin: more than 4 times the units would be needed |
| O09 | PDL1 | thp1 | IRF1 | pooled, and differs | S | one event inside the pool: IRF1 | 0.54 [0.43, 0.65] |  |  | units of the pool |
| O10 | PDL1 | thp1 | pool | pooled, and differs | S | the whole repair across: THP-1 stimulated and melanoma IFN-gamma | 0.62 [0.48, 0.71] |  |  | units of the pool |
| O11 | PDL1 | thp1 | pool | unresolved | E | can the rows be merged across: THP-1 stimulated and melanoma IFN-gamma | 0.14 [0.09, 0.18] | 5 → 9 | O09 | within reach if the pool holds units of the type: a census first |
| O12 | MHC1 | mel_ifn | pool | small response | E | the change of the readout against the half-width of its interval | 3.09 | 5 → see note |  | within reach now: events with states at hand pass the screen: B2M, TAPBP, SSR2, STAT3, CD44 |
| O13 | MHC1 | mel_co | pool | invalid | V | the reading of the input view fails the check: positive control |  |  |  | no new data first; then the owner's decision |
| O14 | MHC1 | mel_co | pool | rule | P | gate 3 is absolute: the change reaches its objective on the cells it was fitted on | 0.28 |  |  | no new data |
| O15 | MHC1 | mel_co | pool | rule | P | the hypothesis had no control arm |  |  |  | no new data |
| O16 | HLAE | mel_ifn | pool | small response | E | the change of the readout against the half-width of its interval | 2.21 | 5 → 10 |  | the pool at hand has no event that moves this readout: not asked, or on the experiment list |
| O17 | HLAE | mel_co | pool | small response | E | the change of the readout against the half-width of its interval | 3.12 | 5 → see note |  | units of the pool, or the owner's decision |
| O18 | all | mel_ifn | pool | mistyped | S | how a cell gets its type | 0.96 |  |  | no new data |
| O19 | all |  | pool | selection | P | events were screened on the cells that were then audited |  |  |  | no new data |
| O20 | all |  | pool | rule | P | a group ready as a pool with no admitted event |  |  |  | no new data |
| O21 | all |  | pool | no final data | E | units of the claim that no run has read | 0 |  |  | units of the pool |

## Counts

| lever | what it is | open items |
|---|---|---|
| M | the model (inner loop) | 0 |
| E | events: which units of the pool are read next, and in which role | 8 |
| V | views: the readers of a stage, the input, the readouts | 3 |
| S | schema: the fields of an event and what is pooled | 5 |
| C | claim: the targets, the grade and the tolerance asked for | 1 |
| P | protocol: the gates and the rules | 4 |

| what the next step needs | open items |
|---|---|
| at the margin | 2 |
| no new data | 7 |
| no new data first; then the owner's decision | 2 |
| the owner's decision | 2 |
| units at hand | 1 |
| units of the pool | 6 |
| units of the pool, or the owner's decision | 1 |
