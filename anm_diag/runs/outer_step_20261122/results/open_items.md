## Targets after the step

| readout | group | status | own output | added beyond the reference |
|---|---|---|---|---|
| PDL1 | mel_ifn | the model's part closed (relative) | 0.04 [−0.03, 0.10] unresolved | 0.01 [−0.02, 0.02] no added loss |
| PDL1 | mel_co | not graded | 0.11 [0.01, 0.18] unresolved | −0.01 [−0.03, 0.01] no added loss |
| PDL1 | thp1 | the model's part closed (relative) | 0.12 [0.00, 0.26] unresolved | −0.01 [−0.06, 0.02] no added loss |
| MHC1 | mel_co | the model's part closed (relative) | 0.36 [0.32, 0.41] reject | 0.04 [0.03, 0.06] no added loss |
| MHC1 | mel_ifn | open before the model | 0.99 [0.97, 1.02] reject | −0.00 [−0.02, 0.01] no added loss |
| HLAE | mel_ifn | not a target |  |  |
| HLAE | mel_co | not a target |  |  |

## Open items after the step

| id | readout | group | event | reading | lever | question | value | the next step needs |
|---|---|---|---|---|---|---|---|---|
| O01 | PDL1 | mel_ifn | CD274 | self event | C | the knockout of the readout's own gene | 0.91 [0.87, 0.95] | the owner's decision |
| O02 | PDL1 | mel_ifn | pool | unresolved | E | is the output of the model inside the margin | 0.04 [−0.03, 0.10] | at the margin: the pool has no further unit that passes admission |
| O03 | PDL1 | mel_ifn | pool | pooled, and differs | S | the whole repair across: melanoma IFN-gamma and melanoma co-culture | 0.24 [0.16, 0.28] (mel_ifn); 0.23 [0.19, 0.28] (mel_co) | units of the pool |
| O04 | PDL1 | mel_ifn | pool | pooled, may be the same | S | the rows across: melanoma IFN-gamma and melanoma co-culture | −0.01 [−0.03, 0.02] (mel_ifn); 0.04 [0.02, 0.06] (mel_co) | units of the pool |
| O05 | PDL1 | mel_co | pool | invalid | V | the reading of the input view fails the check: positive control |  | no new data first; then the owner's decision |
| O06 | PDL1 | mel_co | pool | rule | P | gate 3 as amended (3a) does not pass: on the cells the repair was fitted on, the interval of |u of the output| - |u of the reference| reaches above the step margin | −0.09 | no new data |
| O07 | PDL1 | mel_co | pool | unresolved | E | is anything lost before the model | 0.12 [0.01, 0.19] | at the margin: more than 4 times the units would be needed |
| O08 | PDL1 | thp1 | pool | unresolved | E | is anything lost before the model | 0.13 [0.02, 0.31] | at the margin: more than 4 times the units would be needed |
| O09 | PDL1 | thp1 | IRF1 | upstream | M | the input view reads it and the reference stage does not | 0.21 [0.12, 0.32] | a change of the encoder |
| O10 | PDL1 | thp1 | IRF1 | pooled, and differs | S | one event inside the pool: IRF1 | 0.54 [0.43, 0.65] | units of the pool |
| O11 | PDL1 | thp1 | pool | pooled, and differs | S | the whole repair across: THP-1 stimulated and melanoma IFN-gamma | 0.62 [0.48, 0.71] | units of the pool |
| O12 | PDL1 | thp1 | pool | unresolved | E | can the rows be merged across: THP-1 stimulated and melanoma IFN-gamma | 0.10 [0.01, 0.15] | at the margin: more than 4 times the units would be needed |
| O13 | MHC1 | mel_co | pool | invalid | V | the reading of the input view fails the check: positive control |  | no new data first; then the owner's decision |
| O14 | MHC1 | mel_ifn | B2M | invalid | V | the reading of the input view fails the check: dominated view | −0.07 [−0.10, −0.04] | no new data |
| O15 | MHC1 | mel_ifn | pool | invalid | V | the reading of the input view fails the check: positive control |  | no new data first; then the owner's decision |
| O16 | HLAE | mel_ifn | pool | small response | E | the change of the readout against the half-width of its interval | 2.21 | the pool at hand has no event that moves this readout: not asked, or on the experiment list |
| O17 | HLAE | mel_co | pool | small response | E | the change of the readout against the half-width of its interval | 3.12 | units of the pool, or the owner's decision |

| lever | open items |
|---|---|
| M | 1 |
| E | 6 |
| V | 4 |
| S | 4 |
| C | 1 |
| P | 1 |

## What became of the 21 items the step started from

| id | reading | lever | state after the step | evidence |
|---|---|---|---|---|
| O01 | self event | C | open | not acted on: the owner's decision (a limit for knockouts of the readout's own gene, or a system that is told the perturbation) |
| O02 | unresolved | E | at the margin | bridge reading with the type from the record: 0.037 [-0.034, 0.101], unresolved; on a selection third no further knockout of the library passes for PD-L1, so the pool cannot decide it. Read on development cells |
| O03 | pooled, and differs | S | open | not acted on: needs units the system has not been read on |
| O04 | pooled, may be the same | S | open | not acted on: needs units in which the condition varies; for PD-L1 the only further event of co-culture (NPC1) fell to the final units |
| O05 | invalid | V | decided on development units | richer readers: the best reader of the input view leaves 0.11; from the input view to Z -0.00 [-0.04, 0.03]; dominated: False |
| O06 | invalid | V | open | richer readers: of 60 readings 0 are supported. No positive control: the owner's decision (a tolerance for this target) |
| O07 | unresolved | E | at the margin | unchanged: Z reads 0.12 [0.01, 0.19]; no further event of the type is left for development (NPC1 is final) |
| O08 | unresolved | E | open | read again without IRF1: the reference of the other four is 0.08 [-0.03, 0.15], unresolved; new THP-1 events are in batch 1 |
| O09 | pooled, and differs | S | decided on development units | the split as proposed is not supported: with its own row IRF1 leaves 0.68, with the pooled rows 0.54. Its loss is before the head: the RNA view reads 0.18, Z 0.40; from the RNA view to Z 0.21 [0.12, 0.32]. Read `upstream` (a new item) |
| O10 | pooled, and differs | S | open | not acted on: needs units the system has not been read on |
| O11 | unresolved | E | open | not acted on: needs new THP-1 events; BRD4 and CMTM6 are in batch 1 |
| O12 | small response | E | decided on development units | B2M and TAPBP are admitted on a selection third; the target was read (round2.json). What it showed is in the new items |
| O13 | invalid | V | open | richer readers: of 60 readings 0 are supported. No positive control: the owner's decision (a tolerance for this target) |
| O14 | rule | P | closed | amendment 3a of protocol version 3 is in force (LOOP_PROTOCOL.md) |
| O15 | rule | P | closed | amendment 3b of protocol version 3 is in force (LOOP_PROTOCOL.md) |
| O16 | small response | E | open | on a selection third no knockout of the library passes for HLA-E (the largest is 1.6): the owner's decision (not asked, or the experiment list) |
| O17 | small response | E | open | not acted on: units of the pool without states, or the owner's decision |
| O18 | mistyped | S | closed | amendment 3e of protocol version 3 is in force (LOOP_PROTOCOL.md) |
| O19 | selection | P | closed | amendment 3d of protocol version 3 is in force (LOOP_PROTOCOL.md) |
| O20 | rule | P | closed | amendment 3c of protocol version 3 is in force (LOOP_PROTOCOL.md) |
| O21 | no final data | E | closed | 133 units are set aside as final (census.json); none has states yet |
