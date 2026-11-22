| change | lever | items | hypothesis |
|---|---|---|---|
| 3a gate 3 is stated against the reference stage | P | O14 | P1 |
| 3b the test of a repair by rows has two controls | P | O15 | P2 |
| 3c a group ready as a pool with no admitted event is not a target | P | O20 | P3 |
| 3d a selection split for every new admission | P | O19 | P4 |
| 3e the type of a cell comes from the record | S | O18 O02 | S1 |
| the THP-1 pool is read as two types, provisionally: four genes that signal to STAT1, and IRF1 | S | O09 O08 O11 | S2 |
| reader specification: the input view is read by ridge and boosting on 512 components and by ridge on 4,000 genes; Z by ridge and boosting; a permuted view as control | V | O05 O06 | V1 |
| the same reader specification for MHC-I in co-culture, with the positive control evaluated | V | O13 | V2 |
| the events of MHC-I in melanoma under IFN-gamma are selected from the library on a selection third | E | O12 | E1 |
| a new target: MHC-I in melanoma under IFN-gamma with the admitted events | E | O12 | E2 |
| a census from the record, and a third of the unread units of each stratum set aside as final | E | O21 | E3 |

| items not acted on | why |
|---|---|
| O01 | the owner's decision: a limit for knockouts of the readout's own gene, or a system that is told the perturbation |
| O03 O04 O10 | need units the system has not been read on; no such unit has states |
| O07 | at the margin as it stands; read again after V1 |
| O16 O17 | the owner's decision, or units of the pool without states |

### P1 (rule; items O14)

- **quantity**: |u of the repaired output| - |u of the reference stage| on the cells the repair was fitted on, and its interval
- **forecast**: PDL1 mel_ifn: [-0.02, 0.04]; PDL1 mel_co: [-0.13, -0.07]; PDL1 thp1: [-0.11, -0.03]; MHC1 mel_co: [-0.07, -0.01]; passes: {each of the four: 0.95}
- **rule**: gate 3 passes when the interval lies below the step margin (0.10)
- **rival**: the repair does not reach the reference stage even on its own cells: a failed fit, and the target stays 'not graded'
- **control**: the output as trained fails the same gate on the same cells (forecast: on all four)
- **units**: proposal IFNGR1@mel_co, IFNGR2@mel_co, JAK1@mel_co, JAK2@mel_co, STAT1@mel_co; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### P2 (rule; items O15)

- **quantity**: outcome of the test of the repair by rows with two controls (the repair of another group where there is one; the readout on control cells without rows)
- **forecast**: PDL1 mel_ifn: {unresolved: 0.85, refuted: 0.1, confirmed: 0.05}; PDL1 mel_co: {unresolved: 0.65, confirmed: 0.35}; PDL1 thp1: {unresolved: 0.85, confirmed: 0.1, refuted: 0.05}; MHC1 mel_co: {unresolved: 0.7, confirmed: 0.3}; gain of rows minus gain of controls-only: {PDL1 mel_ifn: [0.05, 0.1], PDL1 mel_co: [0.13, 0.19], PDL1 thp1: [0.3, 0.42], MHC1 mel_co: [0.24, 0.31]}; that difference lies above the step margin: {PDL1 mel_ifn: 0.05, PDL1 mel_co: 0.85, PDL1 thp1: 0.97, MHC1 mel_co: 0.99}
- **rule**: loop.test: confirmed if the rows arm's gain lies above the grid's ceiling minus the step margin and every control arm's below it
- **rival**: a readout fitted on the group's control cells already removes the loss: what was missing is a readout of cells of this type, and the rows add little
- **control**: this amendment is the control; the difference rows minus controls-only is read on paired resamples
- **units**: proposal IFNGR1@mel_co, IFNGR2@mel_co, JAK1@mel_co, JAK2@mel_co, STAT1@mel_co; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### P3 (rule; items O20)

- **quantity**: status of the two groups that pass admission as a pool with no admitted event
- **forecast**: MHC1 mel_ifn (five pathway genes): not a target; HLAE mel_co: not a target
- **rule**: a repair by rows needs an admitted row; without one the group is not a target and the reason is an open item (small response)
- **rival**: more cells per event would let single events pass: then it is a question of resolution (E), not of the rule
- **control**: the pooled response passes A1 while no single event does (the round, replication run)
- **units**: proposal IFNGR1@mel_co, IFNGR1@mel_ifn, IFNGR2@mel_co, IFNGR2@mel_ifn, JAK1@mel_co, JAK1@mel_ifn, JAK2@mel_co, JAK2@mel_ifn, STAT1@mel_co, STAT1@mel_ifn; test none. a definition: no reading changes

### P4 (rule; items O19)

- **quantity**: the events admitted in the round, screened again on a random third of their cells: |change| / half-width
- **forecast**: every event of PDL1 mel_ifn at least 3: 0.8; every event of PDL1 mel_co at least 3: 0.97; IRF1@thp1 at least 3: 0.4; the other four of thp1 at least 3: 0.97; library knockouts other than the readout's own gene and the five that pass by chance (PDL1, of about 200): [0, 1]
- **rule**: from this step on an event is admitted on its selection third only; the events of the round keep their admission and this screen is reported
- **rival**: the selection made on all cells favoured no event: every event passes on a third as well
- **control**: the same screen on the library's other knockouts: how many pass on a third
- **units**: proposal IFNGR1@mel_co, IFNGR1@mel_ifn, IFNGR1@thp1, IFNGR2@mel_co, IFNGR2@mel_ifn, IFNGR2@thp1, IRF1@thp1, JAK1@mel_co, JAK1@mel_ifn, JAK2@mel_co, JAK2@mel_ifn, JAK2@thp1, STAT1@mel_co, STAT1@mel_ifn, STAT1@thp1; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### S1 (rule; items O18 O02)

- **quantity**: u of the repaired output, each cell read by the readout of its recorded type
- **forecast**: PDL1 mel_ifn: [-0.01, 0.05]; PDL1 mel_ifn verdict: {support: 0.6, unresolved: 0.4}; PDL1 mel_co: [0.1, 0.13]; PDL1 thp1 (five events pooled): [0.14, 0.17]; PDL1 mel_ifn status closed (absolute) on development cells: 0.55
- **rule**: the type is taken from the record whatever this reading shows: the reader of types gave 96% of the knockout cells of one group their own type, below the chosen 99%
- **rival**: the 0.04 that was left is not from mistyped cells: u does not move
- **control**: the reader of types as a component: the share of held-out cells it gives their own type, per group
- **units**: proposal IFNGR1@mel_ifn, IFNGR2@mel_ifn, JAK1@mel_ifn, JAK2@mel_ifn, STAT1@mel_ifn; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### S2 (split; items O09 O08 O11)

- **quantity**: the THP-1 pool read as two types: A = IFNGR1, IFNGR2, JAK2, STAT1; B = IRF1
- **forecast**: A reference (other genes of the type, at Z): [0.02, 0.12]; A reference verdict: {unresolved: 0.65, support: 0.2, reject: 0.15}; A repaired output: [0.06, 0.15]; A repaired output verdict: {unresolved: 0.7, support: 0.1, reject: 0.2}; A adds no loss beyond its reference: 0.8; B repaired output (its own row): [0.35, 0.55]; B verdict reject: 0.9; B adds no loss beyond a reader of Z with its own cells: 0.6; B from the RNA view to Z, own cells: [0.1, 0.35]; B reading upstream (interval above the step margin): 0.45
- **rule**: operated provisionally. Adopted only by its test on units not used for the proposal: within each type the repair with the type's rows adds no loss, and the repair with the pooled rows adds a resolved loss on the same units. THP-1 has no such unit with states yet, so the test is not run in this step
- **rival**: the missing field is not the one proposed (the place of the gene in the pathway: signalling to STAT1, or the transcription factor it induces). Rivals from the record: the replicate (IRF1 cells are 56% from replicate 1, the others 50%, controls 40%) and the size of the response, which is the endpoint and may not type an event
- **control**: a split of the same sizes at random does not pass: no other single event of the pool lies outside the pool's interval (IFNGR1 0.10, IFNGR2 0.14, JAK2 0.17, STAT1 0.02; pool 0.16 [0.08, 0.29])
- **units**: proposal IFNGR1@thp1, IFNGR2@thp1, IRF1@thp1, JAK2@thp1, STAT1@thp1; test none. operated in the bridge reading beside the coarser schema; read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### V1 (view; items O05 O06)

- **quantity**: PD-L1 in co-culture: u with the event's own cells at the input view, per reader; the gap from the best reader of the input to Z
- **forecast**: RNA, 512 components, boosting: [0.15, 0.3]; RNA, 4,000 genes, ridge: [0.08, 0.2]; the input view is no longer dominated (the gap's interval reaches zero): 0.6; some event or the pool is supported by some reader at some stage: 0.15
- **rule**: the input view can be read when its gap to Z no longer lies below zero and the positive control holds; until both, no reading of where the loss is
- **rival**: Z carries what no reader of RNA of these classes finds: the view stays dominated, and Z bounds it
- **control**: the same view with its rows permuted across cells reads nothing (u at least 0.9)
- **units**: proposal IFNGR1@mel_co, IFNGR2@mel_co, JAK1@mel_co, JAK2@mel_co, STAT1@mel_co; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### V2 (view; items O13)

- **quantity**: MHC-I in co-culture: u with the event's own cells, per stage and reader (pool of five)
- **forecast**: Z, boosting: [0.2, 0.4]; RNA, 512 components, boosting: [0.22, 0.42]; RNA, 4,000 genes, ridge: [0.22, 0.36]; every reader at every stage rejects: 0.8; some event or the pool is supported by some reader at some stage: 0.1; boosting reads at least 0.10 better than ridge at the same stage: 0.2
- **rule**: `unreadable` is drawn only if the positive control holds, two reader classes of which one is not linear reject at the input view with the event's own cells, and the view is not dominated. Without a positive control the reading stays invalid and the next action is the owner's: a tolerance for this target
- **rival**: the relation between the readout and the state is not linear across the knockout (within the control cells a small change, across the knockout a large one): then boosting reads it
- **control**: the same view with its rows permuted across cells reads nothing (u at least 0.9)
- **units**: proposal IFNGR1@mel_co, IFNGR2@mel_co, JAK1@mel_co, JAK2@mel_co, STAT1@mel_co; test none. read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'

### E1 (resolution; items O12)

- **quantity**: which knockouts of the IFN-gamma library are admitted for MHC-I on a selection third (A1: resolved and at least 3 half-widths; A2: not HLA-A, HLA-B, HLA-C)
- **forecast**: B2M and TAPBP both admitted: 0.95; exactly those two: 0.7; SSR2, STAT3, CD44: each: 0.15; any other: 0.05
- **rule**: an event is admitted on its selection third; those cells, and a third of the control cells, enter no reading of u
- **rival**: the screen on all cells picked its near-threshold events by chance: they fail on a third
- **control**: the count of library knockouts that pass on a third against the count that passed on all cells (7, of which 5 are not the readout's own genes)
- **units**: proposal the screens of 13 and 22 November read the protein response of every library knockout on all its cells; test B2M@mel_ifn, TAPBP@mel_ifn. the system was never read on these units for this readout; their protein response was

### E2 (rows; items O12)

- **quantity**: MHC-I in melanoma under IFN-gamma, the admitted events (expected B2M and TAPBP), on the two thirds of their cells not used for admission
- **forecast**: head as trained, B2M: [0.85, 1.0]; head as trained, TAPBP: [0.85, 1.0]; Z, own cells, B2M: [0.3, 0.9]; Z, own cells, B2M verdict reject: 0.85; Z, own cells, TAPBP: [0.4, 0.95]; Z, own cells, TAPBP verdict reject: 0.85; input view, own cells, B2M (best reader): [0.2, 0.85]; reading upstream, B2M: 0.3; reading upstream, TAPBP: 0.25; the group passes A3 (the head loses what other cells of this kind read): 0.45; if it does: the repair by rows adds no loss beyond the reference: 0.7; if it does: test outcome: {unresolved: 0.5, confirmed: 0.3, refuted: 0.2}; some event supported by some reader at some stage (positive control): 0.15
- **rule**: loop.test with the three arms: rows of the admitted events; rows of the five pathway genes of the same type, which move MHC-I little (cells of a kind already covered); the repair of co-culture applied here (another map). Gate 3 as amended. Read per event: two events give no statement about a family
- **rival**: these knockouts act on the measured complex itself (B2M is its light chain, TAPBP loads it), as CD274 does for PD-L1: no reader of RNA carries them, and the item goes to the claim (a limit, or a system that is told the perturbation), not to rows. The field would be 'relation to the readout'
- **control**: the three arms of the protocol for a coverage claim
- **units**: proposal CD274@mel_ifn, IFNGR1@mel_co, IFNGR1@mel_ifn, IFNGR2@mel_co, IFNGR2@mel_ifn, JAK1@mel_co, JAK1@mel_ifn, JAK2@mel_co, JAK2@mel_ifn, STAT1@mel_co, STAT1@mel_ifn; test B2M@mel_ifn, TAPBP@mel_ifn. forecasts read from the five pathway genes (MHC-I, head as trained 0.88) and from the CD274 probe (0.91 at Z with its own cells)

### E3 (rule; items O21)

- **quantity**: units set aside as final: a third of the units of each stratum that the system has not been read on
- **forecast**: share per stratum: [0.3, 0.37]
- **rule**: drawn once by unit with a fixed seed, stratified by study and condition, before any state is computed for them; a run that touches them is refused unless it is the final audit
- **rival**: a later split along a field that was not stratified leaves one side without final units: then that side has no certified statement
- **control**: the draw is by unit name only: no response and no state enters it
- **units**: proposal none; test none. no reading is made
