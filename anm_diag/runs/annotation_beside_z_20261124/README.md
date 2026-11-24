# The head reads the first annotation position beside Z

24 November 2026. The run before found that TEDDY's annotation positions, a summary the model learned in
pre-training, read IRF1's response better than the mean over the genes, by a margin its rule did not resolve
([`annotation_input_20261124`](../annotation_input_20261124/)). This run tests that lead in the round: the head is
given the first annotation position beside the state Z, for the five targets of the run of 23 November
([`input_beside_z_20261123`](../input_beside_z_20261123/)), with the same cells, halves and resamples. The stages,
the rule, the checks and 19 forecasts were committed before the states of the melanoma cells were computed
([`DESIGN.md`](DESIGN.md)).

Everything here is exploratory. The cells are development cells that earlier runs have read. No final unit was read
and no model was trained.

## What came out

**For IRF1 the repair gains a step when the head reads the annotation position beside Z, and it does not close. No
target loses by it, and every grade is kept. It gives less than the input view gave, and unlike the input view it
costs nothing anywhere.** H1 is unresolved, H2 is refuted, the three controls hold.

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. `ZC` is Z beside
the first annotation position; `ZV` is Z beside the input view (the run of 23 November, made again here).

| target | events | the repaired output, head reads Z | head reads ZC | gain | grade at Z → at ZC | for comparison: head reads ZV |
|---|---|---|---|---|---|---|
| PD-L1, THP-1 | IRF1, on its own | 0.68 [0.57, 0.79] | 0.51 [0.39, 0.61] | **+0.17 [0.12, 0.23]: a step** | not tested (gate 3) → the same | 0.39; gain +0.29, a step |
| MHC-I, melanoma, co-culture | five pathway genes, a pool | 0.36 [0.32, 0.41] | 0.32 [0.28, 0.37] | +0.04 [0.03, 0.05]: no step | relative → relative | 0.26; gain +0.10, unresolved |
| control: PD-L1, THP-1 | four pathway genes | 0.09 [−0.04, 0.16] | 0.07 [−0.00, 0.13] | +0.02 [−0.02, 0.05]: no step | not tested (gate 3) → the same | −0.02 |
| control: PD-L1, co-culture | five pathway genes | 0.11 [0.01, 0.18] | 0.10 [−0.00, 0.16] | +0.01 [−0.01, 0.03]: no step | not tested (gate 3) → the same | 0.07 |
| control: PD-L1, IFN-γ | five pathway genes | 0.04 [−0.03, 0.10] | 0.00 [−0.07, 0.07] | +0.04 [−0.04, 0.05]: no step | **relative → absolute** | −0.05; it lost its grade |

**IRF1 (H1: unresolved).**

- A reader of the joint stage with IRF1's own cells leaves 0.29 [0.16, 0.40], against 0.40 at Z: the 0.29 of the run
  before, now inside the round. From the input view (0.18) that is +0.10 [0.01, 0.20], unresolved: no longer a
  resolved loss, and not closed.
- The round's repair gains a step: 0.68 to 0.51. It is still rejected, it adds 0.22 beyond the reader of the joint
  stage (a loss), and gate 3 does not pass (on the cells it was fitted on it leaves 0.44). No grade is reached.
- Fitted on every cell of the fitted half, the same readout leaves 0.28 [0.15, 0.39], no worse than the reference.
  As on 23 November, what is short at this stage is the record: one record of about 540 cells moves 1,024 weights
  part of the way.
- **Without IRF1's record the annotation position does not hurt, and may help a little.** A readout of control cells
  with no record leaves 0.82 [0.70, 0.96] at ZC against 0.90 at Z (gain +0.09 [0.04, 0.13], unresolved). With the
  input view it was worse than at Z (1.00). The difference between the two wider stages in this arm is 0.19
  [0.11, 0.28] in favour of the annotation position.
- **The input view still gives more where there is a record.** The repaired output is 0.39 with the input view
  beside Z and 0.51 with the annotation position: 0.12 [0.06, 0.20] apart, unresolved by the step reading.
- The annotation position alone, in place of Z, does almost as well as beside it (0.55; gain +0.13 [0.07, 0.20],
  unresolved).

**MHC-I in co-culture (H2: refuted).** The repaired output falls from 0.36 to 0.32; the gain, 0.04 [0.03, 0.05], is
no step. The reader of the joint stage reads 0.25 against 0.32 at Z (+0.07, no step). The grade stays relative.

**The controls, and the grades.** None loses a step (3 of 3), and every target has at ZC at least the grade it had
at Z (5 of 5). For PD-L1 under IFN-γ the grade goes from relative to absolute: the repaired output is supported
(0.00 [−0.07, 0.07]), the reference is supported, gate 3 passes. The change itself is small (0.04, no step); the
grade changes because the interval now lies inside the margin. With the input view on 23 November this was the target
that lost its grade. For the four pathway genes in THP-1 the reader of the joint stage is now supported
(0.03 [−0.05, 0.09]; it was unresolved at Z); gate 3 still does not pass there.

**The checks.** The readings with the head reading Z repeat the first outer step, and the readings with the head
reading ZV repeat the run of 23 November, for all five targets (largest difference of an arm 3e-14). The states of
the melanoma cells without annotations equal the recorded states to 7.2e-7. The joint stage of 23 November was
rebuilt exactly.

## What it means

- **The lead holds as far as it was a lead.** The annotation position carries part of what the mean loses, and a
  repair can use it: a step for IRF1. It does not make the loss go away: the reader of the joint stage is still at
  0.29, the repair at 0.51.
- **It is the first wider stage that could be given to every type.** The input view helped IRF1 more and cost another
  target its grade. The annotation position helps IRF1 less and costs nothing in these five targets; one grade
  improves. In the loop's terms it is a candidate change of what the head reads for all types, where the input view
  was a change for the types with a loss placed before the head. To become one it has to be declared in the
  adjustment between rounds and then read on cells that have not been used.
- **Where the item stands now.** Three things are left between IRF1 and a closed item: the repair made from one
  record stops 0.22 short of what a reader with the event's cells reads (the head, and how
  much one record can move it); that reader itself is at 0.29 where the input view reads 0.18 (the encoder); and
  with no record at all the readout leaves 0.82 (nothing made after training reaches an event it has no record of).
- **For the trial that trains the last layer** (not run; the owner decides): it now has a baseline that costs no
  training (0.51 with the record, 0.82 without) and a target inside the model. What only training can change is the
  last of the three: an event that is readable without its own record.

## What the model says about the melanoma cells

Described, not judged ([`results/stages_manifest.json`](results/stages_manifest.json)). Given the placeholders, the
model's output layer names for the cells of the co-culture stratum: disease `cancer_disease` for 89%, cell type
`epithelial_cell` for 99%, sex `male` for 78%, tissue `respiratory_tissue` for 67%. Under IFN-γ: `cancer_disease`
83% and `infectious_disease` 14% (7% in co-culture).

## The forecasts

Intervals: 8 of 8 held the value. Outcomes: for 9 of 11 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)).

| missed | forecast | observed |
|---|---|---|
| IRF1: the gain of the repair against the step margin | unresolved 0.45 | a step (given 0.20) |
| IRF1, the repaired output: the annotation position against the input view | loss 0.50 (the input view ahead by a step) | unresolved (given 0.40) |

Both misses are on the same side: the annotation position did somewhat better than expected.

## What was not as written

1. **The symbol map was read from another folder.** On the owner's machine the read of the map in the checkout did
   not return, and the script stood still after the counts. The copy of the same file made for the run of
   23 November was given with a new option `--maps`. No state existed yet; the commit of that change follows the
   declaration.
2. **Two grids were also made in the cloud container**, because the owner's machine was under heavy load for part
   of the run. All five grids used for the reading are the owner's machine's. The two made twice (IRF1, and the four
   pathway genes in THP-1) agree to 1.5e-7.

The scripts were tried on made-up values in a separate folder before the grids were made. Nothing had to be changed.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the declaration: the stages, the rule, the checks, the forecasts |
| [`results/tokens_manifest.json`](results/tokens_manifest.json), [`results/stages_manifest.json`](results/stages_manifest.json) | the tokens of the melanoma cells; the stages, their scales, the checks, what the model says at the placeholders |
| [`results/round.json`](results/round.json) | the rounds with the head reading Z, ZC, the annotation position alone, and ZV; the gains; the grades |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/tokens.py` | run on the owner's machine: RNA counts and tokens of the melanoma cells |
| `scripts/build_zc.py`, `round.py`, `score.py` | the arrays and SPEC files; the rounds; the score |

To run again (the arrays of the run of 23 November, the work folder of `encoder_stages_20261124` and the states of
`annotation_input_20261124` are needed):

    python scripts/tokens.py --repo TEDDY_MM --ckpt TEDDY_G_70M --strata STRATA.npz --out TOK
    python ../annotation_input_20261124/scripts/annot.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOK/tokens.npz --recorded TOK/recorded.npz --out ANNOT_MEL --modes N,P
    python scripts/build_zc.py --thp1 THP1_ENC.npz --mel STRATA.npz --annot-thp1 ANNOT_THP1 --annot-mel ANNOT_MEL --ckpt TEDDY_G_70M --work WORK
    for s in WORK/*.yaml; do python -m anm_diag scope $s --out OBSERVE/$(basename $s .yaml); done
    python scripts/round.py --work WORK --observe OBSERVE
    python scripts/score.py

## What is next

- **Declare the stage for a round.** "The head reads Z and the first annotation position" as an adjustment between
  rounds, for every type, and then one reading on units that have not been read. That is the owner's decision, as
  every amendment is.
- The trial that trains the last layer, read by the arms that ask whether an event is readable without its own
  record, against the baseline of this run.
- A second event that is lost before the head, when one exists: this is still one event.
