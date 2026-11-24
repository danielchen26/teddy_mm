# The head reads the first annotation position beside Z: declaration

Written on 24 November 2026, **before the states of the two melanoma strata were computed with the placeholders and
before any grid or reading of this run was made**.

Everything here is exploratory. The cells are development cells that earlier runs have read. No final unit is read
and no model is trained.

## 1. Why this is asked

Two runs stand behind this one.

- On 23 November the head was given 512 components of the model's input beside Z
  ([`input_beside_z_20261123`](../input_beside_z_20261123/)). For IRF1 in THP-1 the repaired output went from 0.68
  to 0.39, a step; it did not close, and given to every type the wider stage cost PD-L1 under IFN-γ its grade.
- Earlier on 24 November the encoder was run with TEDDY's four annotation positions
  ([`annotation_input_20261124`](../annotation_input_20261124/)). The state Z did not change. The annotation
  positions themselves, a summary the model learned in pre-training, read IRF1's response better than the mean over
  the genes: the first position alone leaves 0.29 where Z leaves 0.40. That difference was unresolved by the rule,
  it was one event, and the forecast had gone the other way.

That lead is tested here in the way the plan tests a repair: in the round, with the repair the loop uses, for the
same five targets and with the same cells, halves and resamples as the run of 23 November. The question: **if the
head reads the first annotation position beside Z, does the repaired output for IRF1 gain a step, and does any other
target lose by it?** It needs no training of the encoder, only the encoder run as its official code runs it.

## 2. The stages

For every cell of the three strata of the knockout screens (THP-1, melanoma in co-culture, melanoma under IFN-γ):

| stage | what it is |
|---|---|
| `Z` | the state the head reads today |
| `C` | the output of layer 12 at the first annotation position, with the four placeholders in front of the genes. It is what the official class returns as `cell_emb` |
| `ZC` | Z at unit length beside C, each block divided by the square root of its total variance among the control cells of the stratum: the recipe of the run of 23 November, by the same function |
| `ZV` | Z beside the 512 components of RNA: the joint stage of 23 November, rebuilt, for a paired comparison |
| `V` | the 512 components alone (for the two targets that had them on 23 November) |

The states with placeholders of the THP-1 stratum exist. Those of the two melanoma strata (11,171 cells) are computed
for this run on the owner's machine by the script of the run before, which was checked against the official class.
No readout and no label of an event enters a stage.

## 3. What is read

The five targets of 23 November, unchanged:

| | target | events | read as |
|---|---|---|---|
| T1 | PD-L1, THP-1 | IRF1 | on its own, with its own record |
| T2 | MHC-I, melanoma, co-culture | IFNGR1, IFNGR2, JAK1, JAK2, STAT1 | a pool |
| C1 | PD-L1, THP-1 | IFNGR1, IFNGR2, JAK2, STAT1 | a pool: control |
| C2 | PD-L1, melanoma, co-culture | the five pathway genes | a pool: control |
| C3 | PD-L1, melanoma, IFN-γ | the five pathway genes | a pool: control |

For every target the audit tool reads all stages in one run, with ridge and boosting. Then the round is run with the
general code (`anm_diag.cycle`, the repair `cycle.LinearRows`) with the head reading `Z`, `ZC`, `C` and `ZV`: the same
cells, the same halves (seed 0), the same resamples. Margin 0.10, step margin 0.10, 1,000 resamples. Claim, protocol
and reader specification are those of package 2. Nothing is amended.

## 4. What counts as what

`gain` of an arm = |u| with the head reading Z minus |u| with the head reading the other stage, on paired resamples:
**a step** if its interval lies above 0.10, **no step** if it lies below 0.10, unresolved otherwise. The rule is the
rule of 23 November with `ZC` in the place of `ZV`.

- **H1. IRF1 closes with the head reading ZC.** Confirmed if the grade of the round at ZC is relative or absolute
  and the gain of the arm `rows` is a step. Refuted if that gain is no step. Otherwise unresolved.
- **H2. The same for MHC-I in co-culture**, with the pool.
- **C1 to C3. No harm.** For each control target the arm `rows` at ZC loses no step more than at Z.
- **Grades kept.** Every target has at ZC at least the grade it has at Z. This was not a rule on 23 November and is
  one now, because the input view cost one target its grade without failing "no harm".
- **Side, by the same rule**: the head reads `C` alone.
- **The annotation position against the input view**: |u| at ZC minus |u| at ZV for every arm, by the tool's step
  reading. A loss means the input view beside Z removes a step more than the annotation position does.
- **Checks.** The readings at Z repeat the first outer step, and the readings at ZV repeat the run of 23 November,
  each to 0.005. The states without annotations of the melanoma cells equal the recorded states to 1e-4. If a check
  fails: stop and report.

What each outcome would mean:

| outcome | reading |
|---|---|
| H1 confirmed, controls hold, grades kept | the head should read the annotation position; this loss needs no training |
| the repair gains a step and the round does not close | the same position as the input view on 23 November: part comes back in front of the head, the record is what is short |
| no step | the lead of the run before does not carry into a repair. What the head can be given without training is then the input view, and the item stays with the encoder |
| a control loses, or a grade is not kept | not a stage for every type, as the input view was not |

## 5. The order

1. Commit this file, the forecasts (`results/declared.json`) and the scripts.
2. `scripts/tokens.py`, then `annot.py` of the run before, on the owner's machine: the states of the melanoma cells
   without annotations (the check) and with the placeholders.
3. `scripts/build_zc.py`: the arrays and SPEC files. No readout enters.
4. `anm_diag scope` for the five SPEC files, on the owner's machine.
5. `scripts/round.py`: the rounds and their paired differences, once. A result file that exists is not overwritten.
6. `scripts/score.py`: the forecasts against what was read.

The scripts are tried once on made-up values in a separate folder before step 4.

The tokens of the melanoma cells (step 2, first half) were being made while this file was written. They are the
model's input and the same for any run; no state with annotations existed.

## 6. The forecasts

`results/declared.json`, 19 of them. Briefly, and why:

- **IRF1.** A reader of ZC with IRF1's own cells is expected at 0.24 to 0.36 (the four positions beside the state
  read 0.31 in the run before). That is about half of what the input view gave on 23 November, so the repair is
  expected to gain less than it did then: 0.45 to 0.72 against 0.68 at Z. H1: unresolved 0.50, refuted 0.42,
  confirmed 0.08. The input view beside Z is expected to stay ahead of the annotation position (0.50).
- **MHC-I in co-culture.** No step is the most likely outcome (0.60).
- **Controls.** No harm is expected for all three. Every grade kept: 0.60.

## 7. What would change the plan

- A check of Section 4 fails: stop and report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.
