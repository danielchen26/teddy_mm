# Training by the TEG, against four other kinds of knowledge: a minimal example. Design, written before the fits

Status: exploratory. This file and the scripts were committed before any head of this run was fitted beyond a 20-step
timing check on training cells, and before any audit cell was read for any arm of this run. The commit history shows the
order.

## Why this run

The head was repaired before ([`runs/loop_guard_20261112`](../loop_guard_20261112/)). That repair was guided by the
audit: it changes only the step the audit names (the head's MLP), on the cells the typed events admit. Its loss, however,
is not written by the TEG. It is a squared error to the measured protein, cell by cell, plus two terms that hold the other
outputs. So the question "what if the head is trained by the TEG" has not been asked.

It also left something open. By the rule of that round the loop closed (a readout fitted on the head's output reads the
PD-L1 response again). The head's **own** PD-L1 output was read later
([`runs/scope_20261116`](../scope_20261116/), `results/frozen_arms.json`), on the same cells:

| | head as trained | repaired (`mlp_only`, round 1) |
|---|---|---|
| melanoma | 0.92 [0.90, 0.93] | 0.29 [0.22, 0.35] |
| THP-1 | 0.97 [0.96, 0.99] | 0.43 [0.37, 0.52] |

*u* of the head's own output, margin 0.10. The repair removed most of what was unexplained, and by its own output the head
is still rejected on knockouts of genes it never saw. That is the head "not trained well" of this run, fixed before any
result of this run exists: **the PD-L1 response to knockouts of the IFN-γ pathway, at the head's own output, on genes held
out of the fit, in melanoma and in THP-1.**

Two questions:

1. Trained by the rows of the TEG, is the head's own output enough on held-out events?
2. With the same cells and the same budget, what do four other kinds of knowledge achieve: the class of a cell (an
   ontology), a signed edge (a knowledge graph), the geometry of a graph (a KG embedding), the judgments of a learned model
   (a decision engine)? And the repair already made (the measured protein, cell by cell)?

## What is the same in every arm

- The cells, the split and the audit controls of `runs/loop_20261112` and `runs/loop_guard_20261112` (same code, same seed;
  the audit controls are checked against the record). One head per data set and held-out gene: 5 + 5 fits per arm.
- Only the three layers of the head's MLP change. The decoder stays bit for bit (asserted). TEDDY is frozen.
- Adam, learning rate 3e-4, at most 400 full-batch steps, a check every 20 steps on the set-aside fifth of the training
  cells; kept: the weights with the lowest value of the arm's own loss there; stop after three checks without improvement.
- Two holding terms, weight 1 each, as in `runs/loop_guard_20261112`: the other 133 outputs on the perturbation cells stay
  where they are; all outputs on 3,000 bone-marrow training cells stay where they are.
- **The control cells are fitted to their measured PD-L1, cell by cell, in every arm**, with the target of the earlier
  repair (the measured value, put on the scale of the head's output among the training controls). The arms differ only in
  what they carry about the knockout cells.
- No weight is tuned. No audit cell and no bone-marrow test cell enters a fit, a stopping check or a choice.

## The arms: what the supervision carries about the knockout cells

| arm | stands for | what it is given about the knockout cells of the fit | its term |
|---|---|---|---|
| `data` | the repair already made | the measured PD-L1 of every cell | squared error, cell by cell (with the control side, this is the loss of `runs/loop_guard_20261112`) |
| `class` | an ontology | which class every cell belongs to (knockout of a pathway gene, or control) | a linear probe on the head's latent has to tell the classes apart (cross-entropy) |
| `edge1` | a knowledge graph | a signed edge: the knockout lowers PD-L1. No size, no control | the shift of the own output is at least 1 SD below the controls (hinge, squared); `edge05` and `edge2` use 0.5 and 2 SD, since the graph gives no size |
| `graph` | a KG embedding | the signed edge, and the geometry of the graph: genes with the same edges are held alike | `edge1`, plus: the latent shifts of the events of the fit point the same way (1 − cosine) |
| `judge` | a decision engine | the probability of a learned judgment model that PD-L1 is low in this cell | cross-entropy of the own output, as a judgment, against that probability |
| `rows` | **the TEG** | its rows: for every event of the fit (one gene against the controls), the measured change | *u* of the head's own output for that event, squared, averaged over the events |
| `data_rows` | the repair with the TEG added | both | the terms of `data` and of `rows`, added |

Details, all fixed here:

- **`rows`.** For an event *e*, *u*<sub>e</sub> = 1 − *a* · (mean of the own output on its knockout cells − mean on the
  controls) / (measured mean on its knockout cells − measured mean on the controls), with *a* the slope of the measured
  value on the own output among the control cells of the fit. The own output is the predicted mean (the exponential of the
  head's output), as `anm_diag.scope` reads it. This is the statistic of the audit, written as a loss, on the events of the
  fit. It uses the measured protein of a knockout cell only through the mean of its event.
- **`judge`.** A logistic regression from TEDDY's embedding to the judgment "PD-L1 is below the median of the training
  controls", fitted on the cells of the fit. The head's output becomes a judgment by a logistic curve centred on the median
  of its output among the training controls, with their spread.
- **`graph`.** The graph a knowledge graph would hold gives every one of the five genes the same two edges (it regulates
  PD-L1; it is part of the IFN-γ pathway). An embedding fitted to a graph places entities with the same edges at the same
  point, so the similarity of every pair of genes is 1. No embedding is trained.
- **`class`.** The probe is trained with the head and is not part of it.

The four middle arms are **stand-ins**: the least that carries what that kind of system carries, used as the only thing
known about the knockout cells. They are not the best method of their kind and not a product. Nothing is known here about
how Jev works inside; `judge` is a learned judgment model with probabilities and nothing more specific.

## The reading (made once, after all fits)

For every arm and data set, with `anm_diag.scope` (`frozen` scope, as `runs/scope_20261116/scripts/frozen_arms.py`):

- **Pooled.** Every knockout cell goes through the head fitted with its gene held out; the audit controls are split among
  the five heads as in `runs/loop_20261112`. *u* of the own PD-L1 output, with its interval (1,000 resamples; knockout
  cells resampled by gene, controls by cell). Five genes are few clusters (rule `K1`): the intervals are approximate.
- **Events.** One held-out gene, its head, all audit controls; cells resampled. Ten events.
- **Plain, without the tool.** On the cells of a held-out event: the correlation of the own output with the measured
  protein; the shift of the own output and the measured shift, each in units of the spread among the audit controls.
- **Guard.** The head's own evaluation on the 5,000 bone-marrow test cells; the mean fall of the five heads of a data set
  against the head as trained (0.602).

Readings of *u*, two-sided at the tool's margin 0.10: **supported** (the whole interval inside ±0.10), **rejected** (the
whole interval outside), **unresolved** (otherwise).

## The rule

An arm **trains the head well** if, in both data sets: the pooled reading is supported; the mean fall on bone marrow is at
most 0.02; and the reading is not an artefact of the statistic, meaning the slope among the audit controls is positive and
the shift of the own output has the sign of the measured shift.

- **S1.** `rows` trains the head well, and none of `class`, `edge1`, `graph`, `judge` does.
- **S2.** Against the repair already made: |*u*| of `rows` minus |*u*| of `data`, pooled, on the same resamples, per data
  set. "Lower" needs the interval of the difference below zero in both data sets.

Outcomes: S1 holds or not, S2 holds or not; all four combinations are reported as they come. If `rows` is not supported,
the example has not shown what it set out to show, and it is reported so.

## Forecast

| arm | forecast for the pooled *u* |
|---|---|
| head as trained | 0.92 and 0.97, on record |
| `data` | about 0.3 and 0.4, rejected (on record for the unguarded repair: 0.29 and 0.43) |
| `class` | far from 0 (above 0.5 in size), rejected; its sign is not forecast |
| `edge05`, `edge1`, `edge2` | not supported; `edge05` leaves *u* above 0 and `edge2` drives it below 0; where `edge1` lands depends on a size the graph does not give |
| `graph` | as `edge1`; on THP-1 no better, because IRF1 differs from the other four and the graph holds them alike |
| `judge` | between `data` and the head as trained: 0.3 to 0.7 |
| `rows` | supported on melanoma; on THP-1 within ±0.2 and possibly unresolved, with IRF1 the worst event |
| `data_rows` | close to `rows` |
| guard | kept by every arm |

S2 is forecast to hold, for a reason that has to be said: `rows` is trained on the statistic the reading uses (on other
events), and `data` is not. The plain numbers are reported so that this can be seen from the other side; there `data` is
forecast to be as good as `rows` or better.

## What is done next, and what is not

- If S1 holds, the pooled audits of `runs/loop_guard_20261112` (the full chain, with fitted readouts) are run on the heads
  of `rows` and reported apart. They are not part of the rule.
- No further variant of any arm is fitted against these audit cells after the reading.
- The 73 reserve donors of the second question are not touched.

## Audit of this plan

| problem | what the design does about it |
|---|---|
| The target was chosen knowing that the earlier repair leaves the own output rejected. | Stated; that record is the reason for the run. The target, the events and the arms are fixed here, not after the reading. |
| These audit cells were read before, for the head as trained and for the earlier repairs. | Stated, with the numbers. They have not been read for any arm of this run. |
| `rows` is trained on the statistic it is read by. | It is trained on the events of the fit and read on events of genes it never saw; that is the claim under test. The plain numbers do not use the statistic. S2 names the asymmetry. |
| A head can satisfy the statistic with a wrong prediction (a negative slope among the controls and a shift of the wrong sign). The head as trained has a negative slope on THP-1. | The control cells are fitted cell by cell in every arm; the rule requires a positive slope and the measured sign. |
| The four stand-ins could be made to look weak. | Each is given exactly what its kind carries and nothing is tuned; the margin the graph cannot give is shown at three values; the forecast says in advance which arms cannot succeed and why. They are named stand-ins everywhere. |
| The stopping check of `rows` rests on few knockout cells per event (a fifth of about 200 to 1,300). | The same stopping rule serves every arm; the kept step and the curves are written out. |
| One split of the controls, one fit per gene; five genes per data set as resampling units. | As in the rounds before. Graded exploratory. The ten events are reported one by one. |
| Before this file was committed the scripts were run in `--check` mode. | It compares the head in torch with numpy, reproduces the recorded reading of the head as trained (0.915746 and 0.974302), and runs 20 steps of every arm with one gene held out, printing losses on training cells. No audit cell is read for a fitted head. |
| Only the head is trained. | Said in every summary: this is training a head on a frozen foundation model by the TEG, not pretraining a foundation model by it. |
