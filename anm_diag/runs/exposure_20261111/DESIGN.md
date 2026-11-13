# Does the head lose IRF1 because of what it was trained on? Design, written before the runs

Status: exploratory. This file was committed before any of the runs below were made; the commit history shows the order.

## The question

In THP-1 cells the PD-L1 response to an IRF1 knockout is lost at the head's output (u = 1.07, rejected), far more than for
the four upstream knockouts (0.28 to 0.60). The analysis in `runs/split_20261111` located the shrinkage at the head's last
step and offered a reading that was not tested: the head, trained on unperturbed bone-marrow cells, passes a knockout that
switches the whole IFN-γ programme off and shrinks one that leaves a small footprint. If that is right, a head that has seen
small-footprint knockouts in training should stop losing the IRF1 response, and a head that has only seen the upstream
knockouts should not.

## Data and split

THP-1 cells of Papalexi 2021, stimulated, the same cells and exported features as `runs/split_20261111`
(7,478 cells: controls and guides against IFNGR1, IFNGR2, JAK2, STAT1, IRF1). TEDDY is frozen; only the head changes.

Every group of cells (one guide in one batch; control cells of one batch) is cut in two at random, split seed 0:
a training half and an audit half. The head never sees an audit cell. A fifth of the training half, taken the same way,
is set aside to choose when to stop.

## Three heads

| arm | head | training cells |
|---|---|---|
| `as_trained` | the head as it is | none |
| `without_IRF1` | fine-tuned | controls and the four upstream knockouts |
| `with_IRF1` | fine-tuned | the same plus IRF1 knockout cells |

Fine-tuning: all weights of the head (latent and decoder networks), full batch, Adam with learning rate 3e-4, at most 200
steps, the weights kept are those with the lowest loss on the set-aside fifth, checked every 20 steps. Loss: squared error of
the head's PD-L1 output (CD274) against measured PD-L1, the measurement first rescaled to the mean and spread that this output
has on training control cells; plus, with weight 1, the mean squared change of the other 133 outputs from the head as trained,
so the rest of the output stays where it was. No dropout. The learning rate was chosen between 1e-4 and 3e-4 on the loss of
the training half alone (60 steps); no audit cell was looked at. Nothing else is tuned.

## The audit

Six runs with the branch's `anm_diag`: three heads by two targets, IRF1 and JAK2 (an upstream knockout present in the
training cells of both fine-tuned heads). Each run is the `papalexi2021_target_*` declaration of `runs/split_20261111`
restricted to the audit half, stages Z, head latent, head output, readout PD-L1, resampling by guide, same folds,
bootstrap, margin 0.10 and seed.

Primary outcome: u at the head output for IRF1 and PD-L1.

## What would count as what

Let u0 be the primary outcome for `as_trained` on the audit half.

- **The reading is supported** if `with_IRF1` is no longer rejected at the head output and its u is at most half of u0,
  while `without_IRF1` is still rejected with u above half of u0.
- **The reading is not supported** if `without_IRF1` already meets that mark: then seeing THP-1 cells and this readout was
  enough, and the loss was a matter of the head's training domain, not of the knockout's footprint.
- **Not answered** if `with_IRF1` is still rejected with u above half of u0: this small fine-tune does not recover the
  response, which by itself says nothing about why.
- If u0 is not rejected on the audit half, the question cannot be put on these cells and that is the result.

JAK2 is expected to improve alike under both fine-tuned heads. It is reported and does not enter the rule.

## Also reported, outside the rule

For split seeds 0 to 2 (the audit is run for seed 0 only): the shift of the head's own PD-L1 output between knockout and
control cells of the audit half, in control standard deviations, next to the measured one, for all five targets; and the
group shift over the spread of control cells at the head latent and the head output.

## Limits known in advance

One data set, one split for the audit, four guides per target as resampling units. Training and audit cells share guides,
so this asks whether a head that has seen these knockouts carries their response in new cells, not in new knockouts.
`with_IRF1` trains on more cells than `without_IRF1`. The fine-tuned heads are trained on the audited readout; that is
the point of the test, and it is why the audit cells are kept apart.
