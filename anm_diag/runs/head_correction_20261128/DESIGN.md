# Correcting the one part the loop has never corrected: the head

28 November 2026. Run 52. Asked by the owner on 28 November. This document and the lot
([`results/units.json`](results/units.json)) are committed before any cell of a drawn unit is listed, before any
state of one is computed and before any value is read. The scripts that fit and judge are rehearsed on made-up
values and committed before the reading. The reading is made once; no rule is changed after it.

## 1. Why the head

The loop has corrected the graph and the reading, and the correction held on data never read (run 49). It has never
corrected the model. Two places where the model loses are now measured, and they are not the same place:

| | where | what was measured |
|---|---|---|
| the mean over gene positions | before the state | IRF1 → PD-L1 in THP-1: the state carries 0.52 where the input carries 0.89. Training the pooling: 0.66, unresolved (run 50) |
| **the head** | after the state | of the part of PD-L1 that does not add with IFN-γ, the state carries 0.93 and the head's output 0.09 (run 51) |

The second is the easy one, and this run takes it first: the information is in the state, since a readout taught by
control cells alone reads 0.87 of it from there. The question is not whether it can be recovered but whether the
loop's own rows teach anything beyond a plain refit, and whether a correction made here holds on units never read.
**The author says this plainly before running it: C1 below is the easy half.**

## 2. What is changed, and what is not

The head is an MLP on the unit-length state, six layers, ending in a layer of 134 rows, one per protein of the panel
it was trained on. **Only three rows of that last layer move**: CD274, HLA-A-B-C and HLA-E, which are the three the
loop reads as PD-L1, MHC-I and HLA-E. 3 × 512 + 3 = **1,539 numbers**. The encoder, the pooling, the five layers
before it and the other 131 rows are untouched. Nothing else in the loop changes: version 6 as frozen at `b4f6196`,
the graph as the final audit left it, `account.account` unchanged.

## 3. The cells

**The units.** 93 units drawn by lot (seed 20261128) from the 188 the census never read: 47 in the co-culture, 46
with no signal given. **The other 95 stay closed**, for whatever the loop asks next. No unit of the final audit and
no reserve donor is in the pool.

**The control cells** of each melanoma condition are split in two by a hash of the cell's name,
`sha256("head|" + name)`, even and odd, declared here:

| half | used for |
|---|---|
| teach | fitting the three rows (every arm) |
| judge | the account: the reference the drawn units are read against, and the readouts the account fits |

No arm ever sees a judge cell, and no drawn unit is seen by any arm. The three proteins are read as everywhere here:
the centred log-ratio over the 24 antibody columns.

## 4. The arms

Each arm gives the three rows. The target of a fit is the measured readout, standardized by its mean and spread among
the teach control cells; the head's output is standardized the same way, so the fit is in the space the account
measures in.

| | the rows are fitted on |
|---|---|
| **H0** | nothing: the head as shipped |
| **H1** | the teach control cells alone |
| **H2** | the teach control cells, **and the loop's rows**: for every event already recorded in these two conditions, the shift the head's output shows must equal the shift the loop measured |
| **H3** | the same as H2 with the labels of the events shuffled (described; no pass line) |

The recorded events are the units of batches 1 to 4 and of the final audit in these two conditions, with their
measured shifts; the drawn units are not among them. H1 and H2 are fitted by the same procedure with the same
regularization, chosen on made-up values before the reading; the only difference is the rows.

**The shuffled arm is informative here, and was not in run 50.** In THP-1 every recorded event lowered PD-L1, so
shuffled rows asked for nearly what the true rows asked. In these two conditions the recorded events move the three
readouts up and down and by different amounts, so a shuffled set asks for something else. This is checked before the
arm is run (`results/rehearsal.json`: the correlation of the true and shuffled right-hand sides) and reported.

## 5. What is asked

The judge is not the loss. It is `account.account`, unchanged, on the drawn units: for every unit and each of the
three readouts whose change is at least four errors, the share of that change the arm's output carries
(`as_trained`), with the judge control cells. Pooled over pairs; the interval by resampling the drawn units 2,000
times (seed 20261128).

| | hypothesis | holds if |
|---|---|---|
| **C1** | **A correction made at the head closes the loss, on units never read.** | the best of H1 and H2 carries at least 0.60 where H0 carries at most 0.25. *Refuted* below 0.40, or if H0 is already above 0.40. |
| **C2** | **The loop's rows teach something a plain refit does not.** | H2 − H1 above zero with its 95% interval clear of zero. *Refuted* if the interval lies below zero. Between: unresolved. |
| G1 | The state did not move. | the share the state carries (`carried` at Z) is the same for every arm to within 0.02 — the correction must not change the state, only what is read from it |
| G2 | Nothing got worse where nothing was wrong. | on the drawn pairs whose change is under four errors, no arm's output moves by more than 1.5 spreads; and each arm reads the judge control cells with an error no larger than H0's |

**The run is passed if C1 holds.** C2 is the loop's own question and is reported with its verdict either way; a
refuted C2 is a result, not a failure of the run: it would say the head was simply never fitted to this study and
the loop only pointed at the place.

What the author expects, written before anything is opened: C1 holds 0.85, refuted 0.05. C2 holds 0.2, unresolved
0.5, refuted 0.3. G1 holds 0.95; G2 holds 0.8. H3 within 0.05 of H1 (0.6).

## 6. What can be claimed if it passes

For one study, one model and knockouts: the loop located a loss in the model, a correction of 1,539 numbers aimed at
that place closed it, and the correction held on units no part of the work had seen. Not: that the loop can correct
a loss that lies before the state (that is the open one, and run 50 left it unresolved); not another study, another
model, donors, or time.

## 7. Order

1. This document and `results/units.json`: committed and pushed.
2. On the owner's machine: the cells of the drawn units (names only), their states, and their antibody columns.
3. `scripts/fit.py` and `scripts/read.py`, rehearsed on made-up values (`results/rehearsal.json`), committed and
   pushed.
4. The reading, once.
5. The write-up, as it comes.
