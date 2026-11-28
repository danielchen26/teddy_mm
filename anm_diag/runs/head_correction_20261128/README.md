# Correcting the head — 28 November 2026: the run did not answer its question

Run 52. Asked by the owner. The design, the pass line and the lot are in [`DESIGN.md`](DESIGN.md), committed before
anything was opened (`adfad6d`, amended at `c96fd6d` for how the fit is carried out, before the reading). The arms
were fitted and the judge read once, as written.

## What happened

**Of the 279 pairs the judge had (93 drawn units × 3 readouts), one had a change of at least four errors.** The
median change was 0.77 errors, the largest 4.17. With one pair there is nothing to pool, so **C1 cannot be read**:
the run is neither passed nor refuted. The guards hold (G1, G2).

| readout | pairs | median change, in errors | at least four errors |
|---|---|---|---|
| PD-L1 | 93 | 0.75 | 1 |
| MHC-I | 93 | 0.82 | 0 |
| HLA-E | 93 | 0.76 | 0 |

The one pair: HNRNPA1 in the co-culture, PD-L1, change −0.271 (4.2 errors). The head as shipped carries −0.00 of it,
H1 +0.14, H2 +0.18, H3 +0.10. One pair says nothing.

## Why, and whose fault

The head has an output for three proteins. A knockout drawn by lot does not move them. **This was in the final
audit's own results, which were read before this design was written**: of its 133 melanoma units, 2% reached four
errors on PD-L1, 2% on MHC-I and none on HLA-E. The expected yield of 93 units × 3 readouts was about four pairs.
The design was written without that count. It is a fault of the design, not of the data or the code.

What it cost: the 93 drawn units are opened and cannot serve as sealed data again. **The other 95 are untouched.**

## What was established anyway

The fit uses only the teach control cells and the events recorded before; none of it depends on the judge.

- **A row of the head's last layer can read PD-L1 from the activation.** Fitted on one half of the teach control
  cells and measured on the other, the squared error is 0.597 where the shipped row gives 1.887 — which is what
  "predict the mean" gives. So the information does reach the last layer and the shipped row does not take it. This
  is the premise of C1, and it holds.
- **The loop's rows carry information that generalizes to events the fit was not given.** The weight of the rows was
  chosen by splitting the 197 recorded events in two: the error on the events the fit was not given falls as the
  weight rises (MHC-I 0.346 → 0.096 from weight 0.03 to 3).
- **The shuffled control is informative here**, unlike run 50: the shuffled shifts correlate −0.03 (PD-L1), −0.11
  (MHC-I) and −0.04 (HLA-E) with the true ones.
- **The judge itself measures what it is meant to** ([`results/rehearsal.json`](results/rehearsal.json)): handed the
  readout itself as an output it gives 1.000, handed noise −0.001, handed half of each 0.441.

## What to do instead, and what it costs

The judge needs units that move the three readouts. The one set that moves them and that no part of this fit
touched is the **IFN-γ library**: the fit used control cells and recorded events of the co-culture and the control
condition only, and never a cell under IFN-γ, where PD-L1 moves by 12 to 18 errors for the events of the pathway
(run 51). Judging there costs nothing sealed, and the answer is retrospective: it can show whether a correction at
the head works, not confirm it on data never read. A sealed confirmation afterwards would use the 95 units still
closed, with the units chosen by a written rule on what the frozen system forecasts rather than by lot.

That is a proposal; it waits for the owner.

## The lesson

Before data is drawn, count how many pairs the judge will have. A pass line that cannot be reached by the data it is
written for is not a pass line. The count was available in a run already finished.

## Files

`results/units.json` the lot; `results/rehearsal.json` the made-up worlds and the check of the judge;
`results/reading.json`, `results/table.tsv`, `results/outcome.json` the one reading; `scripts/` the four arms, the
judge, and what ran on the owner's machine. The fitted rows (`rows.npz`, sha256 `e2bf2b9470d6e925…`) and the states
are not in the repository.
