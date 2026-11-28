# A check with a known answer: does a knockout compose with IFN-γ by addition? — 27 November 2026

Run 51. Asked by the owner. Nothing was trained. The design and the pass line are in [`DESIGN.md`](DESIGN.md),
committed before any value was read (`f1578f4`); the scripts were rehearsed on made-up values and committed before the
reading (`bc58988`). Each part was read once and no rule was changed after it.

This is a check of a method against an answer known beforehand. Every cell used was read before (its 20 proteins in
batches 1 to 4 and in the final audit), and the five genes were a declared type of the frozen system from the start;
only the RNA had never been read gene by gene. It is not a finding on data never read and it does not count as a
confirmatory run.

## Outcome: passed

Knockout × IFN-γ, 95 knockouts read before in both conditions, 3,303 genes and 20 proteins, cells in two halves.

| | holds if | came out | |
|---|---|---|---|
| **P1** | all five of IFNGR1, IFNGR2, JAK1, JAK2, STAT1 are refused | 5 of 5 | **holds** |
| **P2** | at least 90% of the other 90 knockouts carry over | 90 of 90 | **holds** |
| **P3** | one class from the witnesses holds all five and at most two others | the five, and no other | **holds** |

Against the known pathway: **recall 5 of 5, precision 5 of 5**, for the refusals and for the class.

The table, one row per knockout: [`results/table_ifn.tsv`](results/table_ifn.tsv) (verdict, the strongest readout
against the end of the interval, the share of resamples in which it is refused, the witnesses, the share of the
signal removed, the class). Everything behind it: [`results/reading.json`](results/reading.json).

| knockout | verdict | strongest readout, over the end *w* | refused in resamples | witnesses: genes (lower / higher), proteins | share of the IFN-γ effect removed: RNA, PD-L1 |
|---|---|---|---|---|---|
| IFNGR2 | refused | IDO1, 3.33 | 1.00 | 56 (55 / 1), PD-L1 | 0.48, 0.63 |
| JAK2 | refused | IDO1, 3.21 | 1.00 | 42 (41 / 1), PD-L1 | 0.50, 0.65 |
| STAT1 | refused | IDO1, 3.13 | 1.00 | 38 (38 / 0), PD-L1 | 0.54, 0.69 |
| IFNGR1 | refused | IDO1, 2.82 | 1.00 | 45 (43 / 2), PD-L1 | 0.52, 0.62 |
| JAK1 | refused | IDO1, 2.44 | 1.00 | 39 (38 / 1), PD-L1 | 0.45, 0.68 |
| the other 90 | carry over | 0.52 to 0.88, median 0.60 | at most 0.18 | none | −0.03 to 0.04 (95%), median 0.00 |

The end of the interval: *w* = 4.00 errors (no readout whose halves disagree in sign reached 4; over the resamples of
knockouts 4.00 to 4.17). Nearest among the others: HLA-H (0.88 of *w*; refused in 3.5% of the resamples) and STAT3
(0.85; 18%).

**The witnesses.** 69 genes are a witness of at least one of the five and 19 of all five: B2M, CD74, CXCL11, ETV7,
GBP1, GBP2, GBP4, HAPLN3, HLA-E, ICAM1, IDO1, IFI16, IRF1, LAP3, PSMB8, RARRES3, RNF213, SECTM1, VAMP5. They are
genes IFN-γ raises, and under the knockout the rise is missing. One gene goes the other way in four of the five
(PMEPA1). The one protein witness is PD-L1, in all five (it stands outside against the isotype antibodies as well).
MHC-I and HLA-E are no witnesses here: IFN-γ changes one protein of the 20 two-fold in these cells, PD-L1 (and 106
genes, 103 of them upward).

**The class.** Grouped by `events.propose_types`, unchanged, on the 70 witness readouts: one class, the five, with
cosines of 0.96 or more between its members (median 0.98). No type was declared by hand.

**The five remove about half of the IFN-γ effect, not all of it.** On the log of a group's mean, over the 106 genes
IFN-γ changes two-fold, the five remove 0.45 to 0.54; on PD-L1, 0.62 to 0.69. Explored after the reading
([`results/explored.json`](results/explored.json), `scripts/explore.py`; it describes): cell by cell on PD-L1, 57% to
81% of the cells of the five look as if no signal had been given (18% of the control cells under IFN-γ do, and 17% of
the cells of the other knockouts). A group of cells that carry a guide holds cells in which the response is gone and
cells in which it is not; the mean of the group shows the mixture.

**The pairs left out, as declared.** Of the 14 pairs the graph declares as direct and of each knockout's own gene, two
stand outside in both halves: CD274 → PD-L1 (−6.1, −7.9 errors) and STAT1's own transcript (−4.4, −6.6). Both are
products IFN-γ raises: there is more to remove under the signal. This is why they were left out and judged apart.

**The sum over all readouts** (*T*; reported, it decides nothing): for the five it is 16 to 20 times its end τ. Two
other knockouts exceed τ with no readout outside: CD274 (1.29) and SEC11C (1.06).

## Secondary: knockout × co-culture

Described; no pass line. 56 knockouts read both with no signal given and in co-culture, the same rule
([`results/table_co.tsv`](results/table_co.tsv)). Refused: **the five and STAT3**; one class holds the six. The
witnesses of the five are CD74, HLA-DRA, B2M, HLA-C and other genes of antigen presentation, and the proteins MHC-I
and PD-L1 (CD47 for two). STAT3 has two witnesses (CD74, B2M), removes 0.07 of the effect of the co-culture and is
refused in 99% of the resamples; under IFN-γ alone it carries over. The five remove 0.40 to 0.48 of the effect of the
co-culture on RNA. The other 50 carry over (at most 0.88 of the end).

## Question 2: what TEDDY carries of the part that does not add

[`results/teddy.json`](results/teddy.json). For the five, by `account.account`, unchanged, in each condition, with the
same halves. With no signal given the five change none of these readouts by four errors (four gene pairs have an
account there, with changes of 0.005 or less), so the part that does not add is the change under IFN-γ.

Pooled over the five (the interval by resampling the five knockouts; the intervals of single pairs are wider and are
in the file):

| the part that does not add | at the input (V) | at the state Z | read from Z by a readout of control cells | at the head's output as trained |
|---|---|---|---|---|
| **PD-L1** (5 pairs) | 0.94 [0.90, 1.00] | 0.93 [0.90, 0.96] | 0.87 [0.83, 0.91] | **0.09 [0.09, 0.10]** |
| the ten strongest RNA witnesses of each (50 pairs) | 0.88 [0.87, 0.90] | 0.83 [0.80, 0.86] | 0.74 [0.70, 0.77] | the head has no output for a gene |

- **The state carries it; the head as trained does not show it.** Nothing is lost before the state for PD-L1 (0.94 at
  the input, 0.93 at Z), and a readout fitted on control cells alone reads 0.87 of it from the state. The head's own
  PD-L1 output moves by 0.09 of what the protein does (single pairs: 0.09 to 0.10, each [0.07, 0.12]).
- A readout taught by the other four members of the class reads the fifth: 0.99 [0.93, 1.05] for PD-L1 and
  0.99 [0.98, 1.00] for the genes. The class the witnesses gave is one a readout can use.
- For a gene, V is a view of the input (512 components); the input itself holds the gene.
- MHC-I: the five change it by 0.06 to 0.21 under IFN-γ, and two of them by four errors (JAK2, STAT1). For those two
  the head's output holds 0.06 and 0.09, the state 1.46 and 1.03. Two pairs give no interval over knockouts (the one
  in the file is not meaningful). HLA-E: no change of four errors, no account.

## The rehearsal, and what it changed

[`results/rehearsal.json`](results/rehearsal.json): six fresh made-up worlds through the same scripts (`rehearse.py`).
In the three with planted knockouts every planted one was refused (the five, one that removes half of the signal, one
with a single gene), none that adds was refused, and the class was the planted one; in the three with nothing planted
nothing was refused. `teddy.py` was run on a planted world for faults of the code.

Two things were changed after the first made-up worlds and before any value was read; both are written into
`DESIGN.md` with their reason. The departure of the genes is centred over the genes of a group (RNA levels are
relative). The verdict is by readouts, and the sum *T* is reported: on worlds with nothing planted it refused 0 to 6
of 40 knockouts, and a refusal through it names no readout. The pass line was not changed.

## What this shows and what it does not

- The rule finds the known answer in this data: the five and nothing else, with a witness for each and a class found
  from the witnesses.
- It was the easy direction. The five were the largest effects in this library and known beforehand; a knockout that
  removes a tenth of a response would stand near the end of the interval, as STAT3 does in the co-culture. The 95
  knockouts hold no second kind of known answer.
- "Carries over" means: no readout departs by more than four errors in both halves, with 46 to 171 cells per
  half. It is not a statement that the departure is zero.
- For TEDDY the check repeats, for the part that does not add, what the first sealed runs found for the whole
  change: the loss is at the head as trained, not before the state.
- One study, one model, knockouts.

Not opened: the 188 units never read and the 73 reserve donors. The unread units with no signal given could carry
the same rule once more for "the others carry over" (the five are used up, so recall cannot be asked again); that
needs the owner's word.

## Reproduce

```
python scripts/units.py     --cells CELLS.tsv --library LIBRARY.npz [--ifn-list OUT.tsv]
python scripts/counts.py    --repo TEDDY_MM --cells OUT.tsv --out DIR            # on the machine with the RNA file
python scripts/rehearse.py  --work DIR                                           # made-up worlds only
python scripts/summarize.py --cells CELLS.tsv --library LIBRARY.npz --counts MEL_COUNTS... --counts-ifn DIR/counts_ifn.npz --packs PACK... --out summary.npz
python scripts/analyze.py   --summary summary.npz                                # once
python scripts/teddy.py     --cells CELLS.tsv --library LIBRARY.npz --states MEL_STATES... --counts MEL_COUNTS... --counts-ifn DIR/counts_ifn.npz --packs PACK...   # once
```

`summary.npz` (sha256 `fda9322ec844dfef…`) holds the sums by group and half and is not in the repository.
