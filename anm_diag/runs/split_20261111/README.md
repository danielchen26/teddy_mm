# The pooled source, separated — 11 November 2026

Exploratory. Every run here was made after the pooled results were seen; the verdicts of the sealed runs stand.

The two sealed head-fix runs ([`frangieh2021_headfix/`](../frangieh2021_headfix/), [`papalexi2021_headfix/`](../papalexi2021_headfix/))
each declare one source that pools several things: five target genes on the melanoma data; five target genes and four
batches on the THP-1 data. Rule `S2` (see the package README, *Typed source events*) flags such pooling: a verdict for the
pool says nothing about each pooled value. This folder tests it on the measured side by running the same audit for each
target alone and, on THP-1, for each batch alone.

## What is here

| folder | what was run |
|---|---|
| `frangieh2021_target_<GENE>/` (5) | melanoma, one target as the treated class, the pooled run's 4,000 controls, guides as the resampling unit |
| `frangieh2021_pooled_quick/`, `frangieh2021_pooled_quick_byguide/` | the pooled source again at Z, head latent and head output (check of the environment; resampling by target and by guide) |
| `papalexi2021_target_<GENE>/` (5) | THP-1, one target, the pooled run's 2,386 controls, guides as the resampling unit |
| `papalexi2021_batch_rep{1,3,4}/` | THP-1, all five targets, one batch with its own controls, targets as the resampling unit (batch 2 has 73 cells and was not run alone) |
| `papalexi2021_pooled_quick/` | the pooled source again at Z, head latent and head output |
| `results/` | everything computed from the records (see below) |
| `scripts/` | the scripts that made the records and the results |

Each run folder holds `declaration.yaml` (the sealed design with only the treated class, the strata or the resampling unit
changed; status exploratory) and `anm_record.json`.

## Results

**The environment reproduces the sealed records.** Rerunning the pooled source differs from the sealed record by at most
0.0004 in *u* on melanoma and 0.00012 on THP-1, with identical verdicts.

**Melanoma, PD-L1: pooling holds.** All five targets respond in the same direction (−0.69 to −0.96); the head output is
rejected for every target (*u* 0.40 to 0.46); at no stage does one target give support and another reject. MHC-I and HLA-E
are weaker: JAK1 alone has no resolved MHC-I response, and only STAT1 alone has a resolved HLA-E response, so those two
pooled verdicts hold for the pool only.

**THP-1, PD-L1: pooling holds, by target and by batch.** The head output is rejected for all five targets and all three
batches, and none of TEDDY's three stages is rejected. Two things the pooled run could not show:

- The unreadable share at the head output differs widely by target: STAT1 0.28; IFNGR1, IFNGR2 and JAK2 0.56 to 0.60;
  IRF1 1.07.
- The batches differ twofold in the size of the response (−0.47 in batch 1, −0.96 and −0.88 in batches 3 and 4) but not
  in the share lost (0.44 to 0.54).

For PD-L2 the head output is rejected for IFNGR2, JAK2 and IRF1 and unresolved for IFNGR1 and STAT1; IRF1 is rejected
already at TEDDY's embedding Z (0.44), the first rejection at one of TEDDY's own stages on this data set (one target, one
readout).

**Why the targets differ** (`results/why_targets.json`, `results/why_gap_over_spread.json`). In TEDDY's embedding the
shifts of IFNGR1, IFNGR2, JAK2 and STAT1 point almost the same way (cosine 0.94 to 1.00); IRF1 has 0.26 to 0.30 with
them and a shift a third as long (1.5 against 4.4 to 6.6). In the RNA, IRF1 changes 293 genes, the receptors and JAK2 742
to 950, STAT1 1,535. The head is not insensitive to IRF1's direction: moved by the same length, its latent reacts 3.1 times
as much as to a random direction (2.3 to 2.9 for the other targets). The difference arises at the last step: from Z to the
head output the gap between knockout and control cells, relative to the spread between cells, grows for the other nine
targets (1.3 to 2.1 times) and shrinks only for IRF1 (0.14 at Z, 0.20 at the latent, 0.09 at the output). A reading
compatible with this, not tested directly: the head carries the large shared change of the whole IFN-γ programme being
switched off, and flattens the small, differently directed change left by a knockout that sits next to the readout. Not
explained: batch 1 has a smaller mean shift and no larger loss. On THP-1 the head's own PD-L1 output moves the wrong way
for all five targets, so "kept" here means that the 134 outputs together still carry information about the knockout, not
that PD-L1 is predicted correctly.

**Ledger** (`results/ledger.json`). The verdicts of 15 runs (both sealed pooled runs, 10 per-target, 3 per-batch), 155 in
all, contain no pair of opposite verdicts at the same readout and stage.

**No leakage** (`results/no_leakage.json`). All 25 real knockout events of the two TEDDY data sets pass on the three
routes (how cells were chosen, shared measurement, the readout's own gene). The two KS-2 declarations do not fully pass:
readout and model input share the RNA counts (declared), and the origin of the labels is not declared. Of the 252 in-silico
events of the forward run, 56 act on a readout's own gene and are flagged.

## Where this departs from the sealed procedure

- The statistics ran on features exported with the real adapter (`scripts/export_features.py`) and served by
  `scripts/precomputed_adapter.py`; the data and the model never left the machine that holds them.
- With one target as the treated class the resampling unit is the guide (3 guides per target on melanoma, 4 on THP-1); with
  one batch it is the target (5). Rule `K1` warns that intervals from so few clusters are approximate.
- The THP-1 runs do not compute the drop-order probe (`SET`), which is off the main chain. The per-batch runs keep the gene
  filter of the pooled selection.
- Batch 4: with the declared seed the fold assignment puts all five targets into one fold. The runner's fallback (package
  README, *Folds when treated clusters are few*) takes the next seed for the folds; the record says so (`fold_seed`). The
  other records were made before the fallback existed and carry no `fold_seed`; their split was the declared one, which
  the fallback leaves unchanged.

## Reproducing

```
cd scripts
python make_declarations.py                      # writes <run>/declaration.yaml from the sealed declarations
# on the machine with the data and the model, once per data set (see the docstring for --head):
python export_features.py --paths ../../../paths/papalexi2021.yaml --declaration ../papalexi2021_pooled_full/declaration.yaml --out $ANM_FEATURES/papalexi2021
python run_one.py papalexi2021_target_IRF1       # one run, about an hour on one core
python analyze_frangieh.py; python analyze_papalexi.py; python ledger.py; python no_leakage_table.py
python join_closure_admission.py; python why_targets.py; python why_gap_over_spread.py
```

`*_pooled_full/` holds only a declaration: the pooled selection with every node, used for the feature export.
`ANM_FEATURES` is the folder holding `frangieh2021/`, `papalexi2021/` and
`head_weights.npz`; the features are not kept in git.
