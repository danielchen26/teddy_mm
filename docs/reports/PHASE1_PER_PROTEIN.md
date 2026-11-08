# Phase 1: per-protein test results and known limitations

Run: `outputs/cite_phase1_official` (frozen TEDDY-G 70M, official preprocessing, gene-mean pooled 512-D embedding,
L2-normalized → negative-binomial head; MLP and latent flow-matching heads; seed 0; size factor = train median 0.92663).
Test set: the 16,750 site4 cells of NeurIPS 2021 BMMC CITE-seq; 134 proteins. Per-protein values:
[`phase1_per_protein.csv`](phase1_per_protein.csv) (copied from the run's `test_per_protein.json`, sorted by MLP Pearson).

## Summary over the 134 proteins

| | MLP head | flow-matching head |
|---|---|---|
| mean Pearson (the headline number) | **0.6027** | 0.5821 |
| median Pearson | 0.6752 | 0.6471 |
| mean R² | 0.0261 | −0.0180 |
| median R² | 0.2938 | 0.2499 |

- **Mean and median come from the same run.** Mean Pearson 0.603 is pulled below the median 0.675 by a long lower tail:
  quartiles of the MLP Pearson are 0.396 / 0.675 / 0.808, range 0.039–0.932; 24 proteins are below 0.3 and 38 are at
  or above 0.8. The flow-matching head beats the MLP on 22 of 134 proteins.
- **Best predicted** (MLP Pearson): CD36 0.932, CD72 0.930, CD19 0.928, CD22 0.911, CD16 0.904 — lineage markers with
  strong mRNA–protein coupling.
- **Worst predicted:** CD42b 0.039, CD223 (LAG-3) 0.049, CD163 0.078, TCRVa7.2 0.086, TCRVd2 0.089, CD169 0.105,
  CD270 0.115, CD23 0.124, CD134 (OX40) 0.135, CD137 (4-1BB) 0.140, CD83 0.152, LOX-1 0.159 — activation markers,
  rare-population markers and TCR variable chains, where the protein level is set after transcription or the cells are
  too rare in training.
- **Proteins used in the ANM diagnosis:** HLA-A-B-C 0.574, CD274 (PD-L1) 0.396, HLA-E 0.322, CD119 (IFNGR1) 0.250.
  The interferon-response proteins sit in the lower half of the panel even in bone marrow.
- **Why mean R² is so low next to Pearson:** 43 of 134 proteins have negative R², a few of them strongly (CD83 −3.45,
  CD141 −3.15, CD23 −2.52, CD154 −1.92, CX3CR1 −1.84). These proteins are ranked reasonably but predicted on the wrong
  scale, and they drag the mean to 0.026 while the median is 0.294. A likely contributor (not yet tested): the head's
  negative-binomial likelihood is fit to the ADT matrix as stored in the h5ad, which is already normalized upstream
  (non-integer values), so the count likelihood is mis-specified for these targets.

## Known limitations of the phase-1 number

These are disclosed elsewhere in the repository (README "Disclosures" and "Honest weak spots",
[`V3_REGISTERED_RESULTS`](V3_REGISTERED_RESULTS.md)); collected here so the 0.60 is read correctly.

1. **Not donor-independent.** The split is by site, not by donor: donor 15078 is in training at sites 1–3 and is 32.6%
   of the site4 test cells. A donor-held-out evaluation (v3 uses donors 13272 and 19593, which appear in no other split)
   is the cleaner estimate.
2. **No raw-RNA baseline.** No head was trained on normalized RNA directly (for example log1p CP10k → PCA → the same
   head), so the 0.60 measures "frozen TEDDY + this head", not what TEDDY adds over the RNA itself.
3. **Likelihood on normalized values.** See the R² note above; an NB head on raw ADT counts, or a Gaussian/CLR-space
   loss, has not been compared.
4. **One seed, one site.** No variance over seeds or over held-out sites is reported for phase 1.
5. **Selection on the test set downstream.** The 9 marker proteins used by the ANM decision-layer experiments were
   chosen by test-set Pearson; numbers that depend on that panel are optimistic.
6. **Registration order.** For v3 E1 and E5 the registration order rests on local commit times (pushed after
   evaluation).
7. **Distribution shift in the ANM tests.** The head was trained on bone-marrow CITE-seq; the ANM diagnosis applies it
   unchanged to melanoma (Frangieh 2021) and THP-1 (Papalexi 2021) cells.

## What would tighten it

A donor-held-out split; a raw-RNA baseline with the same head and budget; NB on raw ADT counts (or a loss matched to the
stored normalization); three seeds; panel selection on validation, not test.
