# Finished runs

Each folder holds the JSON record (`anm_record.json`: declaration, every interval, provenance) and the Markdown report
(`anm_report.md`) written by `python -m anm_diag run`. Results and their reading: [`docs/ANM_DIAGNOSIS.md`](../../docs/ANM_DIAGNOSIS.md).

| folder | data | status | what stands |
|---|---|---|---|
| [`ks2_mitoribo/`](ks2_mitoribo/) | Replogle 2022 K562 CRISPRi, mitoribosome targets (6 treated targets, 1,956 cells) | exploratory | all stages (RAW, TOKENS, Z); read the reread report, see the folder README |
| [`ks2_ribo/`](ks2_ribo/) | same screen, ribosome targets (12 treated targets, 2,196 cells) | exploratory | all stages |
| [`frangieh2021_headfix/`](frangieh2021_headfix/) | Perturb-CITE-seq, melanoma, IFN-γ (5 treated targets, 5,171 cells) | sealed, run once (rerun after the head-input fix) | **all stages**, including the head |
| [`papalexi2021_headfix/`](papalexi2021_headfix/) | ECCITE-seq, THP-1, IFN-γ (5 treated targets, 7,478 cells) | sealed, run once (rerun after the head-input fix) | **all stages**, including the head |
| [`frangieh2021/`](frangieh2021/), [`papalexi2021/`](papalexi2021/) | first sealed run | superseded for the head | RAW, TOKENS, Z only; head stages void |
| [`frangieh2021_headsplit/`](frangieh2021_headsplit/), [`papalexi2021_headsplit/`](papalexi2021_headsplit/) | the head split into its operations (two-sided rule) | sealed, exploratory | all stages; see folder README for code version |
| [`forward_typed_events/`](forward_typed_events/) | typed source events through TEDDY and the head | sealed, exploratory | report and record; see folder README for code version |
| [`split_20261111/`](split_20261111/) | the pooled source of the two head-fix runs separated: melanoma by target (5 runs), THP-1 by target (5) and by batch (3), plus quick pooled reruns | exploratory (made after the pooled results were seen) | per-target and per-batch verdicts, the verdict ledger, the no-leakage table, the comparison of targets; see the folder README |
| [`teddy_mm_head/`](teddy_mm_head/) | singular values of the head's weight matrices (no cell data) | descriptive | — |
