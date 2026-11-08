# Finished runs

Each folder holds the JSON record (`anm_record.json`: declaration, every interval, provenance) and the Markdown report
(`anm_report.md`) written by `python -m anm_diag run`. Results and their reading: [`docs/ANM_DIAGNOSIS.md`](../../docs/ANM_DIAGNOSIS.md).

| folder | data | status | what stands |
|---|---|---|---|
| [`ks2_mitoribo/`](ks2_mitoribo/) | Replogle 2022 K562 CRISPRi, mitoribosome targets (6 treated targets, 1,956 cells) | exploratory | all stages (RAW, TOKENS, Z); read the reread report, see the folder README |
| [`ks2_ribo/`](ks2_ribo/) | same screen, ribosome targets (12 treated targets, 2,196 cells) | exploratory | all stages |
| [`frangieh2021/`](frangieh2021/) | Perturb-CITE-seq, melanoma, IFN-γ (5 treated targets, 5,171 cells) | sealed, run once | **RAW, TOKENS, Z only; head stages void** (see the folder README) |
| [`papalexi2021/`](papalexi2021/) | ECCITE-seq, THP-1, IFN-γ (5 treated targets, 7,478 cells) | sealed, run once | **RAW, TOKENS, Z only; head stages void** |

The head-stage rerun (same declarations, fixed head input, resealed) writes to `frangieh2021_headfix/` and
`papalexi2021_headfix/`; until those folders exist there is no valid head-stage result.
