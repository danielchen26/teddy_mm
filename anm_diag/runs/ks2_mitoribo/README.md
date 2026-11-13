# Status: exploratory; read `anm_report_reread.md`

The design was fixed after looking at these development cells, so every result here is exploratory.

- `anm_report.md` / `anm_record.json`: first pass. Each stage was judged by the readout class with the best point
  estimate. For the MT readout the positive control (READ_MT, the readout itself) was then judged by gradient boosting,
  whose interval reached +0.103 > 0.10, so the root cause was `control_failed` (no verdict on the chain), although ridge
  recovered it exactly (u = 0.00 [−0.00, 0.00]).
- `anm_report_reread.md` / `anm_record_reread.json`: `python -m anm_diag.reanalyze` with the rule the package now uses
  everywhere (a stage is supported if any declared class supports it, rejected only if every class rejects it). No fit
  was redone. MT: RAW −0.01 support, TOKENS +0.01 support, Z +0.47 [0.29, 0.64] reject → `loss_at` TOKENS → Z.
  PSM: RAW +0.19 [0.08, 0.29] unresolved → `reference_unresolved`.

Because the rule was changed after the first pass was seen, the reread is exploratory too. A formal claim needs the
declaration sealed and run once on new mitoribosome targets.
