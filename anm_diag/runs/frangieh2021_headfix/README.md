# Frangieh 2021: sealed rerun after the head-input fix

The record to cite for Frangieh 2021. Declaration: [`../../declarations/teddy_mm_frangieh2021_headfix.sealed.yaml`](../../declarations/teddy_mm_frangieh2021_headfix.sealed.yaml)
(same design as the first sealed run; resealed with the fixed code before any rerun result; code hash `cb8fca78dc47…`,
which is the package in this repository). Run once; the seal was verified at run time.

Provenance notes:

- The rerun ran on a second Mac (Python 3.11.13 vs 3.11.15, identical package versions). RAW, TOKENS and Z match the
  first run to three decimals.
- The seal records, for `teddy_ckpt`, the sha256 of an empty file (`e3b0c442…`): on that machine the checkpoint path
  read as zero bytes. TEDDY was not run again; Z came from the states computed in the first run, whose seal certifies
  the real checkpoint (`ba0a98f9…`). The results are unaffected, but this seal does not certify the TEDDY weights.
- The fix came after the first run's results were seen, so a strict reader may treat the head-stage conclusions as
  exploratory. A post hoc two-sided reread (exploratory) is summarised in [`docs/ANM_DIAGNOSIS.md`](../../../docs/ANM_DIAGNOSIS.md#two-sided-reread-post-hoc-exploratory).
