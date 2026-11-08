# Status: head stages void, TEDDY stages valid

**HEAD_LAT, HEAD_OUT and every step into or out of them (Z → HEAD_LAT, HEAD_LAT → HEAD_OUT) in this report and record
are void.** The adapter fed TEDDY's raw gene-mean embedding to the phase-1 head, but teddy_mm trains and evaluates that
head on embeddings divided by their L2 norm (`scripts/04_train.py` line 101, `scripts/05_eval.py` line 55). Raw
embeddings have norm ≈ 26, so the head was tested on inputs about 26 times larger than any it was trained on; the
head-stage verdicts and root causes (`loss_at` / `relative_loss_at` at Z → HEAD_LAT) describe that wrong operator, not
the head.

**Valid:** RAW, TOKENS, Z, the SET and mRNA probes, the positive controls, the response Δ and the forward AUCs of those
stages. None of them uses the head.

The adapter now applies the same normalization (`PerturbCiteAdapter.head_input`; checked element-wise against
teddy_mm's torch code to within 1e-5, covered by `tests/test_head_operator.py`). The same declaration (design
unchanged) is resealed with the fixed code and run once into `../<dataset>_headfix/`. The error was found after the first
results were seen, so that rerun is recorded as a formal test with a deviation.

This record was sealed against the code before the fix (its `code_hash` no longer matches the package), so it cannot be
rerun from this checkout; it is kept as the record of what was run.
