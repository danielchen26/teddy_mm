# Status: head stages void here; see `../frangieh2021_headfix/`

**HEAD_LAT, HEAD_OUT and every step into or out of them in this report and record are void.** The adapter fed TEDDY's
raw gene-mean embedding to the phase-1 head, but teddy_mm trains and evaluates that head on embeddings divided by their
L2 norm (`scripts/04_train.py` line 101, `scripts/05_eval.py` line 55). Raw embeddings have norm ≈ 26, so the head was
tested on inputs about 26 times larger than any it was trained on.

**Valid:** RAW, TOKENS, Z, the SET and mRNA probes, the positive controls, the response Δ and the forward AUCs of those
stages. None of them uses the head; the rerun reproduces them to three decimals.

The same declaration, resealed with the fixed code, was run once into [`../frangieh2021_headfix/`](../frangieh2021_headfix/): that is
the record to cite for Frangieh 2021. This record was sealed against the code before the fix (its `code_hash` no longer
matches the package), so it cannot be rerun from this checkout; it is kept as the record of what was run.
