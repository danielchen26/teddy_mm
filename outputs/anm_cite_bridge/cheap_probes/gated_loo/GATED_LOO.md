# Gated LOO ±ε / ±ε/2 (Mode B decision sensitivity)

**Not Mode A.** Nudges typed evidence on top LOO protein; counts G1–G4.

- eps=0.1, eps/2=0.05, kept=1500, runtime=0.6s
- existing LOO top1 flip rate (full): **0.2773**
- subsample LOO flip rate: **0.28733333333333333**
- any gated flip (1-G1): **0.0273**
- strong G3+G4: **0.0273**

## Gate counts
- **G1**: 1459 (0.9727) — no decision flip under ±ε or ±ε/2 on top LOO protein
- **G2**: 0 (0.0000) — flips under ±ε but not ±ε/2 (scale-gated)
- **G3**: 0 (0.0000) — flips under ±ε/2 (high sensitivity)
- **G4**: 41 (0.0273) — sign-asymmetric flip at ε and/or ε/2

Do not claim Pearson > ~0.61. Not clinical.
