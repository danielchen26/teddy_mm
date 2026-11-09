# Step 4 (closing the loop): which change to the head, chosen by a rule written before the head-split result

Status: written 2026-11-09, before the head-split declarations were sealed or run. No head-split result exists.
Everything under "Fixed now" stays as written. The items under "Open until the head-split result" are settled after
the head-split result and before any result of this step, in a sealed step-4 declaration. A changed rule is a new
document, recorded as a deviation.

## Why this step (ANM paper, arXiv 2610.14549)
- Table 2, last row (p. 8). Observation: "A task-required contrast is absent despite preserved input records."
  Supported conclusion: "Test how the program or model uses the supplied records." Next experiment: "Compare
  changes to the program or model with controls; remeasure state and readout." Frangieh PD-L1 is this row. The
  sealed headfix run gives bounded support at RAW, TOKENS and Z, and rejection at HEAD_LAT and HEAD_OUT. The
  two-sided rereading gives the same.
- SI S1.3. Joint changes of state and readout establish association under intervention, not mediation, and
  unresolved responses do not localize effects. A state measurement or a known selective intervention separates
  altered propagation from altered readout (paragraph "State measurements and selective interventions").
- p. 34 ("Representation information and intervention benefit") and p. 35 (Section 11). Changing a system is a
  different question from testing a description of it, and a change to what a downstream component receives is a
  system intervention. Establishing a benefit needs independent tasks and comparators with the same information.
- p. 8. This step is designed after the sealed inverse runs, so its evidence is exploratory. It is still sealed
  before it runs.

## Fixed now
1. **What changes.** Only the teddy_mm head (its MLP and decoder) changes. TEDDY, the tokenizer, data, cells, folds,
   readouts, readout classes, margins and the verdict rule all stay as in the head-split declarations. The verdict
   rule is two-sided, and the step between adjacent stages is judged by the change in |u|.
2. **Which result decides.** Frangieh PD-L1 in the sealed head-split record (`runs/frangieh2021_headsplit/anm_record.json`)
   decides. It is the only readout with an absolute conclusion (`loss_at` Z → HEAD_LAT) in the sealed headfix run.
   The other four readouts are reported beside it and do not choose.
3. **The mapping.** The head-split steps fall into three groups:
   - group A: Z → Z_UNIT, the division by the length;
   - group B: Z_UNIT → L1 and L1 → H1, the first linear layer and its SiLU;
   - group C: H1 → L2, L2 → H2 and H2 → HEAD_LAT, the second and third layers.

   Read Frangieh PD-L1's root cause in the record:
   - (a) If it is `loss_at` and its whole span (last supported stage → first rejected stage) lies inside one
     group, choose that group's change.
   - (b) If it is `loss_at` and its span crosses groups, take the first step in the span whose interval for the
     change in |u| has a lower bound above the step margin (0.10). The intervals are the record's
     `steps["<a>-><b>"]["abs_change"]`. Choose that step's group. If no step in the span is resolved, choose C and report that the loss was not localized.
   - (c) If it is anything else, for example a loss outside Z → HEAD_LAT, `unresolved` or `reference_unresolved`,
     choose C and report that the head-split did not localize the loss.
4. **The three changes.** Each is retrained with teddy_mm's own recipe (`scripts/04_train.py`) apart from the
   change itself.
   - A: the head reads z / (|z| + 1e-6) plus log |z| as one extra input coordinate.
   - B: the same architecture, with a declared constraint on the first layer: spectral normalization of net.0, or
     a weight-decay setting.
   - C: the same architecture, trained on teddy_mm's training cells plus cells that carry IFN-γ-driven variation.
5. **Controls** (the "with controls" of Table 2):
   - Each change is compared with a control retrain that differs from it only by the absence of the change: same
     recipe, same data and the same seeds.
   - Each arm is trained with the same number of seeds, at least 3.
   - The current head is reported too. This separates the change from retraining itself.
6. **Disjoint cells.** No cell used for training enters any test, inverse or forward (rule R2: a readout's
   measurement may not feed a candidate state).
7. **Outcomes.**
   - Primary (Frangieh PD-L1): for each seed of each arm, run the head-split inverse declaration as a new sealed
     declaration naming that checkpoint. An arm succeeds if, in every seed, both of these hold:
     - HEAD_LAT has bounded support (two-sided);
     - the step chosen in item 3 adds no resolved loss: the lower bound of its |u|-change interval is at most 0.10.
       If C was chosen under (b) or (c), use the whole step Z → HEAD_LAT.

     The change is credited only if its arm succeeds and the control arm does not. If both succeed, the loss
     belonged to the trained instance, not to the recipe, and that is reported as such.
   - Secondary:
     - the other four readouts at the same steps (Papalexi run with the same arms);
     - the z_direction channel of the forward declaration (version 2), which runs through the head only, for every
       seed of each arm: the readout slope along the measured knockout direction, placed among the random
       directions.
   - Guard: the changed arm's mean per-protein Pearson correlation on teddy_mm's held-out site4 test cells must not
     fall more than 0.02 below the control arm's. Both are seed means. The test cells are set in
     `configs/cite_phase1.yaml` (`test_sites: [site4]`) and the correlation is computed by `scripts/05_eval.py`.
     The 0.02 is a chosen number, not derived. If the guard fails, the change trades one failure for another and is
     not credited.
8. **Reading.** Suppose the loss disappears with the change and not with the control. A recipe change aimed at that
   step then restores the readout, compared with retraining alone. That is evidence by intervention that the step
   matters for this readout. It does not prove that the step alone causes the loss, because retraining changes every
   weight. If the loss does not disappear, the explanation behind the chosen change is rejected for this head, which
   is also a finding. Either way the claim covers this head on these cells.

## Open until the head-split result (settled before any step-4 result, in the sealed step-4 declaration)
- The exact settings of the chosen change: for B, which constraint and its value; for C, which cells and their
  protein targets. Frangieh's antibody panel overlaps only part of the head's 134 outputs.
- The exact form of the control for C: added cells without IFN-γ stimulation, or the original training set alone.
- The number of seeds, if more than 3.
