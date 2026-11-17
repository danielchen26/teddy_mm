# Round 1, gate 5: the one test on a study nobody has looked at — record

Gate 5 of [`DESIGN.md`](DESIGN.md). This file is written in parts, each by a commit made before the step it governs.
Nothing in a part is changed after its commit.

## Part 1. The decision of gate 4 and the frozen recipe

Committed 17 November 2026, before any fit of gate 5 and before the final study was looked at in any way.

**Gate 4.** The go rule holds ([`results/development.json`](results/development.json), `go_to_gate_5`): `coverage` is
admitted in every seed; its gain with the donor's batch not in the fit is 0.053 [0.044, 0.062], above 0; the guard holds
for every head.

**Approval.** The owner of the repository approved, on 17 November 2026, running gate 5 and opening the final study.

**The recipe.** Frozen as [`DESIGN.md`](DESIGN.md) has it ("The change", "Arms"). The one revision that development
allowed is not used. Nothing in `scripts/train_arm.py`, `scripts/plan.py` or `scripts/run_all.py` is changed.

**The heads of the test.**

| arm | heads | fitted |
|---|---|---|
| `as_is`, `amount`, `function` | `s0`, `s1`, `s2` of each | in development; their weights are named by hash in [`results/gate5_heads.json`](results/gate5_heads.json) |
| `coverage` | `s0_coverage_all`, `s1_coverage_all`, `s2_coverage_all`: 47,405 bone-marrow cells and 20,000 PBMC cells of all 8 donors (2,500 each) | now, `scripts/run_all.py --final` |
| `labels` | `s0_labels_all`, `s1_labels_all`, `s2_labels_all`: the same cells, proteins shuffled within the donor | now, the same command |

**Admission of the two new arms** (gate 3, as in development): the validation score does not fall at the end of
training (last 5 epochs not more than 0.005 below the 5 before); `coverage` scores at least 0.035 above `as_is` of its
seed on its own PBMC fit cells; no cell of the final study is in any file a fit reads; the validation donor's score is
no more than 0.02 below `as_is` of the seed. An arm that fails was not tested.

**What is read on the final study F**, fixed here:

- *Units.* The donors of F chosen by the rule of gate 0b, each with its 1,000 cells. They are listed in part 2.
- *Proteins.* Those F measures among the 112 the PBMC study shares with the head, by name. Listed in part 2.
- *Quantity.* For a donor, *u*: 1 minus the mean over those proteins of r² between the head's output and the
  measurement among the donor's cells (`scope.arms`, as `scripts/read.py`). The gain of an arm is *u* of `as_is` of the
  same seed minus *u* of the arm, on the same cells.
- *Over donors.* The mean over the donors of F, with an interval from 1,000 resamples of donors, seeds resampled inside
  each arm (`scope.arms_family`), seed 0, level 0.05.
- *The test.* `loop.test` reads hypothesis `donors_hao_coverage` as it stands in
  `runs/scope_20261116/results/ledger.jsonl`: own arm `coverage`, bound 0.0342; controls `function` and `amount`, each
  below the same bound. It is read once over the three seeds together and once for each seed (the arm of a seed
  against `as_is` of that seed, donors resampled).
- *The outcome of the round* is the outcome over the three seeds together, and it is stated as such only if each seed
  alone reads the same. If the seeds disagree with it or among themselves, the round's outcome is "unresolved".
- *The second hypothesis*, `coverage_not_levels` (own arm `coverage`, control `labels`, the same bound), is written to
  the ledger in part 3 with forecasts from development, and read the same way.
- *The guard* (gate 6): on bone marrow, site 4, the mean score of every head of the test is no more than 0.02 below
  `as_is` of its seed.
- *Also read, and not part of the outcome:* the gain of `coverage` among healthy donors, among patients, and patients
  minus healthy (donors resampled inside each group); the head's ordinary score; and, from the grid of step 6, the
  share of F's own ceiling that was delivered.

**Order.** (1) The new fits. (2) The final study is checked against the criteria of gate 0b from its cell and feature
tables only, its donors and cells are chosen, the pack is written and the cells are embedded (part 2). (3) The second
hypothesis and the declaration with the hashes of every input (part 3). (4) The heads are run on F. (5) The grid on F.

## Part 1b. The rule of gate 0b as it is applied, and the reading script

Committed 17 November 2026, while the file of the first candidate (Stephenson et al. 2021) was being downloaded and
before it was opened.

- [`scripts/final_study.py`](scripts/final_study.py) states the rule of gate 0b in the detail the design left open:
  which sample of a donor, when a donor is eligible, the order in which donors are taken over centres, how cells are
  drawn, and how a protein of the head is matched to a feature of the study. Its first stage reads the cell table and
  the feature table of the file and nothing else.
- [`scripts/read_final.py`](scripts/read_final.py) is the script that reads F, with its tests
  (`tests/test_round1_final.py`). It refuses to run unless the declaration of part 3 is there and every input matches
  it by hash.
- The cells are embedded with `scripts/03_embed_rna.py` and the settings of the PBMC study's manifest (official
  preprocessing, sequence length 2,048, length buckets, batch 32, fp16, gene-mean pooling).

## Part 2. The final study: criteria, donors, cells and proteins

Committed 17 November 2026, before any count of the file was read. What follows was read from the cell table and the
feature table of the file ([`results/final_selection.json`](results/final_selection.json), written by
`scripts/final_study.py select`).

**The study.** Stephenson et al. 2021, Nat Med 27:904, COVID-19 PBMC CITE-seq (ArrayExpress E-MTAB-10026; haniffa21.processed.h5ad): the first candidate of gate 0b. 647,366 cells, 24,929 features. The file was
downloaded on 17 November 2026 from the study's public data portal (7,187,322,881 bytes; its md5 is the one the server
states, add2501947c585ea6b6ef7429e4bd9f2).

**Criteria.**

| criterion | holds | |
|---|---|---|
| open | yes | public file of the study's data portal; no access request |
| raw RNA and protein counts per cell | yes | a layer of raw counts over genes and protein features; that it holds integers is checked when the pack is written |
| a donor and a disease label per cell | yes | `patient_id`, `Status` |
| at least 20 donors, 8 healthy, 8 patients | yes | 23 healthy donors, 90 patients (COVID-19); 17 donors of other groups are not drawn |
| at least 80 of the 112 proteins, by name | yes, **with aliases** | 95 of the 112: 49 by exact name and 46 under another name of the same target (below) |
| never used in this repository | yes | |

**The proteins.** The study names many antibodies by the gene symbol of their target; the head's panel uses CD numbers. By exact name 49 of the 112
match, which is below the 80 the criterion asks. With the names of the same target in the standard nomenclature, written from the study's 192 feature
names and nothing else (`ALIASES` in `scripts/final_study.py`), 95 match. The criterion is taken as met on that basis, and the reading on the
49 exact-name proteins alone is added to `scripts/read_final.py` as a second reading that does not change the outcome. Clones were not compared.

Matched under another name: CD155 = PVR, CD112 = NECTIN2, CD154 = CD40LG, CD11c = ITGAX, CD194 = CCR4, CD279 = PD1, CD335 = NCR1, CD31 = PECAM1, Podoplanin = podoplanin, CD146 = MCAM, CD195 = CCR5, CD196 = CCR6, CD185 = CXCR5, CD103 = ITGAE, CD152 = CTLA4, CD223 = LAG3, CD107a = LAMP1, CD95 = FAS, CD1c = CD1C, CD11b = ITGAM, CD314 = KLRK1, CD35 = CR1, CD57 = B3GAT1, CD272 = BTLA, CD278 = ICOS, CD39 = ENTPD1, CD169 = SIGLEC1, integrinB7 = ITGB7, CD268 = BAFFR, CD54 = ICAM1, CD62P = SELP, CD122 = IL2RB, CD41 = ITGA2B, CD137 = TNFRSF9, CD124 = IL4R, CD29 = ITGB1, CD49b = ITGA2, CD18 = ITGB2, CD127 = IL7R, CD26 = DPP4, CD49d = ITGA4, CD73 = NT5E, TCRVa7.2 = TCR_Va7.2, LOX-1 = OLR1, CD158e1 = KIR3DL1, CD319 = SLAMF7.

Of the 112 not measured: CD270, CD48, CD105, CD134, CD42b, CD119, CD192, CD13, CD115, CD63, CD172a, CD72, CD93, CD49a, CD9, TCRVd2, CD142. `TCRVd2` is left out because the study's `TCR_Vg2` cannot be told from its name to be the same target.

**Donors and cells.** Eligible (one status, one centre, at least 1,000 cells in the donor's first sample): 22 healthy donors and 86 patients.
Chosen by the rule: 20 healthy donors (Cambridge 10, Ncl 10; the third centre has no healthy donor) and 20 patients (Sanger 6, Cambridge 7, Ncl 7); 1,000 cells each, 40,000 cells. The other 73 healthy donors and patients stay unopened as the reserve.

| group | donors (the study's ids) |
|---|---|
| healthy | CV0902, CV0904, CV0911, CV0915, CV0917, CV0926, CV0929, CV0934, CV0939, CV0940, MH8919177, MH8919179, MH8919226, MH8919227, MH8919282, MH8919283, MH8919332, MH8919333, newcastle65, newcastle74 |
| patients | AP10, AP2, AP3, AP5, AP6, AP9, CV0037, CV0100, CV0134, CV0137, CV0178, CV0231, CV0284, MH8919330, MH9143273, MH9143276, MH9143321, MH9143420, MH9143421, MH9179821 |

These 40 donors are the units under `roles.final`.

## Part 3. The pack, the embedding, the second hypothesis and the declaration

Committed 17 November 2026, after the pack and the embedding of the final study were written and before any protein
value of it was scored, summarised or looked at. `scripts/read_final.py` is run after this commit is pushed.

**The pack.** `scripts/final_study.py pack` read the raw counts of the 40,000 chosen rows and wrote them in the layout of
the repository's other packs: 12,052 gene columns of the head's universe, of which 11,747 are in the study (the others
are zero), and 134 protein columns, of which 114 are measured in the study; 95 of those are among the 112 proteins the
head's two panels share. The counts are integers. A protein value is the per-cell CLR over the 192 protein features of the
study, as in the other packs. The script printed the numbers of cells, genes and proteins and no value. One check was
made on the protein values afterwards, and nothing else: none is missing in a measured column and none is present in a
column that is not measured. The names of its genes and proteins are those of the head's packs, in the same order.

**The embedding.** `scripts/03_embed_rna.py` with the settings of part 1b; its configuration hash is `1ac673a8a4ad5204`,
the one in the manifest of the development study's embedding. It ran on the owner's machine (Apple GPU), in 955 seconds.

**Departures from what was written.**

1. Part 2 says it was committed "before any count of the file was read". The pack was started in the same minutes as that
   commit, so the script was reading counts while the commit was made. Nothing of them had been printed or looked at.
2. The embedding ran under Python 3.11.15 with numpy 2.4.6. The manifest of the development study's embedding names
   Python 3.12.9 with numpy 2.5.3. torch is 2.14.0 in both. To measure what that
   changes, the first 1,000 cells of the development study were embedded again in the environment used here: the
   largest difference from the stored embedding is 1e-6 in a coordinate (values are near 1), at most 7e-8 of a
   cell's length. The two environments give the same embedding to that precision.
3. The criterion "at least 80 of the 112 proteins, by name" is met only with the alias list of part 2 (49 by exact name).
   The reading on the 49 alone is in `scripts/read_final.py` and does not change the outcome.
4. [`DESIGN.md`](DESIGN.md) enters the second hypothesis in the ledger before F is opened. Part 1 fixed its own arm, its
   control and its bound before the file was looked at; the ledger entry itself is written here, after the pack and the
   embedding. Its forecasts are numbers of development, and no value of the final study enters it.

**The second hypothesis.** `coverage_not_levels`, entered in [`results/ledger_final.jsonl`](results/ledger_final.jsonl)
by `scripts/declare_final.py`: what the added PBMC cells teach the head is the relation between RNA and protein, not the
protein levels of their study. Its rule: on the donors of the final study the fall of *u* with `coverage` has its interval
at or above 0.0342, and with `labels` (the same cells, the proteins shuffled among the cells of each donor) below it. Its
forecasts are the intervals read in development: `coverage` 0.044 to 0.062, `labels` −0.353 to −0.328. The first
hypothesis, `donors_hao_coverage`, is taken unchanged from `runs/scope_20261116/results/ledger.jsonl`.

**The declaration.** [`results/gate5_declaration.json`](results/gate5_declaration.json), digest `79e1188fd91938664e2ee9861eb68d2971125e51627aca32cf803ede1f45950a`: the two
hypotheses, the five arms and their fifteen heads, the 40 donors, the 95 proteins (and the 49 by exact name), how the
outcome is stated (part 1), and the sha256 of the pack, of the embedding, of every head and of the five scripts the
reading runs. `scripts/read_final.py` stops if any of them differs.

**The six heads fitted for the test** (`coverage_all` and `labels_all`, three seeds) were fitted from the commit of
part 1 and finished on 17 November 2026. Their admission (gate 3) and the guard (gate 6) are read by
`scripts/read_final.py` in the same run as the test, with the code of development (`scripts/read.py`).

**If the reading script stops with an error.** It prints nothing and writes nothing about the final study before it
writes `results/final.json`. If it stops before that for a reason of code, the error is recorded in a part 4, the fix is
committed with a new declaration, and the script is run again. Once `results/final.json` exists it is not run again.

**After the test.** `scripts/grid_final.py` reads the grid on the final study for `as_is` and `coverage`: each donor
against the other donors of the study, and the disease as an event. Seed 0 only, for its cost. The design had not fixed
the number of seeds. It is read after `results/final.json` exists, is for a next round and changes no outcome. Its RNA view
is the projection on the 512 components of `runs/scope_20261116`, refitted on bone marrow and the development study and
checked against the stored projection before it is used.

**What writes the README.** `scripts/write_final.py` (the section about the test and its figure) is committed here,
before the reading. It reads `results/final.json` and decides nothing.

## Part 4. After the reading

Written 17 November 2026, after `results/final.json` existed. Nothing in parts 1 to 3 is changed.

**The reading.** `scripts/read_final.py` was run once, from the commit of part 3 (`4dbca69`), with the inputs of the
declaration. It ran without error and wrote [`results/final.json`](results/final.json) and, to
[`results/ledger_final.jsonl`](results/ledger_final.jsonl), six intervention events (the heads fitted for the test) and
eight test events (two hypotheses; the three seeds together and each seed). It has not been run again.

**The outcome, as part 1 states it.** Both hypotheses: *unresolved*. Over the three seeds together the rule is not
testable, because `s1_coverage_all` is not admitted (gate 3: its validation score over the last 5 epochs is 0.015 below
the 5 before; the limit is 0.005). Seeds 0 and 2, taken alone, read *refuted*: the gain of `coverage` is −0.007 and
−0.009, with intervals below the bound of 0.0342 and below 0. Over the three seeds the gain is −0.013 [−0.022, −0.004].
The guard of part 1 is stated over every head of the test and does not hold: `s2_labels_all`, a control, is 0.022 below
`as_is` of its seed on bone marrow, site 4. The three `coverage` heads are within 0.02. The numbers and their reading
are in [`README.md`](README.md).

**What was not done after the reading.** No head was fitted again and no further seed was fitted. The 73 donors of the
reserve were not opened: the design names them as the step after *unresolved*, but the outcome is unresolved for the
admission of a head, which more donors do not change, and opening them is the owner's decision. No arm, protein or donor
was added or dropped.

**Checks made after the reading (post hoc).**

- *The pack.* `scripts/check_final_pack.py` → [`results/final_pack_check.json`](results/final_pack_check.json): within a
  donor, the correlation between a protein and the RNA of its own gene, for 30 proteins of the final study and, with the
  same code, for the same proteins of the development study. The gene ids were written by hand; the development study
  is the reference for each pair. The pack is sound (CD3 0.63, CD8 0.73, CD20 0.74); over the 26 pairs in both studies
  the mean is 0.38 against 0.54.
- *The matched names.* The reading on the 49 proteins with the same name in both panels was declared in part 2 and is in
  `results/final.json`; it is the same.

**Departures after part 3.**

1. `scripts/grid_final.py` stopped twice at its own check before it ran. The constant it used for the seed of the RNA
   decomposition was 0; the seed of `runs/scope_20261116/scripts/build.py` is 20261116. With the seed corrected the
   decomposition returns the stored view exactly, and the projection of a cell on its components equals the stored view of
   that cell. Between the two attempts the script was also changed to take the RNA view of the training cells from the
   same projection as the cells of the final study; with the right seed this is the stored view. No number had been
   written by the script before it ran to the end.
2. `scripts/write_final.py` was changed after the reading: the axis ranges of the figure are taken from the data, and
   the text, which had been written for a gain, was rewritten for what was read (the admission and the guard are stated
   head by head; the row of the design's table is given both as the rule states it and as the admitted seeds read). No
   number of `results/final.json` is changed by it.
3. While `s2_coverage_all` was being fitted, a trial of the grid script on other data was started on the same machine
   and stopped after two minutes. One epoch of that fit took 135 seconds instead of 19. A fit reads the clock only between
   epochs, for a time budget that was not reached, and that head is admitted. `s1_coverage_all`, the head that is not admitted, had finished a quarter of
   an hour before the trial started.
4. The grid on the final study is of seed 0 only (part 3).
5. The grid on the final study was still running when this part was committed. After the two stops of item 1 it was
   stopped a third time, during its first run and before it had written anything, by the time limit of the job runner,
   and started again. Its reading goes into the README in a later commit, as read after the test.

## Part 5. The grid on the final study (step 6), as it was read

Written 17 November 2026, after the test and after part 4. It is for a next round and changes no outcome.

**What was read.** For each of the 40 donors, *u* of the variation among its cells at the head's own output and under
readouts that were not fitted on the donor: one fitted on the head's training cells, one fitted on the other 39 donors
of the study, and one with the donor's own cells, cross-fitted. Two views of a cell: 512 RNA components, and TEDDY's
state. Two heads of seed 0: `as_is` and `coverage` (`s0_coverage_all`). The donor-by-donor values are in
[`results/grid_final/s0_donors_variation.json`](results/grid_final/s0_donors_variation.json), the means and contrasts in
[`results/what_is_left_final.json`](results/what_is_left_final.json), both written by `scripts/grid_final.py`
(`donors_variation`). The readouts and the resamples are those of a `scope` run of kind `donor` (`scope.Study`,
`scope._u_variation`); the readouts fitted on cells of the final study do not depend on the head and were fitted once.

**Departures from parts 1 and 3 and from the design.**

1. *Not the scope runs that part 3 names.* A scope run with the donors as events also reads the `shift` target, the
   level of a donor against the other donors. That needs a readout for every pair of donors: 780 pairs at each of four
   stages, for each head. Two such runs were stopped by hand after an hour each, with nothing written; a run would have
   taken about ten hours. Only the `variation` target was read, which is the target of the round.
2. *The stages inside the head were not read* on the final study (the head's latent and its output as views for a
   readout). The two views that were read do not depend on the head.
3. *The disease as an event was not read.* The run the script would have made resamples the patients' cells one by one
   and not by donor, which would understate its intervals. The comparison of healthy donors and patients that was
   declared for the test is in `results/final.json`.
4. *Typed events.* Because no scope run was made, the grid on the final study wrote no source or verdict events, and
   `events.check_contrast` was not applied to it. The 40 donors are named under `roles.final` only in the declaration.
5. *How it was run.* Six starts. Two stopped at the script's own check (part 4, item 1). One was stopped by the time
   limit of the job runner. One was killed for memory when the arrays of other runs were prepared beside it. One, as
   two processes side by side, was stopped by hand (item 1 above). The sixth is the reading kept. None of the first
   five wrote a number.
6. *One seed*, as part 3 says: the intervals of this part resample donors only.

**What it read**, in one line: on the final study a readout fitted on the study's own donors leaves 0.006 [0.002, 0.011]
less unexplained than `as_is`, against 0.073 on the development study. The bound of the rule, 0.034, could not have been
met there. The reading is in [`README.md`](README.md); the bound had been set without it, because the final study was
not to be opened before the test.
