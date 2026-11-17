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
