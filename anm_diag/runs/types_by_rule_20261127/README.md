# Types by rule: the system proposes and tests its own types

27 November 2026. Run 48. Exploratory: no unit never read is opened, nothing is trained, and version 6 is not
changed. It is the candidate for one rule of a later version.

## Why

A *type* is a set of events whose records teach a readout for each other. Until now the system had one, the
interferon pathway (five genes), and the author had declared it by hand. Edges grow by rule since run 43 (flag,
propose, test on rows that did not propose them). Types did not. That left the author's hand in the loop.

## The rule

| step | what happens | code |
|---|---|---|
| flag | The law of a library, told of no type, flags an event that moves the state where a readout of control cells misreads what the event does. | `Law.flags` (unchanged) |
| propose | A flagged event seeds a type. Every event that moves the state the same way joins it (cosine of the two shifts at least 0.5): what a readout misreads along a direction, it misreads for every event that moves along it. A group of at least two is proposed. | `events.propose_types` |
| test | Rows that did not propose it (another condition) forecast every pair of its events twice, the event left out: with no type, and with the records of its other events teaching. The type is judged on the pairs it teaches (the two forecasts differ by at least two errors of the measured change). **Admitted** if the forecast with the type is the closer one more often than chance gives (sign test, 2.5%) and the squared error falls beyond the error of the measured changes. **Refused** otherwise. **Undecided** if fewer than two of its events move the state in those rows. | `account.test_types` |

A proposed type enters no forecast until it is admitted.

## On a planted twin, where the answer is known

Real control cells (their states and their 24 antibody columns), made-up events: two types (five and four events),
a group of three units whose cells stand to one side of the state as control cells do and that each do something
else (the rule should propose it and the test should refuse it), three units that each move the state their own
way, six own-gene units, twenty plain units. One pair in 25 holds a large effect the graph does not know.

The rule of the test was settled while looking at worlds 0 to 7. Sixteen more worlds were then run once.

| in 16 worlds run once | |
|---|---|
| each planted type proposed with exactly its events | 16 of 16, and 16 of 16 |
| any other type proposed | 0 |
| the planted types admitted by the other condition | 12 of 16, and 15 of 16 |
| the group that is no type: proposed, then refused | 16 of 16 |
| a made-up type of unrelated events refused | 16 of 16 |
| **a wrong type admitted** | **0** (0 in all 24 worlds) |

The rule errs on one side: it never admitted a type that was not there, and it refused a true one in 5 of 32 tests
(a type's records also pass on the large effects that single events have of their own).

## On the real system, told of no type (rows read before)

| told of no type | proposes | asked of | |
|---|---|---|---|
| the library under interferon (206 events) | **IFNGR1, IFNGR2, JAK1, JAK2, STAT1** | the co-culture | **admitted**: on the 48 pairs it teaches the error falls from 0.45 to 0.09, closer in 44 |
| | | no signal given | undecided: these events do not move the state there |
| the co-culture (38 events) | **IFNGR1, IFNGR2, JAK1, JAK2, STAT1** | the library | refused: it teaches 7 of 99 pairs |
| no signal given (32 events) | nothing | | |

**From either condition the rule proposes exactly the five events the author had declared, and no other.** In the
library two of them are flagged (IFNGR2, JAK1) and three join by direction (cosines 0.94 to 0.97). Five more events
move the state there (among them CCND1 and CDK6, which move it the same way); none is flagged, because the readout
reads them right, so no type is proposed for them.

The library refuses what the co-culture admits, and that is consistent: the library's control cells have seen
interferon, so its readout already reads these events (error 0.11 without the type), and the type has little to
teach there. A type is needed where the readout misreads.

**This is retrospective.** The author's declaration came from the same knowledge of the same data. It shows that
the rule can replace the hand; it is not a test on data never read.

## What this does not show

- A type nobody knew. The rule found the one known type and no other in this system. Whether it finds a new one
  needs units never read that move the state: unread units of the THP-1 screen, or another study.
- Version 6 is unchanged, and the final audit is to be read under it. This rule is a candidate for the version
  after.

## Files

| | |
|---|---|
| `anm_diag/events.py` | `propose_types` |
| `anm_diag/account.py` | `Law.shifts`, `Law.moving`, `test_types` |
| `tests/test_types_by_rule.py` | two tests |
| [`scripts/by_rule.py`](scripts/by_rule.py) | `--part twin`, `--part real` |
| [`results/twin.json`](results/twin.json), [`results/real.json`](results/real.json) | |

```
python scripts/by_rule.py --part twin --library LIBRARY.npz --cells batch4_cells.tsv --states B2 B3 B4 --packs PACK_A PACK_B PACK_4
python scripts/by_rule.py --part real (the same)
```
