# The final audit: read once, passed

27 November 2026. Run 49. The owner gave the word on 27 November. The design, the forecasts and the scripts that
read and score were each committed before what they depend on was opened, and the final cells were exported after
the last of them. Every part was read once. Nothing was trained, and no rule was changed.

| | commit | time (EST) |
|---|---|---|
| the design, the units, the forecast code, the preflight on cells read before | `b6bf10d` | 02:03 |
| the export, the reading and the scoring, as written | `4b70c99` | 02:08 |
| the forecasts of 2,564 pairs | `e79d32f` | 02:28 |
| the rehearsal on made-up values: nothing changed | `62ec0b7` | 02:36 |
| the export on the owner's machine, then the reading | | 02:36 |

## The result

133 final units (64 in the co-culture, 63 with no signal given, 6 in THP-1), 2,564 pairs, under version 6 as frozen
at `b4f6196`. The melanoma units are read against control cells no run had read (621 and 1,333).

| | asked beforehand | measured | |
|---|---|---|---|
| **F1, the intervals hold** | at least 90% of the pairs inside the 95% interval; the median of squared *z* between 0.23 and 0.91 | **2,443 of 2,564 (95.3%); 0.41** | **holds** |
| **F2, the loop corrects** | on the pairs where the systems can differ, version 6 better than the model alone beyond chance, and not worse than version 5 | error **0.336 → 0.140** on 160 pairs; the fall of the summed squared error 14.97, interval [0.54, 35.95]; against version 5: 0.34, [−0.04, 0.89] | **holds** |

**The audit is passed if F1 and F2 hold: passed.** The author had put 0.6 on it.

Reported beside it, as the design said:

| | | measured | |
|---|---|---|---|
| F3 | an edge has the size of its pair | 5 pairs: error 0.48 against 1.77 for the model alone; 5 of 5 inside | holds |
| F4 | the arithmetic of the readout (rule 6a) on data never read | 84 of 94 pairs inside; 85 were asked | **does not hold** |
| F5 | an edge proposed elsewhere holds here | 5 of 8 edges flagged under interferon are admitted | holds |

| stratum | pairs | inside (v6) | median *z*² | error: no change, model alone, v5, v6 |
|---|---|---|---|---|
| co-culture | 1,280 | 1,209 (94.5%) | 0.42 | 0.136, 0.122, 0.090, 0.089 |
| no signal given | 1,260 | 1,212 (96.2%) | 0.40 | 0.134, 0.127, 0.107, 0.106 |
| THP-1 | 24 | 22 (91.7%) | 0.38 | 0.175, 0.152, 0.152, 0.152 |
| all | 2,564 | 2,443 (95.3%) | 0.41 | 0.135, 0.125, 0.099, 0.099 |

## The edges

| pair | the model alone | version 6, with its interval | measured |
|---|---|---|---|
| B2M → MHC-I, co-culture | −0.40 | −4.30 [−6.48, −2.12] | −3.48 |
| CD58 → CD58, no signal | 0.00 | −2.07 [−3.01, −1.14] | −2.11 |
| CD59 → CD59, no signal | −0.04 | −1.48 [−2.80, −0.16] | −1.35 |
| HLA-A → MHC-I, co-culture | −0.04 | −0.57 [−2.75, 1.60] | −0.02 |
| HLA-C → MHC-I, no signal | −0.06 | −0.62 [−2.01, 0.77] | −0.22 |

The size of an edge is its pair's share of the scale, carried from the condition where it was measured. The model
alone says almost nothing changes. Two pairs at the background are inside as well (B2M → HLA-E: −0.21 forecast,
−0.23 measured; KDR → CD309: −0.21 and −0.27).

Edges flagged under interferon and asked here for the first time: admitted SSR2 → CD9 and SSR2 → MHC-I
(co-culture), TAPBP → MHC-I, ITGA3 → CD49f and STAT3 → MHC-I (no signal); refused SSR2 → CD59, SSR2 → CD44
(co-culture) and DAG1 → PD-L1 (no signal).

## How to read it

**What is established.** In one study and for one model, for knockouts, the event graph forecasts protein changes
never read with intervals that hold, and reduces the model's error where it says it can.

**F2 rests on eight units.** The graph has something to say for B2M and HLA-A in the co-culture and for CD58,
CD59, HLA-C, KDR, JAK1 and STAT1 with no signal. Version 6 is closer than the model alone for six of them and
further for two. The gain comes mostly from B2M, CD58 and CD59, and the lower end of the interval is just above
zero.

**On the other 2,404 pairs the three systems give the same values** (error 0.095; 0.103 for "no change"). The
graph says nothing about them. What passed there is the width of the intervals: 95.6% inside.

**Version 5 would have passed as well** (96.1% inside, median 0.36; error 0.147 on the 160 pairs). The audit
establishes the event graph as a way of working. It does not show that version 6 adds much to version 5: its
intervals are closer to exact in the bulk, and its values are a little better where the arithmetic applies.

**F4 does not hold, and it is the rule version 6 added.** Four of the five units that remove a protein of the panel
are as the arithmetic says (74 of 76 pairs inside). B2M in the co-culture is not: 10 of its 18 other pairs are
inside. Its values are still closer with the arithmetic than without (error 0.121 against 0.135 on the 94 pairs).

**THP-1.** Its six units are forecast with no edge, as said beforehand. The two pairs of a unit that removes its own
protein stand outside (CD86: −0.50; PD-L2: −0.49); the other 22 are inside.

## The accounts

14 pairs changed by four errors or more. For 13 the largest part of the change is not in the input: the RNA does
not show it, which is what an edge on the graph is for (seven of the 14 are pairs with an edge: three declared,
four admitted here). One is unresolved (CUL3 → PD-L1 in THP-1). One loses a part before the state (B2M → CD47:
35%). Of their 23 open parts, 7 are on the graph, 9 at the margin and 7 open (6 for the graph, 1 for the model).

## Explored after the reading (no test)

Why the arithmetic fails for B2M. The eight pairs outside split in two: proteins there is much of rise more than
computed (CD47 +0.44 for +0.09 forecast; CD9 +0.54 for +0.28; CD59 +0.51 for +0.26), and proteins at their
background do not rise at all (CD184, CD202b, CD140a, CD172a: +0.01 to +0.13 for +0.21 to +0.42 forecast). Rule 6a
holds the counts of the other antibodies fixed and moves only their mean. The pattern is what a fixed number of
reads per cell would give: when the most abundant protein is removed, the other counts grow in proportion to what
they are, and a column of background does not. That is a candidate for the next version (the arithmetic of the
counts, and not only of the ratio). It was seen after the reading and is no result.

## What this does not show

Another study, another model, donors, disease, time. That the model itself was corrected: it was not. The 73
reserve donors are unread.

## What is left

188 units are unread and not final. No data of this study is left on which nothing was fitted except those.

## Files

| | |
|---|---|
| [`DESIGN.md`](DESIGN.md) | what is asked, written first; what was changed from the draft and why |
| [`results/roles.json`](results/roles.json) | the 133 units |
| [`results/preflight.json`](results/preflight.json) | the code on cells read before |
| [`results/forecasts.json`](results/forecasts.json) | 2,564 pairs by three systems, and the sets of the hypotheses |
| [`results/reading.json`](results/reading.json) | every pair as measured; the edge tests; the records |
| [`results/accounts.json`](results/accounts.json) | the accounts of the 14 pairs with a change |
| [`results/outcome.json`](results/outcome.json) | the hypotheses and the verdict |
| `scripts/` | `roles.py`, `cells.py`, `forecast.py`, `preflight.py`, `pack.py`, `rehearsal.py`, `read.py`, `score.py` |
