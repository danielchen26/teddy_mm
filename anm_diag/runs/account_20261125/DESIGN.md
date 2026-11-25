# The account of a change: why the event graph showed no weight, and the loop rebuilt on what each reading can prove

25 November 2026. The owner asked three things. Why did nothing that is specific to the typed event graph (the TEG) show
any weight? Whatever is missing from first principles in each of the six links of the loop is to be put back, as a
repair of the foundation and not as a patch. And the whole framework is to be looked at again.

This document is the design. Sections 2 and 3 report an exploration on development cells that earlier runs have read;
it is **explored, not counted**. Sections 4 to 6 are the design that follows from it. Section 7 declares what is read
after this document is committed, and Section 8 the forecasts. Nothing is trained. No final unit and no unread unit is
opened.

*u* is the share of a measured change that is left unexplained. A step is 0.10. "Stage" means any array a head could
read: a view of the model's input, a step inside the model, the state its head reads.

## 1. The answer in short

Two causes, one on top of the other.

1. **The measure could not see what a stage holds.** Since the run before this one, a repair was read by letting a
   readout follow the mean change of an event on one half of its cells and reading the other half. That number is near
   one whenever the event's shift of the stage can be found again. It says nothing about whether the stage holds the
   change. On a planted system it calls an event closed whose change no stage holds (Section 2). Read with it, every
   stage looks alike, so a typed state could not show a gain.
2. **The events that need the graph were never targets, and the graph has no edges.** The targets of forty runs were
   knockouts that act through the transcriptome. For those, a readout fitted on control cells alone already reads the
   change, so records and structure have nothing to add. The events with the largest changes that nothing reads are
   knockouts of the gene behind the measured protein itself. Their change is not in the RNA. Only a relation between
   the event and the readout can say anything about them, and the TEG as built is a table of records with labels. The
   one place where that relation was written down, rule S3, treats it as a risk of leakage and keeps such events out.

Read with a measure that can see it, the typed state does carry more (Section 3). Given an edge that was known
beforehand, the events that act on the readout itself are forecast (Section 3). Both were there to be found.

## 2. The first cause: the measure

### 2.1 A planted system

`anm_diag/toy_account.py` plants an event `m` that moves the state along a direction the readout does not follow and
changes the readout directly. No stage holds its change. With three quarters of its cells edited:

| | u |
|---|---|
| a readout of control cells | 0.94 |
| the repair of the run before (a readout made to follow the event's mean, read on the other half) | **0.24** |
| what is left among cells whose state is alike (the evidence inside the arms) | 0.96 |

The repair "closes" three quarters of a change the state does not hold. What it measures is that the event's shift of
the state is found again in the other half.

### 2.2 The same on IRF1 (THP-1, PD-L1; explored)

| the head reads | **carried**: held by the evidence inside the arms | made to follow the mean (run 40's rule) | read by a readout of control cells | forecast from the four pathway events |
|---|---|---|---|---|
| `Z`, the mean state | **52%** [37, 66] | 90% [65, 120] | 10% [−12, 30] | 42% [20, 61] |
| `T`, typed by the TEG (the genes of the other recorded events) | **73%** [59, 90] | 95% [71, 121] | 26% [1, 46] | 49% [25, 70] |
| `R`, genes at random | 57% [42, 70] | 83% | 7% | 35% |
| `P`, programs of the control cells | 50% [39, 64] | 77% | 2% | 31% |
| `V`, the input view | **89%** [71, 112] | 97% [69, 131] | 3% | 64% [34, 98] |

Shares of the measured change (−0.28), by the module of Section 6. On paired resamples:

- the input view carries 37 points more than the mean state [17, 61]: a step;
- the typed state carries 21 points more than the mean state [8, 39] and 17 more than genes at random [7, 28]: above
  zero, not resolved as more than a step;
- made to follow the mean, the typed state and the mean state differ by 5 [−32, +39]: nothing can be seen.

For the four pathway events (the positive control of the measure) the mean state carries 85% to 106%, the input view
no more, and the forecast from the other three is 95% to 107% (STAT1: 154%).

The five last layers trained on 24 November, read the same way: what the state carries of IRF1 changes by 1 to 4
points (no step). One arm (trained with a derived record of IRF1) raises what a readout of control cells reads by 28
points [12, 48] and the forecast from the other events by 27 [13, 44]. Training moved what is easy to read; it did not
add to what a cell's state holds.

What this corrects:

- **Run 40 said** that with its record applied IRF1 closes to 0.17 at the plain mean state, that the repair and not the
  model had been the obstacle, and that the typed state adds 0.05 and never a step. **The first half holds**: the old
  rule did not apply a record. **The second half does not follow**: 0.17 is the number of Section 2.1. By the evidence
  inside the arms the mean state carries about half of IRF1's change and the input view nearly all of it. The loss
  between the input and the mean state, which run 39 had seen (0.40 against 0.18) and run 40 set aside, is there.
- **The typed state does carry more of an event of its type** than the mean state, by about a fifth of the change.
  Run 40's rule could not have shown it. For forecasting the event from the others it adds 7 points [−23, +32]:
  nothing that can be told.

## 3. The second cause: the events and the edges

The melanoma library under IFN-γ: 206 knockout events, 4,000 control cells, 50,229 cells, three protein readouts. All
development units; forty runs used five of these events.

**Three kinds of event** (u, the share left; "own record": fitted on half of the event's cells, read on the other):

| kind | events | no record | with its own record (run 40's rule) |
|---|---|---|---|
| through the state, and control cells teach it | IFNGR1, IFNGR2, JAK1, JAK2, STAT1 on PD-L1 | 0.05–0.18 | 0.02–0.06 |
| **not through the state** | CD274 on PD-L1; B2M, TAPBP, HLA-A, HLA-B, HLA-C on MHC-I | 0.92–1.09 | **0.88–1.02** |
| through the state in part, and it needs a record | IRF1 (THP-1) | 0.90 | (Section 2.2) |

- For the first kind a readout of control cells reads 84% to 100% of the change. Nothing is left for a record, a typed
  state or a trained layer to add. Every target so far was of this kind, or was IRF1.
- For the second kind the change is after the RNA. CD274's knockout lowers PD-L1 by 1.49 while its own transcript
  falls to 54% and the mean state moves by 1.2 times its noise. A readout of control cells reaches 8% from the state
  and 9% with the count of CD274 itself given beside it. With its own record the row is closed to 100% on the fitted
  cells and 0.97 is left on the others.
- Of 173 events whose own transcript can be counted, 122 lower it (median: to 57%). The mean state moves for 24. One
  site in two thousand positions does not survive a mean.

**An edge that was known beforehand.** The declarations already name the gene behind each protein
(`events.FORWARD_READOUT`, used by rule S3 to flag leakage): CD274 for PD-L1, HLA-A, HLA-B and HLA-C for MHC-I, HLA-E
for HLA-E. B2M is the second chain of both MHC molecules. Taken as an edge "the event's site is a gene of the measured
protein", seven pairs of 618 have one. The direct part of each (what a readout of control cells does not reach, in
spreads of the readout), forecast from the other pairs with an edge, each left out in turn:

| | without edges | with the edge |
|---|---|---|
| the forecast | 0 for every pair | −0.7 to −1.1 spreads |
| off by, at worst | 16.8 errors | 3.3 errors (B2M); the other six within one |
| sum of −log density over the seven | 209 | 14 |

The same quantity flags pairs with no edge whose direct part stands out. First on the list: **TAPBP on MHC-I**, 19.5
errors. Tapasin loads MHC-I; the graph did not have that edge. Then SSR2, STAT3 and CD44 on MHC-I.

## 4. First principles

### 4.1 What the loop is for

A model reads cells. Outside it there are readouts that are measured, and a typed event graph: events with records
(measured changes) and relations. The claim to earn is this: *for an event and a readout, the system (the model and
the graph together) says what the change is and how sure it is, and for every part of the change which part of the
system holds it; also for an event that has no record.*

### 4.2 One identity, three kinds of evidence

For an event *e*, a readout *y*, a stage *s* and any function *g* of the stage,

    change  =  [ mean of g(s) in the event's cells − in control cells ]  +  [ mean of y − g(s) in the event's cells − in control cells ].

The first bracket is what *g* reads. Everything depends on which *g*, and on what evidence chose it:

| name | *g* is | the evidence | what it establishes | what it cannot |
|---|---|---|---|---|
| **read** | fitted on control cells | control cells | what a head reads with no record of the event | whether the stage holds more |
| **carried** | the readout given the stage, with the arm given the stage taken out (residuals on residuals, both from the other half of the cells) | the cells of both arms, inside each arm | how much less the two arms differ among cells whose stage is alike: what the stage holds of the change, cell by cell | nothing where no cell of one arm is like a cell of the other (`overlap` near zero) |
| **others** | the control readout, changed as little as possible so that it closes the rows of other events | the rows of other events | what the records of other events say of this one: transfer | — |
| *found again* | the control readout, changed so that it closes this event's own row | this event's own mean | that the event's shift of the stage is found again | **anything about what the stage holds** |

`carried` is the definition the framework lacked. It is a statement of conditional independence (given the stage, does
the event still matter for the readout), it is estimated without ever fitting the difference between the two means,
and its estimator is orthogonal to the errors of its two auxiliary fits. Three limits are part of the design:

- **No overlap.** If the stage tells the event's cells from control cells completely, the stage holds the event's
  identity, and a head with its record reads the change. Whether the stage holds more than the identity cannot be read
  from cells. `overlap` is reported with every reading; below 0.10 the reading says "by identity" and only `others`
  can say more.
- **Per cell.** `carried` is about single cells. A stage that is a noisier measurement of a cell carries less of a
  change per cell, even if its mean over many cells moves by the right amount. That is the property a head needs, and
  it is what differs between a mean state and a typed one.
- **A marker is carried too.** If the stage holds something the event changes that is not on the way to the readout,
  the stage partly tells who was edited and so carries part of the change as information. Whether it would carry to
  another event is the question of `others`.

### 4.3 The chain

Stages in the order of the model: views of the input, steps inside, the state the head reads. Then, exactly,

    1  =  direct  +  lost before the state  +  held, not read  +  read

| part | is | the place it names |
|---|---|---|
| direct | 1 − carried(input) | no stage holds it. Only the event graph can: an edge |
| lost before the state | carried(input) − carried(state) | the model, between its input and its state |
| held, not read | carried(state) − read(state) | the head: a readout with the cells of recorded events |
| read | read(state) | nothing to do |

Every part has an interval from paired resamples. A part is open if its interval lies above a step.

### 4.4 Many events: the law of their rows

With many events the same question is asked of rows (the mean shift of a stage, the mean change of each readout), in
spreads of the readout among control cells:

    change(event, readout) = shift(event) · readout  +  direct(event, readout)  +  error

- **readout**: its prior is the ridge readout of the control cells with what that fit leaves uncertain, widened by one
  factor that the rows choose by their likelihood; the rows of the other events update it.
- **direct**: for a pair with no edge, zero with a small spread and heavy tails (so that one event the stage does not
  hold cannot bend the readout or inflate the widening; in the exploration it had inflated it to 178). For a pair
  with an edge of class *k*, a mean and a spread of its own, learned from the other pairs of that class.
- Every pair is forecast with its event left out. The forecast has three parts (read; from the other events; by its
  edge) and a spread. Its error in units of that spread is what calibration is read from.
- A pair with no edge whose direct part stands out (false discovery rate 0.05) is a **candidate**: if the event moves
  the stage, a direction the law does not cover (its record would teach the readout); if it does not, a missing edge.

## 5. The six links, rebuilt

| link | what was wrong at the foundation | what it is now | how it is checked |
|---|---|---|---|
| **1. Measure** | One number, u, that depended on the class of reader and, since run 40, on a readout made to follow the mean. Neither is a statement about the stage. | Three readings, each named by its evidence (`read`, `carried`, `others`), with `overlap`. `found again` is reported and is not a measure of holding. | The planted suite (Section 7.1): every part against the planted truth; event `m`; a permuted stage must carry nothing; an event made of control cells must have no change to account for. |
| **2. Locate** | A comparison of readers across stages, and "unreadable: no input reads it" as a dead end. | The chain of 4.3: four parts that sum to one, each with a place. "Direct" is a place: the graph. | Planted events `l` (lost before the state), `r1`, `r2` (held, not read), `n*`, `x`, `m` (direct). |
| **3. Repair** | One repair (a readout updated by rows), judged by the number it was fitted to move. | A repair is chosen by the open part and has a **ceiling**: a head cannot be credited with more than the state carries; a state cannot carry more than the input; an edge is credited by forecasts of pairs left out. A candidate state is judged by `carried` against the mean state and against a matched control state. A recorded event can always be answered from its record: that is the graph's memory, reported as such and not as the model holding it. | The gain of a stage on paired resamples; `found again` beside every reading as the ceiling's counter-example. |
| **4. An event without a record** | A rule that derives one row from one type, tested inside that type. | The forecast of 4.4: through the stage from the other events' rows, plus the class of its edge, with an interval. An event that does not move the stage and has no edge is forecast as "no change", with the measured rate at which such forecasts were wrong. | Leave one event out: the share of measured changes inside the interval, against "no change" and "a readout of control cells". |
| **5. What the loop does next** | Nine rule changes waiting for a person; no rule for which record to ask for or when the graph should grow. | The open part of each account names the lever. Candidates of 4.4 are proposals to the graph (an edge, or a record). The widening and the spreads are estimated, not set. | Planted: `x` must be flagged as a missing edge, the events that fill the library must not be. |
| **6. Confirm** | u ≤ 0.10 on items the loop was fitted to. Two of forty runs were confirmatory. | Before a sealed batch is opened: a forecast with an interval for every pair in it, from RNA and edges. After: coverage and error against the two simpler forecasts. The loop is closed for a class of events when its forecasts for events it has not seen are calibrated and beat both. | Section 9: what is prepared and what is the owner's to open. |

What happens to the nine rule changes that were waiting: Section 10.

## 6. The implementation

- `anm_diag/account.py`: `reading` (one event, one readout, one stage: `read`, `carried`, `others`, `found_again`,
  `as_trained`, `overlap`, the shift of the stage against its noise; resamples that are paired across stages and that
  draw the control readout from what its fit leaves uncertain), `gain`, `chain`, `must_fail`, `Law` (`fit`, `forecast`,
  `calibration`, `flags`).
- `anm_diag/toy_account.py`: the planted system of Section 7.1 and its suite.
- `tests/test_account.py`.
- Nothing in `cycle.py`, `outer.py` or `scope.py` is changed. The loop of the protocol runs as before; the account is
  read beside it.

## 7. What is read, and its status

### 7.1 The planted suite (declared)

Four hidden factors; the readout follows three; the state loses one the input holds; one is off in control cells.
Events: `c1 c2` (read), `r1 r2` (held, not read; `r1` forecast from `r2`), `l` (lost before the state), `m` (direct,
though it moves the state), `n1..n4` (direct, with an edge), `x` (direct, no edge), `null`, and thirty that fill the
library. Truth is the same system with very many cells and exact least squares. 73 checks per system.
Read: twenty systems as planted; five each with every cell edited, half edited, noisy cells (measurement noise five
times larger), and few cells. **The suite's first full run finished, and was seen, while this document was being
written. Its forecast (HP) is therefore not counted.** Two things were changed after it: an account now lists the
parts that are unresolved beside those that are open, and the check of the place asks that the account point at the
planted place (first open part, or first unresolved if none is open); how often the place is resolved as a step is
counted separately.

### 7.2 THP-1 (explored, not counted)

The five events on PD-L1 at `Z`, `V`, `T`, `R`, `P`, `ZV`, and for IRF1 at the five last layers trained on 24
November. The cells of Section 2.2 with the module. Reported in full; no hypothesis is counted on them.

### 7.3 The melanoma library under IFN-γ (explored, not counted)

`Law` at `Z` for 206 events and three readouts, with the edge class of Section 3: its forecasts, calibration and
candidates.

### 7.4 Three melanoma targets (declared: read after this document is committed)

The targets of run 40, by event: MHC-I in co-culture (five events), PD-L1 in co-culture (five), PD-L1 under IFN-γ
(five). Stages `Z`, `V`, `T`, `R` (the typed states are those of run 40; their gene sets were chosen on the fitted
half; `T` for an event uses the genes of the other four). These cells were read by run 40 with its repair; they have
not been read with `carried`.

## 8. Hypotheses and forecasts for 7.4, and the forecast for 7.1

| | hypothesis | confirmed if | refuted if | forecast |
|---|---|---|---|---|
| HP | the planted suite | at most 2% of the checks on the twenty systems as planted are wrong | more than 5% | 0.80; **seen before the commit, not counted** (7.1) |
| HA1 | the mean does not see | for at least four of the five MHC-I events, `found again` at `Z` exceeds `carried` at `Z` by more than a step (paired) | for fewer than two | 0.75 |
| HA2 | the typed state carries more | `carried(T) − carried(Z)` is a step for at least three of the five MHC-I events | for none, or the median gain is below 0.03 | 0.35; median gain between 0.02 and 0.20 |
| HA3 | it is the type | `carried(T) − carried(R)` is a step for at least three of the five MHC-I events | for none | 0.40 |
| HA4 | where MHC-I is lost | "lost before the state" is open (a step) for at least three of the five MHC-I events | for none | 0.45; median share between 0.05 and 0.30 |
| HA5 | PD-L1 has nothing open | of the ten PD-L1 accounts, at least seven have no open part | fewer than four | 0.55 |
| HA6 | the controls | at `Z`, for all fifteen: a permuted stage carries less than 0.15 in size, and an event made of control cells has no change to account for | any fails | 0.85 |
| HA7 | the cells can be read | `overlap` at `Z` is at least 0.10 for all fifteen | for fewer than twelve | 0.85 |
| HA8 | transfer inside the type | `others` at `Z` is at least 0.75 for at least ten of fifteen | for fewer than six | 0.50 |

An account is made only where the change stands out (at least four of its errors); a pair below that is counted as
"no change to account for" and leaves the denominators above.

## 9. The last link: what is prepared, and what is the owner's

The forecast of Section 4.4 can be written for any event from its RNA and its edges alone. Two batches would test what
nothing here can:

- **CD274 in THP-1** (and the other THP-1 knockouts no run has read). The edge class learned in melanoma forecasts the
  direct part of CD274 on PD-L1 in another kind of cell: a decrease, of the size of its class, wide. This is the test
  of an edge across cell types.
- **The unread melanoma units.**

Writing those forecasts needs the RNA of unread units. Opening them is the owner's decision and is not taken here.

## 10. The nine rule changes that were waiting

| proposal | under the account |
|---|---|
| the round's reader brings a stage to one scale | done inside every reading |
| the repair applies its rows, with the control fit widened | the widening is chosen by the rows of many events, with a direct part that keeps one event from bending the readout; a single event's own row is no longer a measure |
| gate 3 for a repair that closes its rows | replaced by the ceiling: a repair is not credited beyond what the stage carries |
| a reader class that follows the mean change | withdrawn: that reader is `found again` |
| the five earlier ones | unchanged; they concern the protocol's round and remain the owner's decision |

The protocol of sealed readings is not changed by this run. Whether the account replaces the round's reading in the
protocol is the owner's decision; this run supplies the evidence.

## 11. Threats, named before the reading

| threat | what is done |
|---|---|
| `carried` is biased low by measurement noise in single cells | stated as a property of the definition (4.2); the planted variant with noisy cells reports how far |
| `carried` at a stage with more columns is shrunk more | `R` and `P` have the columns of `T`; `ZV` is reported beside `V` |
| the linear auxiliary fits miss a relation that is not linear | the planted system is linear; on real cells the readings are of linear fits and said to be |
| the typed states were chosen with the RNA of the fitted half | they are read on the other half, as in run 40 |
| the exploration chose the measure on THP-1 | THP-1 is not counted; the melanoma targets are the declared reading; they hold no event unlike its type |
| seven pairs with an edge are few, and one class | reported as explored; the test is Section 9 |
| hyper-parameters of the law are fitted with every event, the forecasts leave one out | stated; the planted suite scores the same procedure |

## 12. Amended after the readings (25 November, after Section 7.4 was read)

The sections above are as they were committed before the reading, except this one. Three things were changed after
it; the hypotheses are scored as declared (`results/melanoma_as_first_read.json`, `scripts/score.py`).

1. **The first part of the chain is not "direct", and its place is not "an edge".** For MHC-I in co-culture the
   input view and the state carry about half of a pathway knockout's change cell by cell, and the other four events
   forecast the fifth in full (91% to 110%). A protein that is measured beside the RNA is not in a cell's RNA in full,
   and its mean change still follows from the mean shift of the stage by a law that holds across events. The chain
   is about single cells; whether the event graph must supply something is read across events. The part is now
   called **not in the input**; an account says beside it what the other events forecast; it needs an edge only if
   they do not (Section 4.4 and `Law` answer that for many events).
2. **The input holds at least what the state holds.** For PD-L1 under IFN-γ the state carries more than the input
   view `V` (97% against 80% to 96%). `V` is one view of the input. The chain now takes for the input the view that
   carries most, and the state itself if it carries more than every view; nothing is then lost before the state.
3. **Section 2.1** was written from a first planted system of the exploration (one seed, another layout). The script
   of this run, on the planted system of the module (five seeds): no record 1.01, run 40's rule 0.43, left among
   cells alike 0.94. The rule "closes" more than half of a change the state does not hold; a row closed exactly
   (`found again`) reads 98% of it.

The numbers of Section 3 are those of the exploration as it was run; `scripts/explore_library.py` repeats it and
differs in the second decimal where halves of cells are drawn (its values are in the README).
