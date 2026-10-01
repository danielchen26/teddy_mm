/* Shared by docs/index.html and docs/anm-loop.html.
 * 1) TERMS: one plain-language definition per technical term.
 * 2) Terms in prose get a dotted underline and show their definition on hover, focus or tap.
 *    Known abbreviations (O0, P_f, z_512, …) are marked automatically; <span class="term" data-term="…"> marks anything else.
 * 3) The guide column: the whole workflow, what you are reading now, and the key concepts.
 *    Each page sets window.GUIDE_CONFIG before loading this file. */
(function () {
  'use strict';

  var TERMS = {
    teddy: ['TEDDY', 'Merck’s single-cell RNA foundation model (TEDDY-G, 70M parameters); it does not predict proteins. Here it turns a cell’s RNA into a 512-number embedding; a small head we train on the training sites’ measured proteins predicts 134 surface proteins from that embedding, and 9 of them (12 in v3) go to ANM. TEDDY is never retrained; the head is trained once and then fixed.'],
    anm: ['ANM · Active Neural Matter', 'A measuring protocol with a decision engine: push a source, follow a declared candidate retained state, read a declared readout, and check linearity and sufficiency (our reading of the paper, not its words: a verifier of declared responses, not a better decision-maker). Here it is used as a decision layer: it reads evidence, applies a written-down question (the observer) and returns a call or an honest “no call”, with checks attached. Every event enters at t = 0, so its calls equal a re-coded declared rule (v3 E1: 0 mismatches).'],
    citeseq: ['CITE-seq', 'A technology that measures RNA and surface proteins in the same cell. The measured proteins are the answer key we hold out.'],
    adt: ['ADT · measured surface protein', 'CITE-seq’s protein readout (antibody-derived tags). Held out and used only to check answers, except where a test says otherwise.'],
    evidence: ['Typed evidence (δu)', 'A value per panel protein, scaled 0–1 and handed to ANM as labelled inputs such as “CD19 is high”. v3: 12 panel proteins (3 each for B, T, NK and myeloid, chosen on the validation donor), each our head’s prediction divided by the 95th percentile of the head’s own training predictions and clipped to [0, 1]; no measured protein enters. v2: 9 proteins scaled by one train-median size factor (0.92663) and divided by each protein’s training 95th percentile; in the withdrawn Experiment 2’s RNA-missing condition they are the protein-only stand-in’s reconstruction (no TEDDY), and in the withdrawn Experiment 6’s complement test the measured proteins. The first run scaled every prediction by a per-cell factor computed from measured protein.'],
    observer: ['Observer · the question', 'The question, written down: which markers count, how strong a signal must be, and when to decline. Changing it needs no labels and no retraining; that holds for any written rule, including the fixed rule. v3 E1: on a changed question a classifier trained only on the old question’s labels was more accurate than the zero-label rule (−0.078, loss).'],
    o0: ['Soft rule (O0)', 'The default lineage call. All 9 markers weigh equally; call the best lineage if its score reaches 0.12.'],
    o1: ['Strict rule (O1)', 'The same question asked more strictly: a weighted mean with the key markers CD19, CD3 and CD16 at weight 2 (redefined in the rerun), and the score must reach 0.28. Almost the same answers where both decide (16 of 14,433 differ); declines more often.'],
    o2: ['B/T-priority rule (O2)', 'A different question: only CD19 (B), CD3 (T) and CD16 (myeloid) count, with B and T weighted up (B ×1.5, T ×1.3, myeloid ×0.5), bar 0.20 as declared (0.1333 on the shared score scale). It changes 12.48% of the expected answers.'],
    abstain: ['Decline to call (abstain)', 'The readout says “no call” because the evidence does not meet the written-down rule: honest silence instead of a guess.'],
    pf: ['Workability (P_f)', '1 if the observer can make a call from the evidence, 0 if it declines. Averaged over cells, it is the share of cells called.'],
    qf: ['Exactness (Q_f)', '1 if a call matches the held-out protein truth. Averaged over labelled cells, a decline counts as a miss.'],
    qdec: ['Accuracy of calls made (Q)', 'Correct calls divided by calls made, checked against the held-out proteins. Used to compare methods like for like.'],
    softp: ['Confidence score (soft_P)', 'How strong ANM’s best option is for a cell. It reaches 1.11, so it is a score, not a probability. Keeping only high-score cells trades coverage for accuracy; the fixed rule’s own score margin does as well at the same coverage, so this is not a gate only ANM has. v3 E1.4a: ANM’s top score (= G(3)·3 × the rule’s top score) orders calls no better than TEDDY’s margin (AURC −0.002 [−0.024, 0.020], inconclusive), and a logistic classifier on the same evidence orders them better (0.942 vs 0.843). On external PBMC data (v3 E6) TEDDY’s margin ranked the calls better in 8 of 8 donors (−0.0087 [−0.0182, −0.0017]), and the classifier’s advantage did not replicate.'],
    coverage: ['Coverage', 'The share of cells that still get a call once a confidence cut is applied.'],
    loo: ['Leave-one-out (LOO)', 'Remove one protein’s evidence and decide again. The marker whose removal moves the score most is the reason for the call.'],
    flip: ['Flip distance', 'How much one marker’s value must change before the call switches lineage. Small means a fragile call.'],
    holdout: ['Held-out cells (site4)', 'v2: 16,750 cells from a site never used in training, not donor-independent: donor 15078 is in training at other sites and is 32.6% of them. v3: the primary test set is the 11,294 site4 cells of donors 13272 (7,365) and 19593 (3,929), which appear in no other split; donor 15078 (5,456 cells) is a secondary split that never changes a verdict. Their measured proteins score the answers; labelled v2 tests, both withdrawn, also used them as input (Experiment 2’s RNA-missing stand-in; Experiment 6’s complement test).'],
    ood: ['Out-of-site check', 'Despite the name, not another site: 2,000 cells of held-out donor 18303 at a training site (val_non_site4), which was also the validation set during training. A second check, but a weak one.'],
    train: ['TEDDY + trained classifier', 'The “just train harder” method: a small classifier fitted on labelled cells. It is a stand-in modelled on Jev’s interface (one probability per option); TypeSafe AI’s Jev itself was never called. v3: a multinomial logistic regression on the 12 evidence values, trained on 58,657 labelled training cells (C = 100 by validation log-loss). v2 inputs depended on the experiment: logistic and MLP heads on the 9 predicted proteins plus the first 32 embedding numbers (Experiments 1 and 5), heads on the 9 values only (Experiment 4, and the withdrawn Experiment 2), and in Experiment 5 also a re-tuned threshold grid on the fixed rule.'],
    fixedrule: ['TEDDY + fixed rule', 'The baseline method: the predicted proteins (our head on TEDDY’s embedding) read by one hard-wired rule. Average each lineage’s markers and call the best lineage if it reaches the bar. Written for the same question as ANM, it makes the same call on every cell (v3 E1.1a: 0 mismatches); v3 calls it the declared rule. On the v2 key its own score margin gated its calls at least as well as ANM’s confidence at every coverage; in v3 the AURC comparison was inconclusive on site4 (E1.4a), and on external data (E6) the margin ranked the calls better than ANM’s top score in 8 of 8 donors.'],
    jev: ['Jev (TypeSafe AI)', 'TypeSafe AI’s “System One” decision model (jev-1.13), closed and reached only through its API. It answers typed questions (a Choice over named options, a Score, or yes/no) with calibrated probabilities and no rationale text. TypeSafe post-trains it (RLCD); it is used zero-shot and cannot be trained or fine-tuned on your labels. This repo never calls it.'],
    z512: ['TEDDY embedding (z_512)', 'TEDDY’s 512-number summary of a cell, from TEDDY-G’s official preprocessing: counts ÷ the cell’s total × 10⁴, ÷ TEDDY’s non-zero gene medians, the top 2,048 genes as rank tokens (no CLS token), then the mean of the last-layer gene tokens, L2-normalized (z_rna.npy, from scripts/03_embed_rna.py --preprocessing official). Every number on the site uses it. The earlier non-official embedding (--preprocessing legacy: mean over real tokens at context 512, no gene-median normalisation) gave the same phase-1 Pearson (0.603). The compact 32-number export (z_keep=32) is not the probe space.'],
    mustpair: ['Must-separate pair', 'Two cells whose TEDDY embeddings are near-identical (cosine ≥ 0.98) but whose measured proteins disagree under the 3-lineage answer key, for example one keyed myeloid and one keyed T. Many of the “myeloid” cells in these pairs are NK cells, which the key has no class for (it has no NK or out-of-scope class).'],
    falseagree: ['False agreement', 'Among must-separate pairs where both cells get a call, the share given the same call even though they differ.'],
    softsep: ['Soft separation', 'The share of must-separate pairs the readout tells apart: different calls, or one cell declined.'],
    corrsep: ['Correct separation', 'The share of must-separate pairs split into the right, different calls.'],
    zplus: ['z⁺ · feature swap', 'TEDDY’s embedding with protein features glued on: the “just swap features” alternative to typed evidence.'],
    phase1: ['Phase 1 · TEDDY + small head (done)', 'Frozen TEDDY embeds each cell’s RNA once (512 numbers, official TEDDY-G preprocessing, up to 2,048 gene tokens). A small head trained in this repo maps that embedding to all 134 surface proteins, fitted with a negative-binomial loss. An MLP head beat a latent flow-matching head: test Pearson 0.603 vs 0.582 with the official preprocessing and the train-median size factor (the earlier non-official preprocessing also gave 0.603; 0.610 vs 0.596 with the old per-cell factor). Every experiment uses its predictions (labelled exceptions, both withdrawn: the RNA-missing stand-in in Experiment 2, measured protein in Experiment 6’s complement test).'],
    phase2: ['Phase 2 · fusion scaffold (unfinished)', 'A multimodal model: TEDDY’s RNA embedding plus a protein encoder, trained with one modality randomly dropped (RNA in about 18% of cells, protein in about 15%, never both) and a latent flow-matching decoder. Retrained on the official embedding after its fixes (ADT transformed once, train-median size factor). Test Pearson (134 proteins), RNA-only / protein-only / joint: 0.613 / 0.767 / 0.763 decoded directly, 0.211 / 0.216 / 0.219 from one flow-matching sample, 0.358 / 0.382 / 0.379 from the mean of 5 (legacy arm, old double transform and measured size factor: 0.622 / 0.745 / 0.738 decoded directly; scripts/phase2_decode_check.py). So the old ≈ 0.25 came mostly from scoring one flow-matching sample; the fix adds about 0.02 when protein is an input; RNA-only is 0.613 against the legacy arm’s 0.622, which used each cell’s own measured protein depth (the leak), and is now slightly above phase 1 (0.603). Protein-only and joint inputs contain the measured protein, so those are reconstructions, not predictions. It contains no ANM; the withdrawn Experiment 2 used it as the RNA-missing stand-in (and half of “both”).'],
    modeb: ['Mode B', 'What this project does: treat TEDDY as an evidence source and study the decision layer on top of it. It runs ANM’s open-loop decision layer: TEDDY is a prescribed evidence source and does not read the field, and every event enters at t = 0, so ANM’s decisions equal a re-coded declared rule and its dynamics (retention over real steps, histories, feedback) are not tested here. Not yet tested is not impossible: ANM’s step can be any real processing order, and TEDDY’s layers are one (a Mode A test).'],
    modea: ['Mode A · looking inside TEDDY', 'Studying TEDDY’s internals: its layers as ANM’s steps, the gene-mean or the per-gene residual stream as a declared candidate retained state, layer Jacobians, gene perturbations. v2 probes (Experiment 6 prototype): an MLP probe on the embedding kept about twice as much of the NK–T gap as our head for CD56, CD94 and CD3, and gene-token states more than their mean. v3 (registered): E5 followed a push through the 12 layers and stopped at its derivative check (float32; float64 passes); E5-M, the sufficiency test of the gene-mean, was inconclusive: for E3’s attention-pooling readout a smaller share of its first-order response flows through what pooling discards on look-alike cells (0.085) than on random cells (0.127), for this one observer. E7 ran the ANM paper’s inverse loop at the layer-11 cut with the consumer fixed (layer 12 + our head): the layer-11 gene-mean was rejected as a sufficient state and no revision was selected. Related: v3 E3 compared readouts of TEDDY’s layer means and token states and rejected “repairable on frozen TEDDY” (both told NK from T less accurately than our head).'],
    sufficiency: ['Sufficiency test', 'A check that a declared candidate retained state holds what a declared readout needs: the same state must give the same readout. If two inputs reach (almost) the same state but need different readouts, the candidate is missing something. Declaring a state does not establish that it is sufficient, and ANM does not “find” the state. v2 prototype: NK and T neighbours in TEDDY’s embedding whose predicted protein differences were compressed to about a quarter to a third of the measured ones. Superseded by v3 E3: on registered within-donor pairs the head’s compression on look-alikes is descriptive only (gap ratio 0.1169 flagged vs 0.3985 unflagged, overlapping intervals; 0.4729 vs 0.5137 without gdT CD158b+ pairs), and both trained readouts of frozen TEDDY lost NK-vs-T accuracy to the head. The ANM test of TEDDY’s gene-mean state (E5-M) is inconclusive: neither “pooling loses the look-alike difference” nor “the gene-mean state is sufficient” is established, for the one observer tested. E7 found exact same-state violations one layer earlier: histories with the same layer-11 gene-mean gave different readouts of TEDDY’s layer 12 + our head (rejected), and no declared revision was selected.'],
    timestep: ['ANM’s step (t)', 'Any declared processing order, not necessarily physical time; memory, feedback and delay are defined over these steps. The order must be real: the system actually processes in that order. TEDDY’s transformer layers are one (h_{l+1} = h_l + f_l(h_l)). The nine marker predictions for one cell are not: they come out of one prediction, so feeding them in panel order made memory decay act as a hidden weight.'],
    rnaonly: ['RNA (TEDDY + head)', 'Evidence from our head’s predictions on TEDDY’s embedding of the cell’s RNA.'],
    adtonly: ['RNA missing (protein-only stand-in)', 'Without RNA, TEDDY cannot run on the cell. A second model from this repo stands in: the phase-2 bidirectional model, trained on TEDDY’s RNA embedding plus measured protein, here run with its RNA input switched off. It reads the cell’s 134 measured proteins and reconstructs the 9 panel proteins: panel Pearson 0.446 in the first run (one flow-matching sample per cell); in the official prototype 0.353 from one sample and 0.945 decoded directly; either way a reconstruction of the measured protein it reads. No TEDDY embedding or head prediction is used for these cells.'],
    joint: ['Both (averaged)', 'The average of the TEDDY + head prediction and the stand-in’s protein-only reconstruction.'],
    lineage: ['Lineage call', 'The decision per cell. v2: B cell, T cell or myeloid, each judged from 3 markers (B: CD19 CD72 CD22 · T: CD3 CD2 CD5 · myeloid: CD16 CD11c CD36). v3: B, T, NK or myeloid, or no call, with 3 panel proteins per class chosen on the validation donor (B: CD20 CD22 CD268 · T: CD3 CD2 CD5 · NK: CD122 CD94 CD56 · myeloid: CD172a CD11c CD62P); out-of-scope cells should get no call.'],
    pearson: ['Pearson ≈ 0.60', 'How well our phase-1 head, reading TEDDY’s embedding, predicts all 134 measured proteins (mean correlation on held-out cells): 0.603 with the train-median size factor, 0.610 with the old per-cell factor. A reference point; we do not try to beat it.'],
    veto: ['Veto', 'The loop’s first check, a step outside ANM’s engine: it flags pairs that TEDDY’s embedding puts together (cosine ≥ 0.98) although their measured proteins disagree, so a readout built on the embedding should not be trusted there.'],
    complement: ['Complement', 'Give the decision the evidence it lacks (here, measured protein as typed evidence, in place of the TEDDY + head predictions) instead of retraining TEDDY or the head.'],
    verify: ['Verify', 'Score the same pairs again after the change: false agreement, soft separation and accuracy.'],
    finitefield: ['Finite field · ANM’s engine', 'ANM’s decision engine: it combines the typed evidence under the observer into a score per lineage, then calls or declines. Markers are linked to their lineage node, and each lineage node has a back-edge of weight 0.15 to each of its markers. All markers enter at once, so each lineage score is a fixed multiple of the fixed rule’s weighted sum (the first run entered them one per step in panel order, a hidden weight).'],
    permutation: ['Permutation test', 'Shuffle the data many times to see how often a pattern this strong appears by chance. Here the evidence values are shuffled within each lineage: the top deciding marker, CD2, decides 17.9% of calls vs 15.7% on average, p ≈ 0.0099 over 100 shuffles (the smallest p they can give).'],
    bootstrap: ['Bootstrap', 'Resample the cells with replacement (200 times) to see how stable a number is; gives a 95% interval. The v3 tests use a two-stage bootstrap instead: donors, then cells within donors, B = 2000.'],
    v3reg: ['Registered v3 tests', 'Experiments E1–E7 and E5-M whose endpoints, margins and decision rules were written down on training and validation data only (registration_v3.json, sha256 e4c8a33e…, amendments A1–A3 and per-experiment addenda, on GitHub branch exp/registered-v3) before the test data were evaluated. E1–E5, E5-M and E7 ran once on the two site4 donors in no other split, E6 once on an external dataset (Hao et al. 2021 PBMC, 8 donors; E7 also ran there as a secondary family), and all are reported whatever they said. For E1 and E5 the governing files were committed before the evaluation but pushed after it; for E2, E3, E4, E5-M, E6 and E7 they were on GitHub before the run (for E7’s selection records, by 1 s on two different clocks).'],
    declrule: ['Declared rule', 'For each class, the mean evidence over its panel; call the highest class if its score reaches the question’s bar, otherwise “no call”. Written down in advance, it needs no labels. Where ANM’s field does not run, or provably computes the same thing (all events at t = 0), this site names the method a declared rule.'],
    outscope: ['OUT · out of scope', 'Cells that are none of B, T, NK or myeloid: erythroid cells, progenitors, plasma cells, pDC, ILC. The correct action is no call. The v2 key had no such class. OUT decline rate: the share of out-of-scope cells that get no call.'],
    primkey: ['Primary key · annotation-only key', 'Primary key: a cell’s class when the annotated cell type and a protein gate on measured CD3, CD19, CD56, … agree; otherwise the cell is unscored (3,561 of the 11,294 test cells). Annotation-only key: the annotated class for every cell. The primary key is weak in donor 13272, so both are reported with equal prominence; the key is not switched, because switching would use test labels to choose.'],
    kappa: ['Key validity (kappa)', 'Agreement between the annotation and the protein gate beyond chance (Cohen’s kappa over five classes). Registered flag: below 0.85 on the test cells, every result is also shown on the annotation-only key. Test cells 0.5781; donor 13272 0.4177, donor 19593 0.8798; validation 0.8905. E6 (external data) registered 0.70 instead: pooled 0.6562, so “key not validated” and the annotation-only key decides.'],
    selacc: ['Coverage · selective accuracy · decision accuracy', 'Coverage: share of all cells that get a call. Selective accuracy: share of called, scored cells whose call is right (a call on an out-of-scope cell is wrong). Decision accuracy: share of scored cells with a right call or a no call on an out-of-scope cell.'],
    matchcov: ['Matched coverage', 'Every method calls its own most confident ⌈c·N⌉ cells, so accuracy is compared at equal coverage; ties are broken by a fixed random permutation (A1.1). This removes the first run’s “declines more, so looks more accurate” effect.'],
    aurc: ['AURC', 'Area under selective accuracy over coverage 0.95 … 0.50, divided by 0.45: how well a score orders calls from safe to risky. Higher is better.'],
    verdict: ['Win · loss · equivalent · inconclusive', 'Registered verdict words. Win: the difference reaches the margin, its 95% interval’s lower bound is above 0, and in each primary donor the difference has the same sign and is at least half the margin. Loss: the mirror image. Equivalent: the interval lies inside ± the margin. Otherwise inconclusive. Intervals: two-stage bootstrap (donors, then cells within donors), B = 2000, seed 1.'],
    gapratio: ['Gap ratio · flagged pair (E3)', 'Within NK–T neighbour pairs, the median predicted protein difference divided by the median measured difference, averaged over CD56, CD94, CD335 and CD3: 1 keeps the whole measured gap, small values compress it. A pair is flagged when its cosine in TEDDY’s embedding is at least 0.949402 (set on the validation donor).'],
    jvp: ['Jacobian-vector product (JVP)', 'The exact first-order change of each layer’s state for a small push at TEDDY’s input, computed by forward-mode derivatives through each of the 12 layers. E5’s registered gate required central finite differences to agree with it (relative error ≤ 0.05, cosine ≥ 0.99) in at least 95% of cases; in float32 at ε = 1e-3 they agreed in 0.854.'],
    candidate: ['Declared candidate retained state', 'In ANM v2 the retained state is a candidate coordinate declared before the run and then tested against a declared readout; its choice does not establish sufficiency. ANM does not “find” the state, and the readout-visible quotient is only a diagnostic. Here the candidate is TEDDY’s gene-mean embedding (inside TEDDY, the per-gene residual stream); E5-M tested it against declared observers, with an inconclusive registered verdict, and E7 rejected the layer-11 gene-mean for TEDDY’s layer 12 + our head, with no revision selected.'],
    replic: ['Replication rule (E6)', 'A v3 conclusion replicates on the external data if its sign holds in at least 5 of the 8 external donors and the pooled 95% interval excludes 0 on the same side. The sign comes from E1’s site4 point; a site4 point of exactly 0 gives “no v3 direction”. E1’s own decision rule (win / loss / equivalent / inconclusive, each external donor as a primary donor) is reported beside it and never replaces it.'],
    inverse: ['Inverse loop (PR #16)', 'The ANM paper’s v2 revision (PR #16, under review) tests a description of a model in three steps, with the consumer fixed. Reject: histories with exactly the same candidate state give readouts that differ beyond the tolerance. Select: choose a revised description from a declared candidate library, by a declared rule, on development data. Validate: test it on fresh comparisons. Selection works only within the library; it does not discover a missing variable. On TEDDY, E7 ran it at the layer-11 cut: the rejection held, no revision was selected, and nothing was validated.'],
    meanpatch: ['Mean-preserving state intervention (E7)', 'An exact change of TEDDY’s layer-11 token states that leaves their mean unchanged: disjoint pairs of token states are each replaced by their average. Two such histories share the gene-mean, so a readout difference beyond the tolerance rejects the gene-mean as a sufficient state. The changed states are not ones a real cell produces (off the data manifold): an intervention on internal states, not biological variation.'],
    clamp: ['Clamp estimator (H)', 'ANM’s estimator of how a declared observer depends on the candidate state: H = Λ·pinv(Γ), computed from clamps of the candidate state alone and never fitted on the tested responses. The residual r = ‖K_O − H·K_z‖ / ‖K_O‖ is the share of the observer’s response that flows through what the candidate discards.'],
    stacker: ['Stacker · learned fusion (E4)', 'F3 stacker: a logistic regression on the 8 class scores of the two channels, trained on the validation donor per noise level. F4 learned fusion: a logistic regression on the 24 evidence values, trained on training cells without noise. Both use labels; the declared averages F1 and F2 do not.'],
    trustavg: ['Trust-weighted average (F2)', 'A declared fusion rule: average the two channels’ class scores weighted by each channel’s validation trust (mean of max(0, Pearson) with the measured panel). ANM’s fusion without its contradiction events equals it up to a constant.']
  };
  window.TEDDY_ANM_TERMS = TERMS;

  // apply the remembered guide width / hidden state now, before the page draws its charts
  try {
    var savedW = parseInt(window.localStorage.getItem('teddyGuide.width'), 10);
    if (savedW) document.documentElement.style.setProperty('--gw', Math.max(260, Math.min(Math.max(260, Math.min(640, window.innerWidth - 640)), savedW)) + 'px');
    if (window.localStorage.getItem('teddyGuide.collapsed') === '1') document.documentElement.classList.add('g-collapsed');
  } catch (e) {}
  // no slide animation while the saved layout is applied on load
  document.documentElement.classList.add('g-notrans');
  window.addEventListener('load', function () { requestAnimationFrame(function () { requestAnimationFrame(function () { document.documentElement.classList.remove('g-notrans'); }); }); });

  // abbreviations marked automatically in prose (headings, buttons, links and code are skipped)
  var AUTO = [
    ['TEDDY + trained classifier', 'train'], ['TEDDY + fixed rule', 'fixedrule'], ['Train + Jev', 'train'], ['Train+Jev', 'train'],   // method names (old names kept for old links): link to the method, not to TypeSafe's Jev
    ['B/T-priority rule', 'o2'],
    ['O0', 'o0'], ['O1', 'o1'], ['O2', 'o2'], ['soft_P', 'softp'], ['z_512', 'z512'], ['z⁺', 'zplus'], ['LOO', 'loo'], ['δu', 'evidence'],
    ['Mode B', 'modeb'], ['Mode A', 'modea'], ['CITE-seq', 'citeseq'], ['ADT', 'adt'], ['P_f', 'pf'], ['Q_f', 'qf'], ['Jev', 'jev'],
    ['site4', 'holdout'], ['adt_only', 'adtonly'], ['rna_only', 'rnaonly'], ['AURC', 'aurc'], ['JVP', 'jvp'],
    ['must-separate', 'mustpair'], ['must-pairs', 'mustpair'], ['must-pair', 'mustpair'], ['false_agree', 'falseagree'], ['soft_sep', 'softsep'], ['correct_sep', 'corrsep']
  ];
  var byToken = {};
  AUTO.forEach(function (a) { byToken[a[0]] = a[1]; });
  function escRx(s) { return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }
  var RX = new RegExp('(^|[^A-Za-z0-9_])(' + AUTO.map(function (a) { return escRx(a[0]); }).join('|') + ')(?![A-Za-z0-9_])', 'g');
  var SKIP = 'a,button,summary,code,pre,svg,script,style,textarea,h1,h2,h3,.term,.topbar,.guide,#tip,#termtip,.yaml,.no-terms,.seg,.ftabs,.bn,.brand';

  function autowrap(root) {
    if (!root) return;
    var walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode: function (n) {
        if (!n.nodeValue || !RX.test(n.nodeValue)) { RX.lastIndex = 0; return NodeFilter.FILTER_REJECT; }
        RX.lastIndex = 0;
        var p = n.parentElement;
        return p && !p.closest(SKIP) ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_REJECT;
      }
    });
    var nodes = [], n;
    while ((n = walker.nextNode())) nodes.push(n);
    nodes.forEach(function (node) {
      var text = node.nodeValue, frag = document.createDocumentFragment(), last = 0, m;
      RX.lastIndex = 0;
      while ((m = RX.exec(text))) {
        var start = m.index + m[1].length;
        if (start > last) frag.appendChild(document.createTextNode(text.slice(last, start)));
        frag.appendChild(makeTerm(m[2], byToken[m[2]]));
        last = start + m[2].length;
      }
      if (last < text.length) frag.appendChild(document.createTextNode(text.slice(last)));
      node.parentNode.replaceChild(frag, node);
    });
  }
  function makeTerm(text, id) {
    var s = document.createElement('span');
    s.className = 'term abbr'; s.textContent = text; s.setAttribute('data-term', id);
    decorate(s);
    return s;
  }
  function decorate(el) {
    if (el.__term) return;
    el.__term = true;
    el.setAttribute('tabindex', '0');
    el.setAttribute('role', 'button');
    var t = TERMS[el.getAttribute('data-term')];
    if (t) el.setAttribute('aria-label', el.textContent + ': ' + t[0] + '. ' + t[1]);
  }
  function decorateAll(root) { [].slice.call((root || document).querySelectorAll('.term[data-term]')).forEach(decorate); }

  // ---------- tooltip ----------
  var tip = document.createElement('div');
  tip.id = 'termtip'; tip.setAttribute('role', 'tooltip');
  document.body.appendChild(tip);
  var cur = null;
  function place(el) {
    var r = el.getBoundingClientRect(), tw = tip.offsetWidth, th = tip.offsetHeight;
    var x = Math.min(Math.max(8, r.left + r.width / 2 - tw / 2), window.innerWidth - tw - 8);
    var y = r.bottom + 8;
    if (y + th > window.innerHeight - 8) y = r.top - th - 8;
    tip.style.left = Math.round(x) + 'px'; tip.style.top = Math.round(Math.max(8, y)) + 'px';
  }
  function show(el) {
    var t = TERMS[el.getAttribute('data-term')];
    if (!t || (el.checkVisibility && !el.checkVisibility())) return;
    if (cur && cur !== el) cur.classList.remove('on');
    tip.innerHTML = '<b>' + t[0] + '</b><span>' + t[1] + '</span>';
    tip.classList.add('on'); el.classList.add('on'); cur = el;
    place(el);
  }
  function hide() { tip.classList.remove('on'); if (cur) cur.classList.remove('on'); cur = null; }
  function termOf(e) { return e.target && e.target.closest ? e.target.closest('.term[data-term]') : null; }
  document.addEventListener('pointerover', function (e) { var t = termOf(e); if (t && e.pointerType !== 'touch') show(t); });
  document.addEventListener('pointerout', function (e) { var t = termOf(e); if (t && e.pointerType !== 'touch' && !(e.relatedTarget && t.contains(e.relatedTarget))) hide(); });
  document.addEventListener('focusin', function (e) { var t = termOf(e); if (t) show(t); });
  document.addEventListener('focusout', function (e) { if (termOf(e)) hide(); });
  // A tap focuses the term (focusin shows the tip) before its click arrives, so the click toggles on the
  // state the tip had when the pointer went down; otherwise a tap on a phone would open and close it at once.
  var downOpen = null;
  document.addEventListener('pointerdown', function (e) { var t = termOf(e); downOpen = t ? cur === t : null; }, true);
  document.addEventListener('click', function (e) {
    var t = termOf(e), was = downOpen;
    downOpen = null;
    if (t) { if (was === null ? cur === t : was) hide(); else show(t); return; }
    if (cur) hide();
  });
  document.addEventListener('keydown', function (e) {
    var t = termOf(e);
    if (e.key === 'Escape') hide();
    else if (t && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); if (cur === t) hide(); else show(t); }
  });
  window.addEventListener('scroll', function () { if (cur) place(cur); }, { passive: true });
  window.addEventListener('resize', function () { if (cur) place(cur); });

  // ---------- guide ----------
  // Per section the column shows: where it sits in the workflow, how ANM decides there (from the page's
  // GUIDE_CONFIG.sections[id].how), one real cell run through ANM's engine, and the terms that section uses.

  // The worked example is the v0 demo cell under the soft rule. Inputs, scores and the CD16 flip come from
  // outputs/anm_cite_bridge/bakeoff/bakeoff_results.json; the per-step levels from scripts/guide_engine_trace.py,
  // which re-runs ANM's engine and checks that the last step equals the published scores.
  var ENGINE = {
    title: 'Inside ANM’s engine',
    sub: 'one real cell, step by step (first run)',
    steps: [
      ['Evidence in', 'Each of the 9 predicted proteins (our head on TEDDY’s embedding) is divided by its 95th-percentile value in training cells: 0 = absent, 1 = as high as it gets. The question then weights them (the soft rule leaves them as they are; the strict rule doubles CD19, CD3 and CD16).'],
      ['Wire it up', 'Every marker becomes a node linked to its lineage: B ← CD19 CD72 CD22 · T ← CD3 CD2 CD5 · myeloid ← CD16 CD11c CD36. Each lineage node also has a back-edge of weight 0.15 to each of its markers, so some of its level flows back to them.'],
      ['Let it flow', 'In this first-run example, markers enter one per step in panel order, each as a pulse of its value. Every step, each node keeps 82% of its level and each lineage node takes in 16% of its markers’ levels. The field runs 4 more steps after the last marker. That entry order was a hidden weight (markers that enter early have decayed more by read-out). The corrected rerun enters all markers at once, so each lineage score is a fixed multiple of the rule’s and ANM makes the rule’s call on every cell.'],
      ['Read out', 'A lineage’s final level is its score. The top score becomes the call if it reaches the question’s bar (soft 0.12 · strict 0.28 · B/T-priority 0.20 as declared); otherwise “no call”.'],
      ['Check', 'The cell’s measured proteins, which the engine never sees, give the answer key. Leave-one-out and flip distance re-run this same engine.']
    ],
    cell: 'cite_site4_73511',
    inputs: [['B', [['CD19', 0.111], ['CD72', 0.153], ['CD22', 0.128]]], ['T', [['CD3', 0.172], ['CD2', 0.477], ['CD5', 0.184]]], ['Myeloid', [['CD16', 0.931], ['CD11c', 0.413], ['CD36', 0.104]]]],
    threshold: 0.12,
    // lineage levels after each step: [B, T, myeloid]
    levels: [[0.0178, 0, 0], [0.0537, 0, 0], [0.0966, 0, 0], [0.1228, 0.0276, 0], [0.1371, 0.1215, 0], [0.1433, 0.2104, 0], [0.1439, 0.2643, 0.149],
      [0.1407, 0.2936, 0.3104], [0.1351, 0.3059, 0.4266], [0.1281, 0.3065, 0.4933], [0.1202, 0.2993, 0.5254], [0.1121, 0.2872, 0.5339], [0.1039, 0.272, 0.5263]],
    ex: {
      threshold: 0.12,
      groups: [
        { label: 'All 9 markers', bars: [['B', 0.104], ['T', 0.272], ['Myeloid', 0.526]], call: 'Myeloid' },
        { label: 'Leave CD16 out', bars: [['B', 0.104], ['T', 0.272], ['Myeloid', 0.189]], call: 'T' }
      ],
      caption: 'Called myeloid, and the measured proteins agree. Leave CD16 out and the myeloid score falls by 0.338, so T wins: CD16 is the deciding marker. Its flip distance is 0.731: CD16 at 0.20 or lower would already flip the call.'
    }
  };

  function esc(s) { return String(s).replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;'); }
  function stepList(steps) {
    return '<ol class="g-howl">' + steps.map(function (st) { return '<li><b>' + st[0] + '</b><span class="d">' + st[1] + '</span></li>'; }).join('') + '</ol>';
  }
  // final scores as small bars, one group per scenario; the winning lineage is filled, the bar is marked
  function exampleBars(ex) {
    if (!ex) return '';
    var thr = (ex.threshold * 100).toFixed(1) + '%';
    var groups = ex.groups.map(function (g) {
      var rows = g.bars.map(function (b) {
        var win = b[0] === g.call;
        return '<div class="g-bar' + (win ? ' win' : '') + '"><span class="l">' + b[0] + '</span><span class="tr"><i style="width:' + (b[1] * 100).toFixed(1) + '%"></i><b class="thr" style="left:' + thr + '"></b></span><span class="v">' + b[1].toFixed(3) + '</span></div>';
      }).join('');
      return '<div class="g-grp-bars"><div class="g-ex-h"><span>' + g.label + '</span><em>→ ' + g.call + '</em></div>' + rows + '</div>';
    }).join('');
    return '<div class="g-ex">' + (ex.title ? '<div class="g-ex-t">' + ex.title + '</div>' : '') + groups +
      '<div class="g-ex-leg"><span><b class="thr"></b>bar ' + ex.threshold.toFixed(2) + '</span><span><i></i>called lineage</span></div>' +
      (ex.caption ? '<p class="g-ex-cap">' + ex.caption + '</p>' : '') + '</div>';
  }
  function inputGrid(inputs) {
    return '<div class="g-inputs">' + inputs.map(function (g) {
      return '<div class="g-in"><span class="lin">' + g[0] + '</span>' + g[1].map(function (m) {
        return '<span class="mk"><span class="nm">' + m[0] + '</span><span class="mv">' + m[1].toFixed(3) + '</span><span class="mb"><i style="width:' + (m[1] * 100).toFixed(1) + '%"></i></span></span>';
      }).join('') + '</div>';
    }).join('') + '</div>';
  }
  // levels of the three lineage nodes over the steps, with the bar; hover or tap a step to read it
  var TR = { W: 284, H: 164, l: 26, r: 50, t: 10, b: 34, max: 0.6 };
  function traceSvg(E) {
    var n = E.levels.length, pw = TR.W - TR.l - TR.r, ph = TR.H - TR.t - TR.b, dx = pw / (n - 1);
    function X(i) { return TR.l + i * dx; }
    function Y(v) { return TR.t + ph * (1 - v / TR.max); }
    var names = ['B', 'T', 'Myeloid'], cls = ['lb', 'lt', 'lm'];
    var g = '';
    [0, 0.2, 0.4, 0.6].forEach(function (v) {
      g += '<line class="gl" x1="' + TR.l + '" x2="' + (TR.l + pw) + '" y1="' + Y(v).toFixed(1) + '" y2="' + Y(v).toFixed(1) + '"/>' +
        '<text class="tk" x="' + (TR.l - 5) + '" y="' + (Y(v) + 3.5).toFixed(1) + '" text-anchor="end">' + (v === 0 ? '0' : v.toFixed(1)) + '</text>';
    });
    g += '<line class="thr" x1="' + TR.l + '" x2="' + (TR.l + pw) + '" y1="' + Y(E.threshold).toFixed(1) + '" y2="' + Y(E.threshold).toFixed(1) + '"/>' +
      '<text class="tl" x="' + (TR.l + 3) + '" y="' + (Y(E.threshold) - 4).toFixed(1) + '">bar ' + E.threshold.toFixed(2) + '</text>';
    // which markers enter when (panel order = event time), then the 4 extra steps
    var phases = [[0, 2, 'B'], [3, 5, 'T'], [6, 8, 'myeloid'], [9, n - 1, '+4 steps']];
    g += '<line class="dv" x1="' + (X(8) + dx / 2).toFixed(1) + '" x2="' + (X(8) + dx / 2).toFixed(1) + '" y1="' + TR.t + '" y2="' + (TR.H - TR.b + 8) + '"/>';
    phases.forEach(function (p) {
      var x1 = X(p[0]) - (p[0] ? dx / 2 : 0) + 1.5, x2 = X(p[1]) + (p[1] < n - 1 ? dx / 2 : 0) - 1.5, y = TR.H - TR.b + 8;
      g += '<line class="ph" x1="' + x1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y1="' + y + '" y2="' + y + '"/>' +
        '<text class="pt" x="' + ((x1 + x2) / 2).toFixed(1) + '" y="' + (y + 12) + '" text-anchor="middle">' + p[2] + '</text>';
    });
    [2, 1, 0].forEach(function (k) {
      var d = E.levels.map(function (row, i) { return (i ? 'L' : 'M') + X(i).toFixed(1) + ' ' + Y(row[k]).toFixed(1); }).join('');
      var last = E.levels[n - 1][k], ly = Y(last);
      g += '<path class="ln ' + cls[k] + '" d="' + d + '"/><circle class="dt ' + cls[k] + '" cx="' + X(n - 1).toFixed(1) + '" cy="' + ly.toFixed(1) + '" r="3.5"/>' +
        '<text class="el" x="' + (X(n - 1) + 7) + '" y="' + (ly - 2).toFixed(1) + '">' + names[k] + '</text>' +
        '<text class="ev" x="' + (X(n - 1) + 7) + '" y="' + (ly + 11).toFixed(1) + '">' + last.toFixed(3) + '</text>';
    });
    g += '<line class="xh" x1="0" x2="0" y1="' + TR.t + '" y2="' + (TR.t + ph) + '" style="display:none"/>' +
      '<rect class="hit" x="' + (TR.l - dx / 2) + '" y="0" width="' + (pw + dx) + '" height="' + (TR.t + ph + 4) + '"/>';
    return '<svg class="g-trace" viewBox="0 0 ' + TR.W + ' ' + TR.H + '" role="img" aria-label="Lineage scores of cell ' + E.cell + ' over ' + n +
      ' steps: myeloid ends at ' + E.levels[n - 1][2].toFixed(3) + ', T at ' + E.levels[n - 1][1].toFixed(3) + ', B at ' + E.levels[n - 1][0].toFixed(3) + '.">' + g + '</svg>';
  }
  function traceTable(E) {
    return '<details class="g-tbl"><summary>Numbers</summary><table><thead><tr><th>Step</th><th>B</th><th>T</th><th>Myeloid</th></tr></thead><tbody>' +
      E.levels.map(function (r, i) { return '<tr><td>' + i + '</td><td>' + r[0].toFixed(3) + '</td><td>' + r[1].toFixed(3) + '</td><td>' + r[2].toFixed(3) + '</td></tr>'; }).join('') +
      '</tbody></table></details>';
  }
  function wireTrace(root, E) {
    var svg = root.querySelector('.g-trace'), out = root.querySelector('.g-trace-read');
    if (!svg || !out) return;
    var xh = svg.querySelector('.xh'), hit = svg.querySelector('.hit'), n = E.levels.length;
    var pw = TR.W - TR.l - TR.r, dx = pw / (n - 1), idle = out.innerHTML;
    function at(e) {
      var r = svg.getBoundingClientRect(), x = (e.clientX - r.left) * TR.W / r.width;
      var i = Math.max(0, Math.min(n - 1, Math.round((x - TR.l) / dx))), row = E.levels[i];
      xh.setAttribute('x1', TR.l + i * dx); xh.setAttribute('x2', TR.l + i * dx); xh.style.display = '';
      out.innerHTML = '<b>Step ' + i + '</b> · myeloid ' + row[2].toFixed(3) + ' · T ' + row[1].toFixed(3) + ' · B ' + row[0].toFixed(3);
    }
    function leave() { xh.style.display = 'none'; out.innerHTML = idle; }
    hit.addEventListener('pointermove', at);
    hit.addEventListener('pointerdown', at);
    hit.addEventListener('pointerleave', leave);
  }
  function engineHtml(E) {
    return '<details class="g-engine" open><summary><span>' + E.title + '</span><em>' + E.sub + '</em></summary><div class="g-eng">' +
      stepList(E.steps) +
      '<div class="g-ex-t">Cell <code>' + E.cell + '</code> · soft rule</div>' +
      '<p class="g-mini">① Evidence in (scaled 0–1)</p>' + inputGrid(E.inputs) +
      '<p class="g-mini">② The field fills up</p>' + traceSvg(E) +
      '<p class="g-trace-read">Hover or tap a step to read its levels.</p>' +
      '<p class="g-ex-cap">Markers enter one per step (B’s three, then T’s, then myeloid’s), then the field runs 4 more steps. Myeloid overtakes T once CD16 arrives. This order was a hidden weight, not a real processing order: the nine predictions come out of one prediction. The corrected rerun enters all markers at once.</p>' + traceTable(E) +
      '<p class="g-mini">③ Read out, ④ check</p>' + exampleBars(E.ex) +
      '</div></details>';
  }
  // workflow map: one row per group (TEDDY's steps, then ANM's), numbered stations on a line
  function mapHtml(cfg) {
    var groups = cfg.groups && cfg.groups.length ? cfg.groups : [{ before: cfg.steps[0].id, cls: '', label: '' }], idx = {};
    cfg.steps.forEach(function (st, i) { idx[st.id] = i; });
    var rows = groups.map(function (g, k) {
      var to = k + 1 < groups.length ? idx[groups[k + 1].before] : cfg.steps.length;
      return { g: g, steps: cfg.steps.slice(idx[g.before], to) };
    });
    var cols = Math.max.apply(null, rows.map(function (r) { return r.steps.length; }));
    function node(st) {
      return '<a class="g-node ' + (st.cls || '') + '" data-step="' + st.id + '" href="' + (st.href || '#') + '" title="' + esc(st.t + ': ' + st.d) + '">' +
        '<span class="n">' + st.n + '</span><span class="s">' + (st.s || st.t) + '</span></a>';
    }
    return '<div class="g-map" style="--n:' + cols + '">' + rows.map(function (r) {
      return '<div class="g-map-grp">' + (r.g.label ? '<span class="g-grp ' + r.g.cls + '">' + r.g.label + '</span>' : '') +
        '<div class="g-map-row" style="--k:' + r.steps.length + '">' + r.steps.map(node).join('') + '</div></div>';
    }).join('') + '</div>';
  }

  // The phase × mode plane: phases build the evidence (x), modes say how ANM uses TEDDY (y).
  // Origin (bottom-left) = phase 1 × Mode B, where this repo's results are.
  var PLANE = {
    x: [['phase1', 'Phase 1', 'TEDDY + small head'], ['phase2', 'Phase 2', 'fusion model']],
    y: [['modeb', 'Mode B', 'on TEDDY + head outputs'], ['modea', 'Mode A', 'inside TEDDY']],
    cells: {
      B1: ['done', 'Done', 'v2: Experiments 1, 3, 4, 5; v3: E1, E2, E4, E6 final'],
      B2: ['part', 'Partial', 'the unfinished phase 2 is only the RNA-missing stand-in (withdrawn Experiment 2)'],
      A1: ['part', 'Started', 'v2 probes; v3: E5, E5-M, E7 final'],
      A2: ['no', 'Not done', 'inside the fusion model']
    }
  };
  function planeHtml(big) {
    function term(t) { return '<span class="term" data-term="' + t[0] + '">' + t[1] + '</span>'; }
    function cell(id) {
      var c = PLANE.cells[id];
      return '<div class="gp-cell ' + c[0] + '" data-cell="' + id + '"><span class="gp-st">' + c[1] + '</span><span class="gp-tx">' + c[2] + '</span><i class="gp-pin" aria-hidden="true"></i></div>';
    }
    var row = function (yi, key) { var y = PLANE.y[yi]; return '<div class="gp-rl">' + term(y) + '<small>' + y[2] + '</small></div>' + cell(key + '1') + cell(key + '2'); };
    return '<div class="g-plane' + (big ? ' big' : '') + '" role="group" aria-label="Phase × mode plane: which part of the project each section belongs to">' +
      '<div class="gp-grid">' + row(1, 'A') + row(0, 'B') +
      '<div class="gp-corner"></div>' + PLANE.x.map(function (x) { return '<div class="gp-cl">' + term(x) + '<small>' + x[2] + '</small></div>'; }).join('') + '</div>' +
      '<div class="gp-axes"><span>Mode →</span><span>Phase →</span></div></div>'  // the Mode label is rotated, so its → reads as ↑;
  }

  // ---------- left timeline: page progress, one marker per section, a short note for where you are ----------
  function shortOf(sec) {
    if (sec.short) return sec.short;
    var m = /^(?:Experiment|Block) (\d)/.exec(sec.label) || /^(\d\d)\b/.exec(sec.label);
    return m ? m[1] : (sec.label === 'Overview' ? '↑' : sec.label.charAt(0));
  }
  // links to other pages (cfg.railLinks, e.g. the ANM framework page) sit at the foot of the rail, below the track
  var FW_ICON = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><circle cx="6" cy="7" r="2.5"/><circle cx="18" cy="7" r="2.5"/><circle cx="12" cy="18" r="2.5"/><path d="M8.2 8.8l2.7 6.9M15.8 8.8l-2.7 6.9M8.6 7h6.8"/></svg>';
  function injectCss() {
    if (document.getElementById('t-ext-css')) return;
    var st = document.createElement('style');
    st.id = 't-ext-css';
    st.textContent = '@media (min-width: 1200px) {' +
      '.t-exts { position: absolute; left: 0; right: 0; bottom: 12px; display: grid; justify-items: center; gap: 6px; }' +
      '.t-ext { width: 32px; height: 32px; border-radius: 10px; display: grid; place-items: center; color: var(--accent); background: var(--bg); border: 1.5px solid color-mix(in srgb, var(--accent) 60%, var(--line)); text-decoration: none; transition: background .2s ease, color .2s ease, transform .15s ease; }' +
      '.t-ext:hover, .t-ext:focus-visible { background: var(--accent); color: var(--bg); transform: scale(1.06); outline: none; } }' +
      '.g-fw { display: block; margin-top: 16px; padding: 10px 12px; border-radius: 12px; border: 1px solid color-mix(in srgb, var(--accent) 35%, var(--line)); background: var(--accent-soft); color: var(--ink-muted); text-decoration: none; font-size: .78rem; line-height: 1.45; }' +
      '.g-fw b { display: block; color: var(--accent); font-size: .84rem; margin-bottom: 2px; }' +
      '.g-fw:hover { border-color: var(--accent); }' +
      '@media (prefers-reduced-motion: reduce) { .t-ext { transition: none !important; } }';
    document.head.appendChild(st);
  }
  function buildRail(cfg, ids) {
    var order = ids.slice().sort(function (a, b) { return document.getElementById(a).offsetTop - document.getElementById(b).offsetTop; });
    var ext = cfg.railLinks || [], exH = ext.length ? 46 + 38 * (ext.length - 1) : 0;
    if (ext.length) injectCss();
    var rail = document.createElement('nav');
    rail.className = 't-rail'; rail.setAttribute('aria-label', 'Page timeline');
    rail.innerHTML = '<div class="t-track"' + (exH ? ' style="bottom:' + (18 + exH) + 'px"' : '') + '><div class="t-fill"></div></div>' + order.map(function (id) {
      var sec = cfg.sections[id];
      return '<a class="t-dot" href="#' + id + '" data-id="' + id + '" aria-label="' + esc(sec.label) + '"><span>' + shortOf(sec) + '</span></a>';
    }).join('') + (ext.length ? '<div class="t-exts">' + ext.map(function (l, i) {
      return '<a class="t-ext" href="' + esc(l.href) + '" data-ext="' + i + '" aria-label="' + esc(l.label + (l.plain ? ': ' + l.plain : '')) + '">' + FW_ICON + '</a>';
    }).join('') + '</div>' : '') + '<div class="t-card" aria-hidden="true"><b></b><span></span></div>';
    var bar = document.createElement('div');
    bar.className = 't-bar'; bar.setAttribute('aria-hidden', 'true'); bar.innerHTML = '<i></i>';
    document.body.appendChild(rail); document.body.appendChild(bar);
    var fill = rail.querySelector('.t-fill'), barFill = bar.querySelector('i'), card = rail.querySelector('.t-card');
    var dots = [].slice.call(rail.querySelectorAll('.t-dot')), tops = {}, active = null, hideT = null, hovering = null;
    function maxScroll() { return Math.max(1, document.documentElement.scrollHeight - window.innerHeight); }
    function layout() {
      var h = rail.clientHeight - 36 - exH, prev = -1e9, pos = [];
      order.forEach(function (id) {
        var el = document.getElementById(id);
        var f = Math.max(0, Math.min(1, (el.getBoundingClientRect().top + window.scrollY - window.innerHeight * 0.3) / maxScroll()));
        var y = Math.max(f * h, prev + 22); pos.push(y); prev = y;
      });
      var over = pos.length ? pos[pos.length - 1] - h : 0;
      if (over > 0) pos = pos.map(function (y, i) { return y * h / (h + over); });
      dots.forEach(function (d, i) { d.style.top = (18 + pos[i]) + 'px'; tops[d.getAttribute('data-id')] = 18 + pos[i]; });
      progress();
    }
    function progress() {
      var f = Math.max(0, Math.min(1, window.scrollY / maxScroll()));
      fill.style.height = (f * 100) + '%'; barFill.style.width = (f * 100) + '%';
    }
    function showText(title, text, at, sticky) {
      card.querySelector('b').textContent = title; card.querySelector('span').textContent = text || '';
      var y = Math.max(8, Math.min(rail.clientHeight - 90, at - 14));
      card.style.top = y + 'px'; card.classList.add('on');
      clearTimeout(hideT);
      if (!sticky) hideT = setTimeout(function () { if (!hovering) card.classList.remove('on'); }, 1600);
    }
    function showCard(id, sticky) {
      var sec = cfg.sections[id]; if (!sec) return;
      showText(sec.label, sec.plain, tops[id] || 0, sticky);
    }
    [].slice.call(rail.querySelectorAll('.t-ext')).forEach(function (a) {
      var l = ext[+a.getAttribute('data-ext')];
      function on() { hovering = 'ext'; showText(l.label, l.plain, a.getBoundingClientRect().top - rail.getBoundingClientRect().top - 64, true); }
      a.addEventListener('pointerenter', on);
      a.addEventListener('focus', on);
      a.addEventListener('pointerleave', function () { hovering = null; if (active) showCard(active, false); else card.classList.remove('on'); });
      a.addEventListener('blur', function () { hovering = null; card.classList.remove('on'); });
    });
    dots.forEach(function (d) {
      var id = d.getAttribute('data-id');
      d.addEventListener('pointerenter', function () { hovering = id; showCard(id, true); });
      d.addEventListener('pointerleave', function () { hovering = null; if (active) showCard(active, false); else card.classList.remove('on'); });
      d.addEventListener('focus', function () { showCard(id, true); });
      d.addEventListener('blur', function () { card.classList.remove('on'); });
    });
    var scrolled = false;
    window.addEventListener('scroll', function () {
      progress();
      if (!scrolled) { scrolled = true; requestAnimationFrame(function () { scrolled = false; if (active && !hovering) showCard(active, false); }); }
    }, { passive: true });
    window.addEventListener('resize', layout);
    if ('ResizeObserver' in window) new ResizeObserver(function () { requestAnimationFrame(layout); }).observe(document.querySelector('main') || document.body);
    window.addEventListener('load', layout);
    layout();
    return {
      setActive: function (id) {
        active = id;
        var passed = true;
        dots.forEach(function (d) {
          var on = d.getAttribute('data-id') === id;
          d.classList.toggle('on', on); d.classList.toggle('past', passed && !on);
          if (on) { d.setAttribute('aria-current', 'location'); passed = false; } else d.removeAttribute('aria-current');
        });
        // after a jump the scroll event can fire before the new section is known: keep the card current
        if (!hovering && card.classList.contains('on')) showCard(id, false);
      }
    };
  }

  function buildGuide(cfg) {
    var aside = document.getElementById('guide');
    if (!aside || !cfg) return;
    if (cfg.railLinks && cfg.railLinks.length) injectCss();
    var E = cfg.engine === false ? null : (cfg.engine || ENGINE);
    var arms = cfg.arms ? '<p class="g-kick">Colours on this page</p><div class="g-arms">' + cfg.arms.map(function (a) { return '<span><i style="background:' + a[0] + '"></i>' + a[1] + '</span>'; }).join('') + '</div>' : '';
    aside.innerHTML = '<div class="g-inner">' +
      '<div class="g-top"><span class="g-kick">Guide · how ANM decides</span><span class="g-btns">' +
      '<button class="g-hide" type="button" aria-label="Hide guide" title="Hide the guide (drag its left edge to resize)"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true"><path d="M9 6l6 6-6 6"/></svg></button>' +
      '<button class="g-close" type="button" aria-label="Close guide">×</button></span></div>' +
      '<div class="g-now" aria-live="polite"><span class="g-now-k">You are reading</span><b id="g-now-t"></b><span id="g-now-d"></span></div>' +
      '<p class="g-kick">Where it sits · phase × mode</p>' + planeHtml(false) +
      '<p class="g-kick">' + (cfg.flowTitle || 'Where this sits in the workflow') + '</p>' + mapHtml(cfg) +
      '<p class="g-kick" id="g-how-k">How ANM decides here</p><div class="g-how" id="g-how"></div>' +
      (E ? engineHtml(E) : '') +
      '<div id="g-terms"><p class="g-kick">Terms in this section</p><div class="g-concepts" id="g-concepts"></div></div>' +
      (cfg.allTerms ? '<a class="g-all" href="' + cfg.allTerms + '">' + (cfg.allTermsLabel || 'All terms') + ' →</a>' : '') +
      (cfg.railLinks || []).map(function (l) { return '<a class="g-fw" href="' + esc(l.href) + '"><b>' + l.label + ' →</b>' + (l.plain || '') + '</a>'; }).join('') + arms +
      '</div>';
    decorateAll(aside);
    if (E) wireTrace(aside, E);

    // a page can place its own opener (e.g. in the top bar); otherwise a floating button is added
    var fab = document.querySelector('.g-open');
    if (!fab) {
      fab = document.createElement('button');
      fab.className = 'g-fab'; fab.type = 'button'; fab.setAttribute('aria-controls', 'guide'); fab.setAttribute('aria-expanded', 'false');
      fab.innerHTML = '<svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true"><circle cx="12" cy="12" r="9"/><path d="M12 11v6M12 7.5v.5"/></svg>Guide';
      document.body.appendChild(fab);
    }
    var scrim = document.createElement('div'); scrim.className = 'g-scrim';
    document.body.appendChild(scrim);
    function open() {
      if (wide()) { setCollapsed(false); return; }
      aside.classList.add('open'); scrim.classList.add('on'); fab.setAttribute('aria-expanded', 'true');
    }
    function close() { aside.classList.remove('open'); scrim.classList.remove('on'); fab.setAttribute('aria-expanded', 'false'); }
    fab.addEventListener('click', open);

    // wide screens: the column floats on the right; it can be hidden to a tab on the edge and resized from its left edge.
    // Both choices are remembered per browser (a convenience only; the page works without storage).
    var root = document.documentElement, W_MIN = 260, W_DEF = 320;
    function wide() { return window.innerWidth >= 1200; }
    function wMax() { return Math.max(W_MIN, Math.min(640, window.innerWidth - 640)); }
    function load(k) { try { return window.localStorage.getItem(k); } catch (e) { return null; } }
    function save(k, v) { try { window.localStorage.setItem(k, v); } catch (e) {} }
    var relayoutT = null;
    function relayout(delay) { clearTimeout(relayoutT); relayoutT = setTimeout(function () { window.dispatchEvent(new Event('resize')); }, delay || 0); }
    var handle = document.createElement('div');
    handle.className = 'g-resize'; handle.setAttribute('role', 'separator'); handle.setAttribute('aria-orientation', 'vertical');
    handle.setAttribute('aria-label', 'Resize the guide'); handle.setAttribute('aria-controls', 'guide'); handle.tabIndex = 0;
    handle.title = 'Drag to resize · double-click to reset';
    var tab = document.createElement('button');
    tab.className = 'g-tab'; tab.type = 'button'; tab.setAttribute('aria-controls', 'guide'); tab.setAttribute('aria-label', 'Show the guide');
    tab.innerHTML = '<svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true"><path d="M15 6l-6 6 6 6"/></svg><span>Guide</span>';
    document.body.appendChild(handle); document.body.appendChild(tab);
    var width = W_DEF;
    function setWidth(w, persist) {
      width = Math.round(Math.max(W_MIN, Math.min(wMax(), w)));
      root.style.setProperty('--gw', width + 'px');
      handle.setAttribute('aria-valuemin', W_MIN); handle.setAttribute('aria-valuemax', wMax()); handle.setAttribute('aria-valuenow', width);
      if (persist) save('teddyGuide.width', String(width));
    }
    function setCollapsed(c) {
      root.classList.toggle('g-collapsed', c);
      tab.setAttribute('aria-expanded', String(!c));
      save('teddyGuide.collapsed', c ? '1' : '0');
      relayout(260);
    }
    setWidth(parseInt(load('teddyGuide.width'), 10) || W_DEF);
    if (load('teddyGuide.collapsed') === '1') root.classList.add('g-collapsed');
    // a non-default layout: let the charts re-measure once the column has its final size
    if (width !== W_DEF || root.classList.contains('g-collapsed')) relayout(0);
    aside.querySelector('.g-hide').addEventListener('click', function () { setCollapsed(true); tab.focus(); });
    tab.addEventListener('click', function () { setCollapsed(false); });
    var dragging = false, lastEmit = 0;
    handle.addEventListener('pointerdown', function (e) {
      if (!wide()) return;
      dragging = true; handle.setPointerCapture(e.pointerId); root.classList.add('g-dragging'); e.preventDefault();
    });
    handle.addEventListener('pointermove', function (e) {
      if (!dragging) return;
      setWidth(window.innerWidth - e.clientX);
      var now = e.timeStamp || 0;
      if (now - lastEmit > 120) { lastEmit = now; relayout(0); }
    });
    function endDrag() { if (!dragging) return; dragging = false; root.classList.remove('g-dragging'); setWidth(width, true); relayout(0); }
    handle.addEventListener('pointerup', endDrag);
    handle.addEventListener('pointercancel', endDrag);
    handle.addEventListener('dblclick', function () { setWidth(W_DEF, true); relayout(0); });
    handle.addEventListener('keydown', function (e) {
      var d = { ArrowLeft: 20, ArrowRight: -20 }[e.key];
      if (d) { setWidth(width + d, true); relayout(150); e.preventDefault(); }
      else if (e.key === 'Home') { setWidth(W_MIN, true); relayout(0); e.preventDefault(); }
      else if (e.key === 'End') { setWidth(wMax(), true); relayout(0); e.preventDefault(); }
    });
    window.addEventListener('resize', function (e) { if (e.isTrusted && width > wMax()) setWidth(width); });
    scrim.addEventListener('click', close);
    aside.querySelector('.g-close').addEventListener('click', close);
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
    aside.addEventListener('click', function (e) { var a = e.target.closest && e.target.closest('a.g-node, a.g-all, a.g-fw'); if (a && window.innerWidth < 1200) close(); });

    var nodes = {};
    [].slice.call(aside.querySelectorAll('.g-node')).forEach(function (a) { nodes[a.getAttribute('data-step')] = a; });
    var nowT = aside.querySelector('#g-now-t'), nowD = aside.querySelector('#g-now-d');
    var howK = aside.querySelector('#g-how-k'), howEl = aside.querySelector('#g-how');
    var terms = aside.querySelector('#g-terms'), conEl = aside.querySelector('#g-concepts');
    var visible = {}, currentId = null, currentKey = null;
    function setNow(id) {
      var key = cfg.resolve ? cfg.resolve(id) : id;
      var sec = cfg.sections[key] || cfg.sections[id];
      if (!sec) return;
      currentId = id;
      if (rail) rail.setActive(id);
      if (key === currentKey) return;
      currentKey = key;
      nowT.textContent = sec.label; nowD.textContent = sec.plain;
      var here = sec.plane || (cfg.planeBySection && cfg.planeBySection[key]) || cfg.planeDefault || ['B1'];
      [].slice.call(aside.querySelectorAll('.gp-cell')).forEach(function (c) { c.classList.toggle('on', here.indexOf(c.getAttribute('data-cell')) > -1); });
      Object.keys(nodes).forEach(function (k) {
        var on = (sec.steps || []).indexOf(k) > -1;
        nodes[k].classList.toggle('on', on);
        if (on) nodes[k].setAttribute('aria-current', 'step'); else nodes[k].removeAttribute('aria-current');
      });
      var how = sec.how;
      howK.textContent = (how && how.k) || 'How ANM decides here';
      howK.hidden = howEl.hidden = !how;
      howEl.innerHTML = how ? '<div class="g-how-t">' + how.t + '</div>' + stepList(how.steps) + exampleBars(how.ex) + (how.note ? '<p class="g-how-note">' + how.note + '</p>' : '') : '';
      var ids = sec.concepts || cfg.concepts || [];
      conEl.innerHTML = ids.map(function (tid) {
        var t = TERMS[tid];
        return t ? '<details><summary>' + t[0] + '</summary><p>' + t[1] + '</p></details>' : '';
      }).join('');
      terms.hidden = !ids.length;
      decorateAll(howEl);
      howEl.classList.remove('fade'); void howEl.offsetWidth; howEl.classList.add('fade');
    }
    window.GUIDE_REFRESH = function () { if (currentId) { currentKey = null; setNow(currentId); } };
    var ids = Object.keys(cfg.sections).filter(function (id) { return document.getElementById(id); });
    var rail = ids.length ? buildRail(cfg, ids) : null;
    // at the very bottom the last section may never reach the trigger band, so it wins there
    function atBottom() { return window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 2; }
    function pick() {
      if (atBottom() && ids.length) { setNow(ids[ids.length - 1]); return; }
      var best = null, bestDepth = -1;
      ids.forEach(function (id) {
        if (!visible[id]) return;
        var d = cfg.sections[id].depth || 1;
        if (d >= bestDepth) { best = id; bestDepth = d; }
      });
      if (best) setNow(best);
    }
    if ('IntersectionObserver' in window) {
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (en) { visible[en.target.id] = en.isIntersecting; });
        pick();
      }, { rootMargin: '-30% 0px -60% 0px' });
      ids.forEach(function (id) { io.observe(document.getElementById(id)); });
      var wasBottom = false;
      window.addEventListener('scroll', function () { var bt = atBottom(); if (bt !== wasBottom) { wasBottom = bt; pick(); } }, { passive: true });
    }
    setNow(ids[0]);
  }

  // ---------- init ----------
  function init() {
    var main = document.querySelector('main') || document.body;
    var pm = document.getElementById('plane-main');
    if (pm) { pm.innerHTML = planeHtml(true); decorateAll(pm); }
    decorateAll(document);
    autowrap(main);
    buildGuide(window.GUIDE_CONFIG);
    // captions, tables and panels are re-rendered by the page; mark their terms too
    if ('MutationObserver' in window) {
      var pending = false;
      new MutationObserver(function () {
        if (pending) return;
        pending = true;
        requestAnimationFrame(function () { pending = false; decorateAll(main); autowrap(main); });
      }).observe(main, { childList: true, subtree: true, characterData: true });
    }
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
