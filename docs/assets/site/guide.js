/* Shared by docs/index.html and docs/anm-loop.html.
 * 1) TERMS: one plain-language definition per technical term.
 * 2) Terms in prose get a dotted underline and show their definition on hover, focus or tap.
 *    Known abbreviations (O0, P_f, z_512, …) are marked automatically; <span class="term" data-term="…"> marks anything else.
 * 3) The guide column: the whole workflow, what you are reading now, and the key concepts.
 *    Each page sets window.GUIDE_CONFIG before loading this file. */
(function () {
  'use strict';

  var TERMS = {
    teddy: ['TEDDY', 'Merck’s single-cell foundation model (TEDDY-G, 70M parameters). Here it reads a cell’s RNA and, with a small head, predicts 9 surface proteins. It is never retrained.'],
    anm: ['ANM · Active Neural Matter', 'A decision layer. It reads evidence, applies a written-down question (the observer) and returns a call or an honest “no call”, with checks attached.'],
    citeseq: ['CITE-seq', 'A technology that measures RNA and surface proteins in the same cell. The measured proteins are the answer key we hold out.'],
    adt: ['ADT · measured surface protein', 'CITE-seq’s protein readout (antibody-derived tags). Held out and used only to check answers, except where a test says otherwise.'],
    evidence: ['Typed evidence (δu)', 'TEDDY’s predicted value for each of the 9 panel proteins, handed to ANM as labelled inputs such as “CD19 is high”.'],
    observer: ['Observer · the question', 'The question, written down: which markers count, how strong a signal must be, and when to decline. Changing it needs no labels and no retraining.'],
    o0: ['Soft rule (O0)', 'The default lineage call. All 9 markers weigh equally; call the best lineage if its score reaches 0.12.'],
    o1: ['Strict rule (O1)', 'The same question asked more strictly: key markers CD19, CD3 and CD16 count double and the score must reach 0.28. Same answers where both decide; declines more often.'],
    o2: ['Key-marker rule (O2)', 'A different question: only CD19 (B), CD3 (T) and CD16 (myeloid) count, with B and T weighted up. It changes 12.48% of the expected answers.'],
    abstain: ['Decline to call (abstain)', 'The readout says “no call” because the evidence does not meet the written-down rule: honest silence instead of a guess.'],
    pf: ['Workability (P_f)', '1 if the observer can make a call from the evidence, 0 if it declines. Averaged over cells, it is the share of cells called.'],
    qf: ['Exactness (Q_f)', '1 if a call matches the held-out protein truth. Averaged over labelled cells, a decline counts as a miss.'],
    qdec: ['Accuracy of calls made (Q)', 'Correct calls divided by calls made, checked against the held-out proteins. Used to compare arms like for like.'],
    softp: ['Confidence score (soft_P)', 'How strong ANM’s best option is for a cell. Keeping only high-score cells trades coverage for accuracy.'],
    coverage: ['Coverage', 'The share of cells that still get a call once a confidence cut is applied.'],
    loo: ['Leave-one-out (LOO)', 'Remove one protein’s evidence and decide again. The marker whose removal moves the score most is the reason for the call.'],
    flip: ['Flip distance', 'How much one marker’s value must change before the call switches lineage. Small means a fragile call.'],
    holdout: ['Held-out cells (site4)', '16,750 cells from a site never used in training. Their measured proteins only score the answers.'],
    ood: ['Out-of-site check', '2,000 cells from other sites (val_non_site4): a second check that the results hold elsewhere.'],
    train: ['Train + Jev · trained head', 'The “just train harder” control: a classifier (logistic, MLP or threshold grid) fitted on labels over the same TEDDY features. A stand-in for a Jev-class learner, not the live Jev API.'],
    z512: ['TEDDY embedding (z_512)', 'TEDDY’s 512-number summary of a cell: the mean of its last-layer tokens at context length 1024.'],
    mustpair: ['Must-separate pair', 'Two cells TEDDY sees as near-identical (cosine ≥ 0.98 in its embedding) whose measured proteins disagree, for example one myeloid and one T cell.'],
    falseagree: ['False agreement', 'Among must-separate pairs where both cells get a call, the share given the same call even though they differ.'],
    softsep: ['Soft separation', 'The share of must-separate pairs the readout tells apart: different calls, or one cell declined.'],
    corrsep: ['Correct separation', 'The share of must-separate pairs split into the right, different calls.'],
    zplus: ['z⁺ · feature swap', 'TEDDY’s embedding with protein features glued on: the “just swap features” alternative to typed evidence.'],
    modeb: ['Mode B', 'What this project does: treat TEDDY as an evidence source and study the decision layer on top of it.'],
    modea: ['Mode A (not done)', 'Studying TEDDY’s internals, such as its residual stream, layer Jacobians or in-silico gene perturbations. Not claimed here.'],
    rnaonly: ['RNA only', 'Evidence from the RNA-based predictions; no protein channel.'],
    adtonly: ['RNA missing (protein only)', 'Evidence only from the protein channel of the phase-2 model, a weak channel (panel Pearson 0.446).'],
    joint: ['Both channels', 'Evidence from the RNA and protein channels together.'],
    lineage: ['Lineage call', 'The decision per cell: B cell, T cell or myeloid, each judged from 3 markers (B: CD19 CD72 CD22 · T: CD3 CD2 CD5 · myeloid: CD16 CD11c CD36).'],
    pearson: ['Pearson ≈ 0.61', 'How well TEDDY’s phase-1 head predicts all measured proteins (correlation on held-out cells). A reference point; we do not try to beat it.'],
    veto: ['Veto', 'ANM flags that TEDDY’s embedding cannot tell a pair apart, so a readout built on it should not be trusted there.'],
    complement: ['Complement', 'Add a missing evidence channel (here, measured protein as typed evidence) instead of retraining TEDDY.'],
    verify: ['Verify', 'Score the same pairs again after the change: false agreement, soft separation and accuracy.'],
    finitefield: ['Finite field · ANM’s engine', 'ANM’s decision engine: it combines the typed evidence under the observer into a score per lineage, then calls or declines.'],
    permutation: ['Permutation test', 'Shuffle the data many times to see how often a pattern this strong appears by chance. Here p ≈ 0.0099 over 100 shuffles.'],
    bootstrap: ['Bootstrap', 'Resample the cells with replacement (200 times) to see how stable a number is; gives a 95% interval.']
  };
  window.TEDDY_ANM_TERMS = TERMS;

  // abbreviations marked automatically in prose (headings, buttons, links and code are skipped)
  var AUTO = [
    ['O0', 'o0'], ['O1', 'o1'], ['O2', 'o2'], ['soft_P', 'softp'], ['z_512', 'z512'], ['z⁺', 'zplus'], ['LOO', 'loo'], ['δu', 'evidence'],
    ['Mode B', 'modeb'], ['Mode A', 'modea'], ['CITE-seq', 'citeseq'], ['ADT', 'adt'], ['P_f', 'pf'], ['Q_f', 'qf'], ['Jev', 'train'],
    ['site4', 'holdout'], ['adt_only', 'adtonly'], ['rna_only', 'rnaonly'],
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
  document.addEventListener('click', function (e) {
    var t = termOf(e);
    if (t) { if (cur === t) hide(); else show(t); return; }
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
  function buildGuide(cfg) {
    var aside = document.getElementById('guide');
    if (!aside || !cfg) return;
    var steps = cfg.steps, groups = cfg.groups || [];
    var flow = '';
    steps.forEach(function (s) {
      groups.forEach(function (g) { if (g.before === s.id) flow += '<li class="g-grp ' + g.cls + '">' + g.label + '</li>'; });
      flow += '<li><a class="g-step ' + (s.cls || '') + '" data-step="' + s.id + '" href="' + (s.href || '#') + '"><span class="n">' + s.n + '</span><span><b>' + s.t + '</b><span class="d">' + s.d + '</span></span></a></li>';
    });
    var concepts = cfg.concepts.map(function (id) {
      var t = TERMS[id];
      return t ? '<details><summary>' + t[0] + '</summary><p>' + t[1] + '</p></details>' : '';
    }).join('');
    var arms = cfg.arms ? '<p class="g-kick">Colours</p><div class="g-arms">' + cfg.arms.map(function (a) { return '<span><i style="background:' + a[0] + '"></i>' + a[1] + '</span>'; }).join('') + '</div>' : '';
    aside.innerHTML = '<div class="g-inner">' +
      '<div class="g-top"><span class="g-kick">Guide</span><button class="g-close" type="button" aria-label="Close guide">×</button></div>' +
      '<p class="g-kick">You are reading</p><div class="g-now" aria-live="polite"><b id="g-now-t"></b><span id="g-now-d"></span></div>' +
      '<p class="g-kick">' + (cfg.flowTitle || 'The whole workflow') + '</p><ol class="g-flow">' + flow + '</ol>' +
      '<p class="g-kick">Key concepts · click to expand</p><div class="g-concepts">' + concepts + '</div>' + arms +
      '</div>';

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
    function open() { aside.classList.add('open'); scrim.classList.add('on'); fab.setAttribute('aria-expanded', 'true'); }
    function close() { aside.classList.remove('open'); scrim.classList.remove('on'); fab.setAttribute('aria-expanded', 'false'); }
    fab.addEventListener('click', open);
    scrim.addEventListener('click', close);
    aside.querySelector('.g-close').addEventListener('click', close);
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
    [].slice.call(aside.querySelectorAll('.g-step')).forEach(function (a) { a.addEventListener('click', function () { if (window.innerWidth < 1200) close(); }); });

    var stepEls = {};
    [].slice.call(aside.querySelectorAll('.g-step')).forEach(function (a) { stepEls[a.getAttribute('data-step')] = a; });
    var nowT = aside.querySelector('#g-now-t'), nowD = aside.querySelector('#g-now-d');
    var visible = {}, currentId = null;
    function setNow(id) {
      var key = cfg.resolve ? cfg.resolve(id) : id;
      var sec = cfg.sections[key] || cfg.sections[id];
      if (!sec) return;
      currentId = id;
      nowT.textContent = sec.label; nowD.textContent = sec.plain;
      Object.keys(stepEls).forEach(function (k) { stepEls[k].classList.toggle('on', (sec.steps || []).indexOf(k) > -1); });
    }
    window.GUIDE_REFRESH = function () { if (currentId) setNow(currentId); };
    var ids = Object.keys(cfg.sections).filter(function (id) { return document.getElementById(id); });
    function pick() {
      var best = null, bestDepth = -1;
      ids.forEach(function (id) {
        if (!visible[id]) return;
        var d = cfg.sections[id].depth || 1;
        if (d >= bestDepth) { best = id; bestDepth = d; }
      });
      if (best && best !== currentId) setNow(best);
    }
    if ('IntersectionObserver' in window) {
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (en) { visible[en.target.id] = en.isIntersecting; });
        pick();
      }, { rootMargin: '-30% 0px -60% 0px' });
      ids.forEach(function (id) { io.observe(document.getElementById(id)); });
    }
    setNow(ids[0]);
  }

  // ---------- init ----------
  function init() {
    var main = document.querySelector('main') || document.body;
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
