/*
 * idx-exp-labels-cost-curve — index.html, #labels (Experiment 5 · zero-label transfer, "Cost of a new question").
 * How the test answers it (the matching idx-q-price-of-question scene shows only the problem).
 * 1 A wall splits the cells: 11,725 label-pool cells and 5,025 scored cells, no overlap (strip at the bottom;
 *   one dot ≈ 167.5 cells, so the two bins hold 70 and 30 dots).
 * 2 ANM writes the B/T-priority question down as three short lines (the rule as declared in the glossary:
 *   only CD19, CD3 and CD16 count, B ×1.5, T ×1.3, myeloid ×0.5). Its answers appear at once: a glowing
 *   "0 labels" tag with the page's numbers (declines 19.0%, exactness 0.835, accuracy 0.963), and a teal dot at
 *   0 on the label axis.
 * 3 Label coins fly up from the pool (never from the scored cells) into three trained heads: threshold grid,
 *   logistic head, MLP head, each tuned to decline as often as ANM. The bars are labels spent on the page's
 *   axis (0, then log scale from 50 to 10,619); the pool dims in proportion to the labels used. Each bar gets
 *   a teal ring where the page says it matches ANM's decline rate and exactness: ≈ 50, ≈ 500, ≈ 2,000.
 *   Only those budgets are measured for the grid and MLP, so no curve shape is drawn between them.
 * 4 The logistic head keeps taking labels up to 10,619 (accuracy 0.972 there, page value), but the dashed slot
 *   on the heads' side stays empty: no written rule.
 * Wide layout 21 × 8 when the stage is at least 2:1 (desktop 21:8), else narrow 16 × 9 (phone 16:9).
 * 12 s seamless loop; the still frame shows everything complete.
 */

/* Every TEDDY-pipeline number this scene shows; all appear in #labels on index.html (text, KPIs, b5 chart and
 * table). Refresh here after the official-preprocessing rerun. ruleWeights is the B/T-priority rule as declared
 * (glossary term o2 on the same page), a declaration rather than a result. */
const NUMBERS = {
  poolCells: 11725,        // label-pool cells (training side of the wall)
  scoredCells: 5025,       // scored cells, no overlap with the pool
  anmDeclines: 0.19,       // ANM, B/T-priority written down, 0 labels: declines 19.0%
  anmExactness: 0.835,     // exactness
  anmAccuracy: 0.963,      // accuracy of calls made
  minLabels: 50,           // smallest training budget
  maxLabels: 10619,        // largest training budget
  matchGrid: 50,           // labels for the re-tuned threshold grid to match ANM's decline rate and exactness
  matchLogistic: 500,      // logistic head
  matchMlp: 2000,          // MLP head
  logisticAccAtMax: 0.972, // logistic head accuracy of calls made at 10,619 labels
  ruleWeights: { B: 1.5, T: 1.3, myeloid: 0.5 },
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const N = NUMBERS;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 8 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const ln = (o = {}) => add(ctx.line([[0, 0], [1, 0]], o));
  const dot = (color, o = {}) => add(ctx.dot([0, 0, 0], Object.assign({ color }, o)));
  const int = (n) => Math.round(n).toLocaleString('en-US');
  const f3 = (v) => v.toFixed(3);
  const pct = (v) => (100 * v).toFixed(1) + '%';

  /* label axis: 0, a break, then log scale from the smallest to the largest budget */
  const QB = 0.1;
  const LOGSPAN = Math.log(N.maxLabels / N.minLabels);
  const qOf = (b) => (b <= N.minLabels ? (QB * b) / N.minLabels : QB + ((1 - QB) * Math.log(b / N.minLabels)) / LOGSPAN);
  const bOf = (q) => (q <= QB ? (N.minLabels * q) / QB : N.minLabels * Math.exp(((q - QB) / (1 - QB)) * LOGSPAN));

  /* timeline (fractions of the loop) */
  const T = { type0: 0.03, typeStep: 0.05, typeLen: 0.045, anm: 0.19, U0: 0.25, U1: 0.78, slot: 0.8, out0: 0.94, out1: 0.985 };
  const uOfQ = (q) => T.U0 + q * (T.U1 - T.U0);

  /* text measurement (CSS px), same fonts as the label sprites */
  const mc = document.createElement('canvas').getContext('2d');
  const fcs = getComputedStyle(ctx.figure);
  const FONT = fcs.fontFamily || 'system-ui, sans-serif';
  const MONO = (fcs.getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace';
  const tw = (text, size, o = {}) => {
    const lines = String(text).split('\n');
    if (!mc) return 4 + Math.max(...lines.map((l) => l.length)) * size * ctx.textScale * 0.6;
    mc.font = `${o.weight || 500} ${size * ctx.textScale}px ${o.mono ? MONO : FONT}`;
    return 4 + 1.05 * Math.max(...lines.map((l) => mc.measureText(l).width));
  };

  /* ── ANM side: the written question and its answers ── */
  const RULE = [`B: CD19 ×${N.ruleWeights.B}`, `T: CD3 ×${N.ruleWeights.T}`, `myeloid: CD16 ×${N.ruleWeights.myeloid}`];
  const hL = lab('ANM · written down', { anchor: 'left', color: 'accent', weight: 600, size: 13 });
  const card = add(ctx.box(1, 1, { color: 'soft', stroke: 'line', strokeWidth: 1.2, radius: 0.16 }));
  const cardT = lab('B/T-priority question', { anchor: 'left', color: 'muted', size: 11 });
  const rule = RULE.map((s) => ({ s, n: -1, l: lab('', { anchor: 'left', mono: true, size: 13 }) }));
  const glow = add(ctx.box(1, 1, { color: 'accent', radius: 0.3, opacity: 0 }));
  const pill = lab('0 labels', { anchor: 'left', bg: 'accent', bgOpacity: 1, color: 'card', weight: 600, size: 13, pad: 5, order: 12 });
  const tgt = lab('', { anchor: 'top-left', color: 'muted', size: 11 });

  /* ── heads side: header, axis, three rows, the empty rule slot ── */
  const HEADS = [
    { name: 'threshold grid', short: 'threshold grid', b: N.matchGrid },
    { name: 'logistic head', short: 'logistic', b: N.matchLogistic, more: true },
    { name: 'MLP head', short: 'MLP', b: N.matchMlp },
  ];
  const hR = lab('trained heads · labels to match ANM', { anchor: 'left', color: 'train', weight: 600, size: 13 });
  const sR = lab('each tuned to decline as often as ANM', { anchor: 'left', color: 'muted', size: 11 });
  const axis = [ln({ color: 'muted', width: 1.2 }), ln({ color: 'muted', width: 1.2 })];
  const brk = [ln({ color: 'muted', width: 1.2 }), ln({ color: 'muted', width: 1.2 })];
  const TICKS = [0, N.minLabels, N.matchLogistic, N.matchMlp, N.maxLabels];
  const ticks = TICKS.map((v, i) => ({
    v, t: ln({ color: 'muted', width: 1.2 }),
    l: lab(i ? int(v) : '0', { anchor: 'top', color: i ? 'muted' : 'accent', weight: i ? 500 : 600, size: 11 }),
  }));
  const axT = lab('labelled cells (log scale)', { anchor: 'top', color: 'muted', size: 11 });
  const anmHalo = dot('accent', { opacity: 0 });
  const anmDot = dot('accent', { opacity: 0 });
  const cursor = ln({ color: 'faint', width: 1.2, dashed: [3, 4] });
  const rows = HEADS.map((h) => Object.assign(h, {
    q: qOf(h.b), y: 0,
    nameL: lab(h.name, { anchor: 'right', size: 12 }),
    track: ln({ color: 'grid', width: 1 }),
    bar: ln({ color: 'train', width: 8, opacity: 0.85 }),
    back: dot('card', { opacity: 0 }),
    tip: dot('train', { hollow: true, ring: 0.36, opacity: 0 }),
    wave: dot('accent', { hollow: true, ring: 0.14, opacity: 0 }),
    mFill: dot('accent', { opacity: 0 }),
    mRing: dot('accent', { hollow: true, ring: 0.34, opacity: 0 }),
    mL: lab('≈ ' + int(h.b), { anchor: 'left', color: 'train', weight: 600, size: 12, bg: 'card', bgOpacity: 1, pad: 2 }),
  }));
  const LOG = rows[1];
  const more = {
    line: ln({ color: 'train', width: 2.2, dashed: [5, 4], opacity: 0.9 }),
    l: lab(`accuracy ${f3(N.logisticAccAtMax)} at ${int(N.maxLabels)}`, { anchor: 'bottom-right', color: 'train', weight: 600, size: 11 }),
  };
  const slot = {
    box: add(ctx.box(1, 1, { color: null, stroke: 'faint', strokeWidth: 1.3, dashed: [4, 4], radius: 0.16 })),
    lines: [0, 1, 2].map(() => ln({ color: 'faint', width: 1.3, dashed: [3, 4] })),
    l: lab('no written rule', { color: 'muted', weight: 600, size: 12 }),
  };

  /* ── the split: scored cells | wall | label pool (one dot ≈ 167.5 cells) ── */
  const NS = 30, NP = 70;
  const binS = add(ctx.box(1, 1, { color: 'soft', stroke: 'line', strokeWidth: 1, radius: 0.12 }));
  const binP = add(ctx.box(1, 1, { color: 'soft', stroke: 'line', strokeWidth: 1, radius: 0.12 }));
  const wall = ln({ color: 'muted', width: 3.5 });
  const wallL = lab('no overlap', { anchor: 'bottom', color: 'muted', size: 11 });
  const binSL = lab(`${int(N.scoredCells)} scored cells`, { anchor: 'top', color: 'muted', size: 11 });
  const binPL = lab('', { anchor: 'top', color: 'muted', size: 11 });
  const inBin = { bg: 'soft', bgOpacity: 1, pad: 2, color: 'muted', size: 10 }; // compact: labels sit inside the bins
  const binSLi = lab(`${int(N.scoredCells)} scored`, inBin);
  const binPLi = lab(`${int(N.poolCells)} label-pool cells · no overlap`, inBin);
  const dS = Array.from({ length: NS }, () => dot('muted', { opacity: 0.55 }));
  const dP = Array.from({ length: NP }, () => dot('muted', { opacity: 0.55 }));
  const rnd = ctx.rand(11);
  const order = dP.map((_, i) => i);
  for (let i = NP - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
  const rankOf = new Array(NP);
  order.forEach((k, r) => { rankOf[k] = r; });
  const usedAt = (u) => (u <= T.U0 ? 0 : (bOf(seg(u, T.U0, T.U1, ease.linear)) / N.poolCells) * NP);

  /* ── label coins: each head draws from the pool while it trains toward its match (the logistic to the end) ── */
  const DU = 0.012, FLY = 0.055;
  const coins = [];
  rows.forEach((h, i) => {
    const end = h.more ? T.U1 : uOfQ(h.q);
    for (let u0 = T.U0 - 0.8 * FLY + (i * DU) / 3; u0 <= end - FLY + 1e-6; u0 += DU) {
      coins.push({ h, u0, jx: rnd() - 0.5, jy: rnd() - 0.5, bow: 0.25 + 0.35 * rnd(), src: 0, d: dot('measured', { px: 3, opacity: 0 }) });
    }
  });
  coins.forEach((c) => { c.src = order[Math.min(NP - 1, Math.floor(usedAt(c.u0 + FLY)))]; });

  /* ── layout ── */
  let L = null;
  function layout() {
    const w = ctx.width || 1100, h = ctx.height || 420, a = w / h;
    const wide = a >= 2;
    const compact = !wide && h < 190; // phones up to 360 px: same message, fewer annotations
    const W = wide ? 21 : 16, H = wide ? 8 : 9;
    camera.userData.animFit = { width: W, height: H, fit: 'contain' };
    const hh = a > W / H ? H / 2 : W / 2 / a, hw = hh * a;
    camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    const ppu = h / (2 * hh), U = (px) => px / ppu, P = (k, lo, hi) => clamp(k * ppu, lo, hi);
    const ts = ctx.textScale;
    const fs = wide
      ? { head: 13, mono: 13, small: 11, pill: 13, name: 12, tick: 11 }
      : { head: 12, mono: compact ? 10 : 11, small: 10, pill: 11, name: 11, tick: 10 };
    const lineH = (size) => U(size * ts * 1.22 + 2);
    const xl = -hw + Math.max(0.3, U(10)), xr = hw - Math.max(0.3, U(10));
    L = { wide, ppu, U };

    const hRs = wide ? fs.head : compact ? 10 : fs.name;
    hL.setSize(fs.head); hR.setSize(hRs);
    rule.forEach((r) => r.l.setSize(fs.mono));
    [cardT, sR, tgt, axT, wallL, binSL, binPL].forEach((l) => l.setSize(fs.small));
    pill.setSize(fs.pill);
    ticks.forEach((k) => k.l.setSize(fs.tick));
    more.l.setSize(fs.small);
    more.l.setText(wide ? `accuracy ${f3(N.logisticAccAtMax)} at ${int(N.maxLabels)}` : `accuracy ${f3(N.logisticAccAtMax)}`);
    const slotTxt = 'no written rule', slotS = compact ? 10 : fs.name;
    rows.forEach((r) => { r.nameL.setSize(fs.name); r.nameL.setText(wide ? r.name : r.short); r.mL.setSize(fs.name); r.track.visible = wide; });
    slot.l.setSize(slotS);

    /* ANM column */
    const pillH = U(fs.pill * ts * 1.22 + 10);
    const HY = hh - (compact ? pillH / 2 + U(5) : Math.max(0.35, lineH(fs.head) / 2 + U(wide ? 10 : 5)));
    const hLt = wide ? 'ANM · written down' : compact ? 'ANM · B/T-priority' : 'ANM · B/T-priority, written down';
    hL.setText(hLt).position.set(xl, HY, 0);
    const pd = wide ? Math.max(0.18, U(9)) : U(compact ? 4 : 6), ls = lineH(fs.mono) + U(wide ? 6 : 1);
    const titleH = wide ? lineH(fs.small) + U(2) : 0;
    cardT.visible = wide;
    const cw = Math.max(wide ? 4.4 : 3.2, U(Math.max(...RULE.map((s) => tw(s, fs.mono, { mono: true })))) + 2 * pd);
    const ch = pd + titleH + 3 * ls + pd * 0.35;
    const pillW = U(tw('0 labels', fs.pill, { weight: 600 }) + 12);
    const cardTop = HY - (compact ? Math.max(lineH(fs.head), pillH) / 2 + U(2) : lineH(fs.head) / 2 + U(wide ? 6 : 3));
    card.setSize(cw, ch, Math.min(0.16, U(8)));
    card.position.set(xl + cw / 2, cardTop - ch / 2, 0);
    cardT.position.set(xl + pd, cardTop - pd - titleH / 2 + U(1), 1);
    rule.forEach((r, i) => r.l.position.set(xl + pd, cardTop - pd - titleH - ls * (i + 0.5), 1));
    let pillX, pillY, leftR, topBot;
    if (wide) {
      tgt.setText(`declines ${pct(N.anmDeclines)} · exactness ${f3(N.anmExactness)}\naccuracy ${f3(N.anmAccuracy)}`);
      pillX = xl; pillY = cardTop - ch - Math.max(0.42, pillH / 2 + U(12));
      tgt.position.set(xl, pillY - pillH / 2 - U(8), 1);
      leftR = xl + Math.max(cw, U(tw(`declines ${pct(N.anmDeclines)} · exactness ${f3(N.anmExactness)}`, fs.small)));
      topBot = cardTop - ch;
    } else if (compact) {
      tgt.setText(`declines ${pct(N.anmDeclines)}\nexactness ${f3(N.anmExactness)}\naccuracy ${f3(N.anmAccuracy)}`);
      pillX = xl + U(tw(hLt, fs.head, { weight: 600 })) + U(6); pillY = HY;
      tgt.position.set(xl + cw + Math.max(0.3, U(9)), cardTop, 1);
      leftR = xl + cw;
      topBot = Math.min(cardTop - ch, cardTop - 3 * U(fs.small * ts * 1.22) - U(4));
    } else {
      tgt.setText(`declines ${pct(N.anmDeclines)}\nexactness ${f3(N.anmExactness)}\naccuracy ${f3(N.anmAccuracy)}`);
      pillX = xl + cw + Math.max(0.3, U(9)); pillY = cardTop - pillH / 2;
      tgt.position.set(pillX, pillY - pillH / 2 - U(3), 1);
      leftR = xl + cw;
      topBot = Math.min(cardTop - ch, pillY - pillH / 2 - U(3) - 3 * U(fs.small * ts * 1.22) - U(3));
    }
    pill.position.set(pillX, pillY, 3);
    const gh = pillH + U(wide ? 12 : compact ? 5 : 8);
    glow.setSize(pillW + U(wide ? 14 : compact ? 6 : 9), gh, gh / 2);
    glow.position.set(pillX + pillW / 2, pillY, 2);

    /* the split, at the bottom */
    const binH = compact ? lineH(fs.small) + U(1) : Math.max(0.36, U(wide ? 16 : 11));
    const binY = compact ? -hh + U(2) + binH / 2 : -hh + Math.max(0.28, lineH(fs.small) + U(wide ? 8 : 4)) + binH / 2 + U(2);
    const gapW = Math.max(0.45, U(16));
    const inner = xr - xl - gapW;
    const sW = (inner * N.scoredCells) / (N.scoredCells + N.poolCells), pW = inner - sW;
    const sX = xl + sW / 2, wX = xl + sW + gapW / 2, pX = xr - pW / 2;
    const wallX = U(wide ? 7 : 4);
    const placeBins = (by) => {
      binS.setSize(sW, binH, Math.min(0.12, binH / 2)); binS.position.set(sX, by, 0);
      binP.setSize(pW, binH, Math.min(0.12, binH / 2)); binP.position.set(pX, by, 0);
      wall.setPoints([[wX, by - binH / 2 - wallX, 1], [wX, by + binH / 2 + wallX, 1]]).setWidth(P(0.07, 2.5, 4));
      wallL.visible = wide;
      wallL.position.set(wX, by + binH / 2 + wallX + U(3), 1);
      binSL.position.set(sX, by - binH / 2 - U(3), 1);
      binPL.setText(`${int(N.poolCells)} label-pool cells` + (wide ? '' : ' · no overlap'));
      binPL.position.set(pX, by - binH / 2 - U(3), 1);
      binSL.visible = binPL.visible = !compact;
      binSLi.visible = binPLi.visible = compact;
      binSLi.position.set(sX, by, 2); binPLi.position.set(pX, by, 2);
      const dpx = P(0.05, 1.5, 3.2);
      const place = (arr, cx, bw, n) => {
        const cols = n / 2, sx = bw / cols;
        arr.forEach((d, i) => {
          d.setPx(dpx);
          d.position.set(cx - bw / 2 + sx * ((i % cols) + 0.5), by + (i < cols ? 1 : -1) * binH * 0.22, 1);
        });
      };
      place(dS, sX, sW, NS);
      place(dP, pX, pW, NP);
    };
    const binTop = binY + binH / 2 + (wide ? wallX + U(3) + lineH(fs.small) : wallX);

    /* heads block */
    const tickW = U(tw(int(N.maxLabels), fs.tick));
    const slotW = wide
      ? Math.max(2.0, U(tw(slotTxt, slotS, { weight: 600 })) + Math.max(0.3, U(12)))
      : U(Math.max(tw('trained heads', hRs, { weight: 600 }), tw(slotTxt, slotS, { weight: 600 })) + 18);
    const nameW = U(Math.max(...rows.map((r) => tw(wide ? r.name : r.short, fs.name))));
    const haloPx = P(0.13, 6, 10);
    const tickOff = Math.max(U(6), U(haloPx + 1));
    let X0, XE, rowTop, rowBot;
    if (wide) {
      X0 = leftR + Math.max(0.5, U(18)) + nameW + Math.max(0.22, U(8));
      XE = xr - slotW - Math.max(0.5, tickW / 2 + U(10));
      const hx = X0 - nameW - Math.max(0.22, U(8));
      hR.setText('trained heads · labels to match ANM').position.set(hx, HY, 0);
      sR.visible = true;
      sR.position.set(hx, HY - lineH(fs.head) / 2 - lineH(fs.small) / 2 - U(1), 0);
      rowTop = sR.position.y - lineH(fs.small) / 2 - Math.max(0.35, U(18));
      rowBot = binTop + U(10) + lineH(fs.small) + lineH(fs.tick) + tickOff + Math.max(0.45, U(18));
    } else {
      X0 = xl + nameW + Math.max(0.2, U(7));
      XE = xr - tickW / 2 - U(2);
      hR.setText('trained heads');
      sR.visible = false;
      rowTop = topBot - lineH(fs.name) / 2 - U(compact ? 4 : 7);
      rowBot = binTop + U(4) + lineH(fs.tick) + tickOff + U(10);
    }
    const rs = clamp((rowTop - rowBot) / 2, U(wide ? 26 : 15), wide ? 1.55 : 1.2);
    rows.forEach((r, i) => { r.y = rowTop - rs * i; });
    placeBins(binY + Math.max(0, rowTop - rowBot - 2 * rs) * 0.5); // spare height: lift the split toward the chart
    const axY = rows[2].y - (compact ? U(8) : Math.max(wide ? 0.45 : 0.3, U(wide ? 18 : 10)));
    const X = (q) => X0 + q * (XE - X0);
    Object.assign(L, { X0, XE, X, axY, rs });

    const aw = P(0.14, 5, 9);
    rows.forEach((r) => {
      r.nameL.position.set(X0 - Math.max(0.15, U(7)), r.y, 1);
      r.track.setPoints([[X0, r.y, 0.5], [XE, r.y, 0.5]]);
      r.bar.setPoints([[X0, r.y, 2], [X(r.q), r.y, 2]]).setWidth(aw);
      r.tip.setPx(aw * 0.75 + 1.5); r.back.setPx(aw * 0.75 + 1.5);
      r.ringPx = aw * 0.8 + 3; r.fillPx = aw * 0.45 + 0.5;
      r.back.position.set(X(r.q), r.y, 3);
      r.mRing.position.set(X(r.q), r.y, 4.2); r.mFill.position.set(X(r.q), r.y, 4.1); r.wave.position.set(X(r.q), r.y, 4);
      r.mL.position.set(X(r.q) + U(r.ringPx + 3), r.y, 5); // right of the ring, on a card-coloured plate
    });
    more.line.setPoints([[X(LOG.q), LOG.y, 1.8], [XE, LOG.y, 1.8]]).setWidth(P(0.045, 1.6, 2.6));
    more.l.position.set(XE + U(aw * 0.75 + 1.5), LOG.y + U(aw * 0.75 + 4), 5);
    const xb = X(QB / 2), bh = U(5), bw2 = U(3);
    axis[0].setPoints([[X0, axY, 1], [xb - U(3), axY, 1]]);
    axis[1].setPoints([[xb + U(3), axY, 1], [XE, axY, 1]]);
    brk[0].setPoints([[xb - U(3) - bw2, axY - bh, 1], [xb - U(3) + bw2, axY + bh, 1]]);
    brk[1].setPoints([[xb + U(3) - bw2, axY - bh, 1], [xb + U(3) + bw2, axY + bh, 1]]);
    ticks.forEach((k, i) => {
      const x = X(qOf(k.v));
      k.t.setPoints([[x, axY - U(4), 1], [x, axY, 1]]);
      k.l.position.set(x, axY - tickOff, 1);
      k.l.visible = wide || i === 0 || i === TICKS.length - 1;
    });
    axT.setText(wide ? 'labelled cells (log scale)' : 'labelled cells (log)');
    if (wide) axT.position.set((X(QB) + XE) / 2, axY - tickOff - lineH(fs.tick), 1);
    else axT.position.set((X(QB) + XE) / 2, axY - tickOff, 1);
    anmDot.setPx(P(0.07, 3, 4.5)); anmHalo.setPx(haloPx);
    anmDot.position.set(X0, axY, 4); anmHalo.position.set(X0, axY, 3.9);
    cursor.setPoints([[0, axY, 1.5], [0, rows[0].y + Math.min(rs * 0.45, U(16)), 1.5]]);
    const cpx = P(0.07, 3, 4.4);
    coins.forEach((c) => c.d.setPx(cpx));

    /* the heads' rule slot: stays empty */
    const lx0 = (x) => x + Math.max(0.2, U(9));
    if (wide) {
      const sx0 = xr - slotW, sy0 = rows[0].y + Math.min(rs * 0.4, U(20)), sy1 = rows[2].y - Math.min(rs * 0.2, U(10)), sh = sy0 - sy1;
      slot.box.setSize(slotW, sh, Math.min(0.16, U(8)));
      slot.box.position.set(sx0 + slotW / 2, (sy0 + sy1) / 2, 0);
      const x1 = sx0 + slotW - Math.max(0.2, U(9));
      slot.lines.forEach((l, i) => {
        const y = sy0 - sh * (0.18 + 0.17 * i);
        l.visible = true;
        l.setPoints([[lx0(sx0), y, 1], [lerp(lx0(sx0), x1, i === 2 ? 0.65 : 1), y, 1]]);
      });
      slot.l.position.set(sx0 + slotW / 2, sy1 + sh * 0.3, 2);
      more.l.visible = true;
    } else {
      /* below the play/pause button, above the logistic head's end label; the heads' header is its title */
      const sx0 = xr - slotW, sy0 = Math.min(cardTop, hh - U(42));
      const labH = lineH(slotS), minH = lineH(hRs) + labH + U(12);
      const accTop = LOG.y + U(aw * 0.75 + 4) + lineH(fs.small);
      const sy1 = compact ? sy0 - minH : Math.min(sy0 - minH, Math.max(accTop + U(4), cardTop - ch)), sh = sy0 - sy1;
      slot.box.setSize(slotW, sh, Math.min(0.16, U(8)));
      slot.box.position.set(sx0 + slotW / 2, (sy0 + sy1) / 2, 0);
      hR.position.set(lx0(sx0), sy0 - U(5) - lineH(hRs) / 2, 1);
      const room = sh - minH;
      slot.lines.forEach((l, i) => {
        l.visible = i === 0 && room > U(6);
        const y = sy0 - U(5) - lineH(hRs) - room / 2;
        l.setPoints([[lx0(sx0), y, 1], [sx0 + slotW - Math.max(0.2, U(9)), y, 1]]);
      });
      slot.l.position.set(sx0 + slotW / 2, sy1 + U(5) + labH / 2, 2);
      /* the logistic head's end note needs its own line above the row; on the smallest phones it is left out
         (the section text states it) */
      more.l.visible = !compact || (rs >= U(19) && accTop + U(3) < sy1);
    }
  }

  const back = { x: 0, y: 0 };
  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    const vis = 1 - seg(u, T.out0, T.out1, ease.inOutSine);
    const { X, U } = L;

    /* 1 ANM writes the question */
    rule.forEach((r, i) => {
      const t0 = T.type0 + T.typeStep * i;
      const n = Math.round(seg(u, t0, t0 + T.typeLen, ease.linear) * r.s.length);
      if (n !== r.n) { r.n = n; r.l.setText(r.s.slice(0, n)); }
      r.l.setOpacity(vis);
    });
    /* 2 its answers, at once */
    const ak = seg(u, T.anm, T.anm + 0.02, ease.outCubic) * vis;
    const flash = seg(u, T.anm, T.anm + 0.02) * (1 - seg(u, T.anm + 0.03, T.anm + 0.12, ease.inOutSine));
    pill.setOpacity(ak); tgt.setOpacity(ak);
    glow.setOpacity(ak * (0.13 + 0.03 * Math.sin(2 * Math.PI * 2 * u)) + 0.22 * flash * vis);
    const gs = 1 + 0.12 * flash;
    glow.scale.set(gs, gs, 1);
    anmDot.setOpacity(ak);
    anmHalo.setOpacity(ak * (0.18 + 0.2 * flash));

    /* 3 labels flow into the trained heads */
    const q = seg(u, T.U0, T.U1, ease.linear);
    const run = seg(u, T.U0 - 0.012, T.U0, ease.linear) * (1 - seg(u, T.U1, T.U1 + 0.03)) * vis;
    cursor.position.x = X(q);
    cursor.setOpacity(0.9 * run);
    const used = usedAt(u);
    dP.forEach((d, k) => {
      const r = rankOf[k];
      d.setOpacity(lerp(0.55, 0.12, clamp(used - r) * vis));
    });
    rows.forEach((r) => {
      const qb = Math.min(q, r.q);
      const on = (u >= T.U0 ? 1 : 0) * vis;
      r.bar.setProgress(r.q > 0 ? qb / r.q : 1).setOpacity(0.85 * on);
      const tipQ = r.more ? q : qb;
      r.tip.position.set(X(tipQ), r.y, 3.2);
      r.tip.setOpacity(on);
      const um = uOfQ(r.q);
      const mk = seg(u, um, um + 0.02, ease.outCubic) * vis;
      const pop = u >= um ? ease.outBack(clamp((u - um) / 0.035)) : 0;
      r.back.setOpacity(mk);
      r.mRing.setOpacity(mk).setPx(r.ringPx * Math.max(0.01, pop));
      r.mFill.setOpacity(mk).setPx(r.fillPx * Math.max(0.01, pop));
      const wk = clamp((u - um) / 0.07);
      r.wave.setPx(r.ringPx * (1 + 1.6 * wk));
      r.wave.setOpacity(u >= um && wk < 1 ? 0.6 * (1 - wk) * vis : 0);
      r.mL.setOpacity(seg(u, um + 0.005, um + 0.03) * vis);
    });
    more.line.setProgress(clamp((q - LOG.q) / (1 - LOG.q))).setOpacity(u >= uOfQ(LOG.q) ? 0.9 * vis : 0);
    more.l.setOpacity(seg(u, T.U1, T.U1 + 0.03) * vis);
    for (const c of coins) {
      const p = (u - c.u0) / FLY;
      if (p < 0 || p > 1 || vis <= 0) { c.d.setOpacity(0); continue; }
      const s = dP[c.src].position;
      const qa = c.h.more ? clamp((c.u0 + FLY - T.U0) / (T.U1 - T.U0)) : Math.min(c.h.q, clamp((c.u0 + FLY - T.U0) / (T.U1 - T.U0)));
      const x0 = s.x + c.jx * U(6), y0 = s.y + c.jy * U(4);
      const x1 = X(qa), y1 = c.h.y;
      const k = ease.inOutSine(p);
      const mx = (x0 + x1) / 2 + c.jx * 0.6, my = Math.max(y0, y1) + c.bow * (L.rs || 1);
      const a1 = 1 - k;
      back.x = a1 * a1 * x0 + 2 * a1 * k * mx + k * k * x1;
      back.y = a1 * a1 * y0 + 2 * a1 * k * my + k * k * y1;
      c.d.position.set(back.x, back.y, 6);
      c.d.setOpacity(0.9 * Math.min(seg(p, 0, 0.12), 1 - seg(p, 0.88, 1)) * vis);
    }

    /* 4 the heads' side never gets a written rule */
    const sk = seg(u, T.slot, T.slot + 0.04) * vis;
    slot.l.setOpacity(sk);
    slot.box.setOpacity(0.55 + 0.45 * sk);
    slot.lines.forEach((l) => l.setOpacity(0.35 + 0.25 * sk));
  }

  let dead = false;
  layout();
  /* label widths drive the layout; lay out again once the page fonts have loaded */
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });
  return {
    scene, camera, period: PERIOD, still: 0.9 * PERIOD,
    update, resize: layout,
    dispose() { dead = true; coins.length = 0; },
  };
}
