/*
 * idx-exp-edit-rewrite — index.html, #edit (Experiment 1 · Change the question), full width between
 * .plain and .bgrid. How the test answers the question, left to right:
 * 1 The written question (ANM teal) rewrites itself in three moves, soft → strict → B/T-priority, and
 *   declares its bar (page guide): soft, every marker counts the same, bar 0.12; strict, CD19, CD3 and
 *   CD16 at weight 2, bar 0.28; B/T-priority, only those three count (×3), then B ×1.5, T ×1.3,
 *   myeloid ×0.5, bar 0.20 as declared.
 * 2 The evidence, each cell's 9 proteins predicted by our head from TEDDY's frozen embedding, stays
 *   unchanged; only the weight column, which belongs to the question, changes. Bar lengths are
 *   illustrative, not data; no value is printed.
 * 3 ANM re-scores the held-out cells (dots; each dot stands for NUMBERS.cells / N cells, so the tray holds
 *   dots in proportion to the page's counts; dot scores and positions are illustrative). Cells whose best
 *   score misses the bar peel off into the “no call” tray: 448 → 2,014 → 3,169, with 0 new labels.
 *   The soft no-calls stay no-calls under strict (page: 2,014 − 1,566 = 448).
 * 4 A dashed orange ghost, the fixed rule re-coded for each question, lies exactly on ANM's bar and
 *   “no call” tray and never separates (same call on 100% of cells). The trained classifier's tray stays
 *   at 40 behind a closed “retrain on new labels” lock: no labels are given, so it never opens.
 * Wide layout 21 × 8 when the stage is at least 2:1 (desktop 21:8), else 16 × 9 (phone 16:9).
 * 12 s seamless loop (soft → strict → B/T-priority → soft); the still frame shows B/T-priority.
 */

/* Every number this scene shows, all on the page (#edit and its guide). Refresh after the
 * official-preprocessing rerun. */
const NUMBERS = {
  cells: 16750,                     // held-out cells
  proteins: 9,                      // predicted proteins ANM reads
  noCall: [272, 1493, 2551],        // ANM “no calls”: soft, strict, B/T-priority
  classifierNoCall: [43, 43, 43],   // trained classifier “no calls”, same three questions
  sameCalls: '100%',                // cells where ANM and the re-coded fixed rule make the same call
  newLabels: 0,
  // the written questions (page guide)
  bar: [0.12, 0.28, 0.20],          // B/T-priority as declared (0.1333 on the shared score scale)
  keyWeight: [1, 2, 3],             // CD19, CD3, CD16
  otherWeight: [1, 1, 0],           // the other six markers
  lineageWeight: { B: 1.5, T: 1.3, myeloid: 0.5 },   // B/T-priority only
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 8 });
  const PERIOD = 12;
  const N = 160;                    // cell dots in the field
  const Q = [0, 1, 2];
  const BAR = NUMBERS.bar;
  const LW = NUMBERS.lineageWeight;
  const fmt = (n) => Number(n).toLocaleString('en-US');
  const f2 = (v) => Number(v).toFixed(2);
  const W = (q, key) => (key ? NUMBERS.keyWeight[q] : NUMBERS.otherWeight[q]);

  const add = (o) => { scene.add(o); return o; };
  const labels = [];
  const lab = (text, o = {}) => {
    const l = ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o));
    l.userData.base = o.size || 12;
    l.userData.nbase = o.nsize || o.size || 12;   // size on the narrow (phone) layout
    l.userData.text = String(text);
    labels.push(l);
    return add(l);
  };
  const setT = (l, text) => { l.userData.text = text; l.setText(text); };
  const ln = (o = {}) => add(ctx.line([[0, 0], [1, 0]], o));
  const dot = (color, o = {}) => add(ctx.dot([0, 0, 0], Object.assign({ color }, o)));
  const rect = (o = {}) => add(ctx.box(1, 1, Object.assign({ radius: 0.16 }, o)));

  /* text width in world units, with the core's label padding */
  const mc = document.createElement('canvas').getContext('2d');
  const P = { ppu: 50, f: 1 };
  function tw(text, size, weight = 500, bg = false, mono = false) {
    const cs = getComputedStyle(ctx.figure);
    const fam = mono ? ((cs.getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace') : (cs.fontFamily || 'system-ui, sans-serif');
    mc.font = `${weight} ${size * (L.wide ? P.f : 1) * ctx.textScale}px ${fam}`;
    let w = 0;
    for (const l of String(text).split('\n')) w = Math.max(w, mc.measureText(l).width);
    return (w + (bg ? 12.8 : 3)) / P.ppu;
  }

  /* ── layouts ── */
  const WIDE = {
    W: 21, H: 8, wide: true, hy: 3.05, cy: 0.27,
    card: [-10.35, -5.75, -1.65, 2.65],
    arrow1: [[-5.68, 0.5], [-5.38, 0.5]],
    table: [-5.3, 0.95, -2.65, 2.65],
    field: [2.05, 5.6, -2.65, 2.65],
    arrow2: [[1.22, 0], [1.93, 0]],
    tray: [6.1, 10.3, -2.65, -0.45], trayHy: -0.07,     // below the classifier: no-calls drop down and right
    cls: [6.1, 10.3, 0.45, 2.65], clsHy: 3.05,
    dotsSide: 'left', ds: 0.26,
  };
  const NARROW = {
    W: 16, H: 9, wide: false, hy: 3.95, cy: -0.1,
    card: [-7.8, -2.75, 0.4, 3.55],
    table: [-2.2, 7.8, 0.4, 3.55],
    field: [-2.2, 2.9, -4.3, -0.9],
    tray: [-7.8, -2.75, -4.3, -0.9],
    cls: [3.45, 7.8, -4.3, -0.9], clsHy: -0.5, botHy: -0.5, trayHy: -0.5,
    dotsSide: 'right', ds: 0.38,
  };
  let L = WIDE;

  /* ── 1 the written question ── */
  const TABS = ['soft', 'strict', 'B/T-priority'];
  const card = rect({ color: 'soft', stroke: 'line' });
  const cardCap = lab('written question', { anchor: 'left', color: 'muted', size: 11.5 });
  const tabs = TABS.map((t) => ({
    on: lab(t, { bg: 'accent', color: 'card', weight: 600, size: 11.5, anchor: 'left', order: 12 }),
    off: lab(t, { bg: 'grid', color: 'muted', size: 11.5, anchor: 'left', order: 11 }),
  }));
  const titles = TABS.map((t) => lab(t, { color: 'accent', weight: 600, size: 13, anchor: 'left' }));
  const steps = TABS.map(() => ({ ring: dot('accent', { hollow: true, ring: 0.34 }), fill: dot('accent') }));
  const RULES = {
    wide: [
      'every marker\ncounts the same',
      `CD19, CD3, CD16\nat weight ${NUMBERS.keyWeight[1]}`,
      `only CD19, CD3, CD16 (×${NUMBERS.keyWeight[2]}),\nthen B ×${LW.B} · T ×${LW.T} · M ×${LW.myeloid}`,
    ],
    narrow: [
      'every marker\ncounts the same',
      `CD19, CD3, CD16\nat weight ${NUMBERS.keyWeight[1]}`,
      'only CD19, CD3,\nCD16 count',
    ],
  };
  const rules = Q.map(() => lab('', { anchor: 'top-left', size: 12 }));
  const pen = ln({ color: 'accent', width: 2 });
  const penTip = dot('accent');
  const cBars = Q.map((q) => lab(`bar ${f2(BAR[q])}`, { anchor: 'left', color: 'accent', weight: 600, size: 16, mono: true }));
  const hint = lab('misses the bar → “no call”', { anchor: 'left', color: 'muted', size: 11 });
  const divider = ln({ color: 'line', width: 1 });
  const note = lab('', { anchor: 'left', color: 'muted', size: 11.5 });
  const arrow1 = add(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.6, head: 8 }));

  /* ── 2 the evidence table: 9 predicted proteins (fixed) and the question's weight column ── */
  const MARK = [['CD19', 'CD72', 'CD22'], ['CD3', 'CD2', 'CD5'], ['CD16', 'CD11c', 'CD36']];
  const LIN = [['B', 'b', LW.B], ['T', 't', LW.T], ['myeloid', 'm', LW.myeloid]];
  const EVID = [0.74, 0.52, 0.6, 0.32, 0.26, 0.38, 0.16, 0.21, 0.12]; // illustrative lengths, never printed
  const stack = [0, 1].map(() => rect({ color: null, stroke: 'line', strokeWidth: 1.2 }));
  const frame = rect({ color: 'soft', stroke: 'line' });
  const wBand = rect({ color: 'accent', opacity: 0.07, radius: 0.1 });
  const tHead = lab(`${NUMBERS.proteins} predicted proteins · unchanged`, { anchor: 'left', weight: 600 });
  const wHead = lab('weight', { color: 'accent', weight: 600, size: 11.5 });
  const rows = [];
  MARK.forEach((ms, g) => ms.forEach((m, j) => {
    const key = j === 0;
    rows.push({
      g, j, key, v: EVID[g * 3 + j],
      bar: ln({ color: LIN[g][1], width: 8 }),
      name: lab(m, { size: 11, anchor: 'right' }),
      nameQ: key ? Q.map((q) => lab(`${m} ×${W(q, true)}`, { size: 11, weight: 600 })) : null,
      wRing: dot('faint', { hollow: true, ring: 0.3 }),
      wDot: dot('accent'),
      wLab: Q.map((q) => lab(`×${W(q, key)}`, { size: 11, color: 'accent', weight: 600, anchor: 'left' })),
    });
  }));
  const lins = LIN.map(([n, c, w]) => ({
    base: ln({ color: c, width: 1.6 }),
    name: lab(n, { color: c, weight: 600, anchor: 'left' }),
    w: lab(`×${w}`, { color: 'accent', weight: 600, size: 11, anchor: 'left' }),
  }));
  const arrow2 = add(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.6, head: 8 }));
  const anmLab = lab('ANM', { color: 'accent', weight: 700, size: 11.5, anchor: 'bottom' });
  const arrowPulse = dot('accent');

  /* ── 3 the held-out cells, the bar plane, the no-call tray; 4 the ghost and the classifier ── */
  const field = rect({ color: 'soft', stroke: 'line' });
  const fHead = lab(`${fmt(NUMBERS.cells)} held-out cells`, { anchor: 'left', weight: 600 });
  const zone = add(ctx.box(1, 1, { color: 'accent', opacity: 0.08, radius: 0 }));
  const barLine = ln({ color: 'accent', width: 2 });
  const ghostBar = ln({ color: 'teddy', width: 2, dashed: [5, 5] });
  const fBars = Q.map((q) => lab(`bar ${f2(BAR[q])}`, { anchor: 'bottom-right', color: 'accent', weight: 600, size: 11, bg: 'soft', bgOpacity: 0.85 }));
  const zoneLab = lab('no call', { anchor: 'bottom-left', color: 'accent', size: 10.5, opacity: 0.8 });
  const recoded = lab('re-coded', { anchor: 'bottom-left', color: 'teddy', weight: 600, size: 10.5, bg: 'soft', bgOpacity: 0.85 });

  const tray = rect({ color: 'soft', stroke: 'accent', strokeWidth: 1.6 });
  const ghostTray = rect({ color: null, stroke: 'teddy', strokeWidth: 1.6, dashed: [5, 5] });
  const aHead1 = lab('ANM', { color: 'accent', weight: 700, anchor: 'left' });
  const aHead2 = lab('', { color: 'teddy', weight: 600, anchor: 'left' });
  const aCap = lab('no call', { color: 'muted', size: 11, anchor: 'left' });
  const aSame = lab('', { color: 'teddy', size: 11, anchor: 'right', weight: 600 });
  const aCounts = Q.map((q) => lab(fmt(NUMBERS.noCall[q]), { size: 22, nsize: 15, weight: 600, color: 'accent', anchor: 'left' }));
  const mkHist = (vals, color) => ({
    e: vals.map((v) => ({ off: lab(fmt(v), { size: 11, color: 'muted', anchor: 'left' }), on: lab(fmt(v), { size: 11, color, weight: 600, anchor: 'left' }) })),
    ar: [0, 1].map(() => lab('→', { size: 11, color: 'faint', anchor: 'left' })),
  });
  const aHist = mkHist(NUMBERS.noCall, 'accent');

  const cls = rect({ color: 'soft', stroke: 'train', strokeWidth: 1.6 });
  const cHead = lab('trained classifier', { color: 'train', weight: 600, anchor: 'left' });
  const cCap = lab('no call', { color: 'muted', size: 11, anchor: 'left' });
  const cCounts = Q.map((q) => lab(fmt(NUMBERS.classifierNoCall[q]), { size: 22, nsize: 15, weight: 600, color: 'train', anchor: 'left' }));
  const cHist = mkHist(NUMBERS.classifierNoCall, 'train');
  const lockBody = add(ctx.box(1, 1, { color: 'train', radius: 0.04, opacity: 0.9 }));
  const shackle = ln({ color: 'train', width: 2.2 });
  const keyhole = dot('card');
  const lockLab = lab('retrain on\nnew labels', { size: 10.5, color: 'muted' });

  /* ── cells: illustrative scores per question; tray membership = score under the bar ── */
  const rnd = ctx.rand(11);
  const rr = (a, b) => a + (b - a) * rnd();
  const nq = NUMBERS.noCall.map((n) => Math.max(1, Math.round((n / NUMBERS.cells) * N)));
  const nS = Math.min(nq[0], nq[1], nq[2]);
  const nStrict = Math.max(0, nq[1] - nS);
  let nStay = nStrict - Math.min(nStrict, Math.round(nStrict * 0.27));
  if (nS + nStay > nq[2]) nStay = Math.max(0, nq[2] - nS);
  const nRet = nStrict - nStay;
  const nNew = Math.max(0, nq[2] - nS - nStay);
  const lo = (q) => BAR[q] - 0.035, hi = (q) => BAR[q] + 0.04;
  const cells = [];
  const mk = (grp, sc) => cells.push({
    grp, sc, xr: rnd(), k: rnd(), col: 'muted', fx: 0,
    fill: dot('muted'), ring: dot('teddy', { hollow: true, ring: 0.3 }),
  });
  for (let i = 0; i < nS; i++) { const v = rr(0.02, Math.min(...BAR) - 0.035); mk(0, [v, v, v]); }
  for (let i = 0; i < nStay + nRet; i++) {
    const s0 = rr(hi(0), Math.max(hi(0), lo(1))), s1 = clamp(s0 - rr(0, 0.04), 0.02, lo(1));
    const s2 = i < nStay ? rr(0.02, lo(2)) : rr(Math.max(hi(2), 0.3), 0.55);
    mk(i < nStay ? 1 : 2, [s0, s1, s2]);
  }
  for (let i = 0; i < nNew; i++) {
    const s0 = rr(Math.max(hi(1) + 0.02, 0.34), 0.5), s1 = Math.max(hi(1), s0 + rr(-0.02, 0.03));
    mk(3, [s0, s1, rr(0.03, lo(2))]);
  }
  while (cells.length < N) {
    const b = rr(0.33, 0.96);
    mk(4, [b, clamp(b + rr(-0.04, 0.04), hi(1) + 0.01, 0.97), clamp(b + rr(-0.07, 0.07), hi(2) + 0.01, 0.97)]);
  }
  cells.forEach((c) => {
    c.tray = Q.map((q) => c.sc[q] < BAR[q]);
    c.slot = [0, 0, 0];
  });
  Q.forEach((q) => { let n = 0; cells.forEach((c) => { if (c.tray[q]) c.slot[q] = n++; }); });
  const NMAX = Math.max(...Q.map((q) => cells.filter((c) => c.tray[q]).length));

  /* ── layout ── */
  const T = { cols: 8, ds: 0.26, dx0: 0, dy0: 0 };
  const F = { x0: 0, x1: 1, y0: 0, y1: 1, pad: 0.2 };
  const yOf = (v) => F.y0 + F.pad + clamp(v) * (F.y1 - F.y0 - 2 * F.pad);
  const slotPos = (j, out) => out.set(T.dx0 + (j % T.cols) * T.ds, T.dy0 - Math.floor(j / T.cols) * T.ds, 0);
  const place = (o, x, y, z = 0) => { o.position.set(x, y, z); return o; };
  const sizeBox = (b, [x0, x1, y0, y1], z = 0, r) => { b.setSize(x1 - x0, y1 - y0, r); place(b, (x0 + x1) / 2, (y0 + y1) / 2, z); };
  const px = (k, a, b) => clamp(k * P.ppu, a, b);
  const LOCK = { x: 0, y: 0, s: 0.3 };

  function layoutTray(b, [x0, x1, y0, y1], cap, counts, hist, countRight = false) {
    place(cap, x0 + 0.18, y1 - 0.3);
    const cy = (y1 - 0.55 + y0 + 0.55) / 2;
    counts.forEach((l) => { l.setAnchor(countRight ? 'right' : 'left'); place(l, countRight ? x1 - 0.2 : x0 + 0.18, cy); });
    let x = x0 + 0.18;
    hist.e.forEach((e, q) => {
      place(e.off, x, y0 + 0.28); place(e.on, x, y0 + 0.28);
      x += tw(e.on.userData.text, 11, 600) + 0.06;
      if (q < 2) { place(hist.ar[q], x, y0 + 0.28); x += tw('→', 11) + 0.06; }
    });
    sizeBox(b, [x0, x1, y0, y1], 0.2);
    return cy;
  }

  function layout() {
    const w = ctx.width || 1100, h = ctx.height || 420, a = w / h;
    L = a >= 2 ? WIDE : NARROW;
    const wide = L.wide;
    camera.userData.animFit = { width: L.W, height: L.H, fit: 'contain' };
    const hh = a > L.W / L.H ? L.H / 2 : L.W / 2 / a;
    camera.left = -hh * a; camera.right = hh * a; camera.top = hh + L.cy; camera.bottom = -hh + L.cy;
    camera.updateProjectionMatrix();
    P.ppu = h / (2 * hh);
    P.f = wide ? clamp(P.ppu / 50, 0.72, 1) : 1;
    labels.forEach((l) => l.setSize(wide ? l.userData.base * P.f : l.userData.nbase));

    /* card */
    const [cx0, cx1, cy0, cy1] = L.card;
    sizeBox(card, L.card, 0);
    place(cardCap, cx0, L.hy);
    Q.forEach((q) => setT(rules[q], (wide ? RULES.wide : RULES.narrow)[q]));
    tabs.forEach((t) => { t.on.visible = t.off.visible = wide; });
    titles.forEach((l) => { l.visible = !wide; });
    steps.forEach((s) => { s.ring.visible = s.fill.visible = !wide; });
    hint.visible = divider.visible = arrow1.visible = wide;
    stack.forEach((s) => { s.visible = wide; });
    arrow2.visible = anmLab.visible = arrowPulse.visible = wide;
    wHead.visible = wide;
    if (wide) {
      let x = cx0 + 0.22;
      tabs.forEach((t) => { place(t.on, x, cy1 - 0.42); place(t.off, x, cy1 - 0.42); x += tw(t.on.userData.text, 11.5, 600, true) + 0.1; });
      L.ruleY = cy1 - 0.85; L.penY = L.ruleY - 0.82; L.barY = L.penY - 0.55;
      place(hint, cx0 + 0.25, L.barY - 0.6);
      divider.setPoints([[cx0 + 0.25, L.barY - 1.15, 1], [cx1 - 0.25, L.barY - 1.15, 1]]);
      setT(note, `${NUMBERS.newLabels} new labels · no retrain`);
      place(note, cx0 + 0.25, L.barY - 1.6);
      arrow1.set(L.arrow1[0], L.arrow1[1]);
    } else {
      titles.forEach((l) => place(l, cx0 + 0.2, cy1 - 0.35));
      steps.forEach((s, k) => {
        const x = cx1 - 0.2 - (2 - k) * 0.42;
        place(s.ring, x, cy1 - 0.35, 1).setPx(4.2); place(s.fill, x, cy1 - 0.35, 1.1).setPx(2.4);
      });
      L.ruleY = cy1 - 0.62; L.penY = L.ruleY - 1.23; L.barY = L.penY - 0.45;
      setT(note, `${NUMBERS.newLabels} new labels`);
      place(note, cx0 + 0.2, cy0 + 0.32);
    }
    rules.forEach((l) => place(l, cx0 + (wide ? 0.25 : 0.2), L.ruleY));
    pen.setPoints([[cx0 + 0.25, L.penY, 1], [cx1 - 0.25, L.penY, 1]]);
    penTip.setPx(3);
    cBars.forEach((l) => place(l, cx0 + (wide ? 0.25 : 0.2), L.barY));

    /* evidence table */
    const [tx0, tx1, ty0, ty1] = L.table;
    sizeBox(frame, L.table, 0);
    setT(tHead, `${NUMBERS.proteins} predicted proteins · ${wide ? 'unchanged' : 'fixed'}`);
    place(tHead, tx0, L.hy);
    if (wide) {
      stack.forEach((s, i) => sizeBox(s, [tx0 + 0.09 * (i + 1), tx1 + 0.09 * (i + 1), ty0 + 0.09 * (i + 1), ty1 + 0.09 * (i + 1)], -0.2 - 0.1 * i));
      const RY = [2.3, 1.8, 1.3, 0.5, 0, -0.5, -1.3, -1.8, -2.3];
      const nameX = tx0 + 2.45, barX0 = tx0 + 2.62, barLen = 2.2, wX = tx1 - 0.85;
      sizeBox(wBand, [wX - 0.4, tx1 - 0.1, ty0 + 0.1, ty1 - 0.1], 0.1, 0.1);
      place(wHead, wX + 0.18, L.hy);
      rows.forEach((r, i) => {
        const y = (r.y = RY[i]);
        r.bar.setPoints([[barX0, y, 1], [barX0 + barLen * r.v, y, 1]]).setWidth(px(0.19, 5, 10));
        place(r.name, nameX, y); r.name.visible = true;
        if (r.nameQ) r.nameQ.forEach((l) => { l.visible = false; });
        place(r.wDot, wX, y, 2); place(r.wRing, wX, y, 2);
        r.wLab.forEach((l) => { place(l, wX + 0.25, y); l.visible = true; });
      });
      lins.forEach((ln0, g) => {
        const ys = [RY[g * 3], RY[g * 3 + 2]];
        ln0.base.setPoints([[barX0, ys[0] + 0.2, 1], [barX0, ys[1] - 0.2, 1]]);
        place(ln0.name, tx0 + 0.2, RY[g * 3 + 1]);
        place(ln0.w, tx0 + 0.2, RY[g * 3 + 2]);
      });
    } else {
      const CX = [-1.35, -0.4, 0.55, 1.95, 2.9, 3.85, 5.25, 6.2, 7.15];
      const yb = 1.8, hMax = 1.45, wY = 1.08, linY = 0.72;
      sizeBox(wBand, [tx0 + 0.2, tx1 - 0.2, wY - 0.2, wY + 0.2], 0.1, 0.1);
      rows.forEach((r, i) => {
        const x = (r.x = CX[i]);
        r.bar.setPoints([[x, yb, 1], [x, yb + hMax * r.v, 1]]).setWidth(px(0.34, 5, 9));
        r.name.visible = false;
        if (r.nameQ) r.nameQ.forEach((l) => { place(l, x, yb - 0.28); l.visible = true; });
        place(r.wDot, x, wY, 2); place(r.wRing, x, wY, 2);
        r.wLab.forEach((l) => { l.visible = false; });
      });
      lins.forEach((ln0, g) => {
        const gx = (CX[g * 3] + CX[g * 3 + 2]) / 2;
        ln0.base.setPoints([[CX[g * 3] - 0.35, yb, 1], [CX[g * 3 + 2] + 0.35, yb, 1]]);
        const nw = tw(ln0.name.userData.text, 12, 600), ww = tw(ln0.w.userData.text, 11, 600);
        place(ln0.name, gx - (nw + ww) / 2, linY);          // name + weight centred on the group
        place(ln0.w, gx - (nw + ww) / 2 + nw + 0.04, linY);
      });
    }
    const wb = px(0.062, 2.3, 3.4);
    rows.forEach((r) => { r.wRing.setPx(wb); r.wb = wb; });
    if (wide) { arrow2.set(L.arrow2[0], L.arrow2[1]); place(anmLab, (L.arrow2[0][0] + L.arrow2[1][0]) / 2, L.arrow2[0][1] + 0.14); arrowPulse.setPx(3.4); }

    /* field */
    [F.x0, F.x1, F.y0, F.y1] = L.field;
    F.pad = wide ? 0.22 : 0.18;
    sizeBox(field, L.field, 0);
    place(fHead, F.x0, L.botHy != null ? L.botHy : L.hy);
    zone.setSize(F.x1 - F.x0 - 0.1, 1, 0);
    zone.position.x = (F.x0 + F.x1) / 2;
    const dp = px(0.055, 1.8, 3.2);
    cells.forEach((c) => {
      c.fx = F.x0 + F.pad + c.xr * (F.x1 - F.x0 - 2 * F.pad);
      c.fill.setPx(dp); c.ring.setPx(dp + px(0.04, 1.3, 2.4));
    });
    place(zoneLab, F.x0 + 0.14, F.y0 + 0.12, 0); zoneLab.visible = wide;

    /* ANM tray + ghost */
    const [ax0, ax1, ay0, ay1] = L.tray;
    const ahy = L.trayHy;
    setT(aHead2, wide ? '= fixed rule, re-coded' : '= re-coded rule');
    place(aHead1, ax0, ahy);
    place(aHead2, ax0 + tw('ANM', 12, 700) + 0.02, ahy);
    setT(aSame, wide ? `${NUMBERS.sameCalls} same calls` : 'same calls');
    place(aSame, ax1 - 0.18, ay1 - 0.3);
    const left = L.dotsSide === 'left';
    layoutTray(tray, L.tray, aCap, aCounts, aHist, left);
    sizeBox(ghostTray, L.tray, 0.3);
    // the dots sit on the tray's side nearest the field, the count on the other side
    const dx1 = left ? ax0 + 2.35 : ax1 - 0.2;
    T.ds = L.ds; T.dx0 = left ? ax0 + 0.25 : ax0 + 2.2; T.dy0 = ay1 - (wide ? 0.62 : 0.82);
    T.cols = Math.max(1, Math.floor((dx1 - T.dx0) / T.ds) + 1);
    while (Math.ceil(NMAX / T.cols) > 1 + Math.floor((T.dy0 - (ay0 + 0.5)) / T.ds) && T.ds > 0.12) {
      T.ds *= 0.92; T.cols = Math.floor((dx1 - T.dx0) / T.ds) + 1;
    }

    /* classifier tray + lock */
    const [kx0, kx1] = L.cls;
    place(cHead, kx0, L.clsHy);
    const kcy = layoutTray(cls, L.cls, cCap, cCounts, cHist);
    LOCK.s = px(0.34, 10, 18) / P.ppu;
    LOCK.x = kx0 + (kx1 - kx0) * 0.72; LOCK.y = kcy + 0.2;
    const s = LOCK.s;
    lockBody.setSize(s, s * 0.78, s * 0.14);
    place(lockBody, LOCK.x, LOCK.y, 1);
    const arc = [[LOCK.x - 0.3 * s, LOCK.y + 0.3 * s, 1.1]];
    for (let i = 0; i <= 12; i++) { const th = Math.PI - (Math.PI * i) / 12; arc.push([LOCK.x + 0.3 * s * Math.cos(th), LOCK.y + 0.62 * s + 0.3 * s * Math.sin(th), 1.1]); }
    arc.push([LOCK.x + 0.3 * s, LOCK.y + 0.3 * s, 1.1]);
    shackle.setPoints(arc).setWidth(px(0.04, 1.6, 2.4));
    place(keyhole, LOCK.x, LOCK.y - 0.02 * s, 1.2).setPx(px(0.03, 1.3, 2));
    place(lockLab, LOCK.x, LOCK.y - 0.39 * s - 0.1, 0).setAnchor('top');
  }

  /* ── timeline: soft (hold) → strict → B/T-priority → soft; each rewrite takes 0.16 of the loop ── */
  const TR = [{ a: 0, b: 1, u0: 0.2, u1: 0.36 }, { a: 1, b: 2, u0: 0.53, u1: 0.69 }, { a: 2, b: 0, u0: 0.84, u1: 1.0 }];
  function stateAt(u) {
    for (const tr of TR) if (u >= tr.u0 && u < tr.u1) return { a: tr.a, b: tr.b, s: (u - tr.u0) / (tr.u1 - tr.u0) };
    const q = u < 0.2 ? 0 : u < 0.53 ? 1 : 2;
    return { a: q, b: q, s: 1 };
  }
  const wts = (a, b, k) => { const w = [0, 0, 0]; if (a === b) w[a] = 1; else { w[a] = 1 - k; w[b] = k; } return w; };
  function fade(arr, a, b, k) {
    const w = arr[a].userData.text === arr[b].userData.text ? wts(a, a, 0) : wts(a, b, k);
    arr.forEach((l, q) => l.setOpacity(w[q]));
  }
  function fadeHist(hist, a, b, k) {
    const w = wts(a, b, k);
    hist.e.forEach((e, q) => { e.on.setOpacity(w[q]); e.off.setOpacity(1 - w[q]); });
  }

  const _p = new THREE.Vector3(), _q = new THREE.Vector3();
  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    const { a, b, s } = stateAt(u);
    const moving = a !== b;

    /* 1 the question rewrites itself */
    const out = 1 - seg(s, 0.02, 0.18), inn = seg(s, 0.16, 0.34);
    const crossText = (arr, dy) => arr.forEach((l, q) => {
      const o = !moving ? (q === a ? 1 : 0) : q === a ? out : q === b ? inn : 0;
      l.setOpacity(o);
      if (dy != null) l.position.y = dy + (moving && q === a ? 0.12 * (1 - out) : moving && q === b ? -0.12 * (1 - inn) : 0);
    });
    crossText(rules, L.ruleY);
    crossText(cBars, L.barY);
    crossText(titles, null);
    const pk = seg(s, 0.06, 0.34, ease.inOutSine), po = moving ? pulse(s, 0.04, 0.5, 0.08) : 0;
    pen.setProgress(pk).setOpacity(po);
    penTip.position.copy(ctx.along(pen, pk)).setZ(1.2);
    penTip.setOpacity(po);
    const tw8 = wts(a, b, seg(s, 0.12, 0.3));
    tabs.forEach((tb, q) => { tb.on.setOpacity(tw8[q]); tb.off.setOpacity(1 - tw8[q]); });
    steps.forEach((st, q) => st.fill.setOpacity(tw8[q]));

    /* 2 the weight column follows the question; the evidence stays */
    const kw = seg(s, 0.2, 0.44);
    rows.forEach((r) => {
      const w = lerp(W(a, r.key), W(b, r.key), kw), on = clamp(w);
      r.bar.setOpacity(0.3 + 0.6 * on);
      r.name.setOpacity(0.45 + 0.55 * on);
      r.wDot.setPx(r.wb * Math.sqrt(Math.max(w, 0.02))).setOpacity(clamp(w * 1.5));
      r.wRing.setOpacity(0.8 * (1 - on));
      if (L.wide) fade(r.wLab, a, b, kw);
      else if (r.nameQ) fade(r.nameQ, a, b, kw);
    });
    const lw = wts(a, b, kw)[2];
    lins.forEach((l) => l.w.setOpacity(lw));
    if (L.wide) {
      const ak = seg(s, 0.34, 0.58, ease.linear), ao = moving ? pulse(s, 0.34, 0.6, 0.04) : 0;
      _p.set(lerp(L.arrow2[0][0], L.arrow2[1][0] - 0.12, ak), L.arrow2[0][1], 2);
      arrowPulse.position.copy(_p);
      arrowPulse.setOpacity(ao);
    }

    /* 3 re-score: the bar plane moves, cells that miss it peel off into the no-call tray */
    const r = seg(s, 0.38, 0.6);
    const yBar = yOf(lerp(BAR[a], BAR[b], r));
    barLine.setPoints([[F.x0 + 0.05, yBar, 1], [F.x1 - 0.05, yBar, 1]]);
    ghostBar.setPoints([[F.x0 + 0.05, yBar, 1.2], [F.x1 - 0.05, yBar, 1.2]]).setDashPhase(u, 18);
    const zh = Math.max(0.001, yBar - (F.y0 + 0.05));
    zone.scale.y = zh; zone.position.y = F.y0 + 0.05 + zh / 2; zone.position.z = 0.5;
    const bw = wts(a, b, r);
    fBars.forEach((l, q) => { place(l, F.x1 - 0.1, yBar + 0.07); l.setOpacity(bw[q]); });
    const gp = moving ? pulse(s, 0.34, 0.86, 0.12) : 0;
    place(recoded, F.x0 + 0.1, yBar + 0.07);
    recoded.setOpacity(gp);
    ghostBar.setWidth(2 + 1.2 * gp);
    ghostTray.outline.setWidth(1.6 + 1.2 * gp);
    ghostTray.outline.setDashPhase(u, 10);

    for (const c of cells) {
      const A = c.tray[a], B = c.tray[b];
      const ya = yOf(c.sc[a]), yb = yOf(c.sc[b]);
      let x, y, tau, hot = 0;
      if (!A && !B) { x = c.fx; y = lerp(ya, yb, r); tau = 0; }
      else if (A && B) { slotPos(c.slot[a], _p); slotPos(c.slot[b], _q); x = lerp(_p.x, _q.x, r); y = lerp(_p.y, _q.y, r); tau = 1; }
      else if (B) {
        const p = seg(s, 0.62 + 0.14 * c.k, 0.78 + 0.14 * c.k);
        slotPos(c.slot[b], _q);
        x = lerp(c.fx, _q.x, p); y = lerp(lerp(ya, yb, r), _q.y, p); tau = p; hot = seg(r, 0.75, 1);
      } else {
        const p = seg(s, 0.56 + 0.1 * c.k, 0.72 + 0.1 * c.k);
        slotPos(c.slot[a], _p);
        x = lerp(_p.x, c.fx, p); y = lerp(_p.y, yb, p); tau = 1 - p;
      }
      const nc = Math.max(tau, hot);
      const col = nc > 0.02 ? 'accent' : 'muted';
      if (col !== c.col) { c.fill.setColor(col); c.col = col; }
      c.fill.position.set(x, y, 3 + tau);
      c.fill.setOpacity(lerp(0.55, 1, nc));
      c.ring.position.set(x, y, 3.1 + tau);
      c.ring.setOpacity(0.9 * tau);
    }
    const kc = seg(s, 0.74, 0.9);
    fade(aCounts, a, b, kc);
    fade(cCounts, a, b, kc);
    fadeHist(aHist, a, b, kc);
    fadeHist(cHist, a, b, kc);
  }

  layout();
  let alive = true;
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (alive) { layout(); ctx.requestRender(); } });
  return {
    scene, camera, period: PERIOD, still: 0.78 * PERIOD,
    update, resize: layout,
    dispose() { alive = false; cells.length = 0; rows.length = 0; },
  };
}
