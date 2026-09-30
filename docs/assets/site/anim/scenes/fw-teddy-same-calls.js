/*
 * fw-teddy-same-calls — anm-framework.html #teddy-same-calls ("Finding: ANM's calls equal the fixed rule").
 * Three beats in a 12 s seamless loop.
 * 1 Scores. Left: the fixed rule's lineage scores (B, T, myeloid); right: ANM's field scores. All nine
 *   markers enter at once, so both charts rise together (three stacked marker slices per bar) and ANM's
 *   starts as a copy. A teal lens sweeps the ANM side: under it every bar and the threshold stretch by the
 *   same factor G·c (S_ANM = G c S_rule, θ_ANM = G c θ_rule); a dashed tick keeps the rule height inside
 *   each ANM bar. Checks light: same ranking, same bar, same call (call T on both sides). The bars then
 *   morph to a second cell whose best lineage stays under the bar: no call on both sides. Bar heights,
 *   slices and the stretch factor are illustrative (the page gives no value for G·c).
 * 2 Every cell. The two charts shrink into two deciders. A stream of cells splits into both; each cell's
 *   two copies leave the deciders in the same colour and drop into matching bins (B, T, myeloid, no call),
 *   tied by a thin line across the bin names. The question tabs (soft, strict, B/T-priority) step on; per
 *   tab a counter runs to the page's 16,750 held-out cells with 0 differences. Which bin a dot takes is
 *   illustrative, not data.
 * 3 What it adds. A grey "better accuracy" plate is struck through and flips into a lit audit record:
 *   written question, per-marker reason, flip distance (the overview page's audit record), with
 *   "not better accuracy" below, as the section says: an auditable record, not better accuracy.
 * Still frame (reduced motion): beat 1 complete, both sides call T, all three checks lit.
 */

/* Page numbers shown by this scene (anm-framework.html #teddy-same-calls). Refresh after the rerun. */
const NUMBERS = {
  heldOutCells: 16750, // "On all 16,750 held-out cells and three questions"
  differences: 0, //      "calls are identical"
};
const QUESTIONS = ['soft', 'strict', 'B/T-priority']; // the three written questions (section 5.2)

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const text = (s, x, y, o = {}) => {
    const l = ctx.label(s, Object.assign({ size: 12, color: 'ink' }, o));
    l.position.set(x, y, 10);
    return add(l);
  };
  const fmt = (n) => Math.round(n).toLocaleString('en-US');
  const colorOf = new Map();
  const tint = (obj, c) => { if (colorOf.get(obj) !== c) { colorOf.set(obj, c); obj.setColor(c); } };
  const sized = new Map();
  const fitBox = (b, w, h, r) => {
    const k = `${w.toFixed(3)} ${h.toFixed(3)}`;
    if (sized.get(b) !== k) { sized.set(b, k); b.setSize(w, h, r); }
  };

  /* ── beat 1: the rule's bars and ANM's bars ─────────────────────────────────────────── */
  const LIN = [{ n: 'B', c: 'b' }, { n: 'T', c: 't' }, { n: 'myeloid', c: 'm' }];
  const BASE = -2.35, BW = 0.95, DX = 1.7, HALF = 2.65, TH = 1.6, K = 1.55, LW = 1.5, LH = 4.75;
  const GX = [-4.7, 4.7];
  const CELLS = [[1.1, 2.7, 0.8], [1.3, 1.0, 1.45]]; // illustrative rule scores, not data
  const SLICES = [[0.42, 0.34], [0.46, 0.3], [0.31, 0.38]]; // three marker slices per lineage (illustrative)
  const CALL = CELLS.map((h) => { let b = 0; h.forEach((v, i) => { if (v > h[b]) b = i; }); return h[b] >= TH ? b : -1; });
  const NTH = 27;
  const groups = [0, 1].map((g) => {
    const grp = add(ctx.group());
    const base = ctx.line([[-HALF + 0.2, 0, 0], [HALF - 0.2, 0, 0]], { color: 'faint', width: 1.2 });
    grp.add(base);
    const bars = LIN.map((l, i) => {
      const b = ctx.box(BW, 1, { color: l.c, radius: 0 });
      b.position.set((i - 1) * DX, 0.5, 0.1);
      grp.add(b);
      const cuts = [0, 1].map(() => { const q = ctx.line([[-BW / 2, 0, 0.3], [BW / 2, 0, 0.3]], { color: 'card', width: 1.4 }); q.position.x = (i - 1) * DX; grp.add(q); return q; });
      let ghost = null;
      if (g) { ghost = ctx.line([[-BW / 2 - 0.16, 0, 0.5], [BW / 2 + 0.16, 0, 0.5]], { color: 'muted', width: 2 }); ghost.position.x = (i - 1) * DX; grp.add(ghost); }
      return { b, cuts, ghost };
    });
    const thPts = Array.from({ length: NTH }, (_, j) => new THREE.Vector3(-HALF + (2 * HALF * j) / (NTH - 1), TH, 0.6));
    const th = ctx.line(thPts, { color: 'ink', width: 1.5, dashed: [6, 4] });
    grp.add(th);
    return { grp, base, bars, thPts, th, key: null };
  });
  const titles = [ // left-aligned over each chart, clear of the loader's play/pause pill (top right)
    text('fixed rule', GX[0] - HALF + 0.2, 3.95, { anchor: 'left', weight: 600, size: 13 }),
    text('ANM', GX[1] - HALF + 0.2, 3.95, { anchor: 'left', weight: 600, size: 13, color: 'accent' }),
  ];
  const names = [];
  GX.forEach((gx) => LIN.forEach((l, i) => names.push(text(l.n, gx + (i - 1) * DX, BASE - 0.16, { anchor: 'top', color: l.c, weight: 600, size: 11.5 }))));
  const thLab = [
    text('θ', -1.93, BASE + TH, { anchor: 'left', size: 13 }),
    text('θ', 1.93, BASE + TH, { anchor: 'right', size: 13 }),
  ];
  const lens = add(ctx.box(LW, LH, { color: 'accent', stroke: 'accent', strokeWidth: 1.3, radius: 0.12 }));
  lens.position.set(0, BASE - 0.1 + LH / 2, 2);
  const lensLab = text('× G·c', 0, BASE - 0.1 + LH + 0.08, { anchor: 'bottom', color: 'accent', weight: 600, size: 12.5 });
  const CY = [2.35, 1.65, 0.95], CT = [0.114, 0.126, 0.158], XL = -1.05;
  const checks = ['same ranking', 'same bar', 'same call'].map((s, i) => ({
    tick: add(ctx.line([[0, 0, 5], [1, 0, 5], [2, 1, 5]], { color: 'good', width: 2 })),
    l: text(s, XL, CY[i], { anchor: 'left', size: 12 }),
  }));
  const chips = GX.map((gx) => text('', gx, -3.62, { weight: 600, size: 12.5 }));
  const eq = text('=', 0, -3.62, { size: 15, color: 'muted', weight: 600 });

  /* ── beat 2: every cell through both deciders ───────────────────────────────────────── */
  const DXC = -4.0, LY = 1.15, DW = 2.0, DH = 1.2, SI = 0.2, RAIL = 0.78, CUPW = 1.4, CUPH = 0.6;
  const BX = [-0.7, 1.35, 3.4, 5.45], SRC = [-6.5, 0];
  const BINS = [...LIN, { n: 'no call', c: 'muted' }];
  const deciders = [0, 1].map((g) => {
    const sg = g ? -1 : 1, y0 = sg * LY;
    const b = add(ctx.box(DW, DH, { color: 'soft', stroke: g ? 'accent' : 'muted', strokeWidth: 1.5, radius: 0.16 }));
    b.position.set(DXC, y0, 0.5);
    const l = text(g ? 'ANM' : 'fixed rule', DXC, y0 + sg * (DH / 2 + 0.12), { anchor: g ? 'top' : 'bottom', weight: 600, color: g ? 'accent' : 'ink' });
    const feed = add(ctx.line([[SRC[0] + 0.3, sg * 0.18, 0.2], [DXC - DW / 2, y0, 0.2]], { color: 'faint', width: 1.2 }));
    const rail = add(ctx.line([[DXC + DW / 2, y0, 0.2], [DXC + DW / 2 + 0.6, y0 + sg * RAIL, 0.2], [BX[3], y0 + sg * RAIL, 0.2]], { color: 'faint', width: 1.2 }));
    const cups = BINS.map((bn, i) => {
      const x = BX[i], yo = y0 + sg * CUPH / 2, yc = y0 - sg * CUPH / 2;
      return add(ctx.line([[x - CUPW / 2, yo, 0.4], [x - CUPW / 2, yc, 0.4], [x + CUPW / 2, yc, 0.4], [x + CUPW / 2, yo, 0.4]], { color: bn.c, width: 1.8 }));
    });
    return { b, l, feed, rail, cups };
  });
  const binLab = BINS.map((bn, i) => text(bn.n, BX[i], 0, { color: bn.c, weight: 600, size: 11.5, bg: 'card', bgOpacity: 1, pad: 3, order: 12 }));
  const srcRing = add(ctx.dot([SRC[0], SRC[1], 3], { px: 6, color: 'muted', hollow: true, ring: 0.34 }));
  const srcLab = text('cells', -7.9, 0, { anchor: 'left', size: 11.5, color: 'muted' });

  // Particles: one per cell, two copies (rule lane, ANM lane) on mirrored paths of equal length.
  const NPT = 34, V = 8.5, E0 = 0.285, E1 = 0.64;
  const rnd = ctx.rand(11);
  const parts = Array.from({ length: NPT }, (_, j) => {
    const r0 = rnd(), bin = r0 < 0.24 ? 0 : r0 < 0.55 ? 1 : r0 < 0.83 ? 2 : 3;
    const x = BX[bin] + (rnd() - 0.5) * 0.8, jy = (rnd() - 0.5) * 0.18;
    const paths = [1, -1].map((sg) => {
      const y0 = sg * LY;
      return [[SRC[0], 0], [DXC - DW / 2, y0], [DXC + DW / 2, y0], [DXC + DW / 2 + 0.6, y0 + sg * RAIL], [x, y0 + sg * RAIL], [x, y0 - sg * (0.13 - jy)]];
    });
    const cum = [0];
    for (let k = 1; k < paths[0].length; k++) cum.push(cum[k - 1] + Math.hypot(paths[0][k][0] - paths[0][k - 1][0], paths[0][k][1] - paths[0][k - 1][1]));
    return {
      bin, paths, cum, total: cum[cum.length - 1],
      e: E0 + ((E1 - E0) * (j + 0.35 * (rnd() - 0.5))) / (NPT - 1),
      dots: [0, 1].map(() => add(ctx.dot([0, 0, 3], { px: 3.4, color: 'ink', opacity: 0 }))),
      tie: add(ctx.line([[x, LY - CUPH / 2, 0.3], [x, -LY + CUPH / 2, 0.3]], { color: 'muted', width: 1, opacity: 0 })),
    };
  });
  const walk = (pts, cum, d, out) => {
    let k = 1;
    while (k < cum.length - 1 && cum[k] < d) k++;
    const f = clamp((d - cum[k - 1]) / Math.max(1e-6, cum[k] - cum[k - 1]));
    out.set(lerp(pts[k - 1][0], pts[k][0], f), lerp(pts[k - 1][1], pts[k][1], f), 3);
  };

  // Question tabs and the counter.
  const TB0 = 0.305, TD = 0.14, QC = -1.0, QY = 3.8, QX = [0, 0, 0]; // row centred left of the play/pause pill
  const qPre = text('question', 0, QY, { anchor: 'left', size: 11.5, color: 'muted' });
  const pills = QUESTIONS.map((q) => ({
    b: add(ctx.box(1, 1, { color: null, stroke: 'line', strokeWidth: 1.2, radius: 0.1 })),
    l: text(q, 0, QY, { size: 12 }),
  }));
  const hi = add(ctx.box(1, 1, { color: 'accent', stroke: 'accent', strokeWidth: 1.4, radius: 0.1 }));
  hi.position.set(0, QY, 0.6);
  const cntL = text('', -0.25, -3.62, { anchor: 'right', mono: true, size: 12 });
  const cntR = text(`differences ${NUMBERS.differences}`, 0.35, -3.62, { anchor: 'left', mono: true, size: 12, weight: 600, color: 'good' });

  /* ── beat 3: not better accuracy, an audit record ───────────────────────────────────── */
  const CYC = 0.35;
  const plate = add(ctx.group());
  plate.position.set(0, CYC, 4);
  const plateBox = ctx.box(1, 1, { color: 'soft', stroke: 'faint', strokeWidth: 1.3, radius: 0.12 });
  plate.add(plateBox);
  const strike = ctx.line([[-1, 0, 0.3], [1, 0, 0.3]], { color: 'muted', width: 2 });
  plate.add(strike);
  const plateLab = text('better accuracy', 0, CYC, { size: 14, color: 'faint' });
  const card = add(ctx.group());
  card.position.set(0, CYC, 4);
  const cardBox = ctx.box(1, 1, { color: 'accent', stroke: 'accent', strokeWidth: 1.6, radius: 0.16 });
  card.add(cardBox);
  const cardTitle = text('audit record', 0, 0, { size: 13.5, weight: 600, color: 'accent' });
  const rows = ['written question', 'per-marker reason', 'flip distance'].map((s) => ({
    dot: add(ctx.dot([0, 0, 6], { px: 3.4, color: 'accent' })),
    l: text(s, 0, 0, { anchor: 'left', size: 12.5 }),
  }));
  const note = text('not better accuracy', 0, 0, { anchor: 'top', size: 11.5, color: 'muted' });

  /* ── pixel-sized details (ticks, dots), refreshed when the scale changes ─────────────── */
  let lastPpu = 0;
  function layoutPx(ppu) {
    const p = (v) => v / ppu;
    checks.forEach((c, i) => {
      const x0 = XL - p(19), y = CY[i];
      c.tick.setPoints([[x0, y + p(0.5), 5], [x0 + p(4), y - p(4), 5], [x0 + p(12), y + p(5), 5]]);
    });
    const dpx = clamp(0.085 * ppu, 2.6, 4.6);
    parts.forEach((P) => P.dots.forEach((d) => d.setPx(dpx)));
    srcRing.setPx(clamp(0.13 * ppu, 5, 8));
  }

  function update(t) {
    const u = ctx.loopT(t, PERIOD), ppu = ctx.ppu(), ts = ctx.textScale, p = (v) => v / ppu;
    if (ppu !== lastPpu) { lastPpu = ppu; layoutPx(ppu); }

    /* beat 1 */
    const aIn = seg(u, 0, 0.03), aLab = aIn * (1 - seg(u, 0.25, 0.27));
    const z = seg(u, 0.255, 0.3), gv = aIn * (1 - seg(u, 0.73, 0.755));
    const rise = seg(u, 0.01, 0.05, ease.outCubic), m = seg(u, 0.19, 0.22);
    const lx = lerp(1.25, 8.5, seg(u, 0.055, 0.135, ease.inOutSine));
    const cover = (x) => ease.smooth((lx - x) / LW + 0.5);
    const on1 = seg(u, 0.145, 0.158) * (1 - seg(u, 0.183, 0.193));
    const on2 = seg(u, 0.222, 0.235) * (1 - seg(u, 0.25, 0.27));
    const ci = u < 0.205 ? 0 : 1, win = CALL[ci], on = ci ? on2 : on1;
    groups.forEach((G, g) => {
      const s = lerp(1, SI, z);
      G.grp.position.set(lerp(GX[g], DXC, z), lerp(BASE, g ? -LY - 0.32 : LY - 0.22, z), 1);
      G.grp.scale.set(s, s, 1);
      G.base.setOpacity(gv);
      G.bars.forEach((B, i) => {
        const k = g ? 1 + (K - 1) * cover(GX[1] + (i - 1) * DX) : 1;
        const h0 = rise * lerp(CELLS[0][i], CELLS[1][i], m), h = Math.max(1e-3, h0 * k);
        B.b.scale.y = h;
        B.b.position.y = h / 2;
        B.b.setOpacity(gv * (i === win ? 1 : lerp(1, win < 0 ? 0.72 : 0.36, on)));
        const f = SLICES[i];
        B.cuts.forEach((q, j) => { q.position.y = h * (j ? f[0] + f[1] : f[0]); q.setOpacity(gv * (1 - seg(z, 0.2, 0.5)) * (h > 0.05 ? 1 : 0)); });
        if (B.ghost) {
          B.ghost.position.y = h0;
          B.ghost.setOpacity(0.9 * gv * clamp(((k - 1) / (K - 1)) * 3) * (1 - seg(u, 0.14, 0.155)));
        }
      });
      if (g) {
        const key = lx.toFixed(4);
        if (G.key !== key) {
          G.key = key;
          G.thPts.forEach((v) => { v.y = TH * (1 + (K - 1) * cover(GX[1] + v.x)); });
          G.th.setPoints(G.thPts);
        }
      }
      G.th.setOpacity(gv).setWidth(win < 0 ? lerp(1.5, 2.4, on) : 1.5);
    });
    const kL = 1 + (K - 1) * cover(GX[1] - HALF);
    thLab[0].setOpacity(aLab);
    thLab[1].setText(kL > 1 + (K - 1) * 0.5 ? 'G·c·θ' : 'θ').setOpacity(aLab);
    thLab[1].position.y = BASE + TH * kL;
    titles.forEach((l) => l.setOpacity(aLab));
    names.forEach((l) => l.setOpacity(aLab));
    const lo = seg(u, 0.05, 0.062) * (1 - seg(u, 0.125, 0.138));
    lens.position.x = lx;
    lens.setOpacity(0.1 * lo);
    lens.outline.setOpacity(0.75 * lo);
    lensLab.position.x = Math.min(lx, 7.85 - lensLab.scale.x / 2);
    lensLab.setOpacity(lo * (1 - seg(lx, 6.4, 7.4)));
    checks.forEach((c, i) => { const o = seg(u, CT[i], CT[i] + 0.012) * aLab; c.tick.setOpacity(o); c.l.setOpacity(o); });
    chips.forEach((c) => {
      c.setText(win >= 0 ? `call ${LIN[win].n}` : 'no call');
      tint(c, win >= 0 ? LIN[win].c : 'muted');
      c.setOpacity(on * aIn);
    });
    eq.setOpacity(on * aIn);

    /* beat 2 */
    const bIn = seg(u, 0.28, 0.315) * (1 - seg(u, 0.73, 0.755));
    deciders.forEach((D) => {
      D.b.setOpacity(bIn);
      D.l.setOpacity(bIn);
      D.feed.setOpacity(bIn);
      D.rail.setOpacity(bIn);
      D.cups.forEach((c) => c.setOpacity(0.9 * bIn));
    });
    binLab.forEach((l) => l.setOpacity(bIn));
    srcRing.setOpacity(bIn);
    srcLab.setOpacity(bIn);
    const sec = u * PERIOD;
    parts.forEach((P) => {
      const d = (sec - P.e * PERIOD) * V, land = (d - P.total) / V;
      const vis = d < 0 || land > 0.75 || bIn <= 0 ? 0 : bIn;
      const inside = d > P.cum[1] && d < P.cum[2];
      const op = vis * (inside ? 0 : 1) * (1 - clamp((land - 0.4) / 0.35)) * clamp(d / 0.25);
      P.dots.forEach((dt, g) => {
        if (op > 0) walk(P.paths[g], P.cum, Math.min(Math.max(d, 0), P.total), dt.position);
        tint(dt, d >= P.cum[2] ? BINS[P.bin].c : 'ink');
        dt.setOpacity(op);
      });
      P.tie.setOpacity(land >= 0 && vis > 0 ? 0.75 * vis * Math.min(1, land / 0.08) * (1 - clamp((land - 0.3) / 0.3)) : 0);
    });
    const qi = clamp(Math.floor((u - TB0) / TD), 0, 2);
    const hPx = 21 * ts, padPx = 9 * ts;
    const pw = pills.map((P) => P.l.scale.x + p(2 * padPx)), gap = p(10 * ts), lead = qPre.scale.x + p(6 * ts);
    let x = QC - (lead + pw.reduce((a, b) => a + b, 0) + 2 * gap) / 2;
    qPre.position.x = x;
    x += lead;
    pw.forEach((w, i) => { QX[i] = x + w / 2; x += w + gap; });
    pills.forEach((P, i) => {
      P.b.position.set(QX[i], QY, 0.5);
      P.l.position.x = QX[i];
      fitBox(P.b, pw[i], p(hPx), p(hPx / 2));
      P.b.setOpacity(bIn);
      tint(P.l, i === qi ? 'accent' : 'muted');
      P.l.setOpacity(bIn);
    });
    const s1 = seg(u, TB0 + TD - 0.012, TB0 + TD + 0.004), s2 = seg(u, TB0 + 2 * TD - 0.012, TB0 + 2 * TD + 0.004);
    hi.position.x = QX[0] + (QX[1] - QX[0]) * s1 + (QX[2] - QX[1]) * s2;
    fitBox(hi, lerp(lerp(pw[0], pw[1], s1), pw[2], s2), p(hPx), p(hPx / 2));
    hi.setOpacity(0.14 * bIn);
    hi.outline.setOpacity(bIn);
    qPre.setOpacity(bIn);
    const n = NUMBERS.heldOutCells * seg(u - (TB0 + qi * TD), 0.005, TD - 0.03, ease.linear);
    cntL.setText(`cells ${fmt(u < TB0 ? 0 : n).padStart(6, ' ')}`).setOpacity(bIn);
    cntR.setOpacity(bIn);

    /* beat 3 */
    const cOut = 1 - seg(u, 0.955, 0.99);
    const pIn = seg(u, 0.75, 0.77), f1 = seg(u, 0.795, 0.81, ease.inQuad), f2 = seg(u, 0.81, 0.825, ease.outQuad);
    const lw = plateLab.scale.x;
    plate.scale.set(Math.max(1e-3, 1 - f1), 1, 1);
    fitBox(plateBox, lw + p(40 * ts), p(34 * ts), p(9));
    plateBox.setOpacity(pIn * (1 - f1));
    const sk = lw.toFixed(4);
    if (strike.userData.key !== sk) { strike.userData.key = sk; strike.setPoints([[-lw / 2 - p(4), 0, 0.3], [lw / 2 + p(4), 0, 0.3]]); }
    strike.setProgress(seg(u, 0.768, 0.785)).setOpacity(pIn);
    plateLab.setOpacity(pIn * (1 - f1));
    const G = 25 * ts, rw = Math.max(cardTitle.scale.x, ...rows.map((r) => r.l.scale.x + p(16)));
    const W = rw + p(46 * ts), H = p(4 * G + 18 * ts);
    card.scale.set(Math.max(1e-3, f2), 1, 1);
    fitBox(cardBox, W, H, p(10));
    const cv = f2 * cOut;
    cardBox.setOpacity(0.07 * cv);
    cardBox.outline.setOpacity(cv);
    cardTitle.position.set(0, CYC + p(1.5 * G), 10);
    cardTitle.setOpacity(cv);
    const xl = -W / 2 + p(23 * ts);
    rows.forEach((r, i) => {
      const y = CYC + p((0.5 - i) * G), o = seg(u, 0.825 + 0.01 * i, 0.837 + 0.01 * i) * cOut;
      r.dot.position.set(xl + p(3), y, 6);
      r.dot.setOpacity(o);
      r.l.position.set(xl + p(13), y, 10);
      r.l.setOpacity(o);
    });
    note.position.set(0, CYC - H / 2 - p(10), 10);
    note.setOpacity(seg(u, 0.86, 0.875) * cOut);
  }

  return {
    scene, camera, period: PERIOD, still: 0.175 * PERIOD,
    update,
    resize() { lastPpu = 0; },
    dispose() { parts.length = 0; groups.length = 0; colorOf.clear(); sized.clear(); },
  };
}
