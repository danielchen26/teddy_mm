/*
 * fw-depth-thinning — anm-framework.html #teddy-depth ("Finding: a response to sequencing depth").
 * Push the source down: keep a fraction f of each cell's RNA counts, re-embed, follow it to the call.
 * 1 A cell's RNA is a cloud of count dots. A dial clicks through the page's f = 100%, 50%, 20%, 10%, 5%.
 *   Each dot carries one fixed random draw, so at every click a random subset blinks out (a nested
 *   drawing of the page's binomial thinning, x_g^(f) ~ Binomial(x_g, f)); dropped counts stay as faint ghosts.
 * 2 After each click the thinned cell goes again through frozen TEDDY (the block never changes) into its
 *   embedding (a ring of beads that shimmers slightly), through our head into the predicted panel
 *   proteins (bars that wobble; heights are illustrative, no values), into ANM, which issues the call
 *   again. The call chip keeps its colour; at 5% it flickers but holds: the pipeline bends, not breaks.
 * 3 Two gauges on a 0–1 scale show only what the page reports, at 100% and at 5% of counts:
 *   panel-protein Pearson 0.937 → 0.805 and lineage accuracy 0.985 → 0.917 (1,500 validation cells,
 *   prototype). The page gives no values for 50%, 20% or 10%, so during those clicks the gauges show
 *   no marker, only the dotted line joining the two endpoints; at 5% the marker slides between them.
 * Every pipeline number shown lives in NUMBERS (refresh it after the official-preprocessing rerun).
 * Layout is sized for a 16:9 phone stage (~358 px) first; labels keep fixed px offsets. 12 s loop.
 */
const NUMBERS = {
  depths: [1, 0.5, 0.2, 0.1, 0.05],     // f: share of each cell's RNA counts kept (page: f ∈ {1, 0.5, 0.2, 0.1, 0.05})
  pearson: { full: 0.937, low: 0.805 }, // panel-protein Pearson, all counts → 5% of counts
  lineage: { full: 0.985, low: 0.917 }, // lineage accuracy, all counts → 5% of counts
  cells: 1500,                          // validation cells (prototype)
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const D = NUMBERS.depths, NL = D.length, LAST = NL - 1;
  const pct = (f) => `${+(f * 100).toFixed(1)}%`;
  const fmt = (v) => String(v);

  /* ── layout (world units, 16 × 9 design frame) ── */
  const T = 1.75;                                   // pipeline row
  const X = { cell: -6.6, teddy: -3.85, ring: -1.55, head: 0.8, anm: 6.45 };
  const CR = 1.05, TW = 1.9, TH = 1.3, RR = 0.6, HW = 2.2, HH = 1.0, AW = 1.5, AH = 0.9;
  const NB = 9, BP = 0.24, BX0 = 2.55, BY = T - 0.55, BHMAX = 1.2, BX1 = BX0 + (NB - 1) * BP;
  const CHIPY = T - 1.3;
  const DIAL = [-6.6, -1.75], KR = 0.48, TR = 0.8;
  const GX0 = -3.6, GX1 = 7.3, GY = [-1.0, -3.05];
  const GX = (v) => lerp(GX0, GX1, v);

  /* ── timeline (fractions of the loop): one cycle per depth, then a reset ── */
  const C = [0.01, 0.17, 0.33, 0.49, 0.65], L = 0.16, R0 = 0.93, R1 = 0.975; // 5% result holds ~1.4 s
  const SLIDE0 = C[LAST] + 0.86 * L, SLIDE1 = C[LAST] + L;

  const add = (o) => { scene.add(o); return o; };
  const text = (s, o) => add(ctx.label(s, Object.assign({ size: 11.5, color: 'ink' }, o)));
  const NAME = { size: 12.5, weight: 600, color: 'ink' };
  const SUB = { size: 11.5, color: 'muted' };
  const circ = (cx, cy, r, n = 48, a0 = 0, a1 = 360) => Array.from({ length: n }, (_, k) => {
    const a = ((a0 + ((a1 - a0) * k) / (a1 - a0 >= 360 ? n : n - 1)) * Math.PI) / 180;
    return [cx + r * Math.cos(a), cy + r * Math.sin(a), 0];
  });
  const block = (x, y, w, h, fill, tint, stroke, r = 0.2) => {
    const b = add(ctx.box(w, h, { color: fill, opacity: tint, stroke, strokeWidth: 1.5, radius: r }));
    b.position.set(x, y, 0);
    if (b.outline) b.outline.setOpacity(1);
    return b;
  };
  const conn = (a, b) => add(ctx.arrow([a[0], a[1], 0.2], [b[0], b[1], 0.2], { color: 'faint', width: 1.5, head: 6 }));
  const mcv = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, size, weight = 500) => { // CSS px, as ctx.label rasterises it
    if (!mcv) return s.length * size * 0.6 * ctx.textScale;
    mcv.font = `${weight} ${(size * ctx.textScale).toFixed(2)}px ${family}`;
    return mcv.measureText(s).width + 3;
  };
  const textH = (size) => size * ctx.textScale * 1.22 + 2;

  /* ── 1 the cell: a cloud of count dots, each with one fixed random draw ── */
  add(ctx.line(circ(X.cell, T, CR), { color: 'ink', width: 1.5, closed: true }));
  const cellLab = text('RNA counts', { anchor: 'bottom' });
  const rnd = ctx.rand(7);
  const dots = [];
  for (let tries = 0; dots.length < 160 && tries < 6000; tries++) {
    const a = rnd() * Math.PI * 2, r = (CR - 0.13) * Math.sqrt(rnd());
    const x = X.cell + r * Math.cos(a), y = T + r * Math.sin(a);
    if (dots.some((q) => (q.x - x) ** 2 + (q.y - y) ** 2 < 0.095 * 0.095)) continue;
    const draw = rnd();
    let lv = -1; // first depth at which this count is not kept (−1 = kept even at 5%)
    for (let k = 1; k < NL; k++) if (draw >= D[k]) { lv = k; break; }
    dots.push({ x, y, lv, j: rnd(), jb: rnd(), d: add(ctx.dot([x, y, 0.4], { r: 0.04, color: 'muted' })) });
  }

  /* ── the depth dial ── */
  const ANG = D.map((_, k) => ((150 - (120 * k) / LAST) * Math.PI) / 180);
  add(ctx.line(circ(DIAL[0], DIAL[1], TR, 24, 158, 22), { color: 'line', width: 1.5 }));
  const ticks = ANG.map((a) => add(ctx.dot([DIAL[0] + TR * Math.cos(a), DIAL[1] + TR * Math.sin(a), 0.6], { px: 2.4, color: 'faint' })));
  add(ctx.dot([DIAL[0], DIAL[1], 0.5], { r: KR, color: 'card' }));
  const knob = add(ctx.line(circ(DIAL[0], DIAL[1], KR, 40), { color: 'ink', width: 1.6, closed: true }));
  const pointer = add(ctx.line([[0, 0], [1, 0]], { color: 'ink', width: 2.2 }));
  add(ctx.dot([DIAL[0], DIAL[1], 0.8], { px: 2.6, color: 'ink' }));
  const dialVal = text(pct(D[0]), { size: 14, weight: 600, anchor: 'top' });
  const dialSub = text('counts kept', Object.assign({ anchor: 'top' }, SUB));
  add(ctx.line([[DIAL[0], DIAL[1] + TR + 0.14, 0.1], [DIAL[0], T - CR - 0.07, 0.1]], { color: 'line', width: 1.5, dashed: [3, 3] }));

  /* ── 2 frozen TEDDY (never animated), with a lock badge ── */
  block(X.teddy, T, TW, TH, 'teddy', 0.12, 'teddy');
  text('TEDDY', Object.assign({}, NAME)).position.set(X.teddy, T, 5);
  const frozen = text('frozen', Object.assign({ anchor: 'top' }, SUB));
  const lock = add(ctx.group());
  lock.position.set(X.teddy - TW / 2 - 0.05, T + TH / 2 + 0.05, 1);
  lock.add(ctx.dot([0, 0, 0], { r: 0.25, color: 'card' }));
  lock.add(ctx.line(circ(0, 0, 0.25, 32), { color: 'teddy', width: 1.5, closed: true }));
  const lockBody = ctx.box(0.2, 0.15, { color: 'teddy', radius: 0.03 });
  lockBody.position.set(0, -0.045, 0.1);
  lock.add(lockBody, ctx.line(circ(0, 0.03, 0.062, 12, 0, 180), { color: 'teddy', width: 1.5 }));

  /* ── the embedding: a ring of beads (count illustrative) ── */
  const rb = ctx.rand(19);
  const beads = Array.from({ length: 36 }, (_, i) => ({
    a: Math.PI / 2 - (i / 36) * Math.PI * 2, g: 2 * rb() - 1, ph: rb() * Math.PI * 2,
    d: add(ctx.dot([0, 0, 0.4], { r: 0.045, color: 'teddy' })),
  }));
  const ringLab = text('embedding', Object.assign({ anchor: 'bottom' }, SUB));

  /* ── our head ── */
  block(X.head, T, HW, HH, 'soft', 1, 'muted');
  text('our head', Object.assign({}, NAME)).position.set(X.head, T, 5);

  /* ── predicted panel proteins: nine bars (illustrative heights, no values) ── */
  const BASE = [0.34, 0.52, 0.28, 0.9, 0.72, 0.8, 0.3, 0.44, 0.24];
  const AMP = [0, 0.04, 0.07, 0.1, 0.15];
  const rh = ctx.rand(31);
  const H = AMP.map((a) => BASE.map((b) => clamp(b + a * (2 * rh() - 1), 0.08, 0.98)));
  const bars = BASE.map((_, j) => ({ ph: rh() * Math.PI * 2, l: add(ctx.line([[0, 0], [0, 1]], { color: 'teddy', width: 6 })) }));
  add(ctx.line([[BX0 - 0.16, BY, 0.3], [BX1 + 0.16, BY, 0.3]], { color: 'faint', width: 1 }));
  const barsLab = text('predicted panel', Object.assign({ anchor: 'top' }, SUB));

  /* ── ANM and its call ── */
  const anm = block(X.anm, T, AW, AH, 'accent', 0.1, 'accent');
  text('ANM', Object.assign({}, NAME)).position.set(X.anm, T, 5);
  const chip = block(X.anm, CHIPY, 1.6, 0.6, 't', 0.16, 't', 0.3);
  const chipLab = text('call: T', { size: 11.5, weight: 600, color: 'ink' });

  /* ── connectors and the tokens that ride them ── */
  conn([X.cell + CR + 0.08, T], [X.teddy - TW / 2 - 0.06, T]);
  conn([X.teddy + TW / 2 + 0.05, T], [X.ring - RR - 0.1, T]);
  conn([X.ring + RR + 0.1, T], [X.head - HW / 2 - 0.06, T]);
  conn([X.head + HW / 2 + 0.05, T], [BX0 - 0.2, T]);
  conn([BX1 + 0.2, T], [X.anm - AW / 2 - 0.06, T]);
  const toChip = conn([X.anm, T - AH / 2 - 0.05], [X.anm, CHIPY + 0.36]);
  const trips = [
    ['muted', [X.cell + CR + 0.1, T], [X.teddy - TW / 2 - 0.06, T], 0.2, 0.32],
    ['teddy', [X.teddy + TW / 2 + 0.06, T], [X.ring - RR - 0.1, T], 0.36, 0.44],
    ['teddy', [X.ring + RR + 0.1, T], [X.head - HW / 2 - 0.06, T], 0.5, 0.57],
    ['teddy', [X.head + HW / 2 + 0.06, T], [BX0 - 0.2, T], 0.6, 0.66],
    ['teddy', [BX1 + 0.2, T], [X.anm - AW / 2 - 0.06, T], 0.72, 0.78],
    ['accent', [X.anm, T - AH / 2 - 0.06], [X.anm, CHIPY + 0.36], 0.8, 0.86],
  ].map(([color, a, b, s0, s1]) => ({ pts: [[a[0], a[1], 3], [b[0], b[1], 3]], s0, s1, d: add(ctx.dot([0, 0, 3], { px: 3.6, color, opacity: 0 })) }));

  /* ── 3 the two gauges: 0–1 track, the page's two endpoints, a dotted line between ── */
  const gauges = [['panel-protein Pearson', 'teddy', NUMBERS.pearson], ['lineage accuracy', 'accent', NUMBERS.lineage]].map(([name, col, v], i) => {
    const y = GY[i], xf = GX(v.full), xl = GX(v.low);
    add(ctx.line([[GX0, y, 0.2], [GX1, y, 0.2]], { color: 'line', width: 2 }));
    const g = {
      y, xf, xl, col,
      title: text(name, { size: 12, weight: 600, anchor: 'bottom-left' }),
      zero: text('0', { size: 10.5, color: 'faint', anchor: 'right' }),
      one: text('1', { size: 10.5, color: 'faint', anchor: 'left' }),
      tf: add(ctx.line([[0, 0], [0, 1]], { color: col, width: 1.5, opacity: 0.75 })),
      tl: add(ctx.line([[0, 0], [0, 1]], { color: col, width: 1.5, opacity: 0.75 })),
      link: add(ctx.line([[xl, y, 0.5], [xf, y, 0.5]], { color: col, width: 2, dashed: [2, 4] })),
      full: text(`${fmt(v.full)} at ${pct(D[0])}`, { size: 11.5, weight: 600, color: 'ink', anchor: 'bottom-right' }),
      low: text(`${fmt(v.low)} at ${pct(D[LAST])}`, { size: 11.5, weight: 600, color: 'muted', anchor: 'top-right' }),
      halo: add(ctx.dot([xf, y, 1], { px: 7.5, color: 'card' })),
      mark: add(ctx.dot([xf, y, 1.1], { px: 5, color: col })),
      act: -2,
    };
    return g;
  });
  const note = text(`only ${pct(D[0])} and ${pct(D[LAST])} reported`, { size: 11, color: 'muted', anchor: 'top-left' });
  const foot = text(`prototype · ${NUMBERS.cells.toLocaleString('en-US')} validation cells`, { size: 11, color: 'muted', anchor: 'bottom-left' });
  gauges.forEach((g, i) => { g.low.userData.t = `${fmt([NUMBERS.pearson, NUMBERS.lineage][i].low)} at ${pct(D[LAST])}`; });
  foot.userData.t = `prototype · ${NUMBERS.cells.toLocaleString('en-US')} validation cells`;

  const P = { dot: 2, bead: 2, ppu: 50 };
  function layout() {
    const ppu = (P.ppu = ctx.ppu()), px = (v) => v / ppu;
    P.dot = Math.max(0.042 * ppu, 1.5);
    P.bead = Math.max(0.045 * ppu, 1.5);
    dots.forEach((q) => q.d.setPx(P.dot));
    beads.forEach((b) => b.d.setPx(P.bead));
    bars.forEach((b) => b.l.setWidth(Math.max(0.14 * ppu, 3)));
    lock.scale.setScalar(clamp(13 / (0.5 * ppu), 1, 1.6));
    cellLab.position.set(X.cell, T + CR + px(5), 5);
    frozen.position.set(X.teddy, T - TH / 2 - px(4), 5);
    ringLab.position.set(X.ring, T + RR + px(6), 5);
    // call chip sized to its text; the bars' label stays clear of it
    const cw = Math.max(1.3, px(textW('call: T', 11.5, 600) + 18)), ch = Math.max(0.55, px(textH(11.5) + 6));
    chip.setSize(cw, ch, ch / 2);
    chipLab.position.set(X.anm, CHIPY, 5);
    toChip.set([X.anm, T - AH / 2 - 0.05, 0.2], [X.anm, CHIPY + ch / 2 + 0.06, 0.2]);
    trips[5].pts[1][1] = CHIPY + ch / 2 + 0.06;
    const bw = px(textW('predicted panel', 11.5)), bcx = (BX0 + BX1) / 2;
    barsLab.position.set(Math.min(bcx, X.anm - cw / 2 - 0.2 - bw / 2), BY - px(5), 5);
    // dial
    dialVal.position.set(DIAL[0], DIAL[1] - KR - px(5), 5);
    dialSub.position.set(DIAL[0], DIAL[1] - KR - px(5 + textH(14)), 5);
    // gauges
    gauges.forEach((g) => {
      g.title.position.set(GX0, g.y + px(7), 5);
      g.zero.position.set(GX0 - px(4), g.y, 5);
      g.one.position.set(GX1 + px(4), g.y, 5);
      g.tf.setPoints([[g.xf, g.y - px(5), 0.6], [g.xf, g.y + px(5), 0.6]]);
      g.tl.setPoints([[g.xl, g.y - px(5), 0.6], [g.xl, g.y + px(5), 0.6]]);
      g.full.position.set(g.xf + px(3), g.y + px(7), 5);
      g.low.position.set(g.xl + px(3), g.y - px(7), 5);
    });
    note.position.set(GX0, GY[0] - px(7), 5);
    // footnote: under the second gauge when it clears that gauge's 5% label, else at the bottom edge
    const g2 = gauges[1], lowL = g2.xl + px(3) - px(textW(g2.low.userData.t, 11.5, 600));
    const fits = GX0 + px(textW(foot.userData.t, 11)) + px(10) < lowL;
    const bottom = camera.position.y + camera.bottom / camera.zoom;
    foot.setAnchor(fits ? 'top-left' : 'bottom-left');
    foot.position.set(GX0, fits ? g2.y - px(7) : bottom + px(5), 5);
  }
  let alive = true;
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (alive) { layout(); ctx.requestRender(); } });

  /* ── state helpers ── */
  const sOf = (u, k) => (u - C[k]) / L;
  const levelAt = (u, lag) => { // index of the depth in force at u (switches `lag` after each click)
    if (u >= R0 + 0.02) return 0;
    let k = 0;
    for (let i = 1; i < NL; i++) if (u >= C[i] + lag) k = i;
    return k;
  };
  const back = (u) => seg(u, R0, R1, ease.inOutSine);
  const MARCH = 12 / (SLIDE0 - C[1]); // whole dash periods over the marching window: seamless
  let shownVal = -1;

  function update(t) {
    const u = ctx.loopT(t, PERIOD), c = ctx.colors, rs = back(u);

    /* dial: clicks one detent per cycle, returns during the reset */
    let dk = -LAST * rs;
    for (let k = 1; k < NL; k++) dk += seg(u, C[k], C[k] + 0.022, ease.outBack);
    const i0 = clamp(Math.floor(dk), 0, LAST - 1), a = lerp(ANG[i0], ANG[i0 + 1], dk - i0);
    pointer.setPoints([[DIAL[0], DIAL[1], 0.7], [DIAL[0] + KR * 0.8 * Math.cos(a), DIAL[1] + KR * 0.8 * Math.sin(a), 0.7]]);
    const lv = levelAt(u, 0.012);
    if (lv !== shownVal) { dialVal.setText(pct(D[lv])); shownVal = lv; }
    ticks.forEach((d, k) => { d.setColor(k === lv ? 'ink' : 'faint'); d.setPx(k === lv ? 3.4 : 2.4); });
    let click = 0;
    for (let k = 1; k < NL; k++) click = Math.max(click, pulse(u, C[k], C[k] + 0.03, 0.01));
    knob.setWidth(1.6 + 1.2 * click);

    /* dots: flash, then fade to a ghost when their depth arrives; all return in the reset */
    for (const q of dots) {
      let flash = 0, gone = 0;
      if (q.lv > 0) {
        const ud = C[q.lv] + L * (0.03 + 0.1 * q.j);
        flash = pulse(u, ud, ud + 0.024, 0.008);
        gone = seg(u, ud + 0.012, ud + 0.036) * (1 - seg(u, R0 + 0.03 * q.jb, R0 + 0.01 + 0.03 * q.jb));
      }
      q.d.material.color.copy(c.muted).lerp(c.ink, 0.7 * flash);
      q.d.setOpacity(lerp(0.85, 0.1, gone));
      q.d.setPx(P.dot * (1 + 0.45 * flash) * lerp(1, 0.7, gone));
    }

    /* tokens ride the connectors once per cycle */
    for (const r of trips) {
      let o = 0, p = 0;
      for (let k = 0; k < NL; k++) {
        const s = sOf(u, k);
        if (s > r.s0 && s < r.s1) { p = seg(s, r.s0, r.s1, ease.inOutSine); o = pulse(s, r.s0, r.s1, 0.015); }
      }
      r.d.position.copy(ctx.along(r.pts, p)); r.d.position.z = 3;
      r.d.setOpacity(o);
    }

    /* embedding: a brief shimmer as each thinned cell arrives, a slight lasting offset */
    let rk = -LAST * rs, sh = 0;
    for (let k = 1; k < NL; k++) rk += seg(u, C[k] + 0.44 * L, C[k] + 0.62 * L);
    for (let k = 0; k < NL; k++) sh = Math.max(sh, pulse(sOf(u, k), 0.42, 0.64, 0.05) * (k ? 0.5 + k / LAST / 2 : 0.25));
    for (const b of beads) {
      const r = RR + 0.018 * b.g * rk + 0.045 * sh * Math.sin(b.ph + 40 * u);
      b.d.position.set(X.ring + r * Math.cos(b.a), T + r * Math.sin(b.a), 0.4);
      b.d.setOpacity(1 - 0.35 * sh * (0.5 + 0.5 * Math.sin(b.ph * 3 + 55 * u)));
    }

    /* bars: wobble and settle to slightly different heights at each depth */
    let kb = 0, e = 1, w0 = 0;
    for (let k = 1; k < NL; k++) if (u >= C[k] + 0.62 * L) kb = k;
    if (kb > 0) {
      const ss = clamp((u - C[kb] - 0.62 * L) / (0.28 * L));
      e = ease.outCubic(ss); w0 = (1 - ss) * (0.04 + (0.1 * kb) / LAST) * Math.sin(Math.PI * Math.min(1, ss * 4));
    }
    const s0 = sOf(u, 0), w00 = s0 > 0.62 && s0 < 0.9 ? 0.03 * Math.sin(Math.PI * ((s0 - 0.62) / 0.28)) : 0;
    bars.forEach((b, j) => {
      let h = kb > 0 ? lerp(H[kb - 1][j], H[kb][j], e) : H[0][j];
      h += (w0 || w00) * Math.sin(b.ph + 38 * u);
      h = lerp(h, H[0][j], rs);
      const x = BX0 + j * BP;
      b.l.setPoints([[x, BY, 0.5], [x, BY + BHMAX * Math.max(0.03, h), 0.5]]);
    });

    /* ANM reads, the call is issued again; at 5% it flickers but holds its colour */
    let glow = 0, pop = 0, dip = 0;
    for (let k = 0; k < NL; k++) {
      const s = sOf(u, k);
      glow = Math.max(glow, pulse(s, 0.76, 0.92, 0.04));
      pop = Math.max(pop, pulse(s, 0.86, 0.98, 0.05));
      if (k === LAST) dip = pulse(s, 0.74, 0.88, 0.035);
    }
    anm.material.opacity = 0.1 + 0.14 * glow;
    const flick = dip * (0.5 + 0.5 * Math.sin(70 * u));
    chip.scale.setScalar(1 + 0.1 * pop - 0.06 * dip);
    chip.material.opacity = 0.16 + 0.12 * pop - 0.12 * flick;
    if (chip.outline) chip.outline.setOpacity(1 - 0.6 * flick);
    chipLab.setOpacity(1 - 0.45 * flick);

    /* gauges: marker only at the page's two endpoints; none in between */
    const mid = u >= C[1] + 0.01 && u < SLIDE1 - 0.01;
    const noteK = pulse(u, C[1], SLIDE0 + 0.01, 0.02);
    gauges.forEach((g) => {
      let x = g.xf, o = 1;
      if (u >= C[1] && u < SLIDE0) o = 1 - seg(u, C[1], C[1] + 0.02);
      else if (u >= SLIDE0 && u < R0) { x = lerp(g.xf, g.xl, seg(u, SLIDE0, SLIDE1, ease.inOutCubic)); o = seg(u, SLIDE0, SLIDE0 + 0.008); }
      else if (u >= R0 && u < R0 + 0.015) { x = g.xl; o = 1 - seg(u, R0, R0 + 0.013); }
      else if (u >= R0 + 0.015) o = seg(u, R0 + 0.017, R1);
      g.mark.position.x = x; g.halo.position.x = x;
      g.mark.setOpacity(o); g.halo.setOpacity(o);
      g.link.setOpacity(0.55 + 0.35 * noteK).setDashPhase(clamp(u - C[1], 0, SLIDE0 - C[1]), MARCH);
      const act = mid ? 0 : u >= SLIDE1 - 0.01 && u < R0 + 0.015 ? -1 : 1; // 1 full, −1 low, 0 neither
      if (act !== g.act) { g.full.setColor(act === 1 ? 'ink' : 'muted'); g.low.setColor(act === -1 ? 'ink' : 'muted'); g.act = act; }
    });
    note.setOpacity(0.6 + 0.4 * noteK);
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: 0.87 * PERIOD,
    update, resize: layout,
    onTheme() { gauges.forEach((g) => { g.act = -2; }); shownVal = -1; },
    dispose() { alive = false; dots.length = 0; beads.length = 0; bars.length = 0; trips.length = 0; gauges.length = 0; },
  };
}
