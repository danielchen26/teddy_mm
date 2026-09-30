/*
 * idx-brings-depth-thinning — index.html #brings, card “3 · Response” (card width, 4:3, above .drop).
 * The card: thin each cell’s RNA, re-embed with TEDDY, let our head predict again. At 5% of counts
 * (median 305 UMIs) panel-protein Pearson goes 0.937 → 0.805 and lineage accuracy 0.985 → 0.917;
 * a declared per-depth trust; 1,500 validation cells (donor 18303), held-out site unused; Experiment 2
 * redesign prototype.
 * 1 A depth dial turns down through its stops, 100% → 5% of each cell’s RNA counts. Each count bead
 *   in the stream into frozen TEDDY has a fixed keep-threshold and shows only while the kept fraction
 *   is above it (binomial thinning), so the stream thins bead by bead.
 * 2 Downstream, TEDDY’s embedding (a ring) wobbles a little and the 9 panel proteins our head predicts
 *   (bars) jitter around their full-depth heights (faint caps). Ring, bar heights and jitter are
 *   illustrative, not data.
 * 3 Two gauges (0 to 1) fall with the dial: panel-protein Pearson and lineage accuracy. Only the page’s
 *   full-depth and 5% values are marked and printed; the needle path between them is illustrative, and
 *   the dial readout is blank between its end stops.
 * 4 A row of trust knobs, one per depth stop: each turns once and lights (ANM teal) as the dial passes
 *   its stop. They show that a trust is declared per depth, not its value (the page gives none).
 * 5 A corner tag carries the card’s caveat.
 * The real stage is 250–340 px wide (4:3) at every viewport, so labels are few and short and the
 * top-right corner stays clear for the play button. Every TEDDY-pipeline number is in NUMBERS below.
 * 10 s seamless loop; the still frame is the 5% hold.
 */
const NUMBERS = {
  depthFull: '100%',          // dial top stop, all counts kept
  depthLow: '5%',             // dial bottom stop
  umisLow: 305,               // median UMIs per cell at 5% of counts
  pearson: [0.937, 0.805],    // panel-protein Pearson, full depth → 5% of counts
  accuracy: [0.985, 0.917],   // lineage accuracy, full depth → 5% of counts
  cells: 1500,                // validation cells
  donor: '18303',
  /* kept fractions of the dial’s stops (anm-framework.html #teddy-depth); never printed here,
     they only set the bead density and the number of stops and knobs */
  stops: [1, 0.5, 0.2, 0.1, 0.05],
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const N = NUMBERS;
  const scene = new THREE.Scene();
  const W = 12, H = 9;
  const camera = ctx.orthoCamera({ width: W, height: H });
  const PERIOD = 10, TAU = Math.PI * 2;
  const NS = N.stops.length;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, x, y, o = {}) => {
    const l = ctx.label(text, Object.assign({ size: 11, color: 'muted' }, o));
    l.position.set(x, y, 10);
    return add(l);
  };
  const arc = (cx, cy, r, a0, a1, n, z = 0) => Array.from({ length: n + 1 }, (_, k) => {
    const a = a0 + ((a1 - a0) * k) / n;
    return [cx + r * Math.cos(a), cy + r * Math.sin(a), z];
  });
  const ring = (cx, cy, r, n, z = 0) => arc(cx, cy, r, 0, TAU, n, z).slice(0, n);
  const rnd = ctx.rand(305);

  /* ── layout (world units, 12 × 9 frame) ── */
  const XL = -5.85;                                  // left margin
  const Y1 = 2.45;                                   // pipeline row
  const D = { x: -5.12, y: Y1, r: 0.6 };             // depth dial
  const T = { x: -2.45, w: 1.55, h: 1.1 };           // frozen TEDDY
  const E = { x: -0.55, r: 0.42 };                   // embedding ring
  const HD = { x: 1.0, w: 1.1, h: 0.8 };             // our head
  const XB0 = 1.85, DXB = 0.3, YB = Y1 - HD.h / 2;   // the 9 bars (end before the play button)
  const HB = [0.5, 0.62, 0.36, 0.55, 0.3, 0.58, 0.45, 0.4, 0.52]; // illustrative; tops stay under the play button
  const YL = Y1 - 0.75;                              // labels under the ring and the bars
  const Y2 = -1.1, RG = 1.45;                        // gauges
  const YK = -2.85, KR = 0.26;                       // trust knobs
  const TH0 = -Math.PI / 4, DTH = (1.5 * Math.PI) / (NS - 1); // dial: 100% at lower right, 5% at lower left

  /* ── 1 · the depth dial and the RNA count stream ── */
  lab('RNA counts kept', XL, D.y + D.r + 0.28, { anchor: 'bottom-left', color: 'ink', weight: 600, size: 11.5 });
  add(ctx.dot([D.x, D.y, 0.2], { r: D.r, color: 'soft' }));
  add(ctx.line(ring(D.x, D.y, D.r, 48, 0.3), { color: 'muted', width: 1.5, closed: true }));
  const pointer = add(ctx.group());
  pointer.position.set(D.x, D.y, 1);
  pointer.add(ctx.line([[D.r * 0.12, 0, 0], [D.r * 0.8, 0, 0]], { color: 'ink', width: 2.4 }));
  pointer.add(ctx.dot([D.r * 0.8, 0, 0.1], { px: 2.8, color: 'ink' }));
  add(ctx.dot([D.x, D.y, 1.2], { px: 3, color: 'ink' }));
  const ticks = N.stops.map((_, k) => {
    const a = TH0 + k * DTH, x = D.x + (D.r + 0.14) * Math.cos(a), y = D.y + (D.r + 0.14) * Math.sin(a);
    return { on: add(ctx.dot([x, y, 0.5], { px: 2.4, color: 'ink', opacity: 0 })), base: add(ctx.dot([x, y, 0.4], { px: 2.2, color: 'faint' })) };
  });
  const RY = D.y - D.r - 0.2;
  const rFull = lab(N.depthFull, D.x, RY, { anchor: 'top', color: 'ink', weight: 600, mono: true, size: 11.5 });
  const rLow = lab(N.depthLow, D.x, RY, { anchor: 'top', color: 'ink', weight: 600, mono: true, size: 11.5, opacity: 0 });
  const umis = lab(`median ${N.umisLow} UMIs`, XL, RY - 0.56, { anchor: 'top-left', color: 'ink', size: 10.5 });

  // Count beads: 3 lanes × 5, each with a fixed keep-threshold r; a bead shows while r < kept fraction.
  const LANES = 3, PER = 5, CYC = 8;
  const XS0 = D.x + D.r + 0.28, XS1 = T.x - T.w / 2 - 0.08;
  const nB = LANES * PER;
  const rs = Array.from({ length: nB }, (_, i) => (i + 0.5) / nB);
  for (let i = nB - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [rs[i], rs[j]] = [rs[j], rs[i]]; }
  const swapTo = (rank, idx) => { const at = rs.indexOf((rank + 0.5) / nB); [rs[at], rs[idx]] = [rs[idx], rs[at]]; };
  swapTo(0, 1 * PER + 4); swapTo(1, 2 * PER + 1);   // the bead kept at 5% sits mid-stream in the still
  const beads = rs.map((r, i) => {
    const l = Math.floor(i / PER), j = i % PER;
    return { r, y: Y1 + (l - 1) * 0.2, ph: (j + 0.37 * l) / PER, d: add(ctx.dot([0, 0, 2], { px: 2, color: 'muted' })) };
  });
  const keptAt = (s) => {
    const i = Math.min(NS - 2, Math.floor(s)), k = s - i;
    return Math.exp(lerp(Math.log(N.stops[i]), Math.log(N.stops[i + 1]), clamp(k)));
  };

  /* ── 2 · frozen TEDDY, the embedding, our head, the 9 predicted panel proteins ── */
  const tb = add(ctx.box(T.w, T.h, { color: 'teddy', opacity: 0.14, stroke: 'teddy', strokeWidth: 1.5, radius: 0.16 }));
  tb.position.set(T.x, Y1, 0.5);
  if (tb.outline) tb.outline.setOpacity(1);
  lab('TEDDY', T.x, Y1, { color: 'ink', weight: 600, size: 11 });
  add(ctx.arrow([T.x + T.w / 2 + 0.05, Y1, 0.3], [E.x - E.r - 0.16, Y1, 0.3], { color: 'faint', width: 1.4, head: 6 }));
  add(ctx.line(ring(E.x, Y1, E.r, 40, 0.3), { color: 'line', width: 1.2, closed: true }));
  const NE = 14;
  const emb = Array.from({ length: NE }, (_, i) => ({
    a: (i / NE) * TAU + 0.1, k1: 2 + (i % 3), k2: 3 + ((i * 7) % 3), p1: rnd(), p2: rnd(),
    d: add(ctx.dot([0, 0, 1], { px: 2, color: 'teddy' })),
  }));
  lab('embedding', E.x, YL, { size: 10.5 });
  add(ctx.arrow([E.x + E.r + 0.14, Y1, 0.3], [HD.x - HD.w / 2 - 0.04, Y1, 0.3], { color: 'faint', width: 1.4, head: 6 }));
  const hb = add(ctx.box(HD.w, HD.h, { color: 'soft', stroke: 'muted', strokeWidth: 1.4, radius: 0.14 }));
  hb.position.set(HD.x, Y1, 0.5);
  lab('head', HD.x, Y1, { color: 'ink', weight: 600, size: 11 });
  const XB1 = XB0 + 8 * DXB;
  add(ctx.line([[HD.x + HD.w / 2, YB, 0.3], [XB1 + 0.16, YB, 0.3]], { color: 'line', width: 1.2 }));
  const bars = HB.map((h, i) => {
    const x = XB0 + i * DXB;
    return {
      x, h, k: 3 + (i % 4), ph: rnd(),
      bar: add(ctx.line([[x, YB, 1], [x, YB + h, 1]], { color: 'teddy', width: 4 })),
      cap: add(ctx.line([[x - 0.11, YB + h, 1.1], [x + 0.11, YB + h, 1.1]], { color: 'faint', width: 1.2, opacity: 0 })),
    };
  });
  lab('9 proteins', (XB0 + XB1) / 2, YL, { size: 10.5 });

  /* ── 3 · two gauges, 0 to 1, with the page’s two values marked ── */
  const thOf = (v) => Math.PI * (1 - v);
  const mono = { mono: true, weight: 600, size: 11.5, color: 'ink' };
  const gauges = [['panel-protein Pearson', N.pearson, -3], ['lineage accuracy', N.accuracy, 3]].map(([name, v, cx]) => {
    const track = add(ctx.line(arc(cx, Y2, RG, 0, Math.PI, 48, 0.3), { color: 'line', width: 5 }));
    const drop = add(ctx.line(arc(cx, Y2, RG, thOf(v[0]), thOf(v[1]), 16, 0.4), { color: 'teddy', width: 5 }));
    v.forEach((val, j) => {
      const a = thOf(val), r0 = RG + 0.1, r1 = RG + 0.32;
      add(ctx.line([[cx + r0 * Math.cos(a), Y2 + r0 * Math.sin(a), 0.5], [cx + r1 * Math.cos(a), Y2 + r1 * Math.sin(a), 0.5]],
        { color: j ? 'ink' : 'faint', width: 1.6 }));
    });
    const needle = add(ctx.group());
    needle.position.set(cx, Y2, 2);
    needle.add(ctx.line([[0, 0, 0], [RG * 0.84, 0, 0]], { color: 'ink', width: 2.2 }));
    add(ctx.dot([cx, Y2, 2.2], { px: 3.6, color: 'ink' }));
    const A = lab(v[0].toFixed(3), cx, Y2 - 0.4, Object.assign({ anchor: 'right' }, mono));
    const B = lab('→ ' + v[1].toFixed(3), cx, Y2 - 0.4, Object.assign({ anchor: 'left', opacity: 0 }, mono));
    lab(name, cx, Y2 - 0.95, { size: 10.5 });
    return { v, cx, track, drop, needle, A, B, wA: 0, wB: 0 };
  });

  /* ── 4 · one declared trust knob per depth stop (ANM) ── */
  const KT = 'declared trust per depth';
  lab(KT, XL, YK, { anchor: 'left', color: 'accent', weight: 600, size: 10.5 });
  const knobs = N.stops.map(() => {
    const g = add(ctx.group());
    const base = ctx.dot([0, 0, 0], { r: KR, hollow: true, ring: 0.2, color: 'faint' });
    const lit = ctx.dot([0, 0, 0.1], { r: KR, hollow: true, ring: 0.24, color: 'accent', opacity: 0 });
    const fill = ctx.dot([0, 0, 0.05], { r: KR * 0.8, color: 'accent', opacity: 0 });
    const spin = ctx.group();
    const notch0 = ctx.line([[0, KR * 0.18, 0.2], [0, KR * 0.7, 0.2]], { color: 'faint', width: 2 });
    const notch1 = ctx.line([[0, KR * 0.18, 0.3], [0, KR * 0.7, 0.3]], { color: 'accent', width: 2.2, opacity: 0 });
    spin.add(notch0, notch1);
    g.add(base, lit, fill, spin);
    return { g, lit, fill, spin, notch0, notch1 };
  });

  /* ── 5 · the card’s caveat ── */
  lab(`${N.cells.toLocaleString('en-US')} validation cells, donor ${N.donor};\nheld-out site not used; Experiment 2 prototype`,
    XL, -H / 2 + 0.12, { anchor: 'bottom-left', size: 10.5, align: 'left' });

  /* ── size-dependent pieces (text is in CSS px, constant on screen; shapes scale) ── */
  const mc = document.createElement('canvas').getContext('2d');
  const fcs = getComputedStyle(ctx.figure);
  const fam = fcs.fontFamily || 'system-ui, sans-serif';
  const monoFam = (fcs.getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace';
  const textW = (s, size, weight, m) => { mc.font = `${weight} ${size * ctx.textScale}px ${m ? monoFam : fam}`; return mc.measureText(s).width + 4; };
  let GAP = 0.1;
  function layout() {
    const w = ctx.width || 320, h = ctx.height || 240;
    const ppu = Math.min(w / W, h / H);
    const px = (k, lo, hi) => clamp(k * ppu, lo, hi);
    beads.forEach((b) => b.d.setPx(px(0.07, 1.6, 3.4)));
    emb.forEach((e) => e.d.setPx(px(0.065, 1.6, 3.2)));
    bars.forEach((b) => b.bar.setWidth(px(0.15, 3, 8)));
    gauges.forEach((g) => {
      const tw = px(0.13, 4, 9);
      g.track.setWidth(tw); g.drop.setWidth(tw);
      g.wA = textW(g.v[0].toFixed(3), mono.size, 600, true) / ppu;
      g.wB = textW('→ ' + g.v[1].toFixed(3), mono.size, 600, true) / ppu;
    });
    GAP = 3 / ppu;
    const x0 = XL + textW(KT, 10.5, 600) / ppu + 0.3 + KR, x1 = W / 2 - 0.2 - KR;
    knobs.forEach((k, i) => k.g.position.set(lerp(x0, x1, i / (NS - 1)), YK, 1));
  }

  /* ── timeline (fractions of the loop) ── */
  const S0 = 0.06, S1 = 0.54, H1 = 0.84, R1 = 0.95;
  const tauK = (k) => S0 + ((S1 - S0) * Math.acos(1 - (2 * k) / (NS - 1))) / Math.PI; // dial reaches stop k
  const progressAt = (u) => (u < H1 ? (u <= S0 ? 0 : u >= S1 ? 1 : ease.inOutSine((u - S0) / (S1 - S0))) : 1 - seg(u, H1, R1, ease.inOutCubic));

  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    const p = progressAt(u), s = p * (NS - 1), f = keptAt(s);
    const hold = seg(u, S1 - 0.01, S1 + 0.04, ease.inOutSine) * (1 - seg(u, H1, H1 + 0.04, ease.inOutSine));
    const full = 1 - seg(u, S0, S0 + 0.04, ease.inOutSine) + seg(u, R1 - 0.03, R1, ease.inOutSine);

    /* dial, its stops and readout (blank between the two end stops), the count stream */
    pointer.rotation.z = TH0 + s * DTH;
    ticks.forEach((tk, k) => { const on = clamp((s - k + 0.12) / 0.12); tk.on.setOpacity(on); tk.base.setOpacity(1 - on); });
    rFull.setOpacity(clamp(full));
    rLow.setOpacity(hold);
    umis.setOpacity(hold);
    for (const b of beads) {
      const q = (((b.ph + u * CYC) % 1) + 1) % 1;
      b.d.position.set(lerp(XS0, XS1, q), b.y, 2);
      const env = Math.min(1, q / 0.14, (1 - q) / 0.14);
      b.d.setOpacity(env * clamp((f - b.r) / (0.3 * b.r + 0.004)));
    }

    /* embedding wobble and head jitter grow as the counts thin (zero at full depth) */
    for (const e of emb) {
      const rr = E.r + p * 0.09 * Math.sin(TAU * (e.k1 * u + e.p1));
      const aa = e.a + p * 0.1 * Math.sin(TAU * (e.k2 * u + e.p2));
      e.d.position.set(E.x + rr * Math.cos(aa), Y1 + rr * Math.sin(aa), 1);
    }
    for (const b of bars) {
      const hh = Math.max(0.08, b.h + p * 0.1 * Math.sin(TAU * (b.k * u + b.ph)));
      b.bar.setPoints([[b.x, YB, 1], [b.x, YB + hh, 1]]);
      b.cap.setOpacity(0.9 * p);
    }

    /* gauges fall with the dial; only the two page values are printed */
    for (const g of gauges) {
      g.needle.rotation.z = thOf(lerp(g.v[0], g.v[1], p));
      g.drop.setProgress(p);
      const ax = lerp(g.cx + g.wA / 2, g.cx - (g.wA + GAP + g.wB) / 2 + g.wA, hold);
      g.A.position.x = ax;
      g.A.setOpacity(lerp(1, 0.55, Math.max(hold, seg(p, 0, 0.12))));   // the full-depth value, not the needle's
      g.B.position.x = ax + GAP;
      g.B.setOpacity(hold);
    }

    /* knobs: each turns once and lights as the dial passes its stop; all reset on the way back */
    const reset = 1 - seg(u, H1, R1, ease.inOutSine);
    knobs.forEach((k, i) => {
      const a = tauK(i), turn = seg(u, a, a + 0.07, ease.inOutCubic), on = turn * reset;
      k.spin.rotation.z = -TAU * turn;
      k.lit.setOpacity(on); k.fill.setOpacity(0.16 * on);
      k.notch1.setOpacity(on); k.notch0.setOpacity(1 - on);
    });
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: 0.7 * PERIOD,
    update, resize: layout,
    dispose() { beads.length = 0; emb.length = 0; bars.length = 0; gauges.length = 0; knobs.length = 0; },
  };
}
