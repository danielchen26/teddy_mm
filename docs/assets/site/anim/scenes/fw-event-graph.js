/*
 * fw-event-graph — anm-framework.html, #ai-graph ("AI step 2: build the typed event graph").
 * Five log records slide one by one through the schema gate: admitted, repaired (missing id filled)
 * or rejected (carries the answer key). Each record that passes becomes a typed event node (colour =
 * type, +/− = sign, size = value) linked to its target answer: a thick edge of weight 1 and a thin
 * return edge of weight 0.15 (W_e→a = 1, W_a→e = 0.15 on the page; a small key shows both). Then steps tick: s' = ρ s + β b (every node keeps the same
 * fraction, events enter at their own step), then one hop along every edge. Answer rings rise with
 * their value; the tallest one above θ is the call. The update below is the page's linear rule with
 * illustrative ρ, κ, β, event values and θ; no values are shown. 12 s seamless loop.
 */
export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;

  // Layout (world units, 16 × 9 frame).
  const XC = -6.85, XG = -4.75, XE = -2.5, TOP = 2.45, BASE = -1.5, HMAX = 3.2, RA = 0.36, YK = -3.95;
  const ROWS = [1.7, 1.0, 0.3, -0.4, -1.1];
  const ANS = [{ n: 'A', x: 1.5 }, { n: 'B', x: 3.4 }, { n: 'C', x: 5.3 }];
  const CARDS = [
    { out: 'admitted', col: 'good', ev: { a: 0, sign: 1, v: 0.6, time: 0, type: 'accent' } },
    { out: 'repaired', why: 'id filled', col: 'warn', ev: { a: 0, sign: 1, v: 0.8, time: 1, type: 'train' } },
    { out: 'rejected', why: 'carries key', col: 'bad' },
    { out: 'admitted', col: 'good', ev: { a: 1, sign: 1, v: 0.4, time: 2, type: 'accent' } },
    { out: 'admitted', col: 'good', ev: { a: 2, sign: -1, v: 0.5, time: 3, type: 'train' } },
  ];
  const EV = CARDS.filter((c) => c.ev).map((c) => Object.assign(c.ev, { y: CARDS.indexOf(c) }));
  const NE = EV.length, NA = ANS.length;

  // The page's update on the event graph: s' = ρ s + β b ; s_next = s' + κ Σ_edges w s'_source.
  const RHO = 0.7, KAP = 0.4, BETA = 1, K = 7;
  const S = [], P = [];
  let s = new Array(NE + NA).fill(0);
  for (let k = 0; k < K; k++) {
    const p = s.map((x) => RHO * x);
    EV.forEach((e, i) => { if (e.time === k) p[i] += BETA * e.sign * e.v; });
    const n = p.slice();
    EV.forEach((e, i) => { n[NE + e.a] += KAP * 1 * p[i]; n[i] += KAP * 0.15 * p[NE + e.a]; });
    S.push(s); P.push(p); s = n;
  }
  S.push(s);
  let VE = 1e-9, VA = 1e-9;
  for (const st of S.concat(P)) st.forEach((x, i) => { if (i < NE) VE = Math.max(VE, Math.abs(x)); else VA = Math.max(VA, Math.abs(x)); });
  const fin = s.slice(NE).map((x, a) => ({ x, a })).sort((q, r) => r.x - q.x);
  const THETA = (fin[0].x + fin[1].x) / 2;
  const WIN = fin[0].x >= THETA ? fin[0].a : -1;
  const yOf = (x) => BASE + (x / VA) * HMAX;

  // Timing (fractions of the loop).
  const C0 = 0.04, CD = 0.066, T0 = 0.405, TD = 0.066, CALL = T0 + K * TD + 0.008, OUT0 = 0.955, OUT1 = 0.995;

  // Quadratic arc from a to b, bowed to the left of a→b by h × length.
  const bez = (a, b, h, n = 16) => {
    const dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy), mx = (a.x + b.x) / 2 - (dy / L) * h * L, my = (a.y + b.y) / 2 + (dx / L) * h * L;
    return Array.from({ length: n + 1 }, (_, j) => { const t = j / n, q = 1 - t; return [q * q * a.x + 2 * q * t * mx + t * t * b.x, q * q * a.y + 2 * q * t * my + t * t * b.y, a.z]; });
  };

  // Static frame: log slots, gate, headings, rails, zero line, θ, step dots.
  const add = (o, z) => { if (z != null) o.position.z = z; scene.add(o); return o; };
  const lab = (text, x, y, o) => { const l = ctx.label(text, Object.assign({ size: 12, color: 'muted' }, o)); l.position.set(x, y, 10); return add(l); };
  ROWS.forEach((y) => { const b = ctx.box(1.2, 0.5, { color: null, stroke: 'line', dashed: [3, 3], radius: 0.08 }); b.position.set(XC, y, 0.2); add(b); });
  const gate = add(ctx.box(0.34, 3.8, { color: 'soft', opacity: 0.92, stroke: 'muted', strokeWidth: 1.6, radius: 0.1 }));
  gate.position.set(XG, 0.3, 3);
  const flash = add(ctx.box(0.34, 3.8, { color: 'good', opacity: 0, radius: 0.1 }));
  flash.position.set(XG, 0.3, 3.1);
  lab('records', XC, TOP, { anchor: 'bottom' });
  lab('schema', XG, TOP, { anchor: 'bottom' });
  lab('events', XE, TOP, { anchor: 'bottom' });
  lab('answers', ANS[1].x, TOP, { anchor: 'bottom' });
  ANS.forEach((A) => {
    add(ctx.line([[A.x, BASE - 1.45, 0.3], [A.x, BASE + HMAX + 0.3, 0.3]], { color: 'line', width: 1.2 }));
    lab(A.n, A.x, BASE - 1.62, { anchor: 'top', color: 'ink', weight: 600 });
  });
  add(ctx.line([[0.8, BASE, 0.4], [6.0, BASE, 0.4]], { color: 'faint', width: 1 }));
  const thLine = add(ctx.line([[0.8, yOf(THETA), 0.5], [6.0, yOf(THETA), 0.5]], { color: 'muted', width: 1.4, dashed: [6, 4] }));
  lab('θ', 6.15, yOf(THETA), { anchor: 'left', size: 14 });
  lab('steps', -2.6, YK, { anchor: 'right', size: 11 });
  const ticks = Array.from({ length: K }, (_, k) => add(ctx.dot([-2.25 + k * 0.45, YK, 1], { px: 3, color: 'faint' })));
  // Edge-weight key: thick edge = 1 (event → answer), thin arc = 0.15 (answer → event).
  add(ctx.line([[-7.6, YK, 1], [-7.1, YK, 1]], { color: 'muted', width: 3, opacity: 0.6 }));
  lab('1', -6.95, YK, { anchor: 'left', size: 11, mono: true });
  add(ctx.line(bez({ x: -6.45, y: YK - 0.05, z: 1 }, { x: -5.95, y: YK - 0.05, z: 1 }, 0.35, 8), { color: 'muted', width: 1.1 }));
  lab('0.15', -5.8, YK, { anchor: 'left', size: 11, mono: true });
  const status = lab('', XG, -2.05, { anchor: 'top', size: 13, weight: 600 });
  const why = lab('', XG, -2.5, { anchor: 'top', size: 11 });

  // Record cards: a plate with two text rows; a repaired card gets a filled patch.
  const cards = CARDS.map((c, i) => {
    const b = ctx.box(1.2, 0.5, { color: 'soft', stroke: 'muted', strokeWidth: 1.3, radius: 0.08 });
    const rows = [[-0.42, 0.28, 0.08], [-0.42, 0.05, -0.1]].map(([x0, x1, y]) => { const l = ctx.line([[x0, y, 0.05], [x1, y, 0.05]], { color: 'faint', width: 1.5 }); b.add(l); return l; });
    const patch = ctx.box(0.14, 0.3, { color: 'warn', radius: 0.03 });
    patch.position.set(0.42, 0, 0.06);
    b.add(patch);
    add(b).position.set(XC, ROWS[i], 1);
    return Object.assign(c, { b, rows, patch, t: C0 + i * CD });
  });

  // Event nodes, their edges (thick forward, thin bowed return), pulses.
  const events = EV.map((e) => {
    const x = XE, y = ROWS[e.y], r = 0.16 + 0.22 * e.v;
    const ring = add(ctx.dot([x, y, 4.2], { r, hollow: true, ring: 0.24, color: e.type }));
    const fill = add(ctx.dot([x, y, 4], { r: r * 0.82, color: e.type, opacity: 0 }));
    const sign = lab(e.sign > 0 ? '+' : '−', x - r - 0.12, y, { anchor: 'right', size: 14, weight: 600, color: 'ink' });
    const fwd = add(ctx.line([[0, 0, 2], [1, 0, 2]], { color: e.type, width: 3, opacity: 0.5 }));
    const back = add(ctx.line([[0, 0, 2.1], [1, 0, 2.1]], { color: 'muted', width: 1.1 }));
    const pf = add(ctx.dot([0, 0, 5], { px: 4, color: e.type, opacity: 0 }));
    const pb = add(ctx.dot([0, 0, 5], { px: 2.4, color: 'muted', opacity: 0 }));
    return Object.assign(e, { x, y0: y, r, ring, fill, sign, fwd, back, pf, pb, ptsF: null, ptsB: null });
  });

  // Answer rings (ANM teal once called).
  const rings = ANS.map((A) => ({
    ring: add(ctx.dot([A.x, BASE, 4.2], { r: RA, hollow: true, ring: 0.2, color: 'muted' })),
    fill: add(ctx.dot([A.x, BASE, 4], { r: RA * 0.8, color: 'accent', opacity: 0 })),
  }));
  const call = lab('call', 0, 0, { anchor: 'bottom', size: 13, weight: 600, color: 'accent' });

  function stateAt(u) {
    if (u < T0) return { v: S[0], k: -1, b: 0 };
    const x = (u - T0) / TD;
    if (x >= K) return { v: S[K], k: K, b: 0 };
    const k = Math.floor(x), f = x - k;
    const a = seg(f, 0, 0.3, ease.inOutSine), b = seg(f, 0.35, 0.95, ease.inOutSine);
    const v = S[k].map((q, i) => (a < 1 ? lerp(q, P[k][i], a) : lerp(P[k][i], S[k + 1][i], b)));
    return { v, k, b: f > 0.35 && f < 0.95 ? b : 0 };
  }

  return {
    scene, camera, period: PERIOD, still: 11.0,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const vin = seg(u, 0, 0.03, ease.inOutSine), vis = vin * (1 - seg(u, OUT0, OUT1, ease.inOutSine));
      const back0 = 1 - seg(u, OUT0, OUT1, ease.inOutSine);

      // Gate: each card slides through (or bounces back) once.
      let st = null;
      for (const c of cards) {
        const p = u - c.t;
        if (c.ev) {
          const m = seg(p, 0, 0.062, ease.inOutSine);
          const sc = lerp(1, 0.3, seg(p, 0.03, 0.062, ease.inOutSine));
          c.b.position.set(lerp(XC, XE, m), ROWS[CARDS.indexOf(c)], 1);
          c.b.scale.set(sc, sc, 1);
          c.o = vin * (1 - seg(p, 0.045, 0.064, ease.inOutSine));
          c.patch.setOpacity(c.why ? c.o * seg(p, 0.02, 0.03) : 0);
        } else {
          const m = seg(p, 0, 0.026, ease.inOutSine) - seg(p, 0.03, 0.062, ease.inOutSine);
          c.b.position.set(lerp(XC, XG - 0.17 - 0.6, m), ROWS[CARDS.indexOf(c)], 1);
          c.b.setStroke(p > 0.024 && u < OUT1 ? 'bad' : 'muted');
          c.o = vis * lerp(1, 0.45, seg(p, 0.03, 0.062));
          c.patch.setOpacity(0);
        }
        c.b.setOpacity(c.o);
        c.rows.forEach((l) => l.setOpacity(c.o));
        if (p > 0.016 && p < CD) st = c;
      }
      const sp = st ? ctx.pulse(u - st.t, 0.016, CD, 0.012) : 0;
      if (st) { status.setText(st.out).setColor(st.col); why.setText(st.why || ''); flash.setColor(st.col); }
      status.setOpacity(sp); why.setOpacity(sp);
      why.position.y = -2.05 - 17 / ctx.ppu();
      flash.setOpacity(st ? 0.55 * ctx.pulse(u - st.t, 0.018, 0.04, 0.01) : 0);

      // Graph state and answer heights.
      const { v, k, b } = stateAt(u);
      ticks.forEach((d, j) => { d.setColor(j <= k && k >= 0 && vis > 0.5 ? 'ink' : 'faint'); d.setOpacity(1); });
      rings.forEach((R, a) => {
        const y = lerp(BASE, yOf(v[NE + a]), back0);
        R.ring.position.y = y; R.fill.position.y = y;
        const won = a === WIN ? seg(u, CALL, CALL + 0.02) * vis : 0;
        R.ring.setColor(won > 0.5 ? 'accent' : 'muted');
        R.fill.setOpacity(0.28 * won);
        R.y = y;
      });
      const cw = WIN >= 0 ? rings[WIN] : null;
      call.setText(WIN >= 0 ? 'call' : 'no call').setOpacity(seg(u, CALL, CALL + 0.02) * vis);
      call.position.set(cw ? ANS[WIN].x : ANS[1].x, (cw ? cw.y : yOf(THETA)) + RA + 0.18, 10);
      thLine.setOpacity(0.6 + 0.4 * seg(u, CALL, CALL + 0.02) * vis);

      events.forEach((e, i) => {
        const c = cards[CARDS.findIndex((q) => q.ev === e)];
        const p = u - c.t, on = seg(p, 0.05, 0.066, ease.outCubic) * vis, grow = seg(p, 0.055, 0.1, ease.inOutSine);
        const E = { x: e.x, y: e.y0, z: 2 }, A = { x: ANS[e.a].x, y: rings[e.a].y, z: 2 };
        const dx = A.x - E.x, dy = A.y - E.y, L = Math.hypot(dx, dy), ux = dx / L, uy = dy / L;
        const e0 = { x: E.x + ux * e.r, y: E.y + uy * e.r, z: 2 }, a0 = { x: A.x - ux * RA, y: A.y - uy * RA, z: 2.1 };
        e.ptsF = [[e0.x, e0.y, 2], [a0.x, a0.y, 2]];
        e.ptsB = bez(a0, { x: e0.x, y: e0.y, z: 2.1 }, -0.14);
        e.fwd.setPoints(e.ptsF).setProgress(grow).setOpacity(0.5 * vis);
        e.back.setPoints(e.ptsB).setProgress(1, 1 - grow).setOpacity(0.8 * vis);
        e.ring.scale.setScalar(e.r * lerp(0.4, 1, on)); e.ring.setOpacity(on);
        e.sign.setOpacity(on);
        e.fill.setOpacity(on * (0.1 + 0.9 * clamp(Math.abs(v[i]) / VE)));
        const env = Math.min(1, b * 5, (1 - b) * 5);
        e.pf.position.copy(ctx.along(e.ptsF, b)); e.pf.position.z = 5;
        e.pf.setOpacity(b > 0 ? env * clamp(Math.abs(P[k][i]) / VE) * vis : 0);
        e.pb.position.copy(ctx.along(e.ptsB, b)); e.pb.position.z = 5;
        e.pb.setOpacity(b > 0 ? 0.7 * env * clamp(Math.abs(P[k][NE + e.a]) / VA) * vis : 0);
      });
    },
    dispose() { events.length = 0; cards.length = 0; rings.length = 0; },
  };
}
