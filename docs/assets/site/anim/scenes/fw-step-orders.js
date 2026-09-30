/*
 * fw-step-orders — anm-framework.html #math-steps, "A step is any real processing order".
 * Three real orders (clock time, iteration, layers) each advance the same state glyph one step at a
 * time: a teal orb trailed by 3 fading copies of its past states, so memory rides inside the state.
 * A fourth track shows the old panel order: the 9 markers (list order B, T, myeloid) enter one step
 * apart and every step multiplies every value by the retention factor, so earlier markers end smaller
 * ('hidden weight'). Then all 9 enter at step 0, shrink equally ('equal weight'). 12 s seamless loop.
 * Marker area ∝ value; the factor is illustrative only and no value is shown.
 */
export default function create(ctx) {
  const { THREE, ease, seg, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12, HALF = 6, STEPS = 8, T0 = 0.45, D = 0.44, MOVE = 0.3, RET = 0.85, RM = 0.33;
  const col = (k) => -3.1 + k * 0.95;
  const TX = -6.85, IX = -7.4, YS = [2.4, 1.2, 0], Y4 = -2.85, HY = 3.5; // top-right stays clear for the pause pill
  const add = (o) => { scene.add(o); return o; };
  const text = (s, x, y, o) => { const l = ctx.label(s, o); l.position.set(x, y, 5); return add(l); };
  const ring = (r, a0, a1, n = 36) => Array.from({ length: n + 1 }, (_, k) => {
    const a = a0 + ((a1 - a0) * k) / n;
    return [r * Math.cos(a), r * Math.sin(a), 0];
  });

  /* ── rails with one tick per step ── */
  const rail = (y) => {
    add(ctx.line([[col(0) - 0.45, y, 0], [col(STEPS) + 0.45, y, 0]], { color: 'line', width: 2 }));
    for (let k = 0; k <= STEPS; k++) add(ctx.line([[col(k), y - 0.02, 0], [col(k), y - 0.16, 0]], { color: 'faint', width: 1.2 }));
  };
  YS.forEach(rail);
  rail(Y4);
  add(ctx.line([[-7.6, -0.95, 0], [7.6, -0.95, 0]], { color: 'line', width: 1, dashed: [4, 4] }));

  /* ── icons + names of the three real orders ── */
  const clock = add(ctx.group());
  clock.position.set(IX, YS[0], 0.2);
  clock.add(ctx.line(ring(0.3, 0, Math.PI * 2, 40), { color: 'muted', width: 1.5, closed: true }));
  const hand = ctx.group();
  hand.add(ctx.line([[0, 0, 0], [0, 0.21, 0]], { color: 'ink', width: 1.8 }));
  clock.add(hand);
  clock.add(ctx.dot([0, 0, 0.1], { px: 1.8, color: 'ink' }));

  const iter = add(ctx.group());
  iter.position.set(IX, YS[1], 0.2);
  const A1 = (320 * Math.PI) / 180, E = [0.27 * Math.cos(A1), 0.27 * Math.sin(A1)], TG = [-Math.sin(A1), Math.cos(A1)];
  iter.add(ctx.line(ring(0.27, (50 * Math.PI) / 180, A1), { color: 'muted', width: 1.6 }));
  iter.add(ctx.arrow([E[0] - TG[0] * 0.02, E[1] - TG[1] * 0.02, 0], [E[0] + TG[0] * 0.1, E[1] + TG[1] * 0.1, 0], { color: 'muted', width: 1.6, head: 6 }));

  [-0.17, 0, 0.17].forEach((dx) => {
    const b = add(ctx.box(0.09, 0.5, { color: 'soft', stroke: 'muted', strokeWidth: 1.2, radius: 0.03 }));
    b.position.set(IX + dx, YS[2], 0.2);
  });
  ['clock time', 'iteration', 'layers'].forEach((s, i) => text(s, TX, YS[i], { size: 13, anchor: 'left' }));

  // network layers: slabs between step columns, drawn over the orb so it passes through them
  for (let k = 0; k < STEPS; k++) {
    const b = add(ctx.box(0.13, 0.78, { color: 'soft', opacity: 0.62, stroke: 'faint', strokeWidth: 1.2, radius: 0.04 }));
    b.position.set(col(k + 0.5), YS[2], 0.8);
  }

  /* ── the state glyph on each real track: orb + 3 past states ── */
  const GH = [0.5, 0.32, 0.18];
  const tracks = YS.map((y) => ({
    y,
    tail: add(ctx.line([[col(0), y, 0.3], [col(0), y, 0.3]], { color: 'accent', width: 3, opacity: 0 })),
    ghosts: GH.map(() => add(ctx.dot([col(0), y, 0.4], { r: 0.19, color: 'accent', opacity: 0 }))),
    orb: add(ctx.dot([col(0), y, 0.5], { r: 0.19, color: 'accent', opacity: 0 })),
  }));
  // header row: what the three tracks are, and a legend for the glyph (past states, then the state)
  text('real orders', TX, HY, { size: 11, color: 'muted', weight: 600, anchor: 'left' });
  [...GH].reverse().concat(1).forEach((o, i) => add(ctx.dot([col(0) + i * 0.34, HY, 0.4], { r: 0.12, color: 'accent', opacity: o })));
  text('state + its past', col(0) + 1.32, HY, { size: 11, color: 'muted', anchor: 'left' });

  // one shared step cursor through all four tracks
  const cursor = add(ctx.line([[0, YS[0] + 0.28, 0.1], [0, Y4 + 2 * RM + 0.12, 0.1]], { color: 'faint', width: 1, dashed: [3, 4], opacity: 0 }));

  /* ── track 4: the 9 markers in list order ── */
  const LIN = ['b', 'b', 'b', 't', 't', 't', 'm', 'm', 'm'];
  const markers = LIN.map((c, k) => add(ctx.dot([col(k), Y4 + RM, 0.6], { r: RM, color: c, opacity: 0 })));
  [['B', 'b', 1], ['T', 't', 4], ['myeloid', 'm', 7]].forEach(([s, c, k]) => text(s, col(k), Y4 - 0.3, { size: 11, color: c, weight: 600, anchor: 'top' }));
  text('9 markers', TX, Y4 + 0.3, { size: 13, anchor: 'left' });
  const sub = text('one per step', TX, Y4 - 0.22, { size: 11, color: 'muted', anchor: 'left' });
  const TAG = { size: 11, color: 'card', weight: 600, anchor: 'left', pad: 5, bgOpacity: 1, opacity: 0 };
  const badTag = text('hidden\nweight', col(STEPS) + 0.62, Y4 + 0.34, Object.assign({ bg: 'bad' }, TAG));
  const goodTag = text('equal\nweight', col(STEPS) + 0.62, Y4 + 0.34, Object.assign({ bg: 'good' }, TAG));
  let shownPhase = 0;

  return {
    scene, camera, period: PERIOD, still: 4.8,
    update(t) {
      const u = ctx.loopT(t, PERIOD), ph = u < 0.5 ? 0 : 1, tau = u * PERIOD - ph * HALF;
      const vis = seg(tau, 0, 0.3, ease.inOutSine) * (1 - seg(tau, 5.55, 5.95, ease.inOutSine));
      const stepStart = (j) => T0 + (j - 1) * D;
      let S = 0;
      for (let j = 1; j <= STEPS; j++) S += seg(tau, stepStart(j), stepStart(j) + MOVE, ease.inOutCubic);

      // real orders: the orb steps on, its past states follow one step behind each other
      for (const tr of tracks) {
        tr.orb.position.x = col(S);
        tr.orb.setOpacity(vis);
        tr.ghosts.forEach((g, i) => {
          g.position.x = col(Math.max(0, S - i - 1));
          g.setOpacity(vis * GH[i] * clamp(S - i));
        });
        tr.tail.setPoints([[col(Math.max(0, S - 3)), tr.y, 0.3], [col(S), tr.y, 0.3]]);
        tr.tail.setOpacity(vis * 0.3 * clamp(S));
      }
      hand.rotation.z = -(ph * STEPS + S) * (Math.PI / 8);
      iter.rotation.z = S * Math.PI * 2;
      cursor.position.x = col(S);
      cursor.setOpacity(vis * 0.8);

      // track 4: every step multiplies every value by the same factor; entry step decides the end size
      markers.forEach((m, k) => {
        const e = ph === 0 && k > 0
          ? seg(tau, stepStart(k), stepStart(k) + MOVE, ease.outCubic)
          : seg(tau, 0.08, 0.4, ease.outCubic);
        const age = ph === 0 ? Math.max(0, S - k) : S;
        const r = RM * Math.sqrt(Math.pow(RET, age));
        m.setRadius(r);
        m.position.y = Y4 + r + 0.85 * (1 - e);
        m.setOpacity(vis * e);
      });
      if (shownPhase !== ph) { shownPhase = ph; sub.setText(ph ? 'all at step 0' : 'one per step'); }
      sub.setOpacity(vis);
      const tag = seg(tau, 3.9, 4.25, ease.inOutSine) * (1 - seg(tau, 5.55, 5.95, ease.inOutSine));
      badTag.setOpacity(ph === 0 ? tag : 0);
      goodTag.setOpacity(ph === 1 ? tag : 0);
    },
    dispose() {
      tracks.length = 0;
      markers.length = 0;
    },
  };
}
