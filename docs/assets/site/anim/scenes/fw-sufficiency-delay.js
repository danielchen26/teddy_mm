/*
 * fw-sufficiency-delay — anm-framework.html, #ai-state ("AI step 3: find and test the state").
 * A time plot seen at a slight slant: time runs right, value runs up, the translucent plane is "now".
 * Both curves solve the page's delay equation ẋ(t) = −a x(t) − k x(t−τ) + u(t) with the same input
 * (u = 0 here), from two histories u₁ (blue) and u₂ (orange) that differ over the past window [t−τ, t)
 * and meet at the present value. Candidate state 1 = the present value: equal for both, yet each future
 * point reads its own value τ earlier (dashed hooks), so the futures split and give two readouts — a
 * same-state violation (red). Candidate state 2 = the last τ of history (teal ribbons, as in fw-delay-ring):
 * the two states differ, so equal states give equal readouts again (green).
 * a, k, the history shapes and the time span are illustrative (Euler steps of τ/200); no values are shown.
 * 12 s seamless loop; the still frame is the red same-state violation.
 */
export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const YAW = -0.35, PITCH = 0.15, TX = -0.3, TY = -0.2;
  camera.position.set(TX + 100 * Math.sin(YAW) * Math.cos(PITCH), TY + 100 * Math.sin(PITCH), 100 * Math.cos(YAW) * Math.cos(PITCH));
  camera.lookAt(TX, TY, 0);
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, x, y, o = {}) => { const l = ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)); l.position.set(x, y, 0); return add(l); };

  /* ── the delay equation, solved once from each history ── */
  const TAU = 1, TEND = 1.15, N = 200, DT = TAU / N, A = 0.3, K = 1.8;
  const HIST = [(s) => -Math.sin((-Math.PI * s) / 2), (s) => 0.9 * Math.sin(-Math.PI * s)];
  const sols = HIST.map((h) => {
    const xs = [];
    for (let i = 0; i <= N; i++) xs.push(h(-TAU + i * DT));
    for (let i = 0, M = Math.round(TEND / DT); i < M; i++) {
      const x = xs[xs.length - 1], xd = xs[xs.length - 1 - N];
      xs.push(x + DT * (-A * x - K * xd));
    }
    return xs;
  });
  const NOWX = -1.2, TW = 4.4, YC = -0.2, VS = 1.9, AY = -3.0, PTOP = 2.45, AX = -6.6;
  const X = (s) => NOWX + s * TW, Y = (v) => YC + VS * v;
  const idx = (s, n) => { const f = (s + TAU) / DT, i = Math.max(0, Math.min(n - 2, Math.floor(f))); return [i, ctx.clamp(f - i)]; };
  const P = (c, s, z = 0) => { const xs = sols[c], [i, k] = idx(s, xs.length); return [X(s), Y(ctx.lerp(xs[i], xs[i + 1], k)), z]; };

  /* ── static frame: past window, now plane, axes ── */
  const win = add(ctx.box(TW, PTOP - AY, { color: 'soft', radius: 0.06, order: 0 }));
  win.position.set(NOWX - TW / 2, (AY + PTOP) / 2, -0.05);
  add(ctx.line([[NOWX - TW, AY, 0], [NOWX - TW, PTOP, 0]], { color: 'faint', width: 1, dashed: [3, 4], order: 1 }));
  const plane = add(ctx.box(2.0, PTOP - AY, { color: 'ink', opacity: 0.06, radius: 0.02, order: 1 }));
  plane.rotation.y = Math.PI / 2;
  plane.position.set(NOWX, (AY + PTOP) / 2, 0);
  add(ctx.line([[NOWX, AY, -1], [NOWX, PTOP, -1], [NOWX, PTOP, 1], [NOWX, AY, 1]], { color: 'faint', width: 1, closed: true, order: 1 }));
  add(ctx.arrow([AX, AY, 0], [X(TEND) + 1.5, AY, 0], { color: 'faint', width: 1.4, head: 7, order: 2 }));
  add(ctx.arrow([AX, AY, 0], [AX, PTOP - 0.1, 0], { color: 'faint', width: 1.4, head: 7, order: 2 }));
  lab('time', X(TEND) + 1.5, AY - 0.28, { anchor: 'top-right', color: 'muted', size: 11 });
  lab('value', AX, PTOP + 0.1, { anchor: 'bottom', color: 'muted', size: 11 });
  lab('past τ', NOWX - TW / 2, AY - 0.28, { anchor: 'top', color: 'muted', size: 11 });
  lab('now', NOWX, AY - 0.42, { anchor: 'top', color: 'ink', size: 11, weight: 600 });

  /* ── the two histories and their futures ── */
  const COL = ['t', 'teddy'], SUB = ['₁', '₂'];
  const curves = sols.map((xs, c) => add(ctx.line(xs.map((v, i) => [X(-TAU + i * DT), v * VS + YC, 0]), { color: COL[c], width: 2.6, order: 4 })));
  const frac = (c, s) => { const ln = curves[c], cum = ln.userData.animCum, [i, k] = idx(s, cum.length); return ctx.lerp(cum[i], cum[i + 1], k) / ln.length; };
  const ribbons = sols.map((xs) => add(ctx.line(xs.slice(0, N + 1).map((v, i) => [X(-TAU + i * DT), v * VS + YC, 0]), { color: 'accent', width: 11, opacity: 0.32, order: 3 })));
  const hooks = COL.map((col) => add(ctx.arrow([0, 0, 0], [1, 0, 0], { color: col, width: 1.4, head: 7, dashed: [4, 4], order: 5 })));
  const tails = COL.map((col) => add(ctx.dot([0, 0, 0], { px: 4, color: col, hollow: true, ring: 0.45, order: 7 })));
  const heads = COL.map((col) => add(ctx.dot([0, 0, 0], { px: 4.2, color: col, order: 7 })));
  const ends = COL.map((col, c) => add(ctx.dot(P(c, TEND), { px: 6.5, color: col, hollow: true, ring: 0.4, order: 7 })));
  const startLabs = COL.map((col, c) => lab('u' + SUB[c], X(-TAU) - 0.22, P(c, -TAU)[1], { anchor: 'right', color: col, weight: 600, size: 13 }));
  const endLabs = COL.map((col, c) => lab('o(u' + SUB[c] + ')', X(TEND) + 0.28, P(c, TEND)[1], { anchor: 'left', color: col, weight: 600 }));

  /* ── the candidate state: a bead at the present, later the last τ ── */
  const bead = add(ctx.dot(P(0, 0, 0.02), { px: 6, color: 'accent', order: 8 }));
  const halo = add(ctx.dot(P(0, 0, 0.02), { px: 8, color: 'accent', hollow: true, ring: 0.2, order: 8 }));
  const link = (color, c) => {
    const a = P(c, 0), e = P(c, TEND), to = [ctx.lerp(a[0], e[0], 0.965), ctx.lerp(a[1], e[1], 0.965), 0];
    return add(ctx.arrow(a, to, { color, width: 2, head: 8, bend: c ? 0.08 : -0.08, order: 6 }));
  };
  const red = [0, 1].map((c) => link('bad', c)), green = [0, 1].map((c) => link('good', c));

  /* ── one status line at the top; the phases never overlap in time ── */
  const status = [
    ['two different pasts', 'ink', 0.02, 0.19],
    ['state = present value', 'accent', 0.19, 0.32],
    ['each reads x(t − τ)', 'ink', 0.32, 0.58],
    ['same state, two readouts', 'bad', 0.58, 0.76],
    ['state = last τ', 'good', 0.76, 0.975],
  ].map(([text, color, a, b]) => ({ l: lab(text, TX + 0.2, 3.72, { size: 13, weight: 600, color }), a, b }));

  return {
    scene, camera, period: PERIOD, still: 0.68 * PERIOD,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const vis = 1 - ctx.seg(u, 0.93, 0.985, ease.inOutSine);
      const past = ctx.seg(u, 0.03, 0.17, ease.inOutSine), fut = ctx.seg(u, 0.32, 0.56, ease.inOutSine);
      const s = u < 0.32 ? -TAU + past * TAU : fut * TEND;
      const drawing = ctx.pulse(u, 0.03, 0.57, 0.02), hook = ctx.pulse(u, 0.32, 0.58, 0.03);
      const endK = ctx.seg(u, 0.54, 0.58) * vis, startK = ctx.seg(u, 0.03, 0.08) * vis;
      const ribK = ctx.seg(u, 0.74, 0.83, ease.inOutCubic);
      for (let c = 0; c < 2; c++) {
        curves[c].setProgress(u < 0.03 ? 0 : frac(c, s)).setOpacity(vis);
        heads[c].position.set(...P(c, s, 0.01)); heads[c].setOpacity(drawing);
        tails[c].position.set(...P(c, s - TAU, 0.01)); tails[c].setOpacity(hook);
        hooks[c].set(P(c, s - TAU), P(c, s), c ? -0.16 : 0.16).setOpacity(0.9 * hook);
        ends[c].setOpacity(endK); endLabs[c].setOpacity(endK); startLabs[c].setOpacity(startK);
        ribbons[c].setProgress(1, 1 - ribK).setOpacity(0.32 * vis);
        red[c].setProgress(ctx.seg(u, 0.58, 0.64, ease.outCubic)).setOpacity(ctx.pulse(u, 0.58, 0.8, 0.04));
        green[c].setProgress(ctx.seg(u, 0.8, 0.86, ease.outCubic)).setOpacity(ctx.seg(u, 0.8, 0.82) * vis);
      }
      const pop = ctx.seg(u, 0.16, 0.21, ease.outBack);
      bead.setPx(6 * Math.max(0, pop)); bead.setOpacity(ctx.seg(u, 0.16, 0.18) * (1 - ribK));
      const hb = ctx.pulse(u, 0.19, 0.33, 0.03), ring = ctx.loopT((u - 0.19) / 0.07, 1);
      halo.setPx(7 + 12 * ease.outCubic(ring)); halo.setOpacity(0.7 * hb * (1 - ring));
      status.forEach((st) => st.l.setOpacity(ctx.pulse(u, st.a, st.b, 0.025)));
    },
    dispose() { curves.length = ribbons.length = hooks.length = red.length = green.length = 0; },
  };
}
