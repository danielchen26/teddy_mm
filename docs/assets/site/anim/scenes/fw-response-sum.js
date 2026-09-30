/*
 * fw-response-sum — anm-framework.html, #math-response ("The response along the real run").
 * Six state slabs s1..s6 sit on the real run (the tube). Pushes enter at steps 2, 4 and 5 and become
 * state changes; each lens between two slabs is that step's A, which turns and rescales every carried
 * change the same way. At the end the carried changes line up head-to-tail and are read on the readout
 * axis, so the output change is one stacked bar with one segment per push. The second half of the loop
 * replays with the step-4 push removed: its segment is missing from the sum (dashed ghost).
 * 12 s seamless loop (two 6 s passes). Angles and lengths are illustrative, not data; no numbers shown.
 */
export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const D2R = Math.PI / 180;

  const N = 6, X0 = -6.3, DX = 1.7, YT = -1.5, XR = 4.6, YTOP = 3.3;
  const xs = Array.from({ length: N }, (_, k) => X0 + k * DX);
  const SUB = '₁₂₃₄₅₆';
  // Entering slab k applies that step's A: turn by R[k] degrees and rescale by S[k] (illustrative).
  const R = [0, 10, 25, -20, 20, -15];
  const S = [1, 1, 1.15, 0.85, 1.1, 0.9];
  // Pushes at steps 2, 4, 5 (slab 1, 3, 4): colour, then angle and length right after entering.
  const P = [
    { s: 1, col: 'train', th: 30, len: 1.19 },
    { s: 3, col: 'teddy', th: 95, len: 1.11 },
    { s: 4, col: 'accent', th: 85, len: 1.33 },
  ];
  const DROPPED = 1; // the second pass replays without the step-4 push

  // Pass timing (p in 0..1 of one 6 s pass): the front reaches slab k at arrive(k), dwells DW, moves on.
  const P0 = 0.06, STEP = 0.1, DW = 0.03;
  const arrive = (k) => P0 + k * STEP;
  function front(p) {
    for (let k = 0; k < N - 1; k++) {
      if (p < arrive(k + 1)) return k + seg(p, arrive(k) + DW, arrive(k + 1), ease.inOutSine);
    }
    return N - 1;
  }
  // A push's state change when the front is at fractional slab f: every later A applied in turn;
  // the turn and rescale happen while the change passes through the lens (middle of the move).
  function carried(q, f) {
    if (f < q.s) return null;
    let ang = q.th, len = q.len;
    const k = Math.floor(f), fr = f - k;
    for (let j = q.s + 1; j <= k; j++) { ang += R[j]; len *= S[j]; }
    if (k + 1 < N && fr > 0) {
      const lam = seg(fr, 0.3, 0.7, ease.inOutSine);
      ang += R[k + 1] * lam;
      len *= Math.pow(S[k + 1], lam);
    }
    return { dx: len * Math.cos(ang * D2R), dy: len * Math.sin(ang * D2R) };
  }

  // State slabs and their names.
  const slabs = xs.map((x, k) => {
    const b = ctx.box(0.9, 5.3, { color: 'soft', opacity: 0.8, stroke: 'muted', strokeWidth: 1.2, radius: 0.12 });
    b.position.set(x, -0.25, 0);
    const lab = ctx.label('s' + SUB[k], { size: 12, color: 'muted', anchor: 'top' });
    lab.position.set(x, -3.2, 5);
    scene.add(b, lab);
    return b;
  });

  // The real run: a soft wide tube with a thin core line.
  const tA = X0 - 1.3, tB = xs[N - 1] + 0.95;
  scene.add(
    ctx.line([[tA, YT, 1], [tB, YT, 1]], { color: 'line', width: 9 }),
    ctx.line([[tA, YT, 1.1], [tB, YT, 1.1]], { color: 'faint', width: 1.6 }),
  );
  const runLab = ctx.label('real\nrun', { size: 12, color: 'muted', anchor: 'bottom-left' });
  runLab.position.set(tA - 0.05, YT + 0.24, 5);
  scene.add(runLab);

  // Lenses between slabs: A of the step being entered.
  const lenses = [];
  for (let k = 1; k < N; k++) {
    const x = (xs[k - 1] + xs[k]) / 2, pts = [];
    for (let i = 0; i < 40; i++) { const a = (i / 40) * Math.PI * 2; pts.push([x + 0.13 * Math.cos(a), YT + 0.42 * Math.sin(a), 2]); }
    const ring = ctx.line(pts, { color: 'muted', width: 1.5, closed: true });
    const lab = ctx.label('A', { size: 12, color: 'muted', weight: 600, anchor: 'top' });
    lab.position.set(x, YT - 0.55, 5);
    scene.add(ring, lab);
    lenses.push({ x, ring, lab });
  }

  // "Now" on the run.
  const dot = ctx.dot([xs[0], YT, 4], { px: 4.5, color: 'ink' });
  scene.add(dot);

  // Readout axis (zero = the output of the real run) and the sum.
  scene.add(
    ctx.line([[XR, YT, 1], [XR, 2.5, 1]], { color: 'muted', width: 1.4 }),
    ctx.line([[XR - 0.22, YT, 1], [XR + 0.22, YT, 1]], { color: 'muted', width: 1.4 }),
  );
  const readLab = ctx.label('readout', { size: 12, color: 'muted', anchor: 'top' });
  readLab.position.set(XR, YT - 0.3, 5);
  const sumLab = ctx.label('sum', { size: 14, weight: 600, color: 'ink', anchor: 'left' });
  const noLab = ctx.label('no push', { size: 12, color: 'muted', anchor: 'bottom' });
  noLab.position.set(xs[P[DROPPED].s], YTOP + 0.15, 5);
  // Second pass: an empty dashed slot the size of the removed push's segment, on top of the shorter sum.
  const ghostLen = carried(P[DROPPED], N - 1).dy;
  const ghost = ctx.box(0.2, ghostLen, { color: null, stroke: P[DROPPED].col, strokeWidth: 1.5, dashed: [4, 3], radius: 0.03 });
  scene.add(readLab, sumLab, noLab, ghost);

  // Per push: the drop (B injecting the push), the carried change, a faint footprint where it entered,
  // its name, a projection line to the readout and its segment of the bar.
  for (const q of P) {
    const X = xs[q.s], a = q.th * D2R;
    q.drop = ctx.arrow([X, YTOP, 3], [X, YTOP - 0.8, 3], { color: q.col, width: 2.2, head: 8 });
    q.vec = ctx.arrow([X, YT, 3.5], [X + 1, YT, 3.5], { color: q.col, width: 2.6, head: 9 });
    q.foot = ctx.arrow([X, YT, 2.5], [X + q.len * Math.cos(a), YT + q.len * Math.sin(a), 2.5], { color: q.col, width: 1.6, head: 7 });
    q.lab = ctx.label('δu' + SUB[q.s], { size: 13, weight: 600, color: q.col, anchor: 'bottom' });
    q.lab.position.set(X, YTOP + 0.15, 5);
    q.proj = ctx.line([[0, 0, 1.5], [1, 0, 1.5]], { color: 'faint', width: 1.2, dashed: [4, 4] });
    q.bar = ctx.line([[XR, 0, 2.2], [XR, 1, 2.2]], { color: q.col, width: 7 });
    scene.add(q.drop, q.vec, q.foot, q.lab, q.proj, q.bar);
  }

  return {
    scene, camera, period: PERIOD, still: 5.3,
    update(t) {
      const u = ctx.loopT(t, PERIOD), passB = u >= 0.5, p = (u % 0.5) * 2;
      const vis = seg(p, 0, 0.04, ease.inOutSine) * (1 - seg(p, 0.93, 0.99, ease.inOutSine));
      const f = front(p), k0 = Math.floor(f), fr = f - k0;
      const xf = lerp(xs[k0], xs[Math.min(k0 + 1, N - 1)], fr);
      dot.position.set(xf, YT, 4);
      dot.setOpacity(vis);
      slabs.forEach((b, k) => b.outline.setOpacity(0.35 + 0.65 * vis * ctx.smooth(clamp(1 - Math.abs(f - k)))));
      for (const L of lenses) {
        const h = vis * ctx.smooth(clamp(1 - Math.abs(xf - L.x) / (DX / 2)));
        L.ring.setOpacity(0.45 + 0.55 * h);
        L.lab.setOpacity(0.7 + 0.3 * h);
      }

      // Readout: changes line up head-to-tail at the last slab, then project onto the axis.
      const chain = seg(p, 0.61, 0.69);
      let cx = xs[N - 1], cy = YT, j = 0;
      P.forEach((q, i) => {
        const on = !(passB && i === DROPPED);
        const a = arrive(q.s), X = xs[q.s];
        const fall = seg(p, a - 0.075, a - 0.005, ease.inOutSine);
        const y = lerp(YTOP, YT + 0.85, fall);
        q.drop.set([X, y, 3], [X, y - 0.8, 3]);
        q.drop.setOpacity(on ? vis * seg(p, a - 0.075, a - 0.06) * (1 - seg(p, a - 0.012, a + 0.01)) : 0);
        q.lab.setOpacity(on ? vis * seg(p, a - 0.08, a - 0.05) : 0);
        q.foot.setOpacity(on ? 0.3 * vis * seg(p, a + 0.02, a + 0.06) : 0);
        const v = on ? carried(q, f) : null;
        if (!v) { q.vec.setProgress(0); q.proj.setProgress(0); q.bar.setProgress(0); return; }
        const tx = lerp(xf, cx, chain), ty = lerp(YT, cy, chain);
        q.vec.set([tx, ty, 3.5], [tx + v.dx, ty + v.dy, 3.5]);
        q.vec.setProgress(seg(p, a, a + 0.03, ease.outCubic));
        q.vec.setOpacity(vis);
        const hx = cx + v.dx, hy = cy + v.dy;
        q.proj.setPoints([[hx, hy, 1.5], [XR, hy, 1.5]]);
        q.proj.setProgress(seg(p, 0.69 + 0.03 * j, 0.74 + 0.03 * j));
        q.proj.setOpacity(vis);
        q.bar.setPoints([[XR, cy, 2.2], [XR, hy, 2.2]]);
        q.bar.setProgress(seg(p, 0.71 + 0.035 * j, 0.77 + 0.035 * j));
        q.bar.setOpacity(vis);
        cx = hx; cy = hy; j++;
      });
      sumLab.position.set(XR + 0.34, (YT + cy) / 2, 5);
      sumLab.setOpacity(vis * seg(p, 0.8, 0.85));
      noLab.setOpacity(passB ? vis * seg(p, arrive(P[DROPPED].s) - 0.08, arrive(P[DROPPED].s) - 0.05) : 0);
      ghost.position.set(XR, cy + ghostLen / 2, 2.1);
      ghost.setOpacity(passB ? 0.85 * vis * seg(p, 0.84, 0.9) : 0);
    },
    dispose() {
      // Every object was made with ctx helpers and added to the scene; the core disposes them all.
      P.forEach((q) => { q.drop = q.vec = q.foot = q.lab = q.proj = q.bar = null; });
      lenses.length = 0;
    },
  };
}
