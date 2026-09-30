/*
 * fw-forward-inverse — anm-framework.html #forward-inverse ("Forward: where lost. Inverse: is it enough?").
 * Split stage, 12 s seamless loop; each half plays for 6 s while the other rests, dimmed.
 * Left, Forward: the push / follow / read chain. The knob turns one tick, the state lattice deforms
 *   (teal), the gauge needle never moves, and a red ring on the lattice→gauge link marks where the
 *   response stops ("lost here"): the information arrived and the readout does not use it.
 * Right, Inverse: a box-shaped state space. Inputs x₁ and x₂ start far apart and land on one state z;
 *   from z two arrows reach different outputs A and B, and z cracks ("not enough"): the state is
 *   missing something the output depends on.
 * Illustrative layout only: no values are shown.
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const tmp = new THREE.Color();
  const mix = (a, b, k) => tmp.copy(ctx.colors[a]).lerp(ctx.colors[b], clamp(k));
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, x, y, o = {}) => {
    const s = ctx.label(text, Object.assign({ size: 12, color: 'muted' }, o));
    s.position.set(x, y, 6);
    return add(s);
  };
  const L = [], R = []; // [object, base opacity]: static pieces dimmed with their half

  add(ctx.line([[0, -3.4], [0, 3.4]], { color: 'line', width: 1 }));
  // Titles sit under each half, so the top-right corner stays free for the page's play/pause button.
  const TY = -2.1;
  const heads = [
    [lab('Forward', -4.08, TY, { size: 14, weight: 600, color: 'ink', anchor: 'right' }), lab('where lost?', -3.92, TY, { anchor: 'left' })],
    [lab('Inverse', 3.92, TY, { size: 14, weight: 600, color: 'ink', anchor: 'right' }), lab('is it enough?', 4.08, TY, { anchor: 'left' })],
  ];

  /* ── left: push → state → output ── */
  const Y0 = 0.45, KX = -6.9, UY = -0.85, deg = Math.PI / 180;
  L.push([add(ctx.dot([KX, Y0, 1], { r: 0.45, hollow: true, ring: 0.2, color: 'ink' })), 1]);
  [130, 100, 70, 40].forEach((a) => {
    const c = Math.cos(a * deg), s = Math.sin(a * deg);
    L.push([add(ctx.line([[KX + 0.56 * c, Y0 + 0.56 * s], [KX + 0.68 * c, Y0 + 0.68 * s]], { color: 'faint', width: 1.5 })), 1]);
  });
  const ind = add(ctx.group());
  ind.position.set(KX, Y0, 1.2);
  const indLine = ctx.line([[0, 0], [0, 0.32]], { color: 'accent', width: 2.5 });
  ind.add(indLine);
  L.push([indLine, 1]);
  const ROT0 = 10 * deg, ROT1 = -20 * deg; // pointer on the 100° tick, then one tick on (70°)
  L.push([add(ctx.arrow([-6.3, Y0, 1], [-5.45, Y0, 1], { color: 'muted', width: 1.5, head: 7 })), 1]);

  const COLS = 5, ROWS = 3, GX = -5.15, GD = 0.55, GY = 0.6;
  const nodes = [];
  for (let j = 0; j < ROWS; j++) for (let i = 0; i < COLS; i++) {
    const x0 = GX + i * GD, y0 = Y0 + (j - 1) * GY;
    nodes.push({ i, j, x0, y0, x: x0, y: y0, d: add(ctx.dot([x0, y0, 1.5], { r: 0.075, color: 'faint' })) });
  }
  const at = (i, j) => nodes[j * COLS + i];
  const rowPts = (j) => Array.from({ length: COLS }, (_, i) => [at(i, j).x, at(i, j).y, 1]);
  const colPts = (i) => Array.from({ length: ROWS }, (_, j) => [at(i, j).x, at(i, j).y, 1]);
  const springs = [];
  for (let j = 0; j < ROWS; j++) springs.push({ ln: add(ctx.line(rowPts(j), { color: 'faint', width: 1.5 })), pts: () => rowPts(j) });
  for (let i = 0; i < COLS; i++) springs.push({ ln: add(ctx.line(colPts(i), { color: 'faint', width: 1.5 })), pts: () => colPts(i) });

  const RX = -2.3;
  L.push([add(ctx.arrow([-2.72, Y0, 1], [-1.86, Y0, 1], { color: 'muted', width: 1.5, head: 7 })), 1]);
  const ring = add(ctx.dot([RX, Y0, 2], { r: 0.2, hollow: true, ring: 0.28, color: 'bad' }));
  const ripple = add(ctx.dot([RX, Y0, 2], { r: 0.2, hollow: true, ring: 0.18, color: 'bad' }));
  const GCX = -1.15, GCY = Y0 - 0.3, GR = 0.62, NA = 118 * deg;
  L.push([add(ctx.line(Array.from({ length: 33 }, (_, k) => [GCX + GR * Math.cos((Math.PI * k) / 32), GCY + GR * Math.sin((Math.PI * k) / 32)]), { color: 'muted', width: 2 })), 1]);
  L.push([add(ctx.line([[GCX, GCY, 1], [GCX + 0.5 * Math.cos(NA), GCY + 0.5 * Math.sin(NA), 1]], { color: 'ink', width: 2.5 })), 1]);
  L.push([add(ctx.dot([GCX, GCY, 1.5], { r: 0.06, color: 'ink' })), 1]);
  const tok1 = add(ctx.dot([0, 0, 3], { px: 4, color: 'accent' }));
  const tok2 = add(ctx.dot([0, 0, 3], { px: 4, color: 'accent' }));
  const lost = lab('lost here', RX, Y0 + 1.2, { color: 'bad', weight: 600, anchor: 'bottom' });
  for (const [s, x] of [['push', KX], ['state', GX + 2 * GD], ['output', GCX]]) L.push([lab(s, x, UY, { anchor: 'top' }), 1]);

  /* ── right: state space → outputs ── */
  const CX = 2.75, CY = 1.0, S = 1.15;
  const P = (X, Y, Z) => [CX + S * (X + 0.42 * Z), CY + S * (Y + 0.3 * Z), 0.5];
  const corner = (k) => P(k & 4 ? 1 : -1, k & 2 ? 1 : -1, k & 1 ? 1 : -1);
  for (let a = 0; a < 8; a++) for (let b = a + 1; b < 8; b++) {
    const x = a ^ b;
    if (x & (x - 1)) continue; // corners that differ in one coordinate share an edge
    const axis = a === 0; // the three edges from the front-bottom-left corner read as axes
    R.push([add(ctx.line([corner(a), corner(b)], { color: axis ? 'faint' : 'line', width: axis ? 1.5 : 1.2 })), 1]);
  }
  const Z3 = [0.15, 0, 0.1], zp = P(...Z3);
  const bez = (a, c, b) => Array.from({ length: 41 }, (_, k) => {
    const s = k / 40, r = 1 - s;
    return P(...[0, 1, 2].map((q) => r * r * a[q] + 2 * r * s * c[q] + s * s * b[q]));
  });
  const OX = 6.7;
  const inputs = [
    { col: 't', name: 'x₁', a: [-0.9, 0.9, 0.6], c: [0.4, 0.9, -0.4], out: [OX - 0.2, CY + 1.0], tag: 'A' },
    { col: 'teddy', name: 'x₂', a: [-0.9, -0.9, -0.7], c: [-0.6, -0.1, 0.9], out: [OX - 0.2, CY - 1.05], tag: 'B' },
  ].map((o) => {
    const pts = bez(o.a, o.c, Z3);
    return Object.assign(o, {
      ln: add(ctx.line(pts, { color: o.col, width: 2.5 })),
      head: add(ctx.dot([0, 0, 3], { px: 4.5, color: o.col })),
      start: add(ctx.dot([pts[0][0], pts[0][1], 2], { px: 3.5, color: o.col })),
      nm: lab(o.name, pts[0][0] - 0.14, pts[0][1], { anchor: 'right', color: o.col, weight: 600, size: 13, bg: 'card' }),
      arr: add(ctx.arrow([zp[0], zp[1], 1], [o.out[0], o.out[1], 1], { color: o.col, width: 2, head: 8 })),
    });
  });
  R.push([add(ctx.line([[OX, CY - 1.35], [OX, CY + 1.55]], { color: 'muted', width: 1.5 })), 1]);
  for (const o of inputs) {
    R.push([add(ctx.line([[OX - 0.14, o.out[1]], [OX + 0.14, o.out[1]]], { color: 'ink', width: 2 })), 1]);
    R.push([lab(o.tag, OX + 0.3, o.out[1], { anchor: 'left', color: 'ink', weight: 600, size: 13 }), 1]);
  }
  R.push([lab('state', CX, UY, { anchor: 'top' }), 1]);
  R.push([lab('output', OX, UY, { anchor: 'top' }), 1]);
  const halo = add(ctx.dot([zp[0], zp[1], 3.5], { r: 0.42, color: 'accent', opacity: 0 }));
  const zDot = add(ctx.dot([zp[0], zp[1], 4], { r: 0.16, color: 'accent' }));
  const zLab = lab('z', zp[0] - 0.28, zp[1] + 0.26, { anchor: 'bottom-right', color: 'ink', weight: 600, size: 13 });
  const crack = add(ctx.line([[0.06, 0.38], [-0.07, 0.13], [0.06, -0.03], [-0.08, -0.21], [0.01, -0.4]].map(([dx, dy]) => [zp[0] + dx, zp[1] + dy, 4.5]), { color: 'bad', width: 1.8 }));
  const notLab = lab('not enough', zp[0] + 0.2, zp[1] - 0.62, { anchor: 'top', color: 'bad', weight: 600, bg: 'card' });

  /* sL, sR: progress (0..1) through each half's 6 s; aL, aR: how active each half is */
  function frame(sL, aL, sR, aR) {
    const dL = lerp(0.3, 1, aL), dR = lerp(0.3, 1, aR);
    heads[0].forEach((h) => h.setOpacity(lerp(0.45, 1, aL)));
    heads[1].forEach((h) => h.setOpacity(lerp(0.45, 1, aR)));
    for (const [o, b] of L) o.setOpacity(b * dL);
    for (const [o, b] of R) o.setOpacity(b * dR);

    // Forward: push one tick, the state follows, the readout stays flat.
    ind.rotation.z = lerp(ROT0, ROT1, seg(sL, 0.06, 0.16) - seg(sL, 0.9, 0.99));
    tok1.position.set(lerp(-6.3, -5.45, seg(sL, 0.1, 0.22, ease.inOutSine)), Y0, 3);
    tok1.setOpacity(pulse(sL, 0.1, 0.24, 0.03) * dL);
    const k = seg(sL, 0.18, 0.32) * (1 - seg(sL, 0.78, 0.92));
    const front = seg(sL, 0.18, 0.46, ease.inOutSine) * (COLS + 1.5);
    const ph = 2 * Math.PI * 2 * sL;
    for (const n of nodes) {
      const a = clamp(front - n.i) * k;
      n.x = n.x0 + 0.16 * a * Math.sin(ph - n.i * 0.9 + n.j * 0.6);
      n.y = n.y0 + 0.16 * a * Math.cos(ph * 1.5 - n.i * 0.7 + n.j * 1.3);
      n.d.position.set(n.x, n.y, 1.5);
      n.d.setColor(mix('faint', 'accent', a));
      n.d.setOpacity(dL);
    }
    mix('faint', 'accent', k);
    for (const s of springs) s.ln.setPoints(s.pts()).setColor(tmp).setOpacity(dL);
    tok2.position.set(lerp(-2.72, RX, seg(sL, 0.42, 0.5, ease.outCubic)), Y0, 3);
    tok2.setOpacity(pulse(sL, 0.42, 0.58, 0.03) * dL);
    const win = pulse(sL, 0.48, 0.9, 0.04);
    ring.setOpacity(win * dL);
    const q = ((clamp(sL, 0.48, 0.9) - 0.48) / 0.21) % 1;
    ripple.setRadius(lerp(0.2, 0.4, ease.outSine(q)));
    ripple.setOpacity(win * (1 - q) * 0.8 * dL);
    lost.setOpacity(pulse(sL, 0.52, 0.9, 0.05) * dL);

    // Inverse: two inputs land on one state; that state leads to two outputs.
    const on = pulse(sR, 0.02, 0.94, 0.05), fade = 1 - seg(sR, 0.84, 0.94);
    const p = seg(sR, 0.06, 0.36, ease.inOutSine);
    for (const o of inputs) {
      o.ln.setProgress(p);
      o.ln.setOpacity(fade * dR);
      o.head.position.copy(ctx.along(o.ln, p)).setZ(3);
      o.head.setOpacity(pulse(sR, 0.05, 0.42, 0.04) * dR);
      o.start.setOpacity(on * dR);
      o.nm.setOpacity(on * dR);
      o.arr.setProgress(seg(sR, 0.46, 0.64));
      o.arr.setOpacity(fade * dR);
    }
    zDot.setOpacity(on * lerp(0.35, 1, seg(sR, 0.3, 0.38)) * dR);
    halo.setOpacity(seg(sR, 0.34, 0.44) * fade * (0.2 + 0.08 * ctx.wave(sR, 3)) * dR);
    zLab.setOpacity(on * dR);
    crack.setProgress(seg(sR, 0.66, 0.74, ease.outCubic));
    crack.setOpacity(pulse(sR, 0.65, 0.94, 0.04) * dR);
    notLab.setOpacity(pulse(sR, 0.7, 0.94, 0.05) * dR);
  }

  return {
    scene, camera, period: PERIOD, still: 3.6,
    update(t) {
      // Reduced motion: one frame with both halves at their key moment.
      if (ctx.reducedMotion) return frame(0.62, 1, 0.8, 1);
      const u = ctx.loopT(t, PERIOD);
      const aR = pulse(u, 0.48, 0.985, 0.035);
      frame(clamp(u / 0.48), 1 - aR, u < 0.5 ? 0 : clamp((u - 0.5) / 0.46), aR);
    },
    dispose() { L.length = 0; R.length = 0; nodes.length = 0; springs.length = 0; inputs.length = 0; },
  };
}
