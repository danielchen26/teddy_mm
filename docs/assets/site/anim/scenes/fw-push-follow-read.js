/*
 * fw-push-follow-read — anm-framework.html #idea, "Push, follow, read, compare".
 * Three stations on a gentle arc. Source: a knob on a tick rail slides by exactly one tick (δu).
 * Retained state: a 7 × 7 lattice of spheres on springs; a ripple enters from the left and moves
 * the spheres (δs), weakening and bending as it spreads. Readout: a semicircular gauge whose needle
 * moves by δO. Before any push, two bands are drawn on the gauge as dashed outlines (the boundary
 * declared in advance): the wide outer band 'still works' and the narrow inner band 'exactly right'.
 * When the needle settles, the band under it fills. Then the knob pushes the other way and it all
 * mirrors. 12 s seamless loop. Geometry is illustrative only: no values are shown or implied.
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9, center: [0, 0] });
  const PERIOD = 12, N = 7, TICK = 0.65, DEG = Math.PI / 180, TH = 31;
  const add = (o) => { scene.add(o); return o; };
  const text = (s, x, y, o) => { const l = ctx.label(s, o); l.position.set(x, y, 5); return add(l); };
  const MUTED = { size: 12, color: 'muted', anchor: 'top' };
  const DELTA = { size: 14, weight: 'italic 600', anchor: 'bottom' };
  const bez = (a, b, c, d, n = 28) => Array.from({ length: n + 1 }, (_, k) => {
    const t = k / n, s = 1 - t, w = [s * s * s, 3 * s * s * t, 3 * s * t * t, t * t * t];
    return [w[0] * a[0] + w[1] * b[0] + w[2] * c[0] + w[3] * d[0], w[0] * a[1] + w[1] * b[1] + w[2] * c[1] + w[3] * d[1], 0.05];
  });

  /* ── source: a knob on a rail with a tick scale ── */
  const KX = -5.9, RY = -2.1, RAIL_END = [KX + 2 * TICK + 0.3, RY];
  add(ctx.line([[KX - 2 * TICK - 0.3, RY, 0.1], [RAIL_END[0], RY, 0.1]], { color: 'muted', width: 2 }));
  for (let k = -2; k <= 2; k++) {
    add(ctx.line([[KX + k * TICK, RY - 0.14, 0.1], [KX + k * TICK, RY - (k ? 0.3 : 0.42), 0.1]], { color: k ? 'faint' : 'muted', width: 1.5 }));
  }
  const knob = add(ctx.dot([KX, RY, 0.6], { r: 0.17, color: 'ink' }));
  const AY = RY + 0.5;
  const pushArrow = add(ctx.arrow([KX, AY, 0.5], [KX + TICK, AY, 0.5], { color: 'accent', width: 2.2, head: 7 }));
  const du = text('δu', KX + TICK / 2, AY + 0.16, Object.assign({ color: 'accent' }, DELTA));
  text('source', KX, RY - 0.6, MUTED);

  /* ── retained state: 7 × 7 spheres joined by springs, drawn as an oblique plane ── */
  const L0 = [-3.95, -0.55], CI = 0.42, CJ = [0.22, 0.34], A0 = 0.55, LAM = 4.2, WIDTH = 1.15;
  const rest = (i, j) => [L0[0] + i * CI + j * CJ[0], L0[1] + j * CJ[1]];
  const ENTRY = rest(0, 3), EXIT = rest(N - 1, 3);
  const nodes = [];
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const p = rest(i, j);
      nodes.push({ p, d: Math.hypot(i, j - 3), h: 0, dot: add(ctx.dot([p[0], p[1], 0.4], { r: 0.08, color: ctx.color('muted') })) });
    }
  }
  const springs = [];
  for (let k = 0; k < N; k++) {
    springs.push(Array.from({ length: N }, (_, i) => k * N + i), Array.from({ length: N }, (_, j) => j * N + k));
  }
  const springPts = (idx) => idx.map((n) => [nodes[n].p[0], nodes[n].p[1] + nodes[n].h, 0.2]);
  const springLines = springs.map((idx) => add(ctx.line(springPts(idx), { color: 'faint', width: 1.2, opacity: 0.75 })));
  const ds = text('δs', -1.95, 2.0, Object.assign({ color: 'accent' }, DELTA));
  text('retained state', -2.2, -1.25, MUTED);

  /* ── readout: semicircular gauge ── */
  const G = { x: 2.45, y: -2.0, R: 2.05 };
  const arc = (r, a0, a1, n = 40, z = 0.3) => Array.from({ length: n + 1 }, (_, k) => {
    const a = (a0 + ((a1 - a0) * k) / n) * DEG;
    return [G.x + r * Math.cos(a), G.y + r * Math.sin(a), z];
  });
  add(ctx.line(arc(G.R, 0, 180, 64), { color: 'muted', width: 1.8, closed: true }));
  // the declared boundary: [inner r, outer r, half-angle, corner the label's leader starts from]
  const BANDS = [
    { r0: 1.6, r1: 1.95, half: 43, name: 'still works', corner: [1.95, 47] },
    { r0: 0.62, r1: 0.95, half: 20, name: 'exactly right', corner: [0.62, 70] },
  ].map((b) => {
    const a0 = 90 - b.half, a1 = 90 + b.half;
    b.outline = add(ctx.line([...arc(b.r1, a0, a1, 32, 0.2), ...arc(b.r0, a1, a0, 32, 0.2)], { color: 'accent', width: 1.5, dashed: [5, 4], closed: true }));
    const mat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, opacity: 0, depthWrite: false })), 'accent');
    b.fill = add(new THREE.Mesh(ctx.track(new THREE.RingGeometry(b.r0, b.r1, 32, 1, a0 * DEG, 2 * b.half * DEG)), mat));
    b.fill.position.set(G.x, G.y, 0.05);
    const cx = G.x + b.corner[0] * Math.cos(b.corner[1] * DEG), cy = G.y + b.corner[0] * Math.sin(b.corner[1] * DEG);
    b.lead = add(ctx.line([[cx, cy, 0.15], [G.x + G.R + 0.06, cy, 0.15]], { color: 'faint', width: 1 }));
    b.lab = text(b.name, G.x + G.R + 0.12, cy, { size: 12, color: 'ink', anchor: 'left' });
    return b;
  });
  const ghost = add(ctx.line([[G.x, G.y, 0.35], [G.x, G.y + 1.85, 0.35]], { color: 'faint', width: 1.3, dashed: [3, 3], opacity: 0 }));
  const sweep = add(ctx.line(arc(0.42, 90, 90, 16, 0.36), { color: 'ink', width: 1.5, opacity: 0 }));
  const needle = add(ctx.group());
  needle.position.set(G.x, G.y, 0.6);
  needle.add(ctx.line([[0, 0, 0], [0, 1.85, 0]], { color: 'ink', width: 2.6 }));
  add(ctx.dot([G.x, G.y, 0.65], { r: 0.1, color: 'ink' }));
  const dO = text('δO', G.x, G.y + 2.3, Object.assign({ color: 'ink' }, DELTA));
  text('readout', G.x, G.y - 0.28, MUTED);

  /* ── traced path: rail → lattice (entry, left middle) and lattice (exit, right middle) → gauge ── */
  const RIM_L = [G.x - G.R, G.y];
  const P1 = bez(RAIL_END, [RAIL_END[0], RY + 1.2], [ENTRY[0] - 0.66, ENTRY[1]], ENTRY);
  const P2 = bez(EXIT, [EXIT[0] + 0.83, EXIT[1]], [RIM_L[0] - 0.5, RIM_L[1]], RIM_L);
  add(ctx.line(P1, { color: 'faint', width: 1.3, dashed: [3, 4] }));
  add(ctx.line(P2, { color: 'faint', width: 1.3, dashed: [3, 4] }));
  const route1 = [1, -1].map((s) => [[KX + s * TICK, RY, 0], [RAIL_END[0], RY, 0], ...P1]);
  const route2 = [...P2, [G.x, G.y, 0]];
  const halo = add(ctx.dot([0, 0, 0.85], { px: 11, color: 'accent', opacity: 0 }));
  const core = add(ctx.dot([0, 0, 0.9], { px: 4.5, color: 'accent', opacity: 0 }));

  /* ── timeline: bands drawn first, then push +, reset, push −, reset, bands fade ── */
  const PH = [{ s: 1, a: 0.06 }, { s: -1, a: 0.5 }];
  let arrowSign = 1, flat = true;

  function setLattice(s, f, on) {
    const { muted, accent } = ctx.colors;
    for (const n of nodes) {
      n.h = on ? s * A0 * Math.exp(-n.d / LAM) * Math.exp(-(((n.d - f) / WIDTH) ** 2)) : 0;
      n.dot.position.y = n.p[1] + n.h;
      n.dot.material.color.copy(muted).lerp(accent, clamp(Math.abs(n.h) / 0.16));
    }
    if (on || !flat) springLines.forEach((ln, k) => ln.setPoints(springPts(springs[k])));
    flat = !on;
  }

  return {
    scene, camera, period: PERIOD, still: 5.0,
    resize() {
      const px = clamp(0.08 * ctx.ppu(), 2.3, 5.2);
      nodes.forEach((n) => n.dot.setPx(px));
    },
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const draw = seg(u, 0, 0.05), bandOp = draw * (1 - seg(u, 0.95, 0.99, ease.inOutSine));
      const ph = u < 0.5 ? PH[0] : PH[1], s = ph.s, r = u - ph.a;
      const off = 1 - seg(r, 0.39, 0.44);

      // push: the knob slides one tick; a fixed-length arrow shows δu
      const push = seg(r, 0, 0.04);
      knob.position.x = KX + s * TICK * push * off;
      if (arrowSign !== s) {
        arrowSign = s;
        pushArrow.set([KX, AY, 0.5], [KX + s * TICK, AY, 0.5]);
        du.position.x = KX + (s * TICK) / 2;
      }
      pushArrow.setProgress(push);
      pushArrow.setOpacity(off);
      du.setOpacity(push * off);

      // follow: a glowing pulse leaves the knob, a ripple crosses the lattice, a weaker pulse goes on
      const second = r >= 0.14;
      const pos = ctx.along(second ? route2 : route1[s > 0 ? 0 : 1], second ? seg(r, 0.215, 0.265, ease.inOutSine) : seg(r, 0.04, 0.095, ease.inOutSine));
      const op = second ? pulse(r, 0.215, 0.265, 0.012) : pulse(r, 0.04, 0.095, 0.012);
      core.position.set(pos.x, pos.y, 0.9);
      halo.position.set(pos.x, pos.y, 0.85);
      core.setPx(second ? 3.6 : 4.5).setOpacity(op);
      halo.setPx(second ? 8 : 11).setOpacity(op * 0.25);
      setLattice(s, lerp(-2.8, 9.8, seg(r, 0.085, 0.225, ease.inOutSine)), r > 0.085 && r < 0.225);
      ds.setOpacity(seg(r, 0.1, 0.13) * off);

      // read: the needle moves by δO from its rest position
      const k = seg(r, 0.26, 0.31, ease.outBack) * off;
      const shown = seg(r, 0.26, 0.285) * off;
      needle.rotation.z = -s * TH * k * DEG;
      ghost.setOpacity(shown * 0.9);
      sweep.setPoints(arc(0.42, 90, 90 - s * TH * k, 16, 0.36));
      sweep.setOpacity(shown);
      const a = (90 - s * TH) * DEG;
      dO.position.set(G.x + 2.2 * Math.cos(a), G.y + 2.2 * Math.sin(a), 5);
      if (dO.userData.s !== s) { dO.userData.s = s; dO.setAnchor(s > 0 ? 'bottom-left' : 'bottom-right'); }
      dO.setOpacity(shown);

      // compare: the band under the settled needle fills (outer 'still works'; the inner one stays open)
      const glow = seg(r, 0.305, 0.345) * off;
      BANDS.forEach((b, i) => {
        b.outline.setProgress(draw);
        b.outline.setOpacity(bandOp);
        b.lead.setOpacity(bandOp * 0.9);
        b.lab.setOpacity(bandOp);
        b.fill.material.opacity = i === 0 ? 0.3 * glow : 0;
      });
    },
    dispose() {
      nodes.length = 0;
      springs.length = 0;
    },
  };
}
