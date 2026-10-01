/*
 * fw-nkt-compression — anm-framework.html, #teddy-nkt ("5.4 · Finding: NK and T nearly merge").
 * The sufficiency test on near-matched cells, drawn with springs.
 * 1 A slowly turning 3D cloud stands for TEDDY's embedding (NK and T dots that touch along a seam;
 *   the layout is schematic, not data). One NK–T pair across the seam is circled: nearly merged.
 * 2 Its CD3 on two rails: measured (the gap is a relaxed spring) and our head's prediction from the
 *   frozen embedding (the same spring, squeezed, still pointing the same way). Bead positions are
 *   illustrative; the squeeze uses the page's CD3 ratio.
 * 3 Six more seam pairs are circled (they stand for the page's 361 pairs) and a rack shows the predicted ÷ measured gap
 *   for the page's four proteins as springs squeezed from 1 (dashed = the whole measured gap).
 *   Then random NK–T pairs far apart: their spring stays about full length, inside the page's band.
 * 4 A fork with two half-open doors, "lost in the embedding" and "unused by our head"; neither is lit,
 *   and a sign says v3 E5-M is running (no result yet).
 * Every number drawn comes from NUMBERS below (section text and its sources). 12 s seamless loop;
 * the still frame shows beat 3 complete. The top-right corner stays empty for the play/pause button.
 */
const NUMBERS = {
  nearPairs: 453,                                                               // NK–T pairs the embedding nearly merges
  kept: [['CD3', 0.27], ['CD56', 0.24], ['CD94', 0.31], ['CD335', 0.88]],       // predicted ÷ measured gap, those pairs
  random: [0.8, 1.5],                                                           // the same ratio, random NK–T pairs
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: 11.5, color: 'ink' }, o)));
  const at = (o, x, y, z = o.position.z) => { o.position.set(x, y, z); return o; };
  const dot = (color, o = {}) => add(ctx.dot([0, 0, 0], Object.assign({ color }, o)));
  const ln = (pts, o) => add(ctx.line(pts, o));

  /* a spring from xa to xb: short straight leads, then 2 × coils zig-zag strokes (constant point count) */
  const spring = (xa, xb, y, z, amp, coils) => {
    const L = Math.max(0.05, xb - xa), lead = Math.min(0.2, L * 0.14), n = 2 * coils, pts = [[xa, y, z], [xa + lead, y, z]];
    for (let k = 1; k < n; k++) pts.push([xa + lead + ((L - 2 * lead) * k) / n, y + (k % 2 ? amp : -amp), z]);
    pts.push([xa + L - lead, y, z], [xa + L, y, z]);
    return pts;
  };
  const springLine = (color, width, coils) => ln(spring(0, 1, 0, 5, 0.1, coils), { color, width });

  /* ── beat 1: the embedding cloud (fixed schematic layout) ── */
  const C = [-4.35, -0.05], SC = 1.7, CP = Math.cos(0.3), SP = Math.sin(0.3);
  let cy = 1, sy = 0;
  const proj = (p, out) => {
    const x = p[0] * cy + p[2] * sy, z = p[2] * cy - p[0] * sy;
    out.x = C[0] + SC * x; out.y = C[1] + SC * (p[1] * CP - z * SP); out.d = p[1] * SP + z * CP;
    return out;
  };
  const rnd = ctx.rand(361);
  const g = () => { let s = 0; for (let i = 0; i < 4; i++) s += rnd(); return clamp((s - 2) / 0.58, -2.2, 2.2); };
  const cloud = [];
  for (const [key, cx, sgn] of [['nk', -0.74, -1], ['t', 0.74, 1]]) {
    for (let i = 0; i < 30; i++) {
      let x = cx + 0.4 * g();
      if (x * sgn < 0.07) x = sgn * (0.14 - x * sgn); // each cluster stays on its side of the seam
      cloud.push({ p: [x, 0.42 * g(), 0.4 * g()], d: dot(key), s: {} });
    }
  }
  /* pairs across the seam; MAIN is the one followed in beat 2 */
  const rs = ctx.rand(985), MAIN = 3;
  const seam = Array.from({ length: 7 }, (_, i) => {
    const y = i === MAIN ? 0.05 : lerp(-0.9, 0.9, i / 6) + 0.06 * (rs() - 0.5);
    const z = i === MAIN ? 0.45 : (i % 2 ? 1 : -1) * (0.3 + 0.3 * rs());
    return {
      a: [-0.085, y, z], b: [0.085, y + 0.05 * (rs() - 0.5), z + 0.08 * (rs() - 0.5)], sa: {}, sb: {},
      back: i === MAIN ? dot('card') : null, da: dot('nk'), db: dot('t'),
      ring: dot('ink', { hollow: true, ring: i === MAIN ? 0.2 : 0.16 }),
    };
  });
  const main = seam[MAIN], mid = { x: 0, y: 0, r: 11 };
  /* random NK–T pairs, far apart */
  const far = [[[-1.25, 0.55, 0.25], [1.2, -0.5, -0.2]], [[-1.15, -0.6, -0.3], [1.3, 0.45, 0.3]], [[-1.45, -0.05, 0.4], [1.05, 0.12, -0.45]]]
    .map(([a, b]) => ({ a, b, sa: {}, sb: {}, da: dot('nk'), db: dot('t'), rod: ln([[0, 0, 2], [1, 0, 2]], { color: 'muted', width: 1.3, dashed: [4, 3] }) }));
  const cTitle = at(lab('TEDDY’s embedding', { color: 'teddy', weight: 600, size: 12.5 }), C[0], 3.85, 10);
  const cSub = at(lab('schematic layout', { color: 'faint', size: 10.5 }), C[0], 3.38, 10);
  const cl = [['NK cells', 'nk', -0.74], ['T cells', 't', 0.74]].map(([s, c, x]) => ({ x, l: lab(s, { color: c, weight: 600, anchor: 'bottom' }) }));
  const merged = lab('nearly merged', { weight: 600, anchor: 'top', bg: 'card', bgOpacity: 0.88, order: 12 });
  /* legend under the cloud */
  const LX = -7.35, LY = [-2.6, -3.3];
  const legRing = at(dot('ink', { hollow: true, ring: 0.16, px: 6.5 }), LX, LY[0], 2);
  const legRod = ln([[LX - 0.3, LY[1], 2], [LX + 0.3, LY[1], 2]], { color: 'muted', width: 1.3, dashed: [4, 3] });
  const leg1 = at(lab(`${NUMBERS.nearPairs} nearly merged pairs`, { anchor: 'left', weight: 600 }), LX + 0.45, LY[0], 10);
  const leg2 = at(lab('random NK–T pairs', { anchor: 'left', color: 'muted', weight: 600 }), LX + 0.45, LY[1], 10);

  /* ── beat 2: CD3 of the main pair, measured vs our head's prediction ── */
  const B = { x0: 0.6, x1: 7.2, ym: 1.45, yp: -0.95, xa: 1.25, xb: 6.55 };
  const bc = (B.xa + B.xb) / 2, bg = (B.xb - B.xa) * NUMBERS.kept[0][1], pa = bc - bg / 2, pb = bc + bg / 2;
  const bTitle = at(lab('CD3 of this pair', { anchor: 'left', weight: 600, size: 12.5 }), 0, 3.8, 10);
  const rails = [B.ym, B.yp].map((y) => [
    ln([[B.x0, y, 2], [B.x1, y, 2]], { color: 'line', width: 1.6 }),
    ln([[B.x0, y - 0.13, 2], [B.x0, y + 0.13, 2]], { color: 'line', width: 1.4 }),
    ln([[B.x1, y - 0.13, 2], [B.x1, y + 0.13, 2]], { color: 'line', width: 1.4 }),
  ]).flat();
  const mLab = at(lab('measured', { anchor: 'bottom-left', weight: 600 }), B.x0, B.ym + 0.42, 10);
  const pLab = at(lab('our head’s prediction', { anchor: 'bottom-left', weight: 600, color: 'teddy' }), B.x0, B.yp + 0.42, 10);
  const note = lab('right direction, squeezed gap', { anchor: 'top', color: 'muted', size: 11 });
  const mSpring = springLine('measured', 1.6, 7), pSpring = springLine('teddy', 1.8, 7);
  const mBeads = [dot('nk'), dot('t')], pBeads = [dot('nk'), dot('t')];

  /* ── beat 3: the rack, predicted ÷ measured gap per protein, then random pairs ── */
  const R = { x0: 1.9, s: 3.3, rows: [1.75, 0.95, 0.15, -0.65], yr: -1.95 };
  const XR = (v) => R.x0 + R.s * v;
  const [LO, HI] = NUMBERS.random;
  const rTitle = at(lab('predicted gap', { anchor: 'left', weight: 600, size: 12.5 }), 0, 3.8, 10);
  const one = ln([[XR(1), 2.35, 3], [XR(1), R.yr - 0.31, 3]], { color: 'muted', width: 1.3, dashed: [4, 4] });
  // right-anchored just past the dashed line: centred on it, the label ran under the corner Play button on phones
  const oneLab = at(lab('1 = measured gap', { anchor: 'bottom-right', color: 'muted' }), XR(1) + 0.25, 2.42, 10);
  const rows = NUMBERS.kept.map(([name, v], i) => ({
    v, y: R.rows[i], spring: springLine('teddy', 1.6, 5), a: dot('nk'), b: dot('t'),
    name: at(lab(name, { anchor: 'right', weight: 600 }), R.x0 - 0.3, R.rows[i], 10),
    val: at(lab(v.toFixed(2), { anchor: 'left', color: 'teddy', weight: 600, mono: true }), XR(1) + 0.3, R.rows[i], 10),
  }));
  const band = at(add(ctx.box(XR(HI) - XR(LO), 0.62, { color: 'faint', radius: 0.08 })), (XR(LO) + XR(HI)) / 2, R.yr, 0.5);
  const bandLab = at(lab(`${LO}–${HI}`, { anchor: 'top', color: 'muted', weight: 600, mono: true }), (XR(LO) + XR(HI)) / 2, R.yr - 0.42, 10);
  const rName = at(lab('random pairs', { anchor: 'right', weight: 600, color: 'muted' }), R.x0 - 0.3, R.yr, 10);
  const rSpring = springLine('teddy', 1.6, 5), rA = dot('nk'), rB = dot('t');

  /* ── beat 4: the fork, two half-open doors and the Mode A sign ── */
  const F = { jx: 1.15, jy: -0.25, dx: 4.55, w: 1.0, h: 1.45 };
  const path = ln([[0, 0, 2], [1, 0, 2]], { color: 'muted', width: 1.4, dashed: [4, 4] });
  const DOORS = [{ y: 1.55, tok: 'teddy', text: 'lost in the\nembedding', col: 'teddy', bend: -0.12 },
    { y: -2.05, tok: 'muted', text: 'unused by\nour head', col: 'ink', bend: 0.12 }];
  const quad = (tok) => {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(12), 3));
    geo.setIndex([0, 1, 2, 0, 2, 3]);
    const mat = new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, side: THREE.DoubleSide });
    ctx.bind(mat, tok);
    const mesh = add(new THREE.Mesh(geo, mat));
    mesh.frustumCulled = false;
    mesh.setPts = (pts) => { const a = geo.attributes.position; pts.forEach((p, i) => a.setXYZ(i, p[0], p[1], p[2])); a.needsUpdate = true; };
    return mesh;
  };
  const doors = DOORS.map((D) => ({
    D,
    branch: add(ctx.arrow([F.jx, F.jy, 2], [F.dx - F.w / 2 - 0.12, D.y, 2], { color: 'muted', width: 1.4, head: 7, bend: D.bend, dashed: [4, 4] })),
    hole: at(add(ctx.box(F.w, F.h, { color: 'faint', radius: 0.03 })), F.dx, D.y, 1),
    frame: at(add(ctx.box(F.w, F.h, { color: null, stroke: D.tok, strokeWidth: 1.6, radius: 0.03 })), F.dx, D.y, 1.2),
    leaf: quad('card'),
    edge: ln([[0, 0, 1.6], [1, 0, 1.6], [1, 1, 1.6], [0, 1, 1.6]], { color: D.tok, width: 1.6, closed: true }),
    knob: dot(D.tok, { px: 2.4 }),
    l: at(lab(D.text, { anchor: 'left', weight: 600, color: D.col, align: 'left' }), F.dx + F.w / 2 + 0.3, D.y, 10),
  }));
  const sign = at(lab('v3 E5-M\nrunning', { color: 'accent', weight: 600, bg: 'accent', bgOpacity: 0.13, pad: 6 }), F.dx, F.jy, 10);

  function update(t) {
    const u = ctx.loopT(t, PERIOD), ppu = ctx.ppu(), px = (p) => p / ppu;
    const th = 0.42 * Math.sin(2 * Math.PI * u);
    cy = Math.cos(th); sy = Math.sin(th);
    const end = 1 - seg(u, 0.95, 0.99, ease.inOutSine);
    const hl = seg(u, 0.03, 0.1) * end;                                   // main pair circled
    const p1 = seg(u, 0.16, 0.2) * (1 - seg(u, 0.45, 0.49));             // beat 2 panel
    const rk = seg(u, 0.49, 0.53) * (1 - seg(u, 0.75, 0.79));            // beat 3 rack
    const rn = seg(u, 0.64, 0.68) * (1 - seg(u, 0.75, 0.79));            // random pairs
    const fk = seg(u, 0.78, 0.82) * end;                                 // beat 4 fork
    const dim = 1 - 0.22 * Math.max(p1, rk, fk);
    const dp = clamp(0.075 * ppu, 2.2, 3.6);
    const shade = (d) => dim * clamp(0.55 + 0.25 * d, 0.25, 0.85);

    /* beat 1: cloud, seam pairs, random pairs */
    for (const q of cloud) { proj(q.p, q.s); at(q.d, q.s.x, q.s.y, 1 + 0.3 * q.s.d).setPx(dp).setOpacity(shade(q.s.d)); }
    cl.forEach((c) => at(c.l, C[0] + SC * c.x * cy, C[1] + SC * 1.3, 10));
    seam.forEach((s, i) => {
      proj(s.a, s.sa); proj(s.b, s.sb);
      const isMain = i === MAIN, k = isMain ? 0 : i < MAIN ? i : i - 1;
      const lit = isMain ? hl : seg(u, 0.47 + 0.012 * k, 0.5 + 0.012 * k) * (1 - seg(u, 0.75, 0.79)) * (1 - 0.55 * rn);
      const r = dp + (isMain ? 1.6 : 0.6) * lit;
      at(s.da, s.sa.x, s.sa.y, 4).setPx(r).setOpacity(Math.max(shade(s.sa.d), lit));
      at(s.db, s.sb.x, s.sb.y, 4).setPx(r).setOpacity(Math.max(shade(s.sb.d), lit));
      const mx = (s.sa.x + s.sb.x) / 2, my = (s.sa.y + s.sb.y) / 2;
      const half = (Math.hypot(s.sb.x - s.sa.x, s.sb.y - s.sa.y) * ppu) / 2 + r;   // the ring always encloses both dots
      const rp = Math.max(isMain ? 11 : 6.5, half + (isMain ? 4 : 2));
      at(s.ring, mx, my, 4.2).setPx(rp).setOpacity(lit * (isMain ? 1 : 0.8));
      if (isMain) { mid.x = mx; mid.y = my; mid.r = rp; at(s.back, mx, my, 3.9).setPx(rp - 1).setOpacity(0.7 * hl); }
    });
    at(merged, mid.x, mid.y - px(mid.r + 4), 10).setOpacity(seg(u, 0.09, 0.13) * (1 - seg(u, 0.45, 0.49)));
    far.forEach((f, i) => {
      proj(f.a, f.sa); proj(f.b, f.sb);
      at(f.da, f.sa.x, f.sa.y, 3).setPx(dp + 0.6 * rn).setOpacity(Math.max(shade(f.sa.d), rn));
      at(f.db, f.sb.x, f.sb.y, 3).setPx(dp + 0.6 * rn).setOpacity(Math.max(shade(f.sb.d), rn));
      f.rod.setPoints([[f.sa.x, f.sa.y, 2], [f.sb.x, f.sb.y, 2]]).setProgress(seg(u, 0.64 + 0.015 * i, 0.68 + 0.015 * i)).setOpacity(0.9 * rn);
    });
    const lk = seg(u, 0.5, 0.53) * (1 - seg(u, 0.75, 0.79));
    legRing.setOpacity(0.8 * lk); leg1.setOpacity(lk);
    legRod.setOpacity(0.9 * rn); leg2.setOpacity(rn);

    /* beat 2: the main pair's CD3, measured (relaxed spring) and predicted (the same spring, squeezed) */
    const bp = clamp(0.17 * ppu, 4, 6);
    bTitle.setOpacity(p1);
    rails.forEach((r) => r.setOpacity(0.9 * p1));
    const fly = seg(u, 0.17, 0.25), mo = p1 * seg(u, 0.17, 0.185);
    [[main.sa, B.xa], [main.sb, B.xb]].forEach(([s, x], j) => at(mBeads[j], lerp(s.x, x, fly), lerp(s.y, B.ym, fly), 6).setPx(bp).setOpacity(mo));
    mSpring.setPoints(spring(B.xa, B.xb, B.ym, 5, px(5), 7)).setOpacity(p1 * seg(u, 0.24, 0.27));
    mLab.setOpacity(p1 * seg(u, 0.2, 0.24));
    const drop = seg(u, 0.28, 0.33), sq = seg(u, 0.34, 0.41), po = p1 * seg(u, 0.28, 0.29);
    const qa = lerp(B.xa, pa, sq), qb = lerp(B.xb, pb, sq), qy = lerp(B.ym, B.yp, drop);
    at(pBeads[0], qa, qy, 6).setPx(bp).setOpacity(po);
    at(pBeads[1], qb, qy, 6).setPx(bp).setOpacity(po);
    pSpring.setPoints(spring(qa, qb, qy, 5.1, px(5), 7)).setOpacity(po);
    pLab.setOpacity(p1 * seg(u, 0.3, 0.33));
    at(note, bc, B.yp - px(12), 10).setOpacity(p1 * seg(u, 0.4, 0.43));

    /* beat 3: springs squeezed from 1 (the measured gap) to the page's ratios; random pairs stay long */
    const sp = clamp(0.13 * ppu, 3.2, 4.8), amp = px(4);
    rTitle.setOpacity(rk); one.setOpacity(0.9 * rk); oneLab.setOpacity(rk);
    rows.forEach((r, i) => {
      const a0 = 0.51 + 0.02 * i, on = seg(u, a0, a0 + 0.02) * rk, v = lerp(1, r.v, seg(u, a0 + 0.03, a0 + 0.08));
      r.spring.setPoints(spring(XR(0), XR(v), r.y, 5, amp, 5)).setOpacity(on);
      at(r.a, XR(0), r.y, 6).setPx(sp).setOpacity(on);
      at(r.b, XR(v), r.y, 6).setPx(sp).setOpacity(on);
      r.name.setOpacity(rk);
      r.val.setOpacity(seg(u, a0 + 0.07, a0 + 0.1) * rk);
    });
    const rr = lerp(LO, HI, 0.5 + 0.38 * Math.sin(2 * Math.PI * 5 * (u - 0.72))); // breathes inside the band; mid-band in the still
    band.setOpacity(0.16 * rn); bandLab.setOpacity(rn); rName.setOpacity(rn);
    rSpring.setPoints(spring(XR(0), XR(rr), R.yr, 5, amp, 5)).setOpacity(rn);
    at(rA, XR(0), R.yr, 6).setPx(sp).setOpacity(rn);
    at(rB, XR(rr), R.yr, 6).setPx(sp).setOpacity(rn);

    /* beat 4: where was the gap lost? two half-open doors, neither lit */
    path.setPoints([[mid.x + px(mid.r + 2), mid.y, 2], [F.jx, F.jy, 2]]).setProgress(seg(u, 0.78, 0.83)).setOpacity(0.9 * fk);
    const dk = seg(u, 0.84, 0.88) * fk, ang = 1.0 * seg(u, 0.86, 0.92, ease.outCubic);
    doors.forEach((d) => {
      const y = d.D.y, xh = F.dx + F.w / 2, xf = xh - F.w * Math.cos(ang), k = 0.11 * Math.sin(ang), h2 = F.h / 2;
      const pts = [[xh, y - h2, 1.5], [xh, y + h2, 1.5], [xf, y + h2 + k, 1.5], [xf, y - h2 - k, 1.5]];
      d.branch.setProgress(seg(u, 0.82, 0.87)).setOpacity(fk);
      d.hole.setOpacity(0.3 * dk); d.frame.setOpacity(dk);
      d.leaf.setPts(pts); d.leaf.material.opacity = dk;
      d.edge.setPoints(pts).setOpacity(dk);
      at(d.knob, xf + 0.13, y, 1.7).setOpacity(dk);
      d.l.setOpacity(seg(u, 0.88, 0.91) * fk);
    });
    sign.setOpacity(seg(u, 0.9, 0.93) * fk);
    cTitle.setOpacity(1); cSub.setOpacity(1);
  }

  return {
    scene, camera, period: PERIOD, still: 0.72 * PERIOD,
    update,
    dispose() { cloud.length = 0; seam.length = 0; far.length = 0; rows.length = 0; doors.length = 0; },
  };
}
