/*
 * idx-brings-sufficiency — index.html, #brings lead card (the sufficiency test).
 * One idea: two cells that are neighbours in TEDDY's embedding get close predicted proteins
 * (TEDDY's frozen embedding + our head), although their measured proteins are far apart.
 * 1 A slowly turning 3D cloud stands for TEDDY's embedding (NK and T dots; the layout is illustrative,
 *   not data). One NK cell and one T cell drift until they almost touch: cosine 0.976 (page value).
 * 2 The same two cells on two CD3 rails: predicted (TEDDY + head) 0.07 / 0.06, measured 0.23 / 0.20.
 * 3 Share of the measured gap the predictions keep (page values): CD3 0.27, CD56 0.24, CD94 0.31,
 *   CD335 0.88; dashed 1 = the whole measured gap; band 0.8–1.5 = NK–T pairs that are not neighbours
 *   (faint pairs that keep their gap).
 * 4 An open fork between the embedding and our head: two unlit lamps (the v2 readings, as questions) and the
 *   Mode A probe, marked with the v3 outcome (E3 rejected the repair reading).
 * Wide layout 21 × 8 when the stage is at least 2:1 (desktop 21:8), else narrow 16 × 10 (phone 16:10).
 * 12 s seamless loop; the still frame shows beats 1–3 complete.
 */

/* Every pipeline number this scene prints (all stated on index.html #brings). Refresh after the
 * official-preprocessing rerun. */
const NUMBERS = {
  cosine: 0.970,                  // the NK–T neighbour pair in TEDDY's embedding
  cd3Predicted: [0.10, 0.28],     // CD3 of the pair, predicted by our head from TEDDY's embedding: NK, T
  cd3Measured: [0.30, 1.00],      // CD3 of the pair, measured: NK, T
  gapKept: [['CD3', 0.27], ['CD56', 0.24], ['CD94', 0.31], ['CD335', 0.88]],   // predicted ÷ measured gap, near pairs
  notNeighbours: [0.8, 1.5],      // the same ratio for NK–T pairs that are not neighbours
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 8 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const ln = (o) => add(ctx.line([[0, 0], [1, 0]], o));
  const dot = (color, o = {}) => add(ctx.dot([0, 0, 0], Object.assign({ color }, o)));

  const LAYOUTS = {
    wide: { W: 21, H: 8, C: [-7.65, -0.55], S: 1.65, top: 3.1, cosRight: true, rt: 'CD3 (T-cell marker)', bt: 'share of gap kept',
      r: { x0: -3.7, x1: 2.0, yp: 1.7, ym: -0.5 },
      b: { x0: 4.35, s: 3.85, tx: 3.3, ty: 3.1, rows: [1.7, 0.6, -0.5, -1.6] },
      f: { hx: 2.2, l1: -0.9, out: true } },
    narrow: { W: 16, H: 10, C: [-4.65, -0.6], S: 1.7, top: 4.45, cosRight: false, rt: 'CD3', bt: 'gap kept',
      r: { x0: -0.6, x1: 7.0, yp: 2.8, ym: 0.85 },
      b: { x0: 0.8, s: 4.45, tx: -0.8, ty: -0.55, rows: [-1.42, -2.14, -2.86, -3.58] },
      f: { hx: 4.6, l1: -0.2, out: false } },
  };

  /* ── beat 1: the embedding cloud (fixed illustrative layout) and the neighbour pair ── */
  const rnd = ctx.rand(23);
  const g = () => { let s = 0; for (let i = 0; i < 4; i++) s += rnd(); return clamp((s - 2) / 0.58, -2.2, 2.2); };
  const pts = [];
  for (const [key, cx] of [['nk', -0.8], ['t', 0.8]]) {
    for (let i = 0; i < 32; i++) pts.push({ p: [cx + 0.42 * g(), 0.36 * g(), 0.42 * g()], d: dot(key) });
  }
  const pair = [['nk', [-1.2, 0.35, 0.3], 1], ['t', [1.15, -0.3, -0.25], -1]].map(([key, home, sgn]) => ({
    home, sgn, w: new THREE.Vector3(), back: dot('card'), d: dot(key),
  }));
  const cl = [['NK cells', 'nk', -0.8], ['T cells', 't', 0.8]].map(([t, c, x]) => ({ x, l: lab(t, { color: c, weight: 600, anchor: 'bottom' }) }));
  const cTitle = lab('TEDDY’s embedding', { color: 'teddy', weight: 600, size: 13 });
  const arc = ln({ color: 'ink', width: 1.3 });
  const arcPts = Array.from({ length: 13 }, () => new THREE.Vector3());
  const cosLab = lab(`cosine ${NUMBERS.cosine}`, { weight: 600, bg: 'card', bgOpacity: 0.85, order: 12 });

  /* ── beat 2: CD3 of the same two cells, predicted vs measured ── */
  const V = [NUMBERS.cd3Predicted, NUMBERS.cd3Measured];
  const rTitle = lab('CD3', { anchor: 'left', weight: 600, size: 13 });
  const rails = ['teddy', 'measured'].map((col, i) => ({
    line: ln({ color: col, width: 2 }), gap: ln({ color: col, width: 6 }),
    caps: [ln({ color: col, width: 1.4 }), ln({ color: col, width: 1.4 })],
    name: lab(i ? 'measured' : 'predicted (TEDDY + head)', { anchor: 'bottom-left', weight: 600, color: i ? 'ink' : 'teddy' }),
    beads: [dot('nk'), dot('t')],
    vals: V[i].map((v, j) => lab((j ? 'T ' : 'NK ') + v.toFixed(2), { weight: 600, color: i ? 'ink' : 'teddy', anchor: i ? 'top' : j ? 'top-left' : 'top-right' })),
  }));

  /* ── beat 3: share of the measured gap kept, 4 proteins, vs pairs that are not neighbours ── */
  const bTitle = lab('gap kept', { anchor: 'left', weight: 600, size: 13 });
  const bars = NUMBERS.gapKept.map(([n, v]) => ({
    v, bar: ln({ color: 'teddy' }), name: lab(n, { anchor: 'right' }),
    val: lab(v.toFixed(2), { anchor: 'left', color: 'teddy', weight: 600 }),
  }));
  const base = ln({ color: 'line', width: 1 });
  const one = ln({ color: 'muted', width: 1.4, dashed: [4, 4] });
  const oneLab = lab('1 = whole gap', { anchor: 'bottom', color: 'muted' });
  const band = add(ctx.box(1, 1, { color: 'faint', radius: 0.08 }));
  const [NB0, NB1] = NUMBERS.notNeighbours;
  const bandLab = lab(`not neighbours: ${NB0}–${NB1}`, { anchor: 'top-right', color: 'muted' });
  const ghosts = [0.26, 0.71, 0.46].map((f) => NB0 + f * (NB1 - NB0)).map((c) => ({ c, y: 0, l: ln({ color: 'faint', width: 1 }), a: dot('nk'), b: dot('t') }));

  /* ── beat 4: the open fork, embedding or head, and the Mode A probe ── */
  const fArrow = add(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.6, head: 8 }));
  const oArrow = add(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.6, head: 8 }));
  const outLab = lab('proteins', { anchor: 'left', color: 'muted' });
  const head = add(ctx.box(2.5, 1, { color: 'soft', stroke: 'muted', radius: 0.2 }));
  const headLab = lab('our head', { weight: 600, size: 13 });
  const lamps = ['our head’s readout?', 'mean pooling? RNA?'].map((t) => ({
    br: ln({ color: 'faint', width: 1.3, dashed: [3, 4] }),
    ring: dot('muted', { hollow: true, ring: 0.3 }), bulb: dot('faint'), l: lab(t, { anchor: 'bottom' }),
  }));
  const probe = { line: ln({ color: 'accent', width: 1.6, dashed: [4, 4] }), ring: dot('accent', { hollow: true, ring: 0.32 }),
    core: dot('accent'), l: lab('v2 probes · v3 E3: no repair', { anchor: 'top', color: 'accent', weight: 600 }) };

  /* ── projection of the turning cloud (yaw oscillates, fixed pitch) ── */
  const CP = Math.cos(0.32), SP = Math.sin(0.32);
  let ct = 1, st = 0, L = LAYOUTS.wide;
  const P = { ppu: 50, pp: 6, mh: 0.1 };
  const proj = (p, out) => {
    const x = p[0] * ct + p[2] * st, z = p[2] * ct - p[0] * st;
    return out.set(L.C[0] + L.S * x, L.C[1] + L.S * (p[1] * CP - z * SP), p[1] * SP + z * CP);
  };
  const RX = (v) => lerp(L.r.x0, L.r.x1, v), XB = (v) => L.b.x0 + L.b.s * v;

  function layout() {
    const w = ctx.width || 1100, h = ctx.height || 420, a = w / h;
    L = LAYOUTS[a >= 2 ? 'wide' : 'narrow'];
    camera.userData.animFit = { width: L.W, height: L.H, fit: 'contain' };
    const hh = a > L.W / L.H ? L.H / 2 : L.W / 2 / a;
    camera.left = -hh * a; camera.right = hh * a; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    const ppu = (P.ppu = h / (2 * hh)), px = (k, lo, hi) => clamp(k * ppu, lo, hi);
    P.pp = px(0.15, 3.6, 7);
    P.mh = (2 * P.pp + 4) / ppu / (2 * L.S * CP); // the pair ends ~4 px apart at any size
    pts.forEach((q) => q.d.setPx(px(0.075, 2.2, 4)));
    pair.forEach((c) => { c.d.setPx(P.pp); c.back.setPx(P.pp + 1.5); });
    cTitle.position.set(L.C[0], L.top, 0);
    cosLab.setAnchor(L.cosRight ? 'left' : 'top');
    rTitle.setText(L.rt).position.set(L.r.x0, L.top, 0);
    rails.forEach((r, i) => {
      const y = (r.y = i ? L.r.ym : L.r.yp);
      r.line.setPoints([[L.r.x0, y], [L.r.x1, y]]);
      r.caps[0].setPoints([[L.r.x0, y - 0.13], [L.r.x0, y + 0.13]]);
      r.caps[1].setPoints([[L.r.x1, y - 0.13], [L.r.x1, y + 0.13]]);
      r.name.position.set(L.r.x0, y + 0.3, 0);
      r.beads.forEach((d) => d.setPx(P.pp));
      r.vals.forEach((l, j) => l.position.set(RX(V[i][j]) + (i ? 0 : (j ? -1 : 1) * (P.pp / ppu)), y - 0.28, 0));
    });
    const rows = L.b.rows, yT = rows[0] + 0.45, yB = rows[3] - 0.45, dy = rows[0] - rows[1];
    bTitle.setText(L.bt).position.set(L.b.tx, L.b.ty, 0);
    bars.forEach((b, i) => {
      b.bar.setPoints([[XB(0), rows[i], 1], [XB(b.v), rows[i], 1]]).setWidth(px(0.36, 7, 15));
      b.name.position.set(XB(0) - 0.15, rows[i], 0);
      b.val.position.set(XB(b.v) + 0.12, rows[i], 0);
    });
    base.setPoints([[XB(0), yT], [XB(0), yB]]);
    one.setPoints([[XB(1), yT, 0.5], [XB(1), yB, 0.5]]);
    oneLab.position.set(XB(1), yT + 0.1, 0);
    band.setSize(XB(NB1) - XB(NB0), yT - yB, 0.08);
    band.position.set((XB(NB0) + XB(NB1)) / 2, (yT + yB) / 2, -0.5);
    bandLab.position.set(XB(NB1), yB - 0.12, 0);
    ghosts.forEach((gh, i) => { gh.y = rows[i] - dy / 2; gh.a.setPx(px(0.06, 2.2, 3.5)); gh.b.setPx(px(0.06, 2.2, 3.5)); });
    const FY = L.C[1] - 1.0, fx = L.C[0] + 1.6 * L.S, hl = L.f.hx - 1.25, r7 = 7 / ppu;
    fArrow.set([fx, FY, 1], [hl - 0.08, FY, 1]);
    head.position.set(L.f.hx, FY, 0.5);
    headLab.position.set(L.f.hx, FY, 2);
    oArrow.set([L.f.hx + 1.33, FY, 1], [L.f.hx + 3.0, FY, 1]);
    outLab.position.set(L.f.hx + 3.15, FY, 0);
    const J = [lerp(fx, hl, 0.45), FY, 1];
    [[J[0] + L.f.l1, FY + 2.7], [L.f.hx, FY + 1.6]].forEach(([x, y], i) => {
      const lp = lamps[i];
      lp.br.setPoints([J, [x, y - r7 - 0.04, 1]]);
      lp.ring.position.set(x, y, 2); lp.ring.setPx(7);
      lp.bulb.position.set(x, y, 2); lp.bulb.setPx(3.2);
      lp.l.position.set(x, y + r7 + 0.12, 0);
    });
    const py = FY - 1.45;
    probe.line.setPoints([[J[0], py + r7 + 0.04, 1], J]);
    probe.ring.position.set(J[0], py, 2); probe.ring.setPx(7);
    probe.core.position.set(J[0], py, 2); probe.core.setPx(3);
    probe.l.position.set(J[0], py - r7 - 0.1, 0);
  }

  const tmp = [0, 0, 0];
  function update(t) {
    const u = ctx.loopT(t, PERIOD), ppu = P.ppu;
    const th = 0.5 * Math.sin(2 * Math.PI * u);
    ct = Math.cos(th); st = Math.sin(th);
    /* beat 1 */
    const dim = 1 - 0.3 * pulse(u, 0.14, 0.97, 0.06);
    for (const q of pts) { proj(q.p, q.d.position); q.d.setOpacity(dim * clamp(0.58 + 0.24 * q.d.position.z, 0.22, 0.9)); }
    cl.forEach((c) => c.l.position.set(L.C[0] + L.S * c.x * ct, L.C[1] + L.S * 1.15, 0));
    const k = seg(u, 0.04, 0.2) * (1 - seg(u, 0.93, 0.99));
    pair.forEach((c) => {
      const m = [0.02, 0.05 + c.sgn * P.mh, 0.05];
      for (let i = 0; i < 3; i++) tmp[i] = lerp(c.home[i], m[i], k);
      proj(tmp, c.w).setZ(5);
      c.d.position.copy(c.w);
      c.back.position.set(c.w.x, c.w.y, 4.9);
    });
    const A = pair[0].w, B = pair[1].w, ak = pulse(u, 0.17, 0.935, 0.03);
    const bow = 0.5 * A.distanceTo(B) + 7 / ppu;
    arcPts.forEach((p, i) => { const s = i / 12; p.set(lerp(A.x, B.x, s) + bow * 4 * s * (1 - s), lerp(A.y, B.y, s), 4.95); });
    arc.setPoints(arcPts).setOpacity(ak);
    if (L.cosRight) cosLab.position.set(Math.max(A.x, B.x) + bow + 6 / ppu, (A.y + B.y) / 2, 6);
    else cosLab.position.set((A.x + B.x) / 2, Math.min(A.y, B.y) - (P.pp + 5) / ppu, 6);
    cosLab.setOpacity(ak);
    /* beat 2: copies of the two cells fly to the rails; the measured pair springs apart */
    const off = 1 - seg(u, 0.665, 0.7, ease.inOutSine);
    const rk = seg(u, 0.22, 0.29);
    rTitle.setOpacity(seg(u, 0.22, 0.26) * off);
    rails.forEach((r, i) => {
      r.line.setProgress(rk).setOpacity(0.5 * off);
      r.caps.forEach((c) => c.setOpacity(0.5 * rk * off));
      r.name.setOpacity(seg(u, 0.24, 0.28) * off);
      const t0 = 0.26 + 0.03 * i, fly = seg(u, t0, t0 + 0.07);
      const sp = i ? ease.outBack(clamp((u - 0.37) / 0.08)) : 0;
      r.beads.forEach((d, j) => {
        const x = RX(lerp(V[0][j], V[i][j], sp));
        d.position.set(lerp(pair[j].w.x, x, fly), lerp(pair[j].w.y, r.y, fly), 5);
        d.setOpacity(seg(u, t0, t0 + 0.015) * off);
      });
      const [p0, p1] = r.beads.map((d) => d.position);
      r.gap.setPoints([[p0.x, p0.y, 4.8], [p1.x, p1.y, 4.8]]).setOpacity(0.3 * seg(u, t0 + 0.05, t0 + 0.08) * off);
      r.vals.forEach((l) => l.setOpacity(seg(u, 0.33 + 0.12 * i, 0.36 + 0.12 * i) * off));
    });
    /* beat 3 */
    const bk = seg(u, 0.44, 0.48) * off;
    bTitle.setOpacity(bk); oneLab.setOpacity(bk); base.setOpacity(bk); one.setOpacity(0.9 * bk);
    bars.forEach((b, i) => {
      const a0 = 0.47 + 0.035 * i;
      b.bar.setProgress(seg(u, a0, a0 + 0.06, ease.outCubic)).setOpacity(off);
      b.name.setOpacity(bk);
      b.val.setOpacity(seg(u, a0 + 0.05, a0 + 0.08) * off);
    });
    band.setOpacity(0.13 * seg(u, 0.56, 0.61) * off);
    bandLab.setOpacity(seg(u, 0.57, 0.61) * off);
    ghosts.forEach((gh, i) => {
      const gk = seg(u, 0.56 + 0.015 * i, 0.6 + 0.015 * i), y = gh.y + 0.3 * (1 - gk), o = 0.65 * gk * off;
      gh.a.position.set(XB(gh.c - 0.12), y, 2); gh.b.position.set(XB(gh.c + 0.12), y, 2);
      gh.l.setPoints([[XB(gh.c - 0.12), y, 1.9], [XB(gh.c + 0.12), y, 1.9]]).setOpacity(o);
      gh.a.setOpacity(o); gh.b.setOpacity(o);
    });
    /* beat 4: the fork stays open; neither lamp lights */
    const fo = 1 - seg(u, 0.92, 0.965, ease.inOutSine);
    fArrow.setProgress(seg(u, 0.7, 0.75)).setOpacity(fo);
    const hk = seg(u, 0.72, 0.76) * fo;
    head.setOpacity(hk); headLab.setOpacity(hk);
    const ok = L.f.out ? fo : 0;
    oArrow.setProgress(seg(u, 0.75, 0.79)).setOpacity(ok);
    outLab.setOpacity(seg(u, 0.77, 0.8) * ok);
    lamps.forEach((lp, i) => {
      const lk = seg(u, 0.77 + 0.02 * i, 0.81 + 0.02 * i);
      lp.br.setProgress(lk).setOpacity(fo);
      lp.ring.setOpacity(lk * fo); lp.bulb.setOpacity(0.3 * lk * fo); lp.l.setOpacity(lk * fo);
    });
    const pk = seg(u, 0.81, 0.85) * fo;
    probe.ring.setOpacity(pk); probe.core.setOpacity(pk); probe.l.setOpacity(pk);
    probe.line.setProgress(seg(u, 0.83, 0.87)).setOpacity(fo).setDashPhase(u, 12);
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: 0.65 * PERIOD,
    update, resize: layout,
    dispose() { pts.length = 0; pair.length = 0; },
  };
}
