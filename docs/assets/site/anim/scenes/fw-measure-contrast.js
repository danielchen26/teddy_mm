/*
 * fw-measure-contrast — anm-framework.html, #ai-response ("AI step 4: measure the response").
 * One event's value sits on a rail at its baseline. The knob is pushed to +ε and −ε (long arrows),
 * then to +ε/2 and −ε/2 (short arrows). Each push is propagated through the graph; the final state
 * is a ghost copy of the graph, and the state moved, Δz = z+ − z−, is the vector between the two
 * ghosts. The fresh model readout of each final state is a score mark (s+, s−); the readout moved is
 * the arrow between the two marks. State and readout are measured separately. The graph update is
 * linear, so the state vector for ε/2 is exactly half the one for ε; the readout need not be.
 * The slopes ĝ(ε) = (s+ − s−)/2ε and ĝ(ε/2) are two rods through one pivot (long = ε, short = ε/2);
 * the grey wedge is the page's 25% band around ĝ(ε). Direction 1 agrees (green); in direction 2 the
 * short rod tilts the other way (sign reversal, red). A tally of 12 directions then fills:
 * state moved 12/12 (imposed by the propagation), readout moved 9/12, slopes agree 2/12.
 * Only those counts and the 25% tolerance are real (page / paper values); the scores, slopes and the
 * order of the tiles are schematic. 12 s seamless loop; the still frame shows direction 2 and the full tally.
 */

/* Numbers shown on the canvas. ANM graph study (paper sections/ai_evidence.tex; stated in the
   #ai-response text and its sources), not TEDDY-pipeline numbers. */
const NUMBERS = {
  directions: 12,    // perturbation directions measured
  stateMoved: 12,    // state change in 12 of 12 (by construction: propagation is imposed)
  readoutMoved: 9,   // fresh model readout changed in 9 of 12
  slopesAgree: 2,    // passed the ε-versus-ε/2 step-halving check in 2 of 12
  tolerancePct: 25,  // |ĝ(ε) − ĝ(ε/2)| / max{|ĝ(ε)|, |ĝ(ε/2)|} ≤ 0.25
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, x, y, o) => {
    const l = ctx.label(text, Object.assign({ size: 11, color: 'muted' }, o));
    l.position.set(x, y, 10);
    return add(l);
  };

  /* ── layout (world units, 16 × 9) ── */
  const Y0 = 0.85, HY = 3.45;                          // panel centre line, heading row
  const RX = -6.35, E = 1.5;                           // rail x, knob travel for ε
  const AXL = RX + 0.4, AXS = RX + 0.7;                // long (ε) and short (ε/2) push arrows
  const KS = [1, 0.5];                                 // step sizes in units of ε
  const XB = [-4.2, -3.0], G = 1.1, GW = 0.8, GH = 0.44; // state sub-columns (ε, ε/2), gain, ghost size
  const XC = [-1.0, 0.2], MW = 0.5;                    // readout sub-columns, mark width
  const PX = 3.55, R1 = 1.5, RW = 1.75, KT = 1.1;       // rod pivot x, long rod half length, band radius, slope → tan
  const D = NUMBERS.directions, TP = 0.62, TX0 = -((D - 1) * TP) / 2 - 0.05;
  const TR = [-2.15, -2.77, -3.39];                    // tally rows
  const TOL = NUMBERS.tolerancePct / 100;

  /* ── two worked directions (schematic scores; ε = 1 rail unit) ── */
  const EX = [
    { tile: 0, a: 0.02, L: 0.37, ds: [1.2, 0.66], mid: [0.04, 0.02] },   // slopes agree
    { tile: 1, a: 0.40, L: 0.31, ds: [0.9, -0.5], mid: [-0.05, 0.0] },  // short step reverses the sign
  ];
  EX.forEach((e) => {
    const g1 = e.ds[0] / 2, g2 = e.ds[1] / (2 * 0.5);                   // ĝ = (s+ − s−) / 2ε
    e.ok = Math.abs(g1 - g2) / Math.max(Math.abs(g1), Math.abs(g2)) <= TOL;
    e.flip = Math.sign(g1) !== Math.sign(g2);
    e.th = [Math.atan(KT * g1), Math.atan(KT * g2)];
    const b0 = Math.atan(KT * g1 * (1 - TOL)), b1 = Math.atan((KT * g1) / (1 - TOL));
    e.band = [Math.min(b0, b1), Math.max(b0, b1)];                    // slopes within 25% of ĝ(ε)
    e.verdict = e.ok ? 'agree' : e.flip ? 'sign flips' : 'slopes differ';
    e.sub = `${e.ok ? 'within' : 'outside'} ${NUMBERS.tolerancePct}%`;
    e.s = KS.map((k, j) => [Y0 + e.mid[j] + e.ds[j] / 2, Y0 + e.mid[j] - e.ds[j] / 2]); // [s+, s−] per step
  });

  /* ── tally pattern: counts from NUMBERS; tile order illustrative ── */
  const changed = new Array(D).fill(true), agree = new Array(D).fill(false);
  for (let i = 0; i < D - NUMBERS.readoutMoved; i++) changed[clamp(6 + i, 2, D - 1)] = false;
  agree[0] = true;
  [4, 10, 3, 9, 5, 11, 2, 7, 8, 6].filter((i) => i < D && changed[i]).slice(0, NUMBERS.slopesAgree - 1).forEach((i) => { agree[i] = true; });
  const FILL0 = 0.72, FDT = 0.0125;
  const tileT = (i, r) => {
    const e = EX.find((q) => q.tile === i);
    return e ? e.a + e.L * (0.9 + 0.016 * r) : FILL0 + (i - 2) * FDT + 0.004 * r;
  };

  /* ── static frame ── */
  lab('event value', RX, HY, { size: 12, weight: 600 });
  lab('state moved', (XB[0] + XB[1]) / 2, HY, { size: 12, weight: 600, color: 'accent' });
  lab('readout moved', (XC[0] + XC[1]) / 2, HY, { size: 12, weight: 600, color: 'ink' });
  // short heading + key: the loader's play/pause pill covers the top-right corner on phones
  lab('slopes', PX, HY, { size: 12, weight: 600 });
  const keyL = add(ctx.line([[0, 0, 1], [1, 0, 1]], { color: 'ink', width: 2.4 })), keyLt = lab('ε', 0, 0, { anchor: 'left' });
  const keyS = add(ctx.line([[0, 0, 1], [1, 0, 1]], { color: 'ink', width: 4 })), keySt = lab('ε/2', 0, 0, { anchor: 'left' });
  add(ctx.line([[RX, Y0 - E - 0.3, 0.2], [RX, Y0 + E + 0.3, 0.2]], { color: 'line', width: 3 }));
  [[1, '+ε'], [0.5, '+ε/2'], [0, 'base'], [-0.5, '−ε/2'], [-1, '−ε']].forEach(([k, s]) => {
    const w = k ? 0.13 : 0.2;
    add(ctx.line([[RX - w, Y0 + E * k, 0.3], [RX + w, Y0 + E * k, 0.3]], { color: k ? 'faint' : 'muted', width: 1.4 }));
    lab(s, RX - 0.3, Y0 + E * k, { anchor: 'right', color: k ? 'muted' : 'ink' });
  });
  ['ε', 'ε/2'].forEach((s, j) => {
    lab(s, XB[j], Y0 - G - GH - 0.3, { anchor: 'top' });
    lab(s, XC[j], Y0 - G - GH - 0.3, { anchor: 'top' });
    add(ctx.line([[XC[j], Y0 - 1.55, 0.2], [XC[j], Y0 + 1.55, 0.2]], { color: 'line', width: 1.2 }));
  });
  add(ctx.line([[PX - RW, Y0, 0.2], [PX + RW, Y0, 0.2]], { color: 'faint', width: 1, dashed: [4, 4] }));
  add(ctx.dot([PX, Y0, 4], { px: 3.5, color: 'ink' }));
  const ROWS = [['state moved', 'accent'], ['readout moved', 'ink'], ['slopes agree', 'good']];
  ROWS.forEach(([s, c], r) => lab(s, TX0 - 0.42, TR[r], { anchor: 'right', color: c }));
  const TY = (TR[0] + TR[2]) / 2, TH = TR[0] - TR[2] + 0.52;
  for (let i = 0; i < D; i++) add(ctx.box(0.5, TH, { color: 'soft', radius: 0.1 })).position.set(TX0 + i * TP, TY, 0.2);
  lab(`${D} directions · order illustrative`, TX0 + ((D - 1) * TP) / 2, TY - TH / 2 - 0.2, { anchor: 'top', color: 'faint', size: 10.5 });

  /* ── event knob and push arrows ── */
  const knobBack = add(ctx.dot([RX, Y0, 3.8], { px: 8.5, color: 'card' }));
  const knobFill = add(ctx.dot([RX, Y0, 3.9], { px: 6.5, color: 'accent', opacity: 0.35 }));
  const knobRing = add(ctx.dot([RX, Y0, 4], { px: 7.5, color: 'accent', hollow: true, ring: 0.28 }));
  const push = KS.map((k, j) => [1, -1].map((sg) => add(ctx.arrow([j ? AXS : AXL, Y0, 1], [j ? AXS : AXL, Y0 + sg * E * k, 1], { color: 'muted', width: j ? 2.4 : 2, head: 8 }))));

  /* ── state: ghost copies of the graph and the Δz vector between them ── */
  const GP = [[-0.22, -0.07], [0.02, 0.09], [0.23, -0.06]];
  const ghost = () => {
    const g = add(ctx.group());
    const box = ctx.box(GW, GH, { color: 'card', stroke: 'accent', strokeWidth: 1.4, radius: 0.08 });
    g.add(box);
    const edges = [[0, 1], [1, 2], [0, 2]].map(([a, b]) => {
      const l = ctx.line([[GP[a][0], GP[a][1], 0.05], [GP[b][0], GP[b][1], 0.05]], { color: 'accent', width: 1 });
      g.add(l);
      return l;
    });
    const nodes = GP.map((p, i) => { const d = ctx.dot([p[0], p[1], 0.1], { px: i ? 2.8 : 2.2, color: 'accent' }); g.add(d); return d; });
    return { g, nodes, set(o) { box.setOpacity(o); edges.forEach((l) => l.setOpacity(0.55 * o)); nodes.forEach((d) => d.setOpacity(o)); } };
  };
  const ghosts = KS.map((k, j) => [1, -1].map((sg) => {
    const gh = ghost();
    gh.g.position.set(XB[j], Y0 + sg * (G * k + GH / 2), 3);
    return gh;
  }));
  const dz = KS.map((k, j) => add(ctx.arrow([XB[j], Y0 - G * k, 2], [XB[j], Y0 + G * k, 2], { color: 'accent', width: 2.2, head: 8 })));
  const dzLab = lab('Δz', XB[0] - 0.14, Y0, { anchor: 'right', color: 'accent', weight: 600 });

  /* ── readout: score marks s+ / s− and the arrow between them ── */
  const marks = KS.map((k, j) => [0, 1].map((s) => ({
    m: add(ctx.line([[XC[j] - MW / 2, Y0, 2], [XC[j] + MW / 2, Y0, 2]], { color: 'ink', width: 2.6 })),
    t: lab(s ? 's−' : 's+', XC[j] + MW / 2 + 0.06, Y0, { anchor: 'left' }),
  })));
  const dsA = KS.map((k, j) => add(ctx.arrow([XC[j], Y0 - 0.5, 2.2], [XC[j], Y0 + 0.5, 2.2], { color: 'ink', width: 2, head: 7 })));

  /* ── slopes: two rods through one pivot and the tolerance band ── */
  const rod = (w) => ({ halo: add(ctx.line([[0, 0, 1.9], [1, 0, 1.9]], { width: 10, opacity: 0 })), core: add(ctx.line([[0, 0, 2], [1, 0, 2]], { width: w })) });
  const rodL = rod(2.4), rodS = rod(4);
  const NW = 12, wPos = new Float32Array(2 * NW * 9);
  const wGeo = ctx.track(new THREE.BufferGeometry());
  wGeo.setAttribute('position', new THREE.BufferAttribute(wPos, 3));
  const wMat = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, opacity: 0, depthWrite: false, side: THREE.DoubleSide }));
  ctx.bind(wMat, 'muted');
  const wedge = new THREE.Mesh(wGeo, wMat);
  wedge.frustumCulled = false;
  wedge.position.z = 0.5;
  add(wedge);
  const wEdges = [0, 1].map(() => add(ctx.line([[0, 0, 0.6], [1, 0, 0.6]], { color: 'muted', width: 1 })));
  const bandLab = lab(`${NUMBERS.tolerancePct}%`, 0, 0, { anchor: 'left' });
  const verdict = lab('', PX, Y0 - 1.45, { size: 13, weight: 600 });
  const vsub = lab('', PX, Y0 - 1.85, { anchor: 'top' });
  const rp = (th, r) => [[PX - r * Math.cos(th), Y0 - r * Math.sin(th), 2], [PX + r * Math.cos(th), Y0 + r * Math.sin(th), 2]];

  /* ── tally cells, active-tile frame, counts ── */
  const cells = Array.from({ length: D }, (_, i) => [true, changed[i], agree[i]].map((on, r) => {
    const col = ROWS[r][1];
    const d = on ? ctx.dot([0, 0, 1], { px: 4.4, color: col }) : ctx.dot([0, 0, 1], { px: 4.4, color: 'faint', hollow: true, ring: 0.3 });
    d.position.set(TX0 + i * TP, TR[r], 1);
    return add(d);
  }));
  const hi = add(ctx.box(0.5, TH, { color: null, stroke: 'muted', strokeWidth: 1.4, radius: 0.1 }));
  const CX = TX0 + (D - 1) * TP + 0.45;
  const counts = [NUMBERS.stateMoved, NUMBERS.readoutMoved, NUMBERS.slopesAgree].map((n, r) =>
    lab(`${n}/${D}`, CX, TR[r], { anchor: 'left', color: ROWS[r][1], mono: true, weight: 600, size: 12 }));
  const imposed = lab('imposed', CX, TR[0], { anchor: 'left', color: 'faint' });
  const MONO = (getComputedStyle(ctx.figure).getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace';
  const meas = document.createElement('canvas').getContext('2d');
  const countW = (px) => { if (!meas) return 3 * px; meas.font = `600 ${px}px ${MONO}`; return meas.measureText(`${NUMBERS.stateMoved}/${D}`).width; };

  /* ── sizes that follow the stage (px-constant marks would crowd a phone-width stage) ── */
  function layout() {
    const ppu = ctx.ppu(), px = (k, lo, hi) => clamp(k * ppu, lo, hi);
    ghosts.forEach((row) => row.forEach((gh) => gh.nodes.forEach((d, i) => d.setPx(px(i ? 0.042 : 0.034, i ? 1.5 : 1.2, i ? 2.8 : 2.2)))));
    cells.forEach((col) => col.forEach((d) => d.setPx(px(0.066, 3, 4.6))));
    const hw = lerp(6, 10, clamp((ctx.width - 360) / 700));
    rodL.halo.setWidth(hw); rodS.halo.setWidth(hw + 1);
    // key under the 'slopes' heading, laid out in CSS px: long thin rod = ε, short thick rod = ε/2
    // (kept left of x = 4.6 so the paused pill cannot cover it on a phone-width stage)
    const k = (v) => v / ppu, ky = HY - k(19 * ctx.textScale), x0 = Math.min(PX - k(37), 4.6 - k(75));
    keyL.setPoints([[x0, ky, 1], [x0 + k(18), ky, 1]]); keyLt.position.set(x0 + k(21), ky, 10);
    keyS.setPoints([[x0 + k(42), ky, 1], [x0 + k(52), ky, 1]]); keySt.position.set(x0 + k(55), ky, 10);
    const kr = px(0.11, 5.5, 7.5);
    knobRing.setPx(kr); knobFill.setPx(kr - 1); knobBack.setPx(kr + 1);
    vsub.position.y = Y0 - 1.45 - (13 * ctx.textScale * 0.62 + 3) / ppu;
  }
  let sizeKey = '';

  /* ── timing inside one direction (p = 0..1) ── */
  const KF = [[0, 0], [0.09, 1], [0.15, 1], [0.27, -1], [0.46, -1], [0.53, 0.5], [0.58, 0.5], [0.66, -0.5], [0.7, -0.5], [0.78, 0]];
  const knobAt = (p) => {
    if (p <= 0 || p >= 0.78) return 0;
    for (let i = 1; i < KF.length; i++) {
      if (p <= KF[i][0]) { const [p0, v0] = KF[i - 1], [p1, v1] = KF[i]; return lerp(v0, v1, ease.inOutSine((p - p0) / (p1 - p0))); }
    }
    return 0;
  };
  const ARR = [[0.09, 0.27], [0.53, 0.66]];           // knob arrives at +, − for ε and ε/2
  const pulses = [add(ctx.dot([0, 0, 6], { px: 4, color: 'accent', opacity: 0 })), add(ctx.dot([0, 0, 6], { px: 4, color: 'ink', opacity: 0 }))];
  const colState = { rodL: '', rodS: '', verdict: -1 };
  let shown = -1;

  function setExample(e) {
    KS.forEach((k, j) => {
      marks[j].forEach((mk, s) => {
        const y = e.s[j][s];
        mk.m.setPoints([[XC[j] - MW / 2, y, 2], [XC[j] + MW / 2, y, 2]]);
        mk.t.position.y = y;
      });
      dsA[j].set([XC[j], e.s[j][1], 2.2], [XC[j], e.s[j][0], 2.2]);
    });
    let o = 0;
    for (const side of [0, Math.PI]) {
      for (let i = 0; i < NW; i++) {
        const b0 = side + lerp(e.band[0], e.band[1], i / NW), b1 = side + lerp(e.band[0], e.band[1], (i + 1) / NW);
        wPos.set([PX, Y0, 0, PX + RW * Math.cos(b0), Y0 + RW * Math.sin(b0), 0, PX + RW * Math.cos(b1), Y0 + RW * Math.sin(b1), 0], o);
        o += 9;
      }
    }
    wGeo.attributes.position.needsUpdate = true;
    e.band.forEach((b, i) => wEdges[i].setPoints(rp(b, RW).map((q) => [q[0], q[1], 0.6])));
    const bm = (e.band[0] + e.band[1]) / 2;
    bandLab.position.set(PX + (RW + 0.1) * Math.cos(bm), Y0 + (RW + 0.1) * Math.sin(bm), 10);
    verdict.setText(e.verdict);
    vsub.setText(e.sub);
    const vc = e.ok ? 'good' : 'bad';
    if (colState.verdict !== vc) { verdict.setColor(vc); colState.verdict = vc; }
  }
  const setRodColor = (key, r, c) => { if (colState[key] !== c) { r.core.setColor(c); r.halo.setColor(c); colState[key] = c; } };

  function update(t) {
    const key = `${ctx.width}x${ctx.height}`;
    if (key !== sizeKey) { sizeKey = key; layout(); }
    const u = ctx.loopT(t, PERIOD);
    const endFade = 1 - seg(u, 0.955, 0.99, ease.inOutSine);
    const ci = u < EX[1].a ? 0 : 1, e = EX[ci];
    if (ci !== shown) { setExample(e); shown = ci; }
    const p = clamp((u - e.a) / e.L, 0, 1);
    const vis = ci === 0 ? 1 - seg(p, 0.955, 0.995, ease.inOutSine) : endFade;

    // knob and push arrows (tips ride the knob)
    const kn = knobAt(u < e.a ? 0 : p);
    const ky = Y0 + E * kn;
    knobBack.position.y = ky; knobFill.position.y = ky; knobRing.position.y = ky;
    const prog = [
      [p < 0.15 ? clamp(kn) : 1, p < 0.15 ? 0 : p < 0.3 ? clamp(-kn) : 1],
      [p < 0.46 ? 0 : p < 0.58 ? clamp(kn / 0.5) : 1, p < 0.58 ? 0 : p < 0.7 ? clamp(-kn / 0.5) : 1],
    ];
    push.forEach((pr, j) => pr.forEach((a, s) => { const q = u < e.a ? 0 : prog[j][s]; a.setProgress(q).setOpacity(q > 0.02 ? vis : 0); }));

    // pushes → ghost states → score marks; the two vectors
    let pa = null, pb = null;
    KS.forEach((k, j) => {
      [0, 1].forEach((s) => {
        const tA = ARR[j][s], sg = s ? -1 : 1;
        const gy = Y0 + sg * (G * k + GH / 2), sy = e.s[j][s];
        ghosts[j][s].set(seg(p, tA + 0.03, tA + 0.05) * vis);
        const mo = seg(p, tA + 0.07, tA + 0.09) * vis;
        marks[j][s].m.setOpacity(mo);
        marks[j][s].t.setOpacity(mo);
        if (p > tA && p < tA + 0.04) pa = { f: (p - tA) / 0.04, x0: (j ? AXS : AXL) + 0.2, y0: Y0 + sg * E * k, x1: XB[j] - GW / 2, y1: gy };
        if (p > tA + 0.04 && p < tA + 0.08) pb = { f: (p - tA - 0.04) / 0.04, x0: XB[j] + GW / 2, y0: gy, x1: XC[j] - MW / 2, y1: sy };
      });
      const tm = ARR[j][1];
      const zq = seg(p, tm + 0.04, tm + 0.09, ease.inOutSine), rq = seg(p, tm + 0.08, tm + 0.13, ease.inOutSine);
      dz[j].setProgress(zq).setOpacity(vis);
      dsA[j].setProgress(rq).setOpacity(vis);
    });
    dzLab.setOpacity(seg(p, 0.34, 0.37) * vis);
    [pa, pb].forEach((q, i) => {
      if (!q || u < e.a) { pulses[i].setOpacity(0); return; }
      pulses[i].position.set(lerp(q.x0, q.x1, ease.inOutSine(q.f)), lerp(q.y0, q.y1, ease.inOutSine(q.f)), 6);
      pulses[i].setOpacity(Math.sin(Math.PI * q.f) * vis);
    });

    // slope rods, band, verdict
    const vk = seg(p, 0.875, 0.9), done = p >= 0.875;
    const thL = e.th[0] * seg(p, 0.4, 0.46, ease.inOutSine);
    rodL.core.setPoints(rp(thL, R1)).setOpacity(seg(p, 0.4, 0.42) * vis);
    rodL.halo.setPoints(rp(thL, R1)).setOpacity(0.2 * vk * vis);
    const thS = lerp(e.th[0], e.th[1], seg(p, 0.8, 0.87, ease.inOutSine));
    rodS.core.setPoints(rp(thS, R1 / 2)).setOpacity(seg(p, 0.79, 0.81) * vis);
    rodS.halo.setPoints(rp(thS, R1 / 2)).setOpacity(0.24 * vk * vis);
    const vc = e.ok ? 'good' : 'bad';
    setRodColor('rodL', rodL, done ? vc : 'muted');
    setRodColor('rodS', rodS, done ? vc : 'ink');
    const bk = seg(p, 0.43, 0.48) * vis;
    wMat.opacity = 0.12 * bk;
    wEdges.forEach((l) => l.setOpacity(0.45 * bk));
    bandLab.setOpacity(bk);
    verdict.setOpacity(vk * vis);
    vsub.setOpacity(vk * vis);

    // tally
    for (let i = 0; i < D; i++) cells[i].forEach((d, r) => { const t0 = tileT(i, r); d.setOpacity(seg(u, t0, t0 + 0.012) * endFade); });
    let hx = null;
    if (u >= e.a && u < e.a + e.L * (ci ? 1 : 0.95)) hx = e.tile;
    else if (u >= FILL0 && u < FILL0 + (D - 2) * FDT) hx = 2 + Math.floor((u - FILL0) / FDT);
    hi.setOpacity(hx == null ? 0 : 0.9);
    if (hx != null) hi.position.set(TX0 + hx * TP, TY, 0.4);
    const ck = seg(u, 0.845, 0.875) * endFade;
    counts.forEach((c) => c.setOpacity(ck));
    imposed.position.x = CX + (countW(12 * ctx.textScale) + 6) / ctx.ppu();
    imposed.setOpacity(ck);
  }

  return {
    scene, camera, period: PERIOD, still: 0.9 * PERIOD,
    update,
    dispose() { cells.length = 0; ghosts.length = 0; },
  };
}
