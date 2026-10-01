/*
 * fw-flock-readouts — anm-framework.html #physics-flock, "Physics II: two flock readouts fail apart".
 * A slightly tilted periodic box holds 256 arrowheads moving at one constant speed. A fixed field arrow
 * points right along the top edge; a translucent wedge (±tolerance) fans out from the box centre around it.
 * Noise η (slider, bottom left) rises and falls once per 12 s; headings spread wider around a direction
 * that stays pinned to the field. Arrows whose heading lies inside the wedge are teal; the Q bar
 * (alignment) is that share and the centre arrow's length is P (polarization), both computed from the
 * spread the arrows are drawn from (the 256 on screen follow it closely, not count for count).
 * A dashed ½ line crosses both bars. Q falls below ½ at the page's η for the tolerance shown
 * (20°: 0.25, 30°: 0.35, 45°: 0.50), P at 0.60. The chips cycle 30° → 45° → 20° (36 s in all), with
 * the page's note that the order can reverse near 45.5°.
 * Illustrative, not data: headings follow a symmetric bounded spread whose median |heading| grows
 * linearly with η and whose shape puts every ½ crossing on the η the page states; the bar heights
 * between crossings, the speed and the arrow motion are schematic. No other numbers are shown.
 */
export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp, smooth } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const LOOP = 12, TOLS = [30, 45, 20], CHIPS = [20, 30, 45], QETA = { 20: '0.25', 30: '0.35', 45: '0.50' };
  const N = 256, G = 0.957, D2R = Math.PI / 180, TAU = Math.PI * 2;
  const rnd = ctx.rand(11);
  const lab = (text, x, y, o) => { const l = ctx.label(text, o); l.position.set(x, y, 8); scene.add(l); return l; };

  /* η over one loop: hold low, rise, hold high, fall back */
  const etaAt = (u) => 0.15 + 0.55 * (seg(u, 0.06, 0.66, ease.inOutSine) - seg(u, 0.76, 0.94, ease.inOutSine));
  // half-width of the heading spread (radians): median |heading| = (100η − 5)°, so Q_tol = ½ at η = (tol + 5)/100
  const halfW = (eta) => Math.max(0.5, 100 * eta - 5) * Math.pow(2, G) * D2R;
  const polar = (w) => { let s = 0; for (let i = 0; i < 48; i++) s += Math.cos(w * Math.pow((i + 0.5) / 48, G)); return s / 48; };
  const share = (tolDeg, w) => Math.min(1, Math.pow((tolDeg * D2R) / w, 1 / G));
  const erf = (x) => {
    const t = 1 / (1 + 0.3275911 * Math.abs(x));
    const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return x < 0 ? -y : y;
  };
  const KS = [2, 3, 5];
  const heading = (ph, u, w) => {
    let g = 0;
    for (let j = 0; j < 3; j++) g += Math.sin(TAU * KS[j] * u + ph[j]);
    const p = erf(g / Math.sqrt(3)); // g has variance 1.5; erf(z/√2) spreads it ~evenly over (−1, 1)
    return w * Math.sign(p) * Math.pow(Math.abs(p), G);
  };

  /* trajectories: integrate each arrow along its heading, pick noise phases whose loop closes in the box */
  const S = 5.9, M = 72, dt = LOOP / M;
  const wm = Array.from({ length: M }, (_, m) => halfW(etaAt((m + 0.5) / M)));
  const V = S / (LOOP * wm.reduce((a, w) => a + polar(w) / M, 0)); // about one box length per loop
  const arrows = [];
  for (let i = 0; i < N; i++) {
    let best = null;
    for (let c = 0; c < 64 && !(best && best.score < 0.3); c++) {
      const ph = [rnd() * TAU, rnd() * TAU, rnd() * TAU];
      const X = new Float32Array(M + 1), Y = new Float32Array(M + 1);
      for (let m = 0; m < M; m++) {
        const th = heading(ph, (m + 0.5) / M, wm[m]);
        X[m + 1] = X[m] + V * dt * Math.cos(th);
        Y[m + 1] = Y[m] + V * dt * Math.sin(th);
      }
      const rx = X[M] - Math.round(X[M] / S) * S, ry = Y[M] - Math.round(Y[M] / S) * S;
      const score = Math.max(Math.abs(rx), Math.abs(ry));
      if (!best || score < best.score) best = { ph, X, Y, rx, ry, score };
    }
    for (let m = 0; m <= M; m++) { best.X[m] -= (best.rx * m) / M; best.Y[m] -= (best.ry * m) / M; }
    best.x0 = (rnd() - 0.5) * S; best.y0 = (rnd() - 0.5) * S;
    arrows.push(best);
  }

  /* box (tilted: y foreshortened, a thin front face) */
  const BX = -3.2, BY = 0.42, K = 0.82, RW = 2.3;
  const floor = ctx.group();
  floor.position.set(BX, BY, 0);
  floor.scale.set(1, K, 1);
  scene.add(floor);
  floor.add(ctx.box(S, S, { color: 'soft', stroke: 'line', radius: 0.12 }));
  const front = ctx.box(S - 0.1, 0.2, { color: 'line', radius: 0.06 });
  front.position.set(BX, BY - (S / 2) * K - 0.09, -0.1);
  scene.add(front);
  const FY = BY + (S / 2) * K + 0.32;
  scene.add(ctx.arrow([BX - S / 2 + 0.3, FY], [BX + S / 2 - 0.3, FY], { color: 'muted', width: 2, head: 9 }));
  lab('field', BX - S / 2 + 0.15, FY, { size: 11, color: 'muted', anchor: 'right' });

  /* wedge: fan + outline in the box plane */
  const FAN = 40;
  const fanPos = new Float32Array((FAN + 2) * 3);
  const fanGeo = ctx.track(new THREE.BufferGeometry());
  fanGeo.setAttribute('position', new THREE.BufferAttribute(fanPos, 3));
  fanGeo.setIndex(Array.from({ length: FAN }, (_, i) => [0, i + 1, i + 2]).flat());
  const fanMat = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, opacity: 0.15, depthWrite: false, side: THREE.DoubleSide }));
  ctx.bind(fanMat, 'accent');
  const fan = new THREE.Mesh(fanGeo, fanMat);
  fan.position.z = 0.5;
  fan.frustumCulled = false;
  floor.add(fan);
  const fanPts = Array.from({ length: FAN + 2 }, () => new THREE.Vector3(0, 0, 0.6));
  const fanEdge = ctx.line(fanPts, { color: 'accent', width: 1.5, opacity: 0.8, closed: true });
  floor.add(fanEdge);
  let fanA = -1;
  const setFan = (a) => {
    if (Math.abs(a - fanA) < 1e-5) return;
    fanA = a;
    for (let i = 0; i <= FAN; i++) {
      const b = -a + (2 * a * i) / FAN;
      fanPos[(i + 1) * 3] = RW * Math.cos(b); fanPos[(i + 1) * 3 + 1] = RW * Math.sin(b);
      fanPts[i + 1].set(RW * Math.cos(b), RW * Math.sin(b), 0.6);
    }
    fanGeo.getAttribute('position').needsUpdate = true;
    fanGeo.computeBoundingSphere();
    fanEdge.setPoints(fanPts);
  };
  const wedgeL = lab('within 30°', 0, 0, { size: 11, color: 'accent', weight: 600, anchor: 'bottom-left', bg: 'card', bgOpacity: 0.85 });

  /* the flock: one instanced arrowhead */
  const geo = ctx.track(new THREE.BufferGeometry());
  geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array([0.55, 0, 0, -0.45, 0.36, 0, -0.22, 0, 0, -0.45, -0.36, 0]), 3));
  geo.setIndex([0, 1, 2, 0, 2, 3]);
  const mat = ctx.track(new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, depthWrite: false, side: THREE.DoubleSide }));
  const flock = new THREE.InstancedMesh(geo, mat, N);
  flock.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  flock.frustumCulled = false;
  flock.position.z = 2;
  floor.add(flock);
  const mtx = new THREE.Matrix4(), col = new THREE.Color(), tmp = new THREE.Color();

  /* mean heading at the centre */
  const halo = ctx.line([[BX, BY, 3.9], [BX + RW, BY, 3.9]], { color: 'card', width: 8, opacity: 0.8 });
  const mean = ctx.arrow([BX, BY, 4], [BX + RW, BY, 4], { color: 'ink', width: 3, head: 11 });
  scene.add(halo, mean);
  lab('mean\nheading', BX - 0.2, BY, { size: 11, color: 'ink', anchor: 'right', bg: 'card', bgOpacity: 0.85 });

  /* bars: Q (alignment) and P (polarization), dashed ½ line */
  const QX = 2.5, PX = 5.4, BW = 0.95, B0 = -1.45, BH = 4.3, I0 = B0 + 0.08, IH = BH - 0.16, HALF = I0 + IH / 2;
  const bar = (x, c) => {
    const track = ctx.box(BW, BH, { color: 'soft', stroke: 'line', radius: 0.08 });
    track.position.set(x, B0 + BH / 2, 0);
    const fill = ctx.box(BW - 0.16, 1, { color: c, radius: 0 });
    fill.position.set(x, I0, 0.2);
    scene.add(track, fill);
    return fill;
  };
  const qFill = bar(QX, 'accent'), pFill = bar(PX, 'ink');
  pFill.setOpacity(0.86);
  scene.add(ctx.line([[QX - BW / 2 - 0.25, HALF, 1], [PX + BW / 2 + 0.25, HALF, 1]], { color: 'ink', width: 1.5, dashed: [5, 4], opacity: 0.75 }));
  lab('½', QX - BW / 2 - 0.36, HALF, { size: 13, color: 'ink', weight: 600, anchor: 'right' });
  const rows = [ // [label, row] under the bars; rows are spaced in screen px (see layout)
    [lab('Q', QX, 0, { size: 13, weight: 700, color: 'accent', anchor: 'top' }), 0],
    [lab('P', PX, 0, { size: 13, weight: 700, color: 'ink', anchor: 'top' }), 0],
    [lab('alignment', QX, 0, { size: 11, color: 'muted', anchor: 'top' }), 1],
    [lab('polarization', PX, 0, { size: 11, color: 'muted', anchor: 'top' }), 1],
    [lab('below ½', QX, 0, { size: 11, weight: 600, color: 'bad', anchor: 'top' }), 2],
    [lab('below ½', PX, 0, { size: 11, weight: 600, color: 'bad', anchor: 'top' }), 2],
    [lab('order can reverse near 45.5°', (QX + PX) / 2, 0, { size: 11, color: 'muted', anchor: 'top' }), 3.1],
  ];
  const qLow = rows[4][0], pLow = rows[5][0], note = rows[6][0];

  /* tolerance chips (top left, clear of the corner play/pause button) */
  const chips = CHIPS.map((d) => {
    const b = ctx.box(1, 1, { color: 'accent', stroke: 'line', strokeWidth: 1.5, radius: 0.2 });
    const l = lab(d + '°', 0, 0, { size: 12, weight: 600, color: 'ink' });
    scene.add(b);
    return { d, b, l };
  });
  const layout = () => {
    const ppu = ctx.ppu(), w = Math.max(0.95, 30 / ppu), h = Math.max(0.42, 14 / ppu), r = Math.max(0.3, 16 / ppu);
    chips.forEach((c, i) => {
      const x = BX - S / 2 + w / 2 + i * (w + 0.12);
      c.b.setSize(w, h, Math.min(0.2, h / 2));
      c.b.position.set(x, 3.95, 0); c.l.position.set(x, 3.95, 8);
    });
    rows.forEach(([l, k]) => { l.position.y = B0 - 0.12 - k * r; });
    const g = Math.max(0.32, 11 / ppu); // slider numbers clear the knob at phone width
    // the number under the slider must not drop off the stage's bottom edge on a phone (it was cut by 2 px)
    const lhN = (11 * ctx.textScale * 1.22 + 2) / ppu;
    qNum.position.y = SY + g; pNum.position.y = SY - Math.min(g, Math.max(7 / ppu, SY + 4.5 - lhN - 1 / ppu));
  };

  /* η slider under the box: only the page's crossing values are labelled */
  const SX0 = -5.5, SX1 = -0.9, SY = -3.45;
  const ex = (eta) => lerp(SX0, SX1, (eta - 0.15) / 0.55);
  scene.add(ctx.line([[SX0, SY], [SX1, SY]], { color: 'line', width: 5 }));
  const sFill = ctx.line([[SX0, SY, 0.1], [SX1, SY, 0.1]], { color: 'muted', width: 5 });
  const knob = ctx.dot([SX0, SY, 0.6], { px: 6.5, color: 'ink' });
  lab('noise η', SX0 - 0.25, SY, { size: 11, color: 'muted', anchor: 'right' });
  const tick = (c) => { const l = ctx.line([[0, SY - 0.24, 0.3], [0, SY + 0.24, 0.3]], { color: c, width: 2 }); scene.add(l); return l; };
  const qTick = tick('accent'), pTick = tick('ink');
  pTick.position.x = ex(0.6);
  const qNum = lab('0.35', 0, SY + 0.32, { size: 11, weight: 600, color: 'accent', anchor: 'bottom' });
  const pNum = lab('0.60', ex(0.6), SY - 0.32, { size: 11, weight: 600, color: 'ink', anchor: 'top' });
  scene.add(sFill, knob);
  layout();

  const wrap = (v) => v - S * Math.floor(v / S + 0.5);
  const low = (v) => smooth((0.5 - v) / 0.016 + 0.5);

  return {
    scene, camera, period: LOOP * 3, still: 4.8,
    update(t) {
      const U = ctx.loopT(t, LOOP * 3), cyc = Math.min(2, Math.floor(U * 3)), u = U * 3 - cyc;
      const kT = seg(u, 0.95, 1, ease.inOutSine), tolA = TOLS[cyc], tolB = TOLS[(cyc + 1) % 3];
      const tol = lerp(tolA, tolB, kT), shown = kT < 0.5 ? tolA : tolB;
      const eta = etaAt(u), w = halfW(eta), Q = share(tol, w), P = polar(w);

      /* flock */
      const ppu = ctx.ppu(), sz = clamp(0.17 * ppu, 7, 13) / ppu, tr = tol * D2R, mm = u * M;
      const m0 = Math.min(M - 1, Math.floor(mm)), fr = mm - m0;
      arrows.forEach((a, i) => {
        const th = heading(a.ph, u, w);
        const x = wrap(a.x0 + lerp(a.X[m0], a.X[m0 + 1], fr)), y = wrap(a.y0 + lerp(a.Y[m0], a.Y[m0 + 1], fr));
        const s = sz * smooth(Math.min(S / 2 - Math.abs(x), S / 2 - Math.abs(y)) / 0.35);
        const c = Math.cos(th) * s, n = Math.sin(th) * s;
        mtx.set(c, -n, 0, x, n, c, 0, y, 0, 0, 1, 0, 0, 0, 0, 1);
        flock.setMatrixAt(i, mtx);
        col.copy(ctx.colors.faint).lerp(ctx.colors.accent, smooth((tr - Math.abs(th)) / (3 * D2R) + 0.5));
        flock.setColorAt(i, col);
      });
      flock.instanceMatrix.needsUpdate = true;
      if (flock.instanceColor) flock.instanceColor.needsUpdate = true;

      /* wedge, mean heading */
      setFan(tr);
      wedgeL.setText('within ' + shown + '°');
      wedgeL.position.set(BX + RW * Math.cos(tr) + 0.05, BY + K * RW * Math.sin(tr) + 0.08, 8);
      const tip = [BX + Math.max(0.02, P) * RW, BY, 4];
      mean.set([BX, BY, 4], tip);
      halo.setPoints([[BX, BY, 3.9], [tip[0] - 0.1, BY, 3.9]]);

      /* bars */
      [[qFill, qLow, Q], [pFill, pLow, P]].forEach(([b, l, v]) => { b.scale.y = Math.max(1e-3, v * IH); b.position.y = I0 + (v * IH) / 2; l.setOpacity(low(v)); });

      /* chips + note */
      chips.forEach((c) => {
        const act = (c.d === tolA ? 1 - kT : 0) + (c.d === tolB ? kT : 0);
        c.b.material.opacity = 0.2 * act;
        tmp.copy(ctx.colors.line).lerp(ctx.colors.accent, act);
        c.b.outline.setColor(tmp);
        c.l.setOpacity(0.5 + 0.5 * act);
        if (c.d === 45) note.setOpacity(act);
      });

      /* slider */
      sFill.setProgress((eta - 0.15) / 0.55);
      knob.position.x = ex(eta);
      const qEta = (tol + 5) / 100;
      qTick.position.x = ex(qEta);
      qTick.setWidth(2 + 2 * (1 - smooth(Math.abs(eta - qEta) / 0.03)));
      pTick.setWidth(2 + 2 * (1 - smooth(Math.abs(eta - 0.6) / 0.03)));
      qNum.setText(QETA[shown]);
      qNum.position.x = ex(qEta);
    },
    resize() { layout(); },
    dispose() {
      ctx.unbind(fanMat);
      flock.dispose();
    },
  };
}
