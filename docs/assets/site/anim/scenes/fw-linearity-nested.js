/*
 * fw-linearity-nested — anm-framework.html #math-checks, "Two checks: linearity and nested readouts".
 * Left (linearity): a response ribbon O(u) around a base point u₀ with two secant rods, one through
 * u₀±ε and one through u₀±ε/2. The curve is odd around u₀, so both rods pivot on the base point.
 * Straight curve: the rods coincide and glow green (same slope). The curve bends (3 s lerp): the short
 * rod tilts less and the rods turn red; it bends further until the short rod tilts the other way.
 * The long rod's slope stays fixed throughout, so only the halved push changes its answer.
 * Right (nested readouts): an 'exactly right' disc inside a 'still works' disc on the floor. As the
 * load slider rises a ball drifts out of the centre, leaves the inner disc first (degraded: works,
 * not exact), then the outer one (broken). Shapes and load are illustrative, not data. 12 s loop.
 */
export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp, smooth } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const tmp = new THREE.Color(), tmp2 = new THREE.Color();
  const lab = (text, x, y, o) => { const l = ctx.label(text, o); l.position.set(x, y, 5); scene.add(l); return l; };

  /* headers + divider */
  const CX = -4.1, CY = -0.1, RX = 3.75, RY = -0.35;
  // titles sit at the left of each half, clear of the loader's play/pause button (top-right corner)
  lab('linearity', -7.55, 3.85, { size: 12, weight: 600, color: 'muted', anchor: 'left' });
  lab('nested readouts', 0.45, 3.85, { size: 12, weight: 600, color: 'muted', anchor: 'left' });
  scene.add(ctx.line([[0, -3.7], [0, 3.7]], { color: 'line', width: 1 }));

  /* ── left: O(u) = b·s + c·s³ with s = (u − u₀)/ε (odd, so every secant passes through u₀) ── */
  const W = 1.7, H = 1.1, S0 = 1.3, N = 64, AX = -6.85, AY = -2.75;
  const SHAPE_A = [0.62, 0], SHAPE_B = [0.2, 0.42], SHAPE_C = [-0.55, 1.17]; // long slope b + c = 0.62 in all three
  scene.add(ctx.arrow([AX, AY], [-1.35, AY], { color: 'faint', width: 1.5, head: 8 }));
  scene.add(ctx.arrow([AX, AY], [AX, 2.0], { color: 'faint', width: 1.5, head: 8 }));
  lab('u', -1.2, AY, { size: 12, color: 'muted', anchor: 'left' });
  lab('O(u)', AX, 2.12, { size: 12, color: 'muted', anchor: 'bottom' });
  lab('u₀', CX, AY - 0.18, { size: 12, color: 'muted', anchor: 'top' });
  [[0, 0.24], [-1, 0.16], [1, 0.16], [-0.5, 0.12], [0.5, 0.12]].forEach(([s, h]) => {
    scene.add(ctx.line([[CX + W * s, AY], [CX + W * s, AY + h]], { color: 'faint', width: 1.5 }));
  });
  scene.add(ctx.line([[CX, AY], [CX, CY]], { color: 'faint', width: 1, dashed: [3, 4] }));

  // ribbon: a faint band behind the curve (offset up-right) gives it depth
  const DX = 0.26, DY = 0.2;
  const ribPos = new Float32Array((N + 1) * 6);
  const ribGeo = ctx.track(new THREE.BufferGeometry());
  ribGeo.setAttribute('position', new THREE.BufferAttribute(ribPos, 3));
  const idx = [];
  for (let i = 0; i < N; i++) { const a = i * 2; idx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2); }
  ribGeo.setIndex(idx);
  const ribMat = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, opacity: 0.13, depthWrite: false, side: THREE.DoubleSide }));
  ctx.bind(ribMat, 'accent');
  const ribbon = new THREE.Mesh(ribGeo, ribMat);
  ribbon.frustumCulled = false;
  scene.add(ribbon);

  const curvePts = Array.from({ length: N + 1 }, () => new THREE.Vector3());
  const curve = ctx.line([[0, 0], [1, 0]], { color: 'accent', width: 3 });
  scene.add(curve);
  const rod = (w) => {
    const halo = ctx.line([[0, 0], [1, 0]], { width: 11, opacity: 0.2 });
    const core = ctx.line([[0, 0], [1, 0]], { width: w });
    scene.add(halo, core);
    return { halo, core };
  };
  const longRod = rod(2.5), shortRod = rod(3.5);
  const hits = [-1, -0.5, 0.5, 1].map((s) => { const d = ctx.dot([0, 0, 3], { px: 4, color: 'good' }); d.userData.s = s; scene.add(d); return d; });
  const base = ctx.dot([CX, CY, 3.2], { px: 4.5, color: 'ink' });
  scene.add(base);
  const epsL = lab('±ε', 0, 0, { size: 12, color: 'muted', anchor: 'left' });
  const halfL = lab('±ε/2', 0, 0, { size: 12, color: 'muted', anchor: 'top-left' });
  const sameL = lab('same slope', CX, 2.9, { size: 13, weight: 600, color: 'good' });
  const diffL = lab('slopes differ', CX, 2.9, { size: 13, weight: 600, color: 'bad' });
  const flipL = lab('sign flips', CX, 2.9, { size: 13, weight: 600, color: 'bad' });

  /* ── right: nested discs on the floor (y squashed by K), ball, load slider ── */
  const RO = 3.05, RI = 1.8, K = 0.52, RMAX = 1.25 * RO, TH = (150 * Math.PI) / 180, BR = 0.24;
  const floor = ctx.group();
  floor.position.set(RX, RY, 0);
  floor.scale.set(1, K, 1);
  floor.add(ctx.dot([0, 0, 0], { r: RO, color: 'accent', opacity: 0.1 }));
  floor.add(ctx.dot([0, 0, 0.1], { r: RI, color: 'accent', opacity: 0.2 }));
  const shadow = ctx.dot([0, 0, 0.3], { r: 0.26, color: 'ink', opacity: 0.16 });
  floor.add(shadow);
  scene.add(floor);
  const ellipse = (r) => Array.from({ length: 72 }, (_, i) => {
    const a = (i / 72) * Math.PI * 2;
    return [RX + r * Math.cos(a), RY + r * K * Math.sin(a), 0.2];
  });
  scene.add(ctx.line(ellipse(RO), { color: 'accent', width: 1.5, opacity: 0.6, closed: true }));
  scene.add(ctx.line(ellipse(RI), { color: 'accent', width: 2, closed: true }));
  lab('still works', RX, RY - RO * K - 0.14, { size: 11, color: 'muted', anchor: 'top' });
  lab('exactly right', RX, RY - 0.45, { size: 11, color: 'ink' });
  const ball = ctx.dot([RX, RY + BR, 4], { r: BR, color: 'good' });
  const shine = ctx.dot([RX, RY, 4.1], { r: 0.065, color: 'card', opacity: 0.45 });
  scene.add(ball, shine);

  const SX = 7.42, SY0 = -1.85, SY1 = 1.55;
  scene.add(ctx.line([[SX, SY0], [SX, SY1]], { color: 'line', width: 6 }));
  const fill = ctx.line([[SX, SY0, 0.1], [SX, SY1, 0.1]], { color: 'muted', width: 6 });
  const knob = ctx.dot([SX, SY0, 0.3], { px: 6.5, color: 'ink' });
  scene.add(fill, knob);
  lab('load', SX, SY0 - 0.25, { size: 11, color: 'muted', anchor: 'top' });

  const exactL = lab('exact', RX, 2.2, { size: 13, weight: 600, color: 'good', anchor: 'bottom' });
  const degL = lab('degraded', RX, 2.2, { size: 13, weight: 600, color: 'warn', anchor: 'bottom' });
  const degL2 = lab('works, not exact', RX, 2.13, { size: 11, color: 'muted', anchor: 'top' });
  const brokeL = lab('broken', RX, 2.2, { size: 13, weight: 600, color: 'bad', anchor: 'bottom' });

  const setRod = (r, m, e, col, glow) => {
    const pts = [[CX - W * e, CY - H * m * e, 2], [CX + W * e, CY + H * m * e, 2]];
    r.core.setPoints(pts).setColor(col);
    r.halo.setPoints(pts).setColor(col).setOpacity(glow).setWidth(lerp(7, 12, clamp((ctx.width - 360) / 700)));
  };

  return {
    scene, camera, period: PERIOD, still: 5.04,
    update(t) {
      const u = ctx.loopT(t, PERIOD);

      /* left: A (straight) → B (bent) → C (short rod flips) → back to A */
      const k1 = seg(u, 0.15, 0.40, ease.inOutSine), k2 = seg(u, 0.43, 0.68, ease.inOutSine), k3 = seg(u, 0.84, 1, ease.inOutSine);
      let b = lerp(SHAPE_A[0], SHAPE_B[0], k1), c = lerp(SHAPE_A[1], SHAPE_B[1], k1);
      b = lerp(b, SHAPE_C[0], k2); c = lerp(c, SHAPE_C[1], k2);
      b = lerp(b, SHAPE_A[0], k3); c = lerp(c, SHAPE_A[1], k3);
      const f = (s) => b * s + c * s * s * s;
      for (let i = 0; i <= N; i++) {
        const s = -S0 + (2 * S0 * i) / N, x = CX + W * s, y = CY + H * f(s);
        curvePts[i].set(x, y, 1);
        ribPos.set([x, y, 0, x + DX, y + DY, 0], i * 6);
      }
      ribGeo.getAttribute('position').needsUpdate = true;
      curve.setPoints(curvePts);

      const mL = b + c, mS = b + c / 4; // secant slopes through u₀±ε and u₀±ε/2
      const d = Math.abs(mL - mS) / Math.abs(mL);
      const kBad = smooth((d - 0.15) / 0.2), kFlip = smooth((0.04 - mS) / 0.08);
      tmp.copy(ctx.colors.good).lerp(ctx.colors.bad, kBad);
      const glow = 0.24 * (1 - kBad); // green rods glow; red rods do not
      setRod(longRod, mL, 1.45, tmp, glow);
      setRod(shortRod, mS, 0.75, tmp, glow);
      hits.forEach((h) => { const s = h.userData.s; h.position.set(CX + W * s, CY + H * f(s), 3); h.setColor(tmp); });
      epsL.position.set(CX + W * 1.45 + 0.12, CY + H * mL * 1.45, 5);
      halfL.position.set(CX + W * 0.75 + 0.08, CY + H * mS * 0.75 - 0.22, 5);
      sameL.setOpacity(1 - kBad);
      diffL.setOpacity(kBad * (1 - kFlip));
      flipL.setOpacity(kBad * kFlip);

      /* right: load rises to 'degraded', holds, rises to 'broken', holds, returns */
      const L = clamp(0.63 * seg(u, 0.08, 0.30, ease.inOutSine) + 0.37 * seg(u, 0.48, 0.64, ease.inOutSine) - seg(u, 0.84, 0.98, ease.inOutSine));
      const R = L * RMAX, fx = R * Math.cos(TH), fy = R * Math.sin(TH);
      shadow.position.set(fx, fy, 0.3);
      ball.position.set(RX + fx, RY + K * fy + BR, 4);
      shine.position.set(ball.position.x - 0.06, ball.position.y + 0.07, 4.1);
      const w = 0.14;
      const a0 = 1 - smooth((R - RI + w) / (2 * w)), a2 = smooth((R - RO + w) / (2 * w)), a1 = clamp(1 - a0 - a2);
      tmp2.copy(ctx.colors.good).lerp(ctx.colors.warn, 1 - a0).lerp(ctx.colors.bad, a2);
      ball.setColor(tmp2);
      exactL.setOpacity(a0);
      degL.setOpacity(a1);
      degL2.setOpacity(a1);
      brokeL.setOpacity(a2);
      fill.setProgress(L);
      knob.position.set(SX, lerp(SY0, SY1, L), 0.3);
    },
    dispose() {
      ctx.unbind(ribMat);
      ribGeo.dispose();
      ribMat.dispose();
    },
  };
}
