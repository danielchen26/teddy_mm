/*
 * idx-setup-mode-b-vs-a — index.html, #setup ("Two ways to study a foundation model"), above the
 * two mode cards. TEDDY is drawn as its 12 transformer layers (TEDDY-G config.json; the count is not
 * labelled) threaded by the residual column, h_{l+1} = h_l + f_l(h_l).
 * Left, lit: Mode B (this repo). RNA goes in at the bottom, a pulse runs up the column and out of
 * the top; the slabs stay dark, a closed box. Our head turns the embedding into predicted proteins;
 * ANM, outside TEDDY, reads the 9 it uses and gives a call or no call.
 * Right, ghosted in greys the whole time and tagged "not done": Mode A. The view swings over and
 * small probe rings slide in between the layers one by one, asking whether NK and T still differ.
 * Nothing is answered there: no values, no result. Wide stages show both halves side by side
 * (the swing is a small pan and a focus shift); narrow stages pan from one half to the other.
 * 12 s seamless loop.
 */
/* Every pipeline number this scene prints (stated on the page). Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  panel: 9,   // predicted proteins ANM reads
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 19.6, height: 6.5 });
  const PERIOD = 12, NL = 12, HALF = 4.8, SW = 1.4;   // layers, half width, slab width
  const XB = -HALF, XA = HALF;                         // half centres (world x)
  const XS_B = -2.3, XS_A = 0;                         // stack x inside each half (local)
  const HEAD_X = -0.1, DOT_X = 1.62, DOT_D = 0.22, ANM_X = 3.35;

  const gB = ctx.group(), gA = ctx.group();
  gB.position.x = XB; gA.position.x = XA;
  scene.add(gB, gA);
  const add = (g, ...o) => { g.add(...o); return o[0]; };

  function stack(g, xs, fill, stroke, fillOp, name, nameColor, sub) {
    const slabs = Array.from({ length: NL }, () =>
      add(g, ctx.box(SW, 0.12, { color: fill, opacity: fillOp, stroke, strokeWidth: 1, radius: 0.03 })));
    const col = add(g, ctx.line([[xs, 0, 0.5], [xs, 1, 0.5]], { color: stroke, width: 2 }));
    const rna = add(g, ctx.label('RNA', { size: 12, color: 'muted', anchor: 'top' }));
    const ted = add(g, ctx.label(name, { size: 13, weight: 600, color: nameColor, anchor: 'bottom-right' }));
    const tsub = add(g, ctx.label(sub, { size: 11, color: 'muted', anchor: 'top-right' }));
    return { xs, slabs, col, rna, ted, tsub, fillOp };
  }
  function header(g, title, titleColor, tag, tagColor) {
    const h = add(g, ctx.label(title, { size: 13, weight: 600, color: titleColor, anchor: 'top-left' }));
    const t = add(g, ctx.label(tag, { size: 11, weight: 600, color: tagColor, bg: 'soft', bgOpacity: 1, anchor: 'top-left' }));
    return { h, t };
  }

  // Mode B (this repo): TEDDY in its own colour, slabs dark; the flow leaves the top.
  const SB = stack(gB, XS_B, 'teddy', 'teddy', 0.1, 'TEDDY', 'teddy', 'closed box');
  const HB = header(gB, 'Mode B', 'ink', 'this repo', 'accent');
  const glow = add(gB, ctx.line([[0, 0], [0, 1]], { color: 'teddy', width: 8, opacity: 0.14 }));
  const pulse = add(gB, ctx.line([[0, 0], [0, 1]], { color: 'teddy', width: 4 }));
  const node = add(gB, ctx.dot([0, 0, 1], { r: 0.1, color: 'teddy' }));
  const head = add(gB, ctx.box(1.7, 0.6, { color: 'soft', stroke: 'muted', strokeWidth: 1.5, radius: 0.14 }));
  const headL = add(gB, ctx.label('our head', { size: 11, weight: 600 }));
  const anm = add(gB, ctx.box(1.5, 0.6, { color: 'accent', opacity: 0.12, stroke: 'accent', strokeWidth: 1.5, radius: 0.14 }));
  const anmL = add(gB, ctx.label('ANM', { size: 12, weight: 600 }));
  const dots = Array.from({ length: 9 }, () => add(gB, ctx.dot([0, 0, 1], { r: 0.07, color: 'teddy' })));
  const protL = add(gB, ctx.label(`${NUMBERS.panel} predicted proteins`, { size: 11, color: 'muted', anchor: 'top' }));
  const callL = add(gB, ctx.label('call / no call', { size: 11, weight: 600, anchor: 'top' }));
  const arrows = [0, 1, 2, 3].map(() => add(gB, ctx.arrow([0, 0, 1], [1, 0, 1], { color: 'muted', width: 1.5, head: 7 })));

  // Mode A (not done): the same stack in greys only; probe rings for the gaps between layers.
  const SA = stack(gA, XS_A, 'faint', 'faint', 0.08, 'TEDDY', 'muted', 'opened up');
  const HA = header(gA, 'Mode A', 'muted', 'first probes', 'muted');
  const rings = Array.from({ length: NL - 1 }, () => add(gA, ctx.dot([0, 0, 3], { r: 0.1, color: 'muted', hollow: true, ring: 0.34 })));
  const askL = add(gA, ctx.label('NK–T differ?', { size: 12, weight: 600, color: 'muted', anchor: 'left' }));

  const divider = add(scene, ctx.line([[0, -1], [0, 1]], { color: 'line', width: 1 }));

  let side = true;
  const L = { y0: 0, y1: 1, yTop: 1, sp: 0.4, gy: [] };

  function layout(w, h) {
    w = w || ctx.width || 640; h = h || ctx.height || 360;
    const a = w / h, ts = ctx.textScale;
    side = w >= 740;                                   // both halves fit; else pan between them
    const DW = side ? 19.6 : 10, hw = DW / 2, hh = hw / a;
    camera.userData.animFit = { width: DW, height: DW / a, fit: 'contain' };
    camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    const ppu = w / DW, P = (px) => px / ppu;
    const half = Math.min(hh, 3.6), T = half - 0.22;
    const hdr1 = P(13 * ts * 1.22 + 5), hdrH = hdr1 + P(11 * ts * 1.22 + 12);
    const clear = side ? 1e9 : half - P(40) - 0.66;   // narrow: keep ANM below the loader's corner button
    L.y1 = Math.min(T - hdrH - 0.12, clear);
    L.yTop = L.y1 + 0.36;
    L.y0 = -half + 0.15 + P(12 * ts * 1.22 + 6) + 0.26;
    L.sp = (L.y1 - L.y0) / (NL - 1);
    L.gy = Array.from({ length: NL - 1 }, (_, i) => L.y0 + (i + 0.5) * L.sp);
    const th = clamp(L.sp * 0.38, 0.06, 0.16), yMid = (L.y0 + L.y1) / 2, yc0 = L.y0 - 0.26;

    for (const S of [SB, SA]) {
      S.slabs.forEach((b, i) => { b.setSize(SW, th, Math.min(0.03, th / 2)); b.position.set(S.xs, L.y0 + i * L.sp, 0); });
      S.col.setPoints([[S.xs, yc0, 0.5], [S.xs, L.yTop, 0.5]]);
      S.rna.position.set(S.xs, yc0 - P(2), 2);
      S.ted.position.set(S.xs - SW / 2 - 0.2, yMid + P(1), 2);
      S.tsub.position.set(S.xs - SW / 2 - 0.2, yMid - P(1), 2);
    }
    for (const [H, x] of [[HB, -4.45], [HA, -4.4]]) {
      H.h.position.set(x, T, 2);
      H.t.position.set(x - P(5), T - hdr1, 2);
    }
    glow.setPoints([[XS_B, yc0, 0.4], [XS_B, L.yTop, 0.4]]);
    pulse.setPoints([[XS_B, yc0, 0.6], [XS_B, L.yTop, 0.6]]);
    const y = L.yTop;
    node.position.set(XS_B, y, 1);
    head.position.set(HEAD_X, y, 0); headL.position.set(HEAD_X, y, 2);
    anm.position.set(ANM_X, y, 0); anmL.position.set(ANM_X, y, 2);
    dots.forEach((d, i) => d.position.set(DOT_X + ((i % 3) - 1) * DOT_D, y + (1 - Math.floor(i / 3)) * DOT_D, 1));
    protL.position.set(DOT_X, y - 0.38, 2);
    arrows[0].set([XS_B + 0.16, y, 1], [HEAD_X - 0.9, y, 1]);
    arrows[1].set([HEAD_X + 0.9, y, 1], [DOT_X - DOT_D - 0.13, y, 1]);
    arrows[2].set([DOT_X + DOT_D + 0.13, y, 1], [ANM_X - 0.8, y, 1]);
    arrows[3].set([ANM_X, y - 0.34, 1], [ANM_X, y - 1.3, 1]);
    callL.position.set(ANM_X, y - 1.36, 2);
    rings.forEach((r) => r.setRadius(clamp(L.sp * 0.4, 0.07, 0.13)));
    divider.setPoints([[0, -half + 0.2, 0], [0, half - 0.2, 0]]);
  }
  layout();

  return {
    scene, camera, period: PERIOD, still: 0.86 * PERIOD,
    resize(w, h) { layout(w, h); },
    update(t) {
      let u = ctx.loopT(t, PERIOD);
      const still = ctx.reducedMotion;
      if (still) u = side ? 0.86 : 0.445;            // narrow still: Mode B, flow complete
      const sw = seg(u, 0.45, 0.55) - seg(u, 0.9, 1.0, ease.inOutSine);   // 0 = Mode B in view, 1 = Mode A
      let fB = 1 - 0.55 * sw, fA = 0.5 + 0.5 * sw;
      let cx = side ? lerp(-0.2, 0.2, sw) : lerp(XB, XA, sw);
      if (still && side) { fB = 1; fA = 1; cx = 0; }
      camera.position.x = cx;
      divider.setOpacity(side ? 1 : 0);

      // Mode B: dark slabs, the column carries the pass; the flow leaves TEDDY at the top.
      SB.slabs.forEach((b) => { b.material.opacity = SB.fillOp * fB; b.outline.setOpacity(0.4 * fB); });
      SB.col.setOpacity(0.75 * fB); glow.setOpacity(0.14 * fB);
      for (const o of [SB.rna, SB.ted, SB.tsub, HB.h, HB.t, headL, anmL]) o.setOpacity(fB);
      const outB = 1 - seg(u, 0.87, 0.94);
      const p = seg(u, 0.04, 0.19, ease.inOutSine) * 1.25;
      pulse.setProgress(clamp(p), clamp(p - 0.25));
      pulse.setOpacity(fB);
      node.setOpacity(seg(u, 0.18, 0.21) * outB * fB);
      const hk = seg(u, 0.23, 0.27) * outB, ak = seg(u, 0.36, 0.4) * outB;
      head.material.opacity = fB; head.outline.setOpacity(fB * (0.45 + 0.55 * hk));
      anm.material.opacity = fB * (0.1 + 0.16 * ak); anm.outline.setOpacity(fB * (0.55 + 0.45 * ak));
      [[0.2, 0.24], [0.25, 0.28], [0.33, 0.36], [0.39, 0.42]].forEach(([a0, a1], i) => {
        const k = seg(u, a0, a1);
        arrows[i].setProgress(k);
        arrows[i].setOpacity(Math.min(1, k * 6) * outB * fB);
      });
      dots.forEach((d, i) => d.setOpacity(seg(u, 0.27 + i * 0.006, 0.3 + i * 0.006) * outB * fB));
      protL.setOpacity(seg(u, 0.3, 0.33) * outB * fB);
      callL.setOpacity(seg(u, 0.41, 0.44) * outB * fB);

      // Mode A: greys only, always tagged; probes slide in between the layers, bottom to top.
      SA.slabs.forEach((b) => { b.material.opacity = SA.fillOp * fA; b.outline.setOpacity(0.55 * fA); });
      SA.col.setOpacity(0.7 * fA);
      for (const o of [SA.rna, SA.ted, SA.tsub, HA.h]) o.setOpacity(fA);
      HA.t.setOpacity(Math.max(fA, 0.85));
      const outA = 1 - seg(u, 0.9, 0.97);
      const q = clamp((u - 0.56) / 0.024, 0, NL - 2);   // gap being asked about
      rings.forEach((r, i) => {
        const k = seg(u, 0.56 + i * 0.024, 0.62 + i * 0.024, ease.outCubic);
        r.position.set(lerp(XS_A + SW / 2 + 0.25, XS_A, k), L.gy[i], 3);
        r.setOpacity(Math.min(1, k * 3) * (1 - 0.5 * clamp(q - i)) * outA * fA);
      });
      askL.position.set(XS_A + SW / 2 + 0.45, L.y0 + (q + 0.5) * L.sp, 2);
      askL.setOpacity(seg(u, 0.55, 0.59) * outA * fA);
    },
    dispose() {
      // Every geometry, material and label texture lives in the scene; the core disposes them.
      SB.slabs.length = 0; SA.slabs.length = 0; dots.length = 0; rings.length = 0; arrows.length = 0;
    },
  };
}
