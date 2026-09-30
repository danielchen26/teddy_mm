/*
 * idx-claims-timestep — index.html, #timestep ("ANM’s step is an order, not necessarily a clock"),
 * above the two .ts2 cards.
 * Left, NOT A REAL ORDER: the nine markers of one prediction, in panel order (CD19 first, CD36 last).
 * Row 1, as the first run did it: they enter ANM’s field one per step, and every step each earlier
 * bead keeps only part of its glow (DECAY, the first run’s per-step retention from the spec; it only
 * shades the beads and is never shown), so at read-out the early beads are dim and the late ones
 * bright: a hidden weight. Row 2, the fix: all nine enter at the same instant and glow equally.
 * Right, A REAL ORDER: TEDDY’s 12 transformer layers (TEDDY-G config; the count is not labelled) on
 * the residual stream, the retained state. A small perturbation δh travels up and is reshaped at
 * each layer by (I + ∂f/∂h); the δh shapes are illustrative, not computed. Footnote: not yet
 * tested; pre-registered Mode A. Wide stages show both halves side by side; narrow stages pan from
 * the left half to the right one. World units are CSS px. 12 s seamless loop.
 */
export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 640, height: 213 });
  const PERIOD = 12, NL = 12, DECAY = 0.82;
  const MARKERS = ['CD19', 'CD72', 'CD22', 'CD3', 'CD2', 'CD5', 'CD16', 'CD11c', 'CD36'];
  const E0 = 0.05, ES = 0.03, EW = 0.02, E2 = 0.4;   // row-1 entry times (one per step), row-2 entry
  const gL = ctx.group(), gR = ctx.group();
  scene.add(gL, gR);
  const add = (g, o) => { g.add(o); return o; };
  const lab = (g, text, o) => add(g, ctx.label(text, o));

  const header = (g, tag, tagColor, title) => ({
    tag: lab(g, tag, { size: 10, mono: true, weight: 600, color: tagColor, anchor: 'top-left' }),
    title: lab(g, title, { size: 13, weight: 600, color: 'ink', anchor: 'top-left' }),
  });
  const HL = header(gL, 'NOT A REAL ORDER', 'bad', 'Panel order');
  const HR = header(gR, 'A REAL ORDER', 'good', 'TEDDY’s layers');

  // Left: two trays in ANM's field, nine slots each.
  const rows = ['one per step', 'all at once'].map((name, r) => ({
    tray: add(gL, ctx.box(10, 10, { color: 'accent', opacity: 0.06, stroke: 'accent', strokeWidth: 1, radius: 5 })),
    name: lab(gL, name, { size: 11, color: 'muted', anchor: 'bottom-left' }),
    verdict: lab(gL, r ? 'equal weight' : 'hidden weight', { size: 12, weight: 600, color: r ? 'good' : 'bad', anchor: 'left' }),
    slots: MARKERS.map(() => add(gL, ctx.dot([0, 0, 1], { r: 5, color: 'faint', hollow: true, ring: 0.22, opacity: 0.7 }))),
    halos: MARKERS.map(() => add(gL, ctx.dot([0, 0, 2], { r: 9, color: 'accent', opacity: 0 }))),
    beads: MARKERS.map(() => add(gL, ctx.dot([0, 0, 3], { r: 5, color: 'accent', opacity: 0 }))),
    xs: [], cy: 0,
  }));
  const anmL = lab(gL, 'ANM field', { size: 11, weight: 600, color: 'accent', anchor: 'bottom-right' });
  const names = MARKERS.map((m) => lab(gL, m, { size: 10, mono: true, color: 'muted', anchor: 'top' }));

  // Right: layer slabs on the residual stream; δh rides up it.
  const slabs = Array.from({ length: NL }, () =>
    add(gR, ctx.box(10, 4, { color: 'teddy', opacity: 0.14, stroke: 'teddy', strokeWidth: 1, radius: 2 })));
  const col = add(gR, ctx.arrow([0, 0, 1], [0, 1, 1], { color: 'teddy', width: 2, head: 7, opacity: 0.35 }));
  const fill = add(gR, ctx.line([[0, 0, 1.5], [0, 1, 1.5]], { color: 'teddy', width: 3.5, opacity: 0.8 }));
  const dh = add(gR, ctx.arrow([0, 0, 4], [1, 0, 4], { color: 'teddy', width: 2.5, head: 8 }));
  const dhL = lab(gR, 'δh', { size: 13, weight: 600, color: 'teddy', anchor: 'left' });
  const resL = lab(gR, 'residual stream', { size: 11, color: 'muted', anchor: 'bottom' });
  const jac = lab(gR, 'each layer\n× (I + ∂f/∂h)', { size: 11, color: 'ink', anchor: 'right' });
  const lead = add(gR, ctx.line([[0, 0, 1], [1, 0, 1]], { color: 'muted', width: 1, dashed: [3, 3] }));
  const foot = lab(gR, 'not yet tested; pre-registered Mode A', { size: 11, color: 'muted', anchor: 'bottom' });
  const divider = add(scene, ctx.line([[0, -1], [0, 1]], { color: 'line', width: 1 }));

  // Illustrative δh direction / length after each layer (fixed seed; not computed from TEDDY).
  const rnd = ctx.rand(11);
  const shapes = Array.from({ length: NL + 1 }, (_, k) => (k ? { a: (rnd() - 0.5) * 1.1, l: 0.55 + rnd() * 0.45 } : { a: 0.25, l: 0.75 }));
  shapes[NL] = { a: -0.15, l: 0.9 };   // last one points level, clear of the label above

  const L = { side: true, w: 640, r: 5, drop: 10, xs: 0, SW: 80, y0: 0, lsp: 10, yIn: 0, yOut: 1, jx: 0, yMid: 0 };

  function layout(w, h) {
    w = w || ctx.width || 640; h = h || ctx.height || 213;
    const ts = ctx.textScale, lh = (s) => s * ts * 1.22;
    L.side = w >= 660; L.w = w;
    const PW = L.side ? w / 2 : w, hw = PW / 2;
    camera.userData.animFit = { width: w, height: h, fit: 'contain' };
    camera.left = -w / 2; camera.right = w / 2; camera.top = h / 2; camera.bottom = -h / 2;
    camera.updateProjectionMatrix();
    gL.position.x = L.side ? -w / 4 : 0; gR.position.x = L.side ? w / 4 : w;
    const pad = clamp(h * 0.06, 12, 22), top = h / 2 - pad, x0 = -hw + pad;
    for (const H of [HL, HR]) { H.tag.position.set(x0, top, 5); H.title.position.set(x0, top - lh(10) - 3, 5); }
    const hdrB = top - lh(10) - 3 - lh(13);

    // Left half: header, then two row blocks centred in what is left.
    const xa = x0, xb = hw - pad - 96 * ts - 10;
    const sp = (xb - xa) / MARKERS.length, r = clamp(sp * 0.19, 4, 9), trayH = 2 * r + 12;
    const allNames = sp >= 34 * ts;
    const blk1 = lh(11) + 6 + trayH + 4 + lh(10), blk2 = lh(11) + 6 + trayH;
    const spanTop = hdrB - 10, spanBot = -h / 2 + pad;
    const gap = clamp((spanTop - spanBot - blk1 - blk2) * 0.5, 10, 48);
    let y = (spanTop + spanBot) / 2 + (blk1 + gap + blk2) / 2;
    rows.forEach((R, ri) => {
      const cy = y - lh(11) - 6 - trayH / 2;
      R.cy = cy;
      R.name.position.set(xa, y - lh(11), 5);
      R.tray.setSize(xb - xa, trayH, trayH / 2);
      R.tray.position.set((xa + xb) / 2, cy, 0);
      R.verdict.position.set(xb + 10, cy, 5);
      R.xs = MARKERS.map((_, i) => xa + (i + 0.5) * sp);
      R.slots.forEach((d, i) => { d.setRadius(r); d.position.set(R.xs[i], cy, 1); });
      R.beads.forEach((d) => d.setRadius(r));
      R.halos.forEach((d) => d.setRadius(r * 1.9));
      if (ri === 0) {
        anmL.position.set(xb, y - lh(11), 5);
        names.forEach((n, i) => { n.position.set(R.xs[i], cy - trayH / 2 - 4, 5); n.visible = allNames || i === 0 || i === 8; });
        y = cy - trayH / 2 - 4 - lh(10) - gap;
      }
    });
    L.r = r; L.drop = Math.min(14, r * 2);

    // Right half: residual column from yIn to yOut, 12 slabs on it, labels around it.
    const SW = clamp(PW * 0.2, 56, 120), xs = PW * 0.04;
    const yOut = top - lh(11) - 4, yIn = -h / 2 + pad + lh(11) + 10;
    const lsp = (yOut - yIn) / (NL - 1 + 1.8), th = clamp(lsp * 0.4, 3, 9);
    Object.assign(L, { xs, SW, lsp, yIn, yOut, y0: yIn + lsp * 0.8 });
    slabs.forEach((b, i) => { b.setSize(SW, th, Math.min(3, th / 2)); b.position.set(xs, L.y0 + i * lsp, 0); });
    col.set([xs, yIn, 1], [xs, yOut, 1]);
    fill.setPoints([[xs, yIn, 1.5], [xs, yOut, 1.5]]);
    resL.position.set(xs, yOut + 4, 5);
    L.jx = xs - SW / 2 - 24; L.yMid = L.y0 + ((NL - 1) * lsp) / 2;
    jac.position.set(L.jx, L.yMid, 5);
    foot.position.set(0, -h / 2 + pad, 5);
    divider.setPoints([[0, -h / 2 + pad, 0], [0, h / 2 - pad, 0]]);
    divider.visible = L.side;
  }
  layout();

  // Row 1: number of later entries that have happened since bead i entered (continuous).
  const stepsAfter = (u, i) => {
    let n = 0;
    for (let j = i + 1; j < MARKERS.length; j++) n += seg(u, E0 + j * ES, E0 + j * ES + EW, ease.inOutSine);
    return n;
  };

  return {
    scene, camera, period: PERIOD, still: 0.85 * PERIOD,
    resize(w, h) { layout(w, h); },
    update(t) {
      let u = ctx.loopT(t, PERIOD);
      if (ctx.reducedMotion) u = L.side ? 0.85 : 0.49;   // side: both halves done; narrow: left half done
      const out = 1 - seg(u, 0.9, 0.97, ease.inOutSine);
      const outL = L.side ? out : 1 - seg(u, 0.6, 0.64);  // narrow: left empties while off screen
      camera.position.x = L.side ? 0 : L.w * seg(u, 0.52, 0.6, ease.inOutCubic) * (1 - seg(u, 0.9, 0.99, ease.inOutCubic));

      rows.forEach((R, ri) => {
        R.beads.forEach((d, i) => {
          const a = ri ? E2 : E0 + i * ES;
          const k = seg(u, a, a + EW, ease.outCubic), vis = k * outL;
          const g = ri ? 1 : Math.pow(DECAY, stepsAfter(u, i));
          const yy = R.cy + (1 - k) * L.drop;
          d.position.set(R.xs[i], yy, 3);
          d.setOpacity(vis * Math.max(0.1, g));
          R.halos[i].position.set(R.xs[i], yy, 2);
          R.halos[i].setOpacity(vis * 0.22 * g);
        });
        R.verdict.setOpacity((ri ? seg(u, 0.45, 0.49) : seg(u, 0.32, 0.36)) * outL);
      });

      // δh from the bottom of the column to just past the top slab, reshaped as it crosses each slab.
      const rv = seg(u, 0.56, 0.6) * out;
      const s = seg(u, 0.58, 0.84, ease.inOutSine);
      const y = lerp(L.yIn, L.y0 + (NL - 0.65) * L.lsp, s), q = (y - L.y0) / L.lsp;
      let kf = 0;
      for (let k = 0; k < NL; k++) kf += ctx.smooth(clamp((q - k) * 2 + 0.5));
      const i0 = Math.min(NL - 1, Math.floor(kf)), f = kf - i0;
      const A = lerp(shapes[i0].a, shapes[i0 + 1].a, f), len = lerp(shapes[i0].l, shapes[i0 + 1].l, f) * L.SW * 0.5;
      const tip = [L.xs + Math.cos(A) * len, y + Math.sin(A) * len, 4];
      dh.set([L.xs, y, 4], tip);
      dh.setOpacity(rv);
      dhL.position.set(tip[0] + 5, tip[1], 5);
      dhL.setOpacity(rv);
      fill.setProgress(clamp((y - L.yIn) / (L.yOut - L.yIn)));
      fill.setOpacity(0.8 * rv);
      slabs.forEach((b, k) => b.setOpacity(0.14 + 0.36 * rv * Math.max(0, 1 - Math.abs(q - k) / 0.7)));
      const ly = clamp(y, L.y0, L.y0 + (NL - 1) * L.lsp);
      lead.setPoints([[L.jx + 5, L.yMid, 1], [L.xs - L.SW / 2 - 4, ly, 1]]);
      lead.setOpacity(0.9 * rv * (1 - seg(u, 0.84, 0.88)));
    },
    // Every geometry, material and label texture lives in the scene; the core disposes them.
    dispose() {},
  };
}
