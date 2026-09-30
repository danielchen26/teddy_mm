/*
 * loop-hero-lookalike (anm-loop.html #overview, between the lede and the four-step strip; the page keeps its
 * "withdrawn · first run" badge over the top-left corner, so nothing is drawn there).
 * A myeloid cell and a T cell drift together in TEDDY's embedding until they almost overlap (cosine 0.982).
 * The readout built on the embedding (our head + the soft rule) stamps both with the same call, and a gauge
 * fills to "falsely agree 85%". The strip's four steps then light in order: Expect (the readout's assumption
 * pulses), Veto (the same call is struck out), Complement (measured protein replaces the embedding as the
 * evidence, strict rule), Verify (the pair is told apart and the gauge falls to 40%, the 85% mark kept).
 * Numbers are the page's first-run values (0.982, 0.85 -> 0.40); positions are illustrative. 12 s loop.
 * Two layouts: the page's 21:7 strip, and 16:9 on narrow screens (the page's --ar-sm).
 */
export default function create(ctx) {
  const { THREE, ease, seg, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 18, height: 6 });
  const PERIOD = 12, BEFORE = 0.85, AFTER = 0.4;
  const BEATS = ['Expect', 'Veto', 'Complement', 'Verify'], BEAT_T = [0.4, 0.48, 0.56, 0.64];
  /* reg / tag = [cx, cy, w, h]; g = gauge [x0, x1, y]; bt = step chips [y, h, w, spacing]; r = cell radius */
  const LAY = {
    wide: { fw: 18, fh: 6, reg: [-5.5, 0.55, 6.2, 3.3], tag: [0.25, 0.72, 2.6, 0.85], g: [2.8, 8.4, 0.72], bt: [-2.1, 0.85, 4.1, 4.4], r: 0.23 },
    narrow: { fw: 16, fh: 9, reg: [-3.3, 0.95, 8.6, 3.8], tag: [4.45, 1.75, 3.1, 1.0], g: [1.6, 7.5, -1.05], bt: [-3.15, 0.9, 3.6, 3.85], r: 0.4 },
  };
  let L = LAY.wide, dead = false;
  const put = (o, z = 0) => { o.position.z = z; scene.add(o); return o; };
  const lab = (s, o) => put(ctx.label(s, Object.assign({ size: 13, weight: 600 }, o)), 5);
  const at = (o, x, y) => { o.position.x = x; o.position.y = y; return o; };
  /* a tinted plate: light fill + solid outline of one token, faded together */
  const plate = (tok, z) => {
    const f = put(ctx.box(1, 1, { color: tok, radius: 0.18 }), z);
    const s = put(ctx.box(1, 1, { color: null, stroke: tok, strokeWidth: 1.5, radius: 0.18 }), z + 0.05);
    return {
      size(w, h, x, y) { for (const b of [f, s]) { b.setSize(w, h, 0.18); at(b, x, y); } },
      set(o, fill = 0.14) { f.setOpacity(fill * o); s.setOpacity(o); },
    };
  };

  /* the panel: TEDDY's embedding, later the measured-protein evidence */
  const reg = put(ctx.box(1, 1, { color: 'soft', radius: 0.22 }), 0);
  const regT = put(ctx.box(1, 1, { color: null, stroke: 'teddy', strokeWidth: 1.5, radius: 0.22 }), 0.1);
  const regM = put(ctx.box(1, 1, { color: null, stroke: 'measured', strokeWidth: 1.5, radius: 0.22 }), 0.1);
  const titleT = lab('TEDDY embedding', { size: 12, color: 'teddy', anchor: 'top-left' });
  const titleM = lab('measured protein', { size: 12, color: 'measured', anchor: 'top-left' });
  const cosine = lab('cosine 0.982', { size: 11.5, mono: true, weight: 500, color: 'muted', anchor: 'top' });
  const cells = [['m', 'myeloid', 'right', -1], ['t', 'T cell', 'left', 1]].map(([tok, name, anchor, side]) => ({
    side, name, start: [0, 0], near: [0, 0], apart: [0, 0], off: 0,
    dot: put(ctx.dot([0, 0], { px: 10, color: tok, opacity: 0.9, order: 2 }), side > 0 ? 1.1 : 1),
    ring: put(ctx.dot([0, 0], { px: 14, hollow: true, ring: 0.2, color: 'teddy', order: 3 }), 1.5),
    lab: lab(name, { anchor }),
  }));

  /* the readout's call on the pair */
  const arrow = put(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.8, head: 8 }), 1);
  const tagT = plate('teddy', 1), tagA = plate('accent', 1.2);
  const same = lab('same call', { color: 'ink' }), apartTxt = lab('told apart', { color: 'ink' });
  const subSoft = lab('head + soft rule', { size: 11.5, weight: 500, color: 'muted', anchor: 'top' });
  const subStrict = lab('strict rule', { size: 11.5, weight: 500, color: 'muted', anchor: 'top' });
  const strike = put(ctx.line([[0, 0], [1, 1]], { color: 'accent', width: 2.5, order: 12 }), 3);

  /* the gauge: share of pairs that falsely agree (a full track = all pairs) */
  const track = put(ctx.line([[0, 0], [1, 0]], { color: 'line', width: 8 }), 1);
  const fillT = put(ctx.line([[0, 0], [1, 0]], { color: 'teddy', width: 8 }), 1.1);
  const fillA = put(ctx.line([[0, 0], [1, 0]], { color: 'accent', width: 8 }), 1.2);
  const ghost = put(ctx.line([[0, 0], [0, 1]], { color: 'teddy', width: 2 }), 1.3);
  const gTitle = lab('falsely agree', { size: 12, color: 'muted', anchor: 'bottom-left' });
  const v85 = lab('85%', { size: 12.5, mono: true, color: 'teddy', anchor: 'top' });
  const v40 = lab('40%', { size: 12.5, mono: true, color: 'accent', anchor: 'top' });

  /* the four steps of the strip below */
  const chips = BEATS.map((name) => ({
    base: put(ctx.box(1, 1, { color: 'card', stroke: 'line', strokeWidth: 1.2, radius: 0.18 }), 0),
    lit: plate('accent', 0.2),
    lab: lab(name, {}),
  }));

  /* text widths in CSS px, so cell labels stay inside the panel at every stage width */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, px) => { if (!mc) return s.length * px * 0.58; mc.font = `600 ${px}px ${family}`; return mc.measureText(s).width; };

  function layout() {
    const ppu = ctx.ppu(), px = (p) => p / ppu;
    const [rx, ry, rw, rh] = L.reg, [tx, ty, tw, th] = L.tag, [g0, g1, gy] = L.g, [by, bh, bw, bdx] = L.bt;
    for (const b of [reg, regT, regM]) { b.setSize(rw, rh, 0.22); at(b, rx, ry); }
    at(titleT, rx - rw / 2 + 0.2, ry + rh / 2 - 0.15); at(titleM, titleT.position.x, titleT.position.y);
    const rp = clamp(L.r * ppu, 7, 14), rr = px(rp), yc = ry + 0.05 * rh;
    const lw = px(Math.max(...cells.map((c) => textW(c.name, 13 * ctx.textScale))));
    const ax = Math.max(1.4 * rr, Math.min(0.3 * rw, rw / 2 - 0.22 - lw - rr - px(8)));
    cells.forEach((c) => {
      c.dot.setPx(rp); c.ring.setPx(rp * 1.45);
      c.start = [rx + c.side * ax, yc + c.side * 0.24 * rh];
      c.near = [rx + c.side * rr * 0.55, yc];
      c.apart = [rx + c.side * ax, yc];
      c.off = rr + px(8);
    });
    at(cosine, rx, yc - rr * 1.45 - px(6));
    tagT.size(tw, th, tx, ty); tagA.size(tw, th, tx, ty);
    at(same, tx, ty); at(apartTxt, tx, ty);
    at(subSoft, tx, ty - th / 2 - px(6)); at(subStrict, tx, ty - th / 2 - px(6));
    arrow.set([rx + rw / 2 + px(6), ty, 0], [tx - tw / 2 - px(6), ty, 0]);
    strike.setPoints([[tx - tw / 2 + 0.14, ty - th / 2 + 0.1], [tx + tw / 2 - 0.14, ty + th / 2 - 0.1]]);
    for (const l of [track, fillT, fillA]) l.setPoints([[g0, gy], [g1, gy]]);
    const gx = (f) => lerp(g0, g1, f);
    ghost.setPoints([[gx(BEFORE), gy - px(9)], [gx(BEFORE), gy + px(9)]]);
    at(gTitle, g0, gy + px(9)); at(v85, gx(BEFORE), gy - px(11)); at(v40, gx(AFTER), gy - px(11));
    chips.forEach((c, i) => {
      const x = (i - 1.5) * bdx;
      c.base.setSize(bw, bh, 0.18); at(c.base, x, by);
      c.lit.size(bw, bh, x, by); at(c.lab, x, by);
    });
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 10.2,
    resize(w, h) {
      L = w / h > 2.3 ? LAY.wide : LAY.narrow;
      const fit = camera.userData.animFit, a = w / h;
      fit.width = L.fw; fit.height = L.fh;
      const hh = a > L.fw / L.fh ? L.fh / 2 : L.fw / 2 / a, hw = hh * a;
      camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
      camera.updateProjectionMatrix();
      layout();
    },
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const out = 1 - seg(u, 0.92, 0.99, ease.inOutSine), kIn = seg(u, 0, 0.05, ease.inOutSine);
      const near = seg(u, 0.04, 0.24), kCall = seg(u, 0.26, 0.32), kV = seg(u, 0.48, 0.52);
      const comp = seg(u, 0.57, 0.63), apart = seg(u, 0.65, 0.77), tint = seg(u, 0.655, 0.675), fall = seg(u, 0.67, 0.78);
      /* texts that share a spot fade out, then in, so they never overlap */
      const tOut = seg(u, 0.57, 0.6), tIn = seg(u, 0.6, 0.63), sOut = seg(u, 0.66, 0.69), sIn = seg(u, 0.695, 0.735);
      const pe = ctx.pulse(u, 0.4, 0.48, 0.03);
      /* cells: together in the embedding, apart once measured protein is the evidence */
      cells.forEach((c) => {
        const x = lerp(lerp(c.start[0], c.near[0], near), c.apart[0], apart);
        const y = lerp(lerp(c.start[1], c.near[1], near), c.apart[1], apart);
        at(c.dot, x, y).setOpacity(0.9 * kIn * out);
        at(c.ring, x, y).setOpacity(kCall * (1 - 0.5 * kV) * (1 - comp) * out * (1 + 0.4 * pe) / 1.4);
        at(c.lab, x + c.side * c.off, y).setOpacity(kIn * out);
      });
      const m = comp * out;
      regT.setOpacity(1 - m); regM.setOpacity(m);
      titleT.setOpacity(clamp(1 - tOut + seg(u, 0.955, 0.99))); titleM.setOpacity(tIn * (1 - seg(u, 0.92, 0.955)));
      cosine.setOpacity(seg(u, 0.2, 0.26) * (1 - comp) * out);
      /* the call: same call (struck by the veto), then told apart */
      arrow.setProgress(seg(u, 0.24, 0.3)); arrow.setOpacity(seg(u, 0.24, 0.26) * out);
      const oT = kCall * (1 - sOut) * out * (1 - 0.45 * kV);
      tagT.set(oT, 0.14 + 0.2 * pe); same.setOpacity(oT);
      tagA.set(sIn * out); apartTxt.setOpacity(sIn * out);
      subSoft.setOpacity(kCall * (1 - tOut) * out); subStrict.setOpacity(kCall * tIn * out);
      strike.setProgress(seg(u, 0.49, 0.54)); strike.setOpacity(1 - seg(u, 0.655, 0.685));
      /* the gauge: fills to 85%, turns teal, falls to 40%; the 85% mark stays as the before value */
      const f = lerp(BEFORE * seg(u, 0.3, 0.38), AFTER, fall);
      fillT.setProgress(f); fillT.setOpacity((1 - tint) * out);
      fillA.setProgress(f); fillA.setOpacity(tint * out);
      ghost.setOpacity(0.85 * seg(u, 0.67, 0.7) * out);
      v85.setOpacity(seg(u, 0.35, 0.39) * out * (1 - 0.3 * fall));
      v40.setOpacity(seg(u, 0.74, 0.78) * out);
      /* the four steps light in order and stay lit until the loop restarts */
      chips.forEach((c, i) => {
        const k = seg(u, BEAT_T[i], BEAT_T[i] + 0.04, ease.inOutSine) * out;
        c.lit.set(k); c.lab.setOpacity(0.5 + 0.5 * k);
      });
    },
    dispose() {
      dead = true;
      cells.length = 0;
      chips.length = 0;
    },
  };
}
