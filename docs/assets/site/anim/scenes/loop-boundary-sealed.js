/*
 * loop-boundary-sealed (anm-loop.html #boundary, "Claim boundary").
 * TEDDY and our head sit in a sealed, padlocked capsule: TEDDY is frozen and the head's checkpoint (best.pt) is
 * not retrained; the only thing leaving it is the head's predicted protein, which flows into the evidence.
 * The decision loop (a teal ring around ANM) runs outside the capsule and touches two things only: the evidence
 * (measured protein swapped in for predicted) and the question (soft rule → strict rule). Measured protein
 * normally flows to the answer key; it reaches the evidence only through one gate, "this test only".
 * A tag on the ring: withdrawn · pending redesign. Illustrative only: no values. 12 s seamless loop.
 * Two layouts: the page's 21:6 strip, and 16:9 on narrow screens (the page's --ar-sm).
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, clamp, lerp } = ctx;
  const C = ctx.colors;
  const PI = Math.PI, PERIOD = 12;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 6 });
  const LAY = {
    /* wide: measured protein → answer key on the right; the gated path drops from it and runs left into the evidence */
    wide: { fw: 21, fh: 6, capL: -10.3, capR: -3.9, capY: -0.3, ring: [0.3, 0.1, 2.6, 1.5], meas: [5.2, 0.1],
      branch: false, eLeg: 'bottom', measText: 'measured\nprotein', tagA: 0.62, tag: [3.0, 1.45], tagSize: 10.5,
      tagText: 'withdrawn · pending redesign', loop: [-2.05, 0.95, 'bottom-right'] },
    /* narrow: measured protein → answer key along the bottom; the gated path branches up into the evidence */
    narrow: { fw: 16, fh: 9, capL: -7.7, capR: -0.6, capY: 1.4, ring: [2.65, 0.9, 2.1, 1.45], meas: [null, -3.5],
      branch: true, eLeg: 'right', measText: 'measured protein', tagA: 0, tag: null, tagSize: 10,
      tagText: 'withdrawn\npending\nredesign', loop: null },
  };
  let L = LAY.wide, dead = false;
  const R = { x: 0, y: 0, rx: 1, ry: 1 };
  let pathA = [[0, 0], [1, 0]], pathB = [[0, 0], [1, 0]];
  const put = (o, x = 0, y = 0, z = 0, p = scene) => { o.position.set(x, y, z); p.add(o); return o; };
  const text = (s, o) => put(ctx.label(s, o), 0, 0, 6);
  const panel = (fill, stroke, z, sw = 1.5) => put(ctx.box(1, 1, { color: fill, stroke, strokeWidth: sw }), 0, 0, z);
  const tint = (color, a, z) => put(ctx.box(1, 1, { color, opacity: a }), 0, 0, z);

  /* the sealed capsule: TEDDY → our head, padlocked; its only output is the head's predicted protein */
  const cap = panel('soft', 'muted', 0, 2), capIn = panel(null, 'line', 0.05, 1.2);
  const tBox = panel('card', 'teddy', 0.1), tTint = tint('teddy', 0.14, 0.12), hBox = panel('card', 'muted', 0.1);
  const tLab = text('TEDDY', { size: 12.5, weight: 600 }), hLab = text('our head', { size: 12.5, weight: 600 });
  const tSub = text('frozen', { size: 10.5, color: 'muted', anchor: 'top' });
  const hSub = text('best.pt', { size: 10.5, color: 'muted', anchor: 'top', mono: true });
  const inArrow = put(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.5, head: 7 }), 0, 0, 0.2);
  const lock = put(ctx.group(), 0, 0, 0.3);
  put(ctx.box(0.34, 0.26, { color: 'muted', radius: 0.05 }), 0, 0, 0, lock);
  const sh = [[0.09, 0.12], [0.09, 0.2]];
  for (let k = 0; k <= 12; k++) { const a = (PI * k) / 12; sh.push([0.09 * Math.cos(a), 0.2 + 0.09 * Math.sin(a)]); }
  sh.push([-0.09, 0.12]);
  put(ctx.line(sh, { color: 'muted', width: 2 }), 0, 0, 0, lock);
  const sealLab = text('sealed', { size: 11, weight: 600, color: 'muted', anchor: 'left' });
  const outArrow = put(ctx.arrow([0, 0], [1, 0], { color: 'teddy', width: 1.8, head: 8, bend: -0.12 }), 0, 0, 0.4);
  const outTok = [0, 1, 2].map(() => put(ctx.dot([0, 0], { px: 3.5, color: 'teddy' }), 0, 0, 3));

  /* the decision loop: a teal ring around ANM, touching only the evidence and the question */
  const ring = put(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'accent', width: 2, dashed: [7, 6], closed: true }), 0, 0, 0.5);
  const tok = put(ctx.dot([0, 0], { px: 5.5, color: 'accent' }), 0, 0, 2);
  const loopLab = text('decision\nloop', { size: 11.5, weight: 600, color: 'accent' });
  const anm = panel('card', 'accent', 1), anmTint = tint('accent', 0.14, 1.05);
  const anmLab = text('ANM', { size: 13, weight: 700 });
  const toAnm = [0, 1].map(() => put(ctx.arrow([0, 0], [0, 1], { color: 'faint', width: 1.5, head: 7 }), 0, 0, 0.9));
  const node = (states, colors, legend) => ({
    base: panel('card', 'accent', 1), tint: tint('accent', 0.07, 1.05),
    hi: put(ctx.box(1, 1, { color: null, stroke: 'accent', strokeWidth: 3, opacity: 0 }), 0, 0, 1.1),
    labs: states.map((s, i) => text(s, { size: 12.5, weight: 600, color: colors[i] })),
    leg: text(legend, { size: 10.5, weight: 600, color: 'muted', bg: 'card', bgOpacity: 1, pad: 3, order: 11 }),
    states, legend, x: 0, y: 0, w: 1, h: 1,
  });
  const card = node(['soft rule', 'strict rule'], ['ink', 'ink'], 'question');
  const plug = node(['predicted', 'measured'], ['teddy', 'ink'], 'evidence');

  /* measured protein: to the answer key, or through the one gate into the evidence */
  const meas = panel('card', 'measured', 1), key = panel('soft', 'measured', 1);
  const measLab = text(L.measText, { size: 12, weight: 600 }), keyLab = text('answer key', { size: 12, weight: 600 });
  const toKey = put(ctx.arrow([0, 0], [1, 0], { color: 'measured', width: 1.8, head: 8 }), 0, 0, 0.9);
  const aTok = [0, 1].map(() => put(ctx.dot([0, 0], { px: 3.5, color: 'measured' }), 0, 0, 3));
  const bPath = put(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 1.5, dashed: [4, 4] }), 0, 0, 0.6);
  const bOn = put(ctx.line([[0, 0], [1, 0]], { color: 'measured', width: 1.8, opacity: 0 }), 0, 0, 0.62);
  const bTok = [0, 1, 2].map(() => put(ctx.dot([0, 0], { px: 3.5, color: 'measured', opacity: 0 }), 0, 0, 3));
  /* the gate: a bar across the path whose two halves slide apart while it is open */
  const gate = put(ctx.group(), 0, 0, 2.5);
  const halves = [1, -1].map(() => put(ctx.line([[0, 0], [0, 1]], { color: 'muted', width: 3.5 }), 0, 0, 0, gate));
  const gateLab = text('this test only', { size: 11, weight: 600, anchor: 'top' });
  let gb = 0.3;

  /* the tag hanging on the ring */
  const str = put(ctx.line([[0, 0], [1, 1]], { color: 'faint', width: 1.2 }), 0, 0, 2.1);
  const knot = put(ctx.dot([0, 0], { px: 2.5, color: 'faint' }), 0, 0, 2.1);
  const tagBox = put(ctx.box(1, 1, { color: 'card', stroke: 'bad', radius: 0.1 }), 0, 0, 2.2);
  const tagLab = text(L.tagText, { size: 10.5, weight: 600, color: 'bad' });

  /* text widths in CSS px, so boxes fit their labels and nothing collides at small stage widths */
  const mc = document.createElement('canvas').getContext('2d');
  const cs = () => getComputedStyle(ctx.figure);
  const tw = (s, size, wt = 600, mono = false) => {
    const px = size * ctx.textScale, ls = s.split('\n');
    if (!mc) return Math.max(...ls.map((l) => l.length)) * px * (mono ? 0.62 : 0.56);
    mc.font = `${wt} ${px}px ${mono ? cs().getPropertyValue('--mono').trim() || 'monospace' : cs().fontFamily || 'sans-serif'}`;
    return Math.max(...ls.map((l) => mc.measureText(l).width));
  };
  const th = (s, size) => s.split('\n').length * size * ctx.textScale * 1.22;
  const fit = (b, x, y, w, h, r = 0.16) => { b.position.set(x, y, b.position.z); b.setSize(w, h, r); };

  function layout() {
    const p = ctx.ppu(), ts = ctx.textScale, W = (px) => px / p;
    [R.x, R.y, R.rx, R.ry] = L.ring;
    /* capsule: two boxes side by side; labels shrink a little if the capsule would not fit */
    const m = Math.max(0.22, W(7)), gap = Math.max(0.6, W(22)), avail = L.capR - L.capL;
    const T = Math.max(tw('TEDDY', 12.5), tw('our head', 12.5));
    const k = clamp((((avail - gap - 2 * m) * p) / 2 - 14) / T, 0.78, 1);
    tLab.setSize(12.5 * k); hLab.setSize(12.5 * k); tSub.setSize(10.5 * k); hSub.setSize(10.5 * k);
    const bw = W(T * k + 14), bh = Math.max(0.85, W(12.5 * ts * k * 1.25 + 14)), subH = W(10.5 * ts * k * 1.25 + 4);
    const cw = Math.max(2 * bw + gap + 2 * m, Math.min(avail, 5.6)), ch = bh + subH + 2 * m + 0.1;
    const cx = L.capR - cw / 2, cy = L.capY, by = cy + ch / 2 - m - bh / 2, cr = Math.min(ch / 2, 0.6);
    fit(cap, cx, cy, cw, ch, cr); fit(capIn, cx, cy, cw - 0.16, ch - 0.16, cr - 0.08);
    const tx = cx - cw / 2 + m + bw / 2, hx = cx + cw / 2 - m - bw / 2;
    fit(tBox, tx, by, bw, bh); fit(tTint, tx, by, bw, bh); fit(hBox, hx, by, bw, bh);
    tLab.position.set(tx, by, 6); hLab.position.set(hx, by, 6);
    tSub.position.set(tx, by - bh / 2 - 0.05, 6); hSub.position.set(hx, by - bh / 2 - 0.05, 6);
    inArrow.set([tx + bw / 2 + 0.06, by], [hx - bw / 2 - 0.06, by]);
    const lk = clamp(11 / (0.26 * p), 1, 2.2), ly = cy - ch / 2 - 0.2 * lk;
    lock.scale.setScalar(lk); lock.position.set(cx, ly, 0.3); sealLab.position.set(cx + 0.26 * lk, ly, 6);
    /* ring, ANM, and the two things the ring touches (tangent at top and bottom) */
    const pts = [];
    for (let i = 0; i < 120; i++) { const a = (i / 120) * 2 * PI; pts.push([R.x + R.rx * Math.cos(a), R.y + R.ry * Math.sin(a)]); }
    ring.setPoints(pts);
    const aw = W(tw('ANM', 13, 700) + 26), ah = Math.max(0.8, W(13 * ts * 1.25 + 14));
    fit(anm, R.x, R.y, aw, ah, 0.18); fit(anmTint, R.x, R.y, aw, ah, 0.18); anmLab.position.set(R.x, R.y, 6);
    loopLab.visible = !!L.loop;
    if (L.loop) { loopLab.position.set(L.loop[0], L.loop[1], 6); loopLab.setAnchor(L.loop[2]); }
    for (const [n, side] of [[card, 1], [plug, -1]]) {
      n.w = W(Math.max(...n.states.map((s) => tw(s, 12.5)), tw(n.legend, 10.5) + 12) + 24);
      n.h = Math.max(0.9, W(12.5 * ts * 1.25 + 14));
      n.x = R.x; n.y = R.y + side * (R.ry + n.h / 2);
      [n.base, n.tint, n.hi].forEach((b) => fit(b, n.x, n.y, n.w, n.h, 0.18));
      n.labs.forEach((l) => l.position.set(n.x, n.y, 6));
    }
    card.leg.position.set(card.x, card.y + card.h / 2, 7);
    if (L.eLeg === 'bottom') { plug.leg.position.set(plug.x, plug.y - plug.h / 2, 7); plug.leg.setAnchor('center'); }
    else { plug.leg.position.set(plug.x + plug.w / 2 + 0.12, plug.y, 7); plug.leg.setAnchor('left'); }
    toAnm[0].set([R.x, plug.y + plug.h / 2 + 0.05], [R.x, R.y - ah / 2 - 0.06]);
    toAnm[1].set([R.x, card.y - card.h / 2 - 0.05], [R.x, R.y + ah / 2 + 0.06]);
    outArrow.set([L.capR + 0.06, by], [plug.x - plug.w / 2 - 0.08, plug.y]);
    /* measured protein, the answer key, the gated path into the evidence */
    measLab.setText(L.measText);
    const mw = W(tw(L.measText, 12) + 22), mh = Math.max(0.9, W(th(L.measText, 12) + 14));
    const kw = W(tw('answer key', 12) + 22), kh = Math.max(0.9, W(th('answer key', 12) + 14));
    const mx = L.meas[0] == null ? L.capL + mw / 2 : L.meas[0], my = L.meas[1], kx = L.fw / 2 - 0.12 - kw / 2;
    fit(meas, mx, my, mw, mh); measLab.position.set(mx, my, 6);
    fit(key, kx, my, kw, kh); keyLab.position.set(kx, my, 6);
    pathA = [[mx + mw / 2 + 0.06, my], [kx - kw / 2 - 0.08, my]];
    toKey.set(pathA[0], pathA[1]);
    aTok[1].visible = (pathA[1][0] - pathA[0][0]) * p > 80;
    const pb = plug.y - plug.h / 2;
    pathB = L.branch ? [pathA[0], [plug.x, my], [plug.x, pb]] : [[mx, my - mh / 2], [mx, plug.y], [plug.x + plug.w / 2, plug.y]];
    const drawn = L.branch ? pathB.slice(1) : pathB;
    bPath.setPoints(drawn); bOn.setPoints(drawn);
    gb = Math.max(0.22, W(9));
    const gx = L.branch ? plug.x : (plug.x + plug.w / 2 + mx) / 2, gy = L.branch ? (my + pb) / 2 : plug.y;
    gate.position.set(gx, gy, 2.5); gate.rotation.z = L.branch ? PI / 2 : 0;
    halves.forEach((h, i) => h.setPoints([[0, 0], [0, (i ? -1 : 1) * gb]]));
    if (L.branch) { gateLab.position.set(gx - 1.7 * gb - 0.1, gy, 6); gateLab.setAnchor('right'); }
    else { gateLab.position.set(gx, gy - 1.7 * gb - 0.06, 6); gateLab.setAnchor('top'); }
    /* the tag, tied to the ring */
    tagLab.setText(L.tagText); tagLab.setSize(L.tagSize);
    const tgw = W(tw(L.tagText, L.tagSize) + 14), tgh = W(th(L.tagText, L.tagSize) + 10);
    const ax = R.x + R.rx * Math.cos(L.tagA), ay = R.y + R.ry * Math.sin(L.tagA);
    const [qx, qy] = L.tag || [ax + 0.3, ay], tcx = qx + tgw / 2;
    fit(tagBox, tcx, qy, tgw, tgh, 0.1); tagLab.position.set(tcx, qy, 6);
    str.setPoints([[ax, ay], [qx, qy]]); knot.position.set(ax, ay, 2.1);
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  /* ring token: left → evidence (dwell) → question (dwell) → left; one lap per loop */
  const angle = (u) => (u < 0.14 ? lerp(PI, 1.5 * PI, ease.inOutSine(u / 0.14))
    : u < 0.4 ? 1.5 * PI
    : u < 0.58 ? lerp(1.5 * PI, 2.5 * PI, ease.inOutSine((u - 0.4) / 0.18))
    : u < 0.8 ? 2.5 * PI
    : lerp(2.5 * PI, 3 * PI, ease.inOutSine((u - 0.8) / 0.14)));

  return {
    scene, camera, period: PERIOD, still: 8.64,
    resize(w, h) {
      L = w / h > 2.4 ? LAY.wide : LAY.narrow;
      const f = camera.userData.animFit, a = w / h;
      f.width = L.fw; f.height = L.fh;
      const hh = a > L.fw / L.fh ? L.fh / 2 : L.fw / 2 / a, hw = hh * a;
      camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
      camera.updateProjectionMatrix();
      layout();
    },
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const th0 = angle(u);
      tok.position.set(R.x + R.rx * Math.cos(th0), R.y + R.ry * Math.sin(th0), 2);
      ring.setDashPhase(u, 6);
      plug.hi.setOpacity(pulse(u, 0.12, 0.42, 0.04));
      card.hi.setOpacity(pulse(u, 0.56, 0.82, 0.04));
      /* the gate opens when the loop reaches the evidence, for this test, and shuts before the loop restarts */
      const open = seg(u, 0.15, 0.21) * (1 - seg(u, 0.86, 0.92));
      halves.forEach((h, i) => {
        h.position.y = (i ? -1 : 1) * 0.7 * gb * open;
        h.material.uniforms.uColor.value.copy(C.muted).lerp(C.accent, open);
      });
      bOn.setOpacity(0.8 * open);
      gateLab.setOpacity(0.72 + 0.28 * open);
      bTok.forEach((d, i) => {
        const a = 0.17 + i * 0.06;
        d.position.copy(ctx.along(pathB, seg(u, a, a + 0.14, ease.inOutSine))); d.position.z = 3;
        d.setOpacity(pulse(u, a, a + 0.14, 0.015));
      });
      /* what the loop changes: the evidence, then the question; both return before the next lap */
      const sw = seg(u, 0.29, 0.35) * (1 - seg(u, 0.88, 0.94));
      plug.labs[0].setOpacity(1 - sw); plug.labs[1].setOpacity(sw);
      const cq = seg(u, 0.62, 0.68) * (1 - seg(u, 0.88, 0.94));
      card.labs[0].setOpacity(1 - cq); card.labs[1].setOpacity(cq);
      anmTint.setOpacity(0.14 + 0.14 * pulse(u, 0.7, 0.84, 0.04));
      /* steady streams: predicted protein out of the capsule, measured protein to the answer key */
      outTok.forEach((d, i) => {
        const s = (u * 4 + i / 3) % 1;
        d.position.copy(ctx.along(outArrow.line, s * 0.92)); d.position.z = 3;
        d.setOpacity(clamp(s / 0.08) * clamp((1 - s) / 0.1));
      });
      aTok.forEach((d, i) => {
        const s = (u * 4 + i / 2) % 1;
        d.position.copy(ctx.along(pathA, s * 0.9)); d.position.z = 3;
        d.setOpacity(clamp(s / 0.1) * clamp((1 - s) / 0.12));
      });
    },
    dispose() {
      dead = true;
      outTok.length = aTok.length = bTok.length = 0;
    },
  };
}
