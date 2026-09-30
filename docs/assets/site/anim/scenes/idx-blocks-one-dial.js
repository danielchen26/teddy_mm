/*
 * idx-blocks-one-dial — index.html, #blocks ("03 — Six experiments"), the strip above the table.
 * ANM's engine sits in the middle, wired to six dials, one per experiment and table row: 1 question,
 * 2 input, 3 leave-one-out, 4 cut-off (on ANM's score), 5 labels (the trained heads' label budget;
 * ANM uses 0), 6 evidence. Everything else stays fixed (section guide), so every dial starts locked.
 * Each experiment in turn unlocks exactly one dial and turns it; the change runs down its wire into
 * the same engine, then the dial turns back and locks again. The other five never move.
 * Dials 2 and 6 are drawn dashed and grey with a "withdrawn" tag, as the page marks them.
 * The dials carry no scale and no values. Wide stages: three dials each side in a row, wired over
 * the top. Narrow stages: three stacked on each side, wired straight in. 12 s seamless loop.
 * World units are CSS px (the camera is refitted to the stage in resize).
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 640, height: 183 });
  const PERIOD = 12, TURN = -2.0;                                 // turn = 115 degrees clockwise
  const DIALS = [['question'], ['input', true], ['leave-one-out'], ['cut-off'], ['labels'], ['evidence', true]];

  const circ = (r, z, n = 56) => Array.from({ length: n }, (_, i) => {
    const a = (i / n) * Math.PI * 2;
    return [r * Math.cos(a), r * Math.sin(a), z];
  });
  const cubic = (p0, p1, p2, p3, z, n = 40) => Array.from({ length: n + 1 }, (_, i) => {
    const t = i / n, a = (1 - t) ** 3, b = 3 * (1 - t) ** 2 * t, c = 3 * (1 - t) * t * t, d = t ** 3;
    return [a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0], a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1], z];
  });
  let mc = null;
  const textW = (s, size) => {
    const px = size * ctx.textScale;
    if (mc === null) mc = document.createElement('canvas').getContext('2d') || false;
    if (!mc) return s.length * px * 0.58;
    mc.font = `500 ${px}px "IBM Plex Sans", system-ui, sans-serif`;
    return mc.measureText(s).width + 2;
  };

  // The engine: the same ANM decision layer in every experiment.
  const engine = ctx.box(120, 60, { color: 'accent', opacity: 0.1, stroke: 'accent', strokeWidth: 1.5, radius: 10 });
  engine.position.z = -0.5;
  const engL = ctx.label('ANM', { size: 14, weight: 600 });
  const engS = ctx.label('engine', { size: 11, color: 'muted' });
  const expL = ctx.label('Experiment 1', { size: 11.5, color: 'muted', anchor: 'top' });
  scene.add(engine, engL, engS, expL);

  const dials = DIALS.map(([name, wd], i) => {
    const col = wd ? 'muted' : 'accent', dash = wd ? [4, 3] : false;
    const g = ctx.group();
    const disc = ctx.dot([0, 0, 0], { r: 20, color: 'soft' });
    const rim = ctx.line(circ(20, 0.2), { color: 'faint', width: 1.5, closed: true, dashed: dash });
    const rimOn = ctx.line(circ(20, 0.3), { color: col, width: 2.2, closed: true, dashed: dash, opacity: 0 });
    const ng = ctx.group();                                        // the pointer turns with this group
    const notch = ctx.line([[0, 10, 0.4], [0, 18, 0.4]], { color: 'muted', width: 2.5 });
    const notchOn = ctx.line([[0, 10, 0.5], [0, 18, 0.5]], { color: col, width: 2.5, opacity: 0 });
    ng.add(notch, notchOn);
    const num = ctx.label(String(i + 1), { size: 13, weight: 600 });
    g.add(disc, rim, rimOn, ng, num);
    const nameL = ctx.label(name, { size: 12 });
    const tag = wd ? ctx.label('withdrawn', { size: 10.5, color: 'muted', bg: 'soft', bgOpacity: 1, pad: 3 }) : null;
    // Padlock: a body and a shackle hinged on one leg; unlocked, it lifts and swings toward the dial.
    const lock = ctx.group();
    const body = ctx.box(10, 8, { color: 'muted', radius: 1.6 });
    body.position.z = 0.2;
    const sg = ctx.group();
    sg.position.set(3.2, 0, 0.1);
    const sp = [[0, 1.5]];
    for (let j = 0; j <= 12; j++) { const a = (Math.PI * j) / 12; sp.push([-3.2 + 3.2 * Math.cos(a), 5.5 + 3.2 * Math.sin(a)]); }
    sp.push([-6.4, 1.5]);
    const sh = ctx.line(sp, { color: 'muted', width: 1.6 });
    sg.add(sh);
    lock.add(body, sg);
    const wire = ctx.line([[0, 0], [1, 0]], { color: 'line', width: 1.5, dashed: wd ? [4, 4] : false });
    const wireOn = ctx.line([[0, 0], [1, 0]], { color: col, width: 2.2, dashed: wd ? [4, 4] : false, opacity: 0 });
    scene.add(g, nameL, lock, wire, wireOn);
    if (tag) scene.add(tag);
    const side = i < 3 ? -1 : 1;                                   // left dials hinge right, right dials left
    sg.position.x = -side * 3.2;
    sh.setPoints(sp.map((q) => [-side * q[0], q[1]]));
    return { i, side, wd: !!wd, g, disc, rim, rimOn, ng, notch, notchOn, num, nameL, tag, lock, body, sg, wire, wireOn };
  });

  function layout(w, h) {
    w = w || ctx.width || 640; h = h || ctx.height || 183;
    camera.userData.animFit = { width: w, height: h, fit: 'contain' };
    camera.left = -w / 2; camera.right = w / 2; camera.top = h / 2; camera.bottom = -h / 2;
    camera.updateProjectionMatrix();
    const ts = ctx.textScale, lh = (s) => s * ts * 1.3;
    const labW = Math.max(...DIALS.map(([n]) => textW(n, 12)), textW('withdrawn', 10.5) + 8);
    const m = 10, gap = 12, EWw = clamp(w * 0.13, 110, 170);
    const s = (w / 2 - m - EWw / 2 - gap) / 3;                    // slot width per dial, wide layout
    const wide = s - 14 >= labW && h >= 150;
    let R, EW, EH, yE, cx = 0;
    const P = [];
    if (wide) {
      R = clamp(Math.min(h * 0.09, s * 0.2), 14, 28);
      EW = EWw; EH = 2 * R + 8;
      const below = 6 + lh(12) + 4 + lh(10.5) + 6, arcZ = clamp(h * 0.4, 50, 130);
      yE = -(arcZ + 2 * R + below + 3) / 2 + below + R;
      const y0 = yE + R + 3, yTop = yE + EH / 2 + 1, r = R + 10;
      dials.forEach((d) => {
        const side = d.i < 3 ? -1 : 1, k = d.i < 3 ? 2 - d.i : d.i - 3;   // k = 0 next to the engine
        const x = side * (EW / 2 + gap + s * (k + 0.5)), e = side * (EW / 2) * (0.78 - 0.28 * k);
        const H = y0 + (arcZ * (0.4 + 0.3 * k)) / 0.75;
        P.push({ x, y: yE, name: [x, yE - R - 6, 'top'], tag: [x, yE - R - 10 - lh(12), 'top'],
          lock: [x + side * 0.707 * r, yE + 0.707 * r - 2],          // padlock on the outer-top diagonal
          wire: cubic([x, y0], [x, H], [e, H], [e, yTop], -1) });
      });
      expL.setAnchor('top'); expL.position.set(0, yE - EH / 2 - 7, 3);
    } else {
      // Three per side, names outside with the padlock inline; top-right corner kept for the pause control.
      const lw = (a, b, c) => Math.max(textW(DIALS[a][0], 12), textW(DIALS[b][0], 12), textW(DIALS[c][0], 12), textW('withdrawn', 10.5) + 8);
      const labL = lw(0, 1, 2), labR = lw(3, 4, 5);
      let mm = 6, lg = 23;
      EW = clamp(w * 0.2, 50, 120);
      cx = (labL - labR) / 2;
      const H0 = () => (w - 2 * mm - labL - labR) / 2 - lg;
      const yBot = -h / 2 + lh(11.5) + 6 + 14;
      const reserved = cx + H0() + lg + labR > w / 2 - 76;
      let yTop = reserved ? h / 2 - 44 : h / 2 - 30;
      R = clamp(((yTop - yBot) / 2) * 0.3, 9, 22);
      yTop = Math.min(yTop, h / 2 - 8 - R);
      const pitch = (yTop - yBot) / 2;
      let H = Math.min(H0(), EW / 2 + 2 * R + 110);
      for (let n = 0; n < 12 && H - EW / 2 - 2 * R < 8; n++) {     // squeeze: engine, then gaps, then dials
        if (EW > 44) EW -= 3; else if (lg > 19) lg -= 1; else if (mm > 4) mm -= 1; else R = Math.max(8, R - 1);
        H = Math.min(H0(), EW / 2 + 2 * R + 110);
      }
      yE = (yTop + yBot) / 2;
      EH = Math.max(pitch * 0.9 + 22, 2 * R + 20);
      dials.forEach((d) => {
        const side = d.i < 3 ? -1 : 1, y = yTop - (d.i % 3) * pitch, x = cx + side * (H - R);
        const x0 = x - side * (R + 3), xe = cx + side * (EW / 2 + 1), ye = yE + (y - yE) * 0.45, dx = (x0 - xe) * 0.5;
        const ny = y + (d.wd ? lh(12) / 2 + 1 : 0), nx = x + side * (R + lg), an = side < 0 ? 'right' : 'left';
        P.push({ x, y, name: [nx, ny, an], tag: [nx, y - lh(10.5) / 2 - 3, an], lock: [x + side * (R + 9), ny - 3],
          wire: cubic([x0, y], [x0 - dx, y], [xe + dx, ye], [xe, ye], -1) });
      });
      expL.setAnchor('bottom'); expL.position.set(cx, -h / 2 + 5, 3);
    }
    engine.setSize(EW, EH, 10);
    engine.position.set(cx, yE, -0.5);
    engL.position.set(cx, yE + 8 * ts, 3);
    engS.position.set(cx, yE - 10 * ts, 3);
    dials.forEach((d, j) => {
      const p = P[j];
      d.g.position.set(p.x, p.y, 0);
      d.disc.setRadius(R);
      d.rim.setPoints(circ(R, 0.2));
      d.rimOn.setPoints(circ(R, 0.3));
      const n0 = Math.max(R * 0.52, 6), n1 = R - 1.5;
      d.notch.setPoints([[0, n0, 0.4], [0, n1, 0.4]]);
      d.notchOn.setPoints([[0, n0, 0.5], [0, n1, 0.5]]);
      d.num.position.set(0, 0, 3);
      d.nameL.position.set(p.name[0], p.name[1], 3); d.nameL.setAnchor(p.name[2]);
      if (d.tag) { d.tag.position.set(p.tag[0], p.tag[1], 3); d.tag.setAnchor(p.tag[2]); }
      d.lock.position.set(p.lock[0], p.lock[1], 1);
      d.wire.setPoints(p.wire);
      d.wireOn.setPoints(p.wire.map((q) => [q[0], q[1], -0.8]));
    });
  }
  layout();

  let cur = -1;
  return {
    scene, camera, period: PERIOD, still: 0.56 * (PERIOD / 6),
    resize(w, h) { layout(w, h); },
    update(t) {
      const u = ctx.loopT(t, PERIOD), a = Math.min(5, Math.floor(u * 6)), v = u * 6 - a;
      if (a !== cur) { cur = a; expL.setText('Experiment ' + (a + 1)); }
      const k = seg(v, 0.02, 0.14) - seg(v, 0.86, 0.98);                               // unlock, relock
      const turn = seg(v, 0.14, 0.44, ease.inOutSine) - seg(v, 0.64, 0.86, ease.inOutSine);
      const wp = seg(v, 0.2, 0.5, ease.inOutSine);                                     // change runs in
      const glow = pulse(v, 0.46, 0.74, 0.08);                                         // engine re-runs
      dials.forEach((d) => {
        const on = d.i === a, kk = on ? k : 0;
        d.rimOn.setOpacity(kk); d.notchOn.setOpacity(kk); d.notch.setOpacity(1 - kk);
        d.ng.rotation.z = on ? TURN * turn : 0;
        d.num.setOpacity(lerp(0.55, 1, kk)); d.nameL.setOpacity(lerp(0.72, 1, kk));
        d.sg.position.y = 3.2 * kk; d.sg.rotation.z = d.side * 0.55 * kk;
        d.body.setOpacity(lerp(1, 0.55, kk));
        d.wire.setOpacity(lerp(1, 0.35, kk));
        d.wireOn.setProgress(on ? wp : 0); d.wireOn.setOpacity(kk);
      });
      engine.material.opacity = 0.08 + (dials[a].wd ? 0.06 : 0.14) * glow;
      expL.setOpacity(pulse(v, 0, 1, 0.05));
    },
    dispose() { mc = null; },
  };
}
