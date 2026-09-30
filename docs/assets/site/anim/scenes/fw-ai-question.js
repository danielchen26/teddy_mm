/*
 * fw-ai-question (anm-framework.html #ai-question, "AI step 1: write the question down").
 * A question card types itself: what is read, which answers are allowed, the bar. Its answers fly out
 * to become the slots of a dial, the bar becomes a dashed threshold ring, and the band past the ring is
 * 'still works' (an allowed answer cleared the bar). Evidence runs in from the left and grows each
 * answer's score; the needle follows the best one. The answer key waits in a locked glass case; a
 * particle from the key's side bounces off the glass (no leakage). Only after the needle stops does the
 * case open: the key's slot past the bar is 'exactly right', and the key stamps the result.
 * Run 1: B clears the bar and matches the key (tick). Run 2: no score reaches the bar, 'no call' (cross).
 * Scores are illustrative, not data. 12 s seamless loop.
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12, STAG = 0.012, FLY = 0.065;
  const G = { x: 0.3, y: -2.0 }, R = 3.4, RT = 2.1, R0 = 0.4, HW = Math.PI / 12;
  const ANG = [150, 90, 30].map((d) => (d * Math.PI) / 180);
  const NAMES = ['A', 'B', 'C'], KEY = 1;
  const CX = 5.7, CY = -1.0, CW = 2.4, CH = 1.8;
  const CL = CX - CW / 2, CR = CX + CW / 2, CT = CY + CH / 2, CB = CY - CH / 2;
  const SX = G.x + 1.55, SY = -3.5;
  const X0 = -7.6, Y0 = 4.05;
  const TEXT = ['Read: answer scores', 'Allowed: A · B · C', 'Bar: score ≥ θ'];
  const TYPE = [[0.015, 0.06], [0.065, 0.11], [0.115, 0.155]];

  const add = (o, x = 0, y = 0, z = 0) => { o.position.set(x, y, z); scene.add(o); return o; };
  const col = (o, c) => { if (o.userData.c !== c) { o.userData.c = c; o.setColor(c); } };
  const arc = (r, a0, a1, n = 48) => Array.from({ length: n + 1 }, (_, i) => {
    const a = a0 + ((a1 - a0) * i) / n;
    return [G.x + r * Math.cos(a), G.y + r * Math.sin(a)];
  });
  const flat = (geo, token, z) => {
    const m = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, side: THREE.DoubleSide }));
    ctx.bind(m, token);
    const o = add(new THREE.Mesh(ctx.track(geo), m), G.x, G.y, z);
    o.setOpacity = (a) => { m.opacity = a; };
    return o;
  };
  const lab = (s, o, x, y) => add(ctx.label(s, Object.assign({ size: 12, color: 'muted' }, o)), x, y, 1);

  /* two readouts: the order in which evidence reaches answers A/B/C, and how much each piece adds */
  const RUNS = [
    { T: 0.2, end: 0.55, order: [0, 1, 1, 0, 1, 2, 1, 1], inc: [0.42, 0.5, 0.4], call: 'call: B', ok: true },
    { T: 0.595, end: 0.945, order: [2, 0, 1, 0, 1, 2, 1], inc: [0.36, 0.33, 0.3], call: 'no call', ok: false },
  ].map((r) => {
    const n = [0, 0, 0];
    r.dots = r.order.map((a, i) => ({ a, s: r.T + i * STAG, rb: R0 + r.inc[a] * n[a]++ }));
    r.done = r.T + (r.order.length - 1) * STAG + FLY + 0.012; // needle has stopped
    r.lid = r.done + 0.022; // case opens only after that
    return r;
  });
  const inRun = (u, r, a, b, fn) => seg(u, a, b, fn) * (1 - seg(u, r.end, r.end + 0.04));

  /* the question card (sized from the measured text so it fits at every stage width) */
  const card = add(ctx.box(1, 1, { color: 'soft', stroke: 'line', radius: 0.14 }), 0, 0, -0.5);
  const lines = TEXT.map(() => add(ctx.label('', { size: 13, anchor: 'left' }), 0, 0, 1));
  const meas = document.createElement('canvas').getContext('2d');
  const letters = NAMES.map(() => [0, 0]);

  /* the dial */
  const rail = add(ctx.line([[-7.8, G.y], [G.x + R, G.y]], { color: 'line', width: 1.5 }));
  const rim = add(ctx.line(arc(R, Math.PI, 0), { color: 'faint', width: 1.5 }));
  const band = flat(new THREE.RingGeometry(RT, R, 64, 1, 0, Math.PI), 'accent', -0.2);
  const sectors = ANG.map((a) => flat(new THREE.CircleGeometry(1, 16, a - HW, 2 * HW), 'accent', 0.1));
  const ring = add(ctx.line(arc(RT, Math.PI, 0), { color: 'ink', width: 2, dashed: [6, 5] }), 0, 0, 0.3);
  const aK = ANG[KEY];
  const wedge = add(ctx.line([...arc(RT, aK - HW, aK + HW, 10), ...arc(R, aK + HW, aK - HW, 10)],
    { color: 'measured', width: 2, closed: true }), 0, 0, 0.4);
  const needle = add(ctx.line([[G.x, G.y], [G.x, G.y + R0]], { color: 'ink', width: 2.5 }), 0, 0, 0.6);
  const tip = add(ctx.dot([0, 0], { px: 3.5, color: 'ink' }), 0, 0, 0.7);
  const hub = add(ctx.dot([0, 0], { px: 4.5, color: 'ink' }), G.x, G.y, 0.7);
  const slots = NAMES.map((n) => add(ctx.label(n, { size: 13, weight: 600 }), 0, 0, 1));
  const evLab = lab('evidence', { anchor: 'bottom-left' }, -7.8, G.y + 0.12);
  const barLab = lab('bar', { anchor: 'top' }, G.x + RT, G.y - 0.16);
  const stillLab = lab('still works', { anchor: 'top', weight: 600 }, G.x - (RT + R) / 2, G.y - 0.16);
  const exactLab = lab('exactly right', { anchor: 'bottom', weight: 600, color: 'measured' }, G.x, 0);
  const chip = lab('call: B', { size: 13, weight: 600, color: 'ink', bg: 'soft', pad: 5, anchor: 'right' }, SX - 0.62, SY);
  const dots = Array.from({ length: 8 }, () => add(ctx.dot([0, 0], { px: 4, color: 'accent' }), 0, 0, 0.8));

  /* the answer key in a locked glass case */
  const caseFill = add(ctx.box(CW, CH, { color: 'soft', radius: 0.06 }), CX, CY, -0.3);
  const walls = add(ctx.line([[CR, CT], [CR, CB], [CL, CB]], { color: 'muted', width: 2 }), 0, 0, 0.5);
  const glass = add(ctx.line([[CL, CB], [CL, CT]], { color: 'muted', width: 2 }), 0, 0, 0.5);
  const glare = [[-1.0, 0.0, -0.7, 0.5], [-1.0, -0.45, -0.45, 0.45]].map(([a, b, c, d]) =>
    add(ctx.line([[CX + a, CY + b], [CX + c, CY + d]], { color: 'faint', width: 1.5 }), 0, 0, 0.45));
  const keyG = add(ctx.group(), CX + 0.3, CY - 0.15, 0.2);
  const keyBox = ctx.box(1.0, 0.64, { color: 'card', stroke: 'measured', strokeWidth: 1.5, radius: 0.1 });
  const keyLab = ctx.label('B', { size: 15, weight: 700, color: 'measured' });
  keyG.add(keyBox, keyLab);
  const lid = add(ctx.group(), CL, CT, 0.55);
  const lidLine = ctx.line([[0, 0], [CW, 0]], { color: 'muted', width: 2.5 });
  lid.add(lidLine);
  const lock = add(ctx.group(), CX, CT, 0.9);
  const lockBody = ctx.box(0.44, 0.34, { color: 'muted', radius: 0.06 });
  const shackle = ctx.line([[0.13, 0.12], ...Array.from({ length: 13 }, (_, i) =>
    [0.13 * Math.cos((Math.PI * i) / 12), 0.2 + 0.13 * Math.sin((Math.PI * i) / 12)]), [-0.13, 0.12]], { color: 'muted', width: 2 });
  lock.add(lockBody, shackle);
  const keyCase = lab('answer key', { anchor: 'top', color: 'ink', weight: 600 }, CX, CB - 0.15);
  const leakLab = lab('no leakage', { anchor: 'top' }, CX, CB - 0.8);
  const leak = add(ctx.dot([0, 0], { px: 4, color: 'measured' }), 0, CY - 0.15, 0.8);

  /* the key grades the result (a dashed arrow from the case, then a stamp) */
  const grade = add(ctx.arrow([CL - 0.15, CB + 0.1], [SX + 0.5, SY + 0.22], { color: 'measured', width: 1.5, head: 7, dashed: [4, 4] }), 0, 0, 1.5);
  const stamp = add(ctx.group(), SX, SY, 2);
  const sFill = ctx.dot([0, 0], { r: 0.4, color: 'card' });
  const sRing = ctx.dot([0, 0, 0.01], { r: 0.4, hollow: true, ring: 0.14, color: 'good' });
  const tick = ctx.line([[-0.19, 0.01], [-0.05, -0.14], [0.2, 0.16]], { color: 'good', width: 3 });
  const crossA = ctx.line([[-0.15, -0.15], [0.15, 0.15]], { color: 'bad', width: 3 });
  const crossB = ctx.line([[-0.15, 0.15], [0.15, -0.15]], { color: 'bad', width: 3 });
  stamp.add(sFill, sRing, tick, crossA, crossB);

  let lay = '', dead = false;
  function layout(force) {
    const ppu = ctx.ppu(), ts = ctx.textScale, k = `${ppu.toFixed(3)}/${ts.toFixed(3)}`;
    if (k === lay && !force) return;
    lay = k;
    meas.font = `500 ${13 * ts}px ${getComputedStyle(ctx.figure).fontFamily}`;
    const w = (s) => meas.measureText(s).width / ppu;
    const pad = (12 * ts) / ppu, gap = (22 * ts) / ppu;
    const W = Math.max(...TEXT.map(w)) + 2 * pad, H = 3 * gap + pad;
    card.setSize(W, H, 0.14);
    card.position.set(X0 + W / 2, Y0 - H / 2, -0.5);
    lines.forEach((l, i) => l.position.set(X0 + pad, Y0 - pad / 2 - gap * (i + 0.5), 1));
    NAMES.forEach((c, i) => {
      const j = TEXT[1].indexOf(c, 1);
      letters[i] = [X0 + pad + w(TEXT[1].slice(0, j)) + w(c) / 2, lines[1].position.y];
    });
    leakLab.position.y = CB - 0.15 - (19 * ts) / ppu;
  }
  if (document.fonts) document.fonts.ready.then(() => { if (!dead) { layout(true); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 0.52 * PERIOD,
    update(t) {
      layout();
      const u = ctx.loopT(t, PERIOD);
      const V = seg(u, 0, 0.02) * (1 - seg(u, 0.945, 0.99));
      const ppu = ctx.ppu(), ts = ctx.textScale;
      const run = u < 0.592 ? RUNS[0] : RUNS[1];

      /* 1. the question is written down first */
      card.setOpacity(V);
      lines.forEach((l, i) => {
        const [a, b] = TYPE[i];
        l.setText(TEXT[i].slice(0, Math.round(seg(u, a, b, ease.linear) * TEXT[i].length)));
        l.setOpacity(V);
        const hot = (i === 1 && pulse(u, 0.105, 0.2, 0.01) > 0.5) || (i === 2 && pulse(u, 0.15, 0.23, 0.01) > 0.5);
        col(l, hot ? 'accent' : 'ink');
      });

      /* 2. it becomes the dial: answers → slots, bar → ring, band past the ring */
      const g0 = seg(u, 0.03, 0.1);
      [rail, rim].forEach((o) => { o.setProgress(g0); o.setOpacity(V); });
      evLab.setOpacity(seg(u, 0.05, 0.09) * V);
      const slotR = R + 0.1 + (10 * ts) / ppu, fly = seg(u, 0.11, 0.17);
      slots.forEach((s, i) => {
        s.position.set(lerp(letters[i][0], G.x + slotR * Math.cos(ANG[i]), fly), lerp(letters[i][1], G.y + slotR * Math.sin(ANG[i]), fly), 1);
        s.setOpacity(seg(u, 0.105, 0.115) * V);
      });
      ring.setProgress(seg(u, 0.155, 0.205));
      ring.setOpacity(V);
      barLab.setOpacity(seg(u, 0.18, 0.21) * V);

      /* 3. evidence grows each answer's score; the needle follows the best one */
      const r = [R0, R0, R0];
      for (const q of RUNS) {
        const back = 1 - seg(u, q.end, q.end + 0.04);
        for (const d of q.dots) r[d.a] += q.inc[d.a] * seg(u, d.s + FLY, d.s + FLY + 0.012, ease.outCubic) * back;
      }
      const grow = seg(u, 0.13, 0.17) * V, rm = Math.max(...r);
      sectors.forEach((s, i) => { s.scale.set(r[i], r[i], 1); s.setOpacity(0.5 * grow); });
      let ws = 0, wa = 0;
      r.forEach((v, i) => { const w = Math.exp(14 * (v - rm)); ws += w; wa += w * ANG[i]; });
      const na = wa / ws, tx = G.x + rm * Math.cos(na), ty = G.y + rm * Math.sin(na);
      needle.setPoints([[G.x, G.y], [tx, ty]]);
      tip.position.set(tx, ty, 0.7);
      [needle, tip, hub].forEach((o) => o.setOpacity(grow));
      const lit = clamp((rm - RT) / 0.12);
      band.setOpacity((0.06 + 0.1 * lit) * seg(u, 0.18, 0.22) * V);
      stillLab.setOpacity(seg(u, 0.19, 0.23) * V);
      col(stillLab, lit > 0.5 ? 'accent' : 'muted');
      dots.forEach((d, i) => {
        const e = run.dots[i], k = e ? (u - e.s) / FLY : -1;
        if (k <= 0 || k >= 1.15) { d.setOpacity(0); return; }
        const p = ctx.along([[-8.6, G.y], [G.x, G.y], [G.x + e.rb * Math.cos(ANG[e.a]), G.y + e.rb * Math.sin(ANG[e.a])]], ease.inOutSine(k));
        d.position.set(p.x, p.y, 0.8);
        d.setOpacity(V * (k < 1 ? 1 : 1 - (k - 1) / 0.15));
      });

      /* 4. the key stays sealed: anything from its side bounces off the glass */
      const cv = seg(u, 0.17, 0.22) * V;
      [caseFill, walls, glass, keyBox, keyCase, lidLine, ...glare].forEach((o) => o.setOpacity(cv));
      const lk = (u - run.T - 0.07) / 0.06;
      if (lk > 0 && lk < 1) {
        const x0 = CX - 0.2, xw = CL + 0.12;
        leak.position.x = lk < 0.5 ? lerp(x0, xw, ease.inQuad(lk * 2)) : lerp(xw, x0 - 0.25, ease.outQuad(lk * 2 - 1));
        leak.setOpacity(V * (lk < 0.5 ? 1 : 2 - 2 * lk));
      } else leak.setOpacity(0);
      const hit = clamp(1 - Math.abs(lk - 0.5) / 0.12);
      glass.setWidth(2 + 2.5 * hit);
      col(glass, hit > 0.05 ? 'ink' : 'muted');
      leakLab.setOpacity(pulse(u, run.T + 0.095, run.T + 0.2, 0.02) * V);

      /* 5. only after the needle stops: unlock, open, reveal the key, stamp the result */
      let open = 0, unlock = 0, reveal = 0, wk = 0;
      for (const q of RUNS) {
        open += inRun(u, q, q.lid, q.lid + 0.035);
        unlock += seg(u, q.lid - 0.014, q.lid) * (1 - seg(u, q.end + 0.03, q.end + 0.045));
        reveal += inRun(u, q, q.lid + 0.015, q.lid + 0.045);
        wk += inRun(u, q, q.lid + 0.03, q.lid + 0.06);
      }
      lid.rotation.z = 1.05 * open;
      lockBody.setOpacity(cv * (1 - unlock));
      shackle.setOpacity(cv * (1 - unlock));
      shackle.position.y = 0.09 * unlock;
      keyLab.setOpacity(cv * (0.2 + 0.8 * reveal));
      keyG.position.y = CY - 0.15 + 0.22 * reveal;
      wedge.setOpacity(wk * V);
      exactLab.position.y = G.y + slotR + (12 * ts) / ppu;
      exactLab.setOpacity(wk * V);

      chip.setText(run.call);
      col(chip, run.ok ? 'ink' : 'muted');
      chip.setOpacity(inRun(u, run, run.done, run.done + 0.02) * V);
      grade.setProgress(seg(u, run.lid + 0.04, run.lid + 0.065));
      grade.setOpacity(inRun(u, run, run.lid + 0.035, run.lid + 0.045) * V);
      const sv = inRun(u, run, run.lid + 0.062, run.lid + 0.072) * V;
      stamp.scale.setScalar(0.4 + 0.6 * seg(u, run.lid + 0.062, run.lid + 0.085, ease.outBack));
      col(sRing, run.ok ? 'good' : 'bad');
      [sFill, sRing, tick, crossA, crossB].forEach((o) => o.setOpacity(sv));
      tick.visible = run.ok;
      crossA.visible = crossB.visible = !run.ok;
    },
    dispose() { dead = true; },
  };
}
