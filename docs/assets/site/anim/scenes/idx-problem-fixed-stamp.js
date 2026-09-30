/*
 * idx-problem-fixed-stamp (index.html #problem, "Read by one fixed rule, the pipeline has six gaps").
 * Cells ride a belt through TEDDY + our head (RNA in, predicted proteins out; both unchanged) and under one
 * stamp, the fixed rule. It comes down every 1.5 s on every cell the same way and leaves the same mark,
 * whatever question is pinned above it (the note cycles through questions the page names; the stamp never
 * changes). Six doors, one per gap on the page (GAP 1–6), light up one by one and are tried: the door shakes,
 * its lock swings and turns red, and none opens. Illustrative only: no values. 12 s seamless loop.
 * Two layouts: the page's 21:6 strip, and 16:9 on narrow screens (the page's --ar-sm).
 */
export default function create(ctx) {
  const { THREE, ease, seg, clamp } = ctx;
  const C = ctx.colors;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 6 });
  const PERIOD = 12, BEAT = 1.5, N = 8, STAMP = 5, MAP = 3, R = 0.27, TRAVEL = 0.5;
  const BW = 1.1, BH = 0.62, HH = 0.4, KH = 0.2;
  const QS = ['stricter?', 'key markers?', 'B/T-priority?'];
  const GAPS = ['change\nquestion', 'missing\ninput', 'explain\na call', 'where to\ntrust', 'price of a\nquestion', 'embedding\nenough?'];
  const LAY = {
    wide: { fw: 21, fh: 6, sp: 1.3, x0: -10.4, belt: -1.32, mapW: 1.9, noteY: 1.88, dx0: -0.35, dsp: 1.95, floor: -1.32, dw: 1.3, dh: 2.7 },
    narrow: { fw: 16, fh: 9, sp: 1.8, x0: -7.2, belt: 0.3, mapW: 2.2, noteY: 3.7, dx0: -6.25, dsp: 2.5, floor: -2.95, dw: 1.7, dh: 1.85 },
  };
  let L = LAY.wide, cellY = L.belt + 0.35, dead = false;
  const put = (o, x, y, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };

  /* belt and cells: a hollow ring per cell, and the stamp's mark (identical on every cell) */
  const belt = put(ctx.line([[0, 0], [1, 0]], { color: 'line', width: 1.5 }), 0, 0);
  const cells = [];
  for (let j = 0; j < N; j++) {
    const g = put(ctx.group(), 0, 0, 1);
    const ring = put(ctx.dot([0, 0], { r: R, hollow: true, ring: 0.3, color: 'faint', order: 1 }), 0, 0, 0, g);
    const mark = put(ctx.box(0.2, 0.2, { color: 'teddy', radius: 0.04, order: 1 }), 0, 0, 0.1, g);
    cells.push({ g, ring, mark });
  }

  /* TEDDY + our head: an opaque box the cells pass behind */
  const map = put(ctx.group(), 0, 0, 2);
  const mapBase = put(ctx.box(1, 1, { color: 'card', stroke: 'teddy', order: 3 }), 0, 0, 0, map);
  const mapTint = put(ctx.box(1, 1, { color: 'teddy', opacity: 0.16, order: 4 }), 0, 0, 0.05, map);
  const mapLab = put(ctx.label('TEDDY + head', { size: 12.5, weight: 600, anchor: 'top' }), 0, 0, 5);
  const stampLab = put(ctx.label('fixed rule', { size: 12.5, weight: 600, anchor: 'top' }), 0, 0, 5);

  /* the stamp: origin at the bottom of its rubber pad */
  const stamp = put(ctx.group(), 0, 0, 3);
  put(ctx.box(BW, BH, { color: 'teddy', opacity: 0.16, radius: 0.1, order: 5 }), 0, BH / 2, 0, stamp);
  put(ctx.box(BW, BH, { color: null, stroke: 'teddy', radius: 0.1, order: 6 }), 0, BH / 2, 0.01, stamp);
  put(ctx.box(BW - 0.1, 0.1, { color: 'teddy', radius: 0.03, order: 6 }), 0, 0.05, 0.02, stamp);
  put(ctx.box(0.22, HH, { color: 'teddy', radius: 0.05, order: 5 }), 0, BH + HH / 2, 0, stamp);
  put(ctx.box(0.6, KH, { color: 'teddy', radius: 0.1, order: 5 }), 0, BH + HH + KH / 2, 0, stamp);

  /* the pinned question: changes; the stamp does not */
  const notes = QS.map((q) => put(ctx.label(q, { size: 12.5, weight: 600, bg: 'grid', bgOpacity: 1, pad: 6 }), 0, 0, 5));
  const tag = put(ctx.label('question', { size: 11.5, color: 'muted', anchor: 'right' }), 0, 0, 5);

  /* the wall of six doors, one per gap */
  const wall = put(ctx.box(1, 1, { color: 'soft', radius: 0.2 }), 0, 0, -1);
  const shacklePts = [[0.075, 0.06], [0.075, 0.14]];
  for (let k = 0; k <= 12; k++) { const a = (Math.PI * k) / 12; shacklePts.push([0.075 * Math.cos(a), 0.14 + 0.075 * Math.sin(a)]); }
  shacklePts.push([-0.075, 0.06]);
  const doors = GAPS.map((name, i) => {
    const g = put(ctx.group(), 0, 0, 0);
    const base = put(ctx.box(1, 1, { color: 'card', stroke: 'line', radius: 0.12, order: 1 }), 0, 0, 0, g);
    const tint = put(ctx.box(1, 1, { color: 'ink', opacity: 0, radius: 0.12, order: 2 }), 0, 0, 0.05, g);
    const ring = put(ctx.box(1, 1, { color: null, stroke: 'ink', strokeWidth: 2, opacity: 0, radius: 0.12, order: 3 }), 0, 0, 0.1, g);
    const num = put(ctx.label(String(i + 1), { size: 11, mono: true, weight: 600, color: 'muted', anchor: 'top' }), 0, 0, 0.2, g);
    const lock = put(ctx.group(), 0, 0, 0.2, g);
    const body = put(ctx.box(0.26, 0.2, { color: 'faint', radius: 0.04, order: 4 }), 0, 0, 0, lock);
    const shackle = put(ctx.line(shacklePts, { color: 'faint', width: 1.8, order: 4 }), 0, 0, 0, lock);
    const lab = put(ctx.label(name, { size: 12, color: 'ink', anchor: 'top' }), 0, 0, 5);
    return { g, base, tint, ring, num, lock, body, shackle, lab, x: 0 };
  });

  /* text widths in CSS px, so labels never collide at small stage widths */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, px, wt = 500) => {
    const lines = s.split('\n');
    if (!mc) return Math.max(...lines.map((l) => l.length)) * px * 0.56;
    mc.font = `${wt} ${px}px ${family}`;
    return Math.max(...lines.map((l) => mc.measureText(l).width));
  };

  function layout() {
    const ppu = ctx.ppu(), ts = ctx.textScale;
    cellY = L.belt + 0.35;
    belt.setPoints([[L.x0, L.belt], [L.x0 + 6.5 * L.sp, L.belt]]);
    const mx = L.x0 + MAP * L.sp, sx = L.x0 + STAMP * L.sp;
    map.position.set(mx, L.belt + 0.5, 2);
    mapBase.setSize(L.mapW, 1);
    mapTint.setSize(L.mapW, 1);
    mapLab.position.set(mx, L.belt - 0.18, 5);
    stampLab.position.set(sx, L.belt - 0.18, 5);
    stamp.position.x = sx;
    notes.forEach((n) => n.position.set(sx, L.noteY, 5));
    const qw = Math.max(...QS.map((q) => textW(q, 12.5 * ts, 600))) + 20;
    tag.position.set(sx - qw / 2 / ppu - 0.18, L.noteY, 5);
    const wx0 = L.dx0 - L.dw / 2 - 0.3, wx1 = L.dx0 + 5 * L.dsp + L.dw / 2 + 0.3, wh = L.dh + 0.2;
    wall.setSize(wx1 - wx0, wh);
    wall.position.set((wx0 + wx1) / 2, L.floor + wh / 2, -1);
    const lw = Math.max(...GAPS.map((s) => textW(s, 12 * ts))), avail = L.dsp * ppu - 10;
    const ls = lw > avail ? Math.max(8.5, (12 * avail) / lw) : 12;
    const lk = clamp(11 / (0.26 * ppu), 1, 2.2);
    doors.forEach((d, i) => {
      d.x = L.dx0 + i * L.dsp;
      d.g.position.set(d.x, L.floor + L.dh / 2, 0);
      [d.base, d.tint, d.ring].forEach((b) => b.setSize(L.dw, L.dh, 0.12));
      d.num.position.set(0, L.dh / 2 - 0.12, 0.2);
      d.lock.position.set(0, -0.2 * L.dh, 0.2);
      d.lock.scale.setScalar(lk);
      d.lab.position.set(d.x, L.floor - 0.16, 5);
      d.lab.setSize(ls);
    });
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 9.15,
    resize(w, h) {
      L = w / h > 2.4 ? LAY.wide : LAY.narrow;
      const fit = camera.userData.animFit, a = w / h;
      fit.width = L.fw; fit.height = L.fh;
      const hh = a > L.fw / L.fh ? L.fh / 2 : L.fw / 2 / a, hw = hh * a;
      camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
      camera.updateProjectionMatrix();
      layout();
    },
    update(t) {
      const u = ctx.loopT(t, PERIOD), tt = u * PERIOD;
      const b = Math.floor(tt / BEAT), ph = tt / BEAT - b;
      /* the stamp: down on each beat, held briefly, up while the belt moves */
      const press = 1 - ease.smooth((Math.min(ph, 1 - ph) - 0.04) / 0.2);
      stamp.position.y = cellY + R + TRAVEL * (1 - press);
      const e = ease.inOutSine((ph - 0.28) / 0.44);
      for (let j = 0; j < N; j++) {
        const c = cells[j], p = (j + b + e) % N;
        c.g.position.set(L.x0 + p * L.sp, cellY, 1);
        const a = clamp(p) * clamp((6.6 - p) / 0.6);
        c.ring.setOpacity(a);
        c.ring.setColor(p >= MAP - 1e-4 ? 'ink' : 'faint');
        const at = Math.abs(p - STAMP) < 1e-4 && ph < 0.5;
        c.mark.visible = at || p > STAMP + 1e-4;
        c.mark.scale.setScalar(at ? Math.max(0.01, ease.outBack(ph / 0.09)) : 1);
        c.mark.setOpacity(a);
      }
      /* the pinned question swaps every 4 s */
      const k3 = Math.min(2, Math.floor(u * 3)), lu = u * 3 - k3;
      const qa = seg(lu, 0, 0.07, ease.inOutSine) * (1 - seg(lu, 0.93, 1, ease.inOutSine));
      notes.forEach((n, i) => n.setOpacity(i === k3 ? qa : 0));
      /* doors light one by one between stamps, are tried, stay shut; all dim before the loop restarts */
      const off = 1 - seg(u, 0.905, 0.975, ease.inOutSine);
      doors.forEach((d, i) => {
        const ti = (i + 0.55) * BEAT, tau = tt - ti;
        const k = seg(tt, ti, ti + 0.4, ease.inOutSine) * off;
        const env = tau > 0 && tau < 0.8 ? Math.sin(2 * Math.PI * 5 * tau) * (1 - tau / 0.8) : 0;
        d.g.position.x = d.x + 0.05 * env;
        d.lock.rotation.z = 0.3 * env;
        d.tint.setOpacity(0.07 * k);
        d.ring.setOpacity(k);
        d.num.setOpacity(0.6 + 0.4 * k);
        d.lab.setOpacity(0.62 + 0.38 * k);
        d.body.material.color.copy(C.faint).lerp(C.bad, k);
        d.shackle.material.uniforms.uColor.value.copy(C.faint).lerp(C.bad, k);
      });
    },
    dispose() {
      dead = true;
      cells.length = 0;
      doors.length = 0;
    },
  };
}
