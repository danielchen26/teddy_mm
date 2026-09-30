/*
 * idx-hero-pipeline — index.html #overview, "One cell, end to end" (hero card, 16:10).
 * One CITE-seq cell splits: its RNA streams into frozen TEDDY (padlocked), which gives a 512-number
 * embedding (a ring of beads); our head turns the embedding into 134 predicted proteins (a strip of
 * thin bars). Nine of them (B, T, myeloid) lift out as evidence cards into ANM; all nine enter at the
 * same instant, three lineage scores rise against the bar, and the one that clears it is the call.
 * The cell's measured proteins drop into a sealed vault at the start and only grade the call at the end.
 * 10 s seamless loop. Bar heights, scores and which lineage wins are illustrative, not data; the only
 * numbers shown (512, 134, 9) are the page's.
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 10 });
  const PERIOD = 10, Y0 = 2.6, SB = 1.8, BASE = -3.85, BAR = -1.95;
  const X = { cell: -6.5, teddy: -3.5, ring: -0.7, head: 1.8, s0: 3.6, s1: 7.5 };
  const CR = 0.45, TW = 2.3, TH = 1.8, RR = 0.65, HW = 2.5, HH = 1.4;
  const P = { x: 4.15, y: -2.625, w: 6.9, h: 4.25 };        // ANM panel
  const V = { x: -6.5, y: -2.4, w: 2.4, h: 2.0 };            // sealed vault of measured proteins
  const B = { x: -1.6, y: -2.4, w: 2.1, h: 0.62 };           // call badge
  const COLS = [
    { k: 'b', name: 'B', x: 2.55, score: -2.65 },
    { k: 't', name: 'T', x: 4.35, score: -1.45 },
    { k: 'm', name: 'myeloid', x: 6.15, score: -3.0 },
  ];
  const WIN = 1;
  const add = (o) => { scene.add(o); return o; };
  const text = (s, x, y, o) => { const l = ctx.label(s, o); l.position.set(x, y, 5); return add(l); };
  const NAME = { size: 13, weight: 600, color: 'ink' };
  const sub = (o) => Object.assign({ size: 11.5, color: 'muted' }, o);
  const block = (x, y, w, h, fill, tint, stroke, o = {}) => {
    const b = add(ctx.box(w, h, Object.assign({ color: fill, opacity: tint, stroke, strokeWidth: 1.5, radius: 0.2 }, o)));
    b.position.set(x, y, 0);
    if (b.outline) b.outline.setOpacity(1);
    return b;
  };
  const conn = (a, b) => add(ctx.arrow([a[0], a[1], 0.2], [b[0], b[1], 0.2], { color: 'faint', width: 1.5, head: 6 }));
  const circ = (cx, cy, r, n = 48, a0 = 0, a1 = 360) => Array.from({ length: n }, (_, k) => {
    const a = ((a0 + ((a1 - a0) * k) / (a1 - a0 >= 360 ? n : n - 1)) * Math.PI) / 180;
    return [cx + r * Math.cos(a), cy + r * Math.sin(a), 0];
  });
  const locks = [];
  const makeLock = (x, y, color) => {
    const g = add(ctx.group());
    g.position.set(x, y, 1);
    g.add(ctx.dot([0, 0, 0], { r: 0.25, color: 'card' }), ctx.line(circ(0, 0, 0.25, 32), { color, width: 1.5, closed: true }));
    const body = ctx.box(0.2, 0.15, { color, radius: 0.03 });
    body.position.set(0, -0.045, 0.1);
    const hinge = ctx.group(); // shackle pivots on its right foot when the lock opens
    hinge.position.set(0.062, 0.03, 0.1);
    hinge.add(ctx.line(circ(-0.062, 0, 0.062, 12, 0, 180), { color, width: 1.5 }));
    g.add(body, hinge);
    locks.push(g);
    return hinge;
  };

  /* ── the cell: RNA beads and measured-protein bars inside one outline ── */
  add(ctx.dot([X.cell, Y0, 0.05], { r: CR, color: 'soft' }));
  add(ctx.line(circ(X.cell, Y0, CR), { color: 'ink', width: 1.5, closed: true }));
  text('cell', X.cell, Y0 + CR + 0.14, sub({ color: 'ink', anchor: 'bottom' }));
  const rna = [[-0.22, 0.14], [-0.02, 0.26], [-0.28, -0.08], [0.02, 0.05], [-0.1, -0.24]].map(([dx, dy], i) => ({
    d: add(ctx.dot([0, 0, 1], { px: 2.6, color: 'muted' })),
    home: [X.cell + dx, Y0 + dy], to: [X.teddy - TW / 2 - 0.02, Y0 + [0.12, -0.12, 0.04, -0.04, 0.2][i]],
  }));
  text('RNA', (X.cell + CR + X.teddy - TW / 2) / 2, Y0 + 0.2, sub({ size: 11, anchor: 'bottom' }));
  const meas = add(ctx.group());
  [[-0.2, 0.5], [0, 0.8], [0.2, 0.36]].forEach(([dx, h]) => meas.add(ctx.line([[dx, 0, 0], [dx, h, 0]], { color: 'measured', width: 2.5 })));
  const MEAS_CELL = [X.cell + 0.18, Y0 - 0.26], MEAS_VAULT = [V.x, V.y - 0.62];

  /* ── frozen TEDDY, the 512 embedding, our head ── */
  block(X.teddy, Y0, TW, TH, 'teddy', 0.14, 'teddy');
  text('TEDDY', X.teddy, Y0, NAME);
  text('frozen', X.teddy, Y0 - TH / 2 - 0.14, sub({ anchor: 'top' }));
  makeLock(X.teddy + TW / 2, Y0 + TH / 2, 'teddy');
  const beads = Array.from({ length: 28 }, (_, i) => {
    const a = Math.PI / 2 - (i / 28) * Math.PI * 2;
    return add(ctx.dot([X.ring + RR * Math.cos(a), Y0 + RR * Math.sin(a), 0.3], { px: 1.5, color: ctx.color('faint') }));
  });
  text('512', X.ring, Y0, { size: 11.5, weight: 600, color: 'ink' });
  text('embedding', X.ring, Y0 - RR - 0.16, sub({ anchor: 'top' }));
  block(X.head, Y0, HW, HH, 'soft', 1, 'muted');
  text('our head', X.head, Y0, Object.assign({}, NAME, { size: 12.5 }));
  const flow = [0, 1, 2].map(() => add(ctx.dot([0, 0, 1], { px: 3.2, color: 'teddy', opacity: 0 })));

  /* ── 134 predicted proteins: one strip of thin bars; 9 of them are ANM's evidence ── */
  const N = 134, sp = (X.s1 - X.s0) / N, rs = ctx.rand(3);
  const hs = Array.from({ length: N }, () => 0.22 + 1.03 * Math.pow(rs(), 1.3));
  const strip = add(ctx.group());
  strip.position.set(X.s0, SB, 0.2);
  const sGeo = ctx.track(new THREE.BufferGeometry());
  const sPos = new THREE.BufferAttribute(new Float32Array(N * 12), 3);
  sGeo.setAttribute('position', sPos);
  sGeo.setIndex(Array.from({ length: N * 6 }, (_, k) => Math.floor(k / 6) * 4 + [0, 1, 2, 0, 2, 3][k % 6]));
  const sMat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, opacity: 0.8 })), 'faint');
  const sMesh = new THREE.Mesh(sGeo, sMat);
  sMesh.frustumCulled = false;
  strip.add(sMesh);
  const buildStrip = (bw) => {
    for (let i = 0; i < N; i++) {
      const x = (i + 0.5) * sp, a = x - bw / 2, b = x + bw / 2, h = hs[i];
      sPos.array.set([a, 0, 0, b, 0, 0, b, h, 0, a, h, 0], i * 12);
    }
    sPos.needsUpdate = true;
  };
  add(ctx.line([[X.s0 - 0.05, SB, 0.1], [X.s1 + 0.05, SB, 0.1]], { color: 'line', width: 1 }));
  // below the strip, on a plate: the top-right corner stays clear for the play/pause button
  text('134 predicted\nproteins', 5.3, SB - 0.14, sub({ color: 'ink', anchor: 'top', bg: 'card', bgOpacity: 1, pad: 2 }));
  const PICKS = [[14, 0], [61, 0], [104, 0], [29, 1], [75, 1], [118, 1], [42, 2], [90, 2], [128, 2]];
  const cards = PICKS.map(([i, c], j) => {
    const slot = PICKS.slice(0, j).filter((p) => p[1] === c).length, x = (i + 0.5) * sp;
    const bar = ctx.line([[x, 0, 0.1], [x, hs[i], 0.1]], { color: COLS[c].k, width: 2.2, opacity: 0 });
    strip.add(bar);
    return {
      bar, c, j, card: add(ctx.box(0.66, 0.28, { color: COLS[c].k, radius: 0.07, opacity: 0 })),
      from: [X.s0 + x, SB + hs[i] * 0.6], to: [COLS[c].x, BASE + 0.35 + slot * 0.38],
    };
  });

  /* ── ANM: three lineage scores against the question's bar ── */
  block(P.x, P.y, P.w, P.h, 'accent', 0.06, 'accent', { radius: 0.25 });
  text('ANM', P.x - P.w / 2 + 0.3, P.y + P.h / 2 - 0.45, Object.assign({}, NAME, { anchor: 'left' }));
  const once = text('all 9 at once', P.x + P.w / 2 - 0.3, P.y + P.h / 2 - 0.45, sub({ size: 11, anchor: 'right', opacity: 0 }));
  add(ctx.line([[COLS[0].x - 0.8, BASE, 0.2], [COLS[2].x + 0.8, BASE, 0.2]], { color: 'line', width: 1 }));
  const zone = add(ctx.box(5.9, 0.75, { color: 'accent', opacity: 0.07, radius: 0.05 }));
  zone.position.set(4.15, BAR + 0.375, 0.1);
  const barPts = [[1.2, BAR, 0.4], [7.1, BAR, 0.4]];
  add(ctx.line(barPts, { color: 'muted', width: 1.5, dashed: [6, 4] }));
  const barHit = add(ctx.line(barPts, { color: 'accent', width: 2.5, opacity: 0 }));
  text('bar', 1.25, BAR + 0.1, sub({ size: 11, anchor: 'bottom-left' }));
  COLS.forEach((c) => {
    c.stem = add(ctx.line([[c.x, BASE, 0.3], [c.x, c.score, 0.3]], { color: c.k, width: 2, opacity: 0 }));
    c.node = add(ctx.dot([c.x, BASE, 0.5], { r: 0.2, color: c.k, opacity: 0 }));
    text(c.name, c.x, BASE - 0.2, sub({ color: 'ink', anchor: 'top' }));
  });

  /* ── the call, and the sealed vault that grades it last ── */
  const slot = block(B.x, B.y, B.w, B.h, null, 1, 'faint', { radius: B.h / 2, dashed: [4, 4] });
  const badge = block(B.x, B.y, B.w, B.h, 'accent', 0, 'accent', { radius: B.h / 2 });
  const badgeLab = text('call: T', B.x, B.y, { size: 12, weight: 600, color: 'ink', opacity: 0 });
  const callDot = add(ctx.dot([0, 0, 3], { px: 3.5, color: 'accent', opacity: 0 }));
  block(V.x, V.y, V.w, V.h, 'soft', 0.7, 'measured', { radius: 0.22 });
  text('measured', V.x, V.y + V.h / 2 - 0.36, { size: 11, weight: 600, color: 'ink' });
  text('sealed until\nthe call', V.x, V.y - V.h / 2 - 0.14, sub({ anchor: 'top' }));
  const vHinge = makeLock(V.x + V.w / 2, V.y, 'measured'); // on the side the grading link leaves from
  const gx0 = V.x + V.w / 2, gx1 = B.x - B.w / 2 - 0.06;
  const link = add(ctx.line([[gx0, V.y, 0.3], [gx1, V.y, 0.3]], { color: 'measured', width: 1.5, dashed: [5, 4] }));
  const stamp = add(ctx.group());
  stamp.position.set((gx0 + 0.3 + gx1) / 2, V.y, 2);
  const stampParts = [
    ctx.dot([0, 0, 0], { r: 0.3, color: 'card' }),
    ctx.line(circ(0, 0, 0.3, 32), { color: 'good', width: 1.5, closed: true }),
    ctx.line([[-0.13, 0.0, 0.1], [-0.03, -0.1, 0.1], [0.14, 0.11, 0.1]], { color: 'good', width: 2 }),
  ];
  stamp.add(...stampParts);

  /* ── connectors ── */
  conn([X.cell + CR + 0.08, Y0], [X.teddy - TW / 2 - 0.06, Y0]);
  conn([X.teddy + TW / 2 + 0.06, Y0], [X.ring - RR - 0.1, Y0]);
  conn([X.ring + RR + 0.08, Y0], [X.head - HW / 2 - 0.06, Y0]);
  conn([X.head + HW / 2 + 0.06, Y0], [X.s0 - 0.08, Y0]);
  conn([7.3, SB - 0.12], [7.3, P.y + P.h / 2 + 0.06]);
  conn([X.cell, Y0 - CR - 0.08], [V.x, V.y + V.h / 2 + 0.06]);
  conn([P.x - P.w / 2 - 0.06, B.y], [B.x + B.w / 2 + 0.08, B.y]);

  let stampScale = 1;
  const setMeas = (p, s, o) => { meas.position.set(p[0], p[1], 1.5); meas.scale.setScalar(s); meas.children.forEach((l) => l.setOpacity(o)); };

  return {
    scene, camera, period: PERIOD, still: 8.95,
    resize() {
      const ppu = ctx.ppu();
      beads.forEach((d) => d.setPx(clamp(0.05 * ppu, 1.2, 3)));
      locks.forEach((g) => g.scale.setScalar(clamp(13 / (0.5 * ppu), 1, 1.6)));
      stampScale = clamp(16 / (0.6 * ppu), 1, 1.7);
      buildStrip(Math.max(sp * 0.55, 1 / ppu));
    },
    update(t) {
      const u = ctx.loopT(t, PERIOD), c = ctx.colors;
      const off = 1 - seg(u, 0.92, 0.98, ease.inOutSine);   // everything that moved fades out
      const back = seg(u, 0.95, 0.995, ease.inOutSine);     // and the cell refills for the next loop
      // 1 · the cell splits: RNA to TEDDY, measured proteins down into the vault
      rna.forEach((r, i) => {
        const a = 0.03 + 0.02 * i, s = seg(u, a, a + 0.1, ease.inOutSine);
        const late = u > 0.9;
        r.d.position.set(late ? r.home[0] : lerp(r.home[0], r.to[0], s), late ? r.home[1] : lerp(r.home[1], r.to[1], s), 1);
        r.d.setOpacity(late ? back : 1 - seg(u, a + 0.09, a + 0.1));
      });
      const drop = seg(u, 0.04, 0.15, ease.inOutCubic);
      if (u > 0.95) setMeas(MEAS_CELL, 0.4, back);
      else setMeas([lerp(MEAS_CELL[0], MEAS_VAULT[0], drop), lerp(MEAS_CELL[1], MEAS_VAULT[1], drop)], lerp(0.4, 1, drop), 1 - seg(u, 0.91, 0.95));
      // 2 · frozen TEDDY gives the embedding: the ring fills clockwise
      beads.forEach((d, i) => {
        const a = 0.16 + (0.09 * i) / beads.length, k = seg(u, a, a + 0.03) * off;
        d.material.color.copy(c.faint).lerp(c.teddy, k);
        d.setOpacity(lerp(0.5, 1, k));
      });
      flow.forEach((d, i) => {
        const a = 0.27 + 0.025 * i, s = seg(u, a, a + 0.05, ease.inOutSine);
        d.position.set(lerp(X.ring + RR + 0.1, X.head - HW / 2 - 0.05, s), Y0, 1);
        d.setOpacity(pulse(u, a, a + 0.05, 0.01));
      });
      // 3 · our head fans out 134 predicted proteins; the 9 evidence proteins light up
      const g = seg(u, 0.34, 0.44, ease.outCubic) * off;
      strip.scale.set(lerp(0.12, 1, g), g, 1);
      strip.visible = g > 0.001;
      const hl = seg(u, 0.44, 0.47) * off;
      // 4 · they lift out together as cards, then all 9 enter ANM at the same instant (no order)
      const merge = seg(u, 0.62, 0.66, ease.inOutCubic);
      const a = 0.48, lift = seg(u, a, a + 0.03, ease.outCubic), s = seg(u, a + 0.02, a + 0.12, ease.inOutCubic);
      for (const cd of cards) {
        cd.bar.setOpacity(hl);
        const x = lerp(cd.from[0], cd.to[0], s), y = lerp(cd.from[1] + 0.45 * lift, cd.to[1], s) + 0.5 * Math.sin(Math.PI * s);
        cd.card.position.set(lerp(x, COLS[cd.c].x, merge), lerp(y, BASE, merge), 2 + cd.j * 0.01);
        cd.card.scale.setScalar(lerp(1, 0.35, merge));
        cd.card.material.opacity = seg(u, a, a + 0.012) * (1 - merge) * 0.9;
      }
      once.setOpacity(pulse(u, 0.6, 0.72, 0.02));
      // 5 · the three lineage scores rise; the one that clears the bar is the call
      const rise = seg(u, 0.66, 0.745, ease.inOutCubic), called = seg(u, 0.74, 0.77);
      COLS.forEach((col, i) => {
        col.node.position.y = lerp(BASE, col.score, rise);
        col.node.setOpacity(merge * off * (i === WIN ? 1 : lerp(1, 0.5, called)));
        col.stem.setProgress(rise).setOpacity(0.45 * merge * off);
      });
      barHit.setOpacity(0.9 * pulse(u, 0.715, 0.8, 0.03));
      const cs = seg(u, 0.725, 0.765, ease.inOutSine);
      callDot.position.set(lerp(P.x - P.w / 2 - 0.05, B.x + B.w / 2 + 0.05, cs), B.y, 3);
      callDot.setOpacity(pulse(u, 0.725, 0.765, 0.008));
      badge.material.opacity = 0.16 * called * off;
      badge.outline.setOpacity(called * off);
      slot.outline.setOpacity(1 - called * off);
      badgeLab.setOpacity(called * off);
      // 6 · only now the vault opens and the measured proteins grade the call
      const open = seg(u, 0.79, 0.83) * off;
      vHinge.rotation.z = 0.75 * open;
      vHinge.position.y = 0.03 + 0.05 * open;
      link.setProgress(seg(u, 0.8, 0.86, ease.inOutSine)).setOpacity(off);
      const sk = seg(u, 0.855, 0.89, ease.outBack) * off;
      stamp.visible = sk > 0.01;
      stamp.scale.setScalar(Math.max(0.01, sk) * stampScale);
    },
    dispose() {
      locks.length = 0; cards.length = 0; beads.length = 0; rna.length = 0; flow.length = 0;
    },
  };
}
