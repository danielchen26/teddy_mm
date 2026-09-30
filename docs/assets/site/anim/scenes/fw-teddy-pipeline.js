/*
 * fw-teddy-pipeline — anm-framework.html #teddy-pipeline, "TEDDY embeds, our head predicts, ANM decides".
 * One CITE-seq cell, left to right. Its RNA enters frozen TEDDY-G 70M (an RNA model, with a lock;
 * the block never changes), which gives a 512-number embedding (a ring of beads). Our head, also fed
 * a size factor, turns the embedding into 134 predicted surface proteins (a grid of 134 dots). Nine of
 * them (B, T, myeloid) light up and slide into ANM, which gives a call or no call. The same cell's
 * measured protein runs on a lower track, under a dashed wall it never crosses, and meets the call
 * only at the grade. 10 s seamless loop. Where the 9 sit inside the grid is illustrative (not the real
 * protein order); no values are shown. Numbers shown (512, 134, 9, 70M) are the page's.
 * layout() sizes chips and the ANM box to their text and keeps the loader's play/pause pill corner
 * clear; ANM's lineage names show only where they fit (the page's panel legend names them anyway).
 */
/* Every pipeline number this scene prints (stated on the page). Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  embedding: 512,   // TEDDY's embedding size
  predicted: 134,   // surface proteins our head predicts
  panel: 9,         // of them, the proteins ANM reads
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9, center: [0, 0.15] });
  const PERIOD = 10, TOP = 1.3, LOW = -2.7, WALL = -1.55, BH = 1.7, CALLY = -1.05, EDGE = 7.75;
  const X = { cell: -7.38, teddy: -4.83, ring: -2.04, head: 0.69, grid: 3.47, anm: 6.35 };
  const TW = 2.7, HW = 2.6, AW = 2.8, GW = 2.3, RR = 0.65, CR = 0.42, GS = 0.145, SFB = TOP + BH / 2 + 0.85;
  const add = (o) => { scene.add(o); return o; };
  const text = (s, x, y, o) => { const l = ctx.label(s, o); l.position.set(x, y, 5); return add(l); };
  const NAME = { size: 12.5, weight: 600, color: 'ink' };
  const SUB = { size: 11.5, color: 'muted' };
  const sub = (s, x, y, anchor) => text(s, x, y, Object.assign({ anchor }, SUB));
  const block = (x, y, w, h, fill, tint, stroke, r = 0.2) => {
    const b = add(ctx.box(w, h, { color: fill, opacity: tint, stroke, strokeWidth: 1.5, radius: r }));
    b.position.set(x, y, 0);
    if (b.outline) b.outline.setOpacity(1);
    return b;
  };
  const conn = (a, b) => add(ctx.arrow([a[0], a[1], 0.2], [b[0], b[1], 0.2], { color: 'faint', width: 1.5, head: 6 }));
  const circ = (cx, cy, r, n = 48, a0 = 0, a1 = 360) => Array.from({ length: n }, (_, k) => {
    const a = ((a0 + ((a1 - a0) * k) / (a1 - a0 >= 360 ? n : n - 1)) * Math.PI) / 180;
    return [cx + r * Math.cos(a), cy + r * Math.sin(a), 0];
  });
  const floors = []; // [dot, world radius, min px]: small dots keep a visible size on phones
  const speck = (x, y, r, color, min, opacity = 1) => {
    const d = add(ctx.dot([x, y, 0.3], { r, color, opacity }));
    floors.push([d, r, min]);
    return d;
  };
  const mcv = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, size, weight = 500) => { // CSS px, as ctx.label rasterises it
    if (!mcv) return s.length * size * 0.6 * ctx.textScale;
    mcv.font = `${weight} ${(size * ctx.textScale).toFixed(2)}px ${family}`;
    return mcv.measureText(s).width + 3;
  };
  const textH = (size) => size * ctx.textScale * 1.22 + 2;

  /* ── the cell: an outline dotted with RNA ── */
  add(ctx.line(circ(X.cell, TOP, CR), { color: 'ink', width: 1.5, closed: true }));
  const rnd = ctx.rand(5);
  for (let i = 0; i < 7; i++) {
    const a = (i / 7) * Math.PI * 2 + rnd(), rr = 0.1 + 0.16 * rnd();
    speck(X.cell + rr * Math.cos(a), TOP + rr * Math.sin(a), 0.045, 'muted', 1.2);
  }
  text('cell', X.cell, TOP + CR + 0.16, Object.assign({}, SUB, { anchor: 'bottom', color: 'ink' }));

  /* ── frozen TEDDY: tinted block + lock badge (never animated) ── */
  block(X.teddy, TOP, TW, BH, 'teddy', 0.12, 'teddy');
  text('TEDDY-G\n70M', X.teddy, TOP, NAME);
  sub('frozen\nRNA model', X.teddy, TOP - BH / 2 - 0.12, 'top');
  const lock = add(ctx.group());
  lock.position.set(X.teddy - TW / 2 - 0.1, TOP + BH / 2 + 0.1, 1);
  lock.add(ctx.dot([0, 0, 0], { r: 0.25, color: 'card' }));
  lock.add(ctx.line(circ(0, 0, 0.25, 32), { color: 'teddy', width: 1.5, closed: true }));
  const lockBody = ctx.box(0.2, 0.15, { color: 'teddy', radius: 0.03 });
  lockBody.position.set(0, -0.045, 0.1);
  lock.add(lockBody, ctx.line(circ(0, 0.03, 0.062, 12, 0, 180), { color: 'teddy', width: 1.5 }));

  /* ── the embedding: a ring of beads around "512" ── */
  const beads = Array.from({ length: 32 }, (_, i) => {
    const a = Math.PI / 2 - (i / 32) * Math.PI * 2;
    return speck(X.ring + RR * Math.cos(a), TOP + RR * Math.sin(a), 0.045, ctx.color('faint'), 1.1);
  });
  text(String(NUMBERS.embedding), X.ring, TOP, { size: 11, weight: 600, color: 'ink' });
  sub('embedding', X.ring, TOP + RR + 0.16, 'bottom');

  /* ── our head, with the size factor dropping in ── */
  block(X.head, TOP, HW, BH, 'soft', 1, 'muted');
  text('our head', X.head, TOP, NAME);
  const sfChip = block(X.head, SFB + 0.28, 1.2, 0.56, null, 1, 'faint', 0.28);
  const sfLabel = sub('size factor', X.head, SFB + 0.28, 'center');

  /* ── 134 predicted proteins: a 10-wide grid; 9 of them are ANM's panel ── */
  const gx0 = X.grid - 4.5 * GS, gy0 = TOP + 6.5 * GS, jit = ctx.rand(11);
  const LIN = ['b', 't', 'm'];
  const PICK = { 12: 0, 57: 0, 96: 0, 25: 1, 70: 1, 113: 1, 38: 2, 83: 2, 127: 2 }; // index → lineage row
  const grid = Array.from({ length: NUMBERS.predicted }, (_, i) => {
    const c = i % 10, r = Math.floor(i / 10);
    const d = add(ctx.dot([gx0 + c * GS, gy0 - r * GS, 0.3], { r: 0.042, color: ctx.color('faint') }));
    return { d, c, j: jit(), row: i in PICK ? PICK[i] : -1 };
  });
  sub(`${NUMBERS.predicted} proteins`, X.grid, gy0 - 13 * GS - 0.2, 'top');

  /* ── ANM: three lineage rows of three slots (placed by layout), then a call ── */
  const anm = block(X.anm, TOP, AW, 2.1, 'accent', 0.1, 'accent');
  const anmTitle = text('ANM', X.anm, TOP + 0.72, NAME);
  const READS = `reads ${NUMBERS.panel}`;
  const reads = sub(READS, X.anm, TOP + 1.2, 'bottom');
  const slotXY = [[], [], []];
  const slots = [], linLabels = [], movers = [];
  ['B', 'T', 'myeloid'].forEach((name, row) => {
    for (let k = 0; k < 3; k++) { slotXY[row].push([0, 0]); slots.push(add(ctx.dot([0, 0, 0.3], { r: 0.085, color: LIN[row], opacity: 0.22 }))); }
    linLabels.push(text(name, X.anm, TOP, { size: 11.5, color: 'ink', anchor: 'left' }));
  });
  grid.forEach((g) => {
    if (g.row < 0) return;
    const k = movers.filter((m) => m.row === g.row).length;
    movers.push({ row: g.row, from: g.d.position.clone(), to: slotXY[g.row][k], d: add(ctx.dot([0, 0, 2], { r: 0.085, color: LIN[g.row], opacity: 0 })) });
  });
  const callChip = block(X.anm, CALLY, 2.4, 0.62, 'accent', 0, 'accent', 0.31);
  const callLabel = text('call / no call', X.anm, CALLY, { size: 11.5, color: 'ink' });
  const grade = block(X.anm, LOW, GW, 0.8, 'measured', 0, 'measured', 0.22);
  text('grade', X.anm, LOW, NAME);

  /* ── connectors (the two under ANM follow layout) ── */
  conn([X.cell + CR + 0.06, TOP], [X.teddy - TW / 2 - 0.06, TOP]);
  conn([X.teddy + TW / 2 + 0.05, TOP], [X.ring - RR - 0.1, TOP]);
  conn([X.ring + RR + 0.1, TOP], [X.head - HW / 2 - 0.06, TOP]);
  conn([X.head + HW / 2 + 0.05, TOP], [gx0 - 0.12, TOP]);
  conn([gx0 + 9 * GS + 0.12, TOP], [X.anm - AW / 2 - 0.06, TOP]);
  conn([X.head, SFB - 0.04], [X.head, TOP + BH / 2 + 0.06]);
  const toCall = conn([X.anm, 0], [X.anm, -0.5]), toGrade = conn([X.anm, -1.5], [X.anm, -2]);

  /* ── measured protein: its own lower track, below a wall it never crosses ── */
  const M0 = [X.cell, TOP - CR - 0.05], M1 = [X.cell, LOW], M2 = [X.anm - GW / 2 - 0.06, LOW];
  add(ctx.line([[M0[0], M0[1], 0.2], [M1[0], M1[1], 0.2]], { color: 'faint', width: 1.5 }));
  conn(M1, M2);
  const mPath = [M0, M1, [M2[0] - 0.4, LOW]].map((p) => [p[0], p[1], 0.6]);
  const bars = add(ctx.group());
  [-0.16, 0, 0.16].forEach((dx) => bars.add(ctx.line([[dx, -0.17, 0], [dx, 0.17, 0]], { color: 'measured', width: 2, opacity: 0 })));
  sub('measured protein', X.cell + 0.3, LOW - 0.3, 'top-left');
  add(ctx.line([[-6.75, WALL, 0.1], [4.0, WALL, 0.1]], { color: 'muted', width: 1.5, dashed: [6, 5] }));
  text('not an input', -1.1, WALL, { size: 11, color: 'muted', bg: 'card', bgOpacity: 1, pad: 5 });

  /* ── travelling tokens: colour, px, from, to, start, end (fractions of the loop) ── */
  const trips = [];
  const trip = (color, px, a, b, t0, t1) => {
    const r = { d: add(ctx.dot([0, 0, 3], { px, color, opacity: 0 })), pts: [[a[0], a[1], 3], [b[0], b[1], 3]], t0, t1 };
    trips.push(r);
    return r;
  };
  [-0.12, 0.1, -0.02].forEach((dy, i) => trip('muted', 2.8, [X.cell + 0.12, TOP + dy], [X.teddy - TW / 2 - 0.05, TOP + dy], 0.03 + i * 0.02, 0.12 + i * 0.02));
  trip('teddy', 4, [X.teddy + TW / 2, TOP], [X.ring - RR - 0.06, TOP], 0.15, 0.21);
  trip('teddy', 4, [X.ring + RR + 0.06, TOP], [X.head - HW / 2, TOP], 0.33, 0.4);
  trip('muted', 3.5, [X.head, SFB], [X.head, TOP + BH / 2], 0.33, 0.4);
  trip('ink', 4, [X.head + HW / 2, TOP], [gx0 - 0.08, TOP], 0.42, 0.47);
  const callTrip = trip('accent', 4, [0, 0], [0, 0], 0.75, 0.79), gradeTrip = trip('accent', 4, [0, 0], [0, 0], 0.81, 0.86);

  const P = { grid: 1.2, lit: 1.7, slot: 2.4 };
  function layout() {
    const ppu = ctx.ppu(), px = (v) => v / ppu;
    P.grid = Math.max(0.042 * ppu, 1.15);
    P.lit = Math.max(0.062 * ppu, 1.7);
    P.slot = Math.max(0.085 * ppu, 2.4);
    floors.forEach(([d, r, min]) => d.setPx(Math.max(r * ppu, min)));
    slots.forEach((d) => d.setPx(P.slot));
    lock.scale.setScalar(clamp(13 / (0.5 * ppu), 1, 1.6));
    // the loader's play/pause pill sits in the top-right corner (≈ 84 × 36 px with its inset)
    const zoneY = camera.position.y + camera.top / camera.zoom - px(42);
    const zoneX = camera.position.x + camera.right / camera.zoom - px(92);
    // chips: at least their text plus a margin
    const sfW = Math.min(3.2, px(textW('size factor', 11.5) + 20)), sfH = Math.max(0.56, px(textH(11.5) + 6));
    sfChip.setSize(sfW, sfH, sfH / 2); sfChip.position.y = SFB + sfH / 2; sfLabel.position.y = SFB + sfH / 2;
    const cw = px(textW('call / no call', 11.5) + 18), ch = Math.max(0.62, px(textH(11.5) + 6)), cx = Math.min(X.anm, EDGE - cw / 2);
    const gh = Math.max(0.8, px(textH(12.5) + 7));
    grade.setSize(GW, gh, 0.22);
    // ANM: title over three rows; lineage names only if they fit the box and stay clear of the pill
    const titleH = px(textH(12.5)), r = px(P.slot), pitch = Math.max(0.26, px(2 * P.slot + 2));
    const labW = px(textW('myeloid', 11.5)), gap = 0.16;
    let rowH = Math.max(0.43, px(textH(11.5))), AH = Math.max(2.1, 0.3 + titleH + 3 * rowH);
    const show = 2 * pitch + 2 * r + gap + labW + 0.2 <= AW && TOP + AH / 2 <= zoneY;
    if (!show) { rowH = Math.max(0.43, px(2 * P.slot + 4)); AH = Math.max(2.1, 0.3 + titleH + 3 * rowH); }
    const ay = Math.min(TOP, zoneY - AH / 2 - 0.02), top = ay + AH / 2, spare = (AH - 0.3 - titleH - 3 * rowH) / 2;
    anm.setSize(AW, AH, 0.2); anm.position.y = ay;
    anmTitle.position.y = top - 0.14 - spare - titleH / 2;
    const gw = 2 * pitch + 2 * r + (show ? gap + labW : 0), x0 = X.anm - gw / 2 + r;
    for (let row = 0; row < 3; row++) {
      const y = top - 0.14 - spare - titleH - 0.06 - rowH / 2 - row * rowH;
      for (let k = 0; k < 3; k++) { slotXY[row][k][0] = x0 + k * pitch; slotXY[row][k][1] = y; slots[row * 3 + k].position.set(x0 + k * pitch, y, 0.3); }
      linLabels[row].position.set(x0 + 2 * pitch + r + gap, y, 5); linLabels[row].visible = show;
    }
    reads.position.y = top + 0.1;
    reads.visible = top + 0.1 + px(textH(11.5)) <= zoneY || X.anm + px(textW(READS, 11.5)) / 2 <= zoneX;
    // call chip under ANM (slides left when its text is wide), then the grade
    const cy = Math.min(CALLY, top - AH - 0.45 - ch / 2);
    callChip.setSize(cw, ch, ch / 2); callChip.position.set(cx, cy, 0); callLabel.position.set(cx, cy, 5);
    const a0 = [X.anm, top - AH - 0.05], a1 = [X.anm, cy + ch / 2 + 0.06], b0 = [X.anm, cy - ch / 2 - 0.05], b1 = [X.anm, LOW + gh / 2 + 0.06];
    toCall.set([a0[0], a0[1], 0.2], [a1[0], a1[1], 0.2]); toGrade.set([b0[0], b0[1], 0.2], [b1[0], b1[1], 0.2]);
    callTrip.pts = [[a0[0], a0[1] + 0.05, 3], [a1[0], a1[1] - 0.06, 3]];
    gradeTrip.pts = [[b0[0], b0[1] + 0.05, 3], [b1[0], b1[1] - 0.06, 3]];
  }
  let alive = true;
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (alive) { layout(); ctx.requestRender(); } });
  const tmp = new THREE.Vector3();

  return {
    scene, camera, period: PERIOD, still: 9.0,
    resize: layout,
    update(t) {
      const u = ctx.loopT(t, PERIOD), c = ctx.colors;
      const off = 1 - seg(u, 0.93, 0.99, ease.inOutSine);
      for (const r of trips) {
        r.d.position.copy(ctx.along(r.pts, seg(u, r.t0, r.t1, ease.inOutSine)));
        r.d.setOpacity(pulse(u, r.t0, r.t1, 0.012));
      }
      // TEDDY's output: the ring fills clockwise
      beads.forEach((d, i) => {
        const k = seg(u, 0.21 + (0.1 * i) / 32, 0.24 + (0.1 * i) / 32) * off;
        d.material.color.copy(c.faint).lerp(c.teddy, k);
        d.setOpacity(lerp(0.55, 1, k));
      });
      // the head's output: 134 dots fill left to right, then the 9 panel proteins take lineage colours
      const hl = seg(u, 0.58, 0.62) * off;
      for (const g of grid) {
        const a = 0.46 + (0.09 * g.c) / 9 + 0.015 * g.j;
        const k = seg(u, a, a + 0.03) * off;
        g.d.material.color.copy(c.faint).lerp(c.ink, k);
        if (g.row >= 0) {
          g.d.material.color.lerp(c[LIN[g.row]], hl);
          g.d.setPx(lerp(P.grid, P.lit, hl)).setOpacity(lerp(lerp(0.45, 0.8, k), 1, hl));
        } else g.d.setPx(P.grid).setOpacity(lerp(0.45, 0.8, k));
      }
      // ANM reads the 9: each copy arcs into its lineage slot
      movers.forEach((m, i) => {
        const s = seg(u, 0.62 + i * 0.01, 0.71 + i * 0.01, ease.inOutCubic);
        tmp.set(lerp(m.from.x, m.to[0], s), lerp(m.from.y, m.to[1], s) + 0.35 * Math.sin(Math.PI * s), 2);
        m.d.position.copy(tmp);
        m.d.setPx(lerp(P.lit, P.slot, s)).setOpacity(s > 0 ? off : 0);
      });
      anm.material.opacity = 0.1 + 0.12 * pulse(u, 0.72, 0.84, 0.04);
      callChip.material.opacity = 0.16 * seg(u, 0.78, 0.81) * off;
      // measured protein: leaves with the cell, waits at the grade, never enters upstream
      bars.position.copy(ctx.along(mPath, seg(u, 0.04, 0.86, ease.inOutSine)));
      const mo = seg(u, 0.03, 0.06) * off;
      bars.children.forEach((b) => b.setOpacity(mo));
      grade.material.opacity = 0.14 * seg(u, 0.86, 0.89) * off;
    },
    dispose() {
      alive = false;
      floors.length = 0; grid.length = 0; movers.length = 0; trips.length = 0; slots.length = 0; linLabels.length = 0;
    },
  };
}
