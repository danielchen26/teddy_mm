/*
 * idx-exp-scope-gate — index.html, #scope (Experiment 4 · scope gate), between .plain and .bgrid.
 * How the test answers “which fixed-rule calls can we act on, and which cells are genuinely hard?”
 * 1 The held-out cells stand as pillars in two mirrored lanes: ranked by ANM’s confidence (teal, up) and
 *   ranked by the fixed rule’s own margin, top score minus runner-up (orange, down). Pillar heights are
 *   illustrative, not data. A coverage slider (100% to 10%) slides a curtain over the least confident
 *   pillars of both lanes: it stops at 80%, where the chip shows ANM’s largest lead (+0.0026), then at 40%.
 *   The accuracy of the kept calls (soft rule) climbs from 0.945 (all cells) to 0.9987 (ANM) and 0.9988
 *   (margin): nearly equal. Arrows show the climb; no in-between values are printed.
 * 2 The hard cells (the call flips under leave-one-out, 4,476) glow red and gather: 73% of them are out of
 *   scope, NK/ILC 3.4× and erythroid 3.7× over-represented; their accuracy 0.818 vs 0.945 for random cells
 *   of the same number (40 draws; the rings re-draw). Dot counts keep the page’s proportions only.
 * Wide layout 21 × 8 when the stage is at least 2:1 (desktop 21:8), else narrow 16 × 9 (phone 16:9).
 * 12 s seamless loop; the still frame shows beat 1 complete (keep 40%, both accuracies).
 */

/* Every pipeline number this scene shows or uses. All are on index.html #scope; refresh after the
 * official-preprocessing rerun (reports/SCOPE_REFINE_PROOF_V2.md, soft rule “O0”). */
const NUMBERS = {
  cells: 16750,        // “16,750 held-out cells”
  covMin: 10,          // “coverage (100% to 10%)”
  leadAt: 80,          // “apart from +0.0026 at 80% (soft) and ≤ 0.0004 elsewhere”
  lead: '+0.0026',
  keep: 40,            // “soft rule, keeping 40%”
  accAll: '0.938',     // “(all cells 0.945)”
  accAnm: '0.9993',    // KPI: ANM confidence, keeping 40%
  accMargin: '0.9994', // KPI: the rule’s margin, keeping 40%
  hard: 4942,          // “calls that flip under leave-one-out (4,476)”
  draws: 40,           // “vs 40 random draws”
  accHard: '0.791',    // “0.818 vs 0.945: accuracy on the 4,476 hard cells vs random cells”
  accRandom: '0.939',
  nkIlc: '3.1×',       // “NK/ILC and erythroid over-represented: 73% of hard cells, out of scope”
  erythroid: '3.4×',
  outOfScope: 73,
};

export default function create(ctx) {
  const { THREE, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 8 });
  const PERIOD = 12;
  const add = (o, z = 0) => { o.position.z = z; scene.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const ln = (o) => add(ctx.line([[0, 0], [1, 0]], o));
  const fmt = (n) => n.toLocaleString('en-US');
  const KEEP = NUMBERS.keep / 100, LEAD = NUMBERS.leadAt / 100, CMIN = NUMBERS.covMin / 100;

  const LAYOUTS = {
    wide: {
      W: 21, H: 8, x0: -9.6, x1: 2.6, yU: 0.12, yD: -0.08, hmax: 2.3, ty: 2.95, above: true, dotPx: 5, handlePx: 6,
      lane: { yU: 1.62, yD: -1.6, yD2: -2.0 }, title: [-9.6, -3.45, 'left'],
      ro: { x: 3.6, colA: 6.35, colB: 8.3, yHead: 2.85, yCols: 2.25, yA: 1.3, yM: -1.2, yChip: 0.45, midApprox: true },
      field: { x0: -9.2, x1: 2.2, rows: [1.8, 0.6, -0.6, -1.8], jit: 0.2, hy: [3.3, 2.82] },
      b2: { cx: 4.55, cy: 1.2, s: 0.4, sy: 0.36, box: 'of hard cells, out of scope', lx: 6.15, ly: [1.6, 1.1, 0.6],
        ax: 3.6, vx: 8.3, ay: [-0.6, -1.3, -2.0, -2.5] },
    },
    narrow: {
      // the slider sits under the lanes here, so the play button's top-right corner stays clear
      W: 16, H: 9, x0: -7.3, x1: 7.3, yU: 1.47, yD: 1.27, hmax: 1.3, ty: -0.5, above: false, dotPx: 3.4, handlePx: 5,
      lane: { yU: 2.12, yD: 0.5, yD2: 0.12 }, title: [7.3, -3.84, 'right'], small: true,
      ro: { x: -7.3, colA: -0.6, colB: 2.7, yHead: -1.55, yCols: -1.55, yA: -2.22, yM: -2.9, yChip: -3.82, midApprox: false },
      field: { x0: -7.0, x1: 7.0, rows: [2.15, 1.4, 0.65, -0.1], jit: 0.13, hy: [3.85, 3.32] },
      b2: { cx: -5.75, cy: -1.85, s: 0.42, sy: 0.36, box: 'out of scope', lx: -4.45, ly: [-1.45, -1.95, -2.45],
        ax: 0.6, vx: 5.3, ay: [-1.1, -1.75, -2.45, -3.05] },
    },
  };
  let L = LAYOUTS.wide, ppu = 50;
  const X = (c) => lerp(L.x0, L.x1, c);
  // rough on-screen width of a label, in world units (for spacing only)
  const est = (text, size, k = 0.56) => (text.length * size * ctx.textScale * k) / ppu;

  /* ── beat 1: two ranked lanes, the coverage slider and its curtain ── */
  const N = 40; // 40 pillars per lane: 80% and 40% cut exactly between pillars
  const rnd = ctx.rand(41);
  const curve = (p, a, b) => Array.from({ length: N }, (_, k) => {
    const r = (k + 0.5) / N;
    return clamp(a + (1 - a) * Math.pow(1 - r, p) + b * (rnd() - 0.5), 0.08, 1);
  }).sort((q, w) => w - q);
  const H = { U: curve(1.4, 0.12, 0.06), D: curve(1.1, 0.12, 0.07) };
  const lanes = [['U', 'accent', 1], ['D', 'teddy', -1]].map(([key, col, sgn]) => ({
    key, sgn,
    base: ln({ color: 'line', width: 1 }),
    bars: H[key].map((h) => ({ h, m: add(ctx.box(1, 1, { color: col, radius: 0 }), 1) })),
  }));
  const curtain = add(ctx.box(1, 1, { color: 'card', radius: 0, opacity: 0.82 }), 2);
  const cutLine = ln({ color: 'ink', width: 1.3, dashed: [4, 4] });
  /* plates keep the lane names readable when the keep cursor sweeps under them (narrow layout) */
  const plate = { bg: 'card', bgOpacity: 0.88, pad: 2 };
  const laneU = lab('ANM confidence', Object.assign({ anchor: 'right', color: 'accent', weight: 600 }, plate));
  const laneD = lab('rule’s margin', Object.assign({ anchor: 'right', color: 'teddy', weight: 600 }, plate));
  const laneD2 = lab('top − runner-up', Object.assign({ anchor: 'right', color: 'muted', size: 11 }, plate));

  const trackBg = ln({ color: 'line', width: 2 });
  const trackOn = ln({ color: 'muted', width: 3 });
  const ticks = [1, LEAD, KEEP, CMIN].map((c) => ({ c, l: ln({ color: 'faint', width: 1.2 }) }));
  const leadTick = ln({ color: 'accent', width: 2.4 });
  const minLab = lab(`${NUMBERS.covMin}%`, { anchor: 'bottom', color: 'muted', size: 11 });
  const handleBack = add(ctx.dot([0, 0], { color: 'card' }), 4);
  const handle = add(ctx.dot([0, 0], { color: 'ink' }), 4.1);
  const handleLab = lab('keep 100%', { anchor: 'bottom-right', weight: 600, size: 12 });

  /* readout: accuracy of the kept calls */
  const roHead = lab('accuracy, soft rule', { anchor: 'left', color: 'muted', size: 11 });
  const colA = lab('all cells', { anchor: 'left', color: 'muted', size: 11 });
  const colB = lab(`keep ${NUMBERS.keep}%`, { anchor: 'left', color: 'muted', size: 11 });
  const rows = [['ANM', 'accent', NUMBERS.accAnm], ['rule’s margin', 'teddy', NUMBERS.accMargin]].map(([n, col, end]) => ({
    name: lab(n, { anchor: 'left', color: col, weight: 600, size: 13 }),
    start: lab(NUMBERS.accAll, { anchor: 'left', weight: 600, size: 15 }),
    arrow: add(ctx.arrow([0, 0], [1, 0], { color: col, width: 1.8, head: 8 }), 3),
    end: lab(end, { anchor: 'left', color: col, weight: 700, size: 17 }),
    endText: end,
  }));
  const approx = lab('≈', { color: 'muted', weight: 600, size: 20 });
  const chip = lab(`${NUMBERS.lead} at ${NUMBERS.leadAt}%: ANM’s largest lead`, {
    anchor: 'left', color: 'accent', weight: 600, size: 11, bg: 'soft', bgOpacity: 1, pad: 4,
  });

  /* ── beat 2: the hard cells ── */
  const COLS = 14, ROWS = 4, ND = COLS * ROWS;
  const NH = Math.round((ND * NUMBERS.hard) / NUMBERS.cells);          // 15 of 56 (page share 4,476 / 16,750)
  const NO = Math.round((NH * NUMBERS.outOfScope) / 100);              // 11 of 15 (page share 73%)
  const rf = ctx.rand(7);
  const cells = Array.from({ length: ND }, (_, i) => ({
    i, c: i % COLS, r: Math.floor(i / COLS), jx: 2 * rf() - 1, jy: 2 * rf() - 1,
    p: new THREE.Vector3(), d: add(ctx.dot([0, 0], { color: 'faint' }), 1),
  }));
  const shuffle = (arr, r) => { for (let i = arr.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [arr[i], arr[j]] = [arr[j], arr[i]]; } return arr; };
  const order = shuffle(cells.map((c) => c.i), ctx.rand(11));
  const hardIdx = order.slice(0, NH);
  const restIdx = order.slice(NH);
  const hard = hardIdx.map((i, k) => ({
    cell: cells[i], oos: k < NO, slot: new THREE.Vector3(),
    halo: add(ctx.dot([0, 0], { color: 'bad' }), 2.5), d: add(ctx.dot([0, 0], { color: 'bad' }), 3),
  }));
  const draws = [ctx.rand(3), ctx.rand(5), ctx.rand(9)].map((r) => shuffle(restIdx.slice(), r).slice(0, NH));
  const rings = Array.from({ length: NH }, () => add(ctx.dot([0, 0], { color: 'ink', hollow: true, ring: 0.3 }), 2));
  const hardLab = lab(`${fmt(NUMBERS.hard)} hard cells`, { anchor: 'left', color: 'bad', weight: 700, size: 14 });
  const hardSub = lab('call flips under leave-one-out', { anchor: 'left', color: 'muted', size: 12 });
  const oosBox = add(ctx.box(1, 1, { color: null, stroke: 'muted', strokeWidth: 1.3, dashed: [4, 3], radius: 0.16 }), 0.5);
  const oosLab = lab('', { anchor: 'bottom-left', weight: 600, size: 12 });
  const who = [`NK/ILC ${NUMBERS.nkIlc}`, `erythroid ${NUMBERS.erythroid}`, 'over-represented'].map((t, i) =>
    lab(t, { anchor: 'left', weight: i < 2 ? 600 : 500, color: i < 2 ? 'ink' : 'muted', size: i < 2 ? 13 : 11 }));
  const accHead = lab('accuracy', { anchor: 'left', color: 'muted', size: 11 });
  const accRows = [['hard cells', NUMBERS.accHard, 'bad'], ['random cells', NUMBERS.accRandom, 'ink']].map(([t, v, col], i) => ({
    mark: add(i ? ctx.dot([0, 0], { color: 'ink', hollow: true, ring: 0.3 }) : ctx.dot([0, 0], { color: 'bad' }), 3),
    name: lab(t, { anchor: 'left', size: 12 }),
    val: lab(v, { anchor: 'left', color: col, weight: 700, size: 16 }),
  }));
  const accNote = lab(`same number, ${NUMBERS.draws} draws`, { anchor: 'left', color: 'muted', size: 11 });

  const title = lab(`${fmt(NUMBERS.cells)} held-out cells`, { anchor: 'left', color: 'muted', size: 12 });

  function layout() {
    const w = ctx.width || 1100, h = ctx.height || 420, a = w / h;
    L = LAYOUTS[a >= 2 ? 'wide' : 'narrow'];
    camera.userData.animFit = { width: L.W, height: L.H, fit: 'contain' };
    const hh = a > L.W / L.H ? L.H / 2 : L.W / 2 / a;
    camera.left = -hh * a; camera.right = hh * a; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    ppu = h / (2 * hh);
    const sp = (L.x1 - L.x0) / N, pw = 0.6 * sp;
    lanes.forEach((ln0) => {
      const y0 = ln0.sgn > 0 ? L.yU : L.yD;
      ln0.base.setPoints([[L.x0 - 0.1, y0, 0.5], [L.x1 + 0.1, y0, 0.5]]);
      ln0.bars.forEach((b, k) => {
        const hgt = L.hmax * b.h;
        b.m.scale.set(pw, hgt, 1);
        b.m.position.set(L.x0 + (k + 0.5) * sp, y0 + ln0.sgn * hgt / 2, 1);
      });
    });
    laneU.position.set(L.x1, L.lane.yU, 0);
    laneD.position.set(L.x1, L.lane.yD, 0);
    laneD2.position.set(L.x1, L.lane.yD2, 0);
    trackBg.setPoints([[L.x0, L.ty, 3], [L.x1, L.ty, 3]]);
    ticks.forEach((tk) => tk.l.setPoints([[X(tk.c), L.ty - 0.12, 3], [X(tk.c), L.ty + 0.12, 3]]));
    leadTick.setPoints([[X(LEAD), L.ty - 0.16, 3.2], [X(LEAD), L.ty + 0.16, 3.2]]);
    minLab.setAnchor(L.above ? 'bottom' : 'top').position.set(X(CMIN), L.ty + (L.above ? 0.2 : -0.2), 0);
    handleLab.setAnchor(L.above ? 'bottom-right' : 'top-right');
    handle.setPx(L.handlePx); handleBack.setPx(L.handlePx + 2);
    const R = L.ro;
    roHead.position.set(R.x, R.yHead, 0);
    colA.position.set(R.colA, R.yCols, 0);
    colB.position.set(R.colB, R.yCols, 0);
    let endR = 0, endW = 0;
    rows.forEach((r, i) => {
      const y = i ? R.yM : R.yA;
      r.name.position.set(R.x, y, 0);
      r.start.position.set(R.colA, y, 0);
      r.end.position.set(R.colB, y, 0);
      r.arrow.set([R.colA + est(NUMBERS.accAll, 15, 0.6) + 0.18, y, 0], [R.colB - 0.18, y, 0]);
      endW = Math.max(endW, est(r.endText, 17, 0.6));
      endR = R.colB + endW;
    });
    approx.position.set(R.midApprox ? R.colB + endW / 2 : endR + 0.3, (R.yA + R.yM) / 2, 0);
    chip.setSize(L.small ? 10 : 11).position.set(R.x, R.yChip, 0);
    title.setSize(L.small ? 11 : 12).setAnchor(L.title[2]).position.set(L.title[0], L.title[1], 0);
    /* beat 2 */
    const F = L.field;
    cells.forEach((c) => {
      c.p.set(lerp(F.x0, F.x1, c.c / (COLS - 1)) + F.jit * c.jx, F.rows[c.r] + F.jit * c.jy, 1);
      c.d.position.copy(c.p); c.d.setPx(L.dotPx);
    });
    rings.forEach((g) => g.setPx(L.dotPx * 1.9));
    hardLab.position.set(L.x0, F.hy[0], 0);
    hardSub.position.set(L.x0, F.hy[1], 0);
    const B = L.b2, s = B.s, sy = B.sy;
    const slots = [[-1.5, 1], [-0.5, 1], [0.5, 1], [1.5, 1], [-1, 0], [0, 0], [1, 0], [-1.5, -1], [-0.5, -1], [0.5, -1], [1.5, -1]];
    let o = 0, q = 0;
    hard.forEach((hc) => {
      if (hc.oos) { const [dx, dy] = slots[o++ % slots.length]; hc.slot.set(B.cx + dx * s, B.cy + dy * sy, 3); }
      else { hc.slot.set(B.cx + (q++ - 1.5) * s, B.cy - sy - 0.62, 3); }
      hc.d.setPx(L.dotPx); hc.halo.setPx(L.dotPx * 2.3);
    });
    const bw = 3 * s + 0.55, bh = 2 * sy + 0.52;
    oosBox.setSize(bw, bh, 0.16);
    oosBox.position.set(B.cx, B.cy, 0.5);
    oosLab.setText(`${NUMBERS.outOfScope}% ${B.box}`).position.set(B.cx - bw / 2, B.cy + bh / 2 + 0.1, 0);
    who.forEach((l, i) => l.position.set(B.lx, B.ly[i], 0));
    accHead.position.set(B.ax, B.ay[0], 0);
    accRows.forEach((r, i) => {
      const y = B.ay[1 + i];
      r.mark.position.set(B.ax + 0.15, y, 3);
      r.mark.setPx(i ? L.dotPx * 1.5 : L.dotPx);
      r.name.position.set(B.ax + 0.5, y, 0);
      r.val.position.set(B.vx, y, 0);
    });
    accNote.position.set(B.ax + 0.5, B.ay[3], 0);
  }

  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    /* visibility of the two beats (beat 1 is back in place before the loop wraps) */
    const b1 = u < 0.7 ? 1 - seg(u, 0.44, 0.475) : seg(u, 0.955, 0.995);
    const b2 = u < 0.7 ? seg(u, 0.475, 0.51) : 1 - seg(u, 0.92, 0.955);

    /* beat 1: coverage 100% → 80% (pause) → 40% */
    const c = u < 0.6 ? 1 - (1 - LEAD) * seg(u, 0.04, 0.1) - (LEAD - KEEP) * seg(u, 0.15, 0.26) : 1;
    const cx = X(c), sweep = (1 - c) / (1 - KEEP), on = clamp((1 - c) / 0.02);
    lanes.forEach((l) => { l.base.setOpacity(b1); l.bars.forEach((b) => b.m.setOpacity(b1)); });
    const top = L.yU + L.hmax + 0.1, bot = L.yD - L.hmax - 0.1, cw = L.x1 + 0.14 - cx;
    const cutTop = L.above ? L.ty - 0.1 : top, cutBot = L.above ? bot : L.ty + 0.1;
    curtain.visible = cw > 0.02 && b1 > 0;
    curtain.scale.set(Math.max(cw, 0.001), top - bot, 1);
    curtain.position.set(cx + cw / 2, (top + bot) / 2, 2);
    curtain.setOpacity(0.82 * b1);
    cutLine.setPoints([[cx, cutTop, 3], [cx, cutBot, 3]]).setOpacity(0.75 * on * b1);
    [laneU, laneD, laneD2].forEach((l) => l.setOpacity(b1));
    trackBg.setOpacity(b1);
    trackOn.setPoints([[L.x0, L.ty, 3.1], [cx, L.ty, 3.1]]).setOpacity(b1);
    ticks.forEach((tk) => tk.l.setOpacity(b1));
    const first = u < 0.6 ? 1 : 0; // the readouts reset while beat 1 is hidden, so the loop wraps cleanly
    const chipK = seg(u, 0.1, 0.125) * first;
    leadTick.setOpacity(chipK * b1);
    minLab.setOpacity(b1);
    handle.position.set(cx, L.ty, 4.1); handleBack.position.set(cx, L.ty, 4);
    handle.setOpacity(b1); handleBack.setOpacity(b1);
    handleLab.setText(`keep ${Math.round(c * 100)}%`).position.set(cx + 0.22, L.ty + (L.above ? 0.2 : -0.2), 0);
    handleLab.setOpacity(b1);
    roHead.setOpacity(b1); colA.setOpacity(b1);
    const endK = seg(u, 0.26, 0.29) * first;
    colB.setOpacity(endK * b1);
    rows.forEach((r) => {
      r.name.setOpacity(b1); r.start.setOpacity(b1);
      r.arrow.setProgress(Math.max(sweep, 0.001)).setOpacity(clamp(sweep / 0.05) * b1);
      r.end.setOpacity(endK * b1);
    });
    approx.setOpacity(seg(u, 0.29, 0.31) * first * b1);
    chip.setOpacity(chipK * b1);

    /* beat 2: hard cells glow, gather, and are named */
    const glow = seg(u, 0.51, 0.54);
    const hk = seg(u, 0.515, 0.545) * b2;
    hardLab.setOpacity(hk); hardSub.setOpacity(hk);
    const isHard = new Set(hardIdx);
    cells.forEach((cl) => cl.d.setOpacity(0.85 * b2 * (isHard.has(cl.i) ? 1 - glow : 1)));
    const breathe = 0.8 + 0.2 * Math.sin(2 * Math.PI * u * 6);
    hard.forEach((hc, i) => {
      const mk = seg(u, 0.555 + 0.0028 * i, 0.625 + 0.0028 * i);
      const p = hc.cell.p;
      const x = lerp(p.x, hc.slot.x, mk), y = lerp(p.y, hc.slot.y, mk) + 0.5 * Math.sin(Math.PI * mk);
      hc.d.position.set(x, y, 3); hc.halo.position.set(x, y, 2.5);
      hc.d.setOpacity(glow * b2);
      hc.halo.setOpacity(0.2 * glow * b2 * breathe * (1 - 0.85 * mk)); // the glow settles once gathered
    });
    const bk = seg(u, 0.63, 0.66) * b2;
    oosBox.setOpacity(bk); oosLab.setOpacity(bk);
    who.forEach((l, i) => l.setOpacity(seg(u, 0.645 + 0.02 * i, 0.675 + 0.02 * i) * b2));
    const ak = seg(u, 0.7, 0.73) * b2, rk = seg(u, 0.72, 0.75) * b2;
    accHead.setOpacity(ak);
    accRows[0].mark.setOpacity(ak); accRows[0].name.setOpacity(ak); accRows[0].val.setOpacity(ak);
    accRows[1].mark.setOpacity(rk); accRows[1].name.setOpacity(rk); accRows[1].val.setOpacity(rk);
    accNote.setOpacity(seg(u, 0.73, 0.76) * b2);
    /* random draws of the same number: rings re-draw twice (the page averages 40 draws) */
    const di = u < 0.8 ? 0 : u < 0.86 ? 1 : 2;
    const dip = Math.min(clamp(Math.abs(u - 0.8) / 0.008), clamp(Math.abs(u - 0.86) / 0.008));
    rings.forEach((g, j) => { g.position.copy(cells[draws[di][j]].p).setZ(2); g.setOpacity(0.8 * rk * dip); });
  }

  layout();
  return { scene, camera, period: PERIOD, still: 0.4 * PERIOD, update, resize: layout };
}
