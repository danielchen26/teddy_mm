/*
 * idx-exp-attr-leave-one-out — index.html, #attr (Experiment 3 · attribution, “Why this call?”), full width
 * between .plain and .bgrid. How the test answers “which marker carried this call, and how near was a flip?”
 * 1 One called cell, cite_site4_75462 (the cell the page quotes). Its 9 evidence bars (the proteins our head
 *   predicts from TEDDY's frozen embedding, each divided by its training 95th percentile and clipped to 0–1;
 *   v2 export; drawn as heights, never printed) stand under three lineage columns: the soft rule's scores,
 *   every marker weighing the same, all 9 entering at once (bridge_anm/lib/lineage_panels.py, v2 scoring).
 *   The call is the tallest column (B for this cell).
 * 2 Leave one out (bridge_anm/run_hard_proof.py, _loo_top_for_instance): each marker is dropped in turn and
 *   the cell is decided again. Dropping CD19, CD72 or CD22 puts B under T and the call flips (a dot under
 *   the marker); T and myeloid markers leave the call alone. The deciding marker is the one whose removal
 *   moves the called lineage's score most: CD22.
 * 3 The deciding marker is scanned from 0 to 1 (x) with the three lineage scores as lines (y). The call
 *   switches where B crosses T (the notch); the gap from this cell's value to the notch is the flip distance,
 *   the smallest change that switches the call (drawn, not printed).
 * 4 The camera pulls back to all 16,302 called cells (page values): how often each protein is the deciding
 *   marker (the page's chart data, drawn to scale), CD5 on top at 19.5%; a dashed ghost settles at 15.7%, the
 *   top marker's share when values are shuffled within each lineage (mean of 100 shuffles), p ≈ 0.0099;
 *   a chip: 27.5% of calls flip without their deciding marker; a note: with all markers at once the deciding
 *   marker is the largest normalised value in the called lineage.
 * Wide layout 21 × 8 when the stage is at least 2:1 (desktop 21:8), else 16 × 9 (phone 16:9).
 * 12 s seamless loop; the still frame shows steps 1–3 complete.
 */

/* Every pipeline number this scene uses. Refresh after the official-preprocessing rerun
 * (evidence: v2 export for the quoted cell; the rest: DATA.b3 in docs/index.html / reports/HARD_PROOF_V2.md). */
const NUMBERS = {
  cell: 'cite_site4_75462',          // page: “Cell cite_site4_75462, annotated CD8+ T naive”
  // v2 export (outputs/anm_cite_bridge_v2/cite_cells_meta.jsonl, adt_pred_panel / norm_p95_train, clipped 0–1).
  // Drawn as bar heights only; never printed.
  evidence: {
    CD19: 0.089, CD72: 0.168, CD22: 0.18,
    CD3: 0.156, CD2: 0.195, CD5: 0.164,
    CD16: 0.08, CD11c: 0.162, CD36: 0.15,
  },
  softBar: 0.12,                     // soft question's bar (page glossary); used to decide, not drawn
  nCalled: 16478,                    // page: “16,302 called cells”
  // page chart data (DATA.b3.top1): [protein, lineage, cells it decides]; drawn to scale, only CD5's share printed
  top1: [['CD2', 'T', 2945], ['CD5', 'T', 2849], ['CD19', 'B', 2749], ['CD36', 'M', 2500], ['CD16', 'M', 2221],
    ['CD11c', 'M', 1049], ['CD3', 'T', 1011], ['CD72', 'B', 869], ['CD22', 'B', 285]],
  permMean: 0.1573,                  // page: “vs 15.7% shuffled” (top marker's share, within-lineage shuffles)
  nPerm: 100,                        // page: “permutation 100×”
  p: 0.0099,                         // page: “p ≈ 0.0099”
  flipRate: 0.2999,                  // page: “27.5% of calls flip without their deciding marker”
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 21, height: 8 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const ln = (o = {}) => add(ctx.line([[0, 0], [1, 0]], o));
  const dot = (color, o = {}) => add(ctx.dot([0, 0, 0], Object.assign({ color }, o)));
  const rect = (color, o = {}) => add(ctx.box(1, 1, Object.assign({ color, radius: 0 }, o)));
  const pct = (x) => (100 * x).toFixed(1) + '%';
  const NCALLED = Number(NUMBERS.nCalled).toLocaleString('en-US');

  /* ── the soft rule, all 9 markers at once: S_l = sum of the lineage's values / c (all weights 1) ── */
  const LIN = [
    { id: 'b', name: 'B', markers: ['CD19', 'CD72', 'CD22'] },
    { id: 't', name: 'T', markers: ['CD3', 'CD2', 'CD5'] },
    { id: 'm', name: 'myeloid', markers: ['CD16', 'CD11c', 'CD36'] },
  ];
  const MK = LIN.flatMap((l) => l.markers);
  const LOF = MK.map((p) => LIN.findIndex((l) => l.markers.includes(p)));
  const POS = MK.map((p, i) => LIN[LOF[i]].markers.indexOf(p));
  const C = Math.max(...LIN.map((l) => l.markers.length));
  const EV0 = MK.map((p) => clamp(+NUMBERS.evidence[p] || 0, 0, 1));
  const score = (ev) => LIN.map((_, li) => ev.reduce((s, v, i) => s + (LOF[i] === li ? v : 0), 0) / C);
  const callOf = (ev) => {
    const s = score(ev);
    let b = 0;
    for (let i = 1; i < s.length; i++) if (s[i] > s[b]) b = i;
    return s[b] >= NUMBERS.softBar ? b : -1;
  };
  const S0 = score(EV0);
  const CALL = Math.max(0, callOf(EV0));

  /* leave one out: drop marker i (its event leaves the field), decide again */
  const LOO = MK.map((_, i) => {
    const drop = (k) => EV0.map((v, j) => (j === i ? v * (1 - k) : v));
    const out = callOf(drop(1));
    const move = Math.abs(score(drop(1))[CALL] - S0[CALL]);
    let kStar = null;
    if (out !== CALL) for (let s = 1; s <= 400; s++) if (callOf(drop(s / 400)) !== CALL) { kStar = s / 400; break; }
    return { out, move, kStar, flips: out !== CALL };
  });
  let DEC = 0;
  LOO.forEach((r, i) => { if (r.move > LOO[DEC].move) DEC = i; });
  const DL = LOF[DEC], V0 = EV0[DEC];

  /* scan the deciding marker from 0 to 1: where does the call switch? */
  const evAt = (v) => EV0.map((x, j) => (j === DEC ? v : x));
  const NS = 1000, calls = [];
  for (let s = 0; s <= NS; s++) calls.push(callOf(evAt(s / NS)));
  const bounds = [];
  for (let s = 1; s <= NS; s++) {
    if (calls[s] === calls[s - 1]) continue;
    let lo = (s - 1) / NS, hi = s / NS;
    for (let it = 0; it < 30; it++) { const m = (lo + hi) / 2; if (callOf(evAt(m)) === calls[s - 1]) lo = m; else hi = m; }
    bounds.push((lo + hi) / 2);
  }
  const cuts = [0, ...bounds, 1];
  const SEGS = cuts.slice(1).map((b, i) => ({ a: cuts[i], b, call: callOf(evAt((cuts[i] + b) / 2)) }));
  const VSTAR = bounds.length ? bounds.reduce((q, x) => (Math.abs(x - V0) < Math.abs(q - V0) ? x : q)) : null;
  const lineAt = (li, v) => S0[li] + (li === DL ? (v - V0) / C : 0);
  const linTok = (c) => (c >= 0 ? LIN[c].id : 'faint');
  const linName = (c) => (c >= 0 ? LIN[c].name : 'no call');

  /* all called cells: deciding-marker counts */
  const TOP = NUMBERS.top1.map(([n, l, c]) => ({ n, l, share: c / NUMBERS.nCalled }));
  const DECROW = Math.max(0, TOP.findIndex((r) => r.n === MK[DEC]));

  /* ── timeline (fractions of the 12 s loop) ── */
  const U = {
    loo0: 0.02, dCall: 0.036, dOther: 0.022,
    dec: [0.262, 0.285], fade: [0.29, 0.31], slide: [0.3, 0.325], axis: [0.298, 0.325], guide: [0.315, 0.33],
    scan: [0.335, 0.435], brk: [0.44, 0.46], out: [0.49, 0.555], back: [0.945, 1.0], reset: 0.6,
  };
  const inCubicInv = (y) => { let lo = 0, hi = 1; for (let i = 0; i < 24; i++) { const m = (lo + hi) / 2; if (ease.inOutCubic(m) < y) lo = m; else hi = m; } return lo; };
  let tw = U.loo0;
  const WIN = MK.map((_, i) => {
    const d = LOF[i] === CALL ? U.dCall : U.dOther, a = tw;
    tw += d;
    const tf = LOO[i].kStar != null ? a + 0.35 * d * inCubicInv(LOO[i].kStar) : null;
    return { a, b: a + d, d, tf };
  });

  /* ── layouts ── */
  const LAYOUTS = {
    wide: {
      W: 21, H: 8, title: 'one called cell · soft rule', tX: -10.1, tY: 3.45,
      pitch: 1.35, ggap: 0.8, bw: 0.62, yE0: -2.35, E: 5.4, lift: 0.72,
      yC0: -0.2, Hc: 11, top: 3.0, nameY: -2.47, flipY: -3.12, decY: -3.58, PW: 13.2, side: true,
      flips: 'removing it flips the call', dec: `deciding marker ${MK[DEC]}: its removal moves ${LIN[CALL].name}’s score most`,
      B: { title: `deciding marker · all ${NCALLED} called cells`, x0: -6.2, len: 35, row0: 2.45, pitch: 0.6, bh: 0.32, wide: true, rx: 4.3 },
    },
    narrow: {
      W: 16, H: 9, title: 'one cell · soft rule', tX: -7.6, tY: 3.95,
      pitch: 1.6, ggap: 0.5, bw: 0.7, yE0: -2.55, E: 6.0, lift: 0.68,
      yC0: -0.35, Hc: 11.6, top: 2.9, nameY: -2.7, flipY: -3.5, decY: -4.0, PW: 11.8, side: false,
      flips: 'flips the call', dec: `deciding marker: moves ${LIN[CALL].name} most`,
      B: { title: `all ${NCALLED} called cells`, x0: -4.6, len: 27, row0: 2.55, pitch: 0.54, bh: 0.3, wide: false },
    },
  };
  let L = LAYOUTS.wide;
  const P = { ppu: 50, XB: 22.5, Z: 1, xs: [], gx: [], cw: 3, P0: 0, PW: 13.2, xDec: 0, VEND: 1, xEnd: 0, yEnd: 0 };

  /* text width in world units (for placing a label right after another) */
  const mctx = (() => { try { return document.createElement('canvas').getContext('2d'); } catch (e) { return null; } })();
  const fam = (() => { try { return getComputedStyle(ctx.stage).fontFamily || 'sans-serif'; } catch (e) { return 'sans-serif'; } })();
  const textW = (s, size, weight = 500) => {
    const px = size * ctx.textScale;
    if (!mctx) return (s.length * px * 0.58 + 6) / P.ppu;
    mctx.font = `${weight} ${px}px ${fam}`;
    return (mctx.measureText(s).width + 6) / P.ppu;
  };
  const lhW = (size) => (size * ctx.textScale * 1.22 + 3) / P.ppu;

  /* ── screen A: one cell ── */
  const titleA = lab('', { anchor: 'left', weight: 600, size: 13 });
  const idA = lab(NUMBERS.cell, { anchor: 'left', mono: true, size: 11, color: 'muted' });
  const sideScore = lab('lineage score', { anchor: 'right', color: 'muted', size: 11 });
  const sideEv = lab('predicted proteins\n(TEDDY + head)', { anchor: 'right', color: 'teddy', size: 11, weight: 600 });
  const cols = LIN.map((l) => ({
    fill: rect(l.id, { opacity: 0.85 }),
    ghost: add(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'muted', width: 1.2, dashed: [4, 3], closed: true })),
    name: lab(l.name, { anchor: 'bottom', weight: 600, size: 12.5, color: l.id }),
    now: dot(l.id),
  }));
  const tag = lab('call', { anchor: 'bottom', bg: 'accent', bgOpacity: 1, color: 'card', weight: 600, size: 11, pad: 3 });
  const bars = MK.map((p) => ({
    fill: rect('teddy'),
    slot: add(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'muted', width: 1.2, dashed: [3, 3], closed: true })),
    name: lab(p, { anchor: 'top', size: 12 }),
    flip: dot('accent'),
  }));
  const decRing = add(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'accent', width: 2.2, closed: true }));
  const decName = lab(MK[DEC], { anchor: 'top', size: 12, weight: 700, color: 'accent' });
  const flipsLeg = lab('', { anchor: 'left', color: 'accent', size: 11.5 });
  const decLeg = lab('', { anchor: 'left', color: 'accent', size: 11.5, weight: 600 });

  // scan plot
  const axisBase = ln({ color: 'faint', width: 1.5 });
  const axisSegs = SEGS.map((s) => ({ s, line: ln({ color: linTok(s.call), width: 4.5 }), l: lab('call ' + linName(s.call), { anchor: 'top', size: 11.5, weight: 600, color: linTok(s.call) }) }));
  const zeroLab = lab('0', { anchor: 'right', color: 'muted', size: 11 });
  const oneLab = lab('1', { anchor: 'left', color: 'muted', size: 11 });
  const axisTitle = lab(MK[DEC] + ' value', { anchor: 'bottom-right', color: 'teddy', size: 11, weight: 600 });
  const guide = ln({ color: 'teddy', width: 1.3, dashed: [3, 3] });
  const v0Ring = dot('teddy', { hollow: true, ring: 0.34 });
  const plotLines = LIN.map((l, li) => (li === DL
    ? add(ctx.arrow([0, 0], [1, 0], { color: l.id, width: 2.2, head: 8 }))
    : ln({ color: l.id, width: 2.2 })));
  const plotLabs = LIN.map((l, li) => lab(l.name, { anchor: li === DL ? 'bottom-right' : 'left', weight: 600, size: 12, color: l.id }));
  const knobBack = dot('card');
  const knob = dot('teddy');
  const notch = ln({ color: 'ink', width: 2 });
  const notchUp = ln({ color: 'muted', width: 1.2, dashed: [3, 3] });
  const cross = dot('ink', { hollow: true, ring: 0.36 });
  const brk = ln({ color: 'ink', width: 1.6 });
  const brkT = [ln({ color: 'ink', width: 1.6 }), ln({ color: 'ink', width: 1.6 })];
  const brkLab = lab('flip distance', { anchor: 'left', weight: 600, size: 12 });

  /* ── screen B: all called cells ── */
  const titleB = lab('', { anchor: 'left', weight: 600, size: 13 });
  const baseB = ln({ color: 'line', width: 1.2 });
  const rows = TOP.map((r, i) => ({
    r, bar: rect('accent', { opacity: i === 0 ? 1 : 0.72 }),
    name: lab(`${r.n} · ${r.l}`, { anchor: 'right', size: 12, weight: i === 0 ? 700 : 500, color: i === 0 ? 'ink' : 'muted' }),
  }));
  const ghost = add(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'ink', width: 1.4, dashed: [4, 3], closed: true }));
  const valA = lab(pct(TOP[0].share), { anchor: 'left', weight: 700, size: 12, color: 'accent' });
  const valB = lab(`${pct(TOP[0].share)} vs ${pct(NUMBERS.permMean)} shuffled`, { anchor: 'left', weight: 700, size: 12, color: 'accent' });
  const pB = lab('', { anchor: 'left', weight: 600, size: 13 });
  const pSub = lab(`${TOP[0].n} vs ${NUMBERS.nPerm} within-lineage shuffles`, { anchor: 'left', color: 'muted', size: 11 });
  const chip = lab('', { anchor: 'left', weight: 600, size: 12.5, color: 'accent', bg: 'soft', bgOpacity: 1, pad: 6 });
  const note = lab('', { anchor: 'left', color: 'muted', size: 11.5 });
  const tracer = dot('accent');

  const rectPts = (x0, y0, x1, y1, z = 0) => [[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]];
  const xOf = (v) => P.P0 + v * P.PW;
  const yOf = (s) => L.yC0 + s * L.Hc;
  const rowY = (i) => L.B.row0 - i * L.B.pitch;

  function layout() {
    const w = ctx.width || 1100, h = ctx.height || 420, a = w / h;
    L = LAYOUTS[a >= 2 ? 'wide' : 'narrow'];
    camera.userData.animFit = { width: L.W, height: L.H, fit: 'contain' };
    const hh = a > L.W / L.H ? L.H / 2 : L.W / 2 / a, hw = hh * a;
    camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    P.ppu = h / (2 * hh);
    P.XB = L.W + 1.5;
    P.Z = Math.max(0, (2 * L.W + 1.5 + 1.2) / (2 * hw) - 1);
    const ppu = P.ppu;

    /* screen A */
    titleA.setText(L.title).position.set(L.tX, L.tY, 0);
    idA.position.set(L.tX + textW(L.title, 13, 600) + 0.35, L.tY, 0);
    const gc = 3 * L.pitch + L.ggap;
    P.gx = [-gc, 0, gc];
    P.xs = MK.map((_, i) => P.gx[LOF[i]] + (POS[i] - 1) * L.pitch);
    P.cw = 2 * L.pitch + L.bw + 0.1;
    P.xDec = P.xs[DEC];
    // the scan runs 0 → 1 with this cell's value under its bar; shorten it if its far end, the "1" and the
    // flat lineages' names after it would leave the frame (with this cell's low value it ran past the right
    // edge at every width)
    const endW = Math.max(textW('1', 11) + 0.15, ...LIN.map((l, li) => (li === DL ? 0 : textW(l.name, 12, 600) + 0.16)));
    P.PW = Math.min(L.PW, (hw - 0.1 - endW - P.xDec) / Math.max(0.05, 1 - V0));
    P.P0 = P.xDec - V0 * P.PW;
    const left = P.gx[0] - P.cw / 2;
    sideScore.position.set(left - 0.3, L.yC0 + 0.6, 0);
    sideEv.position.set(left - 0.3, L.yE0 + 0.55, 0);
    MK.forEach((_, i) => {
      const b = bars[i], x = P.xs[i], hb = EV0[i] * L.E;
      b.slot.setPoints(rectPts(x - L.bw / 2, L.yE0, x + L.bw / 2, L.yE0 + hb, 0.6));
      b.name.position.set(x, L.nameY, 0);
      b.flip.position.set(x, L.flipY, 2); b.flip.setPx(clamp(0.075 * ppu, 3, 4.5));
    });
    const hd = EV0[DEC] * L.E, pad = 4 / ppu;
    decRing.setPoints(rectPts(P.xDec - L.bw / 2 - pad, L.yE0 - pad, P.xDec + L.bw / 2 + pad, L.yE0 + hd + pad, 2));
    decName.position.set(P.xDec, L.nameY, 0);
    // after the rightmost flip dot, so no other marker's dot lands inside the words
    const xFlip = Math.max(P.xDec, ...MK.map((_, i) => (LOO[i].flips ? P.xs[i] : -Infinity)));
    flipsLeg.setText(L.flips).position.set(xFlip + 0.32, L.flipY, 0);
    decLeg.setText(L.dec).position.set(P.xDec - L.bw / 2, L.decY, 0);
    LIN.forEach((_, li) => {
      const c = cols[li], x = P.gx[li], hc = S0[li] * L.Hc;
      c.ghost.setPoints(rectPts(x - P.cw / 2, L.yC0, x + P.cw / 2, L.yC0 + hc, 0.5));
      c.now.setPx(clamp(0.085 * ppu, 3.2, 5));
    });
    // scan plot
    const P1 = P.P0 + P.PW;
    axisBase.setPoints([[P.P0, L.yC0, 1.4], [P1, L.yC0, 1.4]]);
    axisSegs.forEach((q) => {
      q.line.setPoints([[xOf(q.s.a), L.yC0, 1.5], [xOf(q.s.b), L.yC0, 1.5]]);
      q.l.position.set(xOf((q.s.a + q.s.b) / 2), L.yC0 - 0.2, 0);
    });
    zeroLab.position.set(P.P0 - 0.15, L.yC0, 0);
    oneLab.position.set(P1 + 0.15, L.yC0, 0);
    axisTitle.position.set(P1, L.yC0 + 0.14, 0);
    guide.setPoints([[P.xDec, L.yE0 + hd + 0.1, 1.3], [xOf(V0), yOf(S0[DL]), 1.3]]);
    v0Ring.position.set(xOf(V0), L.yC0, 3); v0Ring.setPx(clamp(0.12 * ppu, 5, 7));
    const Smax = (L.top - L.yC0) / L.Hc;
    const s1 = lineAt(DL, 1);
    if (s1 > Smax) { P.VEND = clamp(V0 + (Smax - S0[DL]) * C, 0, 1); P.xEnd = xOf(P.VEND); P.yEnd = L.top; }
    else { P.VEND = 1; P.xEnd = P1; P.yEnd = yOf(s1); }
    LIN.forEach((_, li) => {
      if (li === DL) plotLines[li].set([P.P0, yOf(lineAt(li, 0)), 2], [P.xEnd, P.yEnd, 2]);
      else plotLines[li].setPoints([[P.P0, yOf(S0[li]), 2], [P1, yOf(S0[li]), 2]]);
    });
    knob.setPx(clamp(0.1 * ppu, 4.5, 6)); knobBack.setPx(clamp(0.1 * ppu, 4.5, 6) + 2);
    if (VSTAR != null) {
      const xs = xOf(VSTAR), ys = yOf(lineAt(DL, VSTAR));
      notch.setPoints([[xs, L.yC0 - 0.17, 3], [xs, L.yC0 + 0.17, 3]]);
      notchUp.setPoints([[xs, L.yC0, 1.2], [xs, ys, 1.2]]);
      cross.position.set(xs, ys, 3); cross.setPx(clamp(0.1 * ppu, 4, 6));
      const yb = L.yC0 + 0.36, xa = Math.min(xs, xOf(V0)), xb = Math.max(xs, xOf(V0));
      brk.setPoints([[xa, yb, 3], [xb, yb, 3]]);
      brkT[0].setPoints([[xa, yb - 0.09, 3], [xa, yb + 0.09, 3]]);
      brkT[1].setPoints([[xb, yb - 0.09, 3], [xb, yb + 0.09, 3]]);
      brkLab.position.set(xb + 0.14, yb, 0);
    }

    /* screen B */
    const B = L.B, XB = P.XB, x0 = XB + B.x0;
    titleB.setText(B.title).position.set(XB + L.tX, L.tY, 0);
    baseB.setPoints([[x0, B.row0 + B.pitch / 2, 0.2], [x0, rowY(TOP.length - 1) - B.pitch / 2, 0.2]]);
    rows.forEach((q, i) => q.name.setSize(B.wide ? 12 : 11.5).position.set(x0 - 0.15, rowY(i), 0));
    valA.position.set(x0 + TOP[0].share * B.len + 0.15, B.row0, 0);
    valB.position.set(x0 + TOP[0].share * B.len + 0.15, B.row0, 0);
    const pText = `p ≈ ${NUMBERS.p.toFixed(4)}`;
    if (B.wide) {
      pB.setText(pText).position.set(XB + B.rx, 1.65, 0);
      pSub.position.set(XB + B.rx, 1.25, 0);
      // the right-hand column has hw − rx world units; on mid-size stages (≈ 600–800 px) the long lines ran off
      const room = hw - B.rx - 0.2;
      const pick = (vs, size, wt, extra = 0) => vs.find((v) => Math.max(...v.split('\n').map((t) => textW(t, size, wt))) + extra <= room) || vs[vs.length - 1];
      const chipT = pick([`${pct(NUMBERS.flipRate)} of calls flip\nwithout their deciding marker`, `${pct(NUMBERS.flipRate)} of calls flip\nwithout their\ndeciding marker`], 12.5, 600, 20 / ppu);
      // stack the chip under the subtitle and the note under the chip, whatever their line counts
      const subH = (11 * ctx.textScale * 1.22 + 2) / ppu, chipTop = 1.25 - subH / 2 - 0.22;
      chip.setText(chipT).setAnchor('top-left').position.set(XB + B.rx, chipTop, 0);
      const chipH = (12.5 * ctx.textScale * 1.22 * chipT.split('\n').length + 12) / ppu;
      note.setText(pick(['With all markers at once, the deciding\nmarker is the largest normalised value\nin the called lineage.', 'With all markers at once,\nthe deciding marker is the\nlargest normalised value\nin the called lineage.'], 11.5, 500))
        .setAnchor('top-left').position.set(XB + B.rx, chipTop - chipH - 0.2, 0);
    } else {
      pB.setText(`${pText} · ${NUMBERS.nPerm} shuffles`).position.set(XB + L.tX + textW(B.title, 13, 600) + 0.35, L.tY, 0);
      chip.setText(`${pct(NUMBERS.flipRate)} of calls flip without their deciding marker`).setAnchor('center').position.set(XB, -2.45, 0);
      note.setText('All markers at once: the deciding marker is\nthe largest normalised value in the called lineage.').setAnchor('center').position.set(XB, -3.65, 0);
    }
    const tp = B.wide;
    pSub.visible = tp;
    sideScore.visible = sideEv.visible = L.side;
  }

  const trA = new THREE.Vector3(), trB = new THREE.Vector3();
  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    const B = L.B, XB = P.XB;
    const live = u < U.reset ? 1 : 0;
    const aA = u < 0.5 ? 1 - seg(u, U.out[0] - 0.004, U.out[0] + 0.012) : seg(u, 0.982, 0.998);
    const aB = seg(u, U.out[1] - 0.012, U.out[1] + 0.006) * (1 - seg(u, U.back[0], U.back[0] + 0.015));

    /* camera: A → pull back → B, and back at the loop end */
    let s = 0;
    if (u >= U.back[0]) s = 1 - seg(u, U.back[0], U.back[1], ease.inOutSine);
    else if (u >= U.out[0]) s = seg(u, U.out[0], U.out[1], ease.inOutSine);
    camera.position.x = XB * s;
    camera.zoom = 1 / (1 + P.Z * Math.sin(Math.PI * s));
    camera.updateProjectionMatrix();

    /* ── screen A: leave one out ── */
    titleA.setOpacity(aA); idA.setOpacity(aA);
    sideScore.setOpacity(aA); sideEv.setOpacity(aA);
    const K = MK.map((_, i) => { const w = WIN[i]; return seg(u, w.a, w.a + 0.35 * w.d) * (1 - seg(u, w.b - 0.3 * w.d, w.b)); });
    const S = score(EV0.map((v, i) => v * (1 - K[i])));
    const colOn = u < U.fade[0] || u >= U.reset ? 1 : 1 - seg(u, U.fade[0], U.fade[1]);
    const colTop = S.map((sc) => yOf(sc));
    LIN.forEach((_, li) => {
      const c = cols[li], hc = Math.max(1e-4, S[li] * L.Hc);
      c.fill.scale.set(P.cw, hc, 1);
      c.fill.position.set(P.gx[li], L.yC0 + hc / 2, 0);
      c.fill.setOpacity(0.85 * colOn);
      let kl = 0;
      MK.forEach((_, i) => { if (LOF[i] === li) kl = Math.max(kl, K[i]); });
      c.ghost.setOpacity(0.9 * clamp(kl * 3) * colOn);
      c.name.position.set(P.gx[li], colTop[li] + 0.1, 0);
      c.name.setOpacity(aA * colOn);
    });
    let tTo = CALL, f = 0;
    MK.forEach((_, i) => { if (K[i] > 0 && LOO[i].kStar != null) { tTo = LOO[i].out; f = clamp((K[i] - LOO[i].kStar) / 0.15); } });
    const tDst = tTo >= 0 ? tTo : CALL, lh = lhW(12.5);
    tag.position.set(lerp(P.gx[CALL], P.gx[tDst], f), lerp(colTop[CALL], colTop[tDst], f) + 0.12 + lh, 0);
    tag.setOpacity(aA * colOn * (tTo >= 0 ? 1 : 1 - f));

    const dim = u >= U.fade[0] && u < U.reset ? seg(u, U.fade[0], U.fade[1]) : 0;
    let legK = 0;
    MK.forEach((_, i) => {
      const b = bars[i], x = P.xs[i], hb = Math.max(1e-4, EV0[i] * L.E), k = K[i];
      b.fill.scale.set(L.bw, hb, 1);
      b.fill.position.set(x, L.yE0 + hb / 2 + L.lift * ease.outCubic(k), 1 + k);
      b.fill.setOpacity((1 - 0.62 * k) * (i === DEC ? 1 : 1 - 0.55 * dim));
      b.slot.setOpacity(0.85 * clamp(k * 2));
      b.name.setOpacity(aA * (i === DEC ? 1 - seg(u, U.dec[0], U.dec[1]) * live : 1 - 0.45 * dim));
      const fk = WIN[i].tf != null ? seg(u, WIN[i].tf, WIN[i].tf + 0.012) * live : 0;
      b.flip.setOpacity(fk);
      legK = Math.max(legK, fk);
    });
    flipsLeg.setOpacity(aA * legK);
    const dk = seg(u, U.dec[0], U.dec[1]) * live;
    decRing.setOpacity(dk).setProgress(dk);
    decName.setOpacity(aA * dk);
    decLeg.setOpacity(aA * dk);

    /* columns → three dots at this cell's value → scan plot */
    const nowOn = seg(u, U.fade[0], U.fade[0] + 0.012) * live;
    const sk = seg(u, U.slide[0], U.slide[1]);
    LIN.forEach((_, li) => {
      cols[li].now.position.set(lerp(P.gx[li], xOf(V0), sk), yOf(S0[li]), 3);
      cols[li].now.setOpacity(nowOn);
    });
    const ak = seg(u, U.axis[0], U.axis[1]) * live;
    axisBase.setProgress(ak).setOpacity(live);
    const axL = seg(u, U.axis[1] - 0.01, U.axis[1]) * live * aA;
    zeroLab.setOpacity(axL); oneLab.setOpacity(axL); axisTitle.setOpacity(axL);
    guide.setProgress(seg(u, U.guide[0], U.guide[1])).setOpacity(0.9 * live);
    v0Ring.setOpacity(seg(u, U.guide[0], U.guide[1]) * live);
    const v = seg(u, U.scan[0], U.scan[1], ease.inOutSine);
    const scanOn = (u >= U.scan[0] ? 1 : 0) * live;
    axisSegs.forEach((q) => {
      q.line.setProgress(clamp((v - q.s.a) / Math.max(1e-6, q.s.b - q.s.a))).setOpacity(scanOn);
      const mid = q.s.a + 0.3 * (q.s.b - q.s.a);
      q.l.setOpacity(aA * scanOn * clamp((v - mid) / 0.05));
    });
    // the two flat lineages can score almost the same: spread their end names so they do not overprint
    const flat = LIN.map((_, li) => li).filter((li) => li !== DL).sort((p, q) => S0[q] - S0[p]);
    const labY = {}, nameH = (12 * ctx.textScale * 1.22 + 2) / P.ppu;
    flat.forEach((li) => { labY[li] = yOf(S0[li]); });
    for (let k = 1; k < flat.length; k++) {
      const hi = flat[k - 1], lo = flat[k], d = nameH - (labY[hi] - labY[lo]);
      if (d > 0) { labY[hi] += d / 2; labY[lo] -= d / 2; }
    }
    LIN.forEach((_, li) => {
      const pl = plotLines[li], pl2 = plotLabs[li];
      if (li === DL) {
        const vv = Math.min(v, P.VEND);
        pl.setProgress(clamp(v / Math.max(1e-6, P.VEND))).setOpacity(v > 0.002 ? scanOn : 0);
        const tipX = xOf(vv), tipY = vv >= P.VEND ? P.yEnd : yOf(lineAt(li, vv));
        pl2.position.set(tipX - 0.14, tipY + 0.1, 0);
      } else {
        pl.setProgress(v).setOpacity(scanOn);
        pl2.position.set(xOf(v) + 0.16, labY[li], 0);
      }
      pl2.setOpacity(aA * scanOn * clamp(v / 0.02));
    });
    const kOn = seg(u, U.scan[0] - 0.006, U.scan[0] + 0.004) * (1 - seg(u, U.scan[1] + 0.004, U.scan[1] + 0.014)) * live;
    const segNow = SEGS.find((q) => v <= q.b) || SEGS[SEGS.length - 1];
    knob.setColor(linTok(segNow.call));
    knob.position.set(xOf(v), L.yC0, 4.1); knobBack.position.set(xOf(v), L.yC0, 4);
    knob.setOpacity(kOn); knobBack.setOpacity(kOn);
    const nk = VSTAR != null ? clamp((v - VSTAR) / 0.04) * live * (u >= U.scan[0] ? 1 : 0) : 0;
    notch.setOpacity(nk); notchUp.setOpacity(0.9 * nk); cross.setOpacity(nk);
    const bk = VSTAR != null ? seg(u, U.brk[0], U.brk[1]) * live : 0;
    brk.setProgress(bk).setOpacity(bk);
    brkT.forEach((q) => q.setOpacity(bk));
    brkLab.setOpacity(aA * seg(u, U.brk[1] - 0.008, U.brk[1] + 0.006) * live);

    /* ── screen B: all called cells ── */
    const x0 = XB + B.x0;
    titleB.setOpacity(aB);
    baseB.setOpacity(seg(u, U.out[0] + 0.01, U.out[0] + 0.04));
    rows.forEach((q, i) => {
      const g = seg(u, 0.515 + 0.009 * i, 0.565 + 0.009 * i, ease.outCubic);
      const wv = Math.max(1e-4, q.r.share * B.len * g);
      q.bar.scale.set(wv, B.bh, 1);
      q.bar.position.set(x0 + wv / 2, rowY(i), 0.5);
      q.bar.visible = g > 0;
      q.name.setOpacity(aB * seg(u, 0.552 + 0.006 * i, 0.572 + 0.006 * i));
    });
    const gg = seg(u, 0.64, 0.68, ease.outCubic), gy = B.bh / 2 + 0.07;
    ghost.setPoints(rectPts(x0, B.row0 - gy, x0 + Math.max(1e-3, NUMBERS.permMean * B.len * gg), B.row0 + gy, 1));
    ghost.setOpacity(0.85 * seg(u, 0.635, 0.645));
    valA.setOpacity(aB * seg(u, 0.61, 0.625) * (1 - seg(u, 0.69, 0.7)));
    valB.setOpacity(aB * seg(u, 0.676, 0.69));
    pB.setOpacity(aB * seg(u, 0.688, 0.702));
    pSub.setOpacity(aB * seg(u, 0.688, 0.702));
    chip.setOpacity(aB * seg(u, 0.712, 0.726));
    note.setOpacity(aB * seg(u, 0.74, 0.755));

    /* this cell's deciding marker joins its row */
    const tk = seg(u, U.out[0] + 0.004, U.out[1] + 0.004, ease.inOutSine);
    trA.set(P.xDec, L.yE0 + EV0[DEC] * L.E + 0.12, 5);
    trB.set(x0 + 0.12, rowY(DECROW), 5);
    tracer.position.set(lerp(trA.x, trB.x, tk), lerp(trA.y, trB.y, tk) + 2.2 * Math.sin(Math.PI * tk), 5);
    tracer.setPx(clamp(0.09 * P.ppu, 4, 6));
    tracer.setOpacity(seg(u, U.out[0] - 0.004, U.out[0] + 0.006) * (1 - seg(u, U.out[1] + 0.004, U.out[1] + 0.02)));
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: 0.475 * PERIOD,
    update, resize: layout,
  };
}
