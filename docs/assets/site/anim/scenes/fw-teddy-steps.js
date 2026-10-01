/*
 * fw-teddy-steps — anm-framework.html #teddy-steps ("The same four steps on TEDDY’s workflow").
 * Four pedestals in a row: Question, Sources, Retained state, Readout. Each holds the generic icon of
 * the AI chapter (question card, events entering at their own steps, a state lattice, a gauge). One
 * pedestal at a time the icon turns into its TEDDY version, and that pedestal's top turns TEDDY orange:
 * 1 Question. The card rewrites itself as the written lineage question (B · T · myeloid, call if the
 *   best score clears the bar) and flips through its three tabs: soft, strict, B/T-priority. The nine
 *   bars under the lineage names are the marker weights each question uses, from the page glossary:
 *   soft = all 9 equal; strict = CD19, CD3, CD16 at weight 2; B/T-priority = only those three count,
 *   B ×1.5, T ×1.3, myeloid ×0.5. The weights are drawn, not printed; the bar values are not shown.
 * 2 Sources. The generic events, one per step, give way to the 9 head-predicted proteins (3 B, 3 T,
 *   3 myeloid), which drop in together onto step 0, with no staggering (export_cite_events.py: all
 *   panel events at t = 0). Our head predicts them from TEDDY's embedding; TEDDY predicts no proteins.
 * 3 Retained state. The lattice's nodes slide onto a solid ring of beads: TEDDY's embedding (the bead
 *   ring of #teddy-pipeline; the bead count is illustrative). Mode B works on the head's outputs.
 *   Behind it stand TEDDY-G's 12 layers as dotted grey slabs: Mode A (v3 E5 stopped at its gate; E5-M running).
 * 4 Readout. The gauge folds into a call chip ("call / no call"); only after the chip settles does a
 *   measured-protein stamp arrive (from below, never upstream) to grade it. No result is shown.
 * A thin arrow under the pedestals fills left to right; nothing runs back. The reset is a plain
 * cross-fade. 12 s seamless loop; the still frame shows all four TEDDY versions.
 * layout() sizes everything from the live px-per-unit and text widths (phone to desktop), sits every
 * icon on its pedestal and keeps the loader's play/pause pill corner clear.
 */

/* Numbers shown or drawn. Every one is on anm-framework.html; refresh here after the rerun. */
const NUMBERS = {
  panel: 9,         // #teddy-steps: "the 9 head-predicted proteins" (printed, and drawn as 9 spheres)
  steps: 4,         // #teddy-steps sources: finite_graph_scalar, 4 steps (drawn as ticks, not printed)
  layers: 12,       // TEDDY-G 70M: 12 transformer layers (page glossary; drawn as slabs, not printed)
  strictKey: 2,     // glossary, strict: CD19, CD3 and CD16 at weight 2 (the other six at 1)
  btB: 1.5,         // glossary, B/T-priority: only CD19, CD3, CD16 count; B ×1.5
  btT: 1.3,         //   T ×1.3
  btM: 0.5,         //   myeloid ×0.5
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12, N = NUMBERS;
  const CX = [-5.85, -1.95, 1.95, 5.85], COLW = 3.7, XL = -7.84;
  const LIN = ['b', 't', 'm'], LNAME = ['B', 'T', 'myeloid'];
  const STEPS = ['Question', 'Sources', 'Retained state', 'Readout'];
  const TABS = ['soft', 'strict', 'B/T-priority'];
  // Marker weights per question; order CD19 CD72 CD22 | CD3 CD2 CD5 | CD16 CD11c CD36 (panel order).
  const WTS = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [N.strictKey, 1, 1, N.strictKey, 1, 1, N.strictKey, 1, 1],
    [N.btB, 0, 0, N.btT, 0, 0, N.btM, 0, 0],
  ];
  const WMAX = Math.max(...WTS.flat());
  // When each pedestal turns (fractions of the loop); tabs switch at TAB; reset cross-fade at the end.
  const TURN = [0.04, 0.35, 0.47, 0.67], TAB = [0.22, 0.32];

  const add = (o, z = 0) => { o.position.z = z; scene.add(o); return o; };
  const at = (o, x, y) => { o.position.x = x; o.position.y = y; return o; };
  const lab = (s, o = {}) => add(ctx.label(s, Object.assign({ size: 11.5, color: 'ink' }, o)), 5);
  const tint = (b, a, fill) => { b.material.opacity = fill * a; if (b.outline) b.outline.setOpacity(a); };
  const col = (o, c) => { if (o.userData.c !== c) { o.userData.c = c; o.setColor(c); } };
  const circle = (x, y, r, n = 64) => Array.from({ length: n }, (_, k) => [x + r * Math.cos((k / n) * 2 * Math.PI), y + r * Math.sin((k / n) * 2 * Math.PI)]);
  const arc = (r, n = 40) => Array.from({ length: n + 1 }, (_, k) => [r * Math.cos(Math.PI * (1 - k / n)), r * Math.sin(Math.PI * (1 - k / n))]);

  /* ── pedestals and the one-way arrow ── */
  const peds = STEPS.map((name) => ({
    plaque: add(ctx.box(3.3, 0.72, { color: 'soft', stroke: 'line', strokeWidth: 1.2, radius: 0.12 }), 0.2),
    cap: add(ctx.box(3.5, 0.14, { color: 'line', radius: 0.05 }), 0.22),
    lit: add(ctx.box(3.5, 0.14, { color: 'teddy', radius: 0.05, opacity: 0 }), 0.24),
    name: lab(name, { size: 12, weight: 600 }),
  }));
  const arrowBase = add(ctx.arrow([-7.8, 0], [7.85, 0], { color: 'faint', width: 1.3, head: 7, opacity: 0.7 }), 0.1);
  const arrowLit = add(ctx.line([[-7.8, 0], [7.7, 0]], { color: 'muted', width: 1.6 }), 0.12);
  const arrowTip = add(ctx.dot([0, 0], { px: 3.2, color: 'muted', opacity: 0 }), 0.14);

  /* ── 1 Question: generic card → the lineage question with three tabs ── */
  const card = add(ctx.box(3, 2, { color: 'soft', stroke: 'line', strokeWidth: 1.3, radius: 0.14 }), 0.5);
  const cardAcc = add(ctx.box(3, 2, { color: null, stroke: 'accent', strokeWidth: 1.6, radius: 0.14 }), 0.52);
  const gBars = [0.8, 0.58, 0.7].map((f) => ({ f, l: add(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 3.5 }), 0.6) }));
  const grpLab = LNAME.map((n, j) => lab(n, { size: 11, color: LIN[j], weight: 600 }));
  const wBase = add(ctx.line([[0, 0], [1, 0]], { color: 'line', width: 1 }), 0.62);
  const wBars = WTS[0].map((_, k) => add(ctx.box(1, 1, { color: LIN[Math.floor(k / 3)], radius: 0 }), 0.64));
  const rule = lab('', { size: 11, color: 'muted' });
  const tabBox = add(ctx.box(1, 0.4, { color: 'accent', opacity: 0.14, stroke: 'accent', strokeWidth: 1.3, radius: 0.2 }), 0.5);
  const tabLab = TABS.map((s) => lab(s, { size: 11.5, color: 'muted' }));

  /* ── 2 Sources: events at their own steps → the 9 predicted proteins, all at step 0 ── */
  const axis = add(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 1.3 }), 0.3);
  const ticks = Array.from({ length: N.steps }, () => add(ctx.line([[0, 0], [0, 1]], { color: 'faint', width: 1.3 }), 0.3));
  const gEv = ticks.map((_, k) => ({ s: [0.85, 1.1, 0.75, 0.95][k % 4], d: add(ctx.dot([0, 0], { r: 0.1, color: 'muted' }), 0.8) }));
  const nine = Array.from({ length: N.panel }, (_, k) => ({
    row: Math.floor(k / 3) % 3, c: k % 3, d: add(ctx.dot([0, 0], { r: 0.1, color: LIN[Math.floor(k / 3) % 3], opacity: 0 }), 0.9),
  }));
  const rowLab = LNAME.map((n, j) => lab(n, { size: 11, color: LIN[j], weight: 600, anchor: 'left' }));
  const srcTop = lab('', { size: 11.5, weight: 600, anchor: 'bottom' });
  const srcBot = lab('', { size: 11, color: 'muted', anchor: 'top' });

  /* ── 3 Retained state: lattice → TEDDY's embedding ring; Mode A's layers dotted behind ── */
  const NB = 32;
  // layer slabs: thin dotted boxes where they have room, dotted lines (slabs seen edge-on) on phones
  const slabs = Array.from({ length: N.layers }, () => ({
    box: add(ctx.box(1, 0.1, { color: null, stroke: 'faint', strokeWidth: 1.1, dashed: [1.5, 2.5], radius: 0.03 }), 0.2),
    edge: add(ctx.line([[0, 0], [0, 1]], { color: 'faint', width: 1.6, dashed: [2, 2.5] }), 0.2),
  }));
  const links = Array.from({ length: 6 }, () => add(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 1.3 }), 0.4));
  const disc = add(ctx.dot([0, 0], { r: 1, color: 'card', opacity: 0 }), 0.6);
  const ringLine = add(ctx.line([[0, 0], [1, 0], [1, 1]], { color: 'teddy', width: 1.2, closed: true, opacity: 0 }), 0.62);
  const beads = Array.from({ length: NB }, (_, i) => ({ a: Math.PI / 2 - (i / NB) * 2 * Math.PI, d: add(ctx.dot([0, 0], { r: 0.05, color: 'teddy', opacity: 0 }), 0.8) }));
  const LAT = [];
  for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) LAT.push([c - 1, 1 - r]);
  // the 8 outer nodes slide to the bead at their own angle, the centre node to bead 10 (no crossings)
  const nodes = LAT.map(([x, y]) => ({
    p: [x, y],
    b: x === 0 && y === 0 ? 10 : Math.round(((Math.PI / 2 - Math.atan2(y, x)) / (2 * Math.PI)) * NB + NB) % NB,
    d: add(ctx.dot([0, 0], { r: 0.1, color: ctx.color('muted') }), 0.85), lat: [0, 0], ring: [0, 0],
  }));
  const taken = new Set(nodes.map((n) => n.b));
  const fillOrder = beads.map((_, i) => i).filter((i) => !taken.has(i));
  const emb = lab('', { size: 11.5, weight: 600, color: 'teddy', anchor: 'top' });
  const modeB = lab('', { size: 11, anchor: 'top' });
  const modeA = lab('Mode A: E5, E5-M', { size: 11, color: 'muted', anchor: 'bottom' });

  /* ── 4 Readout: gauge → call chip, then the measured-protein stamp ── */
  const gauge = add(ctx.group(), 1);
  const A55 = (55 * Math.PI) / 180;
  gauge.add(
    ctx.line(arc(1), { color: 'faint', width: 1.5 }),
    ctx.line(arc(0.62), { color: 'muted', width: 1.5, dashed: [4, 3] }),
    ctx.line([[-1, 0], [1, 0]], { color: 'faint', width: 1.2 }),
    ctx.line([[0, 0], [0.85 * Math.cos(A55), 0.85 * Math.sin(A55)]], { color: 'ink', width: 2 }),
    ctx.dot([0, 0, 0.1], { px: 3.5, color: 'ink' }),
  );
  const chip = add(ctx.box(2, 0.6, { color: 'accent', opacity: 0.14, stroke: 'accent', strokeWidth: 1.5, radius: 0.3 }), 1.2);
  const chipLab = lab('call / no call', { size: 11, weight: 600 });
  const stamp = add(ctx.group(), 2);
  const sParts = [
    ctx.dot([0, 0, 0], { r: 1, color: 'card' }),
    ctx.dot([0, 0, 0.01], { r: 1, hollow: true, ring: 0.13, color: 'measured' }),
    ...[-0.3, 0, 0.3].map((x) => ctx.line([[x, -0.3, 0.02], [x, 0.3, 0.02]], { color: 'measured', width: 2 })),
  ];
  stamp.add(...sParts);
  const stampLab = lab('', { size: 11, anchor: 'top' });

  /* ── layout: every size from px per unit and measured text ── */
  const mcv = document.createElement('canvas').getContext('2d');
  const tw = (s, size, weight = 500) => { // CSS px, as ctx.label rasterises it (no plate)
    const ls = String(s).split('\n');
    if (!mcv) return Math.max(...ls.map((l) => l.length)) * size * 0.58 * ctx.textScale + 3;
    mcv.font = `${weight} ${(size * ctx.textScale).toFixed(2)}px ${getComputedStyle(ctx.figure).fontFamily || 'sans-serif'}`;
    return Math.max(...ls.map((l) => mcv.measureText(l).width)) + 3;
  };
  const th = (size, s = '') => size * ctx.textScale * 1.22 * String(s).split('\n').length + 2;
  const fit = (vs, size, weight, maxPx) => vs.find((v) => tw(v, size, weight) <= maxPx) || vs[vs.length - 1];
  const inCol = (i, x, w) => clamp(x, CX[i] - COLW / 2 + w / 2, CX[i] + COLW / 2 - w / 2);

  const Q = {}, S = {}, ST = {}, R = {};
  let AY = -4.1, lay = '';
  function layout(force) {
    const ppu = ctx.ppu(), ts = ctx.textScale, key = `${ppu.toFixed(3)}/${ts.toFixed(3)}`;
    if (key === lay && !force) return;
    lay = key;
    const px = (v) => v / ppu, colPx = COLW * ppu, gap = Math.max(0.1, px(4));

    // Question: lineage names in one row, three bars under each, the rule; tabs listed under the card
    Q.padX = Math.max(0.18, px(5)); Q.padY = Math.max(0.14, px(6));
    Q.bw = Math.max(0.12, px(3)); Q.pb = Math.max(0.24, px(5.5)); Q.hb = Math.max(0.9, px(18));
    const grpW = 2 * Q.pb + Q.bw, wl = LNAME.map((n) => px(tw(n, 11, 600)));
    const gl = Math.max(0.16, px(7)), gmin = grpW + Math.max(0.3, px(5));
    const cT = Math.max(gmin, wl[0] / 2 + gl + wl[1] / 2), cM = cT + Math.max(gmin, wl[1] / 2 + gl + wl[2] / 2);
    const iL = Math.min(-grpW / 2, -wl[0] / 2), iR = cM + Math.max(grpW / 2, wl[2] / 2);
    const QMAX = CX[0] + COLW / 2 - XL;
    Q.rule = fit(['call if best ≥ bar', 'best ≥ bar'], 11, 500, (QMAX - 2 * Q.padX) * ppu);
    Q.cw = Math.min(QMAX, Math.max(3.3, Math.max(iR - iL, px(tw(Q.rule, 11))) + 2 * Q.padX));
    Q.x = Math.min(CX[0], CX[0] + COLW / 2 - Q.cw / 2);
    Q.gx = [0, cT, cM].map((c) => Q.x + c - (iL + iR) / 2);
    Q.hL = px(th(11)); Q.hR = px(th(11)); Q.gap = gap;
    Q.ch = 2 * Q.padY + Q.hL + gap + Q.hb + gap + Q.hR;
    Q.hr = Math.max(0.46, px(th(11.5) + 5)); Q.gapT = Math.max(0.18, px(6));
    Q.tabW = Math.max(...TABS.map((s) => px(tw(s, 11.5)))) + Math.max(0.36, px(14));
    const qRel = [Q.ch / 2, -Q.ch / 2 - Q.gapT - 3 * Q.hr];

    // Sources: 3 × 3 cluster over tick 0, lineage names to its right, step axis under it
    S.rs = Math.max(0.18, px(4.2)); S.ps = Math.max(0.52, px(11.5));
    S.gapL = Math.max(0.18, px(6));
    const totW = 2 * S.ps + 2 * S.rs + S.gapL + px(tw('myeloid', 11, 600));
    S.left = CX[1] - totW / 2; S.t0 = S.left + S.rs + S.ps;
    S.top = fit([`${N.panel} head-predicted proteins`, `${N.panel} head-predicted\nproteins`], 11.5, 600, colPx);
    S.bot = fit(['all entering at once', 'all at once'], 11, 500, colPx);
    S.hTop = px(th(11.5, S.top)); S.hBot = px(th(11, S.bot));
    S.tickH = Math.max(0.18, px(6)); S.gapA = Math.max(0.18, px(6)); S.gap = gap;
    S.axL = S.left - 0.05; S.axR = CX[1] + totW / 2 + 0.05;
    S.drop = Math.max(1.1, px(28));
    const sRel = [S.ps + S.rs + gap + S.hTop, -(S.ps + S.rs) - S.gapA - S.tickH / 2 - gap - S.hBot];

    // Retained state: ring (with a card-coloured disc) in front of the dotted layer stack
    ST.RR = Math.max(0.95, px(20)); ST.dp = Math.max(0.1, px(4));
    ST.bpx = Math.max(1.4, 0.05 * ppu); ST.npx = Math.max(2.6, 0.1 * ppu);
    ST.ls = Math.min(ST.RR * 0.7, Math.max(0.5, px(12)));
    // TEDDY-G's layers as thin upright slabs, layer 1 → 12 left to right, standing behind the ring
    ST.pitch = Math.max(0.2, px(6)); ST.sw = Math.max(0.09, px(2.6)); ST.thin = ST.sw * ppu < 5; ST.sh = 2.4 * ST.RR;
    ST.s0 = -0.55 * ST.RR; ST.sy = 0.35 * ST.RR; ST.sx = ST.s0 + ((N.layers - 1) * ST.pitch) / 2;
    ST.emb = fit(['TEDDY’s embedding', 'TEDDY’s\nembedding'], 11.5, 600, colPx);
    ST.mb = fit(['Mode B: works on the head’s outputs', 'Mode B: works on\nthe head’s outputs', 'Mode B: works on\nthe head’s\noutputs'], 11, 500, colPx);
    ST.hEmb = px(th(11.5, ST.emb)); ST.hMb = px(th(11, ST.mb)); ST.gap = gap;
    const relL = -ST.RR - ST.dp, relR = Math.max(ST.RR + ST.dp, ST.s0 + (N.layers - 1) * ST.pitch + ST.sw / 2);
    ST.x = CX[2] - (relL + relR) / 2;
    const stRel = [Math.max(ST.RR + ST.dp, ST.sy + ST.sh / 2) + gap + px(th(11)), -(ST.RR + ST.dp) - gap - ST.hEmb - ST.hMb];

    // Readout: chip; the stamp lands on its lower-right edge; the grade line under the stamp
    R.cw = Math.min(COLW, Math.max(2.4, px(tw('call / no call', 11, 600) + 16)));
    R.ch = Math.max(0.7, px(th(11) + 8));
    R.sr = Math.max(0.42, px(12)); R.Rg = Math.max(1.05, px(24)); R.gap = gap;
    R.lab = fit(['graded by measured protein', 'graded by\nmeasured protein'], 11, 500, colPx);
    const rRel = [Math.max(R.ch / 2, 0.7 * R.Rg), -R.ch / 2 - 1.75 * R.sr - gap - px(th(11, R.lab))];

    // vertical: arrow, plaque, cap, then each icon sits on its pedestal; the block is centred and
    // pushed down only as far as needed to keep the loader's pill corner (top right) clear
    const PH = Math.max(0.72, px(th(12) + 8)), SH = Math.max(0.12, px(4)), aGap = Math.max(0.24, px(8));
    const gapI = Math.max(0.3, px(10)), mB = Math.max(0.28, px(9));
    const rels = [qRel, sRel, stRel, rRel], H = rels.map(([t, b]) => t - b);
    const tops = [0, 1, 2].map(() => 4.5 - Math.max(0.2, px(8))).concat(4.5 - px(44));
    const lift = aGap + PH + SH + gapI; // arrow line → bottom of the icons
    const hi = Math.min(...H.map((h, i) => tops[i] - h)) - lift;
    AY = Math.max(-4.5 + mB, Math.min(-(lift + Math.max(...H)) / 2 - 0.15, hi));
    const PY = AY + aGap + PH / 2, yLo = AY + lift;
    const Y = rels.map(([, b]) => yLo - b);
    const PW = Math.min(COLW - 0.1, Math.max(3.3, px(tw('Retained state', 12, 600) + 14)));
    peds.forEach((p, i) => {
      p.plaque.setSize(PW, PH, Math.min(0.14, PH / 2)); at(p.plaque, CX[i], PY);
      [p.cap, p.lit].forEach((b) => { b.setSize(PW + 0.16, SH, SH / 2); at(b, CX[i], PY + PH / 2 + SH / 2); });
      at(p.name, CX[i], PY);
    });
    arrowBase.set([-7.8, AY], [7.85, AY]);
    arrowLit.setPoints([[-7.8, AY], [7.7, AY]]);

    // place Question
    const yq = Y[0], qt = yq + Q.ch / 2;
    card.setSize(Q.cw, Q.ch, 0.14); at(card, Q.x, yq);
    cardAcc.setSize(Q.cw, Q.ch, 0.14); at(cardAcc, Q.x, yq);
    Q.yL = qt - Q.padY - Q.hL / 2; Q.base = qt - Q.padY - Q.hL - gap - Q.hb; Q.yR = Q.base - gap - Q.hR / 2;
    grpLab.forEach((l, j) => at(l, Q.gx[j], Q.yL));
    wBase.setPoints([[Q.gx[0] - grpW / 2 - 0.06, Q.base], [Q.gx[2] + grpW / 2 + 0.06, Q.base]]);
    Q.bx = WTS[0].map((_, k) => Q.gx[Math.floor(k / 3)] + ((k % 3) - 1) * Q.pb);
    Q.unit = Q.hb / WMAX;
    at(rule, Q.x, Q.yR);
    const gx0 = Q.x - Q.cw / 2 + Q.padX, giw = Q.cw - 2 * Q.padX, gy = [Q.yL, Q.base + Q.hb / 2, Q.yR];
    gBars.forEach((g, k) => g.l.setPoints([[gx0, gy[k]], [gx0 + g.f * giw, gy[k]]]));
    Q.rows = [0, 1, 2].map((k) => yq - Q.ch / 2 - Q.gapT - Q.hr * (k + 0.5));
    tabLab.forEach((l, k) => at(l, Q.x, Q.rows[k]));
    const tbh = Q.hr - Math.max(0.06, px(3));
    tabBox.setSize(Q.tabW, tbh, tbh / 2);

    // place Sources
    const ys = Y[1];
    S.rowY = [ys + S.ps, ys, ys - S.ps];
    S.ay = ys - S.ps - S.rs - S.gapA;
    axis.setPoints([[S.axL, S.ay], [S.axR, S.ay]]);
    S.tx = ticks.map((_, k) => S.t0 + (k * (S.axR - 0.15 - S.t0)) / Math.max(1, N.steps - 1));
    ticks.forEach((l, k) => l.setPoints([[S.tx[k], S.ay - S.tickH / 2], [S.tx[k], S.ay + S.tickH / 2]]));
    gEv.forEach((e, k) => at(e.d, S.tx[k], ys));
    rowLab.forEach((l, j) => at(l, S.t0 + S.ps + S.rs + S.gapL, S.rowY[j]));
    srcTop.setText(S.top); at(srcTop, CX[1], S.rowY[0] + S.rs + gap);
    srcBot.setText(S.bot); at(srcBot, CX[1], S.ay - S.tickH / 2 - gap);

    // place Retained state
    const yst = Y[2];
    slabs.forEach((b, k) => {
      const x = ST.x + ST.s0 + k * ST.pitch;
      b.box.setSize(ST.sw, ST.sh, Math.min(0.03, ST.sw / 2)); at(b.box, x, yst + ST.sy);
      b.edge.setPoints([[x, yst + ST.sy - ST.sh / 2], [x, yst + ST.sy + ST.sh / 2]]);
      b.box.visible = !ST.thin; b.edge.visible = ST.thin;
    });
    disc.setRadius(ST.RR + ST.dp); at(disc, ST.x, yst);
    ringLine.setPoints(circle(ST.x, yst, ST.RR));
    beads.forEach((b) => { at(b.d, ST.x + ST.RR * Math.cos(b.a), yst + ST.RR * Math.sin(b.a)); b.d.setPx(ST.bpx); });
    nodes.forEach((n) => {
      n.lat = [ST.x + n.p[0] * ST.ls, yst + n.p[1] * ST.ls];
      n.ring = [beads[n.b].d.position.x, beads[n.b].d.position.y];
    });
    const L3 = (r, c) => nodes[r * 3 + c].lat;
    [0, 1, 2].forEach((k) => {
      links[k].setPoints([L3(k, 0), L3(k, 1), L3(k, 2)]);
      links[k + 3].setPoints([L3(0, k), L3(1, k), L3(2, k)]);
    });
    at(modeA, inCol(2, ST.x + ST.sx, px(tw('Mode A: E5, E5-M', 11))), yst + ST.sy + ST.sh / 2 + gap);
    const ey = yst - ST.RR - ST.dp - gap;
    emb.setText(ST.emb); at(emb, inCol(2, ST.x, px(tw(ST.emb, 11.5, 600))), ey);
    modeB.setText(ST.mb); at(modeB, inCol(2, ST.x, px(tw(ST.mb, 11))), ey - ST.hEmb);

    // place Readout
    const yr = Y[3];
    chip.setSize(R.cw, R.ch, R.ch / 2); at(chip, CX[3], yr); at(chipLab, CX[3], yr);
    gauge.position.set(CX[3], yr - 0.3 * R.Rg, 1);
    R.sx = CX[3] + R.cw * 0.3; R.sy = yr - R.ch / 2 - 0.75 * R.sr;
    stampLab.setText(R.lab); at(stampLab, inCol(3, R.sx, px(tw(R.lab, 11))), R.sy - R.sr - gap);
  }
  let alive = true;
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (alive) { layout(true); ctx.requestRender(); } });

  const MARKS = [[0.005, 0.04, -7.8, CX[0]], [0.315, 0.35, CX[0], CX[1]], [0.435, 0.47, CX[1], CX[2]], [0.635, 0.67, CX[2], CX[3]], [0.845, 0.88, CX[3], 7.7]];

  return {
    scene, camera, period: PERIOD, still: 0.89 * PERIOD,
    resize() { layout(true); },
    update(t) {
      layout();
      const u = ctx.loopT(t, PERIOD), c = ctx.colors;
      const out = seg(u, 0.92, 0.95, ease.inOutSine), gin = seg(u, 0.95, 0.995, ease.inOutSine);
      const keep = 1 - out, late = u >= 0.95;
      const G = (a, b) => (late ? gin : 1 - seg(u, a, b, ease.inOutSine)); // generic icon visibility
      const I = (a, b, e) => seg(u, a, b, e || ease.inOutSine) * keep;     // TEDDY piece visibility

      /* pedestal tops turn TEDDY orange as their icon turns; the arrow only ever grows rightwards */
      TURN.forEach((a, i) => peds[i].lit.setOpacity(I(a + 0.03, a + 0.06)));
      let xe = -7.8, tip = 0;
      for (const [a, b, x0, x1] of MARKS) {
        if (u < a) break;
        xe = lerp(x0, x1, seg(u, a, b, ease.inOutSine));
        tip = pulse(u, a, b, 0.008);
      }
      arrowLit.setProgress((xe + 7.8) / 15.5).setOpacity(keep);
      at(arrowTip, xe, AY).setOpacity(tip * keep);

      /* 1 Question: the card rewrites itself, then the tabs flip soft → strict → B/T-priority */
      gBars.forEach((g, k) => {
        g.l.setProgress(late ? 1 : 1 - seg(u, 0.04 + 0.012 * k, 0.07 + 0.012 * k, ease.inOutSine));
        g.l.setOpacity(late ? gin : 1);
      });
      cardAcc.outline.setOpacity(I(0.05, 0.09));
      grpLab.forEach((l, j) => l.setOpacity(I(0.06 + 0.012 * j, 0.075 + 0.012 * j)));
      wBase.setOpacity(I(0.065, 0.09));
      const grow = seg(u, 0.075, 0.11, ease.outCubic);
      const sa = seg(u, TAB[0], TAB[0] + 0.02), sb = seg(u, TAB[1], TAB[1] + 0.02), act = sa + sb;
      wBars.forEach((b, k) => {
        const h = Q.unit * lerp(lerp(WTS[0][k], WTS[1][k], sa), WTS[2][k], sb) * grow;
        b.visible = h > 1e-3;
        b.scale.set(Q.bw, Math.max(h, 1e-3), 1);
        at(b, Q.bx[k], Q.base + h / 2);
        b.setOpacity(keep);
      });
      rule.setText(Q.rule.slice(0, Math.round(seg(u, 0.085, 0.12, ease.linear) * Q.rule.length))).setOpacity(keep);
      const tv = I(0.105, 0.135);
      tabLab.forEach((l, k) => { col(l, Math.abs(act - k) < 0.5 ? 'ink' : 'muted'); l.setOpacity(tv); });
      const ti = Math.min(1, Math.floor(act)), tf = act - ti;
      at(tabBox, Q.x, lerp(Q.rows[ti], Q.rows[Math.min(2, ti + 1)], ease.inOutSine(tf)));
      tint(tabBox, tv, 0.14);

      /* 2 Sources: generic events leave; nine predicted proteins drop in together onto step 0 */
      const gs = G(0.35, 0.375);
      gEv.forEach((e) => { e.d.setRadius(S.rs * e.s * Math.max(late ? 1 : gs, 1e-3)); e.d.setOpacity(gs); });
      const drop = seg(u, 0.38, 0.425, ease.outCubic), nv = seg(u, 0.38, 0.392) * keep;
      nine.forEach((n) => {
        at(n.d, S.t0 + (n.c - 1) * S.ps, S.rowY[n.row] + (1 - drop) * S.drop);
        n.d.setRadius(S.rs).setOpacity(nv);
      });
      const fl = pulse(u, 0.42, 0.47, 0.015) * keep;
      ticks[0].setWidth(1.3 + 1.7 * fl).setColor(fl > 0.3 ? 'ink' : 'faint');
      [srcTop, srcBot, ...rowLab].forEach((l) => l.setOpacity(I(0.415, 0.445)));

      /* 3 Retained state: lattice nodes slide onto the ring; the ring fills; Mode A's stack appears */
      links.forEach((l) => l.setOpacity(G(0.47, 0.5)));
      const mv = late ? 0 : seg(u, 0.48, 0.54, ease.inOutCubic);
      nodes.forEach((n) => {
        at(n.d, lerp(n.lat[0], n.ring[0], mv), lerp(n.lat[1], n.ring[1], mv));
        n.d.material.color.copy(c.muted).lerp(c.teddy, mv);
        n.d.setPx(lerp(ST.npx, ST.bpx, mv)).setOpacity(late ? gin : keep);
      });
      taken.forEach((i) => beads[i].d.setOpacity(0));
      fillOrder.forEach((i, k) => {
        const a = 0.53 + (0.04 * k) / fillOrder.length;
        beads[i].d.setOpacity(I(a, a + 0.012));
      });
      disc.setOpacity(I(0.5, 0.53));
      ringLine.setOpacity(0.55 * I(0.52, 0.56));
      emb.setOpacity(I(0.555, 0.585));
      modeB.setOpacity(I(0.57, 0.6));
      slabs.forEach((b, k) => { const o = I(0.6 + 0.002 * k, 0.625 + 0.002 * k); b.box.setOpacity(o); b.edge.setOpacity(o); });
      modeA.setOpacity(I(0.62, 0.65));

      /* 4 Readout: gauge folds away, the chip pops and settles, then the stamp comes down */
      const gg = G(0.67, 0.7);
      gauge.scale.setScalar(R.Rg * (late ? 1 : lerp(1, 0.45, seg(u, 0.67, 0.7))));
      gauge.children.forEach((o) => o.setOpacity(gg));
      const pop = lerp(0.3, 1, seg(u, 0.69, 0.74, ease.outBack));
      chip.scale.set(pop, pop, 1);
      tint(chip, I(0.69, 0.705), 0.14);
      chipLab.setOpacity(I(0.715, 0.74));
      // the stamp comes in from below (measured protein is never upstream) and never covers the chip text
      const fall = seg(u, 0.765, 0.8, ease.outCubic), press = pulse(u, 0.8, 0.83, 0.012);
      stamp.position.set(R.sx + (1 - fall) * R.sr * 0.3, R.sy - (1 - fall) * R.sr * 0.9, 2);
      stamp.scale.setScalar(R.sr * lerp(1.2, 1, fall) * (1 - 0.06 * press));
      const sv = seg(u, 0.765, 0.78) * keep;
      sParts.forEach((o) => o.setOpacity(sv));
      stampLab.setOpacity(I(0.815, 0.845));
    },
    dispose() { alive = false; },
  };
}
