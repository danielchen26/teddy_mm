/*
 * idx-brings-audit-receipt — index.html, #brings, card “2 · Audit” (“Every call carries an audit record”).
 * ANM prints the audit record of one cell, cite_site4_75462 (annotated CD8+ T naive), line by line:
 * the question, the evidence (the 9 proteins our head predicts from TEDDY's frozen embedding, each divided
 * by its training 95th percentile and clipped to 0–1; drawn as unlabelled bars, no value printed) and the
 * three lineage scores against the question's bar. Scores are computed here with the page's written rules,
 * all 9 markers at once (bridge_anm/lib/lineage_panels.py, v2 scoring):
 *   soft          every marker weighs the same (lineage mean), bar 0.12
 *   strict        weighted mean, key markers CD19 / CD3 / CD16 at weight 2, bar 0.28
 *   B/T-priority  only CD19 / CD3 / CD16 count, B ×1.5, T ×1.3, myeloid ×0.5, bar 0.1333 (shared scale)
 * The call is the best lineage if it reaches the bar, else “no call”. The answer key for this cell is T.
 * 1 Soft: the B column edges past T (B · wrong). A handle drags CD22 down 0.09; B falls under T and the
 *   call flips to T; the record prints “CD22 falls 0.09 → T”. CD22 springs back.
 * 2 Strict: the bar rises; every column stays under it; a grey “no call” stamp lands.
 * 3 B/T-priority: secondary markers drop out; T wins by a hair (T · right). A nudge of CD3 by 0.003, too
 *   small to see, tips the call to B; a “fragile” stamp lands; the record prints “CD3 moves 0.003 → B”.
 * The three record rows end as the table under the figure. The receipt tears off. 12 s seamless loop;
 * the still frame shows the finished record under the B/T-priority question.
 */

/* Every pipeline number this scene uses. Refresh after the official-preprocessing rerun
 * (the deck's examples step writes the same fields; flips and bars must match the page's table/glossary). */
const NUMBERS = {
  cell: 'cite_site4_75462',          // page: “Cell cite_site4_75462, annotated CD8+ T naive”
  cellType: 'CD8+ T naive',
  key: 't',                          // answer key for this cell under all three questions (page: B · wrong, T · right)
  // v2 export (outputs/anm_cite_bridge_v2/cite_cells_meta.jsonl, adt_pred_panel / norm_p95_train, clipped 0–1).
  // Drawn as bar heights only; never printed.
  evidence: {
    CD19: 0.089, CD72: 0.168, CD22: 0.18,
    CD3: 0.156, CD2: 0.195, CD5: 0.164,
    CD16: 0.08, CD11c: 0.162, CD36: 0.15,
  },
  // question bars on the shared score scale (page glossary: 0.12, 0.28, 0.1333)
  bars: { soft: 0.12, strict: 0.28, btp: 0.1333 },
  // page table “Flips if”: soft “CD22 falls 0.09 → T”, B/T-priority “CD3 moves 0.003 → B”
  flips: {
    soft: { marker: 'CD2', verb: 'falls', by: 0.078 },
    btp: { marker: 'CD3', verb: 'falls', by: 0.003 },
  },
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 12, height: 9 });
  const PERIOD = 12;

  /* ── the written questions (declared rules, as on the page) ── */
  const LIN = [
    { id: 'b', name: 'B', markers: ['CD19', 'CD72', 'CD22'] },
    { id: 't', name: 'T', markers: ['CD3', 'CD2', 'CD5'] },
    { id: 'm', name: 'myeloid', markers: ['CD16', 'CD11c', 'CD36'] },
  ];
  const MK = LIN.flatMap((l) => l.markers);
  const LOF = MK.map((p) => LIN.findIndex((l) => l.markers.includes(p)));
  const QS = [
    { id: 'soft', name: 'soft', key: 1, sec: 1, lw: [1, 1, 1] },
    { id: 'strict', name: 'strict', key: 2, sec: 1, lw: [1, 1, 1] },
    { id: 'btp', name: 'B/T-priority', key: 3, sec: 0, lw: [1.5, 1.3, 0.5] },
  ];
  QS.forEach((q) => {
    const raw = MK.map((p, i) => (LIN[LOF[i]].markers[0] === p ? q.key : q.sec) * q.lw[LOF[i]]);
    const mx = Math.max(...raw);
    q.w = raw.map((r) => r / mx);
    q.c = Math.max(...LIN.map((_, li) => q.w.reduce((s, w, i) => s + (LOF[i] === li ? w : 0), 0)));
    q.bar = NUMBERS.bars[q.id];
  });
  const score = (q, ev) => LIN.map((_, li) => ev.reduce((s, v, i) => s + (LOF[i] === li ? q.w[i] * v : 0), 0) / q.c);
  const callOf = (q, ev) => {
    const s = score(q, ev);
    let b = 0;
    for (let i = 1; i < s.length; i++) if (s[i] > s[b]) b = i;
    return s[b] >= q.bar ? b : -1;
  };
  const EV0 = MK.map((p) => clamp(+NUMBERS.evidence[p] || 0, 0, 1));
  const VMAX = Math.max(0.3, ...EV0, ...QS.map((q) => q.bar));

  /* ── timeline (fractions of the loop) ── */
  const U = {
    print: [0, 0.07], grow: [0.05, 0.1], cols: [0.075, 0.12], tear: [0.93, 0.985],
    sw: [[0.37, 0.415], [0.57, 0.615]],
  };
  const DRAGS = [
    { q: 0, f: NUMBERS.flips.soft, show: [0.165, 0.185], go: [0.19, 0.27], back: [0.33, 0.36], hide: [0.355, 0.375] },
    { q: 2, f: NUMBERS.flips.btp, show: [0.665, 0.685], go: [0.69, 0.74], back: [0.8, 0.83], hide: [0.825, 0.845] },
  ];
  DRAGS.forEach((d) => { d.i = MK.indexOf(d.f.marker); });
  const CALLWIN = [[0.125, 0.37, 0], [0.43, 0.57, 1], [0.63, 0.93, 2]];
  const ROWT = [{ call: [0.135, 0.16], flip: [0.272, 0.295] }, { call: [0.445, 0.47], flip: [0.445, 0.47] }, { call: [0.645, 0.67], flip: [0.75, 0.775] }];
  const STAMPS = [{ q: 1, text: 'no call', color: 'muted', in: [0.44, 0.46], out: [0.57, 0.59] }, { q: 2, text: 'fragile', color: 'warn', in: [0.745, 0.765], out: null }];

  const evAt = (u) => {
    const ev = EV0.slice();
    for (const d of DRAGS) {
      if (d.i < 0) continue;
      const k = seg(u, d.go[0], d.go[1], ease.inOutSine) - seg(u, d.back[0], d.back[1], ease.inOutSine);
      ev[d.i] = Math.max(0, ev[d.i] - d.f.by * k);
    }
    return ev;
  };
  const qState = (u) => {
    const [s1, s2] = U.sw;
    if (u < s1[0]) return [0, 0, 0];
    if (u < s1[1]) return [0, 1, seg(u, s1[0], s1[1])];
    if (u < s2[0]) return [1, 1, 0];
    if (u < s2[1]) return [1, 2, seg(u, s2[0], s2[1])];
    return [2, 2, 0];
  };
  const qWeight = (u) => [
    seg(u, 0.07, 0.09) * (1 - seg(u, 0.37, 0.39)),
    seg(u, 0.395, 0.415) * (1 - seg(u, 0.57, 0.59)),
    seg(u, 0.595, 0.615),
  ];
  // Smoothed call indicator: share of a few recent samples that call each lineage (deterministic in u).
  const callW = (u) => {
    const w = [0, 0, 0], N = 6, du = 0.004;
    for (let k = 0; k < N; k++) {
      const uu = u - k * du;
      const win = CALLWIN.find((c) => uu >= c[0] && uu < c[1]);
      if (!win) continue;
      const c = callOf(QS[win[2]], evAt(uu));
      if (c >= 0) w[c] += 1 / N;
    }
    return w;
  };

  /* ── record texts, from the rules above ── */
  const callTxt = (c) => (c < 0 ? 'no call' : `${LIN[c].name} · ${LIN[c].id === NUMBERS.key ? 'right' : 'wrong'}`);
  const callCol = (c) => (c < 0 ? 'muted' : LIN[c].id === NUMBERS.key ? 'good' : 'bad');
  const REC = QS.map((q, qi) => {
    const c0 = callOf(q, EV0);
    const d = DRAGS.find((x) => x.q === qi && x.i >= 0);
    let flip = '–';
    if (d && c0 >= 0) {
      const ev = EV0.slice();
      ev[d.i] = Math.max(0, ev[d.i] - d.f.by);
      const c1 = callOf(q, ev);
      flip = `${d.f.marker} ${d.f.verb} ${d.f.by} → ${c1 < 0 ? 'no call' : LIN[c1].name}`;
    }
    return { name: q.name, call: callTxt(c0), col: callCol(c0), flip };
  });

  /* ── layout constants (world units, 12 × 9 frame) ── */
  const PW = 11, TOP = 3.9, BOT = -3.8, TH = 0.13, NT = 26;
  const XL = -5.1;
  const Y = { anm: 4.2, head: 3.42, q: 2.84, div1: 2.5, ch: 2.16, base: -1.0, glab: -1.28, div2: -1.56, rh: -1.86, rows: [-2.38, -2.9, -3.42] };
  const CH = 2.6, SC = CH / VMAX;
  const BW = 0.34, BSTEP = 0.48, GGAP = 0.42, GW = BW + 2 * BSTEP;
  const BX = MK.map((_, i) => XL + BW / 2 + Math.floor(i / 3) * (GW + GGAP) + (i % 3) * BSTEP);
  const GX = [0, 1, 2].map((g) => BX[g * 3 + 1]);
  const CX = [0.55, 1.55, 2.55], CW = 0.6, LX0 = 0.1, LX1 = 3.0, STX = 1.55, ROT = -0.07;

  const rc = new THREE.Group();
  scene.add(rc);
  const add = (o, parent = rc) => { parent.add(o); return o; };
  const lab = (text, o = {}, parent) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)), parent);
  const ln = (pts, o = {}, parent) => add(ctx.line(pts, o), parent);
  const rectGeo = ctx.track(new THREE.PlaneGeometry(1, 1).translate(0, 0.5, 0));
  const rect = (token, z = 1, order = 2) => {
    const m = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false })), token);
    const mesh = add(new THREE.Mesh(rectGeo, m));
    mesh.position.z = z; mesh.renderOrder = order; mesh.scale.set(1, 1e-4, 1);
    mesh.setOpacity = (a) => { m.opacity = a; return mesh; };
    return mesh;
  };

  /* printer: ANM label + slot (fixed; the receipt hangs from it) */
  const anmL = lab('ANM', { anchor: 'left', color: 'accent', weight: 700, size: 12.5 }, scene);
  const recL = lab('audit record', { anchor: 'left', color: 'muted', size: 12 }, scene);
  const slot = ln([[-PW / 2 - 0.3, TOP, 5], [PW / 2 + 0.3, TOP, 5]], { color: 'ink', width: 4, order: 6 }, scene);

  /* paper: body (top-anchored plane), zigzag teeth, outline */
  const paperMat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false })), 'soft');
  const paper = add(new THREE.Mesh(ctx.track(new THREE.PlaneGeometry(1, 1).translate(0, -0.5, 0)), paperMat));
  paper.position.set(0, TOP, 0);
  const teethShape = new THREE.Shape();
  teethShape.moveTo(-PW / 2, 0);
  for (let j = 1; j <= 2 * NT; j++) teethShape.lineTo(-PW / 2 + (j * PW) / (2 * NT), j % 2 ? -TH : 0);
  teethShape.lineTo(-PW / 2, 0);
  const teethMat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, side: THREE.DoubleSide })), 'soft');
  const teeth = add(new THREE.Mesh(ctx.track(new THREE.ShapeGeometry(teethShape)), teethMat));
  const outlinePts = (bot) => {
    const p = [[-PW / 2, TOP, 0.1], [-PW / 2, bot, 0.1]];
    for (let j = 1; j <= 2 * NT; j++) p.push([-PW / 2 + (j * PW) / (2 * NT), bot - (j % 2 ? TH : 0), 0.1]);
    p.push([PW / 2, TOP, 0.1]);
    return p;
  };
  const outline = ln(outlinePts(BOT), { color: 'line', width: 1.3, closed: true, order: 1 });

  /* header rows */
  const idL = lab(NUMBERS.cell, { anchor: 'left', mono: true, size: 12 });
  const typeL = lab(NUMBERS.cellType, { anchor: 'left', color: 'muted', size: 12 });
  const qL = lab('question', { anchor: 'left', color: 'muted', size: 12 });
  const qV = QS.map((q) => lab(q.name, { anchor: 'left', color: 'accent', weight: 700, size: 13 }));
  const div1 = ln([[XL, Y.div1, 0.5], [-XL, Y.div1, 0.5]], { color: 'line', width: 1.2, dashed: [3, 3] });
  const evH = lab('evidence · TEDDY + head', { anchor: 'left', color: 'teddy', size: 11.5, weight: 600 });
  evH.userData.text = 'evidence · TEDDY + head';
  const scH = lab('lineage score', { anchor: 'left', color: 'muted', size: 11.5 });

  /* evidence bars (TEDDY + head), grouped by lineage */
  const bars = MK.map((_, i) => { const m = rect('teddy', 1, 2); m.position.set(BX[i], Y.base, 1); return m; });
  const evBase = ln([[XL - 0.1, Y.base, 1.2], [BX[8] + BW / 2 + 0.1, Y.base, 1.2]], { color: 'faint', width: 1, order: 3 });
  const gL = LIN.map((l, g) => lab(l.name, { color: l.id, weight: 600, size: 11.5 }));

  /* lineage columns + the question's bar */
  const cols = LIN.map((l, i) => { const m = rect(l.id, 1, 2); m.position.set(CX[i], Y.base, 1); return m; });
  const colBase = ln([[LX0, Y.base, 1.2], [LX1, Y.base, 1.2]], { color: 'faint', width: 1, order: 3 });
  const cL = LIN.map((l) => lab(l.name, { color: l.id, weight: 600, size: 11.5 }));
  const barLine = ln([[LX0, 0, 2], [LX1, 0, 2]], { color: 'ink', width: 1.5, dashed: [5, 4], order: 3 });
  const barL = QS.map((q) => lab('bar ' + q.bar, { anchor: 'left', color: 'ink', size: 11.5 }));
  const callL = LIN.map(() => lab('call', { anchor: 'bottom', weight: 700, size: 11.5 }));

  /* drag handle */
  const ghost = ln([[0, 0, 1.5], [1, 0, 1.5]], { color: 'muted', width: 1.4, dashed: [3, 2], order: 4 });
  const ring = add(ctx.dot([0, 0, 3], { px: 7, color: 'ink', hollow: true, ring: 0.34, order: 5 }));
  const pin = add(ctx.dot([0, 0, 3], { px: 2.2, color: 'ink', order: 5 }));
  const mkL = lab('', { anchor: 'bottom', weight: 700, size: 12 });

  /* record: header + three rows */
  const div2 = ln([[XL, Y.div2, 0.5], [-XL, Y.div2, 0.5]], { color: 'line', width: 1.2, dashed: [3, 3] });
  const rhL = ['question', 'call', 'flips if'].map((t) => lab(t, { anchor: 'left', color: 'faint', size: 11 }));
  const rows = REC.map((r) => ({
    n: lab(r.name, { anchor: 'left', weight: 600, size: 12 }),
    c: lab(r.call, { anchor: 'left', weight: 600, color: r.col, size: 12 }),
    f: lab(r.flip, { anchor: 'left', color: 'ink', size: 12 }),
  }));

  /* stamps */
  const stamps = STAMPS.map((s) => {
    const g = add(new THREE.Group());
    g.rotation.z = ROT;
    const box = ctx.box(1, 0.6, { color: 'soft', stroke: s.color, strokeWidth: 2, radius: 0.1, order: 25 });
    box.position.z = 4;
    g.add(box);
    const l = ctx.label(s.text, { size: 13.5, weight: 700, color: s.color, order: 30 });
    l.material.rotation = ROT;
    l.position.z = 4.5;
    g.add(l);
    return Object.assign({ g, box, l }, s);
  });

  /* text measuring (CSS px incl. the label's padding) for layout that must not collide at phone width */
  const mc = document.createElement('canvas').getContext('2d');
  const fam = (mono) => {
    const cs = getComputedStyle(ctx.figure);
    return mono ? ((cs.getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace') : (cs.fontFamily || 'system-ui, sans-serif');
  };
  const tw = (text, size, weight = 500, mono = false) => {
    if (!mc) return text.length * size * 0.6;
    mc.font = `${weight} ${size * ctx.textScale}px ${fam(mono)}`;
    return mc.measureText(text).width + 3;
  };

  let ppu = 30, fontState = '';
  function layout() {
    ppu = Math.max(1, Math.min((ctx.width || 400) / 12, (ctx.height || 300) / 9));
    fontState = document.fonts ? document.fonts.status : 'x';
    const W = (t, s, w, m) => tw(t, s, w, m) / ppu;
    anmL.position.set(-PW / 2, Y.anm, 6);
    recL.position.set(-PW / 2 + W('ANM', 12.5, 700) + 0.12, Y.anm, 6);
    idL.position.set(XL, Y.head, 2);
    typeL.position.set(XL + W(NUMBERS.cell, 12, 500, true) + 0.25, Y.head, 2);
    qL.position.set(XL, Y.q, 2);
    const qx = XL + W('question', 12) + 0.22;
    qV.forEach((l) => l.position.set(qx, Y.q, 2));
    evH.position.set(XL, Y.ch, 2);
    scH.setText(XL + W(evH.userData.text, 11.5, 600) + 0.3 > LX0 ? 'score' : 'lineage score').position.set(LX0, Y.ch, 2);
    gL.forEach((l, g) => l.position.set(GX[g], Y.glab, 2));
    cL.forEach((l, i) => l.position.set(CX[i], Y.glab, 2));
    // record columns: widest entry of each column + a gap
    const cx = XL + Math.max(...REC.map((r) => W(r.name, 12, 600)), W('question', 11)) + 0.35;
    const fx = cx + Math.max(...REC.map((r) => W(r.call, 12, 600)), W('call', 11)) + 0.35;
    [XL, cx, fx].forEach((x, k) => rhL[k].position.set(x, Y.rh, 2));
    rows.forEach((r, i) => { r.x = [XL, cx, fx]; r.n.position.set(XL, Y.rows[i], 2); });
    const sh = (13.5 * ctx.textScale * 1.55) / ppu;
    // “no call” lands between the tallest column and the raised bar; “fragile” above the near-tied call
    stamps.forEach((s) => {
      s.box.setSize(W(s.text, 13.5, 700) + 0.4, sh, 0.1);
      const q = QS[s.q], top = Y.base + Math.max(...score(q, EV0)) * SC, bar = Y.base + q.bar * SC;
      const y = s.text === 'no call' ? lerp(top, bar, 0.56) : Math.max(top, bar) + (12 * ctx.textScale * 1.3) / ppu + 0.2 + sh / 2;
      s.g.position.set(STX, y, 0);
    });
  }

  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    if (document.fonts && document.fonts.status !== fontState) layout();
    const tear = seg(u, U.tear[0], U.tear[1], ease.inCubic);
    const G = 1 - tear;
    rc.position.y = -1.3 * tear;

    /* printer label stays; paper prints out of the slot */
    const bot = lerp(TOP, BOT, seg(u, U.print[0], U.print[1], ease.inOutSine));
    const len = TOP - bot;
    paper.scale.set(PW, Math.max(len, 1e-4), 1);
    paperMat.opacity = G * clamp(len / 0.2);
    teeth.position.set(0, bot, 0);
    teethMat.opacity = G * clamp(len / 0.4);
    outline.setPoints(outlinePts(bot)).setOpacity(G * clamp(len / 0.2));
    const rev = (y) => G * clamp((y - 0.12 - bot) / 0.3);

    idL.setOpacity(rev(Y.head)); typeL.setOpacity(rev(Y.head));
    qL.setOpacity(rev(Y.q));
    div1.setOpacity(rev(Y.div1)); div2.setOpacity(rev(Y.div2));
    evH.setOpacity(rev(Y.ch)); scH.setOpacity(rev(Y.ch));
    evBase.setOpacity(rev(Y.base)); colBase.setOpacity(rev(Y.base));
    gL.forEach((l) => l.setOpacity(rev(Y.glab))); cL.forEach((l) => l.setOpacity(rev(Y.glab)));
    rhL.forEach((l) => l.setOpacity(rev(Y.rh)));

    /* question: which rule is active (a → b blend during a switch) */
    const [qa, qb, k] = qState(u);
    const qw = qWeight(u);
    qV.forEach((l, i) => l.setOpacity(G * qw[i]));
    const ev = evAt(u);
    const grow = seg(u, U.grow[0], U.grow[1], ease.outCubic);
    const wop = (w) => 0.22 + 0.78 * w;
    bars.forEach((b, i) => {
      b.scale.set(BW, Math.max(1e-4, ev[i] * SC * grow), 1);
      b.setOpacity(rev(Y.base) * lerp(wop(QS[qa].w[i]), wop(QS[qb].w[i]), k));
    });

    const sA = score(QS[qa], ev), sB = score(QS[qb], ev);
    const cg = seg(u, U.cols[0], U.cols[1], ease.outCubic);
    const cw = callW(u);
    const colTop = [];
    cols.forEach((c, i) => {
      const h = lerp(sA[i], sB[i], k) * SC * cg;
      colTop.push(Y.base + h);
      c.scale.set(CW, Math.max(1e-4, h), 1);
      c.setOpacity(G * (cg > 0 ? 1 : 0) * (0.42 + 0.58 * cw[i]));
      callL[i].position.set(CX[i], Y.base + h + 0.1, 2);
      callL[i].setOpacity(G * cw[i]);
    });
    const bar = lerp(QS[qa].bar, QS[qb].bar, k), by = Y.base + bar * SC;
    const bk = G * seg(u, 0.075, 0.1);
    barLine.setPoints([[LX0, by, 2], [LX1, by, 2]]).setOpacity(bk);
    barL.forEach((l, i) => { l.position.set(LX1 + 0.12, by, 2); l.setOpacity(bk * qw[i]); });

    /* the drag handle on the deciding marker */
    let hk = 0, hd = null;
    for (const d of DRAGS) {
      const a = seg(u, d.show[0], d.show[1]) * (1 - seg(u, d.hide[0], d.hide[1]));
      if (a > hk && d.i >= 0) { hk = a; hd = d; }
    }
    if (hd) {
      const x = BX[hd.i], y0 = Y.base + EV0[hd.i] * SC, y = Y.base + ev[hd.i] * SC;
      const moving = seg(u, hd.go[0], hd.go[0] + 0.01) * (1 - seg(u, hd.back[1] - 0.005, hd.back[1] + 0.005));
      ghost.setPoints([[x - BW / 2 - 0.07, y0, 1.5], [x + BW / 2 + 0.07, y0, 1.5]]).setOpacity(G * hk * moving);
      ring.position.set(x, y, 3); pin.position.set(x, y, 3);
      // a nudge too small to see: the ring pulses while it pushes
      const push = Math.sin(Math.PI * clamp((u - hd.go[0]) / (hd.go[1] - hd.go[0])));
      ring.setPx(7 + 2.2 * push);
      ring.setOpacity(G * hk); pin.setOpacity(G * hk);
      mkL.setText(hd.f.marker);
      mkL.position.set(x, Math.max(y0, y) + 9 / ppu + 0.08, 3);
      mkL.setOpacity(G * hk);
    } else {
      ghost.setOpacity(0); ring.setOpacity(0); pin.setOpacity(0); mkL.setOpacity(0);
    }

    /* the record prints row by row */
    rows.forEach((r, i) => {
      const a = seg(u, ROWT[i].call[0], ROWT[i].call[1]), f = seg(u, ROWT[i].flip[0], ROWT[i].flip[1]);
      const y = Y.rows[i];
      r.n.setOpacity(G * a); r.n.position.set(r.x[0], y, 2);
      r.c.setOpacity(G * a); r.c.position.set(r.x[1], y + 0.12 * (1 - a), 2);
      r.f.setOpacity(G * f); r.f.position.set(r.x[2], y + 0.12 * (1 - f), 2);
    });

    /* stamps land on the lineage columns */
    stamps.forEach((s) => {
      const a = seg(u, s.in[0], s.in[1], ease.outCubic) * (s.out ? 1 - seg(u, s.out[0], s.out[1]) : 1);
      const land = seg(u, s.in[0], s.in[1], ease.outCubic);
      const sc = lerp(1.25, 1, land);
      s.box.scale.set(sc, sc, 1);
      s.box.setOpacity(G * a * 0.94);
      s.l.setOpacity(G * a);
    });
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: 0.87 * PERIOD,
    update, resize: layout,
  };
}
