/*
 * idx-q-explain-call — index.html #problem, GAP 3 card "Explain a call" (thumbnail at the top of the card, 16:9).
 * The problem only; idx-exp-attr-leave-one-out (#attr) shows how the test answers it.
 * The card: "Which marker drove this call; how close was a flip?" / TEDDY + fixed rule: "Neither, until written
 * down; then the same leave-one-out runs on it." The block: "Its call is the top lineage score, no reason."
 * The rule and ANM read the same 9 predicted proteins (our head's output on TEDDY's frozen embedding).
 * 1 Nine grey bars, the 9 predicted proteins, each send a grey dot into a closed box, the fixed rule.
 *   A call token "B" slides out of the box. Under it, two fields stay blank: reason, margin.
 * 2 A question mark hovers over the bars, one after another ("which marker drove it?"). Every bar stays the
 *   same grey; none is marked or lit.
 * 3 The question mark moves to the blank margin field ("how close to a flip?"); the field stays blank.
 *   The call fades and the loop starts again.
 * Nothing is removed and nothing flips: that is the experiment's scene. Bar heights are illustrative, not data.
 * The top-right corner stays clear for the play button.
 * 10 s seamless loop; the still frame shows the question mark over the fifth bar.
 */

/* Every pipeline number this scene shows or uses. On index.html: "Rule and ANM read the same 9 predicted
 * proteins", #attr "drop each of the 9 in turn". Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  markers: 9,
};

export default function create(ctx) {
  const { ease, seg, clamp, lerp } = ctx;
  const scene = new ctx.THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 10;
  const N = NUMBERS.markers;

  const put = (o, x, y, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };
  const lab = (s, x, y, o = {}, parent = scene) =>
    put(ctx.label(s, Object.assign({ size: 12, color: 'muted' }, o)), x, y, 10, parent);

  /* ── layout (world units, 16 × 9 frame) ── */
  const YO = -0.65;                                   // the row of box and call
  const YB = -1.75, BW = 0.32, DX = 0.55;             // bar baseline, bar width, bar spacing
  const XB = Array.from({ length: N }, (_, i) => -7.0 + i * DX);
  const HB = [1.5, 2.2, 1.1, 1.8, 2.5, 1.3, 2.0, 0.9, 1.6].slice(0, N); // illustrative heights
  while (HB.length < N) HB.push(1.4);
  const BX = 0.55, BOXW = 3.0, BOXH = 2.3;            // the fixed rule's box
  const TX = 3.95, TS = 1.35;                         // the call token: resting x, size
  const FX = TX - TS / 2, FE = 7.15;                  // the blank fields: left edge, right end
  const YR = [-2.0, -3.0];                            // field rows: reason, margin
  const QY = 2.15;                                    // hover height of the question mark over the bars

  /* questions (top-left; they cross-fade) */
  const qA = lab('which marker drove it?', -7.45, 3.55, { anchor: 'left', size: 13, weight: 600, color: 'ink' });
  const qB = lab('how close to a flip?', -7.45, 3.55, { anchor: 'left', size: 13, weight: 600, color: 'ink', opacity: 0 });

  /* ── the 9 predicted proteins: the same grey, always ── */
  put(ctx.line([[XB[0] - 0.35, YB], [XB[N - 1] + 0.35, YB]], { color: 'line', width: 1.5 }), 0, 0, 0.5);
  HB.forEach((h, i) => {
    const b = put(ctx.box(1, 1, { color: 'faint', opacity: 0.8, radius: 0, order: 1 }), XB[i], YB + h / 2, 1);
    b.scale.set(BW, h, 1);
  });
  lab(`${N} predicted proteins`, (XB[0] + XB[N - 1]) / 2, YB - 0.2, { anchor: 'top' });

  /* one reading per bar travels into the box (all 9 at once) */
  const entry = [BX - 0.5, YO + 0.2];
  const paths = XB.map((x, i) => {
    const p0 = [x, YB + HB[i] + 0.12], p2 = entry, c = [(x + entry[0]) / 2 + 0.4, 2.0];
    return [p0, c, p2];
  });
  const quad = (P, k) => {
    const a = (1 - k) * (1 - k), b = 2 * k * (1 - k), d = k * k;
    return [a * P[0][0] + b * P[1][0] + d * P[2][0], a * P[0][1] + b * P[1][1] + d * P[2][1]];
  };
  const reads = XB.map(() => put(ctx.dot([0, 0], { px: 3, color: 'muted', opacity: 0, order: 2 }), 0, 0, 1.2));
  put(ctx.arrow([XB[N - 1] + 0.55, YO], [BX - BOXW / 2 - 0.1, YO], { color: 'line', width: 1.5, head: 7 }), 0, 0, 0.5);

  /* ── the fixed rule: a closed box (drawn above the token, so the call slides out from inside it) ── */
  const boxG = put(ctx.group(), BX, YO, 3);
  put(ctx.box(BOXW, BOXH, { color: 'card', radius: 0.2, order: 6 }), 0, 0, 0, boxG);
  put(ctx.box(BOXW, BOXH, { color: 'teddy', opacity: 0.13, radius: 0.2, order: 7 }), 0, 0, 0.01, boxG);
  put(ctx.box(BOXW, BOXH, { color: null, stroke: 'teddy', strokeWidth: 1.6, radius: 0.2, order: 8 }), 0, 0, 0.02, boxG);
  put(ctx.line([[-BOXW / 2, BOXH / 2 - 0.5], [BOXW / 2, BOXH / 2 - 0.5]], { color: 'teddy', width: 1.3, order: 8 }), 0, 0, 0.03, boxG);
  lab('fixed rule', 0, -0.3, { size: 13, weight: 600, color: 'teddy' }, boxG);

  /* ── the call token "B", and its blank fields ── */
  const tok = put(ctx.group(), BX, YO, 1);
  const tFill = put(ctx.box(TS, TS, { color: 'card', radius: 0.18, order: 2 }), 0, 0, 0, tok);
  const tTint = put(ctx.box(TS, TS, { color: 'b', opacity: 0.12, radius: 0.18, order: 3 }), 0, 0, 0.01, tok);
  const tStr = put(ctx.box(TS, TS, { color: null, stroke: 'b', strokeWidth: 1.8, radius: 0.18, order: 3 }), 0, 0, 0.02, tok);
  const tLet = lab('B', 0, 0, { size: 24, weight: 700, color: 'b', order: 4 }, tok);
  const callLab = lab('call', TX, YO + TS / 2 + 0.14, { anchor: 'bottom', size: 11.5, opacity: 0 });
  const fields = ['reason', 'margin'].map((w, i) => ({
    w,
    name: lab(w, FX, YR[i], { anchor: 'left', size: 12, opacity: 0 }),
    blank: put(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 1.4, dashed: [4, 3], opacity: 0 }), 0, 0, 0.5),
    x0: 0,
  }));

  /* ── the question mark ── */
  const q = lab('?', XB[0], QY, { size: 20, weight: 700, color: 'ink', opacity: 0 });
  const probe = put(ctx.line([[0, 0], [0, 1]], { color: 'faint', width: 1.2, dashed: [3, 3], opacity: 0 }), 0, 0, 0.8);

  let ppu = 60;
  const est = (text, size, k = 0.56) => (text.length * size * ctx.textScale * k) / ppu;
  function layout() {
    ppu = ctx.ppu() || 60;
    fields.forEach((f, i) => {
      f.x0 = Math.min(FX + est(f.w, 12) + 0.3, FE - 0.9);
      f.blank.setPoints([[f.x0, YR[i] - 0.2, 0.5], [FE, YR[i] - 0.2, 0.5]]);
    });
  }

  /* timing, in seconds */
  const R0 = 0.3, R1 = 1.1;                 // the readings travel into the box
  const P0 = 1.2, P1 = 1.95;                // the call slides out
  const H0 = 2.5, STEP = 0.5, MOVE = 0.2;   // the question mark visits bar i at H0 + i·STEP
  const HE = H0 + (N - 1) * STEP + (STEP - MOVE);   // leaves the last bar
  const F1 = HE + 0.8;                      // lands on the margin field
  const Q0 = 8.9, Q1 = 9.3;                 // the question mark fades
  const C0 = 9.2, C1 = 9.75;                // the call and its fields fade

  function update(t) {
    const u = ctx.loopT(t, PERIOD), tt = u * PERIOD;

    /* 1 · readings: one grey dot per bar, all into the box */
    reads.forEach((d, i) => {
      const k = (tt - (R0 + 0.035 * i)) / (R1 - R0 - 0.035 * (N - 1));
      if (k <= 0 || k >= 1) { d.setOpacity(0); return; }
      const e = ease.inOutSine(k), p = quad(paths[i], e);
      d.position.set(p[0], p[1], 1.2);
      d.setOpacity(Math.min(1, k / 0.15, (1 - k) / 0.1));
    });
    /* the box takes them in (a small squeeze), then lets the call out */
    const sq = Math.sin(Math.PI * clamp((tt - R1) / 0.25));
    boxG.scale.set(1 + 0.02 * sq, 1 - 0.035 * sq, 1);

    /* the call slides out of the box and stays; it fades at the end of the loop */
    const out = seg(tt, P0, P1, ease.outBack);
    const fade = 1 - seg(tt, C0, C1);
    tok.position.set(lerp(BX, TX, out), YO, 1);
    const sc = lerp(0.82, 1, clamp(out));
    tok.scale.set(sc, sc, 1);
    const tv = tt > P0 - 0.01 ? fade : 0;
    tFill.setOpacity(tv); tTint.setOpacity(0.12 * tv); tStr.setOpacity(tv); tLet.setOpacity(tv);
    tok.visible = tv > 0.001;
    const fk = seg(tt, P1 - 0.1, P1 + 0.3) * fade;
    callLab.setOpacity(fk);
    fields.forEach((f, i) => {
      const k = seg(tt, P1 + 0.05 + 0.15 * i, P1 + 0.4 + 0.15 * i) * fade;
      f.name.setOpacity(k); f.blank.setOpacity(k);
    });

    /* 2 · the question mark visits every bar; nothing lights up */
    const qIn = seg(tt, H0 - 0.35, H0);
    const qOut = 1 - seg(tt, Q0, Q1);
    let x, y, pk = 0;
    if (tt < HE) {
      const k = clamp((tt - H0) / STEP, 0, N - 1);
      const i = Math.min(N - 1, Math.floor(k)), f = k - i;
      const hold = (STEP - MOVE) / STEP;
      const m = i < N - 1 && f > hold ? ease.inOutSine((f - hold) / (1 - hold)) : 0;
      x = lerp(XB[i], XB[Math.min(N - 1, i + 1)], m);
      y = QY + 0.22 * Math.sin(Math.PI * m);
      pk = 1 - Math.sin(Math.PI * m);                 // the probe fades while it hops
      const bar = m > 0.5 ? Math.min(N - 1, i + 1) : i;
      probe.setPoints([[x, QY - 0.45, 0.8], [x, YB + HB[bar] + 0.12, 0.8]]);
    } else {
      /* 3 · …then to the blank margin field */
      const f = fields[1];
      const k = ease.inOutSine((tt - HE) / (F1 - HE));
      const mx = (f.x0 + FE) / 2, my = YR[1] + 0.12;
      const P = [[XB[N - 1], QY], [6.5, 2.4], [mx, my]];
      [x, y] = quad(P, k);
    }
    const bob = 0.07 * Math.sin(2 * Math.PI * u * 10);
    q.position.set(x, y + bob, 10);
    q.setOpacity(qIn * qOut);
    probe.setOpacity(0.9 * pk * qIn * (tt < HE ? 1 : 0));

    /* the question on top follows the question mark */
    const kb = seg(tt, HE + 0.1, HE + 0.5) * (1 - seg(tt, Q0 + 0.1, Q1 + 0.3));
    qA.setOpacity(1 - kb);
    qB.setOpacity(kb);
  }

  let dead = false;
  layout();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 4.6,
    update, resize: layout,
    dispose() { dead = true; },
  };
}
