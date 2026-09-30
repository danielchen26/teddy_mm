/*
 * idx-q-missing-input — index.html #problem, GAP 2 card "Missing input" (thumbnail at the top of the card, 16:9).
 * The problem only; Experiment 2 (#missing, withdrawn) is where the page answers it.
 * The card: "RNA missing, only a weaker protein-only model?" / "No RNA, no embedding, no prediction. Withdrawn."
 * The block: without RNA, our head gets no embedding, so a pipeline falls back on another model; the fixed rule
 * is not told this evidence is weaker and keeps answering.
 * 1 RNA beads stream into TEDDY; our head turns TEDDY's embedding into predicted proteins (crisp orange bars),
 *   plugged into one socket; the fixed rule reads them and stamps each cell on the belt.
 * 2 The RNA stream stops and drains ("no RNA"); TEDDY and our head go grey ("no embedding"); the head's plug
 *   leaves the socket and the bars empty ("no prediction": the card's own three phrases). The belt waits:
 *   with no prediction there is nothing to read.
 * 3 A grey, blurred protein-only stand-in slides up and plugs into the same socket. Its bars are grey and
 *   fuzzy: weaker evidence, drawn as blur, not as values. The fixed rule stamps every cell with exactly the
 *   same mark as before ("keeps answering").
 * A "withdrawn" tag sits in the bottom-left corner (the top-right corner stays clear for the play button).
 * No trust knob and no "no call" here: that is the experiment's scene. Bar heights are illustrative, not data.
 * 12 s seamless loop (8 stamps: 4 fed by TEDDY + head, 4 by the stand-in); the still frame is beat 3.
 */

/* Every pipeline number this scene prints (none). Refresh after the official-preprocessing rerun. */
const NUMBERS = {};

export default function create(ctx) {
  const { THREE, ease, seg, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12, BEAT = 1.2;
  const TS = 4.5, TE = 6.9;            // the belt waits TS..TE (stamp up, cells still); 12 − 2.4 = 8 beats
  const CUT = 1.2;                      // the RNA source stops
  const PRESS = [1.2, 2.4, 3.6, 7.2, 8.4, 9.6, 10.8, 12];
  const byStandIn = (p) => p > TE && p < 11;
  void NUMBERS;

  const put = (o, x, y, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };
  const lab = (s, x, y, o = {}, parent = scene) =>
    put(ctx.label(s, Object.assign({ size: 11.5, color: 'muted' }, o)), x, y, 10, parent);

  /* ── layout (world units, 16 × 9 frame) ── */
  const YA = 1.6;                                            // the pipeline row
  const XS = -7.75, TX = -3.75, TW = 2.3, TH = 1.8, XE = TX - TW / 2;
  const HX = -0.95, HW = 1.7, HH = 1.2, HR = HX + HW / 2;
  const PX = 0.62, SOX = 0.8;                                // plug tip, socket
  const BY = 0.7, AX = 1.05, BW = 0.3;                       // bar baseline, chart axis, bar width
  const BX = [1.35, 1.85, 2.35, 2.85, 3.35, 3.85];
  const HT = [1.5, 0.8, 1.85, 1.1, 0.55, 1.35];              // TEDDY + head bars (illustrative)
  const HS = [1.15, 1.0, 1.3, 0.95, 0.85, 1.1];              // stand-in bars (illustrative)
  const WOB = [[3, 0.1], [4, 0.55], [3, 0.8], [5, 0.3], [4, 0.9], [3, 0.45]];
  const SX = -2.3, SY = -1.45, SW = 2.6, SH = 1.1, SY0 = -6.4, SR = SX + SW / 2, PV = 0.3;
  const STX = 4.6, BELT = -2.55, R = 0.27, CY = BELT + 0.35, SP = 1.0, STAMP = 3, X0 = STX - STAMP * SP, N = 7;
  const PAD = 1.1, PH = 0.62, HNDL = 0.4, KNOB = 0.2, TRAVEL = 0.5;

  /* ── 1 · RNA beads flowing into TEDDY ── */
  const NB = 8, GAP = 0.5, V = 1.25, WR = NB * GAP;         // V × PERIOD = 30 gaps: seamless
  const beads = Array.from({ length: NB }, () => put(ctx.dot([0, 0], { r: 0.13, color: 'muted' }), 0, 0, 1));
  const rnaOn = lab('RNA', (XS + XE) / 2, YA + 0.42, { anchor: 'bottom', color: 'ink', weight: 600 });
  const rnaOff = lab('no RNA', (XS + XE) / 2, YA + 0.42, { anchor: 'bottom', color: 'muted', weight: 600 });

  /* TEDDY: an opaque slab, lit (orange) or dark (grey dashed) */
  put(ctx.box(TW, TH, { color: 'card', radius: 0.2, order: 3 }), TX, YA, 2);
  const tFillL = put(ctx.box(TW, TH, { color: 'teddy', opacity: 0.16, radius: 0.2, order: 4 }), TX, YA, 2.1);
  const tFillD = put(ctx.box(TW, TH, { color: 'faint', opacity: 0, radius: 0.2, order: 4 }), TX, YA, 2.1);
  const tStrL = put(ctx.box(TW, TH, { color: null, stroke: 'teddy', radius: 0.2, order: 5 }), TX, YA, 2.2);
  const tStrD = put(ctx.box(TW, TH, { color: null, stroke: 'faint', dashed: [4, 3], radius: 0.2, order: 5 }), TX, YA, 2.2);
  const tLabL = lab('TEDDY', TX, YA, { size: 13, weight: 600, color: 'teddy' });
  const tLabD = lab('TEDDY', TX, YA, { size: 13, weight: 600, color: 'faint' });

  /* the embedding link TEDDY → our head */
  const lkL = put(ctx.arrow([TX + TW / 2 + 0.05, YA], [HX - HW / 2 - 0.05, YA], { color: 'teddy', width: 1.8, head: 7 }), 0, 0, 1);
  const lkD = put(ctx.line([[TX + TW / 2 + 0.05, YA], [HX - HW / 2 - 0.05, YA]], { color: 'faint', width: 1.5, dashed: [3, 3] }), 0, 0, 1);
  const noEmb = lab('no embedding', -2.2, 0.45, { anchor: 'top' });

  /* our head */
  put(ctx.box(HW, HH, { color: 'card', radius: 0.18, order: 3 }), HX, YA, 2);
  const hFill = put(ctx.box(HW, HH, { color: 'soft', radius: 0.18, order: 4 }), HX, YA, 2.1);
  const hStrL = put(ctx.box(HW, HH, { color: null, stroke: 'muted', radius: 0.18, order: 5 }), HX, YA, 2.2);
  const hStrD = put(ctx.box(HW, HH, { color: null, stroke: 'faint', dashed: [4, 3], radius: 0.18, order: 5 }), HX, YA, 2.2);
  const hLabL = lab('our head', HX, YA + HH / 2 + 0.1, { anchor: 'bottom', size: 12.5, weight: 600, color: 'ink' });
  const hLabD = lab('our head', HX, YA + HH / 2 + 0.1, { anchor: 'bottom', size: 12.5, weight: 600, color: 'faint' });

  /* the head's plug, the socket, the protein chart */
  const hPlug = put(ctx.line([[HR, YA], [PX, YA]], { color: 'muted', width: 2 }), 0, 0, 1.5);
  const hTip = put(ctx.box(0.16, 0.34, { color: 'muted', radius: 0.04, order: 2 }), PX, YA, 1.6);
  put(ctx.box(0.14, 0.62, { color: 'card', stroke: 'muted', strokeWidth: 1.5, radius: 0.04, order: 2 }), SOX, YA, 1.7);
  put(ctx.line([[SOX + 0.07, YA], [AX, YA]], { color: 'line', width: 1.5 }), 0, 0, 0.5);
  put(ctx.line([[AX, BY + 2.05], [AX, BY], [BX[5] + BW / 2 + 0.2, BY]], { color: 'line', width: 1.5 }), 0, 0, 0.5);
  lab('predicted\nproteins', AX, BY - 0.12, { anchor: 'top-left' });
  /* inside the emptied chart, clear of its axis on the left and of the feed arrow on the right */
  const noPred = lab('no prediction', AX + 0.2, BY + 0.55, { anchor: 'left', size: 10.5 });
  const tBars = HT.map((h, i) => put(ctx.box(1, 1, { color: 'teddy', radius: 0, order: 3 }), BX[i], BY, 1));
  const HALO = [[0.08, 0.3, 0.16], [0.16, 0.62, 0.09], [0.24, 0.95, 0.05]];   // extra width, height, opacity
  const sBars = HS.map((h, i) => ({
    halos: HALO.map((q, k) => put(ctx.box(1, 1, { color: 'faint', opacity: 0, radius: 0 }), BX[i], BY, 0.9 - 0.05 * k)),
    bar: put(ctx.box(1, 1, { color: 'faint', opacity: 0, radius: 0, order: 2 }), BX[i], BY, 1),
  }));

  /* ── 3 · the protein-only stand-in: grey, blurred (layered soft edges), dashed outline ── */
  const sg = put(ctx.group(), SX, SY0, 3);
  const sLayers = [[0.5, 0.05], [0.26, 0.09], [0, 0.2]].map(([d, o]) => ({
    b: put(ctx.box(SW + d, SH + d, { color: 'faint', opacity: o, radius: 0.22 + d / 2, order: 3 }), 0, 0, 0.1 - d / 10, sg), o,
  }));
  const sStr = put(ctx.box(SW, SH, { color: null, stroke: 'faint', dashed: [4, 3], radius: 0.22, order: 4 }), 0, 0, 0.2, sg);
  const sName = lab('stand-in', 0, 0, { size: 12.5, weight: 600, color: 'muted' }, sg);
  const sSub = lab('protein-only', 0, -SH / 2 - 0.12, { anchor: 'top' }, sg);
  const sPts = [[SR, SY], [PV, SY], [PV, YA], [PX, YA]];
  const sPlug = put(ctx.line(sPts, { color: 'faint', width: 2 }), 0, 0, 1.4);
  const sTip = put(ctx.box(0.16, 0.34, { color: 'faint', radius: 0.04, order: 2 }), PX, YA, 1.6);

  /* ── the fixed rule: reads the bars, stamps the cell under it ── */
  const feed = put(ctx.arrow([BX[5] + 0.4, YA + 0.15], [STX, BELT + 2.57], { color: 'faint', width: 1.4, head: 7, bend: 0.22, opacity: 0.6 }), 0, 0, 0.5);
  const ev = put(ctx.dot([0, 0], { px: 3.6, color: 'teddy', opacity: 0 }), 0, 0, 2);
  put(ctx.line([[X0 + 0.75 * SP, BELT], [X0 + 6.8 * SP, BELT]], { color: 'line', width: 1.5 }), 0, 0, 0);
  const cells = Array.from({ length: N }, () => {
    const g = put(ctx.group(), 0, 0, 1);
    const ring = put(ctx.dot([0, 0], { r: R, hollow: true, ring: 0.3, color: 'faint', order: 1 }), 0, 0, 0, g);
    const mark = put(ctx.box(0.2, 0.2, { color: 'teddy', radius: 0.04, order: 1 }), 0, 0, 0.1, g);
    return { g, ring, mark, done: null };
  });
  const stamp = put(ctx.group(), STX, CY + R + TRAVEL, 3);
  put(ctx.box(PAD, PH, { color: 'teddy', opacity: 0.16, radius: 0.1, order: 5 }), 0, PH / 2, 0, stamp);
  put(ctx.box(PAD, PH, { color: null, stroke: 'teddy', radius: 0.1, order: 6 }), 0, PH / 2, 0.01, stamp);
  put(ctx.box(PAD - 0.1, 0.1, { color: 'teddy', radius: 0.03, order: 6 }), 0, 0.05, 0.02, stamp);
  put(ctx.box(0.22, HNDL, { color: 'teddy', radius: 0.05, order: 5 }), 0, PH + HNDL / 2, 0, stamp);
  put(ctx.box(0.6, KNOB, { color: 'teddy', radius: 0.1, order: 5 }), 0, PH + HNDL + KNOB / 2, 0, stamp);
  lab('fixed rule', STX, BELT - 0.15, { anchor: 'top', size: 12.5, weight: 600, color: 'ink' });
  const keeps = lab('keeps answering', STX, BELT - 0.7, { anchor: 'top', size: 12, weight: 600, color: 'teddy' });

  lab('withdrawn', XS, -4.2, { anchor: 'bottom-left', size: 10.5, mono: true, weight: 600, color: 'bad', bg: 'soft', bgOpacity: 1, pad: 4 });

  function layout() {
    const ppu = ctx.ppu();
    beads.forEach((d) => d.setPx(clamp(0.12 * ppu, 2.4, 5)));
    keeps.position.y = BELT - 0.15 - (12.5 * ctx.textScale * 1.3 + 2) / ppu;
  }

  const tauOf = (tt) => (tt < TS ? tt : tt < TE ? TS : tt - (TE - TS));
  let evTok = '';

  function update(t) {
    const u = ctx.loopT(t, PERIOD), tt = u * PERIOD;

    /* RNA: after CUT no bead is born, so the stream's tail recedes into TEDDY; it fades back in for the loop */
    const back = seg(tt, 10.9, 11.4, ease.inOutSine);
    beads.forEach((d, k) => {
      const x = XS + ((((V * tt + k * GAP) % WR) + WR) % WR);
      const born = tt - (x - XS) / V;
      const on = born <= CUT ? 1 : Math.max(back, 1 - ease.smooth((born - CUT) / 0.3));
      d.position.set(x, YA + 0.13 * Math.sin(2.4 * x), 1);
      d.setOpacity(on * clamp((x - XS) / 0.4) * clamp((XE - 0.05 - x) / 0.3));
    });
    /* "RNA" → "no RNA": one label out, then the other in (no overlap) */
    rnaOn.setOpacity(1 - seg(tt, 1.3, 1.55) + seg(tt, 11.15, 11.4));
    rnaOff.setOpacity(seg(tt, 1.55, 1.8) * (1 - seg(tt, 10.95, 11.15)));

    /* TEDDY and our head go dark once the last RNA is in; they relight for the loop */
    const kt = seg(tt, 3.0, 3.7) * (1 - seg(tt, 10.95, 11.4));
    tFillL.setOpacity(0.16 * (1 - kt)); tStrL.setOpacity(1 - kt); tLabL.setOpacity(1 - kt);
    tFillD.setOpacity(0.1 * kt); tStrD.setOpacity(kt); tLabD.setOpacity(kt);
    lkL.setOpacity(1 - kt); lkD.setOpacity(kt);
    const kh = seg(tt, 3.5, 4.1) * (1 - seg(tt, 11.0, 11.4));
    hFill.setOpacity(1 - 0.7 * kh); hStrL.setOpacity(1 - kh); hStrD.setOpacity(kh);
    hLabL.setOpacity(1 - kh); hLabD.setOpacity(kh);
    noEmb.setOpacity(seg(tt, 3.9, 4.3) * (1 - seg(tt, 10.9, 11.2)));
    noPred.setOpacity(seg(tt, 4.25, 4.6) * (1 - seg(tt, 5.5, 5.8)));

    /* the head's plug comes out; its bars empty; later they come back */
    const kp = clamp(1 - seg(tt, 4.1, 4.5) + seg(tt, 11.05, 11.35));
    hPlug.setProgress(kp).setOpacity(kp > 0.01 ? 1 : 0);
    hTip.position.x = lerp(HR + 0.08, PX, kp);
    hTip.setOpacity(kp > 0.01 ? 1 : 0);
    const kb = clamp(1 - seg(tt, 3.75, 4.25) + seg(tt, 11.15, 11.5, ease.outCubic));
    tBars.forEach((b, i) => {
      const h = Math.max(1e-3, HT[i] * kb);
      b.scale.set(BW, h, 1);
      b.position.y = BY + h / 2;
      b.visible = kb > 0.002;
    });

    /* the stand-in slides up, plugs into the same socket, and gives grey, fuzzy bars */
    const ks = seg(tt, 4.5, 5.3, ease.outCubic) * (1 - seg(tt, 10.95, 11.45, ease.inCubic));
    sg.position.y = lerp(SY0, SY, ks);
    sLayers.forEach((l, i) => {
      l.b.setOpacity(l.o * ks);
      l.b.position.x = 0.05 * Math.sin(2 * Math.PI * (u * (3 + i) + 0.3 * i));
    });
    sStr.setOpacity(ks); sName.setOpacity(ks); sSub.setOpacity(ks);
    const kq = seg(tt, 5.3, 5.75) * (1 - seg(tt, 10.9, 11.1));
    sPlug.setProgress(kq).setOpacity(kq > 0.01 ? 1 : 0);
    sTip.position.copy(ctx.along(sPts, kq)).setZ(1.6);
    sTip.setOpacity(ease.smooth((kq - 0.85) / 0.15));
    const kg = seg(tt, 5.7, 6.3) * (1 - seg(tt, 10.95, 11.25));
    sBars.forEach((s, i) => {
      const [c, ph] = WOB[i];
      const h = Math.max(1e-3, kg * HS[i] * (1 + 0.14 * Math.sin(2 * Math.PI * (u * c + ph))));
      s.bar.scale.set(BW, h, 1); s.bar.position.y = BY + h / 2; s.bar.setOpacity(0.75 * kg);
      s.halos.forEach((b, k) => {
        const [dw, hh, o] = HALO[k];
        b.scale.set(BW + dw, Math.min(hh, 2 * h + 0.02), 1); b.position.y = BY + h; b.setOpacity(o * kg);
      });
    });

    /* one reading per stamp travels from the bars to the fixed rule: orange from TEDDY + head, grey from the stand-in */
    let evA = 0;
    for (const P of PRESS) {
      const k = (tt - (P - 0.62)) / 0.56;
      if (k <= 0 || k >= 1) continue;
      const tok = byStandIn(P) ? 'faint' : 'teddy';
      if (tok !== evTok) { ev.setColor(tok); evTok = tok; }
      ev.position.copy(ctx.along(feed.line, ease.inOutSine(k))).setZ(2);
      evA = Math.min(1, k / 0.12, (1 - k) / 0.12);
    }
    ev.setOpacity(evA);

    /* the fixed rule: the same stamp on every cell, whatever fed it; the belt waits while nothing does */
    const tau = tauOf(tt), b = Math.floor(tau / BEAT + 1e-9), ph = tau / BEAT - b;
    const press = 1 - ease.smooth((Math.min(ph, 1 - ph) - 0.04) / 0.2);
    stamp.position.y = CY + R + TRAVEL * (1 - press);
    const e = ease.inOutSine((ph - 0.28) / 0.44);
    cells.forEach((c, j) => {
      const p = (j + b + e) % N;
      c.g.position.set(X0 + p * SP, CY, 1);
      const a = clamp((p - 0.45) / 0.6) * clamp((6.3 - p) / 0.5);
      const at = Math.abs(p - STAMP) < 1e-4 && ph < 0.5, done = p > STAMP + 1e-4;
      c.ring.setOpacity(a);
      const tok = done || at ? 'ink' : 'faint';
      if (tok !== c.done) { c.ring.setColor(tok); c.done = tok; }
      c.mark.visible = at || done;
      c.mark.scale.setScalar(at ? Math.max(0.01, ease.outBack(ph / 0.09)) : 1);
      c.mark.setOpacity(a);
    });
    keeps.setOpacity(seg(tt, 7.3, 7.7) * (1 - seg(tt, 10.85, 11.1)));
  }

  let dead = false;
  layout();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 8.1,
    update, resize: layout,
    dispose() { dead = true; beads.length = 0; cells.length = 0; },
  };
}
