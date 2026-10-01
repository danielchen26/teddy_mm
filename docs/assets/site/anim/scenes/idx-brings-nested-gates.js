/*
 * idx-brings-nested-gates — index.html #brings, card “4 · Coverage vs accuracy” (“Nested readouts make the
 * trade explicit”), card width, 4:3, above table.mini.
 * The card: soft → strict → B/T-priority on the same evidence; “no calls” 448 / 2,014 / 3,169; accuracy of
 * calls 0.953 / 0.978 / 0.964; 64–91% of out-of-scope cells still get a call: the key needs NK and
 * out-of-scope classes.
 * 1 The same stream of cells (identical in the three lanes: the same evidence) falls onto three gates, one
 *   per written question. The gates stand side by side, not in series: every question reads every cell
 *   (the audit card’s cell gets “no call” under strict yet T under B/T-priority, so the questions do not
 *   filter one after another). The more a question holds back, the narrower its slot is drawn.
 * 2 A held cell stops on the closed gate and rolls into that lane’s “no call” tray; a called cell drops
 *   through and turns ANM teal. Each tray then prints the page’s “no calls”, each pile its accuracy of calls.
 * 3 Grey cells are out of scope (NK, ILC, erythroid, progenitor; legend top left). Most of them slip through
 *   every gate and get a call anyway (teal ring). A meter shows 64–91% still called and points at the answer
 *   key (B, T, myeloid), which “needs” NK and out-of-scope classes (dashed chips).
 * The top-right corner stays clear for the loader's play/pause pill; layout is sized for a ~292 px card.
 * The particles are illustrative: 20 cells a lane, 4 of them grey; lanes hold 1, 3 and 4 cells in the order
 * of the page’s counts, and 3 of 4 grey cells are called (inside 64–91%). Only NUMBERS are printed.
 * 10 s seamless loop; the still frame is the finished state.
 */

/* Every pipeline number this scene prints. Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  noCalls: [272, 1493, 2551],       // page table “No calls”: soft, strict, B/T-priority
  accuracy: [0.943, 0.967, 0.955],  // page table “Accuracy of calls”: soft, strict, B/T-priority
  oosCalled: [73, 95],              // “64–91% of out-of-scope cells still get a call” (%)
};

export default function create(ctx) {
  const { THREE, ease, seg, pulse, lerp, clamp } = ctx;
  const N = NUMBERS;
  const scene = new THREE.Scene();
  const W = 12, H = 9;
  const camera = ctx.orthoCamera({ width: W, height: H });
  const PERIOD = 10;
  const add = (o) => { scene.add(o); return o; };
  const lab = (text, x, y, o = {}) => {
    const l = ctx.label(text, Object.assign({ size: 11, color: 'muted' }, o));
    l.position.set(x, y, 10);
    return add(l);
  };

  /* ── layout (world units, 12 × 9 frame). Sized for the smallest card stage (~292 × 219 px, 24 px a unit):
     the lane names sit below the loader's play/pause pill (top-right, ≈ 84 × 36 px); the legend takes the
     free top-left corner. ── */
  const QN = ['soft', 'strict', 'B/T-priority'];
  const LX = [-4, 0, 4];                        // lane centres
  const YLEG = 3.95, YN = 2.6, Y0 = 2.17, G = 0.92; // legend, lane names, spawn height, gate
  const HB = 1.0, PKL = 1.85, PKD = 0.6;        // gate half-length; tray left edge and depth (from lane centre)
  const PB = -1.63, DP = 0.36, YACC = -1.87;    // pile base row, pile pitch, accuracy row
  const YSEP = -2.5, YB = -3.17, YC = -4.02;    // strip: separator, meter, answer key
  const XL = -5.8, XR = 5.8, XT0 = -2.9, XT1 = 5.6;

  /* ── the same 20 cells fall in every lane; 4 of them are out of scope ── */
  const NC = 20, GREY = [3, 8, 12, 17];
  const rank = N.noCalls.map((v) => N.noCalls.filter((w) => w < v).length); // 0 = holds back fewest
  const SLOT = [1.1, 0.7, 0.42];                                              // slot width by rank
  const HELD_IN = [[], [5, 14], [1, 9, 14]];                                  // in-scope cells held, by rank
  const GREY_ORDER = [[12, 3, 17, 8], [17, 8, 3, 12], [8, 17, 12, 3]];        // which grey cell a lane holds first
  const oosMid = (N.oosCalled[0] + N.oosCalled[1]) / 200;
  const greyHeld = GREY.length - clamp(Math.round(GREY.length * oosMid), 0, GREY.length);
  const S0 = 0.04, DS = 0.026, F = 0.075, GR = (Y0 - G) / (F * F);            // spawn, spacing, fall time, gravity
  const rnd = ctx.rand(2014);
  const JX = Array.from({ length: NC }, () => (rnd() - 0.5) * 0.12);

  /* ── lanes: name, gate shutters, tray, counts, cells ── */
  const lanes = LX.map((lc, q) => {
    const r = rank[q];
    const held = new Set(HELD_IN[Math.min(r, 2)].concat(GREY_ORDER[q].slice(0, greyHeld)));
    const L = { lc, slot: SLOT[Math.min(r, 2)], cells: [] };
    lab(QN[q], lc, YN, { size: 12, weight: 600, color: 'ink' });
    add(ctx.line([[lc, Y0 + 0.05, 0.4], [lc, PB + 3 * DP + 0.22, 0.4]], { color: 'line', width: 1, dashed: [2, 4] }));
    L.shL = add(ctx.line([[lc - HB, G, 1], [lc, G, 1]], { color: 'ink', width: 3.5 }));
    L.shR = add(ctx.line([[lc, G, 1], [lc + HB, G, 1]], { color: 'ink', width: 3.5 }));
    add(ctx.line([[lc - PKL, G + 0.12, 0.8], [lc - PKL, G - PKD, 0.8], [lc - HB, G - PKD, 0.8], [lc - HB, G, 0.8]],
      { color: 'muted', width: 1.5 }));
    lab('no call', lc - PKL, G + 0.34, { anchor: 'bottom-left' });
    L.cnt = lab(N.noCalls[q].toLocaleString('en-US'), lc - PKL, G - PKD - 0.1, { anchor: 'top-left', size: 12, weight: 600, color: 'ink' });
    L.accL = lab('accuracy', lc + 0.02, YACC, { anchor: 'top-right' });
    L.accV = lab(N.accuracy[q].toFixed(3), lc + 0.12, YACC, { anchor: 'top-left', size: 12, weight: 600, color: 'ink' });
    let kp = 0, kh = 0;
    for (let i = 0; i < NC; i++) {
      const grey = GREY.includes(i), h = held.has(i);
      const c = { grey, held: h, s: S0 + i * DS, x0: lc + JX[i], k: h ? kh++ : kp++, tx: 0, ty: 0, T: 0, col: grey ? 'faint' : 'ink' };
      c.d = add(ctx.dot([c.x0, Y0, 5], { px: 3, color: c.col, opacity: 0 }));
      if (grey && !h) c.ring = add(ctx.dot([c.x0, Y0, 5.2], { px: 5.5, color: 'accent', hollow: true, ring: 0.42, opacity: 0 }));
      L.cells.push(c);
    }
    L.heldCells = L.cells.filter((c) => c.held);
    return L;
  });

  /* ── legend (top left), then the strip: “still called” meter and answer key ── */
  const legDot = add(ctx.dot([XL + 0.12, YLEG, 5], { px: 3, color: 'faint' }));
  lab('out of scope:\nNK, ILC, erythroid, progenitor', XL + 0.4, YLEG, { anchor: 'top-left' });
  add(ctx.line([[XL, YSEP], [XR, YSEP]], { color: 'line', width: 1 }));
  const icoDot = add(ctx.dot([XL + 0.12, YB, 5], { px: 3, color: 'faint' }));
  const icoRing = add(ctx.dot([XL + 0.12, YB, 5.2], { px: 5.5, color: 'accent', hollow: true, ring: 0.42 }));
  lab('still called', XL + 0.4, YB, { anchor: 'left' });
  add(ctx.line([[XT0, YB, 1], [XT1, YB, 1]], { color: 'faint', width: 1.5 }));
  add(ctx.line([[XT0, YB - 0.13, 1], [XT0, YB + 0.13, 1]], { color: 'faint', width: 1.5 }));
  add(ctx.line([[XT1, YB - 0.13, 1], [XT1, YB + 0.13, 1]], { color: 'faint', width: 1.5 }));
  const X0 = lerp(XT0, XT1, N.oosCalled[0] / 100), X1 = lerp(XT0, XT1, N.oosCalled[1] / 100);
  const band = add(ctx.line([[X0, YB, 2], [X1, YB, 2]], { color: 'accent', width: 9, opacity: 0 }));
  const bandLab = lab(`${N.oosCalled[0]}–${N.oosCalled[1]}%`, X0 - 0.1, YB,
    { anchor: 'right', size: 12, weight: 600, color: 'ink', bg: 'card', bgOpacity: 1, pad: 2, opacity: 0 });
  lab('answer key', XL, YC, { anchor: 'left' });
  const chips = [['B', 'b', false], ['T', 't', false], ['myeloid', 'm', false], ['NK', 'muted', true], ['out of scope', 'muted', true]]
    .map(([t, col, miss]) => ({
      t, miss,
      box: add(ctx.box(1, 0.5, { color: null, stroke: col, strokeWidth: 1.3, dashed: miss ? [3, 3] : false, radius: 0.12 })),
      l: lab(t, 0, YC, { color: col, weight: 600 }),
    }));
  const missLab = lab('needs', 0, YC, { anchor: 'right', color: 'warn', weight: 600, opacity: 0 }); // page: “the key needs NK and out-of-scope classes”
  const pointer = add(ctx.arrow([0, 0, 3], [1, 0, 3], { color: 'muted', width: 1.4, head: 7 }));

  /* ── text measuring for the chip row (label size × textScale, page font) ── */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'system-ui, sans-serif';
  const tw = (text, size, weight) => {
    if (!mc) return text.length * size * 0.55 * ctx.textScale;
    mc.font = `${weight} ${size * ctx.textScale}px ${family}`;
    return mc.measureText(text).width;
  };

  const P = { ppu: 27, rw: 0.11, lw: 0.07 };
  function layout() {
    const w = ctx.width || 400, h = ctx.height || 300;
    const ppu = (P.ppu = Math.min(w / W, h / H));
    // show the whole 12 x 9 frame (ppu above assumes it): in the wide 21:9 stage the page gives this card
    // between 681 and 1040 px, the default camera fit filled the width and cut the lane names and answer key
    const a = w / h, hh = a > W / H ? H / 2 : W / 2 / a;
    camera.left = -hh * a; camera.right = hh * a; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    const px = clamp(0.11 * ppu, 2.6, 5);
    P.rw = px / ppu; P.lw = 1.75 / ppu;
    lanes.forEach((L) => {
      L.cells.forEach((c) => {
        c.d.setPx(px);
        if (c.ring) c.ring.setPx(px + 2.4);
        if (c.held) {
          c.tx = L.lc - PKL + 0.22 + (c.k % 2) * 0.4;
          c.ty = G - PKD + P.rw + 0.05 + Math.floor(c.k / 2) * (2 * P.rw + 0.07);
        } else {
          c.tx = L.lc + ((c.k % 5) - 2) * DP;
          c.ty = PB + Math.floor(c.k / 5) * DP;
          c.T = c.s + Math.sqrt((Y0 - c.ty) / GR);
        }
      });
    });
    legDot.setPx(px); icoDot.setPx(px); icoRing.setPx(px + 2.4);
    legDot.position.y = YLEG - (11 * ctx.textScale * 1.22 * 0.5 + 1) / ppu; // centred on the legend's first line
    /* chip row: present classes after the label, missing ones flush right under the meter band */
    const u = (v) => v / ppu, padX = u(5), chipH = u(16), gap = u(4);
    let x = XL + u(tw('answer key', 11, 500)) + u(8);
    const place = (c, x0) => {
      const cw = u(tw(c.t, 11, 600)) + 2 * padX;
      c.box.setSize(cw, chipH, Math.min(0.12, chipH / 2));
      c.box.position.set(x0 + cw / 2, YC, 1);
      c.l.position.set(x0 + cw / 2, YC, 10);
      return cw;
    };
    chips.filter((c) => !c.miss).forEach((c) => { x += place(c, x) + gap; });
    let xr = XR;
    const miss = chips.filter((c) => c.miss);
    const wm = miss.map((c) => u(tw(c.t, 11, 600)) + 2 * padX);
    for (let i = miss.length - 1; i >= 0; i--) { xr -= wm[i]; place(miss[i], xr); if (i) xr -= gap; }
    missLab.position.set(xr - u(6), YC, 10);
    P.miss = xr - u(6) - u(tw('needs', 11, 600)) > x + u(4); // drop the word, keep the dashed chips, if the row is too tight
    const px0 = clamp((xr + XR) / 2, X0 + 0.1, X1 - 0.1); // straight down from the band onto the missing classes
    pointer.set([px0, YB - u(7), 3], [px0, YC + chipH / 2 + u(3), 3]);
  }

  function update(t) {
    const u = ctx.loopT(t, PERIOD);
    const fo = 1 - seg(u, 0.93, 0.99, ease.inOutSine);
    const yRest = G + P.rw + P.lw;
    const glow = pulse(u, 0.74, 0.88, 0.04);
    lanes.forEach((L) => {
      /* the gate closes around each held cell */
      let cl = 0;
      for (const c of L.heldCells) { const a = c.s + F; cl = Math.max(cl, pulse(u, a - 0.01, a + 0.022, 0.007)); }
      const hs = (L.slot / 2) * (1 - cl);
      L.shL.setPoints([[L.lc - HB, G, 1], [L.lc - hs, G, 1]]);
      L.shR.setPoints([[L.lc + hs, G, 1], [L.lc + HB, G, 1]]);
      for (const c of L.cells) {
        const dt = u - c.s;
        const o = seg(dt, 0, 0.012) * fo;
        let x = c.x0, y = Y0 - GR * Math.max(0, dt) * Math.max(0, dt), col = c.grey ? 'faint' : 'ink';
        if (c.held) {
          const a = c.s + F;
          y = Math.max(y, yRest);
          if (u > a) {
            const k1 = seg(u, a + 0.002, a + 0.032, ease.inOutSine), k2 = seg(u, a + 0.032, a + 0.054, ease.linear);
            x = lerp(c.x0, L.lc - HB, k1);
            if (k2 > 0) { x = lerp(L.lc - HB, c.tx, ease.outSine(k2)); y = lerp(yRest, c.ty, ease.inQuad(k2)); }
          }
        } else {
          y = Math.max(y, c.ty);
          if (dt > F) {
            x = lerp(c.x0, c.tx, ease.inOutSine((u - c.s - F) / Math.max(1e-4, c.T - c.s - F)));
            if (!c.grey) col = 'accent';
          }
        }
        if (col !== c.col) { c.d.setColor(col); c.col = col; }
        c.d.position.set(x, y, 5);
        c.d.setOpacity(o);
        if (c.ring) {
          c.ring.position.set(x, y, 5.2);
          c.ring.setOpacity(seg(dt, F, F + 0.012) * fo);
          c.ring.setPx(P.ppu ? clamp(0.11 * P.ppu, 2.6, 5) + 2.4 + 2.2 * glow : 5.5);
        }
      }
      const k = seg(u, 0.68, 0.72) * fo, ka = seg(u, 0.7, 0.74) * fo;
      L.cnt.setOpacity(k);
      L.accL.setOpacity(ka); L.accV.setOpacity(ka);
    });
    /* the meter, then the pointer at the answer key’s missing classes */
    band.setProgress(seg(u, 0.73, 0.8, ease.outCubic)).setOpacity(0.55 * fo);
    bandLab.setOpacity(seg(u, 0.74, 0.78) * fo);
    icoRing.setPx(clamp(0.11 * P.ppu, 2.6, 5) + 2.4 + 2.2 * glow);
    pointer.setProgress(seg(u, 0.8, 0.85)).setOpacity(fo);
    missLab.setOpacity(P.miss ? seg(u, 0.83, 0.87) * fo : 0);
  }

  layout();
  let dead = false;
  if (document.fonts && document.fonts.ready) {
    document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });
  }
  return {
    scene, camera, period: PERIOD, still: 0.9 * PERIOD,
    update, resize: layout,
    dispose() { dead = true; lanes.length = 0; chips.length = 0; },
  };
}
