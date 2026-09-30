/*
 * idx-q-embedding-enough — index.html #problem, GAP 6 card "Is the embedding enough?" (thumbnail at the top of
 * the card, 16:9). The problem only; Experiment 6 (#anm-loop, withdrawn) and anm-loop.html answer it.
 * The card: "Can it flag cells TEDDY's embedding nearly merges?" / TEDDY + fixed rule: "Some myeloid and T cells
 * sit at cosine 0.982. Withdrawn." The block (#anm-loop, Problem): "TEDDY's embedding puts some myeloid and T
 * cells at cosine ≈ 0.982. Any readout that trusts it gives both the same call, and nothing in TEDDY + fixed rule
 * flags that the embedding is insufficient there."
 * 1 Two cells whose measured protein differs (core colour = measured lineage, myeloid and T; "≠" between them)
 *   float from the left into TEDDY's embedding (their rings turn TEDDY orange; names hide in flight) and drift
 *   into almost the same point: "cosine 0.982" (page value). The faint dots are other cells; positions are illustrative, not data.
 * 2 The readout built on the embedding (our head's predicted proteins, read by the fixed rule) reads that point
 *   and stamps both cells with the same call: two identical orange marks (no lineage is named).
 * 3 Beside the two stamped cells a small flagpole stays empty: "no flag".
 * No veto, no measured-protein swap, no pulling apart: that is the experiment's scene.
 * A "withdrawn" tag sits in the bottom-left corner (the top-right corner stays clear for the play button).
 * 10 s seamless loop; the still frame shows all three beats.
 */

/* Every pipeline number this scene shows. On index.html: GAP 6 card "Some myeloid and T cells sit at cosine
 * 0.982", #anm-loop "cosine ≈ 0.982" (first run, withdrawn). Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  cosine: 0.982,
};

export default function create(ctx) {
  const { THREE, ease, seg, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 10;

  const put = (o, x = 0, y = 0, z = 0) => { o.position.set(x, y, z); scene.add(o); return o; };
  const lab = (s, x, y, o = {}) => put(ctx.label(s, Object.assign({ size: 12, color: 'muted' }, o)), x, y, 10);

  /* ── layout (world units, 16 × 9 frame) ── */
  const Y = 0;                                         // the row of the pair, the readout and the calls
  const SX = -6.8, SDY = 0.95;                         // start column (measured) and half its gap
  const PX = -2.5, PY = -0.15, PW = 4.6, PH = 5.1;     // TEDDY's embedding
  const PL = PX - PW / 2;                              // its left edge
  const MX = -2.3;                                     // where the pair ends up
  const HY = PY + PH / 2 + 0.15;                       // headings row
  const RX = 1.65, RW = 2.75, RH = 1.15;               // the readout
  const CX = 4.2, ROWS = [Y + 0.75, Y - 0.75];         // the two called cells
  const FX = 6.6, FB = Y - 0.7, FT = Y + 1.35;         // the flagpole: base, top
  const FW = 0.95, FH = 0.6;                           // the empty flag's outline

  /* text widths in CSS px, so headings stay inside the frame on a small card */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, px, wt = 600) => { if (!mc) return s.length * px * 0.58; mc.font = `${wt} ${px}px ${family}`; return mc.measureText(s).width; };

  /* ── headings ── */
  const mHead = lab('measured', SX, HY, { anchor: 'bottom', size: 11.5, weight: 600, color: 'ink' });
  lab('TEDDY’s embedding', PL + 0.05, HY, { anchor: 'bottom-left', size: 12, weight: 600, color: 'teddy' });

  /* ── TEDDY's embedding: a soft panel with a few other cells (illustrative) ── */
  put(ctx.box(PW, PH, { color: 'soft', radius: 0.22 }), PX, PY, 0);
  put(ctx.box(PW, PH, { color: null, stroke: 'teddy', strokeWidth: 1.5, radius: 0.22 }), PX, PY, 0.1);
  const BG = [[-4.3, 1.75, 1, 0.1], [-3.25, 1.95, 2, 0.6], [-1.25, 1.85, 1, 0.35], [-0.6, 1.1, 2, 0.8],
    [-4.35, 0.45, 1, 0.55], [-4.2, -0.9, 2, 0.2], [-4.3, -2.35, 1, 0.9], [-3.0, -2.45, 2, 0.45],
    [-1.4, -2.45, 1, 0.7], [-0.5, -0.85, 2, 0.15]];
  const bg = BG.map(([x, y, c, ph]) => ({ x, y, c, ph, d: put(ctx.dot([0, 0], { px: 3, color: 'faint', opacity: 0.55 }), x, y, 0.5) }));

  /* ── the two cells: core = measured lineage, ring = measured (ink) outside, TEDDY orange inside ── */
  const cells = [['m', 'myeloid', 1], ['t', 'T cell', -1]].map(([tok, name, side]) => ({
    tok, side,
    core: put(ctx.dot([0, 0], { px: 6, color: tok }), 0, 0, side > 0 ? 5 : 5.1),
    ringM: put(ctx.dot([0, 0], { px: 9, hollow: true, ring: 0.24, color: 'measured' }), 0, 0, 5.3),
    ringT: put(ctx.dot([0, 0], { px: 9, hollow: true, ring: 0.24, color: 'teddy' }), 0, 0, 5.35),
    name: lab(name, 0, 0, { anchor: side > 0 ? 'bottom' : 'top', weight: 600, color: tok }),
  }));
  const neq = lab('≠', SX, Y, { size: 17, weight: 500, color: 'muted' });
  const cosLab = lab(`cosine ${NUMBERS.cosine}`, MX, Y - 1, { anchor: 'top', size: 11.5, weight: 600, color: 'ink', bg: 'card', bgOpacity: 0.92, pad: 3 });

  /* ── the readout built on the embedding: our head + the fixed rule ── */
  const rFill = put(ctx.box(RW, RH, { color: 'teddy', opacity: 0.14, radius: 0.18 }), RX, Y, 1);
  put(ctx.box(RW, RH, { color: null, stroke: 'teddy', strokeWidth: 1.6, radius: 0.18 }), RX, Y, 1.1);
  lab('readout', RX, Y, { size: 13, weight: 600, color: 'teddy' });
  const rSub = lab('our head +\nfixed rule', RX, Y - RH / 2 - 0.1, { anchor: 'top', size: 11 });
  const read = put(ctx.arrow([0, 0], [1, 0], { color: 'muted', width: 1.6, head: 7 }), 0, 0, 2);
  const readDot = put(ctx.dot([0, 0], { px: 3.5, color: 'teddy', opacity: 0 }), 0, 0, 3);

  /* ── the calls: one row per cell, the same orange mark on both ── */
  const rows = cells.map((c, i) => ({
    y: ROWS[i],
    link: put(ctx.line([[0, 0], [1, 0]], { color: 'faint', width: 1.4 }), 0, 0, 0.8),
    core: put(ctx.dot([0, 0], { px: 5, color: c.tok }), CX, ROWS[i], 4),
    ring: put(ctx.dot([0, 0], { px: 7, hollow: true, ring: 0.26, color: 'measured' }), CX, ROWS[i], 4.1),
    tag: put(ctx.box(0.5, 0.5, { color: 'teddy', radius: 0.08 }), 0, ROWS[i], 4.2),
  }));
  const same = lab('same call', CX, ROWS[0] + 0.6, { anchor: 'bottom', size: 12.5, weight: 600, color: 'ink' });

  /* ── the flagpole: stands, stays empty ── */
  put(ctx.line([[FX - 0.35, FB], [FX + 0.35, FB]], { color: 'muted', width: 2 }), 0, 0, 1);
  put(ctx.line([[FX, FB], [FX, FT]], { color: 'muted', width: 2 }), 0, 0, 1);
  put(ctx.dot([0, 0], { px: 3, color: 'muted' }), FX, FT, 1.1);
  const ghost = put(ctx.line([[FX, FT - 0.08], [FX + FW, FT - 0.08], [FX + FW, FT - 0.08 - FH], [FX, FT - 0.08 - FH]],
    { color: 'faint', width: 1.3, dashed: [3, 3] }), 0, 0, 0.9);
  const noFlag = lab('no flag', FX, FB - 0.15, { anchor: 'top', weight: 600 });

  lab('withdrawn', -7.75, -4.2, { anchor: 'bottom-left', size: 10.5, mono: true, weight: 600, color: 'bad', bg: 'soft', bgOpacity: 1, pad: 4 });

  /* ── sizes that depend on the stage (px ↔ world) ── */
  const P = { ppu: 20, ring: 0.45, d: 0.17, g3: 0.15 };
  function layout() {
    const ppu = (P.ppu = ctx.ppu()), w = (px) => px / ppu, ts = ctx.textScale;
    const cp = clamp(0.2 * ppu, 5, 13), rp = cp * 1.5;
    P.ring = w(rp); P.d = 0.55 * w(cp); P.g3 = w(3);
    cells.forEach((c) => { c.core.setPx(cp); c.ringM.setPx(rp); c.ringT.setPx(rp); });
    bg.forEach((b) => b.d.setPx(clamp(0.07 * ppu, 2.2, 4)));
    /* "measured" centred on the column, but never past the left edge */
    const hw = w(textW('measured', 11.5 * ts) / 2 + 2);
    mHead.position.x = Math.max(SX, -7.9 + hw);
    /* the cosine chip sits under the T cell's label once the pair has merged */
    cosLab.position.set(MX, Y - P.ring - P.g3 - w(12 * ts * 1.22 + 2) - w(4), 10);
    rSub.position.y = Y - RH / 2 - w(5);
    read.set([MX + P.d + P.ring + w(6), Y, 2], [RX - RW / 2 - w(4), Y, 2]);
    const mcp = cp * 0.8, mr = w(rp * 0.8), s = w(clamp(0.55 * ppu, 10, 26));
    let right = CX;
    rows.forEach((r) => {
      r.core.setPx(mcp); r.ring.setPx(rp * 0.8);
      r.tag.setSize(s, s, s * 0.18);
      r.tag.position.x = CX + mr + w(5) + s / 2;
      right = r.tag.position.x + s / 2;
      r.link.setPoints([[RX + RW / 2 + w(3), Y, 0.8], [CX - mr - w(4), r.y, 0.8]]);
    });
    same.position.set((CX - mr + right) / 2, ROWS[0] + Math.max(mr, s / 2) + w(5), 10);
    noFlag.position.y = FB - w(5);
  }

  const bez = (a, c, b, k) => {
    const q = 1 - k;
    return [q * q * a[0] + 2 * q * k * c[0] + k * k * b[0], q * q * a[1] + 2 * q * k * c[1] + k * k * b[1]];
  };

  function update(t) {
    const u = ctx.loopT(t, PERIOD), TAU = 2 * Math.PI;
    const kIn = seg(u, 0, 0.04, ease.inOutSine), kOut = 1 - seg(u, 0.9, 0.96, ease.inOutSine), vis = kIn * kOut;

    /* other cells idle in place (integer cycles: seamless) */
    bg.forEach((b) => b.d.position.set(b.x + 0.05 * Math.sin(TAU * (u * b.c + b.ph)), b.y + 0.05 * Math.cos(TAU * (u * b.c + b.ph)), 0.5));

    /* 1 · measured apart (≠), then into TEDDY's embedding, where they drift into almost the same point */
    mHead.setOpacity(Math.max(1 - seg(u, 0.1, 0.16), seg(u, 0.95, 0.99)));
    neq.setOpacity(vis * (1 - seg(u, 0.07, 0.12)));
    const k = seg(u, 0.06, 0.34, ease.inOutCubic);
    cells.forEach((c) => {
      const [x, y] = bez([SX, Y + c.side * SDY], [-4.3, Y + c.side * 1.25], [MX - c.side * P.d, Y], k);
      const kT = clamp((x - (PL - 0.3)) / 0.6);
      c.core.position.set(x, y, c.core.position.z); c.core.setOpacity(vis);
      c.ringM.position.set(x, y, 5.3); c.ringM.setOpacity(vis * (1 - kT));
      c.ringT.position.set(x, y, 5.35); c.ringT.setOpacity(vis * kT);
      c.name.position.set(x, y + c.side * (P.ring + P.g3), 10);
      c.name.setOpacity(vis * Math.max(1 - seg(u, 0.07, 0.11), seg(u, 0.3, 0.34))); // hidden in flight
    });
    cosLab.setOpacity(seg(u, 0.33, 0.37) * kOut);

    /* 2 · the readout reads that one point and stamps both cells the same */
    read.setProgress(seg(u, 0.38, 0.44)).setOpacity(kOut);
    const kd = (u - 0.43) / 0.07;
    if (kd > 0 && kd < 1) {
      readDot.position.copy(ctx.along(read.line, ease.inOutSine(kd))).setZ(3);
      readDot.setOpacity(Math.min(1, kd / 0.15, (1 - kd) / 0.15) * kOut);
    } else readDot.setOpacity(0);
    rFill.setOpacity(0.14 + 0.22 * ctx.pulse(u, 0.49, 0.57, 0.03));
    rows.forEach((r, i) => {
      r.link.setProgress(seg(u, 0.53, 0.58)).setOpacity(kOut);
      const kc = seg(u, 0.55, 0.59) * kOut;
      r.core.setOpacity(kc); r.ring.setOpacity(kc);
      const kt = clamp((u - (0.59 + 0.035 * i)) / 0.04);
      r.tag.visible = kt > 0;
      r.tag.scale.set(Math.max(1e-3, ease.outBack(kt)), Math.max(1e-3, ease.outBack(kt)), 1);
      r.tag.setOpacity(kOut);
    });
    same.setOpacity(seg(u, 0.66, 0.7) * kOut);

    /* 3 · nothing flags it: the pole stays empty */
    const kf = seg(u, 0.72, 0.76) * kOut;
    ghost.setOpacity(0.9 * kf);
    noFlag.setOpacity(kf);
  }

  let dead = false;
  layout();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  return {
    scene, camera, period: PERIOD, still: 8.2,
    update, resize: layout,
    dispose() { dead = true; cells.length = 0; rows.length = 0; bg.length = 0; },
  };
}
