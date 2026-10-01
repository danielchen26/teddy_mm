/*
 * fw-mode-a-layers — anm-framework.html, #mode-a ("Next: TEDDY's layers as the steps").
 * Mode A (planned) treats TEDDY-G's 12 transformer layers as the steps. Two lanes, an NK cell and a
 * T cell, each carry a small sheet of gene tokens: the per-gene residual stream, the retained state.
 * Crossing a layer slab adds a small update to every token (h + f(h)), thin arcs flash between tokens
 * (attention couples genes) and the sheet briefly squeezes and restores (TEDDY-G blocks are post-norm:
 * LayerNorm closes each block). After every layer a dashed probe bar appears between the lanes: that
 * is where NK versus T will be compared. At the end the tokens are averaged into the embedding z
 * (mean of tokens) and probed once more. All bars are dashed, equal and carry no values, and the
 * lanes keep a fixed gap: no layer-wise result exists (v3 E5 stopped at its derivative gate; E5-M
 * tested only the layer-12 gene-mean). Updates and arc pairs are illustrative, not data.
 * The top-right corner stays empty for the play/pause button. 12 s seamless loop.
 */
/* Every pipeline number this scene prints (stated on the page). Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  layers: 12,   // TEDDY-G transformer layers (config.json)
};

export default function create(ctx) {
  const { THREE, ease, seg, lerp, clamp, smooth } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const rnd = ctx.rand(12);

  const L = NUMBERS.layers, SX0 = -4.4, SX1 = 2.8, DX = (SX1 - SX0) / (L - 1);
  const slabX = Array.from({ length: L }, (_, i) => SX0 + i * DX);
  const barX = slabX.map((x) => x + DX / 2);          // one probe after each layer
  const XS = -5.2, XZ = 5.6;                           // sheets start here; the embedding z sits here
  const WAY = [XS, ...slabX, XZ];                      // 14 way points, 13 moves
  const T0 = 0.05, T1 = 0.67;                          // travel window (fraction of the loop)
  const NT = 5, GAP = 0.34, YC = 1.4;                  // tokens per sheet, spacing, |y| of the lanes
  const HALF = ((NT - 1) * GAP) / 2 + 0.18;            // half-height of a sheet frame
  const PAIRS = [[0, 2], [1, 3], [2, 4], [0, 3], [1, 4], [0, 4], [0, 1], [3, 4], [1, 2]];

  // Sheet position: steady glide with a gentle slow-down at each slab (one layer = one step).
  function sheetX(u) {
    const n = WAY.length - 1;
    const p = clamp((u - T0) / (T1 - T0)) * n;
    const i = Math.min(n - 1, Math.floor(p)), f = p - i;
    return lerp(WAY[i], WAY[i + 1], 0.55 * f + 0.45 * ease.inOutSine(f));
  }
  function pick3() {
    const a = PAIRS.slice();
    for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
    return a.slice(0, 3);
  }
  // Attention arc between two tokens of one sheet, bowing sideways (side = +1 right, -1 left).
  function arcPts(a, b, side) {
    const out = [], bow = side * (0.2 + 0.32 * Math.abs(b[1] - a[1]));
    for (let i = 0; i <= 12; i++) {
      const s = i / 12, r = 1 - s;
      out.push([r * a[0] + s * b[0] + 2 * r * s * bow, r * a[1] + s * b[1], 2]);
    }
    return out;
  }
  function text(str, x, y, opts) {
    const s = ctx.label(str, Object.assign({ size: 13 }, opts));
    s.position.set(x, y, 5);
    scene.add(s);
    return s;
  }

  // TEDDY's 12 layer slabs, spanning both lanes.
  const SLAB_H = 2 * (YC + HALF) + 0.3;
  const slabs = slabX.map((x) => {
    const b = ctx.box(0.2, SLAB_H, { color: 'teddy', opacity: 0.1, stroke: 'teddy', strokeWidth: 1, radius: 0.08 });
    b.position.set(x, 0, 0);
    scene.add(b);
    return b;
  });

  // Static labels: legend top-left, equation top-centre, lane names left, axis and probe legend bottom.
  const lo = -SLAB_H / 2 - 0.12;
  text('dot = gene token', -7.6, 3.95, { color: 'muted', anchor: 'left' });
  text('arc = attention', -7.6, 3.35, { color: 'teddy', anchor: 'left' });
  text('each layer: h + f(h)', 0, 3.65, { color: 'ink', weight: 600 });
  text('NK cell', -7.6, YC, { size: 14, weight: 600, color: 'nk', anchor: 'left' });
  text('T cell', -7.6, -YC, { size: 14, weight: 600, color: 't', anchor: 'left' });
  text('layer 1', slabX[0] - 0.12, lo, { weight: 600, color: 'teddy', anchor: 'top-left' });
  text(`layer ${NUMBERS.layers}`, slabX[L - 1] + 0.12, lo, { weight: 600, color: 'teddy', anchor: 'top-right' });
  text('embedding z\nmean of tokens', XZ, lo, { color: 'ink', anchor: 'top', align: 'center' });
  text('NK vs T probe · first result', -7.1, -4.0, { weight: 600, color: 'accent', anchor: 'left' });
  scene.add(ctx.line([[-7.42, -3.74, 1], [-7.42, -4.26, 1]], { color: 'accent', width: 1.5, dashed: [3, 3] }));

  // Dashed probe bars: one after each layer, one on the embedding z. Equal length, no values.
  const BAR = YC - HALF - 0.1;
  const bars = barX.map((x) => {
    const b = ctx.line([[x, BAR, 1], [x, -BAR, 1]], { color: 'accent', width: 1.5, dashed: [3, 3] });
    scene.add(b);
    return b;
  });
  const zBar = ctx.line([[XZ, YC - 0.32, 1], [XZ, -YC + 0.32, 1]], { color: 'accent', width: 2, dashed: [4, 3] });
  scene.add(zBar);

  // Two sheets of gene tokens (the per-gene residual stream of one cell each).
  const sheets = [{ key: 'nk', y: YC }, { key: 't', y: -YC }].map((s) => {
    const frame = ctx.box(0.5, 2 * HALF, { color: null, stroke: s.key, strokeWidth: 1, radius: 0.14 });
    const dots = Array.from({ length: NT }, () => ctx.dot([0, 0, 3], { r: 0.1, color: s.key }));
    const ring = ctx.dot([XZ, s.y, 2.5], { r: 0.2, hollow: true, ring: 0.25, color: s.key, opacity: 0.55 });
    const z = ctx.dot([XZ, s.y, 3.2], { r: 0.2, color: s.key });
    const arcs = [0, 1, 2].map(() => ctx.line([[0, 0, 2], [0, 0.1, 2]], { color: 'teddy', width: 1.3 }));
    // Residual update f_l(h_l) per layer and token (illustrative), and which genes attend per layer.
    const D = slabX.map(() => Array.from({ length: NT }, () => [(rnd() - 0.5) * 0.1, (rnd() - 0.5) * 0.1]));
    const A = slabX.map(() => pick3());
    scene.add(frame, ring, z, ...dots, ...arcs);
    return Object.assign(s, { frame, dots, ring, z, arcs, D, A });
  });

  return {
    scene, camera, period: PERIOD,
    still: PERIOD * (T0 + (9 / 13) * (T1 - T0)),   // sheets inside layer 9, arcs lit, 8 bars drawn
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const vis = seg(u, 0, 0.04, ease.inOutSine) * (1 - seg(u, 0.93, 0.99, ease.inOutSine));
      const sx = sheetX(u);
      const m = seg(u, 0.63, 0.71);           // tokens averaged into z while arriving

      let li = 0;
      for (let i = 1; i < L; i++) if (Math.abs(sx - slabX[i]) < Math.abs(sx - slabX[li])) li = i;
      const near = 1 - smooth(Math.abs(sx - slabX[li]) / 0.3);
      // Post-norm: just after leaving a slab the sheet squeezes and restores (LayerNorm).
      const ln = 1 - smooth(Math.abs(sx - slabX[li] - 0.24) / 0.16);
      const squash = 1 - 0.3 * ln;
      slabs.forEach((b, i) => {
        const h = i === li ? near : 0;
        b.material.opacity = 0.08 + 0.14 * h;
        b.outline.setOpacity(0.35 + 0.5 * h);
      });

      for (const s of sheets) {
        const pos = [];
        for (let j = 0; j < NT; j++) {
          let ox = 0, oy = (j - (NT - 1) / 2) * GAP;
          for (let i = 0; i < L; i++) {        // h_{l+1} = h_l + f_l(h_l): updates add up layer by layer
            const k = smooth((sx - slabX[i] + 0.12) / 0.24);
            if (k <= 0) break;
            ox += k * s.D[i][j][0];
            oy += k * s.D[i][j][1];
          }
          const x = sx + ox * squash * (1 - m), y = s.y + oy * squash * (1 - m);
          s.dots[j].position.set(x, y, 3);
          s.dots[j].setOpacity(vis * (1 - m));
          pos.push([x, y]);
        }
        s.frame.position.set(sx, s.y, 1.5);
        s.frame.scale.set(1, lerp(squash, 0.2, m), 1);
        s.frame.setOpacity(vis * 0.5 * (1 - m));
        s.z.position.set(sx, s.y, 3.2);
        s.z.setOpacity(vis * m);
        const pairs = s.A[li];
        s.arcs.forEach((a, k) => {
          const [p, q] = pairs[k];
          a.setPoints(arcPts(pos[p], pos[q], k % 2 ? -1 : 1));
          a.setOpacity(vis * near * (1 - m) * 0.9);
        });
      }

      bars.forEach((b, i) => {
        b.setProgress(smooth((sx - barX[i]) / 0.25));
        b.setOpacity(vis * 0.85);
      });
      zBar.setProgress(seg(u, 0.71, 0.79));
      zBar.setOpacity(vis);
    },
    dispose() {
      // Geometries, materials and label textures all live in the scene; the core disposes them.
      slabs.length = 0; bars.length = 0; sheets.length = 0;
    },
  };
}
