/*
 * fw-inventory-join — docs/anm-framework.html #ai-inventory (4.5 · AI workflows, exploratory).
 * The page: swapping two requests' warehouses should move the answer 14→9; the model stayed 9→9.
 * A request–return join, fixed in advance, restored 14→9; three controls stayed 9→9. The snapshot
 * (ALB-C, 10 units, rev 7) is unchanged by the swap. Only values the page states are drawn; C's
 * return is left blank because the page gives no value for it. 12 s seamless loop.
 * Step 1: A's and B's warehouse chips trade places; "required" draws 14→9, "model" stays 9→9 (red).
 * Step 2: each return token docks with its own request; "with join" draws 14→9, controls stay 9→9.
 * The caveat sits bottom-right because the loader's play/pause button owns the top-right corner.
 */
export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const add = (o) => { scene.add(o); return o; };
  const fadeBox = (b, o, tint = 1) => { b.setOpacity(o * tint); if (b.outline) b.outline.setOpacity(o); };
  const put = (o, x, y, z = 5) => { o.position.set(x, y, z); return o; };

  /* snapshot: never changes */
  put(add(ctx.label('snapshot · unchanged', { size: 12, color: 'muted', anchor: 'left' })), -7.6, 3.95);
  put(add(ctx.box(6.8, 0.8, { color: null, stroke: 'ink', strokeWidth: 1.2, radius: 0.14 })), -4.2, 3.1, 0);
  put(add(ctx.label('ALB-C · 10 units · rev 7', { size: 13, weight: 600 })), -4.2, 3.1);
  /* caveat bottom-right: the loader's play/pause button owns the top-right corner */
  put(add(ctx.label('exploratory · one known case\ncalibration missed', { size: 11, color: 'warn', anchor: 'bottom-right', align: 'right' })), 7.7, -4.3);

  /* request cards A, B, C */
  const CX = -6.1, CHX = -5.65;
  const Y = { A: 1.35, B: -0.2, C: -1.75 };
  for (const id of ['A', 'B', 'C']) {
    put(add(ctx.box(2.9, 1.0, { color: 'soft', stroke: 'line', strokeWidth: 1.2, radius: 0.16 })), CX, Y[id], 0);
    put(add(ctx.label(id, { size: 14, weight: 700, color: id === 'C' ? 'faint' : 'ink' })), CX - 1.05, Y[id]);
  }

  /* warehouse chips on A and B (the two that are swapped) */
  const chip = (text, target) => {
    const g = add(ctx.group());
    const b = ctx.box(1.75, 0.58, { color: target ? 'accent' : null, opacity: 0.14, stroke: target ? 'accent' : 'faint', strokeWidth: 1.4, radius: 0.12 });
    const l = ctx.label(text, { size: 12, weight: 600, color: target ? 'accent' : 'muted' });
    l.position.z = 5;
    g.add(b, l);
    return { g, b, l, tint: target ? 0.14 : 1 };
  };
  const chipA = chip('ALB-C', true), chipB = chip('other', false);

  /* join bands (behind cards), one per request */
  const bands = ['A', 'B', 'C'].map((id) => put(add(ctx.box(5.4, 1.3, { color: 'accent', opacity: 0.07, stroke: 'accent', strokeWidth: 1.2, radius: 0.2 })), -5.05, Y[id], -1));

  /* loose return tokens, arriving B, A, C; linked to requests only by id tags */
  const DOCK = -3.95;
  const tokens = [
    { id: 'B', val: '9', x0: -1.5, y0: Y.A, t0: 0.52, yFirst: true },
    { id: 'A', val: '14', x0: -1.9, y0: Y.B, t0: 0.46, yFirst: false },
    { id: 'C', val: '', x0: -1.6, y0: Y.C, t0: 0.58, yFirst: false },
  ].map((k, i) => {
    const g = add(ctx.group());
    const b = ctx.box(1.0, 0.78, { color: 'card', stroke: k.id === 'C' ? 'faint' : 'muted', strokeWidth: 1.4, radius: 0.36 });
    const v = k.val ? ctx.label(k.val, { size: 13, weight: 600 }) : null;
    if (v) { v.position.z = 5; g.add(v); }
    const tag = ctx.label('id ' + k.id, { size: 11, color: k.id === 'C' ? 'faint' : 'muted', anchor: 'left' });
    tag.position.set(0.65, 0, 5);
    g.add(b, tag);
    g.position.z = 1;
    return Object.assign(k, { g, b, v, tag, arrive: 0.02 + i * 0.03, y1: Y[k.id] });
  });

  const byId = Object.fromEntries(tokens.map((k) => [k.id, k]));

  /* answer rows (right): value before the swap → value after it */
  const VX0 = 4.2, VX1 = 6.5, AX0 = 4.65, AX1 = 6.05;
  put(add(ctx.label('before', { size: 11, color: 'muted' })), VX0, 2.35);
  put(add(ctx.label('after', { size: 11, color: 'muted' })), VX1, 2.35);
  const rows = [
    { name: 'required', y: 1.35, from: '14', to: '9', color: 'ink', t0: 0.27 },
    { name: 'model', y: 0.05, from: '9', to: '9', color: 'bad', t0: 0.33, dimAt: 0.62 },
    { name: 'with join', y: -1.25, from: '14', to: '9', color: 'good', t0: 0.66 },
    { name: '3 controls', y: -2.55, from: '9', to: '9', color: 'faint', t0: 0.72, dashed: true },
  ].map((r) => {
    const drop = r.from !== r.to ? 0.17 : 0;
    const lab = put(add(ctx.label(r.name, { size: 12, weight: 600, color: r.color === 'faint' ? 'muted' : r.color, anchor: 'left' })), 0.9, r.y);
    const v0 = put(add(ctx.label(r.from, { size: 13, weight: 600 })), VX0, r.y);
    const v1 = put(add(ctx.label(r.to, { size: 13, weight: 600 })), VX1, r.y);
    const arrow = add(ctx.arrow([AX0, r.y + drop, 1], [AX1, r.y - drop, 1],
      { color: r.color, width: 2.2, head: 8, dashed: r.dashed ? [5, 4] : false }));
    return Object.assign(r, { lab, v0, v1, arrow });
  });

  /* step labels (bottom left), one at a time */
  const step1 = put(add(ctx.label('1 · swap warehouses', { size: 12, weight: 600, anchor: 'left' })), -7.6, -3.6);
  const step2 = put(add(ctx.label('2 · request–return join', { size: 12, weight: 600, color: 'accent', anchor: 'left' })), -7.6, -3.6);

  const placeChip = (c, y, x) => { c.g.position.set(x, y, 2); };

  return {
    scene, camera, period: PERIOD, still: 0.86 * PERIOD,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const out = 1 - ctx.seg(u, 0.91, 0.98, ease.inOutSine);

      /* step 1: A's and B's chips trade places, bowing apart so they never cross */
      const k = ctx.seg(u, 0.14, 0.27, ease.inOutCubic);
      const bow = Math.sin(Math.PI * k);
      placeChip(chipA, ctx.lerp(Y.A, Y.B, k), CHX + 1.9 * bow);
      placeChip(chipB, ctx.lerp(Y.B, Y.A, k), CHX - 0.4 * bow);
      const chipVis = ctx.seg(u, 0.0, 0.05, ease.inOutSine) * out;
      for (const c of [chipA, chipB]) { fadeBox(c.b, chipVis, c.tint); c.l.setOpacity(chipVis); }

      /* tokens arrive loose, then (step 2) dock with their own request */
      for (const tk of tokens) {
        const a = ctx.seg(u, tk.arrive, tk.arrive + 0.07, ease.outCubic);
        const s = ctx.seg(u, tk.t0, tk.t0 + 0.08, ease.linear);
        const sx = ease.inOutCubic(tk.yFirst ? (s - 0.4) / 0.6 : s / 0.6);
        const sy = ease.inOutCubic(tk.yFirst ? s / 0.6 : (s - 0.4) / 0.6);
        tk.g.position.x = ctx.lerp(tk.x0 + (1 - a) * 1.0, DOCK, sx);
        tk.g.position.y = ctx.lerp(tk.y0, tk.y1, sy);
        const o = a * out;
        fadeBox(tk.b, o);
        if (tk.v) tk.v.setOpacity(o);
        tk.tag.setOpacity(o * 0.95);
      }
      bands.forEach((b, i) => { const t0 = byId['ABC'[i]].t0; fadeBox(b, ctx.seg(u, t0 + 0.07, t0 + 0.11) * out, 0.07); });

      /* answer rows draw on in order; the step-1 model row dims once the join row is coming */
      for (const r of rows) {
        const show = ctx.seg(u, r.t0, r.t0 + 0.03, ease.inOutSine) * out;
        const o = show * (r.dimAt ? 1 - 0.35 * ctx.seg(u, r.dimAt, r.dimAt + 0.04) : 1);
        r.lab.setOpacity(o);
        r.v0.setOpacity(o);
        r.v1.setOpacity(o * ctx.seg(u, r.t0 + 0.04, r.t0 + 0.06));
        r.arrow.setProgress(ctx.seg(u, r.t0 + 0.01, r.t0 + 0.05, ease.inOutCubic));
        r.arrow.setOpacity(o);
      }

      step1.setOpacity(ctx.pulse(u, 0.1, 0.43, 0.03));
      step2.setOpacity(ctx.pulse(u, 0.44, 0.93, 0.03));
    },
  };
}
