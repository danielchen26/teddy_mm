/*
 * idx-q-where-to-trust — index.html #problem, GAP 4 card "Know where to trust it" (thumbnail, 16:9).
 * The problem only; idx-exp-scope-gate (#scope) shows how Experiment 4 answers it.
 * The card: "Which calls are safe to act on?" / "Its score margin is a graded confidence it doesn't report."
 * The block (#scope, Problem · TEDDY + fixed rule): "It reports a call, not how sure it is."
 * A queue of calls (B, T, myeloid: ring colour + letter) leaves the fixed rule and rides the belt to
 * "act on". Every call wears the same grey check badge, so a sure call and a shaky one look alike. Above the queue a follow-up basket
 * with room for only a few cells (three empty slots) hovers and swings back and forth. Three times it
 * lowers a probe onto a call ("sure?", "shaky?", "which?"), finds the same badge, shakes and moves on;
 * its slots stay empty. No ranking, no confidence, no gate here: that is the experiment's scene.
 * Which calls are probed and the lineage mix are illustrative; the scene prints no pipeline number.
 * 12 s seamless loop (each call crosses the belt once per loop); the still frame shows the "shaky?" probe.
 */

/* Every pipeline number this scene prints (none). Refresh after the official-preprocessing rerun. */
const NUMBERS = {};

export default function create(ctx) {
  const { THREE, ease, seg, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  void NUMBERS;

  /* ── timing and layout (world units, 16 × 9 frame) ── */
  const PERIOD = 12, STILL = 5.15;
  const N = 7, SP = 2.05, SPAN = N * SP, X0 = -7.1, V = SPAN / PERIOD;  // each call crosses once per loop
  const CY = -1.7, R = 0.6, BELT = CY - R - 0.06;                      // calls, belt
  const BDX = 0.5, BDY = 0.5, BR = 0.27;                             // badge offset and radius
  const EDGE = 5.75;                                                   // inner edges of the two end boxes
  const YB = 0.8;                                                     // basket height
  const W = 2.2;                                                       // one probe's length (s)
  const LIN = [['t', 'T'], ['b', 'B'], ['m', 'M'], ['t', 'T'], ['b', 'B'], ['t', 'T'], ['m', 'M']];

  const mod = (a, b) => ((a % b) + b) % b;
  const cellX = (j, tt) => X0 + mod(j * SP + V * tt, SPAN);
  /* three probes: which call, where the basket meets it, what it asks (left, then right, then middle) */
  const PROBES = [[1, -4.6, 'sure?'], [2, 1.6, 'shaky?'], [5, -1.8, 'which?']].map(([j, x, text]) => {
    const a = mod(x - X0 - SP * j, SPAN) / V;
    return { j, x, text, a, b: a + W };
  });

  const put = (o, x = 0, y = 0, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };
  const lab = (text, o = {}) => put(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const circle = (r, n = 36) => Array.from({ length: n }, (_, i) => [r * Math.cos((2 * Math.PI * i) / n), r * Math.sin((2 * Math.PI * i) / n)]);

  /* ── the belt between the fixed rule and "act on" ── */
  put(ctx.line([[-EDGE, BELT], [EDGE, BELT]], { color: 'line', width: 1.5 }));

  /* ── the calls: a lineage letter in its colour, and the same grey check badge on every one ── */
  const cells = LIN.map(([tok, letter]) => {
    const g = put(ctx.group(), 0, CY, 1);
    const badge = put(ctx.group(), BDX, BDY, 0.4, g);
    const c = {
      g, badge,
      fill: put(ctx.dot([0, 0], { r: R, color: tok, opacity: 0.14 }), 0, 0, 0, g),
      ring: put(ctx.dot([0, 0], { r: R, hollow: true, ring: 0.13, color: tok }), 0, 0, 0.1, g),
      letter: put(ctx.label(letter, { size: 12, weight: 700, color: 'ink' }), 0, 0, 0.2, g),   // ink: legible on dark T/B
      disc: put(ctx.dot([0, 0], { r: BR, color: 'faint' }), 0, 0, 0, badge),
      check: put(ctx.line([[-0.11, 0.0], [-0.03, -0.09], [0.12, 0.085]], { color: 'card', width: 1.6 }), 0, 0, 0.05, badge),
    };
    return c;
  });

  /* ── the two ends: the fixed rule emits calls, "act on" takes every one ── */
  const endBox = (x, stroke, tint) => {
    put(ctx.box(2.0, 1.7, { color: 'card', radius: 0.2 }), x, CY + 0.05, 3);
    if (tint) put(ctx.box(2.0, 1.7, { color: tint, opacity: 0.14, radius: 0.2 }), x, CY + 0.05, 3.05);
    put(ctx.box(2.0, 1.7, { color: null, stroke, strokeWidth: 1.5, radius: 0.2 }), x, CY + 0.05, 3.1);
  };
  endBox(-EDGE - 1.0, 'teddy', 'teddy');
  endBox(EDGE + 1.0, 'muted', null);
  const ruleLab = lab('fixed\nrule', { color: 'teddy', weight: 600 });
  const actLab = lab('act on', { weight: 600, size: 11.5 });
  const note1 = lab('same badge on every call', { anchor: 'top-left', weight: 600 });
  const note2 = lab('a call, not how sure it is', { anchor: 'top-left', color: 'muted', size: 11.5 });

  /* ── the follow-up basket: room for a few cells, all slots empty ── */
  const basket = put(ctx.group(), 0, YB, 4);
  const TW = 1.5, BW = 1.15, BH = 0.48;
  const trap = [[-TW, BH], [-BW, -BH], [BW, -BH], [TW, BH]];
  const bMat = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, side: THREE.DoubleSide }));
  ctx.bind(bMat, 'soft');
  basket.add(new THREE.Mesh(ctx.track(new THREE.ShapeGeometry(new THREE.Shape(trap.map(([x, y]) => new THREE.Vector2(x, y))))), bMat));
  put(ctx.line(trap, { color: 'ink', width: 1.8 }), 0, 0, 0.1, basket);
  put(ctx.line([[-TW - 0.12, BH], [TW + 0.12, BH]], { color: 'ink', width: 2.4 }), 0, 0, 0.12, basket);
  const slots = [-0.8, 0, 0.8].map((x) => put(ctx.line(circle(0.34), { color: 'faint', width: 1.2, closed: true }), x, 0.0, 0.2, basket));
  let slotTok = 'faint';
  const bLab = lab('follow-up', { anchor: 'bottom', weight: 600 });
  const qMark = lab('?', { size: 20, weight: 700, color: 'muted', anchor: 'left' });

  /* ── the probe it lowers onto a call, and what it asks ── */
  const probe = put(ctx.line([[0, 0], [0, -1]], { color: 'ink', width: 1.4, dashed: [3, 3] }), 0, 0, 2.5);
  const hit = put(ctx.dot([0, 0], { r: 0.46, hollow: true, ring: 0.15, color: 'ink' }), 0, 0, 2.6);
  const asks = PROBES.map((p) => lab(p.text, { anchor: 'right', weight: 600, size: 12.5 }));

  /* basket x at time tt: follows the probed call inside a window, eases to the next call in between */
  function basketX(tt) {
    tt = mod(tt, PERIOD);
    for (let i = 0; i < PROBES.length; i++) {
      const p = PROBES[i];
      const d = mod(tt - p.a, PERIOD);
      if (d <= W) return p.x + V * d + BDX;
    }
    for (let i = 0; i < PROBES.length; i++) {
      const p = PROBES[i], n = PROBES[(i + 1) % PROBES.length];
      const gap = mod(n.a - p.b, PERIOD), d = mod(tt - p.b, PERIOD);
      if (d <= gap) {
        const from = p.x + V * W + BDX, to = n.x - V * (gap - d) + BDX;  // the next call, unwrapped
        return lerp(from, to, ease.inOutSine(d / gap));
      }
    }
    return 0;
  }

  /* label widths in world units: letters grow with the calls on large stages */
  function layout() {
    const ppu = ctx.ppu(), ts = ctx.textScale;
    const lp = clamp((0.42 * ppu) / ts, 12, 24);
    cells.forEach((c) => c.letter.setSize(lp));
    ruleLab.position.set(-EDGE - 1.0, CY + 0.05, 5);
    actLab.position.set(EDGE + 1.0, CY + 0.05, 5);
    const lh = (px) => (px * ts * 1.22 + 2) / ppu;
    note1.position.set(-EDGE - 2.0, CY - 0.95, 5);
    note2.position.set(-EDGE - 2.0, CY - 0.95 - lh(12) - 0.02, 5);
    const cw = clamp(ppu * 0.03, 1.4, 3);
    cells.forEach((c) => c.check.setWidth(cw));
    slots.forEach((sl) => sl.setDashed(ppu > 40 ? [3, 3] : false));   // dashes only where they read
  }

  function update(t) {
    const tt = ctx.loopT(t, PERIOD) * PERIOD;

    /* calls ride the belt; letters fade under the end boxes (labels draw on top of meshes) */
    cells.forEach((c, j) => {
      const x = cellX(j, tt);
      c.g.position.x = x;
      c.g.visible = Math.abs(x) < EDGE + 1.2;   // hidden once fully inside an end box
      c.letter.setOpacity(clamp((x + EDGE - 0.4) / 0.35) * clamp((EDGE - 0.4 - x) / 0.35));
      c.badge.scale.setScalar(1);
    });

    /* the basket: hover bob, lean with its swing, settle on arrival, shake when it gives up */
    const bx = basketX(tt);
    const vx = (basketX(tt + 0.04) - basketX(tt - 0.04)) / 0.08;
    const bob = 0.07 * Math.sin((2 * Math.PI * tt * 3) / PERIOD);
    let tilt = -0.11 * clamp(vx / 4, -1, 1);
    let probeK = 0, hitK = 0, ask = -1, askK = 0, qK = 0.55;
    PROBES.forEach((p, i) => {
      const d = mod(tt - p.a, PERIOD);
      if (d > W + 0.6) return;
      if (d <= 0.9) tilt += 0.05 * Math.sin(2 * Math.PI * 1.8 * d) * Math.exp(-3 * d);
      const s = (d - (W - 0.55)) / 0.5;
      if (s > 0 && s < 1) tilt += 0.1 * Math.sin(4 * Math.PI * s) * Math.sin(Math.PI * s);
      if (d > W) return;
      probeK = seg(d, 0.2, 0.6, ease.outCubic) * (1 - seg(d, W - 0.65, W - 0.3, ease.inCubic));
      hitK = seg(d, 0.55, 0.75) * (1 - seg(d, W - 0.7, W - 0.45));
      ask = i;
      askK = seg(d, 0.45, 0.75) * (1 - seg(d, W - 0.35, W - 0.05));
      qK = 0.55 + 0.45 * seg(d, W - 0.6, W - 0.4) * (1 - seg(d, W - 0.1, W + 0.3));
      const c = cells[p.j];
      c.badge.scale.setScalar(1 + 0.3 * hitK);
    });
    const by = YB + bob;
    basket.position.set(bx, by, 4);
    basket.rotation.z = tilt;
    bLab.position.set(bx, by + BH + 0.2, 6);
    qMark.position.set(bx + TW + 0.4, by + 0.05, 6);
    qMark.setOpacity(qK);

    /* the probe drops from the basket to the call's badge */
    const y0 = by - BH - 0.06, y1 = CY + BDY + BR + 0.04;
    probe.setPoints([[bx, y0, 2.5], [bx, y1, 2.5]]);
    probe.setProgress(Math.max(probeK, 0.0001));
    probe.setOpacity(probeK > 0.001 ? 1 : 0);
    hit.position.set(bx, CY + BDY, 2.6);
    hit.setOpacity(hitK);
    asks.forEach((a, i) => {
      a.position.set(bx - 0.24, (y0 + y1) / 2 + 0.05, 6);
      a.setOpacity(i === ask ? askK : 0);
    });
    /* a slot lights as if to take the call, then stays empty */
    const tok = hitK > 0.5 ? 'ink' : 'faint';
    if (tok !== slotTok) { slots[0].setColor(tok); slotTok = tok; }
  }

  layout();
  return { scene, camera, period: PERIOD, still: STILL, update, resize: layout };
}
