/*
 * idx-protocol-fair-test — index.html #protocol ("One fair test, reused in every experiment").
 * One tray of the same 9 predicted proteins is copied onto three lanes: the fixed rule (hard-wired
 * average), ANM (written question) and a trained classifier that also gets 32 embedding numbers
 * ("can only favour it"). Each lane gives calls; the answer key (measured proteins) stays sealed and
 * opens only at the end to score all three. help = ANM − fixed rule then reads "no better calls",
 * while ANM's checks and reasons stay lit. Layout is illustrative; no values are drawn. 12 s loop.
 */
/* Every pipeline number this scene prints (stated on the page). Refresh after the official-preprocessing rerun. */
const NUMBERS = {
  panel: 9,              // predicted proteins every lane reads
  cells: 16750,          // held-out cells
  embeddingNumbers: 32,  // the classifier's extra input: the first 32 numbers of TEDDY's embedding
};

export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const ROWS = [1.85, -0.3, -2.45], TY = ROWS[1], BH = 1.4, PH = 0.66, KY = 3.62, CY = -3.72;
  const add = (o) => { scene.add(o); return o; };
  const L = (text, o) => add(ctx.label(text, Object.assign({ size: 12, color: 'ink' }, o)));
  const tint = (b, fill, stroke) => { b.setOpacity(fill); if (b.outline) b.outline.setOpacity(stroke); };
  const fadeAll = (g, k) => g.children.forEach((c) => { if (c.setOpacity) c.setOpacity(k); fadeAll(c, k); });
  const grid = (s, gap, z) => {
    const g = ctx.group();
    for (let i = 0; i < 9; i++) {
      const c = ctx.box(s, s, { color: 'teddy', radius: s * 0.2 });
      c.position.set(((i % 3) - 1) * gap, (1 - Math.floor(i / 3)) * gap, z);
      g.add(c);
    }
    return add(g);
  };

  // the evidence tray: the same 9 predicted proteins per cell
  const tray = add(ctx.box(1.6, 1.6, { color: 'soft', stroke: 'line', radius: 0.16 }));
  const trayCells = grid(0.36, 0.46, 0.2);
  const trayTop = L('same evidence', { weight: 600, anchor: 'bottom' });
  const trayBot = L(`${NUMBERS.panel} predicted\nproteins\n${NUMBERS.cells.toLocaleString('en-US')} cells`, { size: 11, color: 'muted', anchor: 'top' });

  // three deciders, each ending in its calls
  const lanes = [
    { c: 'teddy', name: 'fixed rule', sub: 'hard-wired average' },
    { c: 'accent', name: 'ANM', sub: 'written question' },
    { c: 'train', name: 'classifier', sub: 'trained on labels' },
  ].map((ln, i) => {
    ln.y = ROWS[i];
    ln.path = add(ctx.line([[0, 0], [1, 0]], { color: 'line', width: 1.5 }));
    ln.box = add(ctx.box(4.4, BH, { color: ln.c, stroke: ln.c, strokeWidth: 1.5, radius: 0.2 }));
    ln.nameL = L(ln.name, { size: 13, weight: 600 });
    ln.subL = L(ln.sub, { size: 11, color: i === 1 ? 'accent' : 'muted', weight: i === 1 ? 600 : 500 });
    ln.arrow = add(ctx.arrow([0, 0], [1, 0], { color: ln.c, width: 2, head: 8 }));
    ln.pill = add(ctx.box(1.5, PH, { color: ln.c, stroke: ln.c, strokeWidth: 1.5, radius: PH / 2 }));
    ln.calls = L('calls', { weight: 600 });
    ln.ring = add(ctx.box(1.9, PH + 0.3, { color: null, stroke: 'measured', strokeWidth: 1.2, radius: (PH + 0.2) / 2 }));
    ln.copy = grid(0.15, 0.2, 3);
    return ln;
  });

  // hard-wired averaging gear, a badge on the fixed rule's corner (8 teeth: a 2-tooth turn loops seamlessly)
  const gear = add(ctx.group());
  gear.add(ctx.dot([0, 0, 0], { r: 0.4, color: 'card' }));
  const cog = ctx.group();
  gear.add(cog);
  cog.add(ctx.dot([0, 0, 0.1], { r: 0.2, hollow: true, ring: 0.45, color: 'teddy' }));
  for (let k = 0; k < 8; k++) {
    const a = (k * Math.PI) / 4, th = ctx.box(0.1, 0.13, { color: 'teddy', radius: 0.02 });
    th.position.set(Math.cos(a) * 0.25, Math.sin(a) * 0.25, 0.1);
    th.rotation.z = a - Math.PI / 2;
    cog.add(th);
  }

  // the classifier's extra input: 32 numbers of TEDDY's embedding
  const cable = add(ctx.line([[0, 0], [1, 0]], { color: 'teddy', width: 1.2 }));
  const cableLab = L(`+${NUMBERS.embeddingNumbers} embedding numbers`, { size: 11, color: 'muted', anchor: 'top-left' });
  const favour = L('can only favour it', { size: 11, color: 'muted', anchor: 'left' });
  const spark = add(ctx.dot([0, 0, 0], { px: 3.5, color: 'teddy' }));

  // the answer key: measured proteins, sealed until the end
  const key = add(ctx.box(3.5, 0.84, { color: 'soft', stroke: 'measured', strokeWidth: 1.5, dashed: [5, 4], radius: 0.2 }));
  const keyLab = L('answer key', { weight: 600 });
  const lockBody = add(ctx.box(0.34, 0.26, { color: 'measured', radius: 0.05 }));
  const arc = [[-0.11, -0.06]];
  for (let k = 0; k <= 12; k++) arc.push([-0.11 * Math.cos((k * Math.PI) / 12), 0.11 * Math.sin((k * Math.PI) / 12)]);
  arc.push([0.11, -0.06]);
  const shackle = add(ctx.line(arc, { color: 'measured', width: 1.6 }));
  const sealed = L('sealed', { size: 11, color: 'muted', anchor: 'right' });
  const scores = L('scores all three', { size: 11, weight: 600, anchor: 'right' });
  const measured = L('measured proteins', { size: 11, color: 'muted' });

  // help = ANM − fixed rule, between lanes 1 and 2
  const minus = add(ctx.group());
  minus.add(ctx.dot([0, 0, 0], { r: 0.28, hollow: true, ring: 0.2, color: 'ink' }));
  minus.add(ctx.line([[-0.13, 0, 0.1], [0.13, 0, 0.1]], { color: 'ink', width: 2 }));
  const help1 = L('help = ANM − fixed rule', { size: 11, color: 'muted', anchor: 'bottom-left' });
  const help2 = L('no better calls', { size: 12, weight: 600, anchor: 'top-left' });
  const extras = L('checks · reasons', { weight: 600, color: 'accent', anchor: 'left' });

  let cablePts = [[0, 0], [1, 0]], keyX = 0, dashedNow = true;
  function place(aspect) {
    // true "contain" of the 16 × 9 design (set here, like idx-setup-mode-b-vs-a, so wide stages add width)
    const a = aspect > 0 ? aspect : 16 / 9, hh = a > 16 / 9 ? 4.5 : 8 / a, hw = hh * a;
    camera.userData.animFit = { width: 2 * hw, height: 2 * hh, fit: 'contain' };
    camera.left = -hw; camera.right = hw; camera.top = hh; camera.bottom = -hh;
    camera.updateProjectionMatrix();
    const sx = ctx.clamp(hw / 8, 1, 1.6), X = (x) => x * sx + (sx - 1) * 3;
    const bw = 4.4 * (1 + (sx - 1) * 0.3), pw = 1.5 * (1 + (sx - 1) * 0.3);
    const tx = X(-6.3), bus = X(-4.65), bx = X(-2.0), px = X(1.75);
    tray.position.set(tx, TY, 0);
    trayCells.position.set(tx, TY, 0);
    trayTop.position.set(tx, TY + 0.95, 5);
    trayBot.position.set(tx, TY - 0.95, 5);
    lanes.forEach((ln, i) => {
      const y = ln.y, left = bx - bw / 2, right = bx + bw / 2;
      ln.path.setPoints(i === 1 ? [[tx + 0.8, y], [left, y]] : [[tx + 0.8, TY], [bus, TY], [bus, y], [left, y]]);
      ln.travel = i === 1 ? [[tx, y], [left + 0.45, y]] : [[tx, TY], [bus, TY], [bus, y], [left + 0.45, y]];
      ln.box.setSize(bw, BH, 0.2);
      ln.box.position.set(bx, y, 0);
      ln.nameL.position.set(bx, y + 0.24, 5);
      ln.subL.position.set(bx, y - 0.3, 5);
      ln.arrow.set([right + 0.1, y, 1], [px - pw / 2 - 0.12, y, 1]);
      ln.pill.setSize(pw, PH, PH / 2);
      ln.pill.position.set(px, y, 0);
      ln.calls.position.set(px, y, 5);
      ln.ring.setSize(pw + 0.26, PH + 0.2, (PH + 0.2) / 2);
      ln.ring.position.set(px, y, 0.5);
    });
    gear.position.set(bx - bw / 2 + 0.06, ROWS[0] + BH / 2 - 0.06, 2);
    const up = bx - bw * 0.2;
    cablePts = [[X(-7.8), CY], [up, CY], [up, ROWS[2] - BH / 2]];
    cable.setPoints(cablePts);
    cableLab.position.set(X(-7.8), CY - 0.16, 5);
    favour.position.set(up + 0.3, CY, 5);
    // key left of the loader's corner button; its state reads to the key's left
    keyX = px - 0.2;
    key.position.set(keyX, KY, 0);
    keyLab.position.set(keyX + 0.2, KY, 5);
    lockBody.position.set(keyX - 1.25, KY - 0.07, 1);
    sealed.position.set(keyX - 2.05, KY, 5);
    scores.position.set(keyX - 2.05, KY, 5);
    measured.position.set(keyX, 2.8, 5);
    // narrow stages: text starts right of the scored rings, so the two lines never touch them
    const my = (ROWS[0] + ROWS[1]) / 2, hx = ctx.width / (2 * hw) < 40 ? px + pw / 2 + 0.25 : px + 0.45;
    minus.position.set(px, my, 2);
    help1.position.set(hx, my + 0.03, 5);
    help2.position.set(hx, my - 0.03, 5);
    extras.position.set(px + pw / 2 + 0.38, ROWS[1], 5);
  }
  place(ctx.aspect);

  return {
    scene, camera, period: PERIOD, still: 0.88 * PERIOD,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const F = 1 - ctx.seg(u, 0.93, 0.99, ease.inOutSine);
      const travel = ctx.seg(u, 0.06, 0.27, ease.inOutSine), carry = ctx.pulse(u, 0.06, 0.3, 0.03);
      const work = ctx.pulse(u, 0.27, 0.42, 0.05), out = ctx.seg(u, 0.42, 0.48) * F;
      lanes.forEach((ln, i) => {
        ln.copy.position.copy(ctx.along(ln.travel, travel));
        ln.copy.position.z = 3;
        fadeAll(ln.copy, carry);
        tint(ln.box, 0.08 + 0.14 * work, 1);
        ln.arrow.setProgress(ctx.seg(u, 0.36, 0.44));
        ln.arrow.setOpacity(F);
        tint(ln.pill, 0.05 + 0.13 * out, 0.35 + 0.65 * out);
        ln.calls.setOpacity(out);
        ln.ring.outline.setOpacity(ctx.seg(u, 0.62 + i * 0.035, 0.68 + i * 0.035) * F);
      });
      cog.rotation.z = -(Math.PI / 2) * ctx.seg(u, 0.26, 0.42, ease.inOutSine);
      spark.position.copy(ctx.along(cablePts, ctx.seg(u, 0.1, 0.28, ease.inOutSine)));
      spark.position.z = 3;
      spark.setOpacity(ctx.pulse(u, 0.1, 0.29, 0.03));
      // the key opens only after every lane has answered, then re-seals for the next loop
      const open = ctx.seg(u, 0.54, 0.6) * F;
      shackle.position.set(keyX - 1.25, KY + 0.06 + 0.13 * open, 1);
      if ((open < 0.5) !== dashedNow) { dashedNow = open < 0.5; key.outline.setDashed(dashedNow ? [5, 4] : false); }
      // one state word at a time: "sealed" leaves before "scores all three" arrives, and returns after it left
      sealed.setOpacity(1 - ctx.seg(u, 0.52, 0.555) + ctx.seg(u, 0.965, 0.995));
      scores.setOpacity(ctx.seg(u, 0.565, 0.6) * (1 - ctx.seg(u, 0.93, 0.96)));
      extras.setOpacity(ctx.seg(u, 0.46, 0.52) * F);
      const m = ctx.seg(u, 0.72, 0.77) * F;
      fadeAll(minus, m);
      help1.setOpacity(m);
      help2.setOpacity(ctx.seg(u, 0.78, 0.83) * F);
    },
    resize(w, h) { if (w > 0 && h > 0) place(w / h); },
    dispose() { lanes.length = 0; cablePts = []; },
  };
}
