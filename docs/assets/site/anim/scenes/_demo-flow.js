/*
 * Demo scene (not used on a page): the teddy_mm pipeline as a flow.
 * RNA → frozen TEDDY → our head → ANM (+ written question) → call.
 * Links draw on, a few tokens travel along the path, then everything fades for a seamless 12 s loop.
 */
export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;

  const W = 2.5, H = 1.35, Y = -0.6;
  const steps = [
    { x: -6.35, text: 'RNA', fill: 'soft', stroke: 'line', sub: 'one cell' },
    { x: -3.18, text: 'TEDDY', fill: 'teddy', stroke: 'teddy', sub: 'frozen', tint: 0.14 },
    { x: 0, text: 'our head', fill: 'soft', stroke: 'muted', sub: 'trained' },
    { x: 3.18, text: 'ANM', fill: 'accent', stroke: 'accent', sub: 'typed evidence', tint: 0.14 },
    { x: 6.35, text: 'call', fill: 'soft', stroke: 'ink', sub: 'call / no call' },
  ].map((s) => {
    const b = ctx.box(W, H, { color: s.fill, opacity: s.tint || 1, stroke: s.stroke, strokeWidth: 1.5, radius: 0.22 });
    b.position.set(s.x, Y, 0);
    b.userData.tint = s.tint || 1;
    const lab = ctx.label(s.text, { size: 14, weight: 600, color: 'ink' });
    lab.position.set(s.x, Y, 2);
    const sub = ctx.label(s.sub, { size: 11, color: 'muted', anchor: 'top' });
    sub.position.set(s.x, Y - H / 2 - 0.22, 2);
    scene.add(b, lab, sub);
    return Object.assign(s, { b, lab, sub });
  });

  const links = [];
  for (let i = 0; i < steps.length - 1; i++) {
    const a = ctx.arrow([steps[i].x + W / 2 + 0.08, Y, 1], [steps[i + 1].x - W / 2 - 0.1, Y, 1], { color: 'muted', width: 2, head: 8 });
    scene.add(a);
    links.push(a);
  }

  const q = ctx.box(2.3, 0.9, { color: null, stroke: 'accent', strokeWidth: 1.5, dashed: [5, 4], radius: 0.2 });
  q.position.set(3.18, 2.45, 0);
  const qLab = ctx.label('question', { size: 12, color: 'ink' });
  qLab.position.set(3.18, 2.45, 2);
  const qArrow = ctx.arrow([3.18, 1.95, 1], [3.18, Y + H / 2 + 0.1, 1], { color: 'accent', width: 2, head: 8, dashed: true });
  scene.add(q, qLab, qArrow);

  const path = [[steps[0].x, Y, 3], [steps[4].x, Y, 3]];
  const tokens = [0, 1, 2].map(() => ctx.dot([0, 0, 3], { px: 4.5, color: 'ink' }));
  tokens.forEach((d) => scene.add(d));

  return {
    scene, camera, period: PERIOD, still: 7.2,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const fade = 1 - ctx.seg(u, 0.9, 0.99, ease.inOutSine);
      const fadeIn = ctx.seg(u, 0.0, 0.05, ease.inOutSine);
      const vis = Math.min(fade, fadeIn);
      steps.forEach((s, i) => {
        const k = ctx.seg(u, 0.02 + i * 0.07, 0.1 + i * 0.07, ease.outCubic) * vis;
        s.b.setOpacity(k * s.b.userData.tint);
        if (s.b.outline) s.b.outline.setOpacity(k);
        s.lab.setOpacity(k);
        s.sub.setOpacity(k * 0.95);
      });
      links.forEach((a, i) => {
        a.setProgress(ctx.seg(u, 0.08 + i * 0.07, 0.15 + i * 0.07, ease.inOutCubic));
        a.setOpacity(vis);
      });
      const qk = ctx.seg(u, 0.36, 0.44) * vis;
      q.setOpacity(qk); qLab.setOpacity(qk);
      qArrow.setProgress(ctx.seg(u, 0.4, 0.48)); qArrow.setOpacity(vis);
      qArrow.setDashPhase(u, 14);
      tokens.forEach((d, i) => {
        const s = ctx.seg(u, 0.48 + i * 0.08, 0.78 + i * 0.04, ease.inOutSine);
        d.position.copy(ctx.along(path, s));
        d.position.z = 3;
        d.setOpacity(ctx.pulse(u, 0.48 + i * 0.08, 0.8 + i * 0.04, 0.03) * vis);
        d.setColor(s > 0.84 ? 'accent' : s > 0.18 ? 'teddy' : 'ink');
      });
    },
  };
}
