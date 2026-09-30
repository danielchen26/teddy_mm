/*
 * Demo scene (not used on a page): four lineage clusters; NK and T drift together until they
 * nearly merge, a teal dashed ring marks the merged zone, then they part. 10 s seamless loop.
 * Positions are a fixed illustrative layout, not data.
 */
export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 10;
  const rnd = ctx.rand(7);
  const gauss = () => { let s = 0; for (let i = 0; i < 4; i++) s += rnd(); return (s - 2) / 1.15; };

  const clusters = [
    { key: 'b', name: 'B', home: [-5.3, 1.7], meet: [-5.3, 1.7] },
    { key: 'm', name: 'myeloid', home: [5.3, 1.9], meet: [5.3, 1.9] },
    { key: 'nk', name: 'NK', home: [-2.6, -1.4], meet: [-0.35, -1.25] },
    { key: 't', name: 'T', home: [2.6, -1.4], meet: [0.35, -1.25] },
  ].map((c) => {
    const g = ctx.group();
    const dots = [];
    for (let i = 0; i < 18; i++) {
      const off = [gauss() * 0.75, gauss() * 0.6];
      const d = ctx.dot([off[0], off[1], 0.1], { r: 0.13, color: c.key, opacity: 0.9 });
      d.userData.off = off;
      d.userData.phase = rnd();
      dots.push(d);
      g.add(d);
    }
    const lab = ctx.label(c.name, { size: 13, color: 'ink', weight: 600, anchor: 'bottom' });
    lab.position.set(0, 1.25, 1);
    g.add(lab);
    scene.add(g);
    return Object.assign(c, { g, dots, lab });
  });

  const ringPts = [];
  for (let i = 0; i < 64; i++) { const a = (i / 64) * Math.PI * 2; ringPts.push([Math.cos(a) * 2.1, Math.sin(a) * 1.35, 0]); }
  const ring = ctx.line(ringPts, { color: 'accent', width: 2, dashed: true, closed: true });
  ring.position.set(0, -1.25, 0.5);
  scene.add(ring);
  const ringLab = ctx.label('nearly merged', { size: 12, color: 'muted', anchor: 'top' });
  ringLab.position.set(0, -2.75, 1);
  scene.add(ringLab);

  const base = ctx.line([[-7.2, -3.6], [7.2, -3.6]], { color: 'line', width: 1 });
  scene.add(base);

  return {
    scene, camera, period: PERIOD, still: 4.6,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const k = ctx.seg(u, 0.1, 0.4, ease.inOutCubic) - ctx.seg(u, 0.62, 0.92, ease.inOutCubic);
      for (const c of clusters) {
        c.g.position.set(ctx.lerp(c.home[0], c.meet[0], k), ctx.lerp(c.home[1], c.meet[1], k), 0);
        for (const d of c.dots) {
          const w = 2 * Math.PI * (u * 2 + d.userData.phase);
          d.position.set(d.userData.off[0] + 0.05 * Math.cos(w), d.userData.off[1] + 0.05 * Math.sin(w), 0.1);
        }
      }
      const show = ctx.pulse(u, 0.36, 0.66, 0.06);
      ring.setOpacity(show);
      ring.setDashPhase(u, 8);
      ringLab.setOpacity(show);
    },
  };
}
