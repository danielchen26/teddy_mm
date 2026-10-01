/*
 * fw-delay-ring — anm-framework.html, #physics-ring ("Physics I: a delayed ring, measured vs theory").
 * 64 sites on a ring (seen at a slant); each sphere's height is its activity c_i(t). Every site is linked
 * to 2 neighbours on each side (degree four). Beads carry each value along every link and arrive a delay
 * τ later, so a site is pushed by where its neighbours were τ ago. The teal ribbon behind each sphere is
 * the last τ of its height history: the retained state. One front site and its four incoming links are
 * highlighted. Above: the drive (grey, dashed) and the ring mean (teal) over one drive cycle; the teal
 * disc at the centre bobs with the mean, which lags the drive. The lag is the linear mean-mode response
 * with the page's parameters (a = 1, K = 0.40, Ω = 0.35) at τ = 0.445, the middle of the tested band
 * τ ∈ [0.43, 0.46]: φ = arg(a − K e^(−iΩτ) + iΩ).
 * τ is drawn longer than to scale (1/16 of the cycle) so the beads can be seen; the site jitter stands in
 * for the noise term and is illustrative. Side panel: the analytic amplitude F/|…| and phase φ against τ,
 * the tested band τ ∈ [0.43, 0.46] and the page's bounds. No per-delay measured points are drawn (not in
 * the released data here). 10 s seamless loop = one drive cycle.
 */
export default function create(ctx) {
  const { THREE, ease } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 10, N = 64, BURSTS = 16, TAU = 1 / BURSTS, TP = Math.PI * 2;
  const A = 1, K = 0.4, F = 0.25, OM = 0.35, TAU0 = (0.43 + 0.46) / 2;
  const resp = (tau) => {
    const re = A - K * Math.cos(OM * tau), im = OM + K * Math.sin(OM * tau);
    return { amp: F / Math.hypot(re, im), ph: Math.atan2(im, re) };
  };
  const LAG = resp(TAU0).ph / TP; // fraction of one drive cycle
  const add = (o) => { scene.add(o); return o; };
  const seg2 = () => [[0, 0], [1, 0]];
  const lab = (text, x, y, o = {}) => { const l = ctx.label(text, Object.assign({ size: 11, color: 'ink' }, o)); l.position.set(x, y, 6); return add(l); };

  /* ── ring ── */
  const CX = -2.55, CY = -1.45, RX = 4.2, RY = 1.2, H = 0.66, ARCH = 0.22, RL = 0.26, RS = 6, LS = 6;
  const rnd = ctx.rand(11);
  const sites = [];
  for (let i = 0; i < N; i++) {
    const th = -Math.PI / 2 + (TP * i) / N;
    const d = add(ctx.dot([0, 0, 0], { r: 0.09, color: new THREE.Color() }));
    sites.push({ bx: CX + RX * Math.cos(th), by: CY + RY * Math.sin(th), depth: (1 - Math.sin(th)) / 2, d, p: [rnd(), rnd(), rnd()], y: 0 });
  }
  const heightAt = (s, u) => H * (Math.cos(TP * (u - LAG)) + 0.06 * Math.sin(TP * (2 * u + s.p[0]))
    + 0.05 * Math.sin(TP * (3 * u + s.p[1])) + 0.035 * Math.sin(TP * (5 * u + s.p[2])));
  const tintSites = () => sites.forEach((s) => s.d.material.color.copy(ctx.colors.ink).lerp(ctx.colors.card, 0.5 * (1 - s.depth)));
  tintSites();

  const ell = [];
  for (let k = 0; k < 96; k++) { const a = (k / 96) * TP; ell.push([CX + RX * Math.cos(a), CY + RY * Math.sin(a), 0]); }
  add(ctx.line(ell, { color: 'faint', width: 1, dashed: [3, 5], closed: true, opacity: 0.7 }));

  // Links i→i+1 (straight) and i→i+2 (arched): one LineSegments, re-posed every frame.
  const pairs = [];
  for (let i = 0; i < N; i++) for (const d of [1, 2]) pairs.push([i, (i + d) % N, d]);
  const linkArr = new Float32Array(pairs.length * LS * 6);
  const linkGeo = ctx.track(new THREE.BufferGeometry());
  linkGeo.setAttribute('position', new THREE.BufferAttribute(linkArr, 3));
  const linkMat = ctx.bind(ctx.track(new THREE.LineBasicMaterial({ transparent: true, opacity: 0.55, depthWrite: false })), 'faint');
  const links = add(new THREE.LineSegments(linkGeo, linkMat));
  links.position.z = 1.5; links.frustumCulled = false;
  const P = new THREE.Vector3(), Q = new THREE.Vector3(), M = new THREE.Matrix4();
  const bez = (a, b, d, s, out) => {
    const cx = (a.bx + b.bx) / 2, cy = (a.y + b.y) / 2 + (d === 2 ? 2 * ARCH : 0), r = 1 - s;
    return out.set(r * r * a.bx + 2 * r * s * cx + s * s * b.bx, r * r * a.y + 2 * r * s * cy + s * s * b.y, 0);
  };

  // One bead per directed link (256), all leaving together every τ and arriving τ later.
  const beadGeo = ctx.track(new THREE.CircleGeometry(1, 12));
  const beadMat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false })), 'accent');
  const beads = add(new THREE.InstancedMesh(beadGeo, beadMat, pairs.length * 2));
  beads.position.z = 1.9; beads.frustumCulled = false;

  // Ribbons: the last τ of each site's height, trailing to the left (older = further left).
  const ribArr = new Float32Array(N * (RS + 1) * 6);
  const ribGeo = ctx.track(new THREE.BufferGeometry());
  ribGeo.setAttribute('position', new THREE.BufferAttribute(ribArr, 3));
  const idx = [];
  for (let i = 0; i < N; i++) for (let k = 0; k < RS; k++) { const v = (i * (RS + 1) + k) * 2; idx.push(v, v + 1, v + 2, v + 1, v + 3, v + 2); }
  ribGeo.setIndex(idx);
  const ribMat = ctx.bind(ctx.track(new THREE.MeshBasicMaterial({ transparent: true, opacity: 0.17, depthWrite: false, side: THREE.DoubleSide })), 'accent');
  const ribs = add(new THREE.Mesh(ribGeo, ribMat));
  ribs.position.z = 1.7; ribs.frustumCulled = false;

  // The highlighted front site: its four incoming links, beads, ribbon and arrival ring.
  const F0 = sites[0], FN = [1, 2, N - 1, N - 2];
  const fLinks = FN.map(() => add(ctx.line([[0, 0], [1, 0]], { color: 'accent', width: 1.6, opacity: 0.85 })));
  const fBeads = FN.map(() => add(ctx.dot([0, 0, 0], { px: 3.6, color: 'accent' })));
  const fRib = add(ctx.line([[0, 0], [1, 0]], { color: 'accent', width: 3.5, opacity: 0.85 }));
  const fRing = add(ctx.dot([0, 0, 0], { r: 0.1, color: 'accent', hollow: true, ring: 0.25 }));
  const retLab = lab('retained state', 0, 0, { anchor: 'top-right' });
  const tauLab = lab('delay τ', 0, 0, { anchor: 'top-left' });
  lab('64 sites\n4 neighbours each', -7.5, 1.25, { anchor: 'left', color: 'muted' });

  // Ring mean: a teal disc at the centre on a faint vertical track.
  add(ctx.line([[CX, CY - H - 0.25, 0.4], [CX, CY + H + 0.25, 0.4]], { color: 'faint', width: 1, dashed: [2, 4] }));
  const disc = add(ctx.dot([CX, CY, 3], { r: 0.34, color: 'accent' }));
  disc.scale.y = 0.12;
  const discLab = lab('mean', 0, 0, { anchor: 'left' });

  /* ── drive vs mean, one cycle ── */
  const SX0 = -7.5, SW = 9.9, SY = 2.7, SA = 0.7, U0 = -0.3;
  const sx = (f) => SX0 + f * SW;
  const dPts = [], mPts = [];
  for (let k = 0; k <= 120; k++) {
    const f = k / 120, u = U0 + f;
    dPts.push([sx(f), SY + SA * Math.cos(TP * u), 1]);
    mPts.push([sx(f), SY + SA * Math.cos(TP * (u - LAG)), 1.1]);
  }
  add(ctx.line(dPts, { color: 'muted', width: 1.5, dashed: [5, 4], opacity: 0.8 }));
  add(ctx.line(mPts, { color: 'accent', width: 2.4 }));
  const xp = sx(-U0), xm = sx(-U0 + LAG), ya = SY + SA + 0.16;
  add(ctx.arrow([xp, ya, 1.2], [xm, ya, 1.2], { color: 'ink', width: 1.5, head: 7, bend: 0.55 }));
  const lagLab = lab('lag', (xp + xm) / 2, ya + 0.34, { anchor: 'bottom' });
  lab('drive', SX0, SY + 0.45, { anchor: 'bottom-left', color: 'muted' });
  lab('ring mean', sx(0.8 + LAG), SY - SA - 0.3, { anchor: 'top' });
  const head = add(ctx.line(seg2(), { color: 'faint', width: 1.2 })); // drive dot to mean dot: the gap now
  const hD = add(ctx.dot([0, 0, 1.3], { px: 3.2, color: 'muted' }));
  const hM = add(ctx.dot([0, 0, 1.35], { px: 4.2, color: 'accent' }));

  /* ── panel: analytic amplitude and phase against τ ── */
  const PL = 3.2, PR = 7.85, PX0 = 4.05, PX1 = 7.55, PW = PX1 - PX0, BAND = TAU0;
  const panel = add(ctx.box(PR - PL, 8, { color: 'soft', stroke: 'line', strokeWidth: 1, radius: 0.22 }));
  panel.position.set((PL + PR) / 2, 0, -1);
  const taus = Array.from({ length: 61 }, (_, k) => k / 60), rs = taus.map(resp), BX = PX0 + BAND * PW;
  const rows = [['amplitude', 'within 2.51×10⁻³ rel.', 'amp'], ['phase', 'within 0.74°', 'ph']].map(([name, bound, key]) => {
    const v = rs.map((r) => r[key]);
    return {
      v, lo: Math.min(...v), hi: Math.max(...v),
      axis: add(ctx.line(seg2(), { color: 'faint', width: 1 })),
      band: add(ctx.line(seg2(), { color: 'accent', width: 4, opacity: 0.4 })),
      curve: add(ctx.line(seg2(), { color: 'muted', width: 2 })),
      name: lab(name, PX0, 0, { anchor: 'bottom', color: 'muted' }), // turned into a y-axis title below
      bound: lab(bound, PL, 0, { anchor: 'top-left', weight: 600 }),
    };
  });
  const theoryLab = lab('theory', PX1, 0, { anchor: 'bottom-right', color: 'muted' });
  rows.forEach((r) => { r.name.material.rotation = Math.PI / 2; });
  const bandLab = lab('13 delays', BX, 0, { anchor: 'top' });
  const tauAx = lab('τ', PX1, 0, { anchor: 'top-right', color: 'muted' });
  // Stack the panel in CSS px so the text rows never collide and the corner play button stays clear.
  let laidOut = 0;
  function layout(ppu) {
    const g = 4 / ppu, pad = 7 / ppu, lh = (11 * ctx.textScale * 1.22 + 2) / ppu;
    // on a phone the label grows in world units: keep its top inside the 9-unit frame (it was cut at 360 px)
    lagLab.position.y = Math.min(ya + 0.34, 4.5 - lh - 2 / ppu);
    const PT = 4.5 - 40 / ppu, PB = -4.5 + 8 / ppu, ph = (PT - PB - 2 * pad - 3 * lh - 5 * g) / 2;
    panel.setSize(PR - PL, PT - PB, 0.22);
    panel.position.y = (PT + PB) / 2;
    let y = PT - pad;
    rows.forEach((r, n) => {
      r.bound.position.set(PL + pad, y, 6);
      const y1 = y - lh - g, y0 = y1 - ph;
      r.axis.setPoints([[PX0, y1, 0], [PX0, y0, 0], [PX1, y0, 0]]);
      r.band.setPoints([[BX, y0, 0.1], [BX, y1, 0.1]]).setWidth(Math.max(3, 0.03 * PW * ppu));
      r.curve.setPoints(taus.map((tau, k) => [PX0 + tau * PW, y0 + (0.15 + (0.7 * (r.v[k] - r.lo)) / (r.hi - r.lo)) * ph, 0.2]));
      r.name.position.set(PX0 - g, (y0 + y1) / 2, 6);
      if (n === 1) theoryLab.position.set(PX1 - g, y0 + g, 6);
      bandLab.position.y = tauAx.position.y = y0 - g;
      y = y0 - 2 * g;
    });
  }

  return {
    scene, camera, period: PERIOD, still: 1.0,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const ppu = ctx.ppu();
      const sr = Math.max(0.085, 2.3 / ppu);
      let sum = 0;
      for (const s of sites) {
        const h = heightAt(s, u);
        sum += h; s.y = s.by + h;
        s.d.position.set(s.bx, s.y, 2 + 2 * s.depth);
        s.d.setRadius(s === F0 ? sr * 1.35 : sr);
      }
      const hm = sum / N;

      // Links and beads: b = how far the current wave of beads has travelled (0 = sent, 1 = arrived).
      const b = ctx.loopT(u * BURSTS, 1);
      const fade = Math.min(ctx.seg(b, 0, 0.12, ease.inOutSine), 1 - ctx.seg(b, 0.85, 1, ease.inOutSine));
      beadMat.opacity = 0.55 * fade;
      const br = Math.max(0.035, 1.5 / ppu);
      let li = 0, bi = 0;
      for (const [i, j, d] of pairs) {
        const a = sites[i], c = sites[j];
        for (let k = 0; k < LS; k++) {
          bez(a, c, d, k / LS, P); bez(a, c, d, (k + 1) / LS, Q);
          linkArr[li++] = P.x; linkArr[li++] = P.y; linkArr[li++] = 0;
          linkArr[li++] = Q.x; linkArr[li++] = Q.y; linkArr[li++] = 0;
        }
        bez(a, c, d, b, P); beads.setMatrixAt(bi++, M.makeScale(br, br, 1).setPosition(P.x, P.y, 0));
        bez(c, a, d, b, P); beads.setMatrixAt(bi++, M.makeScale(br, br, 1).setPosition(P.x, P.y, 0));
      }
      linkGeo.attributes.position.needsUpdate = true;
      beads.instanceMatrix.needsUpdate = true;

      let ri = 0;
      const rt = Math.max(0.07, 2.2 / ppu) / 2;
      for (const s of sites) {
        for (let k = 0; k <= RS; k++) {
          const x = s.bx - (k / RS) * RL, y = s.by + heightAt(s, u - (k / RS) * TAU);
          ribArr[ri++] = x; ribArr[ri++] = y - rt; ribArr[ri++] = 0;
          ribArr[ri++] = x; ribArr[ri++] = y + rt; ribArr[ri++] = 0;
        }
      }
      ribGeo.attributes.position.needsUpdate = true;

      FN.forEach((j, n) => {
        const s = sites[j], d = j === 1 || j === N - 1 ? 1 : 2, pts = [];
        for (let k = 0; k <= LS; k++) { bez(s, F0, d, k / LS, P); pts.push([P.x, P.y, 1.6]); }
        fLinks[n].setPoints(pts);
        bez(s, F0, d, b, P);
        fBeads[n].position.set(P.x, P.y, 1.95);
        fBeads[n].setOpacity(fade);
      });
      const rp = [];
      for (let k = 0; k <= RS; k++) rp.push([F0.bx - (k / RS) * RL, F0.by + heightAt(F0, u - (k / RS) * TAU), 1.75]);
      fRib.setPoints(rp);
      const arr = ctx.clamp(b / 0.4);
      fRing.position.set(F0.bx, F0.y, 4.2);
      fRing.setRadius(sr * 1.35 + 0.22 * ease.outCubic(arr));
      fRing.setOpacity(0.7 * (1 - arr));
      retLab.position.set(F0.bx - 0.05, F0.y - 0.3, 6);
      tauLab.position.set(F0.bx + 0.18, F0.y - 0.3, 6);

      disc.position.y = CY + hm;
      discLab.position.set(CX + 0.48, CY + hm, 6);

      // Playhead over one drive cycle; it fades at the wrap so the loop has no jump.
      const f = ctx.loopT(u - U0, 1), x = sx(f);
      const op = Math.min(ctx.seg(f, 0, 0.03), 1 - ctx.seg(f, 0.97, 1));
      const yD = SY + SA * Math.cos(TP * u), yM = SY + SA * Math.cos(TP * (u - LAG));
      head.setPoints([[x, yD, 0.9], [x, yM, 0.9]]).setOpacity(0.9 * op);
      hD.position.set(x, yD, 1.3); hD.setOpacity(op);
      hM.position.set(x, yM, 1.35); hM.setOpacity(op);

      if (Math.abs(ppu - laidOut) > 0.01) { layout(ppu); laidOut = ppu; }
    },
    onTheme() { tintSites(); },
    dispose() {
      [linkGeo, linkMat, beadGeo, beadMat, ribGeo, ribMat].forEach((o) => o.dispose());
      beads.dispose();
    },
  };
}
