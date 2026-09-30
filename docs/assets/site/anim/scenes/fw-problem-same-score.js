/*
 * fw-problem-same-score (anm-framework.html #problem, "A score says it failed, not where").
 * Two identical pipelines: three inputs → carried state (teal box, a cloud of points) → output → score.
 * Run A: input 2 enters the state and fades out there; its band of the state stays dark (lost in state).
 * Run B: input 2 crosses the state and lights it up, then bounces off the output (ignored at output).
 * Both score tiles flip to the same red cross. Then "?" marks appear at the input, the state and the
 * output, and a dial wobbles between "degraded" and "broken" without settling: the score alone cannot
 * say where it failed, or how badly. Illustrative only: a fixed layout, no values. 12 s seamless loop.
 */
export default function create(ctx) {
  const { THREE, ease, seg, pulse, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });
  const PERIOD = 12;
  const rnd = ctx.rand(23);
  const C = ctx.colors;
  const EM = -5.35, BOX = -1.6, BW = 3.8, BH = 1.9, PLUG = 2.0, TILE = 4.7, TS = 1.3;
  const boxL = BOX - BW / 2, boxR = BOX + BW / 2, face = PLUG - 0.2, start = EM + 0.25;
  const BAND = [0.6, 0, -0.6]; // inputs 1, 2, 3 (input 2 is the one followed)
  const SPEED = 28, E0 = 0.03, E1 = 0.56, GAP = 0.022; // world units per loop; emission window, spacing
  const NP = Math.floor((E1 - E0) / GAP) + 1;
  const N_CLOUD = 200;
  const own = [];
  const tmp = new THREE.Color();

  /* one draw call per point set: round points, per-point colour, alpha and size (world units, floored on small stages) */
  const uPx = { value: 1 }, uFloor = { value: 38 };
  const VS = `attribute vec3 aColor; attribute float aAlpha; attribute float aSize; uniform float uPx; uniform float uFloor;
varying vec3 vColor; varying float vAlpha;
void main() {
  vColor = aColor; vAlpha = aAlpha;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  gl_PointSize = max(aSize * uPx, aSize * uFloor);
}`;
  const FS = `varying vec3 vColor; varying float vAlpha;
void main() {
  vec2 c = gl_PointCoord * 2.0 - 1.0;
  float d = dot(c, c);
  if (d > 1.0 || vAlpha < 0.004) discard;
  gl_FragColor = vec4(vColor, vAlpha * (1.0 - smoothstep(0.7, 1.0, d)));
  #include <colorspace_fragment>
}`;
  function points(n, order) {
    const geo = new THREE.BufferGeometry();
    const pos = new Float32Array(n * 3), col = new Float32Array(n * 3), al = new Float32Array(n), sz = new Float32Array(n);
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('aColor', new THREE.BufferAttribute(col, 3));
    geo.setAttribute('aAlpha', new THREE.BufferAttribute(al, 1));
    geo.setAttribute('aSize', new THREE.BufferAttribute(sz, 1));
    const mat = new THREE.ShaderMaterial({ uniforms: { uPx, uFloor }, vertexShader: VS, fragmentShader: FS, transparent: true, depthTest: false, depthWrite: false });
    own.push(ctx.track(geo), ctx.track(mat));
    const obj = new THREE.Points(geo, mat);
    obj.frustumCulled = false;
    obj.renderOrder = order;
    scene.add(obj);
    return {
      set(i, x, y, z, c, a, s) {
        pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z;
        col[i * 3] = c.r; col[i * 3 + 1] = c.g; col[i * 3 + 2] = c.b;
        al[i] = a; sz[i] = s;
      },
      flush() { for (const k in geo.attributes) geo.attributes[k].needsUpdate = true; },
    };
  }
  const add = (o, x, y) => { o.position.set(x, y, 0); scene.add(o); return o; };

  [['inputs', EM], ['state', BOX], ['output', PLUG]].forEach(([s, x]) => add(ctx.label(s, { size: 12, color: 'muted' }), x, 3.75));
  /* 'score' sits under the first tile: the stage's top-right corner is kept clear for the play/pause button */
  add(ctx.label('score', { size: 12, color: 'muted', anchor: 'top' }), TILE, 2.0 - TS / 2 - 0.16);

  const rows = [{ y: 2.0, name: 'run A', lost: true }, { y: -1.3, name: 'run B', lost: false }].map((r) => {
    const y = r.y;
    add(ctx.label(r.name, { size: 12.5, weight: 600, anchor: 'left' }), -7.75, y);
    BAND.forEach((o, b) => {
      scene.add(ctx.line([[EM, y + o], [boxL, y + o]], { color: 'line', width: 1 }));
      add(ctx.dot([0, 0], { r: b === 1 ? 0.15 : 0.12, color: b === 1 ? 'ink' : 'muted', order: 2 }), EM, y + o);
    });
    add(ctx.box(BW, BH, { color: 'accent', opacity: 0.07, radius: 0.28, order: -1 }), BOX, y);
    add(ctx.box(BW, BH, { color: null, stroke: 'accent', strokeWidth: 1.5, opacity: 0.7, radius: 0.28 }), BOX, y);
    scene.add(ctx.line([[boxR, y], [face, y]], { color: 'faint', width: 1.5 }));
    add(ctx.box(0.36, 0.9, { color: 'soft', stroke: 'muted', radius: 0.08, order: 1 }), PLUG, y);
    scene.add(ctx.arrow([PLUG + 0.24, y], [TILE - TS / 2 - 0.1, y], { color: 'muted', width: 1.5, head: 8 }));

    const tile = add(ctx.group(), TILE, y);
    const red = ctx.box(TS, TS, { color: 'bad', opacity: 0.12, radius: 0.2, order: 2 });
    const redEdge = ctx.box(TS, TS, { color: null, stroke: 'bad', strokeWidth: 2, radius: 0.2, order: 3 });
    const dash = ctx.line([[-0.22, 0], [0.22, 0]], { color: 'faint', width: 2.5, order: 4 });
    const k = 0.3;
    const x1 = ctx.line([[-k, -k], [k, k]], { color: 'bad', width: 3.5, order: 4 });
    const x2 = ctx.line([[-k, k], [k, -k]], { color: 'bad', width: 3.5, order: 4 });
    tile.add(ctx.box(TS, TS, { color: 'card', stroke: 'line', radius: 0.2, order: 1 }), red, redEdge, dash, x1, x2);

    const note = add(ctx.label(r.lost ? 'lost in state' : 'ignored at output', { size: 12, weight: 600, color: 'bad', anchor: 'top' }),
      r.lost ? BOX : PLUG, y - BH / 2 - 0.16);
    const qs = [-4.4, BOX, PLUG].map((x) => add(ctx.label('?', { size: 16, weight: 700, color: 'ink', bg: 'card', bgOpacity: 1, pad: 5, order: 12 }), x, y));

    const cloud = points(N_CLOUD, 3), cp = [];
    for (let i = 0; i < N_CLOUD; i++) {
      const px = boxL + 0.2 + rnd() * (BW - 0.4), py = (rnd() * 2 - 1) * (BH / 2 - 0.17), pz = (rnd() * 2 - 1) * 0.7;
      const band = py > 0.27 ? 0 : py < -0.27 ? 2 : 1;
      cp.push({ x: px, y: py, z: pz, arrive: E0 + (px - start) / SPEED, lit: !(r.lost && band === 1), ph: rnd() });
    }
    return Object.assign(r, { tile, red, redEdge, dash, x1, x2, note, qs, cloud, cp, parts: points(3 * NP, 5) });
  });

  /* dial: degraded ↔ broken */
  const DX = 4.9, DY = -4.15, DR = 0.8;
  const dial = add(ctx.group(), DX, DY);
  const arcPts = [];
  for (let i = 0; i <= 40; i++) { const a = Math.PI * (1 - i / 40); arcPts.push([Math.cos(a) * DR, Math.sin(a) * DR]); }
  const dialParts = [
    ctx.line(arcPts, { color: 'faint', width: 2 }),
    ctx.dot([-DR, 0], { r: 0.08, color: 'warn', order: 2 }),
    ctx.dot([DR, 0], { r: 0.08, color: 'bad', order: 2 }),
    ctx.dot([0, 0], { r: 0.07, color: 'ink', order: 3 }),
  ];
  const needle = ctx.group();
  const needleLine = ctx.line([[0, 0], [0, DR * 0.84]], { color: 'ink', width: 2.5, order: 2 });
  needle.add(needleLine);
  const dLab = [ctx.label('degraded', { size: 12, color: 'muted', anchor: 'right' }), ctx.label('broken', { size: 12, color: 'muted', anchor: 'left' })];
  dLab[0].position.set(-DR - 0.18, 0.02, 0);
  dLab[1].position.set(DR + 0.18, 0.02, 0);
  dial.add(...dialParts, needle, ...dLab);
  dialParts.push(needleLine, ...dLab);

  return {
    scene, camera, period: PERIOD, still: 0.61 * PERIOD,
    update(t) {
      const u = ctx.loopT(t, PERIOD);
      const w = 2 * Math.PI * u;
      /* slow, seamless camera sway: the flat layout stays put (z = 0), only the cloud shows depth */
      const ay = 0.1 * Math.sin(w), ax = 0.05 * Math.cos(w);
      camera.position.set(100 * Math.sin(ay) * Math.cos(ax), 100 * Math.sin(ax), 100 * Math.cos(ay) * Math.cos(ax));
      camera.lookAt(0, 0, 0);
      uPx.value = ctx.ppu() * ctx.dpr;
      uFloor.value = 38 * ctx.dpr;

      const fadeCloud = 1 - seg(u, 0.8, 0.88);
      const flip = seg(u, 0.38, 0.46, ease.inOutSine) - seg(u, 0.88, 0.96, ease.inOutSine);
      const cross = flip >= 0.5;
      const qOut = 1 - seg(u, 0.84, 0.9);
      const dim = 1 - 0.45 * seg(u, 0.5, 0.58); // the trail dims once only the score is left to read
      for (const r of rows) {
        r.cp.forEach((p, i) => {
          const l = p.lit ? seg(u, p.arrive, p.arrive + 0.05, ease.outCubic) * fadeCloud : 0;
          const j = 2 * Math.PI * (u * 2 + p.ph);
          tmp.copy(C.faint).lerp(C.accent, l);
          r.cloud.set(i, p.x + 0.03 * Math.cos(j), r.y + p.y + 0.03 * Math.sin(j), p.z, tmp, 0.5 + 0.45 * l, 0.075 + 0.015 * l);
        });
        r.cloud.flush();

        let n = 0;
        for (let b = 0; b < 3; b++) {
          for (let i = 0; i < NP; i++) {
            const s = (u - (E0 + i * GAP)) * SPEED;
            let x = start + s, y = r.y + BAND[b];
            let a = s > 0 ? clamp(s / 0.5) : 0;
            let c = b === 1 ? C.ink : C.muted;
            if (b !== 1) { // inputs 1, 3: absorbed into the state
              a *= 1 - clamp((x - (boxL + 0.1)) / 1.3);
              c = tmp.copy(C.muted).lerp(C.accent, clamp((x - boxL) / 0.8));
            }
            else if (r.lost) a *= 1 - clamp((x - (boxL + 0.15)) / 1.2); // run A: fades out inside the state
            else if (x > face) { // run B: reaches the output and bounces off
              const q = (x - face) / 0.9;
              x = face - 0.55 * ease.outQuad(q);
              y += (i % 2 ? 0.32 : -0.32) * ease.outQuad(q);
              a *= 1 - clamp(q);
            }
            r.parts.set(n++, x, y, 0, c, a * dim, b === 1 ? 0.14 : 0.115);
          }
        }
        r.parts.flush();

        r.tile.scale.x = Math.max(0.02, Math.abs(Math.cos(Math.PI * flip)));
        r.red.visible = r.redEdge.visible = r.x1.visible = r.x2.visible = cross;
        r.dash.visible = !cross;
        r.note.setOpacity(pulse(u, r.lost ? 0.13 : 0.28, 0.86, 0.06));
        r.qs.forEach((q, i) => {
          const a0 = 0.5 + 0.03 * i;
          q.setOpacity(seg(u, a0, a0 + 0.04) * qOut);
          q.position.y = r.y + 0.25 * (1 - ease.outBack((u - a0) / 0.06));
        });
      }

      const dv = seg(u, 0.52, 0.58) * qOut;
      dialParts.forEach((o) => o.setOpacity(dv));
      needle.rotation.z = 0.55 * Math.sin(w * 3) + 0.3 * Math.sin(w * 7 + 0.7);
    },
    dispose() {
      own.forEach((o) => { try { o.dispose(); } catch (e) { /* core disposes too */ } });
      own.length = 0;
    },
  };
}
