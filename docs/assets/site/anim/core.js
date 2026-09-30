/*
 * teddy_mm animation core — three.js r186 (vendored at docs/assets/vendor/three/).
 *
 * Pages never import this file directly. They include, in this order:
 *   <script type="importmap">{"imports":{"three":"./assets/vendor/three/three.module.min.js"}}</script>
 *   <link rel="stylesheet" href="assets/site/anim/anim.css">
 *   <script type="module" src="assets/site/anim/loader.js"></script>
 * and one placeholder per animated section or question:
 *   <figure class="anim" data-anim="NAME">
 *     <div class="anim-stage" role="img" aria-label="ALT TEXT"></div>
 *     <figcaption>ONE SHORT CAPTION</figcaption>
 *   </figure>
 * Optional: data-aspect="4:3" (or "1:1", "21:9", "3/2", "1.6") on the figure or the stage; 16:9 by default.
 *
 * loader.js finds the figures, mounts a scene when its stage comes near the viewport, disposes it
 * (renderer.dispose + forceContextLoss) when it is far away, keeps at most 6 live WebGL contexts,
 * pauses offscreen stages and hidden tabs, follows theme and reduced-motion changes, and toggles
 * pause on click/tap. This file builds one stage: the renderer, the ctx object, the helpers.
 *
 * ─── Scene module ────────────────────────────────────────────────────────────────────────────
 * docs/assets/site/anim/scenes/NAME.js   (NAME = kebab-case, unique site-wide)
 *
 *   export default function create(ctx) {
 *     const scene = new ctx.THREE.Scene();
 *     const camera = ctx.orthoCamera({ width: 16, height: 9 });   // or any THREE camera
 *     ...build objects with the helpers below...
 *     return {
 *       scene, camera,
 *       update(t, dt) { ... },   // t = seconds since start (paused time excluded), dt = seconds (0..0.1, may be 0)
 *       resize(w, h) { ... },    // optional; CSS px of the stage. Cameras made by ctx.orthoCamera refit on their own.
 *       onTheme() { ... },       // optional; after ctx.colors changed (token-coloured helpers already updated)
 *       dispose() { ... },       // optional; core also disposes every geometry/material/texture in the scene
 *       period: 10,              // optional; loop length in seconds (documentation + default still frame)
 *       still: 7.5,              // optional; t of the representative frame shown under reduced motion
 *     };
 *   }
 *
 * Write update() as a function of t (use ctx.loopT(t, period) for seamless loops): the core calls
 * update(t, 0) again after a theme change or resize while paused, and update(still, 0) for reduced motion.
 * Create objects once in create(); change them in update() (setText, setProgress, position, opacity...).
 *
 * ─── ctx API ─────────────────────────────────────────────────────────────────────────────────
 *   ctx.THREE                 the three.js module namespace (r186)
 *   ctx.VERSION               this core's version string
 *   ctx.name                  scene NAME
 *   ctx.figure, ctx.stage     the <figure class="anim"> and its .anim-stage element
 *   ctx.canvas, ctx.renderer  the WebGL canvas and THREE.WebGLRenderer (owned by the core; do not dispose)
 *   ctx.width, ctx.height     stage size in CSS px (live getters)
 *   ctx.aspect                width / height
 *   ctx.dpr                   device pixel ratio in use (capped at 2)
 *   ctx.theme                 'light' | 'dark' (from the card colour's luminance)
 *   ctx.reducedMotion         true when a single still frame is shown instead of the loop
 *   ctx.textScale             label size multiplier: clamp(width / 640, 0.82, 1.12)
 *
 *   ctx.colors                { ink, muted, faint, line, card, paper, accent, teddy, measured,
 *                               nk, t, b, m, good, bad, warn, soft, grid, train } as THREE.Color.
 *                             Read from the page's CSS variables (--anim-* tokens in anim.css, which map to
 *                             --ink, --ink-muted, --ink-faint, --line, --bg-elev, --bg, --s-anm, --s-teddy, …)
 *                             and updated IN PLACE on theme change. Treat as read-only; use ctx.color() for a copy.
 *                             accent = ANM teal, teddy = TEDDY orange, measured = ink (measured protein),
 *                             nk / t / b / m = NK, T, B, myeloid lineage colours.
 *   ctx.color(token)          -> new THREE.Color copy of a token
 *   ctx.hex(token)            -> '#rrggbb' string of a token
 *   ctx.bind(target, token, prop = 'color')
 *                             keep target coloured by token across theme changes. target = a material
 *                             (target[prop] is a THREE.Color), a uniform { value: THREE.Color }, or a THREE.Color.
 *                             Returns target. ctx.unbind(target) stops it.
 *
 *   Colour arguments of every helper: a token name ('ink', 'accent', 'teddy', 'nk', …) follows the theme
 *   automatically; anything else ('#hex', 0xhex, THREE.Color) is fixed.
 *
 *   ctx.label(text, { size = 13, color = 'ink', anchor = 'center', weight = 500, mono = false,
 *                     bg = null, bgOpacity = 0.92, pad = 4, align, opacity = 1, order = 10 })
 *       -> THREE.Sprite: a crisp text sprite. size is CSS px (× ctx.textScale), constant on screen whatever
 *          the camera; rasterised at device resolution with the page font (IBM Plex Sans / --mono).
 *          anchor: which point of the text sits at label.position: 'center' | 'left' | 'right' | 'top' |
 *          'bottom' | 'top-left' | 'top-right' | 'bottom-left' | 'bottom-right'. '\n' makes lines.
 *          bg: token for a rounded plate behind the text. Drawn on top (depthTest off, renderOrder = order).
 *          Methods: setText(s), setColor(c), setSize(px), setOpacity(o), setAnchor(a).
 *   ctx.line(points, { color = 'muted', width = 2, dashed = false, dash = 6, gap = 5, opacity = 1,
 *                      closed = false, order = 0 })
 *       -> THREE.Mesh: a screen-space polyline. points = [[x,y],[x,y,z],Vector3,{x,y,z},…] (≥ 2).
 *          width, dash, gap in CSS px (constant on screen). dashed: true or [dash, gap].
 *          Methods: setPoints(points), setProgress(p, p0 = 0) (draw only the part between fractions p0..p of
 *          the length), setDashOffset(px), setDashPhase(u, cycles = 1) (marching dashes that stay seamless:
 *          offset = u × cycles × (dash + gap); pass u = ctx.loopT(t, period)), setColor(c),
 *          setOpacity(o), setWidth(px), setDashed(bool | [dash, gap]). Property: length (world units).
 *   ctx.arrow(from, to, { color = 'muted', width = 2, head = 9, bend = 0, dashed = false, opacity = 1, order = 0 })
 *       -> THREE.Group (line + triangular head; head = length in CSS px). bend = sideways bow as a fraction
 *          of the length (+ = to the left of from→to). Methods: set(from, to, bend?), setProgress(p)
 *          (grows from `from`; the head rides the tip), setColor(c), setOpacity(o), setDashOffset(px),
 *          setDashPhase(u, cycles).
 *          Properties: line, head.
 *   ctx.dot(pos, { r = 0.12, px = null, color = 'ink', opacity = 1, hollow = false, ring = 0.3, order = 0 })
 *       -> THREE.Mesh: a flat disc. r in world units; or px = radius in CSS px (constant on screen).
 *          hollow: ring of thickness ring × r. Methods: setColor(c), setOpacity(o), setRadius(r), setPx(px).
 *   ctx.box(w, h, { color = 'soft', radius = 0.18, opacity = 1, stroke = null, strokeWidth = 1.5,
 *                   dashed = false, order = 0 })
 *       -> THREE.Mesh: a flat rounded rectangle centred on its position (world units). color = null → no fill.
 *          stroke: token/colour for an outline (strokeWidth CSS px). Methods: setColor(c), setStroke(c),
 *          setOpacity(o), setSize(w, h, radius?). Property: outline (the stroke line or null).
 *   ctx.group()               -> new THREE.Group()
 *   ctx.orthoCamera({ width = 16, height = 9, fit = 'contain', center = [0, 0], zoom = 1 })
 *       -> THREE.OrthographicCamera looking down -z at the XY plane that always shows the width × height
 *          design rectangle ('contain') or fills the stage with it ('cover'), refitted on resize.
 *          Layering: larger z draws on top (all helper materials are transparent, depthWrite off).
 *   ctx.v3(x, y, z = 0)       -> THREE.Vector3; ctx.toV3(point) converts any accepted point form.
 *   ctx.along(pointsOrLine, u, outTangent?) -> THREE.Vector3 at fraction u (0..1) of a polyline's length.
 *   ctx.ppu(pos?)             CSS px per world unit at pos (orthographic: the same everywhere).
 *   ctx.requestRender()       draw one frame even when paused/still (after changing something outside update).
 *   ctx.track(obj)            register an extra geometry/material/texture for disposal. Returns obj.
 *
 *   Timing and easing (all pure functions):
 *   ctx.loopT(t, periodSeconds)    -> u in [0, 1): position inside a seamless loop
 *   ctx.seg(u, a, b, easeFn?)      -> 0 before a, 1 after b, eased ramp between (default ease.inOutCubic)
 *   ctx.pulse(u, a, b, fade = 0.06, easeFn?) -> 0 outside [a, b], 1 inside, eased fade in/out
 *   ctx.wave(u, cycles = 1)        -> 0 → 1 → 0 cosine, seamless over u in [0, 1]
 *   ctx.clamp(x, lo = 0, hi = 1), ctx.lerp(a, b, k), ctx.smooth(x) (smoothstep)
 *   ctx.ease.{ linear, inQuad, outQuad, inOutQuad, inCubic, outCubic, inOutCubic, inSine, outSine,
 *              inOutSine, outBack, smooth }  (inputs clamped to 0..1)
 *   ctx.rand(seed)                 -> deterministic PRNG function () => [0, 1) (mulberry32), for fixed layouts
 *
 * Visual language: calm, flat, mostly orthographic, 6–12 s seamless loops, few objects, labels of 1–3
 * words, the site's colours only (TEDDY = 'teddy', ANM = 'accent', measured protein = 'measured',
 * NK / T / B / myeloid = 'nk' / 't' / 'b' / 'm'). Animations illustrate an idea; never show made-up values
 * as if measured, and show a real number only if the page already states it.
 * Size check: on a 320 px phone the stage is ~288 px wide, so in a 16 × 9 frame 1 world unit ≈ 18 px
 * (≈ 69 px on a 1100 px desktop stage). Labels stay ~11–15 px on screen at every width; shapes scale.
 * Test a scene in docs/assets/site/anim/_harness.html?name=NAME (served over http, not file://).
 */
import * as THREE from 'three';

export const VERSION = '1.0.0';
export const MAX_DPR = 2;
export const COLOR_TOKENS = Object.freeze([
  'ink', 'muted', 'faint', 'line', 'card', 'paper', 'accent', 'teddy', 'measured',
  'nk', 't', 'b', 'm', 'good', 'bad', 'warn', 'soft', 'grid', 'train',
]);
/* Light-theme values of index.html, used only if a token cannot be read at all. */
const FALLBACK = {
  ink: '#16181a', muted: '#565b62', faint: '#868b92', line: '#dcd9d0', card: '#ffffff', paper: '#f7f6f2',
  accent: '#008c7e', teddy: '#c8641f', measured: '#16181a', nk: '#c2458f', t: '#2a78d6', b: '#4a3aa7',
  m: '#a88600', good: '#1f6b3a', bad: '#a13a3a', warn: '#8a4b12', soft: '#f0efe9', grid: '#ece9e1',
  train: '#6554c9',
};

/* ───────────────────────── colour reading ───────────────────────── */
let pc = null;
function cssToHex(str) {
  if (!str) return null;
  if (!pc) {
    const c = document.createElement('canvas');
    c.width = c.height = 1;
    pc = c.getContext('2d', { willReadFrequently: true });
    if (!pc) return null;
  }
  pc.fillStyle = '#010203';
  pc.fillStyle = str;
  if (pc.fillStyle === '#010203' && !/^#010203$/i.test(str.trim())) return null;
  pc.clearRect(0, 0, 1, 1);
  pc.fillRect(0, 0, 1, 1);
  const d = pc.getImageData(0, 0, 1, 1).data;
  return '#' + [d[0], d[1], d[2]].map((v) => v.toString(16).padStart(2, '0')).join('');
}

/** Resolve every --anim-* token as seen from `el` into '#rrggbb' strings. */
export function readColors(el) {
  const probe = document.createElement('i');
  probe.setAttribute('aria-hidden', 'true');
  probe.style.cssText = 'position:absolute;width:0;height:0;overflow:hidden;visibility:hidden;pointer-events:none';
  el.appendChild(probe);
  const out = {};
  const cs = getComputedStyle(probe);
  for (const k of COLOR_TOKENS) {
    probe.style.color = `var(--anim-${k}, ${FALLBACK[k]})`;
    out[k] = cssToHex(cs.color) || FALLBACK[k];
  }
  probe.remove();
  return out;
}

function luminance(hex) {
  const n = parseInt(hex.slice(1), 16);
  const f = (v) => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
  return 0.2126 * f((n >> 16) & 255) + 0.7152 * f((n >> 8) & 255) + 0.0722 * f(n & 255);
}

/* ───────────────────────── easing and timing ───────────────────────── */
const c01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);
export const ease = Object.freeze({
  linear: (t) => c01(t),
  inQuad: (t) => { t = c01(t); return t * t; },
  outQuad: (t) => { t = c01(t); return t * (2 - t); },
  inOutQuad: (t) => { t = c01(t); return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2; },
  inCubic: (t) => { t = c01(t); return t * t * t; },
  outCubic: (t) => { t = c01(t); return 1 - Math.pow(1 - t, 3); },
  inOutCubic: (t) => { t = c01(t); return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2; },
  inSine: (t) => { t = c01(t); return 1 - Math.cos((t * Math.PI) / 2); },
  outSine: (t) => { t = c01(t); return Math.sin((t * Math.PI) / 2); },
  inOutSine: (t) => { t = c01(t); return -(Math.cos(Math.PI * t) - 1) / 2; },
  outBack: (t) => { t = c01(t); const a = 1.70158, b = a + 1; return 1 + b * Math.pow(t - 1, 3) + a * Math.pow(t - 1, 2); },
  smooth: (t) => { t = c01(t); return t * t * (3 - 2 * t); },
});
export const clamp = (x, lo = 0, hi = 1) => (x < lo ? lo : x > hi ? hi : x);
export const lerp = (a, b, k) => a + (b - a) * k;
export const smooth = ease.smooth;
export function loopT(t, period) {
  const p = period > 0 ? period : 1;
  return (((t % p) + p) % p) / p;
}
export function seg(u, a, b, fn = ease.inOutCubic) {
  if (b <= a) return u >= b ? 1 : 0;
  return fn((u - a) / (b - a));
}
export function pulse(u, a, b, fade = 0.06, fn = ease.inOutSine) {
  if (u <= a || u >= b) return 0;
  const f = Math.min(fade, (b - a) / 2);
  return Math.min(seg(u, a, a + f, fn), 1 - seg(u, b - f, b, fn));
}
export const wave = (u, cycles = 1) => 0.5 - 0.5 * Math.cos(2 * Math.PI * u * cycles);
export function rand(seed = 1) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* ───────────────────────── geometry utils ───────────────────────── */
export function toV3(p, out = new THREE.Vector3()) {
  if (p == null) return out.set(0, 0, 0);
  if (Array.isArray(p)) return out.set(+p[0] || 0, +p[1] || 0, +p[2] || 0);
  return out.set(+p.x || 0, +p.y || 0, +p.z || 0);
}
function polyAlong(pts, cum, u, out, tan) {
  const n = pts.length;
  out = out || new THREE.Vector3();
  if (n === 0) return out.set(0, 0, 0);
  if (n === 1) { if (tan) tan.set(1, 0, 0); return out.copy(pts[0]); }
  const total = cum[n - 1];
  const d = clamp(u) * total;
  let i = 1;
  while (i < n - 1 && cum[i] < d) i++;
  const segLen = cum[i] - cum[i - 1];
  const k = segLen > 0 ? (d - cum[i - 1]) / segLen : 0;
  out.lerpVectors(pts[i - 1], pts[i], clamp(k));
  if (tan) {
    tan.subVectors(pts[i], pts[i - 1]);
    if (tan.lengthSq() < 1e-12) tan.set(1, 0, 0); else tan.normalize();
  }
  return out;
}
function cumulative(pts) {
  const cum = new Float32Array(pts.length);
  for (let i = 1; i < pts.length; i++) cum[i] = cum[i - 1] + pts[i].distanceTo(pts[i - 1]);
  return cum;
}
function bentPoints(a, b, bend) {
  if (!bend) return [a.clone(), b.clone()];
  const mid = new THREE.Vector3().addVectors(a, b).multiplyScalar(0.5);
  const d = new THREE.Vector3().subVectors(b, a);
  const nrm = new THREE.Vector3(-d.y, d.x, 0).multiplyScalar(bend);
  const c = mid.add(nrm);
  const out = [];
  const N = 28;
  for (let i = 0; i <= N; i++) {
    const s = i / N, r = 1 - s;
    out.push(new THREE.Vector3(
      r * r * a.x + 2 * r * s * c.x + s * s * b.x,
      r * r * a.y + 2 * r * s * c.y + s * s * b.y,
      r * r * a.z + 2 * r * s * c.z + s * s * b.z,
    ));
  }
  return out;
}

/* ───────────────────────── screen-space line shader ───────────────────────── */
const LINE_VS = /* glsl */ `
uniform vec2 uRes;
uniform float uWidth;
uniform float uDpr;
attribute vec3 aPrev;
attribute vec3 aNext;
attribute float aSide;
attribute float aDist;
varying float vDist;
void main() {
  mat4 mvp = projectionMatrix * modelViewMatrix;
  vec4 cur = mvp * vec4(position, 1.0);
  vec4 prv = mvp * vec4(aPrev, 1.0);
  vec4 nxt = mvp * vec4(aNext, 1.0);
  vec2 asp = vec2(uRes.x / uRes.y, 1.0);
  vec2 c = cur.xy / cur.w * asp;
  vec2 p = prv.xy / prv.w * asp;
  vec2 n = nxt.xy / nxt.w * asp;
  vec2 d1 = c - p;
  vec2 d2 = n - c;
  float l1 = length(d1);
  float l2 = length(d2);
  vec2 dir = vec2(1.0, 0.0);
  float miter = 1.0;
  if (l1 > 1e-7 && l2 > 1e-7) {
    vec2 t1 = d1 / l1;
    vec2 t2 = d2 / l2;
    vec2 s = t1 + t2;
    float ls = length(s);
    dir = ls < 1e-5 ? t1 : s / ls;
    miter = 1.0 / max(dot(vec2(-dir.y, dir.x), vec2(-t1.y, t1.x)), 0.5);
  } else if (l2 > 1e-7) {
    dir = d2 / l2;
  } else if (l1 > 1e-7) {
    dir = d1 / l1;
  }
  vec2 nrm = vec2(-dir.y, dir.x);
  float hw = 0.5 * uWidth * uDpr;
  vec2 off = nrm * aSide * hw * miter * (2.0 / uRes.y);
  off.x /= asp.x;
  cur.xy += off * cur.w;
  vDist = aDist;
  gl_Position = cur;
}`;
const LINE_FS = /* glsl */ `
uniform vec3 uColor;
uniform float uOpacity;
uniform float uDash;
uniform float uGap;
uniform float uDashOffset;
uniform float uStart;
uniform float uEnd;
uniform float uPpu;
varying float vDist;
void main() {
  if (vDist < uStart || vDist > uEnd) discard;
  if (uDash > 0.0) {
    float ppu = max(uPpu, 1e-4);
    float period = (uDash + uGap) / ppu;
    float m = mod(vDist - uDashOffset / ppu, period);
    if (m > uDash / ppu) discard;
  }
  gl_FragColor = vec4(uColor, uOpacity);
  #include <colorspace_fragment>
}`;

function fillLineGeometry(geo, pts, closed) {
  const src = closed && pts.length > 2 ? pts.concat([pts[0]]) : pts;
  const n = src.length;
  const pos = new Float32Array(n * 6), prv = new Float32Array(n * 6), nxt = new Float32Array(n * 6);
  const side = new Float32Array(n * 2), dist = new Float32Array(n * 2);
  let d = 0;
  for (let i = 0; i < n; i++) {
    const p = src[i];
    let a = src[Math.max(i - 1, 0)], b = src[Math.min(i + 1, n - 1)];
    if (closed && pts.length > 2) {
      if (i === 0) a = src[n - 2];
      if (i === n - 1) b = src[1];
    }
    if (i > 0) d += p.distanceTo(src[i - 1]);
    for (let s = 0; s < 2; s++) {
      const k = i * 2 + s;
      pos[k * 3] = p.x; pos[k * 3 + 1] = p.y; pos[k * 3 + 2] = p.z;
      prv[k * 3] = a.x; prv[k * 3 + 1] = a.y; prv[k * 3 + 2] = a.z;
      nxt[k * 3] = b.x; nxt[k * 3 + 1] = b.y; nxt[k * 3 + 2] = b.z;
      side[k] = s ? 1 : -1;
      dist[k] = d;
    }
  }
  const same = geo.getAttribute('position') && geo.getAttribute('position').count === n * 2;
  if (same) {
    geo.getAttribute('position').array.set(pos);
    geo.getAttribute('aPrev').array.set(prv);
    geo.getAttribute('aNext').array.set(nxt);
    geo.getAttribute('aDist').array.set(dist);
    for (const k of ['position', 'aPrev', 'aNext', 'aDist']) geo.getAttribute(k).needsUpdate = true;
  } else {
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('aPrev', new THREE.BufferAttribute(prv, 3));
    geo.setAttribute('aNext', new THREE.BufferAttribute(nxt, 3));
    geo.setAttribute('aSide', new THREE.BufferAttribute(side, 1));
    geo.setAttribute('aDist', new THREE.BufferAttribute(dist, 1));
    const idx = [];
    for (let i = 0; i < n - 1; i++) { const a = i * 2; idx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2); }
    geo.setIndex(idx);
  }
  geo.computeBoundingSphere();
  return { total: d, pts: src, cum: cumulative(src) };
}

const ANCHORS = {
  center: [0.5, 0.5], left: [0, 0.5], right: [1, 0.5], top: [0.5, 1], bottom: [0.5, 0],
  'top-left': [0, 1], 'top-right': [1, 1], 'bottom-left': [0, 0], 'bottom-right': [1, 0],
};

/* ───────────────────────── one stage ───────────────────────── */
const _v = new THREE.Vector3(), _s = new THREE.Vector3(), _tan = new THREE.Vector3();

class AnimStage {
  constructor(o) {
    this.figure = o.figure;
    this.el = o.stage;
    this.name = o.name;
    this.onFail = o.onFail || (() => {});
    this.t = typeof o.t === 'number' ? o.t : null;
    this.last = null;
    this.running = false;
    this.reduced = !!o.reducedMotion;
    this.dead = false;
    this.screen = new Set();
    this.bound = new Map();
    this.owned = new Set();
    this.labels = new Set();
    this.renderQueued = false;
    this.width = Math.max(1, this.el.clientWidth);
    this.height = Math.max(1, this.el.clientHeight);
    this.dpr = Math.min(window.devicePixelRatio || 1, MAX_DPR);
    this.hexes = readColors(this.figure);
    this.colors = {};
    for (const k of COLOR_TOKENS) this.colors[k] = new THREE.Color(this.hexes[k]);
    const fcs = getComputedStyle(this.figure);
    this.font = fcs.fontFamily || 'system-ui, sans-serif';
    this.mono = (fcs.getPropertyValue('--mono') || '').trim() || 'ui-monospace, Menlo, monospace';

    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'low-power' });
    } catch (e) {
      this.dead = true;
      this.onFail(e);
      return;
    }
    this.renderer.setPixelRatio(this.dpr);
    this.renderer.setSize(this.width, this.height, false);
    this.renderer.setClearColor(this.colors.card, 1);
    this.canvas = this.renderer.domElement;
    this.canvas.setAttribute('aria-hidden', 'true');
    this.canvas.className = 'anim-canvas';
    this.el.appendChild(this.canvas);
    this.canvas.addEventListener('webglcontextlost', this._onLost = (e) => {
      if (this.dead) return;
      e.preventDefault();
      this.fail(new Error('WebGL context lost'), true);
    });

    this.ctx = buildCtx(this);
    let api;
    try {
      api = o.factory(this.ctx);
    } catch (e) {
      this.fail(e);
      return;
    }
    if (!api || !api.scene || !api.camera) {
      this.fail(new Error('create(ctx) must return { scene, camera, update }'));
      return;
    }
    this.api = api;
    this.scene = api.scene;
    this.camera = api.camera;
    this.still = typeof api.still === 'number' ? api.still : (api.period > 0 ? api.period * 0.6 : 0);
    if (this.t == null) this.t = o.startStill ? this.still : 0;

    this.ro = new ResizeObserver(() => this.resize());
    this.ro.observe(this.el);
    this._fonts = () => { this.labels.forEach((it) => { it.dirty = true; }); this.requestRender(); };
    if (document.fonts) {
      document.fonts.ready.then(() => { if (!this.dead) this._fonts(); });
      document.fonts.addEventListener && document.fonts.addEventListener('loadingdone', this._fonts);
    }
    this.fitCamera();
    if (this.api.resize) this.safe(() => this.api.resize(this.width, this.height));
    this.draw(0);
  }

  get textScale() { return clamp(this.width / 640, 0.82, 1.12); }
  get theme() { return luminance(this.hexes.card) < 0.3 ? 'dark' : 'light'; }

  safe(fn) {
    try { fn(); return true; } catch (e) { this.fail(e); return false; }
  }

  fail(err, quiet) {
    if (this.dead) return;
    if (!quiet) console.warn(`[anim] ${this.name}: ${err && err.message ? err.message : err}`);
    this.dispose();
    this.onFail(err);
  }

  ppuAt(pos) {
    const cam = this.camera;
    if (!cam) return 1;
    if (cam.isOrthographicCamera) return (this.height * cam.zoom) / Math.max(1e-6, cam.top - cam.bottom);
    if (cam.isPerspectiveCamera) {
      _v.copy(pos || _s.set(0, 0, 0)).applyMatrix4(cam.matrixWorldInverse);
      const dist = Math.max(1e-3, -_v.z);
      return (this.height * cam.zoom) / (2 * dist * Math.tan(THREE.MathUtils.degToRad(cam.fov) / 2));
    }
    return 1;
  }

  fitCamera() {
    const cam = this.camera, a = this.width / this.height;
    if (!cam) return;
    const fit = cam.userData && cam.userData.animFit;
    if (cam.isOrthographicCamera && fit) {
      let hw, hh;
      const da = fit.width / fit.height;
      const wider = a > da;
      if ((fit.fit === 'cover') !== wider) { hw = fit.width / 2; hh = hw / a; } else { hh = fit.height / 2; hw = hh * a; }
      cam.left = -hw; cam.right = hw; cam.top = hh; cam.bottom = -hh;
      cam.updateProjectionMatrix();
    } else if (cam.isPerspectiveCamera) {
      cam.aspect = a;
      cam.updateProjectionMatrix();
    } else if (cam.isOrthographicCamera) {
      const hh = (cam.top - cam.bottom) / 2, cx = (cam.left + cam.right) / 2;
      cam.left = cx - hh * a; cam.right = cx + hh * a;
      cam.updateProjectionMatrix();
    }
  }

  resize() {
    if (this.dead) return;
    const w = this.el.clientWidth, h = this.el.clientHeight;
    const dpr = Math.min(window.devicePixelRatio || 1, MAX_DPR);
    if (!w || !h) return;
    if (w === this.width && h === this.height && dpr === this.dpr) return;
    this.width = w; this.height = h;
    if (dpr !== this.dpr) { this.dpr = dpr; this.renderer.setPixelRatio(dpr); }
    this.renderer.setSize(w, h, false);
    this.fitCamera();
    if (this.api.resize && !this.safe(() => this.api.resize(w, h))) return;
    if (!this.running) this.draw(0);
  }

  refreshTheme() {
    if (this.dead) return;
    const next = readColors(this.figure);
    let changed = false;
    for (const k of COLOR_TOKENS) if (next[k] !== this.hexes[k]) { changed = true; break; }
    if (!changed) return;
    this.hexes = next;
    for (const k of COLOR_TOKENS) this.colors[k].set(next[k]);
    this.bound.forEach((b, target) => applyBinding(this, target, b));
    this.labels.forEach((it) => { it.dirty = true; });
    this.renderer.setClearColor(this.colors.card, 1);
    if (this.api.onTheme && !this.safe(() => this.api.onTheme())) return;
    if (!this.running) this.draw(0);
  }

  setReduced(on) {
    this.reduced = !!on;
  }

  showStill() {
    this.t = this.still;
    this.draw(0);
  }

  setRunning(on) {
    if (on && !this.running) this.last = null;
    this.running = !!on;
  }

  tick(now) {
    if (this.dead) return;
    const dt = this.last == null ? 0 : Math.min(0.1, Math.max(0, (now - this.last) / 1000));
    this.last = now;
    this.t += dt;
    this.draw(dt);
  }

  updateScreen() {
    this.scene.updateMatrixWorld();
    this.camera.updateMatrixWorld();
    this.screen.forEach((it) => it.update());
  }

  draw(dt) {
    if (this.dead) return;
    if (this.api.update && !this.safe(() => this.api.update(this.t, dt))) return;
    if (this.dead) return;
    this.updateScreen();
    this.renderer.render(this.scene, this.camera);
  }

  requestRender() {
    if (this.dead || this.running || this.renderQueued) return;
    this.renderQueued = true;
    requestAnimationFrame(() => { this.renderQueued = false; if (!this.running) this.draw(0); });
  }

  dispose() {
    if (this.dead) return;
    this.dead = true;
    this.running = false;
    if (this.ro) this.ro.disconnect();
    if (document.fonts && document.fonts.removeEventListener && this._fonts) document.fonts.removeEventListener('loadingdone', this._fonts);
    if (this.api && this.api.dispose) { try { this.api.dispose(); } catch (e) { /* scene cleanup is best effort */ } }
    const seen = this.owned;
    if (this.scene) {
      this.scene.traverse((o) => {
        if (o.geometry) seen.add(o.geometry);
        const mats = Array.isArray(o.material) ? o.material : o.material ? [o.material] : [];
        mats.forEach((m) => {
          seen.add(m);
          for (const key in m) { const v = m[key]; if (v && v.isTexture) seen.add(v); }
          if (m.uniforms) for (const key in m.uniforms) { const v = m.uniforms[key].value; if (v && v.isTexture) seen.add(v); }
        });
      });
    }
    seen.forEach((d) => { try { d.dispose(); } catch (e) { /* ignore */ } });
    this.owned.clear(); this.screen.clear(); this.bound.clear(); this.labels.clear();
    try { this.renderer.dispose(); } catch (e) { /* ignore */ }
    try { this.renderer.forceContextLoss(); } catch (e) { /* ignore */ }
    if (this.canvas && this.canvas.parentNode) this.canvas.parentNode.removeChild(this.canvas);
    this.scene = null; this.camera = null;
  }
}

/* ───────────────────────── colour binding ───────────────────────── */
function colorTarget(target, prop) {
  if (!target) return null;
  if (target.isColor) return target;
  if (prop && target[prop] && target[prop].isColor) return target[prop];
  if (target.value && target.value.isColor) return target.value;
  return null;
}
function applyBinding(S, target, b) {
  const dst = colorTarget(target, b.prop);
  const src = S.colors[b.token];
  if (dst && src) dst.copy(src);
}
function isToken(S, c) {
  return typeof c === 'string' && Object.prototype.hasOwnProperty.call(S.colors, c);
}
function setColorOf(S, target, c, prop = 'color') {
  if (c == null) return;
  if (isToken(S, c)) {
    const b = { token: c, prop };
    S.bound.set(target, b);
    applyBinding(S, target, b);
  } else {
    S.bound.delete(target);
    const dst = colorTarget(target, prop);
    if (dst) dst.set(c);
  }
}
function cssColorOf(S, c) {
  if (c == null) return null;
  if (isToken(S, c)) return S.hexes[c];
  if (c.isColor) return '#' + c.getHexString();
  if (typeof c === 'number') return '#' + new THREE.Color(c).getHexString();
  return String(c);
}

/* ───────────────────────── ctx + helpers ───────────────────────── */
function buildCtx(S) {
  const track = (o) => { if (o) S.owned.add(o); return o; };
  let circleGeo = null;
  const ringGeos = new Map();
  const circle = () => circleGeo || (circleGeo = track(new THREE.CircleGeometry(1, 48)));
  const ring = (k) => {
    const key = Math.round(clamp(k, 0.02, 0.98) * 100);
    if (!ringGeos.has(key)) ringGeos.set(key, track(new THREE.RingGeometry(1 - key / 100, 1, 48)));
    return ringGeos.get(key);
  };
  const flatMat = (c, opacity) => {
    const m = track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, opacity, side: THREE.DoubleSide }));
    setColorOf(S, m, c);
    return m;
  };

  function label(text, opts = {}) {
    const o = Object.assign({ size: 13, color: 'ink', anchor: 'center', weight: 500, mono: false, bg: null,
      bgOpacity: 0.92, pad: 4, align: null, opacity: 1, order: 10 }, opts);
    const mat = track(new THREE.SpriteMaterial({ transparent: true, depthTest: false, depthWrite: false, opacity: o.opacity }));
    const spr = new THREE.Sprite(mat);
    spr.renderOrder = o.order;
    const it = { spr, text: String(text), o, dirty: true, px: 0, dpr: 0, cssW: 1, cssH: 1, tex: null, canvas: document.createElement('canvas') };
    const setAnchor = (a) => { const v = ANCHORS[a] || ANCHORS.center; spr.center.set(v[0], v[1]); o.anchor = a; };
    setAnchor(o.anchor);
    it.raster = () => {
      const px = o.size * S.textScale;
      const dpr = S.dpr;
      const fpx = px * dpr;
      const c = it.canvas, g = c.getContext('2d');
      const family = o.mono ? S.mono : S.font;
      const font = `${o.weight} ${fpx.toFixed(2)}px ${family}`;
      g.font = font;
      const lines = it.text.split('\n');
      const lh = fpx * 1.22;
      let maxW = 0;
      for (const l of lines) maxW = Math.max(maxW, g.measureText(l).width);
      const padX = (o.bg ? o.pad * 1.6 : 1.5) * dpr, padY = (o.bg ? o.pad : 1) * dpr;
      const W = Math.max(2, Math.ceil(maxW + padX * 2)), H = Math.max(2, Math.ceil(lh * lines.length + padY * 2));
      c.width = W; c.height = H;
      g.clearRect(0, 0, W, H);
      if (o.bg) {
        g.globalAlpha = o.bgOpacity;
        g.fillStyle = cssColorOf(S, o.bg);
        const r = Math.min(H / 2, 6 * dpr);
        g.beginPath();
        if (g.roundRect) g.roundRect(0, 0, W, H, r); else g.rect(0, 0, W, H);
        g.fill();
        g.globalAlpha = 1;
      }
      g.font = font;
      g.fillStyle = cssColorOf(S, o.color);
      g.textBaseline = 'middle';
      const align = o.align || (/left/.test(o.anchor) ? 'left' : /right/.test(o.anchor) ? 'right' : 'center');
      g.textAlign = align;
      const x = align === 'left' ? padX : align === 'right' ? W - padX : W / 2;
      lines.forEach((l, i) => g.fillText(l, x, padY + lh * (i + 0.5) + fpx * 0.04));
      if (it.tex) it.tex.dispose();
      const tex = new THREE.CanvasTexture(c);
      tex.colorSpace = THREE.SRGBColorSpace;
      tex.minFilter = THREE.LinearFilter;
      tex.magFilter = THREE.LinearFilter;
      tex.generateMipmaps = false;
      it.tex = track(tex);
      mat.map = tex;
      mat.needsUpdate = true;
      it.cssW = W / dpr; it.cssH = H / dpr; it.px = px; it.dpr = dpr; it.dirty = false;
    };
    it.update = () => {
      if (it.dirty || Math.abs(o.size * S.textScale - it.px) > 0.2 || it.dpr !== S.dpr) it.raster();
      spr.getWorldPosition(_v);
      const ppu = S.ppuAt(_v);
      if (spr.parent) spr.parent.getWorldScale(_s); else _s.set(1, 1, 1);
      spr.scale.set(it.cssW / ppu / (_s.x || 1), it.cssH / ppu / (_s.y || 1), 1);
    };
    it.raster();
    S.screen.add(it);
    S.labels.add(it);
    spr.setText = (t) => { t = String(t); if (t !== it.text) { it.text = t; it.dirty = true; S.requestRender(); } return spr; };
    spr.setColor = (c) => { o.color = c; it.dirty = true; S.requestRender(); return spr; };
    spr.setSize = (px) => { if (px !== o.size) { o.size = px; it.dirty = true; S.requestRender(); } return spr; };
    spr.setOpacity = (a) => { mat.opacity = a; return spr; };
    spr.setAnchor = (a) => { setAnchor(a); it.dirty = true; return spr; };
    return spr;
  }

  function line(points, opts = {}) {
    const o = Object.assign({ color: 'muted', width: 2, dashed: false, dash: 6, gap: 5, opacity: 1, closed: false, order: 0 }, opts);
    if (Array.isArray(o.dashed)) { o.dash = o.dashed[0]; o.gap = o.dashed[1]; o.dashed = true; }
    const uniforms = {
      uColor: { value: new THREE.Color() }, uOpacity: { value: o.opacity }, uWidth: { value: o.width },
      uDpr: { value: S.dpr }, uRes: { value: new THREE.Vector2(1, 1) },
      uDash: { value: o.dashed ? o.dash : 0 }, uGap: { value: o.gap }, uDashOffset: { value: 0 },
      uStart: { value: 0 }, uEnd: { value: 1e9 }, uPpu: { value: 1 },
    };
    const mat = track(new THREE.ShaderMaterial({ uniforms, vertexShader: LINE_VS, fragmentShader: LINE_FS, transparent: true, depthWrite: false, side: THREE.DoubleSide }));
    setColorOf(S, uniforms.uColor, o.color);
    const geo = track(new THREE.BufferGeometry());
    const mesh = new THREE.Mesh(geo, mat);
    mesh.frustumCulled = false;
    mesh.renderOrder = o.order;
    let info = null, p0 = 0, p1 = 1;
    const norm = (pts) => pts.map((p) => toV3(p));
    const applyProgress = () => {
      uniforms.uStart.value = info.total * p0 - 1e-5;
      uniforms.uEnd.value = p1 >= 1 ? 1e9 : info.total * p1;
      mesh.visible = p1 > p0;
    };
    mesh.setPoints = (pts) => {
      const v = norm(pts);
      if (v.length < 2) v.push(v[0] ? v[0].clone() : new THREE.Vector3());
      info = fillLineGeometry(geo, v, o.closed);
      mesh.userData.animPts = info.pts; mesh.userData.animCum = info.cum;
      applyProgress();
      return mesh;
    };
    mesh.setProgress = (p, from = 0) => { p1 = clamp(p); p0 = clamp(from); applyProgress(); return mesh; };
    mesh.setDashOffset = (px) => { uniforms.uDashOffset.value = px; return mesh; };
    mesh.setDashPhase = (u, cycles = 1) => { uniforms.uDashOffset.value = u * cycles * (uniforms.uDash.value + uniforms.uGap.value); return mesh; };
    mesh.setColor = (c) => { setColorOf(S, uniforms.uColor, c); return mesh; };
    mesh.setOpacity = (a) => { uniforms.uOpacity.value = a; return mesh; };
    mesh.setWidth = (px) => { uniforms.uWidth.value = px; return mesh; };
    mesh.setDashed = (d) => {
      if (Array.isArray(d)) { uniforms.uDash.value = d[0]; uniforms.uGap.value = d[1]; } else uniforms.uDash.value = d ? o.dash : 0;
      return mesh;
    };
    Object.defineProperty(mesh, 'length', { get: () => info.total });
    mesh._lineEnd = (end) => { uniforms.uEnd.value = end; };
    mesh.setPoints(points);
    S.screen.add({
      update() {
        uniforms.uDpr.value = S.dpr;
        uniforms.uRes.value.set(S.width * S.dpr, S.height * S.dpr);
        mesh.getWorldPosition(_v);
        uniforms.uPpu.value = S.ppuAt(_v);
      },
    });
    return mesh;
  }

  let headGeo = null;
  function arrow(from, to, opts = {}) {
    const o = Object.assign({ color: 'muted', width: 2, head: 9, bend: 0, dashed: false, opacity: 1, order: 0 }, opts);
    const g = new THREE.Group();
    let a = toV3(from), b = toV3(to), bend = o.bend, p = 1;
    const ln = line(bentPoints(a, b, bend), { color: o.color, width: o.width, dashed: o.dashed, opacity: o.opacity, order: o.order });
    if (!headGeo) {
      headGeo = track(new THREE.BufferGeometry());
      headGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array([0, 0, 0, -1, 0.46, 0, -1, -0.46, 0, -0.78, 0, 0]), 3));
      headGeo.setIndex([0, 1, 3, 0, 3, 2]);
    }
    const hm = flatMat(o.color, o.opacity);
    const head = new THREE.Mesh(headGeo, hm);
    head.renderOrder = o.order;
    g.add(ln, head);
    g.line = ln; g.head = head;
    g.set = (f, t2, bnd) => { a = toV3(f); b = toV3(t2); if (bnd !== undefined) bend = bnd; ln.setPoints(bentPoints(a, b, bend)); return g; };
    g.setProgress = (v) => { p = clamp(v); return g; };
    g.setColor = (c) => { ln.setColor(c); setColorOf(S, hm, c); return g; };
    g.setOpacity = (v) => { ln.setOpacity(v); hm.opacity = v; return g; };
    g.setDashOffset = (px) => { ln.setDashOffset(px); return g; };
    g.setDashPhase = (u, cycles = 1) => { ln.setDashPhase(u, cycles); return g; };
    S.screen.add({
      update() {
        const pts = ln.userData.animPts, cum = ln.userData.animCum;
        polyAlong(pts, cum, p, head.position, _tan);
        head.rotation.set(0, 0, Math.atan2(_tan.y, _tan.x));
        head.getWorldPosition(_v);
        const ppu = S.ppuAt(_v);
        g.getWorldScale(_s);
        const L = o.head / ppu / (_s.x || 1);
        head.scale.set(L, L, 1);
        head.visible = p > 0.001 && ln.length > 0;
        ln.visible = p > 0.001;
        ln._lineEnd(Math.max(0, ln.length * p - L * 0.7));
      },
    });
    return g;
  }

  function dot(pos, opts = {}) {
    const o = Object.assign({ r: 0.12, px: null, color: 'ink', opacity: 1, hollow: false, ring: 0.3, order: 0 }, opts);
    const m = flatMat(o.color, o.opacity);
    const mesh = new THREE.Mesh(o.hollow ? ring(o.ring) : circle(), m);
    mesh.renderOrder = o.order;
    toV3(pos, mesh.position);
    let r = o.r, px = o.px;
    mesh.scale.set(r, r, 1);
    mesh.setColor = (c) => { setColorOf(S, m, c); return mesh; };
    mesh.setOpacity = (a) => { m.opacity = a; return mesh; };
    mesh.setRadius = (v) => { r = v; px = null; mesh.scale.set(v, v, 1); return mesh; };
    mesh.setPx = (v) => { px = v; return mesh; };
    S.screen.add({
      update() {
        if (px == null) return;
        mesh.getWorldPosition(_v);
        const s = px / S.ppuAt(_v);
        mesh.scale.set(s, s, 1);
      },
    });
    return mesh;
  }

  function roundRectPoints(w, h, r) {
    r = Math.max(0, Math.min(r, w / 2, h / 2));
    const s = new THREE.Shape();
    const x = -w / 2, y = -h / 2;
    s.moveTo(x + r, y);
    s.lineTo(x + w - r, y);
    if (r) s.absarc(x + w - r, y + r, r, -Math.PI / 2, 0, false);
    s.lineTo(x + w, y + h - r);
    if (r) s.absarc(x + w - r, y + h - r, r, 0, Math.PI / 2, false);
    s.lineTo(x + r, y + h);
    if (r) s.absarc(x + r, y + h - r, r, Math.PI / 2, Math.PI, false);
    s.lineTo(x, y + r);
    if (r) s.absarc(x + r, y + r, r, Math.PI, Math.PI * 1.5, false);
    return s;
  }

  function box(w, h, opts = {}) {
    const o = Object.assign({ color: 'soft', radius: 0.18, opacity: 1, stroke: null, strokeWidth: 1.5, dashed: false, order: 0 }, opts);
    const m = flatMat(o.color == null ? 'card' : o.color, o.opacity);
    if (o.color == null) m.visible = false;
    let shape = roundRectPoints(w, h, o.radius);
    const mesh = new THREE.Mesh(track(new THREE.ShapeGeometry(shape, 10)), m);
    mesh.renderOrder = o.order;
    mesh.outline = null;
    const outlinePts = () => {
      const pts = shape.getPoints(10);
      const clean = [];
      for (const p of pts) { if (!clean.length || clean[clean.length - 1].distanceTo(p) > 1e-6) clean.push(p); }
      if (clean.length > 2 && clean[0].distanceTo(clean[clean.length - 1]) < 1e-6) clean.pop();
      return clean;
    };
    const addOutline = (c) => {
      mesh.outline = line(outlinePts(), { color: c, width: o.strokeWidth, dashed: o.dashed, closed: true, opacity: o.opacity, order: o.order });
      mesh.add(mesh.outline);
    };
    if (o.stroke != null) addOutline(o.stroke);
    mesh.setColor = (c) => { if (c == null) { m.visible = false; } else { m.visible = true; setColorOf(S, m, c); } return mesh; };
    mesh.setStroke = (c) => { if (!mesh.outline) addOutline(c); else mesh.outline.setColor(c); return mesh; };
    mesh.setOpacity = (a) => { m.opacity = a; if (mesh.outline) mesh.outline.setOpacity(a); return mesh; };
    mesh.setSize = (nw, nh, nr = o.radius) => {
      o.radius = nr;
      shape = roundRectPoints(nw, nh, nr);
      const old = mesh.geometry;
      mesh.geometry = track(new THREE.ShapeGeometry(shape, 10));
      old.dispose(); S.owned.delete(old);
      if (mesh.outline) mesh.outline.setPoints(outlinePts());
      return mesh;
    };
    return mesh;
  }

  function orthoCamera(opts = {}) {
    const o = Object.assign({ width: 16, height: 9, fit: 'contain', center: [0, 0], zoom: 1 }, opts);
    const cam = new THREE.OrthographicCamera(-o.width / 2, o.width / 2, o.height / 2, -o.height / 2, 0.1, 1000);
    const c = toV3(o.center);
    cam.position.set(c.x, c.y, 100);
    cam.lookAt(c.x, c.y, 0);
    cam.zoom = o.zoom;
    cam.userData.animFit = { width: o.width, height: o.height, fit: o.fit };
    return cam;
  }

  function along(src, u, tan) {
    let pts, cum;
    if (src && src.isObject3D && src.userData.animPts) { pts = src.userData.animPts; cum = src.userData.animCum; }
    else { pts = src.map((p) => toV3(p)); cum = cumulative(pts); }
    return polyAlong(pts, cum, u, new THREE.Vector3(), tan);
  }

  const ctx = {
    THREE, VERSION,
    name: S.name, figure: S.figure, stage: S.el, canvas: S.canvas, renderer: S.renderer,
    get width() { return S.width; },
    get height() { return S.height; },
    get aspect() { return S.width / S.height; },
    get dpr() { return S.dpr; },
    get theme() { return S.theme; },
    get reducedMotion() { return S.reduced && !S.running; },
    get textScale() { return S.textScale; },
    colors: S.colors,
    color: (tok) => (S.colors[tok] ? S.colors[tok].clone() : new THREE.Color(tok)),
    hex: (tok) => S.hexes[tok] || cssColorOf(S, tok),
    bind: (target, token, prop = 'color') => { setColorOf(S, target, token, prop); return target; },
    unbind: (target) => { S.bound.delete(target); return target; },
    label, line, arrow, dot, box,
    group: () => new THREE.Group(),
    orthoCamera,
    v3: (x = 0, y = 0, z = 0) => new THREE.Vector3(x, y, z),
    toV3: (p) => toV3(p),
    along,
    ppu: (pos) => S.ppuAt(pos ? toV3(pos) : undefined),
    requestRender: () => S.requestRender(),
    track,
    ease, clamp, lerp, smooth, loopT, seg, pulse, wave, rand,
  };
  return ctx;
}

/**
 * Build one live stage. Returns the stage, or null when the scene failed (onFail is then called).
 * Used by loader.js only.
 */
export function createStage(opts) {
  const S = new AnimStage(opts);
  return S.dead ? null : S;
}
