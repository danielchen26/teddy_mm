/*
 * teddy_mm animation loader. Include once per page, after the import map for "three":
 *   <script type="importmap">{"imports":{"three":"./assets/vendor/three/three.module.min.js"}}</script>
 *   <link rel="stylesheet" href="assets/site/anim/anim.css">
 *   <script type="module" src="assets/site/anim/loader.js"></script>
 *
 * It finds every <figure class="anim" data-anim="NAME"> (also ones added later), and:
 *   - mounts scenes/NAME.js when the stage is within ~0.75 viewport of the screen, disposes it
 *     (renderer.dispose + forceContextLoss) when it moves away; at most 6 live WebGL contexts
 *     (override with <html data-anim-max-live="N">), nearest-to-centre stages win;
 *   - runs one requestAnimationFrame loop for all stages, only for stages on screen, not paused,
 *     and only while the tab is visible;
 *   - click/tap on the stage, or the corner button, toggles pause;
 *   - prefers-reduced-motion (or <html data-motion="reduce">) shows one still frame per scene;
 *     the button then plays it on request. <html data-motion="full"> forces motion on;
 *   - follows the theme: <html data-theme="light|dark"> and prefers-color-scheme changes;
 *   - without WebGL2 (or if a scene fails) the figure shows its caption only.
 * Loading core.js (and three.js) is deferred until a figure first needs it.
 * Figure attributes set here: data-state = playing | paused | static | fallback, data-live = 1 | 0;
 * events on the figure (bubbling): anim:mount, anim:dispose, anim:fallback.
 * window.teddyAnim (harness and tests): version, scan(), maxLive, setMaxLive(n), liveCount(), stats,
 * records(), toggle(el), pauseAll(), playAll(), refreshTheme(), seek(el, t).
 */
const VERSION = '1.0.0';
const BASE = new URL('./', import.meta.url);
const SCENE_DIR = new URL('scenes/', BASE);
const NAME_RE = /^[a-z0-9_][a-z0-9_-]{0,80}$/;
const root = document.documentElement;
const recs = new Map();
const stats = { mounts: 0, disposes: 0, fallbacks: 0, maxLiveSeen: 0 };
let maxLive = parseInt(root.getAttribute('data-anim-max-live') || window.TEDDY_ANIM_MAX_LIVE || '6', 10) || 6;
let overCap = false;

/* ───────────── lazy modules ───────────── */
let coreP = null;
const loadCore = () => coreP || (coreP = import('./core.js'));
const sceneP = new Map();
function loadScene(name) {
  if (!sceneP.has(name)) sceneP.set(name, import(new URL(name + '.js', SCENE_DIR).href));
  return sceneP.get(name);
}

let glOK = null;
function webglAvailable() {
  if (glOK !== null) return glOK;
  try {
    const c = document.createElement('canvas');
    const gl = c.getContext('webgl2');
    glOK = !!gl;
    if (gl) { const ext = gl.getExtension('WEBGL_lose_context'); if (ext) ext.loseContext(); }
  } catch (e) {
    glOK = false;
  }
  return glOK;
}

/* ───────────── motion preference ───────────── */
const mqReduce = window.matchMedia ? window.matchMedia('(prefers-reduced-motion: reduce)') : { matches: false };
function reduced() {
  const o = root.getAttribute('data-motion');
  if (o === 'reduce') return true;
  if (o === 'full' || o === 'no-preference') return false;
  return !!mqReduce.matches;
}
const playing = (r) => (r.intent ? r.intent === 'play' : !reduced());
const running = (r) => !!r.live && r.visible && playing(r) && !document.hidden;

/* ───────────── per-figure state ───────────── */
const ICONS =
  '<svg class="anim-i anim-i-pause" viewBox="0 0 16 16" aria-hidden="true"><rect x="4" y="3" width="2.8" height="10" rx="1"/><rect x="9.2" y="3" width="2.8" height="10" rx="1"/></svg>' +
  '<svg class="anim-i anim-i-play" viewBox="0 0 16 16" aria-hidden="true"><path d="M5 3.3v9.4c0 .5.5.8.9.5l7.3-4.7c.4-.2.4-.8 0-1L5.9 2.8c-.4-.3-.9 0-.9.5z"/></svg>' +
  '<span class="anim-btn-t"></span>';

function applyAspect(r) {
  const a = r.stage.getAttribute('data-aspect') || r.fig.getAttribute('data-aspect');
  if (!a) return;
  const m = String(a).trim().match(/^(\d+(?:\.\d+)?)\s*(?:[:/x×]\s*(\d+(?:\.\d+)?))?$/);
  if (!m) return;
  const w = +m[1], h = m[2] ? +m[2] : 1;
  if (w > 0 && h > 0) r.stage.style.aspectRatio = `${w} / ${h}`;
}

function setState(r) {
  const disp = r.failed ? 'fallback' : playing(r) ? 'playing' : r.intent === 'pause' ? 'paused' : 'static';
  r.fig.setAttribute('data-state', disp);
  r.fig.setAttribute('data-live', r.live ? '1' : '0');
  if (r.btn) {
    const t = r.btn.querySelector('.anim-btn-t');
    if (disp === 'playing') { r.btn.setAttribute('aria-label', 'Pause animation'); t.textContent = ''; }
    else if (disp === 'paused') { r.btn.setAttribute('aria-label', 'Play animation'); t.textContent = 'Paused'; }
    else { r.btn.setAttribute('aria-label', 'Play animation'); t.textContent = 'Play'; }
    r.btn.setAttribute('aria-pressed', disp === 'playing' ? 'false' : 'true');
  }
  if (r.live) {
    r.live.setReduced(!playing(r) && r.intent == null);
    r.live.setRunning(running(r));
  }
  kick();
}

function toggle(r) {
  if (!r || r.failed) return;
  r.intent = playing(r) ? 'pause' : 'play';
  setState(r);
}

function fallback(r) {
  if (r.failed) return;
  r.failed = true;
  r.gen++;
  if (r.live) { const s = r.live; r.live = null; s.dispose(); }
  r.fig.classList.add('anim--fallback');
  stats.fallbacks++;
  if (nearIO) nearIO.unobserve(r.stage);
  if (visIO) visIO.unobserve(r.stage);
  setState(r);
  r.fig.dispatchEvent(new CustomEvent('anim:fallback', { bubbles: true }));
}

function setup(fig) {
  const stage = fig.querySelector(':scope > .anim-stage') || fig.querySelector('.anim-stage');
  if (!stage || recs.has(stage)) return;
  const name = (fig.getAttribute('data-anim') || '').trim();
  const r = { fig, stage, name, near: false, visible: false, want: false, live: null, mounting: false,
    gen: 0, t: null, intent: null, failed: false, btn: null, rank: 0 };
  recs.set(stage, r);
  applyAspect(r);
  if (!NAME_RE.test(name)) { console.warn(`[anim] invalid data-anim name "${name}"`); fallback(r); return; }
  const btn = document.createElement('button');
  btn.type = 'button';
  btn.className = 'anim-btn';
  btn.innerHTML = ICONS;
  btn.addEventListener('click', () => toggle(r));
  stage.insertAdjacentElement('afterend', btn);
  r.btn = btn;
  stage.addEventListener('click', () => toggle(r));
  if (nearIO) { nearIO.observe(stage); visIO.observe(stage); } else { r.near = r.visible = true; }
  setState(r);
}

/* ───────────── mount / dispose ───────────── */
async function mount(r) {
  if (r.live || r.mounting || r.failed) return;
  r.mounting = true;
  const gen = ++r.gen;
  if (!webglAvailable()) { r.mounting = false; fallback(r); return; }
  let core, mod;
  try {
    [core, mod] = await Promise.all([loadCore(), loadScene(r.name)]);
  } catch (e) {
    r.mounting = false;
    if (gen === r.gen) { console.warn(`[anim] ${r.name}: scene did not load (${e && e.message ? e.message : e})`); fallback(r); }
    return;
  }
  if (gen !== r.gen) return;
  r.mounting = false;
  if (!r.want || !r.stage.isConnected || r.failed) return;
  if (typeof mod.default !== 'function') { console.warn(`[anim] ${r.name}: no default export create(ctx)`); fallback(r); return; }
  const live = core.createStage({
    figure: r.fig, stage: r.stage, name: r.name, factory: mod.default,
    t: r.t, startStill: !playing(r), reducedMotion: !playing(r) && r.intent == null,
    onFail: () => { r.live = null; fallback(r); },
  });
  if (!live) return;
  r.live = live;
  r.t = live.t;
  stats.mounts++;
  stats.maxLiveSeen = Math.max(stats.maxLiveSeen, liveCount());
  setState(r);
  r.fig.dispatchEvent(new CustomEvent('anim:mount', { bubbles: true }));
}

function unmount(r) {
  r.gen++;
  r.mounting = false;
  if (!r.live) return;
  const s = r.live;
  r.t = s.t;
  r.live = null;
  s.dispose();
  stats.disposes++;
  setState(r);
  r.fig.dispatchEvent(new CustomEvent('anim:dispose', { bubbles: true }));
}

function liveCount() {
  let n = 0;
  recs.forEach((r) => { if (r.live) n++; });
  return n;
}

let schedQueued = false;
function schedule() {
  if (schedQueued) return;
  schedQueued = true;
  requestAnimationFrame(runSchedule);
}
function runSchedule() {
  schedQueued = false;
  const vh = window.innerHeight || 800;
  const cands = [];
  recs.forEach((r) => {
    r.want = false;
    if (!r.near || r.failed || !r.stage.isConnected) return;
    const b = r.stage.getBoundingClientRect();
    if (!b.width || !b.height) return;
    r.rank = Math.abs((b.top + b.bottom) / 2 - vh / 2) - (r.visible ? 1e7 : 0);
    cands.push(r);
  });
  cands.sort((a, b) => a.rank - b.rank);
  for (let i = 0; i < cands.length && i < maxLive; i++) cands[i].want = true;
  overCap = cands.length > maxLive;
  recs.forEach((r) => { if (!r.want && (r.live || r.mounting)) unmount(r); });
  recs.forEach((r) => { if (r.want) mount(r); });
}

/* ───────────── observers ───────────── */
const hasIO = 'IntersectionObserver' in window;
const nearIO = hasIO ? new IntersectionObserver((es) => {
  es.forEach((e) => { const r = recs.get(e.target); if (r) r.near = e.isIntersecting; });
  schedule();
}, { rootMargin: '75% 0px 75% 0px' }) : null;
const visIO = hasIO ? new IntersectionObserver((es) => {
  es.forEach((e) => { const r = recs.get(e.target); if (r) { r.visible = e.isIntersecting; setState(r); } });
  schedule();
}, { threshold: 0 }) : null;

let scrollQueued = false;
window.addEventListener('scroll', () => {
  if (!overCap || scrollQueued) return;
  scrollQueued = true;
  setTimeout(() => { scrollQueued = false; schedule(); }, 120);
}, { passive: true });

/* ───────────── loop ───────────── */
let raf = 0;
function kick() {
  if (!raf && !document.hidden) raf = requestAnimationFrame(loop);
}
function loop(now) {
  raf = 0;
  let any = false;
  recs.forEach((r) => {
    const s = r.live;
    if (s && s.running) {
      s.tick(now);
      if (r.live) { r.t = r.live.t; any = true; }
    }
  });
  if (any) kick();
}
document.addEventListener('visibilitychange', () => { recs.forEach(setState); });

/* ───────────── theme + motion ───────────── */
let themeQueued = false;
function onTheme() {
  if (themeQueued) return;
  themeQueued = true;
  requestAnimationFrame(() => {
    themeQueued = false;
    recs.forEach((r) => { if (r.live) r.live.refreshTheme(); });
  });
}
function onMotion() {
  const red = reduced();
  recs.forEach((r) => {
    r.intent = null;
    if (r.live) {
      r.live.setReduced(red);
      if (red) r.live.showStill();
    }
    setState(r);
  });
}
new MutationObserver((muts) => {
  let theme = false, motion = false;
  for (const m of muts) { if (m.attributeName === 'data-motion') motion = true; else theme = true; }
  if (theme) onTheme();
  if (motion) onMotion();
}).observe(root, { attributes: true, attributeFilter: ['data-theme', 'class', 'data-motion'] });
if (window.matchMedia) {
  const mqDark = window.matchMedia('(prefers-color-scheme: dark)');
  if (mqDark.addEventListener) mqDark.addEventListener('change', onTheme);
  if (mqReduce.addEventListener) mqReduce.addEventListener('change', onMotion);
}

/* ───────────── DOM discovery ───────────── */
function prune() {
  recs.forEach((r, st) => {
    if (st.isConnected) return;
    unmount(r);
    if (nearIO) { nearIO.unobserve(st); visIO.unobserve(st); }
    recs.delete(st);
  });
}
function scan() {
  document.querySelectorAll('.anim[data-anim]').forEach(setup);
  prune();
  schedule();
}
const touchesAnim = (n) => n.nodeType === 1 && (n.matches('.anim, .anim-stage') || !!n.querySelector('.anim'));
function watchDom() {
  new MutationObserver((muts) => {
    for (const m of muts) {
      for (const n of m.addedNodes) if (touchesAnim(n)) { scan(); return; }
      for (const n of m.removedNodes) if (touchesAnim(n)) { scan(); return; }
    }
  }).observe(document.body, { childList: true, subtree: true });
}
function init() {
  watchDom();
  new MutationObserver(onTheme).observe(document.body, { attributes: true, attributeFilter: ['data-theme', 'class'] });
  scan();
}
if (document.body) init();
else document.addEventListener('DOMContentLoaded', init);

/* ───────────── inspection API ───────────── */
window.teddyAnim = {
  version: VERSION,
  scan,
  get maxLive() { return maxLive; },
  setMaxLive(n) { maxLive = Math.max(1, n | 0); schedule(); },
  liveCount,
  stats,
  records: () => Array.from(recs.values()).map((r) => ({
    name: r.name, state: r.fig.getAttribute('data-state'), live: !!r.live, near: r.near,
    visible: r.visible, running: !!(r.live && r.live.running), t: r.live ? r.live.t : r.t, failed: r.failed,
  })),
  toggle: (el) => { const st = el && (el.classList.contains('anim-stage') ? el : el.querySelector('.anim-stage')); toggle(recs.get(st)); },
  pauseAll: () => recs.forEach((r) => { if (!r.failed) { r.intent = 'pause'; setState(r); } }),
  playAll: () => recs.forEach((r) => { if (!r.failed) { r.intent = 'play'; setState(r); } }),
  refreshTheme: onTheme,
  /** Jump a stage (element inside the figure, or the figure) to scene time t seconds; redraws when not running. */
  seek: (el, t) => {
    const st = el && (el.classList.contains('anim-stage') ? el : el.querySelector('.anim-stage'));
    const r = recs.get(st);
    if (!r) return false;
    r.t = t;
    if (r.live) { r.live.t = t; r.live.last = null; if (!r.live.running) r.live.draw(0); }
    return true;
  },
};
