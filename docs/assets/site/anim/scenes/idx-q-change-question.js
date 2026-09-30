/*
 * idx-q-change-question — index.html #problem, GAP 1 card "Change the question" (thumbnail, 16:9).
 * The problem only (how Experiment 1 answers it is idx-exp-edit-rewrite).
 * A question card sits over a gear: the fixed rule, hard-wired to the soft rule ("soft" is written on it).
 * The card is swapped from "soft" to "B/T-priority". The link to the gear breaks, the gear tries to turn
 * and does not, and its stamp keeps landing on every cell on the belt the same way. Cells the new question
 * would hold back (dashed ring, "hold back") are stamped anyway; each sends a dot to the tally, which shows only
 * 2,279, the page's count of cells the kept soft rule answers that B/T-priority holds back (no made-up partial
 * counts on the way). A note gives
 * the page's way out under the belt: until re-coded, or re-tuned on labels.
 * Which belt cells are flagged is illustrative (their share is not to scale); the only number is 2,279.
 * 12 s seamless loop; the still frame (reduced motion) shows the counter at 2,279.
 */
const NUMBERS = {
  // index.html #problem, GAP 1: "until then it answers 2,279 cells B/T-priority holds back"
  heldBackAnswered: 2279,
};

export default function create(ctx) {
  const { THREE, ease, seg, clamp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });

  /* timing (seconds; one stamp per beat) */
  const PERIOD = 12, N = 12, STILL = 9.25;
  const HITS = [4, 6, 8];                 // beats at which a cell the new question holds back is stamped
  const SWAP = [2.3, 2.8], BACK = [11.35, 11.85], FADE = [11.3, 11.9];
  const RAMP = [4.95, 8.95], WAY = [8.9, 9.2];

  /* layout (world units, 16 × 9) */
  const S = 5, SP = 1.4, SX = -1.6, X0 = SX - S * SP;          // stamp slot, spacing, stamp x, belt start
  const BY = -2.75, CY = BY + 0.35, R = 0.3, FR = 0.55, TRAVEL = 0.5;
  const BW = 1.0, BH = 0.46, HH = 0.34, KH = 0.18;             // stamp pad, handle, knob
  const QY = 3.5, QH = 1.2, GY = 1.2, RO = 1.15, RT = 0.93;    // question card, gear (tip / root radius)
  const COL = 0.7, NY = 1.65, RIGHT = 7.75;                    // right column: left edge, counter row, limit
  const TXT = 13.5, QTXT = 15, NTXT = 22;

  const add = (o, x = 0, y = 0, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: TXT, color: 'ink' }, o)));
  const circle = (r, n = 44) => Array.from({ length: n }, (_, i) => [r * Math.cos((2 * Math.PI * i) / n), r * Math.sin((2 * Math.PI * i) / n)]);

  /* belt and cells: a ring per cell; the stamp leaves the same teddy mark on every one */
  add(ctx.line([[-8.3, BY], [8.3, BY]], { color: 'line', width: 1.5 }));
  const cells = [];
  for (let j = 0; j < N; j++) {
    const g = add(ctx.group(), 0, CY, 1);
    const c = {
      g,
      ring: add(ctx.dot([0, 0], { r: R, hollow: true, ring: 0.3, color: 'muted' }), 0, 0, 0, g),
      mark: add(ctx.box(0.26, 0.26, { color: 'teddy', radius: 0.05 }), 0, 0, 0.1, g),
      flag: null, hold: null, tok: 'ink',
    };
    if (HITS.includes((((S - j) % N) + N) % N)) {
      c.flag = add(ctx.line(circle(FR), { color: 'ink', width: 1.5, dashed: [3, 3], closed: true }), 0, 0, 0.2, g);
      c.hold = add(ctx.label('hold back', { size: TXT, color: 'ink', anchor: 'top' }), 0, -(FR + 0.08), 0.3, g);
    }
    cells.push(c);
  }

  /* the question card */
  const card = add(ctx.box(4.4, QH, { color: 'soft', stroke: 'ink', strokeWidth: 1.5, radius: 0.2 }), SX, QY, 2);
  const qTag = lab('question', { color: 'muted', anchor: 'right' });
  const qSoft = lab('soft', { size: QTXT, weight: 600 });
  const qNew = lab('B/T-priority', { size: QTXT, weight: 600 });

  /* link card → gear: solid while they agree, broken once they do not */
  const linkPts = [[SX, QY - QH / 2, 1], [SX, GY + RO + 0.02, 1]];
  const link = add(ctx.line(linkPts, { color: 'ink', width: 1.6 }));
  const linkBad = add(ctx.line(linkPts, { color: 'bad', width: 1.6, dashed: [3, 3], opacity: 0 }));

  /* the gear: the fixed rule, wired to the soft rule */
  const gear = add(ctx.group(), SX, GY, 2);
  const gp = [];
  const TEETH = 10, dl = Math.PI / TEETH;
  for (let i = 0; i < TEETH; i++) {
    const c = i * 2 * dl + Math.PI / 2;
    for (const [r, a] of [[RT, c - 0.62 * dl], [RO, c - 0.3 * dl], [RO, c + 0.3 * dl], [RT, c + 0.62 * dl], [RT, c + dl]]) {
      gp.push(new THREE.Vector2(r * Math.cos(a), r * Math.sin(a)));
    }
  }
  const gMat = ctx.track(new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, opacity: 0.16, side: THREE.DoubleSide }));
  ctx.bind(gMat, 'teddy');
  gear.add(new THREE.Mesh(ctx.track(new THREE.ShapeGeometry(new THREE.Shape(gp))), gMat));
  add(ctx.line(gp.map((p) => [p.x, p.y]), { color: 'teddy', width: 1.6, closed: true }), 0, 0, 0.05, gear);
  add(ctx.label('soft', { size: TXT, weight: 600 }), SX, GY, 6);
  const ruleLab = lab('fixed rule', { color: 'teddy', weight: 600, anchor: 'right' });
  const wiredLab = lab('hard-wired', { color: 'muted', anchor: 'right' });
  /* the turn the new question would need; the gear only shivers */
  const AR = RO + 0.32, A0 = -0.95, A1 = 0.95;
  const turn = add(ctx.arrow([SX + AR * Math.cos(A0), GY + AR * Math.sin(A0), 2], [SX + AR * Math.cos(A1), GY + AR * Math.sin(A1), 2],
    { color: 'muted', width: 1.5, head: 7, bend: -0.515, dashed: [3, 3] }));

  /* the stamp, driven by a rod from the gear (origin at the bottom of its pad) */
  const rod = add(ctx.line([[SX, 0], [SX, -1]], { color: 'teddy', width: 3 }));
  const stamp = add(ctx.group(), SX, 0, 3);
  add(ctx.box(BW, BH, { color: 'teddy', opacity: 0.16, radius: 0.1 }), 0, BH / 2, 0, stamp);
  add(ctx.box(BW, BH, { color: null, stroke: 'teddy', radius: 0.1 }), 0, BH / 2, 0.01, stamp);
  add(ctx.box(BW - 0.1, 0.1, { color: 'teddy', radius: 0.03 }), 0, 0.05, 0.02, stamp);
  add(ctx.box(0.2, HH, { color: 'teddy', radius: 0.05 }), 0, BH + HH / 2, 0, stamp);
  add(ctx.box(0.56, KH, { color: 'teddy', radius: 0.09 }), 0, BH + HH + KH / 2, 0, stamp);

  /* the counter: a held-back cell that was stamped, and how many there are */
  const icon = add(ctx.group(), COL + 0.55, NY, 2);
  const iRing = add(ctx.dot([0, 0], { r: R, hollow: true, ring: 0.3, color: 'muted' }), 0, 0, 0, icon);
  const iMark = add(ctx.box(0.26, 0.26, { color: 'teddy', radius: 0.05 }), 0, 0, 0.1, icon);
  const iFlag = add(ctx.line(circle(FR), { color: 'bad', width: 1.5, dashed: [3, 3], closed: true }), 0, 0, 0.2, icon);
  const count = lab(NUMBERS.heldBackAnswered.toLocaleString('en-US'), { size: NTXT, weight: 700, color: 'bad', anchor: 'left' });
  const countLab = lab('answered anyway', { color: 'muted', anchor: 'top-left' });
  const WAYTXT = 'until re-coded, or re-tuned on labels';
  const way = lab(WAYTXT, { anchor: 'top' });   // under the belt, once the last 'hold back' has gone
  const fly = add(ctx.dot([0, 0], { r: 0.16, color: 'bad' }), 0, 0, 7);
  const P0 = [SX + 0.15, CY], P2 = [COL + 0.55, NY], P1 = [P2[0], CY + 0.3];

  /* text widths in CSS px, so the labels fit the ~250–420 px thumbnail */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, px, wt = 500) => {
    const lines = s.split('\n');
    if (!mc) return Math.max(...lines.map((l) => l.length)) * px * 0.56;
    mc.font = `${wt} ${px}px ${family}`;
    return Math.max(...lines.map((l) => mc.measureText(l).width));
  };
  const fit = (spr, text, span, base, wt = 500) => {   // span = world units the text may take
    const ppu = ctx.ppu(), ts = ctx.textScale, w = textW(text, base * ts, wt) + 4;
    spr.setSize(Math.max(12, Math.min(base, (base * span * ppu) / w)));
  };
  let dead = false;

  function layout() {
    const ppu = ctx.ppu(), ts = ctx.textScale, lh = (px) => (px * ts * 1.22 + 2) / ppu;
    const qw = Math.max(textW('soft', QTXT * ts, 600), textW('B/T-priority', QTXT * ts, 600)) / ppu;
    const cw = clamp(qw + 0.9, 3.6, 6);
    card.setSize(cw, QH, 0.2);
    qTag.position.set(SX - cw / 2 - 0.22, QY, 5);
    qSoft.position.x = SX; qNew.position.x = SX;
    const d = lh(TXT) / 2 + 0.02;
    ruleLab.position.set(SX - RO - 0.28, GY + d, 5);
    wiredLab.position.set(SX - RO - 0.28, GY - d, 5);
    const nx = COL + 1.3;
    count.position.set(nx, NY, 5);
    countLab.position.set(nx, NY - lh(NTXT) / 2 + 0.02, 5);
    fit(countLab, 'answered anyway', RIGHT - nx, TXT);
    way.position.set(0, CY - FR - 0.08, 5);
    fit(way, WAYTXT, 2 * RIGHT, TXT);
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  function update(t) {
    const u = ctx.loopT(t, PERIOD), tt = u * PERIOD;
    const b = Math.floor(tt), ph = tt - b;
    const endK = 1 - seg(tt, FADE[0], FADE[1], ease.inOutSine);

    /* the stamp comes down on every beat, whatever the question */
    const press = 1 - ease.smooth((Math.min(ph, 1 - ph) - 0.04) / 0.2);
    const padY = CY + R + TRAVEL * (1 - press);
    stamp.position.y = padY;
    rod.setPoints([[SX, GY - RT + 0.05, 1.5], [SX, padY + BH + HH + KH - 0.02, 1.5]]);

    /* the belt moves one slot between stamps */
    const e = ease.inOutSine((ph - 0.28) / 0.44);
    for (let j = 0; j < N; j++) {
      const c = cells[j], k = (j + b) % N, p = (j + b + e) % N;
      c.g.position.x = X0 + p * SP;
      const a = clamp(p / 0.7) * clamp((11.5 - p) / 0.6);
      c.ring.setOpacity(a);
      const at = k === S && ph < 0.5;
      c.mark.visible = at || p > S + 1e-4;
      c.mark.scale.setScalar(at ? Math.max(0.01, ease.outBack(ph / 0.09)) : 1);
      c.mark.setOpacity(a);
      if (c.flag) {
        const tHit = b - k + S;   // the beat this trip is stamped; outside 0..N-1 = another loop
        const on = tHit >= 0 && tHit < N ? seg(tt, tHit - 1.25, tHit - 0.95) * endK : 0;
        const tok = tt >= tHit ? 'bad' : 'ink';
        if (tok !== c.tok) { c.flag.setColor(tok); c.tok = tok; }
        c.flag.setOpacity(on * a);
        c.hold.setOpacity(on * a * (1 - seg(tt, tHit + 0.25, tHit + 0.55)));
      }
    }

    /* the question changes; the gear does not */
    const qk = seg(tt, SWAP[0], SWAP[1], ease.inOutSine) * (1 - seg(tt, BACK[0], BACK[1], ease.inOutSine));
    qSoft.setOpacity(1 - seg(qk, 0, 0.55)); qSoft.position.y = QY + 0.2 * qk;
    qNew.setOpacity(seg(qk, 0.45, 1)); qNew.position.y = QY - 0.2 * (1 - qk);
    const lk = seg(tt, 2.7, 3.0) * endK;
    link.setOpacity(1 - lk); linkBad.setOpacity(lk);
    turn.setProgress(seg(tt, 2.85, 3.3, ease.outCubic));
    turn.setOpacity(1 - seg(tt, 3.7, 4.1));
    const tau = tt - 3.05;
    gear.rotation.z = tau > 0 && tau < 0.75 ? 0.06 * Math.sin(2 * Math.PI * 4 * tau) * (1 - tau / 0.75) : 0;

    /* counter: every held-back cell the stamp answers sends a dot; the page's count appears with the first one */
    const ck = seg(tt, 2.95, 3.35) * endK;
    count.setOpacity(seg(tt, RAMP[0], RAMP[0] + 0.45) * endK); countLab.setOpacity(ck);
    iRing.setOpacity(ck); iMark.setOpacity(ck); iFlag.setOpacity(ck);
    let fo = 0, bump = 0;
    for (const h of HITS) {
      const f = (tt - (h + 0.3)) / 0.65;
      if (f >= 0 && f <= 1) {
        const s = ease.inOutSine(f), r = 1 - s;
        fly.position.set(r * r * P0[0] + 2 * r * s * P1[0] + s * s * P2[0], r * r * P0[1] + 2 * r * s * P1[1] + s * s * P2[1], 7);
        fo = seg(f, 0, 0.12) * (1 - seg(f, 0.85, 1));
      }
      const x = (tt - (h + 0.9)) / 0.35;
      if (x > 0 && x < 1) bump = Math.sin(Math.PI * x);
    }
    fly.setOpacity(fo);
    icon.scale.setScalar(1 + 0.25 * bump);
    way.setOpacity(seg(tt, WAY[0], WAY[1]) * endK);
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: STILL,
    update, resize: layout,
    dispose() { dead = true; cells.length = 0; },
  };
}
