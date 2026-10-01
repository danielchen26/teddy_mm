/*
 * idx-q-price-of-question — index.html #problem, GAP 5 card "Price of a new question" (thumbnail, 16:9).
 * The problem only (how Experiment 5 answers it is idx-exp-labels-cost-curve).
 * Question cards queue on a lane in front of a trained readout, whose door stays shut until the question is
 * paid for in labelled cells. Stacks of label coins slide from a pile into the coin jar; once it is full the
 * link to the door lights, the booth runs the page's two follow-on steps ("retraining", "recalibrating"),
 * the door swings open and the card goes through. Then the jar empties, the next card steps up and the
 * readout asks for labels again ("pay again").
 * The price tag carries the card's only number, "50–2,000 labels", with the card's "again at every change";
 * stack sizes and the jar's fill are illustrative, not counts. First card "B/T-priority" (the question the
 * card asks about), then "next question".
 * 12 s seamless loop (two payments, then the lane clears); the still frame shows the second payment under way
 * ("next question" at the door, "pay again", a stack sliding into the jar).
 */
const NUMBERS = {
  // index.html #problem, GAP 5: "Trained readout — v2: 50–1,000 new labelled cells, again at every change." (v3 E1.5: n* = 0)
  labelsLow: 50,
  labelsHigh: 1000,
};

export default function create(ctx) {
  const { THREE, ease, seg, clamp, lerp } = ctx;
  const scene = new THREE.Scene();
  const camera = ctx.orthoCamera({ width: 16, height: 9 });

  /* timing: seconds inside one payment (s = t mod 6); the loop holds two payments */
  const PERIOD = 12, HALF = 6, STILL = 7.55;
  const STACKS = [[0.6, 1.2], [1.2, 1.9], [1.9, 2.7]];   // each stack's slide into the jar
  const COINS = [3, 5, 7], TOTAL = 15;
  const LINK = [2.7, 3.1], OPEN = [4.1, 4.5], PASS = [4.4, 5.5], CLOSE = [5.45, 5.85];
  const DRAIN = [5.5, 5.95], UNLIT = [5.5, 5.9], STEP = [4.8, 5.8];
  const STATUS = [[0, 'closed', 'muted'], [3.1, 'retraining', 'ink'], [3.6, 'recalibrating', 'ink'],
    [4.1, 'open', 'train'], [5.45, 'closed', 'muted']];

  /* layout (world units, 16 × 9) */
  const GX = 1.2;                                          // door, link and jar share this x
  const FLOOR = 0.05, CH = 1.0, CY = FLOOR + 0.08 + CH / 2;
  const HINGE = 2.05, FLAP = HINGE - (FLOOR + 0.05);
  const BOOTH_W = 5.4, BOOTH_H = 1.6, BOOTH_Y = HINGE + BOOTH_H / 2;
  const BX = GX + FLAP + 0.15 - BOOTH_W / 2;               // the open door folds under the booth's right side
  const JW = 2.0, JH = 1.75, JY = -2.775, JT = JY + JH / 2, JB = JY - JH / 2;
  const CW = 0.9, COIN = 0.19, PITCH = 0.23, RAIL = JB + 0.07;
  const SX0 = -5.0, PILE = [[-7.1, 4], [-6.1, 6]];
  const TAG_X = GX + JW / 2 + 0.3, RIGHT = 7.85, D = 6.2;
  const TXT = 13, QTXT = 12.5;

  const add = (o, x = 0, y = 0, z = 0, parent = scene) => { o.position.set(x, y, z); parent.add(o); return o; };
  const lab = (text, o = {}) => add(ctx.label(text, Object.assign({ size: TXT, color: 'ink' }, o)));

  /* the lane the questions queue on */
  add(ctx.line([[-8.3, FLOOR], [8.3, FLOOR]], { color: 'line', width: 1.5 }));
  const queueLab = lab('questions', { color: 'muted', anchor: 'bottom' });

  /* the trained readout: a booth that shows what it is doing, and a door hinged under it */
  add(ctx.box(BOOTH_W, BOOTH_H, { color: 'soft', stroke: 'train', strokeWidth: 1.6, radius: 0.2 }), BX, BOOTH_Y, 5);
  const title = lab('trained readout', { color: 'train', weight: 600 });
  const status = lab('closed', { color: 'muted' });
  const door = add(ctx.group(), GX, HINGE, 4);
  add(ctx.line([[0, 0], [0, -FLAP]], { color: 'train', width: 4 }), 0, 0, 0, door);
  add(ctx.dot([0, 0], { px: 3.5, color: 'train' }), 0, 0, 0.1, door);

  /* question cards: K0 "B/T-priority" pays first, K1 "next question" pays second, K2 waits behind */
  const TEXTS = ['B/T-priority', 'next question', 'next question'];
  const cards = TEXTS.map((text) => {
    const g = add(ctx.group(), 0, CY, 3);
    const b = add(ctx.box(4, CH, { color: 'soft', stroke: 'ink', strokeWidth: 1.4, radius: 0.18 }), 0, 0, 0, g);
    const l = add(ctx.label(text, { size: QTXT, weight: 600 }), 0, 0, 0.5, g);
    return { g, b, l };
  });

  /* the coin jar the labels go into, and the link from it to the door */
  add(ctx.line([[GX, JT, 0.5], [GX, FLOOR, 0.5]], { color: 'faint', width: 1.4, dashed: [3, 4] }));
  const link = add(ctx.line([[GX, JT, 0.6], [GX, FLOOR, 0.6]], { color: 'train', width: 2.4 }));
  add(ctx.box(JW, JH, { color: 'card', stroke: 'train', strokeWidth: 1.6, radius: 0.16 }), GX, JY, 2);
  const fill = add(ctx.box(JW - 0.24, JH - 0.24, { color: 'train', opacity: 0.3, radius: 0.06 }), GX, JY, 2.2);
  const FH = JH - 0.24, FB = JB + 0.12;
  const TAG1 = `${NUMBERS.labelsLow}–${NUMBERS.labelsHigh.toLocaleString('en-US')} labels`, TAG2 = 'v2: again each change';
  const tag = lab(TAG1, { weight: 600, anchor: 'left' });
  const tag2 = lab(TAG2, { color: 'muted', anchor: 'left' });
  const prompt = lab('pay to pass', { color: 'bad', weight: 600, anchor: 'right' });

  /* labelled cells as coins: a pile, and three stacks that slide from it into the jar */
  const coin = (parent, y, token) => {   // a coin with a dot for its label (no lineage colours: nothing is called here)
    add(ctx.box(CW, COIN, { color: 'card', stroke: token, strokeWidth: 1.2, radius: 0.06 }), 0, y, 0, parent);
    add(ctx.dot([0, 0], { r: 0.05, color: token }), -CW / 2 + 0.17, y, 0.05, parent);
  };
  for (const [x, n] of PILE) {
    const g = add(ctx.group(), x, RAIL + COIN / 2, 1);
    for (let i = 0; i < n; i++) coin(g, i * PITCH, 'muted');
  }
  const pileLab = lab('labelled cells', { color: 'muted', anchor: 'bottom-left' });
  const stacks = COINS.map((n) => {
    const g = add(ctx.group(), SX0, RAIL + COIN / 2, 1);
    for (let i = 0; i < n; i++) coin(g, i * PITCH, 'measured');
    const parts = [];
    g.traverse((o) => { if (o !== g && o.setOpacity) parts.push(o); });
    return { g, parts, o: -1 };
  });

  /* text widths in CSS px, so the labels fit the ~300–360 px thumbnail */
  const mc = document.createElement('canvas').getContext('2d');
  const family = getComputedStyle(ctx.figure).fontFamily || 'sans-serif';
  const textW = (s, px, wt = 500) => {
    if (!mc) return s.length * px * 0.56;
    mc.font = `${wt} ${px}px ${family}`;
    return mc.measureText(s).width;
  };
  const fit = (spr, text, span, base, wt = 500, min = 11) => {   // span = world units the text may take
    const w = textW(text, base * ctx.textScale, wt) + 4;
    spr.setSize(Math.max(min, Math.min(base, (base * span * ctx.ppu()) / w)));
  };
  const L = { W: 4, XG: -1.0, XW: -5.3 };
  let dead = false;

  function layout() {
    const ppu = ctx.ppu(), ts = ctx.textScale, lh = (px) => (px * ts * 1.22 + 2) / ppu;
    const tw = Math.max(...TEXTS.map((s) => textW(s, QTXT * ts, 600))) / ppu;
    L.W = clamp(tw + 0.5, 1.8, 4.3);
    L.XG = GX - 0.2 - L.W / 2;
    L.XW = L.XG - L.W - 0.3;
    cards.forEach((c) => c.b.setSize(L.W, CH, 0.18));
    queueLab.position.set(L.XW, CY + CH / 2 + 0.1, 6);
    const d = lh(TXT) / 2 + 0.02;
    title.position.set(BX, BOOTH_Y + d, 6);
    status.position.set(BX, BOOTH_Y - d, 6);
    fill.setOpacity(ctx.theme === 'dark' ? 0.5 : 0.3);
    tag.position.set(TAG_X, JY + d + 0.02, 6);
    tag2.position.set(TAG_X, JY - d - 0.02, 6);
    fit(tag, TAG1, RIGHT - TAG_X, TXT, 600);
    fit(tag2, TAG2, RIGHT - TAG_X, TXT, 500, 10); // the muted second line may go a size smaller: at 11 it ran off the ~270 px thumbnail
    prompt.position.set(GX - 0.3, (JT + FLOOR) / 2, 6);
    pileLab.position.set(PILE[0][0] - CW / 2, RAIL + (PILE[1][1] - 1) * PITCH + COIN + 0.16, 6);
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { if (!dead) { layout(); ctx.requestRender(); } });

  const place = (c, x, o, dy = 0) => {
    c.g.position.set(x, CY + dy, 3);
    c.b.setOpacity(o); c.l.setOpacity(o);
    c.g.visible = o > 0.002;
  };
  let stTok = 'muted';

  function update(t) {
    const u = ctx.loopT(t, PERIOD), tt = u * PERIOD, s = tt % HALF;

    /* the stacks slide in and vanish behind the jar; the jar fills, then empties for the next card */
    let lvl = 0;
    STACKS.forEach(([a, b], i) => {
      const st = stacks[i];
      st.g.position.x = lerp(SX0, GX, seg(s, a, b, ease.inOutSine));
      const o = s < a || s > b + 0.05 ? 0 : seg(s, a, a + 0.12);
      if (o !== st.o) { st.parts.forEach((p) => p.setOpacity(o)); st.o = o; st.g.visible = o > 0; }
      lvl += (COINS[i] / TOTAL) * seg(s, lerp(a, b, 0.62), b + 0.08, ease.outCubic);
    });
    lvl *= 1 - seg(s, DRAIN[0], DRAIN[1], ease.inOutSine);
    fill.visible = lvl > 0.004;
    fill.scale.y = Math.max(0.004, lvl);
    fill.position.y = FB + (lvl * FH) / 2;

    /* once paid: the link lights, the booth retrains and recalibrates, and only then the door opens */
    link.setProgress(seg(s, LINK[0], LINK[1])).setOpacity(1 - seg(s, UNLIT[0], UNLIT[1]));
    let k = 0;
    while (k + 1 < STATUS.length && s >= STATUS[k + 1][0]) k++;
    const [, sText, sTok] = STATUS[k];
    status.setText(sText);
    if (sTok !== stTok) { status.setColor(sTok); stTok = sTok; }
    let dip = 1;
    for (let i = 1; i < STATUS.length; i++) dip = Math.min(dip, clamp(Math.abs(s - STATUS[i][0]) / 0.1));
    status.setOpacity(ease.smooth(dip));
    door.rotation.z = (Math.PI / 2) * seg(s, OPEN[0], OPEN[1], ease.outCubic) * (1 - seg(s, CLOSE[0], CLOSE[1], ease.inOutSine));

    /* the price is asked for when a card reaches the door: first "pay to pass", then "pay again" */
    prompt.setText(tt >= 5.5 && tt < 11.5 ? 'pay again' : 'pay to pass');
    prompt.setOpacity(Math.max(seg(tt, 0.3, 0.6) * (1 - seg(tt, 2.7, 3.0)), seg(tt, 5.75, 6.05) * (1 - seg(tt, 8.7, 9.0))));

    /* the cards: K0 pays and passes, K1 steps up and pays again, K2 joins the queue; the lane clears */
    const pass = (a) => D * seg(tt, a, a + PASS[1] - PASS[0], ease.inOutSine);
    const intro = 0.4 * (1 - seg(tt, 0, 0.45, ease.outCubic));   // the first two cards settle onto the lane
    const step = seg(tt, STEP[0], STEP[1], ease.inOutSine), gap = L.XG - L.XW;
    place(cards[0], L.XG + pass(PASS[0]), seg(tt, 0, 0.4) * (1 - seg(tt, 5.1, 5.5)), intro);
    place(cards[1], L.XW + gap * step + pass(PASS[0] + HALF), seg(tt, 0, 0.4) * (1 - seg(tt, 11.1, 11.5)), intro);
    const k2 = seg(tt, STEP[1] - 0.15, STEP[1] + 0.35, ease.outCubic);   // joins the queue once K1 has moved up
    place(cards[2], L.XW - 0.6 * (1 - k2), k2 * (1 - seg(tt, 11.3, 11.8)));
  }

  layout();
  return {
    scene, camera, period: PERIOD, still: STILL,
    update, resize: layout, onTheme: layout,
    dispose() { dead = true; cards.length = 0; stacks.length = 0; },
  };
}
