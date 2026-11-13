"""Charts for the site (scripts/build_site.py): the sealed chain, u at Z, head latent and head output for each head (one
panel per cell type), and the bone-marrow score of each head against the guard. Plain SVG strings; colours are CSS tokens
of the page. `tx` writes one label in both languages."""
ARMS = [("as_trained", "第 0 轮：原来的头", "round 0: head as trained", "var(--muted)", "6 4", "c"), ("decoder_only", "对照：只修解码器（定位排除的）", "control: decoder repaired (ruled out)", "var(--s2)", "", "c"),
        ("mlp_only", "第 1 轮：只修 MLP（按定位）", "round 1: MLP repaired (as located)", "var(--s1)", "", "c"), ("mlp_guarded", "第 2 轮：再按住骨髓上的输出", "round 2: bone-marrow outputs held", "var(--accent)", "", "s"),
        ("joint_s0", "第 3 轮：两种细胞、全部实测蛋白一起修", "round 3: both cell types, all measured proteins", "var(--s3)", "", "d")]
STAGES = [("Z", "Z", "Z"), ("HEAD_LAT", "头的中间层", "head latent"), ("HEAD_OUT", "头的输出", "head output")]
DS = [("frangieh2021", "黑色素瘤", "Melanoma"), ("papalexi2021", "THP-1", "THP-1")]
def mark(shape, x, y, c, r=4.6):
    if shape == "d":
        return f'<path d="M{x:.1f} {y - r - 1:.1f} L{x + r + 1:.1f} {y:.1f} L{x:.1f} {y + r + 1:.1f} L{x - r - 1:.1f} {y:.1f} Z" style="fill:{c};stroke:var(--bg)" stroke-width="1.4"/>'
    if shape == "s":
        return f'<rect x="{x - r:.1f}" y="{y - r:.1f}" width="{2 * r}" height="{2 * r}" style="fill:{c};stroke:var(--bg)" stroke-width="1.4"/>'
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" style="fill:{c};stroke:var(--bg)" stroke-width="1.4"/>'
def chart(tx, A, readout="PDL1", kind="pooled", colors=None):
    """A: audits by name (loop.json['audits'], plus the guarded ones if present). Lines join the three stages; whiskers are the intervals."""
    col = dict((a[0], a[3]) for a in ARMS); col.update(colors or {})
    arms = [a for a in ARMS if all(f"{kind}_{ds}_{a[0]}" in A for ds, _, _ in DS)]
    rows_l = (len(arms) + 1) // 2 if len(arms) > 3 else 1
    W, top, bot, pl, pw, gap = 760, 62, 52 + 26 * rows_l, 56, 290, 68; H = top + 200 + bot
    vals = [A[f"{kind}_{ds}_{a[0]}"]["readouts"][readout]["stages"][s[0]] for ds, _, _ in DS for a in arms for s in STAGES]
    ymax = max(0.6, min(1.6, max(v["ci"][1] for v in vals) + 0.05)); ymin = min(-0.1, max(-0.4, min(v["ci"][0] for v in vals) - 0.05))
    Y = lambda v: top + (H - top - bot) * (1 - (min(max(v, ymin), ymax) - ymin) / (ymax - ymin))
    desc_zh, desc_en, s = [], [], []
    for pi, (ds, zh, en) in enumerate(DS):
        x0 = pl + pi * (pw + gap); X = lambda i: x0 + 40 + i * (pw - 80) / 2
        s.append(f'<rect x="{x0}" y="{Y(0.10):.1f}" width="{pw}" height="{Y(0) - Y(0.10):.1f}" style="fill:var(--accent-tint)"/>')
        v = round(round(ymin / 0.2) * 0.2, 1)
        while v <= ymax + 1e-9:
            if v >= ymin - 1e-9:
                s.append(f'<line x1="{x0}" y1="{Y(v):.1f}" x2="{x0 + pw}" y2="{Y(v):.1f}" stroke="var(--line)" stroke-width="{1.2 if abs(v) < 1e-9 else .7}"/><text x="{x0 - 8}" y="{Y(v) + 4:.1f}" font-size="11" text-anchor="end" style="fill:var(--muted);font-family:var(--mono)">{v:.1f}</text>')
            v = round(v + 0.2, 1)
        s.append(tx(x0, top - 14, zh, en, "start", 13, weight=600))
        for i, (sid, szh, sen) in enumerate(STAGES):
            s.append(tx(X(i), H - bot + 20, szh, sen, "middle", 11.5))
        labels = []
        for ai, (arm, azh, aen, _, dash, shape) in enumerate(arms):
            st = A[f"{kind}_{ds}_{arm}"]["readouts"][readout]["stages"]; c = col[arm]; off = (ai - (len(arms) - 1) / 2) * 7
            pts = " ".join(f"{X(i) + off:.1f},{Y(st[sid]['u']):.1f}" for i, (sid, _, _) in enumerate(STAGES))
            s.append(f'<polyline points="{pts}" fill="none" style="stroke:{c}" stroke-width="2.2"' + (f' stroke-dasharray="{dash}"' if dash else "") + '/>')
            for i, (sid, _, _) in enumerate(STAGES):
                v = st[sid]; x = X(i) + off
                s.append(f'<line x1="{x:.1f}" y1="{Y(v["ci"][0]):.1f}" x2="{x:.1f}" y2="{Y(v["ci"][1]):.1f}" style="stroke:{c}" stroke-width="1.4" opacity=".75"/>' + mark(shape, x, Y(v["u"]), c))
            labels.append([Y(st["HEAD_OUT"]["u"]) + 4, st["HEAD_OUT"]["u"]])
            desc_zh.append(f"{zh}，{azh}：" + "、".join(f"{st[sid]['u']:.2f}" for sid, _, _ in STAGES)); desc_en.append(f"{en}, {aen}: " + ", ".join(f"{st[sid]['u']:.2f}" for sid, _, _ in STAGES))
        labels.sort()                                                          # value labels at the right edge, pushed apart so they never overlap
        for k in range(1, len(labels)):
            labels[k][0] = max(labels[k][0], labels[k - 1][0] + 12)
        for y, u in labels:
            s.append(f'<text x="{X(2) + (len(arms) - 1) * 3.5 + 10:.1f}" y="{y:.1f}" font-size="11" style="fill:var(--fg, currentColor);font-family:var(--mono)">{u:.2f}</text>')
    for k, (arm, azh, aen, _, dash, shape) in enumerate(arms):
        c = col[arm]; lx = pl + (k % 2) * (pw + gap) if len(arms) > 3 else pl + k * 232; ly = H - 22 - 26 * (rows_l - 1 - k // 2) if len(arms) > 3 else H - 22
        s.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 26}" y2="{ly}" style="stroke:{c}" stroke-width="2.4"' + (f' stroke-dasharray="{dash}"' if dash else "") + '/>' + mark(shape, lx + 13, ly, c, 4.2) + tx(lx + 33, ly + 4, azh, aen, "start", 11.5))
    s.append(tx(pl - 44, 14, "读不出的比例 u（线段是区间，浅色带是容差 0.10）", "Unreadable share u (bars are intervals; the pale band is the 0.10 margin)", "start", 11.5))
    head = f'<svg style="min-width:640px" viewBox="0 0 {W} {H}" role="img" data-al-zh="各个头在 Z、头的中间层、头的输出三环上读不出的比例。{"；".join(desc_zh)}。" data-al-en="Unreadable share at Z, head latent and head output for each head. {"; ".join(desc_en)}.">'
    return head + "".join(s) + "</svg>"
def cost(tx, rows, base, guard=0.02, colors=None):
    """rows: [(arm, zh, en, lo, hi)] bone-marrow scores of a group of heads; base: the head as trained. A line marks the guard."""
    col = dict((a[0], a[3]) for a in ARMS); col.update(colors or {}); shp = dict((a[0], a[5]) for a in ARMS)
    W, left, right, top, rowh = 760, 300, 46, 34, 26; H = top + rowh * (len(rows) + 1) + 46
    x0 = min(0.34, min(r[3] for r in rows) - 0.02); x1 = base + 0.085; X = lambda v: left + (W - left - right) * (v - x0) / (x1 - x0)
    d_zh = "；".join(f"{r[1]} {r[3]:.3f} 到 {r[4]:.3f}" for r in rows); d_en = "; ".join(f"{r[2]} {r[3]:.3f} to {r[4]:.3f}" for r in rows)
    s = [f'<svg style="min-width:620px" viewBox="0 0 {W} {H}" role="img" data-al-zh="各个头在骨髓测试细胞上的得分。原来的头 {base:.3f}，界限 {base - guard:.3f}。{d_zh}。" data-al-en="Score of each head on bone-marrow test cells. Head as trained {base:.3f}, guard {base - guard:.3f}. {d_en}.">']
    yb = top + rowh * (len(rows) + 1)
    s.append(f'<rect x="{X(base - guard):.1f}" y="{top - 6}" width="{X(base) - X(base - guard):.1f}" height="{yb - top + 6}" style="fill:var(--accent-tint)"/>')
    v = round((int(x0 / 0.05) + 1) * 0.05, 2)
    while v <= x1:
        s.append(f'<line x1="{X(v):.1f}" y1="{top - 6}" x2="{X(v):.1f}" y2="{yb}" stroke="var(--line)" stroke-width=".7"/><text x="{X(v):.1f}" y="{yb + 15}" font-size="11" text-anchor="middle" style="fill:var(--muted);font-family:var(--mono)">{v:.2f}</text>'); v = round(v + 0.05, 2)
    s.append(f'<line x1="{X(base - guard):.1f}" y1="{top - 6}" x2="{X(base - guard):.1f}" y2="{yb}" style="stroke:var(--accent)" stroke-width="1.6"/>' + tx(X(base - guard) - 6, top - 12, f"界限 {base - guard:.3f}", f"guard {base - guard:.3f}", "end", 10.5))
    cy = top + rowh / 2
    s.append(tx(left - 12, cy + 4, "原来的头", "head as trained", "end", 11.5) + f'<line x1="{X(base):.1f}" y1="{top - 6}" x2="{X(base):.1f}" y2="{yb}" style="stroke:var(--muted)" stroke-width="1.4" stroke-dasharray="4 3"/>' + mark("c", X(base), cy, "var(--muted)", 5) + f'<text x="{X(base) + 9:.1f}" y="{cy + 4}" font-size="11" style="fill:var(--fg, currentColor);font-family:var(--mono)">{base:.3f}</text>')
    for k, (arm, zh, en, lo, hi) in enumerate(rows):
        cy = top + rowh * (k + 1.5); c = col[arm]
        s.append(tx(left - 12, cy + 4, zh, en, "end", 11.5) + f'<line x1="{X(lo):.1f}" y1="{cy}" x2="{X(hi):.1f}" y2="{cy}" style="stroke:{c}" stroke-width="3" stroke-linecap="round"/>' + mark(shp[arm], X(lo), cy, c, 4.4) + mark(shp[arm], X(hi), cy, c, 4.4)
                 + (f'<text x="{X(lo) - 9:.1f}" y="{cy + 4}" font-size="11" text-anchor="end" style="fill:var(--fg, currentColor);font-family:var(--mono)">{lo:.3f}–{hi:.3f}</text>' if X(lo) - 9 - 84 > left + 6 and lo < base - guard
                    else f'<text x="{(max(X(hi), X(base)) if lo >= base - guard else X(hi)) + 9:.1f}" y="{cy + 4}" font-size="11" style="fill:var(--fg, currentColor);font-family:var(--mono)">{lo:.3f}–{hi:.3f}</text>'))
    s.append(tx((left + W - right) / 2, yb + 34, "骨髓测试细胞上的得分（每行是一组头的最低到最高）", "Score on bone-marrow test cells (each row: lowest to highest of a group of heads)", "middle", 11.5))
    return "".join(s) + "</svg>"


def share(tx, rows, arms=(("mlp_only", "修 MLP（定位到的那一步）", "MLP repaired (the located step)"), ("decoder_only", "只修解码器（对照）", "decoder repaired (control)")), colors=None):
    """rows: [(zh, en, {arm: (point, lo, hi)})]: the share of the baseline's added loss that each change removes, with intervals."""
    col = dict((a[0], a[3]) for a in ARMS); col.update(colors or {}); shp = dict((a[0], a[5]) for a in ARMS)
    W, left, right, top, band = 760, 168, 150, 30, 50; H = top + band * len(rows) + 74
    x0, x1 = -0.25, 1.3; X = lambda v: left + (W - left - right) * (min(max(v, x0), x1) - x0) / (x1 - x0); yb = top + band * len(rows)
    pc = lambda v: f"{100 * v:.0f}%".replace("-0%", "0%").replace("-", "−")
    d_zh = "；".join(f"{r[0]}：" + "，".join(f"{a[1]} {pc(r[2][a[0]][0])}（{pc(r[2][a[0]][1])} 到 {pc(r[2][a[0]][2])}）" for a in arms) for r in rows)
    d_en = "; ".join(f"{r[1]}: " + ", ".join(f"{a[2]} {pc(r[2][a[0]][0])} ({pc(r[2][a[0]][1])} to {pc(r[2][a[0]][2])})" for a in arms) for r in rows)
    s = [f'<svg style="min-width:620px" viewBox="0 0 {W} {H}" role="img" data-al-zh="头多丢的那部分响应被每种修法去掉的比例。{d_zh}。" data-al-en="Share of the loss the head adds that each repair removes. {d_en}.">']
    for v in (0, 0.25, 0.5, 0.75, 1.0):
        s.append(f'<line x1="{X(v):.1f}" y1="{top - 4}" x2="{X(v):.1f}" y2="{yb}" stroke="var(--line)" stroke-width="{1.3 if v in (0, 1.0) else .7}"/><text x="{X(v):.1f}" y="{yb + 15}" font-size="11" text-anchor="middle" style="fill:var(--muted);font-family:var(--mono)">{pc(v)}</text>')
    for k, (zh, en, vals) in enumerate(rows):
        cy = top + band * (k + 0.5)
        if k:
            s.append(f'<line x1="{left - 150}" y1="{top + band * k}" x2="{W - 12}" y2="{top + band * k}" stroke="var(--line)" stroke-width=".7"/>')
        s.append(tx(left - 14, cy + 4, zh, en, "end", 12, weight=600))
        for j, (arm, _, _) in enumerate(arms):
            pt, lo, hi = vals[arm]; y = cy + (j - (len(arms) - 1) / 2) * 17; c = col[arm]
            s.append(f'<line x1="{X(lo):.1f}" y1="{y:.1f}" x2="{X(hi):.1f}" y2="{y:.1f}" style="stroke:{c}" stroke-width="2.6" stroke-linecap="round"/>' + mark(shp[arm], X(pt), y, c, 4.6)
                     + f'<text x="{X(hi) + 10:.1f}" y="{y + 4:.1f}" font-size="11" style="fill:var(--fg, currentColor);font-family:var(--mono)">{pc(pt)} [{pc(lo)}, {pc(hi)}]</text>')
    s.append(tx((left + W - right) / 2, yb + 34, "头在 Z 到输出之间多丢的部分，被去掉的比例", "Share removed of the loss the head adds between Z and its output", "middle", 11.5))
    for j, (arm, zh, en) in enumerate(arms):
        lx = left + j * 270; ly = yb + 58; c = col[arm]
        s.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 26}" y2="{ly}" style="stroke:{c}" stroke-width="2.6"/>' + mark(shp[arm], lx + 13, ly, c, 4.2) + tx(lx + 33, ly + 4, zh, en, "start", 11.5))
    return "".join(s) + "</svg>"

CHAIN = [("RAW", "归一化的 RNA", "normalized RNA"), ("TOKENS", "TEDDY 的词元", "TEDDY's tokens"), ("Z", "嵌入 Z", "embedding Z"), ("HEAD_LAT", "头的中间层", "head latent"), ("HEAD_OUT", "头的输出", "head output")]
VC = {"support": "var(--good)", "reject": "var(--bad)", "unresolved": "var(--warn)"}
def chain(tx, panels, readout="PDL1"):
    """panels: [(zh, en, record)] sealed records; u with its interval at each of the five chain stages; the marker colour is the verdict."""
    W, top, bot, pl, pw, gap = 760, 60, 86, 56, 300, 48; H = top + 210 + bot
    vals = [r["readouts"][readout]["states"][s[0]] for _, _, r in panels for s in CHAIN]
    ymax = max(0.6, max(v["ci"][1] for v in vals) + 0.05); ymin = min(-0.1, min(v["ci"][0] for v in vals) - 0.05)
    Y = lambda v: top + (H - top - bot) * (1 - (v - ymin) / (ymax - ymin)); s, dz, de = [], [], []
    for pi, (zh, en, r) in enumerate(panels):
        st = r["readouts"][readout]["states"]; x0 = pl + pi * (pw + gap); X = lambda i: x0 + 26 + i * (pw - 52) / 4
        s.append(f'<rect x="{x0}" y="{Y(0.10):.1f}" width="{pw}" height="{Y(-0.10) - Y(0.10) if ymin <= -0.10 else Y(ymin) - Y(0.10):.1f}" style="fill:var(--accent-tint)"/>')
        v = round(round(ymin / 0.2) * 0.2, 1)
        while v <= ymax + 1e-9:
            if v >= ymin - 1e-9:
                s.append(f'<line x1="{x0}" y1="{Y(v):.1f}" x2="{x0 + pw}" y2="{Y(v):.1f}" stroke="var(--line)" stroke-width="{1.2 if abs(v) < 1e-9 else .7}"/><text x="{x0 - 8}" y="{Y(v) + 4:.1f}" font-size="11" text-anchor="end" style="fill:var(--muted);font-family:var(--mono)">{v:.1f}</text>')
            v = round(v + 0.2, 1)
        s.append(tx(x0, top - 14, zh, en, "start", 13, weight=600))
        s.append(f'<line x1="{X(2.5):.1f}" y1="{top}" x2="{X(2.5):.1f}" y2="{H - bot}" style="stroke:var(--muted)" stroke-width="1" stroke-dasharray="3 4"/>' + tx(X(2.5) - 6, top + 12, "TEDDY", "TEDDY", "end", 10.5) + tx(X(2.5) + 6, top + 12, "头", "head", "start", 10.5))
        s.append('<polyline points="' + " ".join(f"{X(i):.1f},{Y(st[k]['u']):.1f}" for i, (k, _, _) in enumerate(CHAIN)) + '" fill="none" style="stroke:var(--fg, currentColor)" stroke-width="1.8"/>')
        for i, (k, szh, sen) in enumerate(CHAIN):
            v = st[k]; c = VC[v["verdict"]]
            s.append(f'<line x1="{X(i):.1f}" y1="{Y(v["ci"][0]):.1f}" x2="{X(i):.1f}" y2="{Y(v["ci"][1]):.1f}" style="stroke:{c}" stroke-width="2"/>' + mark("c", X(i), Y(v["u"]), c, 5)
                     + f'<text x="{X(i):.1f}" y="{Y(v["ci"][1]) - 7:.1f}" font-size="11" text-anchor="middle" style="fill:var(--fg, currentColor);font-family:var(--mono)">{f"{v['u']:.2f}".replace("-", "−")}</text>')
            for cls, label in (("zh", szh), ("en", sen)):                      # one or two lines per language
                for li, part in enumerate(label.split(" ", 1)):
                    s.append(f'<text class="{cls}" x="{X(i):.1f}" y="{H - bot + 18 + 13 * li}" font-size="10.5" text-anchor="middle">{part}</text>')
        dz.append(f"{zh}：" + "、".join(f"{k} {st[k]['u']:.2f}" for k, _, _ in CHAIN)); de.append(f"{en}: " + ", ".join(f"{k} {st[k]['u']:.2f}" for k, _, _ in CHAIN))
    for k, (v, zh, en) in enumerate((("support", "支持", "support"), ("unresolved", "未决", "unresolved"), ("reject", "否定", "reject"))):
        lx = pl + k * 130; s.append(mark("c", lx + 6, H - 20, VC[v], 5) + tx(lx + 18, H - 16, zh, en, "start", 11.5))
    s.append(tx(W - 12, H - 16, "竖线是区间，浅色带是容差 ±0.10", "bars are intervals; the pale band is the ±0.10 margin", "end", 11.5))
    s.append(tx(pl - 44, 14, "读不出的比例 u", "Unreadable share u", "start", 11.5))
    return (f'<svg style="min-width:660px" viewBox="0 0 {W} {H}" role="img" data-al-zh="封存检验里 PD-L1 的响应在五环上读不出的比例。{"；".join(dz)}。" data-al-en="Unreadable share of the PD-L1 response at the five stages of the sealed tests. {"; ".join(de)}.">' + "".join(s) + "</svg>")
