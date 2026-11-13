#!/usr/bin/env python3
"""Build the site: docs/index.html, docs/anm-framework.html and docs/anm-loop.html.

    python scripts/build_site.py

Every number on the three pages is read here from a record committed under anm_diag/runs/. The figures are plain SVG
written by scripts/site_charts.py; the two explanatory figures and the formal definitions of the framework page are kept
as fragments in docs/_src/. Each page carries English and Chinese text and shows one of them (English first).
The pages of October 2026 (the registered decision-layer tests) are kept as docs/v3*.html.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import site_charts as SC  # noqa: E402

RUNS, DOCS, SRC = ROOT / "anm_diag" / "runs", ROOT / "docs", ROOT / "docs" / "_src"
GH = "https://github.com/danielchen26/teddy_mm"
UPDATED = ("2026 年 11 月 12 日", "12 November 2026")
U = '<span class="m">u</span>'
DS = [("frangieh2021", "黑色素瘤细胞（IFN-γ）", "Melanoma cells (IFN-γ)", "黑色素瘤", "melanoma"), ("papalexi2021", "THP-1 细胞", "THP-1 cells", "THP-1", "THP-1")]
GENES = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}
RO = {"PDL1": "PD-L1", "PDL2": "PD-L2", "HLAE": "HLA-E", "MHC1": "MHC-I"}
ST = {"located in the head": ("no", "在头里丢", "lost in the head"), "no loss in the head": ("ok", "头不再丢", "no loss in the head"), "not resolved": ("un", "没定清", "not resolved"), "no reading": ("un", "无读数", "no reading")}


def J(*p):
    return json.load(open(RUNS.joinpath(*p)))


def bi(zh, en):
    return f'<span class="zh">{zh}</span><span class="en">{en}</span>'


def tx(x, y, zh, en, anchor="start", size=12.5, cls="", weight=None):
    a = f' x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}"' + (f' font-weight="{weight}"' if weight else "")
    return f'<text class="zh {cls}"{a}>{zh}</text><text class="en {cls}"{a}>{en}</text>'


def gh(path, label=None):
    return f'<a href="{GH}/tree/main/{path}"><span class="mono">{label or path}</span></a>'


# ---------------------------------------------------------------- records
SEAL = {"frangieh2021": J("frangieh2021_headfix", "anm_record.json"), "papalexi2021": J("papalexi2021_headfix", "anm_record.json")}
L, LC, LF = J("loop_20261112", "results", "loop.json"), J("loop_20261112", "results", "cost_bone_marrow.json"), J("loop_20261112", "results", "fits.json")
G, GB = J("loop_guard_20261112", "results", "loop_guard.json"), J("loop_guard_20261112", "results", "bone_marrow.json")
HAVE3 = (RUNS / "repair_20261112" / "results" / "repair.json").exists()
R = J("repair_20261112", "results", "repair.json") if HAVE3 else None
A = {**L["audits"], **G["audits"], **(R["audits"] if HAVE3 else {})}
BASE = LC["as_trained"]["mean_pearson"]
ARM_OF_ROUND = ["as_trained", "mlp_only", "mlp_guarded"] + (["joint_s0"] if HAVE3 else [])
PASSED3 = bool(HAVE3 and R["one_head_for_both"])        # round 3 by its own rule
BEST = "joint_s0" if PASSED3 else "mlp_guarded"         # the last round that passed its rule
FE = J("repair_20261112", "results", "full_test_eval.json") if HAVE3 else None
TC = J("repair_20261112", "results", "torch_check.json") if HAVE3 else None
OS = J("repair_20261112", "results", "own_shift.json") if HAVE3 else None
assert all(SEAL[d]["seal_verified"] for d in SEAL) and all(L["rule"][d]["outcome"] == "closed at the located step" for d, *_ in DS) and all(G["rule"][d]["outcome"] == "closed within the guard" for d, *_ in DS)


def po(ds, arm, ro="PDL1"):
    return A[f"pooled_{ds}_{arm}"]["readouts"][ro]


def u_out(ds, arm, ro="PDL1"):
    return po(ds, arm, ro)["stages"]["HEAD_OUT"]["u"]


def bm_range(arm):
    """Lowest and highest bone-marrow test score of the heads of one arm, over both data sets."""
    if arm == "as_trained":
        return BASE, BASE
    if arm in ("mlp_only", "decoder_only"):
        v = [x["mean_pearson"] for k, x in LC["heads"].items() if k.endswith(arm)]
    elif arm == "mlp_guarded":
        v = [x["mean_pearson"] for x in GB["heads"].values()]
    else:
        v = R["bone_marrow"]["by_seed"]["0"]["heads"]
    return min(v), max(v)


def bm_ds(ds, arm):
    if arm in ("mlp_only", "decoder_only"):
        v = [x["mean_pearson"] for k, x in LC["heads"].items() if k.startswith(ds) and k.endswith(arm)]
    else:
        v = [x["mean_pearson"] for k, x in GB["heads"].items() if k.startswith(ds)]
    return min(v), max(v)


def cnt(arm, status):
    suf = "_joint_s0" if arm == "joint_s0" else "_" + arm
    return sum(1 for ds, *_ in DS for g in GENES[ds] if A[f"target_{ds}_{g}{suf}"]["readouts"]["PDL1"]["status"] == status)


def own(ds, arm):
    v = [x["own_output_shift_sd"] / x["measured_shift_sd"] for x in LF["side"][f"{ds}|{arm}"].values()]
    return min(v), max(v)


def stat(status):
    c, zh, en = ST[status]
    return f'<span class="st {c}">{bi(zh, en)}</span>'


def sealed(ds, ro, stage):
    return SEAL[ds]["readouts"][ro]["states"][stage]


def f2(v):
    return f"{v:.2f}".replace("-", "−")


def rng3(lo, hi):
    return f"{lo:.3f}" if abs(hi - lo) < 5e-4 else f"{lo:.3f}–{hi:.3f}"


# ---------------------------------------------------------------- shell
CSS = """
:root{ --bg:#fbfbf9; --panel:#f1f3f2; --line:#d3d9d6; --fg:#14212b; --muted:#55626b; --accent:#0f6b6e; --accent-tint:#dcecea;
  --s1:#2a78d6; --s2:#c8641f; --s3:#7a4fd0; --good:#12873a; --warn:#a86e00; --bad:#c93535;
  --sans:"IBM Plex Sans","PingFang SC","Hiragino Sans GB","Noto Sans SC","Microsoft YaHei",system-ui,sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,"SF Mono",Menlo,Consolas,monospace; --math:"STIX Two Text","Times New Roman",serif; color-scheme:light }
@media (prefers-color-scheme: dark){ :root{ --bg:#121618; --panel:#1b2125; --line:#2f393e; --fg:#e8ecea; --muted:#a6b1b6; --accent:#5fbdb8; --accent-tint:#17302f;
  --s1:#5a9bea; --s2:#e58a4a; --s3:#a98bf0; --good:#3fbf6b; --warn:#f0b33a; --bad:#ee7b7b; color-scheme:dark } }
*{ box-sizing:border-box }
html{ scroll-padding-top:4.5rem }
body{ margin:0; background:var(--bg); color:var(--fg); font-family:var(--sans); font-size:16px; line-height:1.65 }
:root[data-ui="zh"] body{ line-height:1.75 }
:root:not([data-ui="zh"]) .zh, :root[data-ui="zh"] .en{ display:none !important }
a{ color:var(--accent) } a:focus-visible, button:focus-visible, summary:focus-visible{ outline:2px solid var(--accent); outline-offset:2px }
header.site{ position:sticky; top:0; z-index:20; background:var(--bg); border-bottom:1px solid var(--line) }
header.site .in{ max-width:62rem; margin:0 auto; padding:.6rem max(16px, env(safe-area-inset-right, 0px)) .6rem max(16px, env(safe-area-inset-left, 0px)); display:flex; flex-wrap:wrap; gap:.4rem 1.2rem; align-items:center }
.brand{ font-weight:600; color:var(--fg); text-decoration:none; margin-right:auto; white-space:nowrap }
nav.pages{ display:flex; flex-wrap:wrap; gap:.2rem 1rem; font-size:.94rem }
nav.pages a{ color:var(--fg); text-decoration:none; padding:.15rem 0; border-bottom:2px solid transparent }
nav.pages a[aria-current="page"]{ border-color:var(--accent); font-weight:600 } nav.pages a:hover{ border-color:var(--line) }
.lang{ display:flex; gap:.3rem }
.lang button{ font:inherit; font-size:.86rem; line-height:1.3; padding:.3rem .6rem; border:1px solid var(--line); background:var(--bg); color:var(--fg); cursor:pointer; border-radius:3px }
.lang button[aria-pressed="true"]{ background:var(--accent-tint); border-color:var(--accent); font-weight:600 }
main{ max-width:62rem; margin:0 auto; padding:2.2rem max(16px, env(safe-area-inset-right, 0px)) 3rem max(16px, env(safe-area-inset-left, 0px)); display:flex; flex-direction:column; gap:2.75rem }
section{ display:flex; flex-direction:column; gap:1rem; min-width:0; scroll-margin-top:4.5rem }
.col{ max-width:46rem }
h1{ font-size:2.05rem; line-height:1.2; font-weight:600; margin:0; text-wrap:balance; max-width:50rem }
h2{ font-size:1.32rem; line-height:1.3; font-weight:600; margin:0; text-wrap:balance }
h3{ font-size:1.02rem; font-weight:600; margin:.4rem 0 0 }
p{ margin:0 }
.eyebrow{ font-family:var(--mono); font-size:.8rem; color:var(--muted) }
.lede{ font-size:1.1rem }
.src{ color:var(--muted); font-size:.9rem }
code,.mono{ font-family:var(--mono); font-size:.88em }
.m{ font-family:var(--math); font-style:italic; font-size:1.08em; white-space:nowrap }
.sh{ display:flex; gap:.7rem; align-items:baseline }
.sn{ font-family:var(--mono); font-size:.85rem; font-weight:500; color:var(--accent); border:1px solid var(--accent); min-width:1.55rem; text-align:center; border-radius:3px; flex:none; line-height:1.5 }
figure{ margin:0; display:flex; flex-direction:column; gap:.6rem; min-width:0 }
figcaption{ color:var(--muted); font-size:.9rem; max-width:46rem }
.scroll{ overflow-x:auto; border:1px solid var(--line); background:var(--panel) }
.scroll svg{ display:block; width:100%; height:auto }
svg text{ font-family:var(--sans); fill:currentColor }
svg .t-m{ font-family:var(--mono) } svg .t-mu{ fill:var(--muted) }
svg .bx{ fill:var(--bg); stroke:currentColor; stroke-width:1.2 }
svg .sys{ fill:none; stroke:var(--muted); stroke-width:1; stroke-dasharray:5 4 }
svg .ro{ fill:var(--accent-tint); stroke:var(--accent); stroke-width:1.4 }
svg .ar{ fill:none; stroke:currentColor; stroke-width:1.3 } svg .ar-d{ fill:none; stroke:var(--accent); stroke-width:1.5; stroke-dasharray:6 4 }
svg .tick{ stroke:var(--muted); stroke-width:1; stroke-dasharray:2 3 } svg .band{ fill:var(--accent-tint) } svg .ax{ stroke:var(--muted); stroke-width:1 }
.tw{ overflow-x:auto }
table{ border-collapse:collapse; width:100%; font-size:.92rem; font-variant-numeric:tabular-nums }
th,td{ text-align:left; padding:.5rem .7rem; border-bottom:1px solid var(--line); vertical-align:top }
th{ font-weight:600; color:var(--muted); font-size:.82rem }
td.n{ white-space:nowrap; font-family:var(--mono); font-size:.84rem }
td.k{ font-weight:600; white-space:nowrap }
.st{ white-space:nowrap; font-weight:500 }
.st::before{ content:""; display:inline-block; width:.6em; height:.6em; margin-right:.4em; border-radius:50% }
.st.ok::before{ background:var(--good) } .st.no::before{ background:var(--bad); border-radius:0 } .st.un::before{ background:var(--warn); border-radius:0; transform:rotate(45deg) }
.note{ border-left:2px solid var(--accent); padding:.2rem 0 .2rem 1rem; max-width:46rem }
.glance{ display:grid; grid-template-columns:repeat(4, minmax(0, 1fr)); gap:1px; background:var(--line); border:1px solid var(--line) }
.gcard{ background:var(--bg); padding:.9rem 1rem 1rem; display:flex; flex-direction:column; gap:.4rem; min-width:0 }
.gcard .k{ font-family:var(--mono); font-size:.76rem; color:var(--accent); letter-spacing:.04em }
.gcard .big{ font-size:1.16rem; font-weight:600; line-height:1.3; text-wrap:balance; font-variant-numeric:tabular-nums }
.gcard p{ font-size:.9rem; line-height:1.55; color:var(--muted) }
.gcard a{ font-size:.86rem; margin-top:auto; padding-top:.2rem }
ul.pts, ol.pts{ margin:0; padding-left:1.2rem; display:flex; flex-direction:column; gap:.45rem; max-width:46rem }
ul.pts b, ol.pts b, dl.roles b{ font-weight:600 }
dl.roles{ display:grid; grid-template-columns:minmax(7rem,11rem) 1fr; gap:.55rem 1.25rem; margin:0; max-width:46rem }
dl.roles dt{ font-weight:600 } dl.roles dd{ margin:0; min-width:0 }
ol.defs{ margin:0; padding-left:1.4rem; display:flex; flex-direction:column; gap:.65rem; max-width:46rem }
details.more{ max-width:46rem; border-top:1px solid var(--line); border-bottom:1px solid var(--line) }
details.more.wide{ max-width:none }
details.more > summary{ cursor:pointer; padding:.6rem 0; color:var(--accent); font-size:.94rem; font-weight:500 }
details.more > .db{ display:flex; flex-direction:column; gap:.9rem; padding:.2rem 0 1.1rem; min-width:0 }
a.next{ display:flex; flex-wrap:wrap; gap:.2rem 1rem; align-items:baseline; max-width:46rem; padding:.8rem 1rem; border:1px solid var(--accent); background:var(--accent-tint); color:var(--fg); text-decoration:none }
a.next .t{ flex:1 1 16rem; min-width:0; font-weight:500 } a.next .go{ color:var(--accent); font-weight:600 } a.next:hover .t{ text-decoration:underline }
.steps{ list-style:none; margin:0; padding:0; display:grid; grid-template-columns:repeat(auto-fit, minmax(13.5rem, 1fr)); gap:.6rem }
.steps li{ border:1px solid var(--line); background:var(--panel); padding:.65rem .8rem .75rem; min-width:0; display:flex; flex-direction:column; gap:.3rem; font-size:.92rem; line-height:1.6 }
.steps .h{ display:flex; gap:.5rem; align-items:baseline; font-weight:600; font-size:1rem } .steps .h .mono{ color:var(--accent) }
footer{ border-top:1px solid var(--line); color:var(--muted); font-size:.88rem }
footer .in{ max-width:62rem; margin:0 auto; padding:1.2rem max(16px, env(safe-area-inset-right, 0px)) 2rem max(16px, env(safe-area-inset-left, 0px)); display:flex; flex-direction:column; gap:.4rem }
@media (max-width:900px){ .glance{ grid-template-columns:1fr 1fr } }
@media (max-width:560px){ .glance{ grid-template-columns:1fr } h1{ font-size:1.6rem } dl.roles{ grid-template-columns:1fr; gap:.1rem 0 } dl.roles dd{ margin-bottom:.6rem } }
"""
JS = """
(function(){
  var root = document.documentElement, saved = null;
  try { saved = localStorage.getItem('teddy-anm-lang'); } catch (e) {}
  function set(l){
    root.setAttribute('data-ui', l); root.lang = l === 'zh' ? 'zh-Hans' : 'en';
    document.querySelectorAll('.lang button').forEach(function(b){ b.setAttribute('aria-pressed', String(b.dataset.l === l)); });
    document.querySelectorAll('[data-al-' + l + ']').forEach(function(el){ el.setAttribute('aria-label', el.getAttribute('data-al-' + l)); });
    try { localStorage.setItem('teddy-anm-lang', l); } catch (e) {}
  }
  document.querySelectorAll('.lang button').forEach(function(b){ b.addEventListener('click', function(){ set(b.dataset.l); }); });
  set(saved === 'zh' ? 'zh' : 'en');
})();
"""
PAGES = [("index.html", "总览", "Overview"), ("anm-framework.html", "审计怎么做", "How the audit works"), ("anm-loop.html", "闭环", "The loop")]


def page(fname, title, desc, body):
    nav = "".join(f'<a href="{f}"' + (' aria-current="page"' if f == fname else "") + f'>{bi(zh, en)}</a>' for f, zh, en in PAGES)
    html_ = f"""<!doctype html>
<html lang="en" data-ui="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{title}</title>
<meta name="description" content="{desc}">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=STIX+Two+Text:ital@1&display=swap">
<style>{CSS}</style>
</head>
<body>
<header class="site"><div class="in">
  <a class="brand" href="index.html">TEDDY × ANM</a>
  <nav class="pages" aria-label="Pages">{nav}<a href="{GH}/blob/main/docs/ANM_DIAGNOSIS.md">{bi("诊断全文", "Full write-up")}</a><a href="{GH}">GitHub</a></nav>
  <div class="lang" role="group" aria-label="Language"><button type="button" data-l="en" aria-pressed="true">English</button><button type="button" data-l="zh" aria-pressed="false">中文</button></div>
</div></header>
<main>
{body}
</main>
<footer><div class="in">
  <p>{bi(f"这三页上的每个数字都由 <span class='mono'>scripts/build_site.py</span> 从仓库里提交的记录读出（<span class='mono'>anm_diag/runs/</span>）。更新于 {UPDATED[0]}。", f"Every number on these three pages is read by <span class='mono'>scripts/build_site.py</span> from records committed in the repository (<span class='mono'>anm_diag/runs/</span>). Updated {UPDATED[1]}.")}</p>
  <p>{bi("方法：Active Neural Matter（arXiv 2610.14549）。把模型的中间量当作候选状态，是对论文定义的套用。", "Method: Active Neural Matter (arXiv 2610.14549). Treating a model's intermediate quantities as candidate states is an application of the paper's definitions.")} · <a href="v3.html">{bi("2026 年 10 月的旧页面（已存档）", "The pages of October 2026 (archived)")}</a></p>
</div></footer>
<script>{JS}</script>
</body>
</html>
"""
    (DOCS / fname).write_text(html_, encoding="utf-8")
    print("wrote", fname, len(html_))


def sh(n, zh, en):
    return f'<div class="sh"><span class="sn">{n}</span><h2>{bi(zh, en)}</h2></div>'


# ---------------------------------------------------------------- shared blocks
def rounds_table():
    """The loop so far, one row per round."""
    m, t = "frangieh2021", "papalexi2021"
    left2 = sorted(RO[r] for r, v in G["audits"][f"pooled_{m}_mlp_guarded"]["readouts"].items() if v["status"] == "located in the head")
    rows = [
        (("0", "不改：原来的头", "none: the head as trained"), ("第一次审计", "the first audit"), "as_trained",
         ("PD-L1 的响应在头里丢，从嵌入到头的中间层（头的 MLP）开始", "the PD-L1 response is lost in the head, from the embedding to the head latent (the head's MLP)")),
        (("1", "只修 MLP，训 PD-L1", "repair the MLP on PD-L1"), ("审计指到的那一步", "the step the audit names"), "mlp_only",
         ("骨髓原任务掉分超过 0.02（由准确率的护栏发现）", "the original bone-marrow task falls by more than 0.02 (found by the accuracy guard)")),
        (("2", "同样的修法，按住骨髓上的输出", "the same, with the bone-marrow outputs held"), ("上一轮的代价", "the cost of round 1"), "mlp_guarded",
         (f"{'、'.join(left2)} 在黑色素瘤上仍在头里丢；在一种细胞上修的头，到另一种细胞上仍然丢", f"{', '.join(left2)} is still lost in the head on melanoma; a head repaired on one cell type still loses it on the other")),
    ]
    if HAVE3:
        s_zh, s_en = r3_left()
        rows.append((("3", "一个头，两种细胞和全部实测蛋白一起修", "one head: both cell types and every measured protein in one repair"), ("上一轮审计还指出的两处", "the two findings the round-2 audit still located"), "joint_s0", (s_zh, s_en)))
    body = []
    for (n, czh, cen), (wzh, wen), arm, (lzh, len_) in rows:
        lo, hi = bm_range(arm)
        body.append(f'<tr><td class="k">{n}</td><td>{bi(czh, cen)}<br><span class="src">{bi("依据：" + wzh, "prompted by " + wen)}</span></td><td class="n">{f2(u_out(m, arm))}<br>{stat(po(m, arm)["status"])}</td><td class="n">{f2(u_out(t, arm))}<br>{stat(po(t, arm)["status"])}</td><td class="n">{rng3(lo, hi)}</td><td>{bi(lzh, len_)}</td></tr>')
    return ('<div class="tw"><table><thead><tr><th>' + bi("轮", "Round") + '</th><th>' + bi("改了什么", "What was changed") + '</th><th>' + bi("PD-L1，黑色素瘤", "PD-L1, melanoma") + '</th><th>' + bi("PD-L1，THP-1", "PD-L1, THP-1")
            + '</th><th>' + bi("骨髓测试得分", "Bone-marrow test score") + '</th><th>' + bi("再审之后还剩什么", "What the next audit still found") + '</th></tr></thead><tbody>' + "".join(body) + '</tbody></table></div>'
            + '<p class="src col">' + bi(f"PD-L1 两列是头的输出上读不出的比例 {U}，在修复时留出的基因的敲除细胞上审计。骨髓得分是这一轮各个头的最低到最高；事先定的界限是比原来的头低不超过 0.02。", f"The PD-L1 columns give the unreadable share {U} at the head output, audited on knockout cells of genes held out of the repair. The bone-marrow score is the lowest to highest of that round's heads; the guard set in advance is at most 0.02 below the head as trained.") + '</p>')


def r3_left():
    """What the audit of round 3 found, per data set, in both languages."""
    zh, en = [], []
    for d in DS:
        ru = R["rule"][d[0]]; g = po(d[0], "joint_s0"); loc = [RO[x] for x in ru["located_in_the_head"]]
        if ru["outcome"] == "regressed":
            if g["status"] != "no loss in the head":
                zh.append(f"<b>{d[3]}：没有通过规则。</b>PD-L1 在输出上比 Z 高 {g['added_by_head']:.3f}，过了 0.10 的界限")
                en.append(f"<b>{d[4]}: did not pass its rule.</b> PD-L1 at the output is {g['added_by_head']:.3f} above Z, over the margin of 0.10")
            if not ru["guard_kept"]:
                zh.append(f"{d[3]}：骨髓得分掉出界限"); en.append(f"{d[4]}: the bone-marrow score falls beyond the guard")
        if loc:
            zh.append(f"{d[3]}：{'、'.join(loc)} 仍在头里丢"); en.append(f"{d[4]}: {', '.join(loc)} still lost in the head")
    return ("；".join(zh) or "头里不再有被指出的丢失"), ("; ".join(en) or "no loss is located in the head any more")


def loop_chart():
    return SC.chart(tx, A)


def cost_rows():
    rows = []
    for arm, zh, en in (("decoder_only", "对照：只修解码器", "control: decoder repaired"), ("mlp_only", "第 1 轮：只修 MLP", "round 1: MLP repaired")):
        for ds, _, _, dzh, den in DS:
            rows.append((arm, f"{zh}，{dzh}", f"{en}, {den}", *bm_ds(ds, arm)))
    for ds, _, _, dzh, den in DS:
        rows.append(("mlp_guarded", f"第 2 轮：按住骨髓输出，{dzh}", f"round 2: bone marrow held, {den}", *bm_ds(ds, "mlp_guarded")))
    if HAVE3:
        rows.append(("joint_s0", "第 3 轮：一个头，两种细胞", "round 3: one head, both cell types", *bm_range("joint_s0")))
    return rows


# ---------------------------------------------------------------- index
def index():
    m, t = "frangieh2021", "papalexi2021"
    z, lat, out = sealed(m, "PDL1", "Z"), sealed(m, "PDL1", "HEAD_LAT"), sealed(m, "PDL1", "HEAD_OUT")
    assert z["verdict"] == "support" and lat["verdict"] == "reject" and SEAL[m]["readouts"]["PDL1"]["root_cause"]["cause"] == "loss_at"
    last = BEST; lo, hi = bm_range(last)
    cards = [
        (("定位", "LOCATE"), f'{U} {f2(z["u"])} → {f2(lat["u"])} → {f2(out["u"])}', ("封存、只跑一次。IFN-γ 通路敲除后 PD-L1 的响应，在 TEDDY 的嵌入里读得出，到头的中间层开始读不出，到输出丢了近一半。", "Sealed and run once. The PD-L1 response to IFN-γ-pathway knockouts is readable from TEDDY's embedding, starts to go at the head latent and is nearly half gone at the head's output."), ("#where", "看第 1 节", "See section 1")),
        (("拆开", "SEPARATE"), bi(f"{cnt('as_trained', 'located in the head')} / 10 个靶基因", f"{cnt('as_trained', 'located in the head')} of 10 target genes"), ("合并的源拆成单个靶基因后，PD-L1 的响应在每一个上都被指出在头里丢；TEDDY 自己的三环没有一环被否定。", "With the pooled source separated into single target genes, the PD-L1 response is located in the head for each one; none of TEDDY's own three stages is rejected."), ("anm-framework.html#events", "看类型化事件", "See typed events")),
        (("修", "REPAIR"), f'{U} {f2(u_out(m, "as_trained"))} → {f2(u_out(m, last))}', (f"只修审计指到的那一步（头的 MLP），在没见过的基因上丢失消失；只修解码器，仍是 {f2(u_out(m, 'decoder_only'))}。THP-1 上 {f2(u_out(t, 'as_trained'))} → {f2(u_out(t, last))}。", f"Repairing only the step the audit names (the head's MLP) removes the loss on unseen genes; repairing the decoder leaves {f2(u_out(m, 'decoder_only'))}. On THP-1: {f2(u_out(t, 'as_trained'))} → {f2(u_out(t, last))}."), ("anm-loop.html", "看闭环", "See the loop")),
        (("守住原任务", "KEEP THE OLD TASK"), f"{BASE:.3f} → {rng3(lo, hi)}", ("修的时候把头在骨髓细胞上的输出按住，原任务的得分留在事先定的 0.02 以内。不按住时掉到 0.50 到 0.57。", "With the head's outputs held on bone-marrow cells during the repair, the score on its original task stays within the 0.02 set in advance. Without that it falls to 0.50–0.57."), ("anm-loop.html#cost", "看代价", "See the cost")),
    ]
    glance = '<div class="glance">' + "".join(f'<div class="gcard"><span class="k">{bi(*k)}</span><span class="big">{big}</span><p>{bi(*p)}</p><a href="{a[0]}">{bi(a[1], a[2])}</a></div>' for k, big, p, a in cards) + '</div>'
    other = []
    for ds, zh, en, *_ in DS:
        for ro, v in SEAL[ds]["readouts"].items():
            st = v["states"]; rc = v["root_cause"]
            other.append(f'<tr><td>{bi(zh, en)}</td><td class="k">{RO[ro]}</td><td class="n">{f2(v["delta"])}</td>' + "".join(f'<td class="n">{f2(st[k]["u"])} <span class="src">{st[k]["verdict"]}</span></td>' for k in ("RAW", "TOKENS", "Z", "HEAD_LAT", "HEAD_OUT")) + f'<td><span class="mono">{rc["cause"]}</span>' + (f' {rc["from"]} → {rc["to"]}' if "from" in rc else "") + '</td></tr>')
    repo = [
        ("frangieh2021_headfix", "papalexi2021_headfix", ("封存，只跑一次", "sealed, run once"), ("丢失在哪一步：五个蛋白读出，两套数据", "where the loss is: five protein readouts, two data sets")),
        ("frangieh2021_headsplit", "papalexi2021_headsplit", ("探索性，封存", "exploratory, sealed"), ("把头拆成七环：丢失在后面的层和解码器，不在输入归一化和第一层", "the head split into seven stages: the loss is in its later layers and decoder, not in the input normalization or the first layer")),
        ("split_20261111", "exposure_20261111", ("探索性", "exploratory"), ("合并的源按靶基因和批次拆开，源写成类型化事件；对 IRF1 的一个解释做了直接检验，按规则没有回答", "the pooled source separated by target gene and batch, with the source written as typed events; a direct test of one explanation for IRF1, not answered by its rule")),
        ("loop_20261112", "loop_guard_20261112", ("探索性，规则先提交", "exploratory, rule committed first"), ("闭环第 1、2 轮：修 MLP 有效，修解码器无效；按住骨髓输出后在界限内", "rounds 1 and 2 of the loop: repairing the MLP works and the decoder does not; within the guard once the bone-marrow outputs are held")),
    ] + ([("repair_20261112", None, ("探索性，规则先提交", "exploratory, rule committed first"), ("闭环第 3 轮：一个头同时修两种细胞和全部实测蛋白，另在两次重新划分上重复" + ("" if PASSED3 else "；按它自己的规则没有通过"), "round 3 of the loop: one head repaired for both cell types and every measured protein, repeated on two fresh splits" + ("" if PASSED3 else "; it did not pass its own rule")))] if HAVE3 else [])
    repo_rows = "".join(f'<tr><td>{gh("anm_diag/runs/" + a, a)}' + (f'<br>{gh("anm_diag/runs/" + b, b)}' if b else "") + f'</td><td>{bi(*g)}</td><td>{bi(*w)}</td></tr>' for a, b, g, w in repo)
    body = f"""<section>
  <p class="eyebrow">TEDDY-G 70M · teddy_mm phase-1 head · {bi("ANM 审计", "ANM audit")} · {bi(UPDATED[0], UPDATED[1])}</p>
  <h1>{bi("敲除引起的响应，TEDDY 留住了；接在它后面的蛋白头丢了近一半。", "TEDDY keeps the response to a knockout. The protein head built on it loses nearly half.")}</h1>
  <p class="lede col">{bi("ANM 审计沿着模型自己的计算，一环一环地跟踪一个实测的响应，指出它从哪一步开始读不出来。在两套 Perturb-CITE-seq 数据上，这一步是 teddy_mm 蛋白头的第一段网络。只修这一步，丢失就在修复时没见过的敲除上消失；修别处没有用。", "An ANM audit follows one measured response through each stage of a model's own computation and names the step where it stops being readable. On two Perturb-CITE-seq data sets that step is the first network of the teddy_mm protein head. Repairing that step removes the loss on knockouts the repair never saw; repairing elsewhere does not.")}</p>
  {glance}
  <p class="src col">{bi("现在的状态：修复是在头的微调副本上验证的，<span class='mono'>scripts/04_train.py</span> 训出来的头还是没修的那个。" + ("修法已经写成脚本 <span class='mono'>scripts/07_repair_head.py</span>。" + ("" if PASSED3 else "通过了规则的是每种细胞一个头的修法（第 2 轮）；想修成一个头管两种细胞的第 3 轮没有通过。") + "见“闭环”一页。" if HAVE3 else ""), "Where this stands: the repair is shown on fine-tuned copies of the head; the head that <span class='mono'>scripts/04_train.py</span> trains is still the unrepaired one." + (" The repair exists as a script, <span class='mono'>scripts/07_repair_head.py</span>." + ("" if PASSED3 else " The form that passed its rule gives one head per cell type (round 2); round 3, which tried one head for both, did not pass.") + " See the page on the loop." if HAVE3 else ""))}</p>
</section>

<section id="where">
  {sh(1, "响应在哪一步丢", "Where the response is lost")}
  <p class="col">{bi(f"系统是 TEDDY-G 70M 加上 teddy_mm 第一阶段从 RNA 预测 134 个蛋白的头。源是 IFN-γ 通路基因的 CRISPR 敲除，带非靶向对照。读出是 CITE-seq 实测的 PD-L1 蛋白，模型从没读过它。每一环问同一个问题：只凭这一环，实测的响应还有多少解释不了（{U}；容差 0.10）。", f"The system is TEDDY-G 70M followed by the teddy_mm phase-1 head, which predicts 134 proteins from RNA. The source is CRISPR knockout of IFN-γ-pathway genes, with non-targeting controls. The readout is PD-L1 protein measured by CITE-seq, which the model never reads. Each stage is asked the same question: from this stage alone, how much of the measured response is left unexplained ({U}; margin 0.10).")}</p>
  <figure>
    <div class="scroll">{SC.chain(tx, [(d[1], d[2], SEAL[d[0]]) for d in DS])}</div>
    <figcaption>{bi("封存的检验，只跑一次。虚线左边是 TEDDY 自己的三环，右边是头的两环。", "The sealed tests, run once. Left of the dashed line are TEDDY's own three stages, right of it the two stages of the head.")}</figcaption>
  </figure>
  <p class="note">{bi(f"黑色素瘤：嵌入得到支持（{f2(z['u'])}，区间 {f2(z['ci'][0])} 到 {f2(z['ci'][1])}），头的中间层被否定（{f2(lat['u'])}），输出 {f2(out['u'])}。工具给出的定位是嵌入到头的中间层这一步。THP-1 上是同一步，但输入端本身没定清（5 个靶基因），所以算相对丢失。", f"Melanoma: the embedding is supported ({f2(z['u'])}, interval {f2(z['ci'][0])} to {f2(z['ci'][1])}), the head latent is rejected ({f2(lat['u'])}) and the output is at {f2(out['u'])}. The tool locates the loss at the step from the embedding to the head latent. On THP-1 it is the same step, as a relative loss, because the input view itself is unresolved with five target genes.")}</p>
  <details class="more wide"><summary>{bi("五个蛋白读出的全部数值", "All five protein readouts")}</summary><div class="db"><div class="tw"><table><thead><tr><th>{bi("数据", "Data")}</th><th>{bi("读出", "Readout")}</th><th>{bi("实测变化", "Measured change")}</th><th>RAW</th><th>TOKENS</th><th>Z</th><th>{bi("头的中间层", "Head latent")}</th><th>{bi("头的输出", "Head output")}</th><th>{bi("工具给出的根因", "Root cause printed by the tool")}</th></tr></thead><tbody>{"".join(other)}</tbody></table></div>
  <p class="src col">{bi("TEDDY 的三环在五个读出上都没有被否定。MHC-I 和 HLA-E 的实测变化很小，区间很宽。", "None of TEDDY's three stages is rejected for any of the five readouts. MHC-I and HLA-E change little and their intervals are wide.")}</p></div></details>
</section>

<section id="loop">
  {sh(2, "闭环：定位，修，再审", "The loop: locate, repair, audit again")}
  <p class="col">{bi("定位只是观察。要知道它对不对，就照它去改，再用同一把尺子量一遍。每一轮的设计和判定规则都在跑之前提交进仓库。", "A location is an observation. To learn whether it is right, change the model where it points and measure again with the same instrument. The design and the decision rule of every round were committed to the repository before the round was run.")}</p>
  {rounds_table()}
  <a class="next" href="anm-loop.html"><span class="t">{bi("每一轮的图、对照、代价和限制", "Figures, the control, the cost and the limits of each round")}</span><span class="go" aria-hidden="true">→</span></a>
</section>

<section id="repo">
  {sh(3, "仓库里有什么", "What the repository holds")}
  <div class="tw"><table><thead><tr><th>{bi("运行记录", "Run records")}</th><th>{bi("证据等级", "Evidence grade")}</th><th>{bi("内容", "What it is")}</th></tr></thead><tbody>{repo_rows}</tbody></table></div>
  <ul class="pts">
    <li>{bi(f"审计的代码和测试在 {gh('anm_diag')}；写法和全部修正记录在 <a href='{GH}/blob/main/docs/ANM_DIAGNOSIS.md'>诊断全文</a>。", f"The audit's code and tests are in {gh('anm_diag')}; the write-up with every correction is the <a href='{GH}/blob/main/docs/ANM_DIAGNOSIS.md'>full diagnosis</a>.")}</li>
    <li>{bi("封存的第 4 步规则（三个种子、带对照的重训）在仓库里，但按它做的重训没有记录在仓库里，这里不报告它的任何数字。", "The sealed step-4 rule (three seeds, a retrain with a control arm) is in the repository; no retrain made under it is recorded there, and no number from one is reported here.")}</li>
  </ul>
</section>

<section id="claims">
  {sh(4, "能说什么，不能说什么", "What can and cannot be claimed")}
  <ul class="pts">
    <li>{bi("<b>能说（封存的检验）。</b>在 Frangieh 2021 的 IFN-γ 通路敲除上，TEDDY-G 70M 的嵌入保留了解释 PD-L1 蛋白响应所需的信息（容差 0.10 以内）；五个蛋白读出上 TEDDY 没有一环被否定；teddy_mm 的头在嵌入和它的中间层之间丢掉了 PD-L1 的响应。", "<b>Can claim (sealed tests).</b> For the IFN-γ-pathway knockouts of Frangieh 2021, TEDDY-G 70M's embedding keeps the information that explains the PD-L1 protein response, within a 0.10 margin; no TEDDY stage is rejected for any of five protein readouts; the teddy_mm head loses the PD-L1 response between the embedding and its latent layer.")}</li>
    <li>{bi("<b>探索性。</b>把头拆开、把源拆开、类型化事件，以及闭环的各轮：修头的 MLP 消除了留出基因上的丢失，修解码器没有；按住骨髓输出后代价在 0.02 以内。", "<b>Exploratory.</b> The head split, the source separated, typed events, and the rounds of the loop: repairing the head's MLP removes the loss on held-out genes and repairing its decoder does not; with the bone-marrow outputs held the cost stays within 0.02.")}</li>
    <li>{bi("<b>不能说。</b>训练脚本训出的头已经修好（修好的是微调副本，一个种子）；修复能带到第三种细胞；头为什么丢；TEDDY 对一般的扰动响应都够用；被否定的环“坏了”（否定只是在声明的读出模型下指出一步）。", "<b>Cannot claim.</b> That the head the training script produces is repaired (the repaired heads are fine-tuned copies, one seed); that the repair carries to a third cell type; why the head loses the response; that TEDDY is sufficient for perturbation responses in general; that a rejected stage is “broken” (a rejection names a step, under the declared readout classes).")}</li>
  </ul>
</section>

<section id="earlier">
  {sh(5, "这个站点以前的内容", "Earlier work on this site")}
  <p class="col">{bi("2026 年 10 月的页面问的是另一个问题：把 ANM 当作 TEDDY 之上的决策层，有没有增益。注册的 v3 检验的结论是没有：ANM 的引擎算出的正好是一条声明的规则。那些页面原样存档：", "The pages of October 2026 asked a different question: does ANM add anything as a decision layer on top of TEDDY. The registered v3 tests found that it does not: ANM's engine computes exactly a declared rule. Those pages are kept as they were:")} <a href="v3.html">{bi("注册的决策层检验", "registered decision-layer tests")}</a> · <a href="v3-anm-framework.html">{bi("当时的框架页", "the framework page of that time")}</a> · <a href="v3-anm-loop.html">{bi("决策回路的首次记录", "the first decision-loop record")}</a>.</p>
</section>"""
    page("index.html", "TEDDY × ANM: where a knockout response is lost, and the loop that repairs it",
         "An ANM audit of TEDDY-G 70M and the teddy_mm protein head on two Perturb-CITE-seq data sets: TEDDY's stages keep the PD-L1 response to IFN-γ-pathway knockouts, the head loses nearly half of it at its first network, and repairing that step removes the loss on held-out genes.", body)


# ---------------------------------------------------------------- framework
CHECKS = [
    ("S1", "必须有带对照的受控源扰动", "a controlled source perturbation with controls is required", "没有源的“决策层”分析", "a decision layer with no source"),
    ("R1", "读出必须是测量，不能是模型的节点", "a readout must be a measurement, not a node", "拿模型自己的预测给模型打分", "scoring the model against its own prediction"),
    ("R2", "读出的测量不能流进候选状态；如果它也是系统的输入，要声明共享方式", "a readout's measurement must not flow into a candidate state; if it is also a system input, the sharing must be declared", "把实测蛋白喂进去，再和实测蛋白比", "feeding measured protein in and scoring against measured protein"),
    ("C1", "链上每一环都由上一环算出", "each chain stage must be computed from the previous one", "比较不在同一条计算路径上的量", "comparing quantities that are not on one computation path"),
    ("C2", "参照环之后的每一环都是系统自己的运算", "stages after the reference must be the system's own operators", "路径上混进近似的分词器", "an approximate tokenizer on the path"),
    ("C3", "分析者不加压缩（PCA、k-means 等）", "no analyst compression (PCA, k-means, …)", "没声明的降维造成的假象", "artefacts of an undeclared dimension reduction"),
    ("C4", "不遮输入，不改权重，不替换状态", "no input masking, weight edits or state replacement", "把对系统的干预当成充分性检验", "taking an intervention on the system for a sufficiency test"),
    ("E1", "检验量必须是连续的条件检验", "the statistic must be a continuous conditional test", "粗糙的聚类划分", "a coarse partition into clusters"),
    ("K1", "容差必须声明；处理的靶点少于 10 个时给出警告", "the tolerance must be declared; a warning below 10 treated targets", "事后挑容差；把宽区间当结论", "choosing the tolerance afterwards; reading wide intervals as results"),
    ("seal", "运行前核对封存的摘要、数据哈希和代码哈希", "the run verifies the seal digest, the data hashes and the code hash", "封存之后改设计、数据或代码", "changing design, data or code after sealing"),
]
CAUSES = [
    ("no_effect", "读出没有可解释的响应", "the readout has no response to explain"),
    ("control_failed", "阳性对照没读出来，这次检验本身不可信", "the positive control is not recovered, so the test itself is not trusted"),
    ("loss_at", "参照环得到支持，后面某一环被否定：丢失在两者之间", "the reference stage is supported and a later one is rejected: the loss lies between them"),
    ("relative_loss_at", "参照环没在容差内，但有一步明确多丢了", "the reference stage is not within tolerance, but one step adds a resolved loss"),
    ("reference_unresolved", "参照环没定清，也没有哪一步明确多丢", "the reference stage is unresolved and no step adds a resolved loss"),
    ("unresolved", "区间跨过容差，下不了结论", "the intervals cross the margin; no conclusion"),
    ("no_deficit", "每一环都在容差内", "every stage is within the margin"),
]


def loop_diagram():
    steps = [("审计", "Audit", ("跟踪实测的响应", "一环一环"), ("follow the measured", "response, stage by stage")), ("定位", "Locate", ("工具指出", "丢在哪一步"), ("the tool names the", "step that loses it")),
             ("开方", "Prescribe", ("改那一步；再写下", "不该有用的改法"), ("change that step; name", "a change that should fail")), ("修", "Repair", ("留出一部分基因", "和对照细胞"), ("hold out genes", "and control cells")),
             ("再审", "Audit again", ("只用修复时", "没见过的细胞"), ("only on cells the", "repair never saw")), ("护栏", "Guard", ("原任务", "不能掉"), ("the original task", "must hold"))]
    W, H, bw, bh, gap, x0, y0 = 960, 166, 138, 80, 22, 12, 34
    s = [f'<svg style="min-width:820px" viewBox="0 0 {W} {H}" role="img" data-al-zh="闭环的六步：审计，定位，开方，修，再审，护栏；然后回到审计。" data-al-en="The six steps of the loop: audit, locate, prescribe, repair, audit again, guard; then back to the audit."><defs><marker id="la" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0 L10 5 L0 10 z" fill="currentColor"/></marker></defs>']
    for i, (zh, en, szh, sen) in enumerate(steps):
        x = x0 + i * (bw + gap)
        s.append(f'<rect class="{"ro" if i in (1, 4) else "bx"}" x="{x}" y="{y0}" width="{bw}" height="{bh}" rx="3"/>' + tx(x + bw / 2, y0 + 26, f"{i + 1} · {zh}", f"{i + 1} · {en}", "middle", 13, weight=600) + tx(x + bw / 2, y0 + 47, szh[0], sen[0], "middle", 10.5, cls="t-mu") + tx(x + bw / 2, y0 + 62, szh[1], sen[1], "middle", 10.5, cls="t-mu"))
        if i < 5:
            s.append(f'<path class="ar" d="M{x + bw} {y0 + bh / 2} H{x + bw + gap - 2}" marker-end="url(#la)"/>')
    xe = x0 + 5 * (bw + gap) + bw / 2; xs = x0 + bw / 2
    s.append(f'<path class="ar-d" d="M{xe} {y0 + bh} V{y0 + bh + 26} H{xs} V{y0 + bh + 3}" marker-end="url(#la)"/>' + tx((xs + xe) / 2, y0 + bh + 44, "再审还指出的丢失，进下一轮", "what the next audit still locates goes into the next round", "middle", 11.5))
    s.append(tx(x0, 18, "每一轮的设计和判定规则先提交，再运行", "The design and the decision rule of each round are committed before it runs", "start", 11.5, cls="t-mu"))
    return "".join(s) + "</svg>"


def framework():
    frag = lambda n: (SRC / n).read_text(encoding="utf-8")
    checks = "".join(f'<tr><td class="k mono">{c}</td><td>{bi(rz, re_)}</td><td>{bi(pz, pe)}</td></tr>' for c, rz, re_, pz, pe in CHECKS)
    causes = "".join(f'<tr><td class="mono">{c}</td><td>{bi(z, e)}</td></tr>' for c, z, e in CAUSES)
    body = f"""<section>
  <p class="eyebrow">ANM · arXiv 2610.14549 · <span class="mono">anm_diag</span></p>
  <h1>{bi("审计怎么做", "How the audit works")}</h1>
  <p class="lede col">{bi("给定一个模型和一次受控的扰动，审计回答三个问题：实测的响应在模型内部哪一步还读得出来，从哪一步开始读不出来，改了模型以后有没有变好。", "Given a model and one controlled perturbation, the audit answers three questions: at which internal step the measured response can still be read, where it stops being readable, and whether changing the model made it better.")}</p>
</section>

<section id="roles">
  {sh(1, "先声明五样东西", "Five things are declared first")}
  <p class="col">{bi("审计前先声明五样东西，连同数据和代码的指纹（sha256）一起封存，再只跑一次。", "Five things are declared before the audit and sealed together with fingerprints (sha256) of the data and code. The run is then made once.")}</p>
  <figure><div class="scroll">{frag("roles.svg").replace('<svg ', '<svg style="min-width:720px" ', 1)}</div>
  <figcaption>{bi("扰动同时改变系统的输入和系统外独立测量的读出。审计沿着系统自己的计算，一环一环地问这个响应还在不在。", "The perturbation changes both the system's input and a readout measured independently, outside the system. The audit follows the system's own computation and asks at each stage whether the response is still there.")}</figcaption></figure>
  {frag("roles_list.html")}
</section>

<section id="verdict">
  {sh(2, "检验量和判定", "The statistic and the verdict")}
  <p class="col">{bi(f"对每一环 <span class='m'>z</span>，用声明的读出模型（线性的 ridge 和树模型 HGB，按靶点交叉拟合）拟合读出和处理，再看扣掉这一环知道的一切之后，处理还能解释读出多少：<span class='m'>θ(z) = E[(O − m(z))(A − e(z))] / E[(A − e(z))²]</span>，{U} <span class='m'>= θ(z) / Δ</span>。{U} 是读不出的那部分占实测响应 <span class='m'>Δ</span> 的比例：0 是全读得出，1 是一点读不出。", f"For each stage <span class='m'>z</span>, the declared readout models (linear ridge and tree-based HGB, cross-fitted by target) fit the readout and the treatment; what counts is how much of the readout the treatment still explains once everything the stage knows has been used: <span class='m'>θ(z) = E[(O − m(z))(A − e(z))] / E[(A − e(z))²]</span>, {U} <span class='m'>= θ(z) / Δ</span>. {U} is the unreadable part as a share of the measured response <span class='m'>Δ</span>: 0 means all of it is readable, 1 none.")}</p>
  <figure><div class="scroll">{frag("verdict.svg").replace('<svg ', '<svg style="min-width:720px" ', 1)}</div>
  <figcaption>{bi("判定规则，用封存检验里的真实区间画出。支持只需要有一类读出模型读得出；否定要求每一类都读不出，所以否定总是相对于声明的那几类读出模型。", "The verdict rule, drawn with real intervals from the sealed test. Support needs only one readout class to read the response; rejection needs every class to fail, so a rejection is always relative to the declared classes.")}</figcaption></figure>
  <p class="col">{bi("沿链把判定连起来，工具给出一个根因（按顺序取第一个成立的）：", "Reading the verdicts along the chain, the tool prints one root cause (the first that applies, in this order):")}</p>
  <div class="tw col"><table><thead><tr><th>{bi("根因", "Root cause")}</th><th>{bi("含义", "Meaning")}</th></tr></thead><tbody>{causes}</tbody></table></div>
</section>

<section id="checks">
  {sh(3, "不合法的设计，不给跑", "Illegal designs refuse to run")}
  <p class="col">{bi("统计之前先查设计。下面每一条都对应一种我们自己或别人犯过的错。", "The design is checked before any statistic. Each rule below stands for a mistake that we or others have made.")}</p>
  <div class="tw"><table><thead><tr><th>{bi("代号", "Code")}</th><th>{bi("规则", "Rule")}</th><th>{bi("防的是什么", "What it prevents")}</th></tr></thead><tbody>{checks}</tbody></table></div>
</section>

<section id="events">
  {sh(4, "类型化的源事件：事件图怎么建", "Typed source events: how the event graph is built")}
  <p class="col">{bi("一次运行的源不是一个标签，而是一组事件。每个事件写明：什么时候、打在哪个位点（哪个基因）、走哪个通道、方向和幅度，再加上角色、层、来源、合并方式和准入。<span class='mono'>python -m anm_diag events</span> 只读声明和已有的记录，不碰细胞数据，写出事件并报告三件事。", "The source of a run is not one label but a set of events. Each event states when, at which site (which gene), through which channel, with what sign and amplitude, plus its role, stratum, provenance, pooling and admission. <span class='mono'>python -m anm_diag events</span> reads declarations and finished records only, touches no cell data, writes the events and reports three things.")}</p>
  <ul class="pts">
    <li>{bi("<b>合并丢了哪些字段（规则 S2）。</b>把几个靶基因、几个批次并成一个处理组，等于丢掉“位点”或“层”这个字段。合并的判定不等于每个被合并的值的判定，所以合并要声明，否则给出警告。", "<b>Which fields a pooled source drops (rule S2).</b> Merging several target genes or batches into one treated arm drops the field “site” or “stratum”. A verdict for the pool is not a verdict for each pooled value, so pooling must be declared or a warning is raised.")}</li>
    <li>{bi("<b>防泄漏（规则 S3、S4、R2）。</b>每个事件带三条检查：细胞是怎么选出来的（按向导或设计，不能按实测的响应）；读出的测量是不是也进了系统的输入；事件打的是不是读出自己的基因。", "<b>No leakage (rules S3, S4, R2).</b> Each event carries three checks: how its cells were chosen (by guide or design, never by the measured response); whether the readout's measurement is also a system input; whether the event acts on the readout's own gene.")}</li>
    <li>{bi("<b>判定台账。</b>每个判定也写成事件，再列出不同运行之间相反的判定，以及它们在哪个字段上不同。", "<b>The verdict ledger.</b> Every verdict is written as an event too, and opposite verdicts between runs are listed with the field in which the runs differ.")}</li>
  </ul>
  <p class="note">{bi(f"用在封存的两次运行上：合并的源被标出来，于是对每个靶基因、每个批次单独重审（{gh('anm_diag/runs/split_20261111', 'runs/split_20261111')}）。PD-L1 的定位对每个靶基因都成立；15 次运行的 155 条判定里没有一对相反；两套数据的 25 个真实敲除事件都通过了三条防泄漏检查。", f"Applied to the two sealed runs: the pooled source is flagged, so the audit was repeated for each target gene and each batch alone ({gh('anm_diag/runs/split_20261111', 'runs/split_20261111')}). The PD-L1 location holds for every target gene; the 155 verdicts of 15 runs contain no opposite pair; all 25 real knockout events of the two data sets pass the three leakage checks.")}</p>
</section>

<section id="loop">
  {sh(5, "闭环", "The loop")}
  <p class="col">{bi("定位是观察，闭环是检验这个观察：照定位去改系统，再用同一个审计量一遍。", "A location is an observation; the loop tests it: change the system where the location points and measure again with the same audit.")}</p>
  <figure><div class="scroll">{loop_diagram()}</div>
  <figcaption>{bi("第 3 步同时写下一个按定位不该有用的改法，和该有用的那个用同样的数据、损失和步数各做一遍。第 5 步只用修复时留出的基因和没参与拟合的对照细胞。", "Step 3 also names a change that, by the location, should not help; it is made with the same data, loss and budget as the one that should. Step 5 uses only genes held out of the repair and control cells that no fit used.")}</figcaption></figure>
  <a class="next" href="anm-loop.html"><span class="t">{bi("在 teddy_mm 的头上走过的几轮", "The rounds run on the teddy_mm head")}</span><span class="go" aria-hidden="true">→</span></a>
</section>

<section id="grades">
  {sh(6, "证据等级", "Evidence grades")}
  <ul class="pts">
    <li>{bi("<b>封存、只跑一次：确认性。</b>声明连同数据和代码的哈希先封存；封存后改过设计、数据或代码，运行会被拒绝。", "<b>Sealed and run once: confirmatory.</b> The declaration is sealed with hashes of data and code; a run is refused if design, data or code changed after sealing.")}</li>
    <li>{bi("<b>设计先提交，但在看过别的结果之后才定：探索性。</b>闭环的各轮属于这一类，提交历史能看出先后。", "<b>Design committed first, but fixed after other results were seen: exploratory.</b> The rounds of the loop are of this kind; the commit history shows the order.")}</li>
    <li>{bi("<b>事后重读：探索性，原来的判定不变。</b>", "<b>Post hoc rereading: exploratory, and the original verdicts stand.</b>")}</li>
  </ul>
  <details class="more"><summary>{bi("形式定义（10 条）", "Formal definitions (10 items)")}</summary><div class="db"><ol class="defs">{frag("definitions.html")}</ol></div></details>
</section>"""
    page("anm-framework.html", "How the audit works — TEDDY × ANM",
         "How an ANM audit is declared, sealed and read: roles, the statistic u and its verdicts, the checks that refuse illegal designs, typed source events, the loop from location to repair, and evidence grades.", body)


# ---------------------------------------------------------------- loop
def readout_rows():
    out = []
    for ds, zh, en, *_ in DS:
        ros = [r for r in ("PDL1", "PDL2", "HLAE", "MHC1") if r in A[f"pooled_{ds}_as_trained"]["readouts"]]
        for ro in ros:
            cells = "".join(f'<td class="n">{f2(u_out(ds, arm, ro))}<br>{stat(po(ds, arm, ro)["status"])}</td>' for arm in ARM_OF_ROUND)
            out.append(f'<tr><td>{bi(zh, en)}</td><td class="k">{RO[ro]}</td><td class="n">{f2(po(ds, "as_trained", ro)["delta"])}</td><td class="n">{f2(po(ds, "as_trained", ro)["stages"]["Z"]["u"])}</td>{cells}</tr>')
    return "".join(out)


def gene_rows():
    out = []
    for ds, _, _, zh, en in DS:
        for g in GENES[ds]:
            def c(arm):
                r = A[f"target_{ds}_{g}_{arm}"]["readouts"]["PDL1"]
                return f'<td class="n">{f2(r["stages"]["HEAD_OUT"]["u"])}<br>{stat(r["status"])}</td>'
            r0 = A[f"target_{ds}_{g}_as_trained"]
            out.append(f'<tr><td>{bi(zh, en)} <b>{g}</b></td><td class="n">{r0["n_treated"]:,}</td><td class="n">{f2(r0["readouts"]["PDL1"]["stages"]["Z"]["u"])}</td>' + "".join(c(a) for a in ["as_trained", "decoder_only", "mlp_only", "mlp_guarded"] + (["joint_s0"] if HAVE3 else [])) + '</tr>')
    return "".join(out)


def hla_note():
    """HLA-E on melanoma: the head's own output against the measured protein, on held-out and on fitted genes."""
    if "HLAE" not in R["rule"]["frangieh2021"]["located_in_the_head"]:
        return "", ""
    v = [OS["frangieh2021"][g]["HLAE"] for g in GENES["frangieh2021"]]; lo = lambda k, sub=None: min((x[sub][k] if sub else x[k]) for x in v); hi = lambda k, sub=None: max((x[sub][k] if sub else x[k]) for x in v)
    assert lo("repaired") > 0 and lo("repaired", "fitted_genes") > 0 and hi("measured") < 0
    sd = lambda a, b: f"{a:+.2f} 到 {b:+.2f}".replace("-", "−"); se = lambda a, b: f"{a:+.2f} to {b:+.2f}".replace("-", "−")
    return (f"头自己的 HLA-E 输出在这些敲除上是升的（修复后 {sd(lo('repaired'), hi('repaired'))} 个对照标准差），实测蛋白略降（{sd(lo('measured'), hi('measured'))}）；在拟合用过的基因上也一样是升的（{sd(lo('repaired', 'fitted_genes'), hi('repaired', 'fitted_genes'))}）。也就是说，这次拟合在 400 步之内连训练过的基因上都没有学到 HLA-E 的响应，所以它说明的是这次拟合修不好，不是 MLP 一定带不动；原因没有查。",
            f"The head's own HLA-E output rises on these knockouts ({se(lo('repaired'), hi('repaired'))} control standard deviations after the repair) while the measured protein falls slightly ({se(lo('measured'), hi('measured'))}); it rises just as much on the genes the fit was trained on ({se(lo('repaired', 'fitted_genes'), hi('repaired', 'fitted_genes'))}). So within its 400 steps this fit did not learn the HLA-E response even where it was trained on it: the round shows that this fit does not repair HLA-E, not that the MLP cannot; the reason was not examined. ")


def round3_text():
    """Round 3 in words, in both languages: the main paragraph, what the single-gene audits leave open, and which head to use."""
    m, t = "frangieh2021", "papalexi2021"
    rule, rep = R["rule"], R["replication"]
    still = [(d, rule[d[0]]["located_in_the_head"]) for d in DS if rule[d[0]]["located_in_the_head"]]
    left = [(d, x) for d in DS for x in rule[d[0]]["left_outside_the_head"]]
    zh = ("<b>第 3 轮。</b>第 2 轮之后审计还指出两处，都在头里、都指向同一步：黑色素瘤上的 HLA-E，以及换一种细胞就失效。所以这一轮还是修 MLP，但把缺的补进拟合：两种细胞放在一起修成一个头，凡是实测了、头又有对应输出的蛋白都训（黑色素瘤上 PD-L1、MHC-I、HLA-E，THP-1 上 PD-L1）。六个靶基因轮流留出，在两套数据里同时留出。")
    en = ("<b>Round 3.</b> After round 2 the audit still located two things, both in the head and both at the same step: HLA-E on melanoma, and the failure on the other cell type. So this round repairs the MLP again and adds what the fit lacked: both cell types in one fit, giving one head, and every measured protein the head has an output for is trained (PD-L1, MHC-I and HLA-E on melanoma; PD-L1 on THP-1). The six target genes are held out in turn, in both data sets at once. ")
    if PASSED3:
        zh += f"按事先写定的规则：一个头同时守住了两种细胞。PD-L1 在头的输出上是 {f2(u_out(m, 'joint_s0'))}（黑色素瘤）和 {f2(u_out(t, 'joint_s0'))}（THP-1），工具不指出头里的任何一步；六个头在骨髓上是 {rng3(*bm_range('joint_s0'))}。"
        en += f"By the rule written in advance, one head holds for both cell types. PD-L1 at the head output is {f2(u_out(m, 'joint_s0'))} (melanoma) and {f2(u_out(t, 'joint_s0'))} (THP-1) and the tool names no step in the head; the six heads score {rng3(*bm_range('joint_s0'))} on bone marrow. "
    else:
        zh += "<b>按事先写定的规则，这一轮没有通过。</b>"; en += "<b>By the rule written in advance, this round did not pass.</b> "
        for d in DS:
            g = po(d[0], "joint_s0"); z = g["stages"]["Z"]["u"]; o = g["stages"]["HEAD_OUT"]["u"]; o2 = u_out(d[0], "mlp_guarded")
            if rule[d[0]]["outcome"] == "regressed":
                zh += f"{d[3]}：PD-L1 在头的输出上是 {f2(o)}，Z 处 {f2(z)}，头多丢了 {g['added_by_head']:.3f}，过了 0.10 的界限，读作“没定清”，按规则记为退步（第 2 轮是 {f2(o2)}）。"
                en += f"On {d[4]}, PD-L1 is {f2(o)} at the head output against {f2(z)} at Z: the head adds {g['added_by_head']:.3f}, over the margin of 0.10, which reads “not resolved” and counts as regressed under the rule (round 2: {f2(o2)}). "
            else:
                zh += f"{d[3]}：通过，PD-L1 是 {f2(o)}（Z 处 {f2(z)}），头里没有被指出的丢失。"
                en += f"On {d[4]} it passes: PD-L1 is {f2(o)} ({f2(z)} at Z) and no loss is located in the head. "
        zh += f"六个头在骨髓上是 {rng3(*bm_range('joint_s0'))}，在界限内。"; en += f"The six heads score {rng3(*bm_range('joint_s0'))} on bone marrow, within the guard. "
    if still:
        zh += "还被指出在头里丢的：" + "；".join(f"{d[3]}的 " + "、".join(f"{RO[x]}（输出 {f2(u_out(d[0], 'joint_s0', x))}，第 2 轮 {f2(u_out(d[0], 'mlp_guarded', x))}）" for x in v) for d, v in still) + "。按停止规则，这算作“改 MLP、用这些细胞修不好”，不再换别的拟合去试。" + hla_note()[0]
        en += "Still located in the head: " + "; ".join(", ".join(f"{RO[x]} (output {f2(u_out(d[0], 'joint_s0', x))}; round 2: {f2(u_out(d[0], 'mlp_guarded', x))})" for x in v) + f" on {d[4]}" for d, v in still) + ". By the stopping rule this counts as not repaired by changing the MLP with these cells, and no further variant of the fit is tried. " + hla_note()[1]
    else:
        zh += "审计不再指出头里有任何丢失。"; en += "The audit no longer locates any loss in the head. "
    if any(rep[d[0]] for d in DS):
        ok = lambda d: sum(1 for x in rep[d[0]].values() if x["outcome"] != "regressed")
        zh += "换两次对照细胞的划分重做，同一条规则：" + "；".join(f"{d[3]} {len(rep[d[0]])} 次里 {ok(d)} 次通过" for d in DS) + "。"
        en += "Repeated on two fresh splits of the control cells and read by the same rule: " + "; ".join(f"{d[4]} passes in {ok(d)} of {len(rep[d[0]])}" for d in DS) + ". "
    if not PASSED3:
        zh += "所以保留第 2 轮的修法，每种细胞一个头；“一个头管两种细胞”没有成立。"
        en += "So the repair of round 2 stands, with one head per cell type; one head for both cell types is not established."
    lzh = "；".join(f"{d[3]} {x['target']}（输出 {f2(x['u_output'])}，Z 处 {f2(x['u_Z'])}，{ST[x['status']][1]}）" for d, x in left); len_ = "; ".join(f"{d[4]} {x['target']} (output {f2(x['u_output'])}, {f2(x['u_Z'])} at Z, {ST[x['status']][2]})" for d, x in left)
    # which head to use
    fin = R["bone_marrow"]["final_head"]; r2 = TC.get("round_2_recipe") if TC else None
    if PASSED3:
        uzh = f"留出基因的那些头是证据。真正拿去用的头在全部六个基因上拟合，骨髓测试得分 {fin['mean_pearson']:.3f}（原来 {BASE:.3f}）。它没有在这些敲除细胞上再审计，因为它就是在这些细胞上修的。同样的修法用 teddy_mm 自己的模型代码写成了 <span class='mono'>scripts/07_repair_head.py</span>。"
        uen = f"The heads with a gene held out are the evidence. The head to use is fitted on all six genes; it scores {fin['mean_pearson']:.3f} on the bone-marrow test cells ({BASE:.3f} as trained). It is not audited on these knockout cells, because it was repaired on them. The same repair is written with teddy_mm's own model code as <span class='mono'>scripts/07_repair_head.py</span>."
    else:
        diff = max(r2[d[0]]["perturbation_cells"]["outputs_rms_difference_numpy_torch"] for d in DS); ch = [r2[d[0]]["perturbation_cells"]["outputs_rms_change_from_as_trained"] for d in DS]
        assert diff < 1e-3 and all(r2[d[0]]["decoder_identical_to_as_trained"] for d in DS)
        a, b, c, j0 = FE["as_trained"]["mean_pearson_mlp"], FE["melanoma_pdl1"]["mean_pearson_mlp"], FE["thp1_pdl1"]["mean_pearson_mlp"], FE["joint"]["mean_pearson_mlp"]
        uzh = (f"通过了规则的是第 2 轮的修法：每种细胞一个头，只训 PD-L1。<span class='mono'>scripts/07_repair_head.py</span> 用 teddy_mm 自己的模型代码做同样的修复；它做出的头和这里审计过的头对过，输出相差不到 0.001（修复本身带来的变化是 {f2(min(ch))} 和 {f2(max(ch))}）。"
               f"用它在全部基因上各修一个头，在完整的骨髓测试集（{FE['as_trained']['n_test']:,} 个细胞）上，teddy_mm 自己的评估是 {b:.3f}（黑色素瘤的头）和 {c:.3f}（THP-1 的头），原来 {a:.3f}。这两个头本身没有审计过，因为它们就是在这些敲除细胞上修的；证据是留出基因的那些头。"
               f"一个头只对修它时用的那种细胞有依据。这一轮的联合头（{j0:.3f}）没有通过规则，不建议用。")
        uen = (f"The repair that passed its rule is that of round 2: one head per cell type, trained on PD-L1 only. <span class='mono'>scripts/07_repair_head.py</span> makes the same repair with teddy_mm's own model code; its heads were compared with the ones audited here and their outputs differ by less than 0.001 (the repair itself moves them by {f2(min(ch))} and {f2(max(ch))}). "
               f"Fitted with it on all genes, one head per cell type, they score {b:.3f} (the melanoma head) and {c:.3f} (the THP-1 head) in teddy_mm's own evaluation on the full bone-marrow test split ({FE['as_trained']['n_test']:,} cells), against {a:.3f} as trained. These two heads are not audited themselves, because they were repaired on these knockout cells; the evidence is the heads with a gene held out. "
               f"A head is supported only for the cell type it was repaired on. The joint head of this round ({j0:.3f}) did not pass its rule and is not recommended.")
    return zh, en, lzh, len_, uzh, uen


def loop():
    m, t = "frangieh2021", "papalexi2021"
    n = {ds: (A[f"pooled_{ds}_as_trained"]["n_treated"], A[f"pooled_{ds}_as_trained"]["n_cells"] - A[f"pooled_{ds}_as_trained"]["n_treated"]) for ds, *_ in DS}
    assert cnt("as_trained", "located in the head") == 10 and cnt("decoder_only", "located in the head") == 10 and cnt("mlp_only", "located in the head") == 0 and cnt("mlp_guarded", "located in the head") == 0
    om, oa = own(m, "decoder_only"), own(m, "as_trained")
    round3 = ""
    if HAVE3:
        zh, en, lzh, len_, uzh, uen = round3_text()
        round3 = f"""<section id="round3">
  {sh(3, "第 3 轮：一个头，两种细胞" + ("" if PASSED3 else "；没有通过"), "Round 3: one head for both cell types" + ("" if PASSED3 else " — not passed"))}
  <p class="col">{bi(zh, en)}</p>
  <div class="tw"><table><thead><tr><th>{bi("细胞", "Cells")}</th><th>{bi("读出", "Readout")}</th><th>{bi("实测变化", "Measured change")}</th><th>Z</th>{"".join(f"<th>{bi('第 ' + str(i) + ' 轮', 'Round ' + str(i))}</th>" for i in range(len(ARM_OF_ROUND)))}</tr></thead><tbody>{readout_rows()}</tbody></table></div>
  <p class="src col">{bi(f"每格是头的输出上的 {U} 和读法。PD-L2 没有对应的头输出，从没进过损失。MHC-I 和 HLA-E 的实测变化很小，区间很宽。", f"Each cell gives {U} at the head output and the reading. PD-L2 has no output in the head and was never in a loss. MHC-I and HLA-E change little and their intervals are wide.")}</p>
  <p class="col">{bi("<b>逐个基因看，还没解决的：</b>" + (lzh or "没有") + "。工具在这些审计里没有指出头里的哪一步；其中 Z 处本身就高的，问题在头之前。", "<b>Gene by gene, what is not settled:</b> " + (len_ or "nothing") + ". In these audits the tool names no step in the head; where the value at Z is itself high, the problem lies before the head.")}</p>
  <h3>{bi("可以拿去用的头" if PASSED3 else "现在该用哪个头", "The head to use" if PASSED3 else "Which head to use now")}</h3>
  <p class="col">{bi(uzh, uen)}</p>
</section>"""
    body = f"""<section>
  <p class="eyebrow">{bi("闭环", "The loop")} · teddy_mm phase-1 head · {bi(UPDATED[0], UPDATED[1])}</p>
  <h1>{bi("照审计的定位去修，丢失就没了", "Repair where the audit points, and the loss goes")}</h1>
  <p class="lede col">{bi("审计说 PD-L1 的响应丢在头的第一段网络（MLP）。这一页是对这句话的检验：照它去修，再审；按它的说法不该有用的修法，也做一遍；然后看修复有没有弄坏头原来的任务。" + ("" if (not HAVE3 or PASSED3) else "前两轮通过了各自事先写定的规则；第 3 轮想更进一步，没有通过。"), "The audit says the PD-L1 response is lost in the first network of the head (its MLP). This page is the test of that statement: repair there and audit again; make the repair that should not help as well; then check whether the repair damaged the head's original task." + ("" if (not HAVE3 or PASSED3) else " The first two rounds passed the rules written for them in advance; the third, which went further, did not."))}</p>
  {rounds_table()}
</section>

<section id="rounds12">
  {sh(1, "第 1、2 轮：修审计指到的那一步", "Rounds 1 and 2: repair the step the audit names")}
  <p class="col">{bi(f"每套数据 5 个靶基因轮流留出一个：头在另外 4 个基因的敲除细胞上修，在留出基因的敲除细胞上审计；一半对照细胞不参加任何拟合，只用来审计（黑色素瘤 {n[m][0]:,} 个敲除细胞、{n[m][1]:,} 个对照；THP-1 {n[t][0]:,} 个、{n[t][1]:,} 个）。修法是在现有的头上微调：只动 MLP 的三层，解码器逐位不变；目标是实测的 PD-L1。", f"Each data set's 5 target genes are held out in turn: the head is repaired on knockout cells of the other 4 and audited on those of the held-out gene; half of the control cells take part in no fit and are used only for the audit (melanoma {n[m][0]:,} knockout cells and {n[m][1]:,} controls; THP-1 {n[t][0]:,} and {n[t][1]:,}). The repair fine-tunes the existing head: only the three MLP layers move, the decoder stays bit for bit; the target is measured PD-L1.")}</p>
  <figure><div class="scroll">{loop_chart()}</div>
  <figcaption>{bi(f"PD-L1，合并留出基因的审计。每条线是一个头在三环上的 {U}，竖线是区间。只修解码器的头在前两环就是原来的头，到输出仍在高处；修了 MLP 的头从中间层起落回 Z 的高度。" + ("" if (not HAVE3 or PASSED3) else "第 3 轮的头在黑色素瘤上到输出略高，过了规则的界限（见第 3 节）。"), f"PD-L1, the audit pooled over the held-out genes. Each line is one head's {U} at the three stages; bars are intervals. The head with only the decoder repaired is the head as trained at the first two stages and stays high at the output; the heads with the MLP repaired are back at the level of Z from the latent on." + ("" if (not HAVE3 or PASSED3) else " The head of round 3 ends slightly higher at the output on melanoma, over the rule's margin (section 3)."))}</figcaption></figure>
  <ul class="pts">
    <li>{bi(f"<b>按规则，两种细胞上都在定位到的那一步闭环。</b>输出上的 {U}：黑色素瘤 {f2(u_out(m, 'as_trained'))} → {f2(u_out(m, 'mlp_only'))}，只修解码器是 {f2(u_out(m, 'decoder_only'))}；THP-1 {f2(u_out(t, 'as_trained'))} → {f2(u_out(t, 'mlp_only'))}，只修解码器是 {f2(u_out(t, 'decoder_only'))}。", f"<b>By the rule, the loop closes at the located step on both cell types.</b> {U} at the output: melanoma {f2(u_out(m, 'as_trained'))} → {f2(u_out(m, 'mlp_only'))}, against {f2(u_out(m, 'decoder_only'))} with only the decoder repaired; THP-1 {f2(u_out(t, 'as_trained'))} → {f2(u_out(t, 'mlp_only'))}, against {f2(u_out(t, 'decoder_only'))}.")}</li>
    <li>{bi(f"<b>逐个基因。</b>原来的头，10 个留出基因都被指出在头里丢；只修解码器后还是 10 个；修 MLP 后 0 个（{cnt('mlp_only', 'no loss in the head')} 个头不再丢，{cnt('mlp_only', 'not resolved')} 个没定清）。", f"<b>Gene by gene.</b> As trained, the tool names a head step for all 10 held-out genes; with only the decoder repaired, still 10; with the MLP repaired, 0 ({cnt('mlp_only', 'no loss in the head')} read “no loss in the head” and {cnt('mlp_only', 'not resolved')} “not resolved”).")}</li>
    <li>{bi(f"<b>只看终点的平均值会选错。</b>只修解码器后，头自己的 PD-L1 输出在没见过的敲除上平均下降了实测的 {om[0]:.2f} 到 {om[1]:.2f} 倍（黑色素瘤；原来只有 {oa[0]:.2f} 到 {oa[1]:.2f}），看上去修好了。审计仍然否定它，因为它读的中间层没有变。", f"<b>The average at the endpoint would pick the wrong repair.</b> With only the decoder repaired, the head's own PD-L1 output falls for unseen knockouts by {om[0]:.2f} to {om[1]:.2f} times the measured fall (melanoma; {oa[0]:.2f} to {oa[1]:.2f} as trained), which looks repaired. The audit still rejects it, because the latent it reads from has not changed.")}</li>
  </ul>
  <details class="more wide"><summary>{bi("逐个留出基因的结果", "Held-out genes one by one")}</summary><div class="db"><div class="tw"><table><thead><tr><th>{bi("留出的基因", "Held-out gene")}</th><th>{bi("敲除细胞", "Knockout cells")}</th><th>Z</th><th>{bi("第 0 轮", "Round 0")}</th><th>{bi("只修解码器", "Decoder only")}</th><th>{bi("第 1 轮", "Round 1")}</th><th>{bi("第 2 轮", "Round 2")}</th>{("<th>" + bi("第 3 轮", "Round 3") + "</th>") if HAVE3 else ""}</tr></thead><tbody>{gene_rows()}</tbody></table></div>
  <p class="src col">{bi(f"每格是头的输出上的 {U} 和读法。“没定清”：工具没有指出头里的哪一步，但输出比 Z 高出 0.10 以上。“头不再丢”是和 Z 比：THP-1 的 JAK2 单看输出仍被否定，因为 Z 处已经是 0.15。", f"Each cell gives {U} at the head output and the reading. “Not resolved”: the tool names no step in the head, but the output is more than 0.10 above Z. “No loss in the head” compares with Z: for JAK2 on THP-1 the output alone is still rejected, with Z already at 0.15.")}</p></div></details>
</section>

<section id="cost">
  {sh(2, "代价：头原来的任务", "The cost: the head's original task")}
  <p class="col">{bi(f"头原来的任务是在骨髓细胞上从 RNA 预测 134 个蛋白。第 1 轮的修复把这项得分从 {BASE:.3f} 拉到 {rng3(*bm_range('mlp_only'))}，超过事先定的 0.02。第 2 轮在损失里加一项，把头在 3,000 个骨髓训练细胞上的输出按在原处；这一项的权重只在骨髓验证细胞上选。丢失照样去掉，得分是 {rng3(*bm_range('mlp_guarded'))}。", f"The head's original task is predicting 134 proteins from RNA on bone-marrow cells. The repair of round 1 pulled that score from {BASE:.3f} to {rng3(*bm_range('mlp_only'))}, beyond the 0.02 set in advance. Round 2 adds one term to the loss that holds the head's outputs on 3,000 bone-marrow training cells where they are; its weight is chosen on bone-marrow validation cells only. The loss still goes, and the score is {rng3(*bm_range('mlp_guarded'))}.")}</p>
  <figure><div class="scroll">{SC.cost(tx, cost_rows(), BASE)}</div>
  <figcaption>{bi("每个蛋白的预测和实测的相关，再取平均；5,000 个骨髓测试细胞。浅色带是事先定的界限。", "The correlation between predicted and measured, per protein, averaged; 5,000 bone-marrow test cells. The pale band is the guard set in advance.")}</figcaption></figure>
</section>

{round3}

<section id="limits">
  {sh(4 if HAVE3 else 3, "这些结果没有说明什么", "What these results do not show")}
  <ul class="pts">
    <li>{bi(f"<b>修好的头没有被判“支持”。</b>它只是不比 Z 多丢：这批细胞上 Z 自己就没定清（{f2(po(m, 'as_trained')['stages']['Z']['u'])} 和 {f2(po(t, 'as_trained')['stages']['Z']['u'])}）。", f"<b>The repaired head is not given “support”.</b> It only adds no loss beyond Z: on these cells Z itself is unresolved ({f2(po(m, 'as_trained')['stages']['Z']['u'])} and {f2(po(t, 'as_trained')['stages']['Z']['u'])}).")}</li>
    <li>{bi(f"<b>修复不会自己带到别的细胞类型。</b>第 1 轮里，在一种细胞上修的头到另一种细胞上仍被指出在头里丢（{f2(A['transfer_papalexi2021_mlp_only']['readouts']['PDL1']['stages']['HEAD_OUT']['u'])} 和 {f2(A['transfer_frangieh2021_mlp_only']['readouts']['PDL1']['stages']['HEAD_OUT']['u'])}）。没有审计过第三种细胞。" + ("" if (not HAVE3 or PASSED3) else "第 3 轮把两种细胞放进同一个拟合，也没有得到两边都通过的一个头。"), f"<b>A repair does not carry to another cell type by itself.</b> In round 1 a head repaired on one cell type is still located in the head on the other ({f2(A['transfer_papalexi2021_mlp_only']['readouts']['PDL1']['stages']['HEAD_OUT']['u'])} and {f2(A['transfer_frangieh2021_mlp_only']['readouts']['PDL1']['stages']['HEAD_OUT']['u'])}). No third cell type has been audited." + ("" if (not HAVE3 or PASSED3) else " Round 3 put both cell types into one fit and still did not give one head that passes on both."))}</li>
    <li>{bi("<b>这还不是封存的第 4 步。</b>那条规则要三个训练种子、一组只差这处改动的对照重训，并要求头的中间层得到“支持”。这里是在现有的头上微调，比的是原来的头。", "<b>This is not yet the sealed step 4.</b> That rule asks for three training seeds, a control retrain that differs only by the change, and “support” at the head latent. Here the existing head is fine-tuned and compared with the head as trained.")}</li>
    <li>{bi("<b>证据等级：探索性。</b>每套数据 5 个基因作重抽样单位，区间只能看方向。每一轮都是看过上一轮的结果之后才设计的；设计和规则在各自运行之前提交。", "<b>Evidence grade: exploratory.</b> Five genes per data set are the resampling units, so intervals show direction only. Each round was designed after the previous one was read; its design and rule were committed before it ran.")}</li>
  </ul>
</section>

<section id="records">
  {sh(5 if HAVE3 else 4, "记录", "Records")}
  <ul class="pts">
    <li>{gh("anm_diag/runs/loop_20261112")} · {gh("anm_diag/runs/loop_guard_20261112")}{(" · " + gh("anm_diag/runs/repair_20261112")) if HAVE3 else ""}{bi("：每一轮的设计（DESIGN.md）、每次审计的声明和记录、结果和脚本。", ": for each round the design (DESIGN.md), the declaration and record of every audit, the results and the scripts.")}</li>
    <li>{bi("细胞特征、头的权重和骨髓细胞不在仓库里，脚本通过环境变量读取。", "Cell features, head weights and bone-marrow cells are not in the repository; the scripts read them through environment variables.")}</li>
  </ul>
</section>"""
    page("anm-loop.html", "The loop on the teddy_mm head — TEDDY × ANM",
         "Closing the loop on the teddy_mm protein head: the ANM audit locates the loss of the PD-L1 response in the head's MLP; repairing that step removes it on held-out genes, repairing the decoder does not, and with the bone-marrow outputs held the original task is kept." + ("" if (not HAVE3 or PASSED3) else " A third round, one head for both cell types, did not pass its rule."), body)


if __name__ == "__main__":
    index(); framework(); loop()
