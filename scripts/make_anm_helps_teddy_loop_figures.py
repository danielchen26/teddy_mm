#!/usr/bin/env python3
"""NotebookLM-style poster panels: how the ANM Mode-B reverse decision loop helps TEDDY.

Large visual storytelling · icon-heavy · minimal text · verified number chips only.
Sources: docs/reports/reverse_loop_small_stats.json (and REVERSE_LOOP_SMALL.md).

Usage:
  /usr/bin/python3 scripts/make_anm_helps_teddy_loop_figures.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import (
    Arc,
    Circle,
    FancyBboxPatch,
    FancyArrowPatch,
    RegularPolygon,
    Rectangle,
    Wedge,
)
from matplotlib.colors import to_rgb, to_hex
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "assets" / "infographics"
FONT_DIR = ROOT / "docs" / "assets" / "fonts"

W, H, DPI = 1600, 900, 140

THEMES = {
    "light": dict(
        bg="#f7f6f2", panel="#ffffff", soft="#efeee8", ink="#1a1c1e", muted="#5a5f66",
        line="#d8d5cc", accent="#0d6e6e", accent_soft="#d8efef", warn="#8a4b12",
        warn_bg="#f5e6d4", good="#1f6b3a", good_bg="#dcefe3", bad="#8b2e2e",
        bad_bg="#f3e4e4", anm="#2a78d6", anm_soft="#dce8f8", jev="#eb6834",
        gray="#8c959f", rna="#0d6e6e", adt="#8a4b12", cell="#3d5a80",
        myeloid="#c45c26", tlin="#2a78d6",
    ),
    "dark": dict(
        bg="#121416", panel="#1b1e22", soft="#23272c", ink="#eceae4", muted="#a3a8b0",
        line="#2e333a", accent="#5ec4c4", accent_soft="#1a3333", warn="#e0a86a",
        warn_bg="#2f2418", good="#7dcea0", good_bg="#1a2e22", bad="#e08a8a",
        bad_bg="#2a1c1c", anm="#3987e5", anm_soft="#1a2a40", jev="#d95926",
        gray="#6e7681", rna="#5ec4c4", adt="#e0a86a", cell="#7ba3c9",
        myeloid="#e08a5a", tlin="#5aa0ef",
    ),
}


def mix(c1, c2, t):
    a, b = to_rgb(c1), to_rgb(c2)
    return to_hex(tuple((1 - t) * x + t * y for x, y in zip(a, b)))


def register_fonts():
    props = {}
    mapping = {
        "regular": "IBMPlexSans-Regular.ttf",
        "medium": "IBMPlexSans-Medium.ttf",
        "semibold": "IBMPlexSans-SemiBold.ttf",
        "bold": "IBMPlexSans-Bold.ttf",
        "mono": "IBMPlexMono-Regular.ttf",
        "mono_med": "IBMPlexMono-Medium.ttf",
    }
    arial = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
    arial_bold = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")
    for key, name in mapping.items():
        path = FONT_DIR / name
        if path.exists():
            fm.fontManager.addfont(str(path))
            props[key] = fm.FontProperties(fname=str(path))
        elif key in ("bold", "semibold") and arial_bold.exists():
            props[key] = fm.FontProperties(fname=str(arial_bold))
        elif arial.exists():
            props[key] = fm.FontProperties(fname=str(arial))
        else:
            props[key] = fm.FontProperties(family="DejaVu Sans")
    return props


FP = register_fonts()


def fp(weight="regular", size=12):
    p = FP.get(weight, FP["regular"]).copy()
    p.set_size(size)
    return p


def box(ax, x, y, w, h, fc, ec=None, lw=1.0, r=14, z=1, alpha=1.0):
    p = FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec or "none", linewidth=lw, alpha=alpha,
        zorder=z, mutation_aspect=1,
    )
    ax.add_patch(p)
    return p


def arrow(ax, x1, y1, x2, y2, color, lw=2.4, ms=16, z=5, style="-|>"):
    ax.annotate(
        "",
        xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=ms),
        zorder=z,
    )


def new_fig(theme):
    tk = THEMES[theme]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=tk["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ax.set_facecolor(tk["bg"])
    return fig, ax, tk


def header(ax, tk, eyebrow, title, subtitle=None):
    hh = 108 if subtitle else 82
    box(ax, 36, 24, W - 72, hh, tk["panel"], tk["line"], lw=1.0, r=18, z=1)
    box(ax, 36, 24, 10, hh, tk["accent"], r=5, z=2)
    ax.text(64, 50, eyebrow, fontproperties=fp("semibold", 12), color=tk["accent"], va="center", zorder=3)
    ax.text(64, 82, title, fontproperties=fp("bold", 26), color=tk["ink"], va="center", zorder=3)
    if subtitle:
        ax.text(64, 110, subtitle, fontproperties=fp("regular", 13), color=tk["muted"], va="center", zorder=3)


def footer(ax, tk, left, right="TEDDY × ANM · Mode B"):
    ax.text(36, H - 26, left, fontproperties=fp("regular", 11), color=tk["muted"], va="center", zorder=3)
    ax.text(W - 36, H - 26, right, fontproperties=fp("semibold", 11),
            color=tk["accent"], ha="right", va="center", zorder=3)


def save(fig, tk, slug, theme):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{slug}_{theme}.png"
    fig.savefig(path, dpi=DPI, facecolor=tk["bg"], bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return path


def chip(ax, x, y, w, h, value, label, tk, accent=None, z=3):
    accent = accent or tk["accent"]
    box(ax, x, y, w, h, tk["panel"], tk["line"], lw=1.0, r=14, z=z)
    ax.text(x + w / 2, y + h * 0.38, value, fontproperties=fp("bold", 22),
            color=accent, ha="center", va="center", zorder=z + 1)
    ax.text(x + w / 2, y + h * 0.72, label, fontproperties=fp("medium", 11),
            color=tk["muted"], ha="center", va="center", zorder=z + 1)


def stage_label(ax, x, y, text, color, z=4):
    ax.text(x, y, text, fontproperties=fp("semibold", 12), color=color,
            ha="center", va="center", zorder=z)


def draw_cell(ax, cx, cy, r, tk, fill=None, ec=None, z=3):
    fill = fill or mix(tk["cell"], tk["panel"], 0.55)
    ec = ec or tk["cell"]
    ax.add_patch(Circle((cx, cy), r, facecolor=fill, edgecolor=ec, linewidth=2.4, zorder=z))
    ax.add_patch(Circle((cx - r * 0.12, cy - r * 0.08), r * 0.38,
                         facecolor=mix(ec, tk["panel"], 0.45), edgecolor="none", zorder=z + 1))


def draw_lock(ax, cx, cy, color, size=26, z=4):
    s = size
    box(ax, cx - s * 0.45, cy - s * 0.1, s * 0.9, s * 0.7, color, r=4, z=z)
    ax.add_patch(Arc((cx, cy - s * 0.15), s * 0.7, s * 0.7, theta1=0, theta2=180,
                      color=color, lw=2.6, zorder=z))


def draw_yaml(ax, cx, cy, tk, z=3):
    box(ax, cx - 30, cy - 38, 60, 76, tk["panel"], tk["accent"], lw=1.8, r=7, z=z)
    for i, y in enumerate((cy - 18, cy - 2, cy + 14)):
        ax.plot([cx - 18, cx + 18], [y, y], color=tk["accent"] if i == 0 else tk["line"],
                lw=2.2 if i == 0 else 1.5, zorder=z + 1, solid_capstyle="round")


def draw_x(ax, cx, cy, color, size=16, z=5, lw=3.0):
    s = size
    ax.plot([cx - s, cx + s], [cy - s, cy + s], color=color, lw=lw, zorder=z, solid_capstyle="round")
    ax.plot([cx - s, cx + s], [cy + s, cy - s], color=color, lw=lw, zorder=z, solid_capstyle="round")


def draw_check(ax, cx, cy, color, size=18, z=5):
    ax.plot([cx - size, cx - size * 0.2, cx + size],
            [cy, cy + size * 0.7, cy - size * 0.7],
            color=color, lw=3.4, zorder=z, solid_capstyle="round", solid_joinstyle="round")


def draw_protein_dot(ax, cx, cy, r, color, z=3):
    ax.add_patch(Circle((cx, cy), r, facecolor=color, edgecolor="none", zorder=z, alpha=0.95))
    ax.add_patch(Circle((cx, cy), r * 0.4, facecolor="white", edgecolor="none", zorder=z + 1, alpha=0.35))


# ── PANEL 1 · LOOP POSTER ──────────────────────────────────────────────────

def draw_loop_poster(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(
        ax, tk,
        "ANM HELPS TEDDY  ·  MODE B DECISION LOOP",
        "Veto → Complement → Verify",
        "z_512 insufficient on must-pairs  ·  typed evidence + editable observer  ·  not Mode A",
    )

    # Big circular loop center
    cx, cy, R = 800, 520, 250
    # soft ring
    ax.add_patch(Circle((cx, cy), R + 18, facecolor=mix(tk["accent"], tk["bg"], 0.88),
                         edgecolor="none", zorder=1))
    ax.add_patch(Circle((cx, cy), R - 70, facecolor=tk["bg"], edgecolor="none", zorder=2))

    # three arc segments with arrows
    # VETO (top) · COMPLEMENT (bottom-left) · VERIFY (bottom-right)
    nodes = [
        ("VETO", 90, tk["bad"], "z_512\ncan't separate"),
        ("COMPLEMENT", 330, tk["anm"], "true ADT\n+ O1 YAML"),
        ("VERIFY", 210, tk["good"], "false_agree\n↓  ·  Q ↑"),
    ]
    for name, ang, color, sub in nodes:
        rad = np.deg2rad(ang)
        # node center on ring
        nx = cx + R * np.cos(rad)
        ny = cy - R * np.sin(rad)
        ax.add_patch(Circle((nx, ny), 78, facecolor=tk["panel"], edgecolor=color, linewidth=3.0, zorder=4))
        ax.add_patch(Circle((nx, ny), 78, facecolor=mix(color, tk["panel"], 0.82), edgecolor="none", zorder=3))
        ax.text(nx, ny - 14, name, fontproperties=fp("bold", 15), color=color, ha="center", va="center", zorder=5)
        ax.text(nx, ny + 22, sub, fontproperties=fp("medium", 11), color=tk["ink"], ha="center", va="center", zorder=5,
                linespacing=1.25)

    # curved direction arrows between nodes (approximate with FancyArrowPatch)
    def arc_arrow(a1, a2, color):
        # mid angle
        mid = (a1 + a2) / 2 if a2 > a1 else (a1 + a2 + 360) / 2
        if mid >= 360:
            mid -= 360
        r = R
        # start/end slightly inside arc
        s = np.deg2rad(a1 - 28)
        e = np.deg2rad(a2 + 28)
        x1, y1 = cx + r * np.cos(s), cy - r * np.sin(s)
        x2, y2 = cx + r * np.cos(e), cy - r * np.sin(e)
        # control outward
        mrad = np.deg2rad(mid)
        xm, ym = cx + (r + 55) * np.cos(mrad), cy - (r + 55) * np.sin(mrad)
        path = FancyArrowPatch(
            (x1, y1), (x2, y2),
            connectionstyle=f"arc3,rad={0.35 if a2 > a1 else -0.35}",
            arrowstyle="-|>", mutation_scale=22, lw=3.2, color=color, zorder=3,
        )
        # simpler radial chevrons instead
        ax.add_patch(path)

    # simpler chevrons along the ring
    for ang, color in ((30, tk["bad"]), (270, tk["anm"]), (150, tk["good"])):
        rad = np.deg2rad(ang)
        x = cx + (R - 10) * np.cos(rad)
        y = cy - (R - 10) * np.sin(rad)
        # tangent direction (clockwise for loop: decreasing angle in our y-flip)
        tang = rad - np.pi / 2
        dx, dy = 28 * np.cos(tang), -28 * np.sin(tang)
        arrow(ax, x - dx, y - dy, x + dx, y + dy, color, lw=3.0, ms=18, z=3)

    # center badge
    box(ax, cx - 110, cy - 55, 220, 110, tk["panel"], tk["accent"], lw=2.0, r=18, z=4)
    ax.text(cx, cy - 18, "TEDDY", fontproperties=fp("bold", 20), color=tk["ink"], ha="center", zorder=5)
    ax.text(cx, cy + 18, "frozen z_512", fontproperties=fp("medium", 13), color=tk["muted"], ha="center", zorder=5)
    draw_lock(ax, cx, cy + 48, tk["accent"], size=18, z=5)

    # side chips — verified only
    chip(ax, 48, 170, 200, 100, "5,000", "must-pairs", tk, tk["warn"])
    chip(ax, 48, 290, 200, 100, "1,490", "myeloid ↔ T", tk, tk["myeloid"])
    chip(ax, 48, 410, 200, 100, "0.982", "median cos", tk, tk["bad"])
    chip(ax, 48, 530, 200, 100, "VETO", "z_512 alone", tk, tk["bad"])

    chip(ax, 1352, 170, 210, 100, "0.40", "false_agree ↓", tk, tk["good"])
    chip(ax, 1352, 290, 210, 100, "0.99", "Q decided", tk, tk["good"])
    chip(ax, 1352, 410, 210, 100, "0.63", "soft_sep ↑", tk, tk["anm"])
    chip(ax, 1352, 530, 210, 100, "ADT+O1", "decision win", tk, tk["anm"])

    # bottom honesty strip
    box(ax, 280, 790, 1040, 58, tk["warn_bg"], mix(tk["warn"], tk["line"], 0.5), lw=1.2, r=12, z=3)
    ax.text(800, 819, "Mode B decision loop  ·  not Mode A residual/Jacobian/gene-perturb  ·  not Pearson > 0.61",
            fontproperties=fp("semibold", 12.5), color=tk["warn"], ha="center", va="center", zorder=4)

    footer(ax, tk, "reverse_loop_small_stats.json · myeloid↔T focus")
    return save(fig, tk, "06_anm_loop_poster", theme)


# ── PANEL 2 · MYELOID ↔ T BEFORE / AFTER ───────────────────────────────────

def draw_myeloid_t(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(
        ax, tk,
        "MYELOID ↔ T  ·  MUST-SEPARATE STORY",
        "Same near-identical z. Different true lineage.",
        "Baseline pred+O0 silently agrees  ·  true-ADT typed + O1 separates",
    )

    # LEFT: before
    box(ax, 48, 155, 720, 620, tk["bad_bg"], tk["bad"], lw=2.0, r=20, z=1)
    stage_label(ax, 408, 195, "BEFORE  ·  TEDDY PRED + O0", tk["bad"])

    # two near cells overlapping
    draw_cell(ax, 280, 360, 70, tk, fill=mix(tk["myeloid"], tk["panel"], 0.55), ec=tk["myeloid"])
    draw_cell(ax, 380, 360, 70, tk, fill=mix(tk["tlin"], tk["panel"], 0.55), ec=tk["tlin"])
    # blur / merge hint
    ax.add_patch(Circle((330, 360), 95, facecolor=mix(tk["bad"], tk["panel"], 0.75),
                         edgecolor="none", alpha=0.35, zorder=2))
    ax.text(330, 470, "cos ≈ 0.982", fontproperties=fp("bold", 18), color=tk["bad"], ha="center", zorder=4)
    ax.text(330, 500, "near-identical z_512", fontproperties=fp("medium", 13), color=tk["ink"], ha="center", zorder=4)

    # false agree icon
    box(ax, 120, 540, 580, 200, tk["panel"], tk["bad"], lw=1.6, r=16, z=3)
    ax.text(410, 590, "FALSE AGREE", fontproperties=fp("bold", 22), color=tk["bad"], ha="center", zorder=4)
    ax.text(410, 635, "0.85", fontproperties=fp("bold", 48), color=tk["bad"], ha="center", zorder=4)
    ax.text(410, 695, "soft_sep 0.17   ·   Q 0.74", fontproperties=fp("medium", 14),
            color=tk["muted"], ha="center", zorder=4)

    # RIGHT: after
    box(ax, 832, 155, 720, 620, tk["good_bg"], tk["good"], lw=2.0, r=20, z=1)
    stage_label(ax, 1192, 195, "AFTER  ·  TRUE ADT + O1", tk["good"])

    draw_cell(ax, 1040, 340, 62, tk, fill=mix(tk["myeloid"], tk["panel"], 0.45), ec=tk["myeloid"])
    draw_cell(ax, 1340, 340, 62, tk, fill=mix(tk["tlin"], tk["panel"], 0.45), ec=tk["tlin"])
    # protein markers between
    for i, (lab, col) in enumerate((("CD16", tk["myeloid"]), ("CD3", tk["tlin"]), ("CD5", tk["tlin"]))):
        px = 1100 + i * 70
        draw_protein_dot(ax, px, 340, 16, col)
        ax.text(px, 372, lab, fontproperties=fp("semibold", 10), color=col, ha="center", zorder=5)
    # separation arrow
    arrow(ax, 1110, 340, 1270, 340, tk["good"], lw=3.2, ms=20, z=4)
    ax.text(1190, 300, "typed ADT evidence", fontproperties=fp("semibold", 13),
            color=tk["good"], ha="center", zorder=5)
    ax.text(1190, 420, "observer O1 (abstain gate)", fontproperties=fp("medium", 12),
            color=tk["ink"], ha="center", zorder=5)

    box(ax, 904, 540, 580, 200, tk["panel"], tk["good"], lw=1.6, r=16, z=3)
    ax.text(1194, 590, "FALSE AGREE ↓", fontproperties=fp("bold", 22), color=tk["good"], ha="center", zorder=4)
    ax.text(1194, 635, "0.40", fontproperties=fp("bold", 48), color=tk["good"], ha="center", zorder=4)
    ax.text(1194, 695, "soft_sep 0.63   ·   Q 0.99", fontproperties=fp("medium", 14),
            color=tk["muted"], ha="center", zorder=4)

    # tiny chips under header area inside panels for n
    ax.text(408, 230, "n = 1,490 pairs", fontproperties=fp("medium", 12), color=tk["muted"], ha="center", zorder=4)
    ax.text(1192, 230, "winner complement_b_true_ADT_O1", fontproperties=fp("medium", 12),
            color=tk["muted"], ha="center", zorder=4)

    footer(ax, tk, "VERIFY myeloid↔T · reverse_loop_small_stats.json")
    return save(fig, tk, "07_myeloid_t_before_after", theme)


# ── PANEL 3 · MISSING-LINK REDIRECT (ANM vs feature swap) ─────────────────

def draw_missing_link(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(
        ax, tk,
        "WHAT ANM UNIQUELY PROVIDES",
        "Not just swapping features.",
        "Typed evidence · editable observers · abstain · verify  —  vs z+ geometry alone",
    )

    # Four unique ANM pillars
    pillars = [
        ("TYPED\nEVIDENCE", "δu events\ninto finite_field", tk["anm"], 0),
        ("EDITABLE\nOBSERVER", "O0->O1 YAML\nzero retrain", tk["accent"], 1),
        ("ABSTAIN", "silence when\nz insufficient", tk["warn"], 2),
        ("VERIFY", "false_agree · Q\nsoft_sep chips", tk["good"], 3),
    ]
    for title, sub, color, i in pillars:
        x = 56 + i * 380
        box(ax, x, 155, 350, 280, tk["panel"], color, lw=2.2, r=18, z=2)
        ax.add_patch(Circle((x + 175, 230), 36, facecolor=mix(color, tk["panel"], 0.75),
                             edgecolor=color, linewidth=2.4, zorder=3))
        ax.text(x + 175, 230, str(i + 1), fontproperties=fp("bold", 20), color=color, ha="center", va="center", zorder=4)
        ax.text(x + 175, 300, title, fontproperties=fp("bold", 16), color=tk["ink"], ha="center", va="center",
                zorder=4, linespacing=1.15)
        ax.text(x + 175, 370, sub, fontproperties=fp("medium", 13), color=tk["muted"], ha="center", va="center",
                zorder=4, linespacing=1.25)

    # Bottom comparison: feature swap FAIL vs ANM path WIN
    box(ax, 56, 470, 720, 310, tk["bad_bg"], tk["bad"], lw=2.0, r=18, z=2)
    stage_label(ax, 416, 510, "FEATURE SWAP ALONE  ·  z+", tk["bad"])
    draw_x(ax, 200, 600, tk["bad"], size=28, lw=4.0)
    ax.text(300, 600, "pred z+ still cos≥0.98\non 83% of myeloid↔T", fontproperties=fp("semibold", 15),
            color=tk["ink"], va="center", zorder=4)
    ax.text(416, 670, "global relative_must_reduction = −0.27", fontproperties=fp("bold", 16),
            color=tk["bad"], ha="center", zorder=4)
    ax.text(416, 710, "z+ alone worsens must-rate among near", fontproperties=fp("medium", 13),
            color=tk["muted"], ha="center", zorder=4)
    ax.text(416, 745, "geometry ≠ decision ownership", fontproperties=fp("semibold", 13),
            color=tk["bad"], ha="center", zorder=4)

    box(ax, 824, 470, 720, 310, tk["good_bg"], tk["good"], lw=2.0, r=18, z=2)
    stage_label(ax, 1184, 510, "ANM DECISION PATH", tk["good"])
    draw_check(ax, 980, 600, tk["good"], size=26)
    ax.text(1080, 600, "true ADT typed + O1\n→ decision winner", fontproperties=fp("semibold", 15),
            color=tk["ink"], va="center", zorder=4, linespacing=1.2)
    ax.text(1184, 670, "false_agree 0.85 → 0.40    Q 0.74 → 0.99", fontproperties=fp("bold", 15),
            color=tk["good"], ha="center", zorder=4)
    ax.text(1184, 710, "soft_sep 0.17 → 0.63   ·   Mode B only", fontproperties=fp("medium", 13),
            color=tk["muted"], ha="center", zorder=4)
    ax.text(1184, 745, "holdout ADT = complement probe (not production leak)", fontproperties=fp("semibold", 12),
            color=tk["warn"], ha="center", zorder=4)

    footer(ax, tk, "ANM_HELPS_TEDDY_LOOP · complements a/b/c from reverse_loop_small")
    return save(fig, tk, "08_anm_vs_feature_swap", theme)


def main():
    drawers = [draw_loop_poster, draw_myeloid_t, draw_missing_link]
    written = []
    for draw in drawers:
        for theme in ("light", "dark"):
            p = draw(theme)
            written.append(p)
            print("wrote", p, p.stat().st_size)
    print(f"done: {len(written)} PNGs → {OUT}")


if __name__ == "__main__":
    main()
