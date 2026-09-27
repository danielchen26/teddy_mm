#!/usr/bin/env python3
"""Figure-first story infographics for TEDDY×ANM docs.

Icons / arrows / cell·protein glyphs carry the story. Minimal labels.
Verified numbers only (HARD_PROOF / SCOPE_REFINE / MISSING_MODALITY).

Usage:
  python scripts/make_story_infographics.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Arc, Wedge
from matplotlib.colors import to_rgb, to_hex
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "assets" / "infographics"
FONT_DIR = ROOT / "docs" / "assets" / "fonts"
for _cand in (FONT_DIR, Path("/workspace/teddy_fonts"), Path("/tmp/teddy_fonts")):
    if _cand.exists() and any(_cand.glob("*.ttf")):
        FONT_DIR = _cand
        break

W, H, DPI = 1600, 900, 140

THEMES = {
    "light": dict(
        bg="#f7f6f2", panel="#ffffff", soft="#efeee8", ink="#1a1c1e", muted="#5a5f66",
        line="#d8d5cc", accent="#0d6e6e", accent_soft="#d8efef", warn="#8a4b12",
        warn_bg="#f5e6d4", good="#1f6b3a", good_bg="#dcefe3", bad="#8b2e2e",
        bad_bg="#f3e4e4", anm="#2a78d6", jev="#eb6834", gray="#8c959f",
        rna="#0d6e6e", adt="#8a4b12", cell="#3d5a80",
    ),
    "dark": dict(
        bg="#121416", panel="#1b1e22", soft="#23272c", ink="#eceae4", muted="#a3a8b0",
        line="#2e333a", accent="#5ec4c4", accent_soft="#1a3333", warn="#e0a86a",
        warn_bg="#2f2418", good="#7dcea0", good_bg="#1a2e22", bad="#e08a8a",
        bad_bg="#2a1c1c", anm="#3987e5", jev="#d95926", gray="#6e7681",
        rna="#5ec4c4", adt="#e0a86a", cell="#7ba3c9",
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
    for key, name in mapping.items():
        path = FONT_DIR / name
        if path.exists():
            fm.fontManager.addfont(str(path))
            props[key] = fm.FontProperties(fname=str(path))
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


def arrow(ax, x1, y1, x2, y2, color, lw=2.0, ms=14, z=5):
    ax.annotate(
        "",
        xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=ms),
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
    box(ax, 36, 28, W - 72, 100 if subtitle else 78, tk["panel"], tk["line"], lw=1.0, r=16, z=1)
    box(ax, 36, 28, 8, 100 if subtitle else 78, tk["accent"], r=4, z=2)
    ax.text(64, 52, eyebrow, fontproperties=fp("semibold", 11), color=tk["accent"], va="center", zorder=3)
    ax.text(64, 86, title, fontproperties=fp("bold", 24), color=tk["ink"], va="center", zorder=3)
    if subtitle:
        ax.text(64, 112, subtitle, fontproperties=fp("regular", 12.5), color=tk["muted"], va="center", zorder=3)


def footer(ax, tk, left, right="TEDDY × ANM"):
    ax.text(36, H - 28, left, fontproperties=fp("regular", 11), color=tk["muted"], va="center", zorder=3)
    ax.text(W - 36, H - 28, right, fontproperties=fp("semibold", 11),
            color=tk["accent"], ha="right", va="center", zorder=3)


def save(fig, tk, slug, theme):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{slug}_{theme}.png"
    fig.savefig(path, dpi=DPI, facecolor=tk["bg"], bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return path


# ── glyph helpers ──────────────────────────────────────────────────────────

def draw_cell(ax, cx, cy, r, tk, fill=None, ec=None, label=None, z=3):
    fill = fill or mix(tk["cell"], tk["panel"], 0.55)
    ec = ec or tk["cell"]
    ax.add_patch(Circle((cx, cy), r, facecolor=fill, edgecolor=ec, linewidth=2.2, zorder=z))
    # nucleus
    ax.add_patch(Circle((cx - r * 0.12, cy - r * 0.08), r * 0.38,
                         facecolor=mix(ec, tk["panel"], 0.45), edgecolor="none", zorder=z + 1))
    if label:
        ax.text(cx, cy + r + 18, label, fontproperties=fp("semibold", 11),
                color=tk["ink"], ha="center", va="center", zorder=z + 2)


def draw_protein(ax, cx, cy, r, color, label=None, z=3):
    ax.add_patch(Circle((cx, cy), r, facecolor=color, edgecolor="none", zorder=z, alpha=0.92))
    ax.add_patch(Circle((cx, cy), r * 0.45, facecolor="white", edgecolor="none", zorder=z + 1, alpha=0.35))
    if label:
        ax.text(cx, cy + r + 16, label, fontproperties=fp("medium", 10),
                color=color, ha="center", va="center", zorder=z + 2)


def draw_rna(ax, cx, cy, color, scale=1.0, z=3):
    """Simple helix-ish RNA glyph."""
    t = np.linspace(0, 2.2 * np.pi, 60)
    x = cx + 7 * scale * np.sin(t)
    y = cy - 38 * scale + (t / (2.2 * np.pi)) * 76 * scale
    ax.plot(x, y, color=color, lw=2.4 * scale, zorder=z, solid_capstyle="round")
    ax.plot(cx - 7 * scale * np.sin(t), y, color=mix(color, "#ffffff", 0.35),
            lw=1.6 * scale, zorder=z, solid_capstyle="round", alpha=0.85)


def draw_lock(ax, cx, cy, color, size=22, z=4):
    """Frozen / lock glyph."""
    s = size
    box(ax, cx - s * 0.45, cy - s * 0.1, s * 0.9, s * 0.7, color, r=4, z=z)
    arc = Arc((cx, cy - s * 0.15), s * 0.7, s * 0.7, theta1=0, theta2=180,
              color=color, lw=2.4, zorder=z)
    ax.add_patch(arc)
    ax.add_patch(Circle((cx, cy + s * 0.15), s * 0.12, facecolor=mix(color, "#ffffff", 0.5),
                         edgecolor="none", zorder=z + 1))


def draw_yaml(ax, cx, cy, tk, z=3):
    """YAML doc glyph."""
    box(ax, cx - 28, cy - 36, 56, 72, tk["panel"], tk["accent"], lw=1.6, r=6, z=z)
    for i, y in enumerate((cy - 18, cy - 2, cy + 14)):
        ax.plot([cx - 16, cx + 16], [y, y], color=tk["accent"] if i == 0 else tk["line"],
                lw=2.0 if i == 0 else 1.4, zorder=z + 1, solid_capstyle="round")


def draw_x(ax, cx, cy, color, size=18, z=5):
    s = size
    ax.plot([cx - s, cx + s], [cy - s, cy + s], color=color, lw=3.2, zorder=z, solid_capstyle="round")
    ax.plot([cx - s, cx + s], [cy + s, cy - s], color=color, lw=3.2, zorder=z, solid_capstyle="round")


def metric_chip(ax, x, y, w, h, value, label, tk, accent=None, z=3):
    accent = accent or tk["accent"]
    box(ax, x, y, w, h, tk["panel"], tk["line"], lw=1.0, r=12, z=z)
    ax.text(x + w / 2, y + h * 0.38, value, fontproperties=fp("bold", 20),
            color=accent, ha="center", va="center", zorder=z + 1)
    ax.text(x + w / 2, y + h * 0.72, label, fontproperties=fp("medium", 10),
            color=tk["muted"], ha="center", va="center", zorder=z + 1)


def stage_label(ax, x, y, text, color, z=4):
    ax.text(x, y, text, fontproperties=fp("semibold", 11), color=color,
            ha="center", va="center", zorder=z)


# ── CASE 00 · GLOBAL FRAME ─────────────────────────────────────────────────

def draw_00(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "GLOBAL FRAME  ·  CITE / TEDDY FROZEN MAP",
           "RNA → protein evidence. Editable decisions on top.",
           "site4/test n=16,750  ·  phase-1 full-ADT Pearson ~0.61")

    # Three visual stages
    stages = [
        (180, "CITE cell", tk["rna"]),
        (560, "Frozen TEDDY map", tk["gray"]),
        (980, "ANM decision field", tk["good"]),
        (1380, "Edit / abstain", tk["anm"]),
    ]
    top = 200
    card_h = 520

    # Stage 1: cell with RNA → proteins
    box(ax, 60, top, 280, card_h, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 200, top + 36, "1 · CARE", tk["accent"])
    draw_cell(ax, 200, top + 200, 70, tk, label="CITE cell")
    draw_rna(ax, 130, top + 200, tk["rna"], scale=0.85)
    for i, (dx, lab, c) in enumerate([(-55, "CD", tk["adt"]), (0, "CD", mix(tk["adt"], tk["accent"], 0.4)),
                                       (55, "CD", tk["accent"])]):
        draw_protein(ax, 200 + dx, top + 340, 16, c)
    ax.text(200, top + 390, "RNA → protein / state", fontproperties=fp("medium", 12),
            color=tk["ink"], ha="center", zorder=4)
    ax.text(200, top + 420, "What biologists want", fontproperties=fp("regular", 11),
            color=tk["muted"], ha="center", zorder=4)

    # Stage 2: frozen map
    box(ax, 380, top, 300, card_h, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 530, top + 36, "2 · TEDDY GIVES", tk["muted"])
    # pipeline pills
    pills = [("RNA", tk["rna"]), ("TEDDY-G", tk["gray"]), ("z₅₁₂", tk["cell"]), ("ADT", tk["adt"])]
    py = top + 160
    for i, (lab, c) in enumerate(pills):
        px = 420 + i * 60
        box(ax, px, py, 52, 36, mix(c, tk["panel"], 0.7), c, lw=1.2, r=10, z=2)
        ax.text(px + 26, py + 18, lab, fontproperties=fp("semibold", 9),
                color=tk["ink"], ha="center", va="center", zorder=3)
        if i < len(pills) - 1:
            arrow(ax, px + 54, py + 18, px + 58, py + 18, tk["line"], lw=1.4, ms=10)
    draw_lock(ax, 530, top + 280, tk["gray"], size=36)
    ax.text(530, top + 340, "FROZEN", fontproperties=fp("bold", 14),
            color=tk["gray"], ha="center", zorder=4)
    metric_chip(ax, 430, top + 390, 200, 90, "~0.61", "Pearson (phase-1)", tk, tk["accent"])

    # Stage 3: ANM field
    box(ax, 720, top, 300, card_h, tk["good_bg"], mix(tk["good"], tk["line"], 0.4), lw=1.4, r=18, z=1)
    stage_label(ax, 870, top + 36, "3 · ANM ADDS", tk["good"])
    draw_yaml(ax, 870, top + 180, tk)
    ax.text(870, top + 260, "YAML observers", fontproperties=fp("semibold", 13),
            color=tk["ink"], ha="center", zorder=4)
    # decision chips
    for i, (lab, c) in enumerate([("call", tk["good"]), ("abstain", tk["warn"]),
                                   ("markers", tk["anm"]), ("scope", tk["accent"])]):
        bx = 750 + (i % 2) * 120
        by = top + 310 + (i // 2) * 70
        box(ax, bx, by, 108, 52, tk["panel"], c, lw=1.4, r=12, z=2)
        ax.text(bx + 54, by + 26, lab, fontproperties=fp("semibold", 13),
                color=c, ha="center", va="center", zorder=3)

    # Stage 4: no retrain
    box(ax, 1060, top, 480, card_h, tk["accent_soft"], mix(tk["accent"], tk["line"], 0.35), lw=1.4, r=18, z=1)
    stage_label(ax, 1300, top + 36, "4 · WITHOUT", tk["accent"])
    # crossed-out retrain + big check for edit
    ax.text(1300, top + 160, "NO  TEDDY retrain", fontproperties=fp("bold", 18),
            color=tk["bad"], ha="center", zorder=4)
    ax.text(1300, top + 210, "NO  Pearson > ~0.61", fontproperties=fp("bold", 18),
            color=tk["bad"], ha="center", zorder=4)
    ax.text(1300, top + 260, "NO  Clinical claim", fontproperties=fp("bold", 18),
            color=tk["bad"], ha="center", zorder=4)
    box(ax, 1140, top + 320, 320, 120, tk["panel"], tk["good"], lw=2.0, r=14, z=2)
    ax.text(1300, top + 360, "YES  Declaration edit", fontproperties=fp("bold", 18),
            color=tk["good"], ha="center", zorder=3)
    ax.text(1300, top + 400, "same frozen δu", fontproperties=fp("regular", 13),
            color=tk["muted"], ha="center", zorder=3)

    # flow arrows between cards
    for x in (340, 680, 1020):
        arrow(ax, x, top + card_h / 2, x + 36, top + card_h / 2, tk["line"], lw=2.2, ms=16)

    footer(ax, tk, "Not clinical · not a Pearson win · holdout adt_true = verifier only")
    return save(fig, tk, "00_global_frame", theme)


# ── CASE 01 · EDIT QUESTION ────────────────────────────────────────────────

def draw_01(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "CASE 1  ·  EDIT THE QUESTION  ·  HARD_PROOF",
           "Change the criterion — keep the frozen map.",
           "Same δu  ·  O0 equal-panel → O1 tighter abstain → O2 key-marker")

    # Left: frozen TEDDY + cell
    box(ax, 48, 160, 340, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 218, 196, "FROZEN EVIDENCE", tk["muted"])
    draw_cell(ax, 218, 340, 55, tk)
    draw_lock(ax, 218, 440, tk["gray"], size=28)
    ax.text(218, 490, "TEDDY δu", fontproperties=fp("semibold", 14),
            color=tk["ink"], ha="center", zorder=4)
    ax.text(218, 520, "no retrain", fontproperties=fp("regular", 12),
            color=tk["muted"], ha="center", zorder=4)
    # lineage glyphs
    for i, (lab, c) in enumerate([("B", tk["anm"]), ("T", tk["accent"]), ("M", tk["adt"])]):
        bx = 100 + i * 80
        box(ax, bx, 580, 64, 64, mix(c, tk["panel"], 0.7), c, lw=1.6, r=14, z=2)
        ax.text(bx + 32, bx and 612, lab, fontproperties=fp("bold", 18),
                color=c, ha="center", va="center", zorder=3)
        ax.text(bx + 32, 612, lab, fontproperties=fp("bold", 18),
                color=c, ha="center", va="center", zorder=3)
    ax.text(218, 700, "lineage call", fontproperties=fp("medium", 12),
            color=tk["muted"], ha="center", zorder=4)

    # Center: three observers as dials
    box(ax, 420, 160, 720, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 780, 196, "YAML OBSERVERS", tk["accent"])

    observers = [
        ("O0", "equal weight", 317, tk["accent"], 0.25),
        ("O1", "tighter abstain", 834, tk["anm"], 0.45),
        ("O2", "key markers", 3373, tk["good"], 0.95),
    ]
    max_a = 3373
    for i, (name, desc, abstain, color, fill) in enumerate(observers):
        cx = 520 + i * 220
        # dial circle
        ax.add_patch(Circle((cx, 380), 70, facecolor=mix(color, tk["panel"], 0.75),
                             edgecolor=color, linewidth=3, zorder=2))
        # fill wedge showing abstain intensity
        wedge = Wedge((cx, 380), 70, 90, 90 - 360 * fill,
                      facecolor=mix(color, tk["panel"], 0.4), edgecolor="none", zorder=3)
        ax.add_patch(wedge)
        ax.text(cx, 370, name, fontproperties=fp("bold", 22),
                color=color, ha="center", va="center", zorder=4)
        ax.text(cx, 400, desc, fontproperties=fp("regular", 10),
                color=tk["muted"], ha="center", va="center", zorder=4)
        # abstain bar
        bar_h = 180 * (abstain / max_a)
        box(ax, cx - 28, 520 + (180 - bar_h), 56, bar_h, color, r=8, z=2)
        ax.text(cx, 720, f"{abstain:,}", fontproperties=fp("bold", 16),
                color=color, ha="center", zorder=4)
        ax.text(cx, 745, "abstain", fontproperties=fp("medium", 10),
                color=tk["muted"], ha="center", zorder=4)
        if i < 2:
            arrow(ax, cx + 85, 380, cx + 125, 380, tk["line"], lw=2.0, ms=12)

    # Right: contrast
    box(ax, 1170, 160, 390, 640, tk["soft"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 1365, 196, "VS OTHER ARMS", tk["muted"])

    metric_chip(ax, 1210, 250, 310, 100, "0.1248", "O0→O2 rewrite rate", tk, tk["bad"])
    metric_chip(ax, 1210, 380, 310, 100, "185", "Train abstain (stuck)", tk, tk["jev"])
    metric_chip(ax, 1210, 510, 310, 100, "1,752", "TEDDY over-answers @O2", tk, tk["bad"])
    box(ax, 1210, 640, 310, 100, tk["good_bg"], tk["good"], lw=1.6, r=12, z=2)
    ax.text(1365, 675, "Declaration only", fontproperties=fp("bold", 16),
            color=tk["good"], ha="center", va="center", zorder=3)
    ax.text(1365, 710, "YAML edit · 0 labels", fontproperties=fp("regular", 12),
            color=tk["ink"], ha="center", va="center", zorder=3)

    footer(ax, tk, "HARD_PROOF · n=16,750 · Q secondary · no retrain")
    return save(fig, tk, "01_edit_question", theme)


# ── CASE 02 · MISSING MODALITY ─────────────────────────────────────────────

def draw_02(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "CASE 2  ·  MISSING MODALITY  ·  MISSING_MODALITY",
           "Weak channel → honest silence. Don’t silently over-answer.",
           "Same site4 cells  ·  rna_only / adt_only / joint")

    masks = [
        ("rna_only", True, False, 233, 76, 0.9815, tk["rna"]),
        ("adt_only", False, True, 1664, 0, 0.8675, tk["adt"]),
        ("joint", True, True, 0, 0, 1.0000, tk["good"]),
    ]

    for i, (name, has_rna, has_adt, anm_a, teddy_a, pf, color) in enumerate(masks):
        x = 60 + i * 510
        is_bad = name == "adt_only"
        bg = tk["bad_bg"] if is_bad else tk["panel"]
        ec = tk["bad"] if is_bad else tk["line"]
        box(ax, x, 160, 480, 640, bg, ec, lw=1.6 if is_bad else 1.2, r=18, z=1)
        stage_label(ax, x + 240, 200, name, color)

        # channel glyphs
        # RNA
        rx, ry = x + 140, 320
        box(ax, rx - 60, ry - 50, 120, 100, tk["panel"], tk["rna"] if has_rna else tk["bad"], lw=2, r=12, z=2)
        if has_rna:
            draw_rna(ax, rx, ry, tk["rna"], scale=0.7)
        else:
            draw_x(ax, rx, ry, tk["bad"], size=22)
        ax.text(rx, ry + 70, "RNA", fontproperties=fp("semibold", 12),
                color=tk["rna"] if has_rna else tk["bad"], ha="center", zorder=4)

        # ADT
        ax_x, ay = x + 340, 320
        box(ax, ax_x - 60, ay - 50, 120, 100, tk["panel"], tk["adt"] if has_adt else tk["bad"], lw=2, r=12, z=2)
        if has_adt:
            for j, dx in enumerate((-28, 0, 28)):
                draw_protein(ax, ax_x + dx, ay, 14, mix(tk["adt"], tk["accent"], j * 0.3))
        else:
            draw_x(ax, ax_x, ay, tk["bad"], size=22)
        ax.text(ax_x, ay + 70, "ADT", fontproperties=fp("semibold", 12),
                color=tk["adt"] if has_adt else tk["bad"], ha="center", zorder=4)

        # comparison bars
        ax.text(x + 240, 460, "ABSTAIN", fontproperties=fp("semibold", 11),
                color=tk["muted"], ha="center", zorder=4)

        # TEDDY bar
        max_a = 1664
        tw = max(8, 360 * (teddy_a / max_a)) if teddy_a else 8
        aw = max(8, 360 * (anm_a / max_a)) if anm_a else 8
        ax.text(x + 40, 510, "TEDDY", fontproperties=fp("medium", 11), color=tk["gray"], va="center", zorder=4)
        box(ax, x + 110, 495, tw, 30, tk["gray"], r=6, z=2)
        ax.text(x + 110 + tw + 12, 510, str(teddy_a), fontproperties=fp("bold", 13),
                color=tk["bad"] if teddy_a == 0 and is_bad else tk["ink"], va="center", zorder=4)

        ax.text(x + 40, 570, "ANM", fontproperties=fp("medium", 11), color=tk["anm"], va="center", zorder=4)
        box(ax, x + 110, 555, aw, 30, tk["anm"], r=6, z=2)
        ax.text(x + 110 + aw + 12, 570, f"{anm_a:,}", fontproperties=fp("bold", 13),
                color=tk["anm"], va="center", zorder=4)

        # P_f
        box(ax, x + 80, 640, 320, 90, tk["panel"], tk["line"], lw=1, r=12, z=2)
        ax.text(x + 240, 670, f"P_f  {pf:.4f}", fontproperties=fp("bold", 18),
                color=tk["accent"], ha="center", va="center", zorder=3)
        ax.text(x + 240, 705, "ANM mean workability", fontproperties=fp("regular", 11),
                color=tk["muted"], ha="center", va="center", zorder=3)

        if is_bad:
            ax.text(x + 240, 760, "silent over-answer → honest abstain",
                    fontproperties=fp("semibold", 12), color=tk["bad"], ha="center", zorder=4)

    footer(ax, tk, "MISSING_MODALITY · O2 under adt_only abstains on all 12,563")
    return save(fig, tk, "02_missing_modality", theme)


# ── CASE 03 · ATTRIBUTION ──────────────────────────────────────────────────

def draw_03(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "CASE 3  ·  ATTRIBUTION  ·  HARD_PROOF",
           "Which marker flipped this call?",
           "Closed-form LOO on the declared field — auditable per cell")

    # Left: question glyph
    box(ax, 48, 160, 420, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 258, 200, "REVIEWER ASKS", tk["accent"])
    draw_cell(ax, 258, 340, 60, tk, label=None)
    ax.text(258, 430, "T-lineage call", fontproperties=fp("semibold", 14),
            color=tk["ink"], ha="center", zorder=4)

    # markers around cell
    markers = [("CD5", tk["bad"], True), ("CD2", tk["anm"], False), ("CD36", tk["adt"], False)]
    for i, (lab, c, drop) in enumerate(markers):
        mx = 120 + i * 140
        my = 520
        if drop:
            box(ax, mx - 40, my - 30, 80, 60, tk["bad_bg"], tk["bad"], lw=2, r=12, z=2)
            draw_x(ax, mx, my, tk["bad"], size=14)
            ax.text(mx, my + 50, f"drop {lab}?", fontproperties=fp("semibold", 12),
                    color=tk["bad"], ha="center", zorder=4)
        else:
            draw_protein(ax, mx, my, 22, c, label=lab)

    ax.text(258, 700, "Does the call flip?", fontproperties=fp("bold", 16),
            color=tk["ink"], ha="center", zorder=4)
    ax.text(258, 735, "How far was the margin?", fontproperties=fp("regular", 13),
            color=tk["muted"], ha="center", zorder=4)

    # Center: LOO flow
    box(ax, 500, 160, 520, 640, tk["accent_soft"], mix(tk["accent"], tk["line"], 0.35), lw=1.4, r=18, z=1)
    stage_label(ax, 760, 200, "ANM LOO", tk["accent"])

    # visual: panel → drop one → flip?
    steps = [
        (280, "panel", "all markers"),
        (400, "−1", "leave one out"),
        (520, "flip?", "yes / no + Δ"),
    ]
    for i, (yy, big, small) in enumerate(steps):
        box(ax, 580, yy, 360, 90, tk["panel"], tk["accent"] if i == 2 else tk["line"], lw=1.4, r=14, z=2)
        ax.text(760, yy + 35, big, fontproperties=fp("bold", 22),
                color=tk["accent"] if i == 2 else tk["ink"], ha="center", va="center", zorder=3)
        ax.text(760, yy + 65, small, fontproperties=fp("regular", 12),
                color=tk["muted"], ha="center", va="center", zorder=3)
        if i < 2:
            arrow(ax, 760, yy + 95, 760, steps[i + 1][0] - 4, tk["line"], lw=1.8, ms=12)

    ax.text(760, 700, "declared field only", fontproperties=fp("semibold", 13),
            color=tk["good"], ha="center", zorder=4)
    ax.text(760, 735, "not SHAP / not residuals", fontproperties=fp("regular", 12),
            color=tk["muted"], ha="center", zorder=4)

    # Right: metrics
    box(ax, 1050, 160, 510, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 1305, 200, "VERIFIED", tk["good"])

    metrics = [
        ("16,433", "n_attr (full-n LOO)", tk["accent"]),
        ("0.277", "top-1 flip rate", tk["anm"]),
        ("0.75", "flip-distance q50", tk["warn"]),
        ("≈0.0099", "permutation p", tk["good"]),
    ]
    for i, (val, lab, c) in enumerate(metrics):
        my = 250 + i * 120
        metric_chip(ax, 1090, my, 430, 100, val, lab, tk, c)

    footer(ax, tk, "HARD_PROOF · top-1 flips: CD5 3,882 · CD2 2,516 · CD36 2,394")
    return save(fig, tk, "03_attribution", theme)


# ── CASE 04 · SCOPE GATE ───────────────────────────────────────────────────

def draw_04(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "CASE 4  ·  SCOPE GATE  ·  SCOPE_REFINE_PROOF",
           "Spend budget on trustworthy cells.",
           "soft_P = max(action_scores) gate  ·  ANM LOO flip-sensitive scope")

    # Funnel visual
    box(ax, 48, 160, 700, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 398, 200, "COVERAGE → Q", tk["accent"])

    # funnel trapezoids approximated as stacked boxes narrowing
    widths = [560, 420, 280, 160]
    qs = ["Q 0.945", "Q ↑", "Q ↑↑", "Q 1.00"]
    covs = ["cov = 1", "raise τ", "soft_P gate", "peak"]
    colors_f = [tk["gray"], tk["anm"], tk["accent"], tk["good"]]
    for i, (ww, q, cov, c) in enumerate(zip(widths, qs, covs, colors_f)):
        yy = 260 + i * 110
        xx = 398 - ww / 2
        box(ax, xx, yy, ww, 80, mix(c, tk["panel"], 0.65), c, lw=1.6, r=12, z=2)
        ax.text(398, yy + 30, q, fontproperties=fp("bold", 16),
                color=c, ha="center", va="center", zorder=3)
        ax.text(398, yy + 55, cov, fontproperties=fp("regular", 11),
                color=tk["ink"], ha="center", va="center", zorder=3)

    # Right: metric cards
    box(ax, 780, 160, 780, 640, tk["soft"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 1170, 200, "VERIFIED ΔQ", tk["good"])

    cards = [
        ("O0 soft_P", "+0.0547", "0.9453 → 1.0000", tk["accent"]),
        ("O2 soft_P", "+0.0554", "0.9149 → 0.9703", tk["anm"]),
        ("adt_only O0", "+0.3206", "0.6288 → 0.9494", tk["warn"]),
        ("flip − random", "−0.1359", "Q 0.8099 vs 0.9458", tk["bad"]),
    ]
    for i, (title, delta, detail, c) in enumerate(cards):
        col = i % 2
        row = i // 2
        cx = 820 + col * 360
        cy = 260 + row * 240
        box(ax, cx, cy, 330, 200, tk["panel"], c, lw=1.6, r=14, z=2)
        ax.text(cx + 165, cy + 45, title, fontproperties=fp("semibold", 13),
                color=tk["muted"], ha="center", zorder=3)
        ax.text(cx + 165, cy + 100, delta, fontproperties=fp("bold", 32),
                color=c, ha="center", zorder=3)
        ax.text(cx + 165, cy + 155, detail, fontproperties=fp("regular", 12),
                color=tk["ink"], ha="center", zorder=3)

    footer(ax, tk, "SCOPE_REFINE_PROOF · n=50 hard-scope logreg ΔQ ≈ +0.034")
    return save(fig, tk, "04_scope_gate", theme)


# ── CASE 05 · LABEL COST ───────────────────────────────────────────────────

def draw_05(theme="light"):
    fig, ax, tk = new_fig(theme)
    header(ax, tk, "CASE 5  ·  LABEL COST  ·  HARD_PROOF",
           "New criterion. Zero new labels.",
           "O0→O2 mid-project  ·  numerical match ≠ criterion ownership")

    # Left: ANM at 0
    box(ax, 48, 160, 520, 640, tk["good_bg"], tk["good"], lw=2.0, r=18, z=1)
    stage_label(ax, 308, 200, "ANM DECLARATION", tk["good"])
    draw_yaml(ax, 308, 340, tk)
    ax.text(308, 430, "O2 YAML edit", fontproperties=fp("bold", 18),
            color=tk["ink"], ha="center", zorder=4)

    box(ax, 120, 480, 360, 240, tk["panel"], tk["good"], lw=1.6, r=14, z=2)
    ax.text(308, 540, "0", fontproperties=fp("bold", 64),
            color=tk["good"], ha="center", va="center", zorder=3)
    ax.text(308, 620, "endpoint labels", fontproperties=fp("semibold", 16),
            color=tk["ink"], ha="center", zorder=3)
    ax.text(308, 660, "abstain ≈0.201  ·  strict ≈0.821", fontproperties=fp("regular", 12),
            color=tk["muted"], ha="center", zorder=3)
    ax.text(308, 690, "editable field: YES", fontproperties=fp("bold", 13),
            color=tk["good"], ha="center", zorder=3)

    # Right: Train staircase
    box(ax, 600, 160, 960, 640, tk["panel"], tk["line"], lw=1.2, r=18, z=1)
    stage_label(ax, 1080, 200, "TRAIN+JEV TO MATCH NUMBERS", tk["jev"])

    steps = [
        ("grid", 50, 0.8227, 0.2063),
        ("MLP", 200, 0.8020, 0.1976),
        ("logistic", 500, 0.8129, 0.2025),
    ]
    max_n = 500
    for i, (name, n, strict, abst) in enumerate(steps):
        yy = 260 + i * 160
        bar_w = 120 + 520 * (n / max_n)
        box(ax, 660, yy, bar_w, 110, mix(tk["jev"], tk["panel"], 0.7), tk["jev"], lw=1.6, r=12, z=2)
        ax.text(680, yy + 40, name, fontproperties=fp("bold", 16),
                color=tk["ink"], va="center", zorder=3)
        ax.text(680, yy + 75, f"≈{n} O2 labels", fontproperties=fp("semibold", 13),
                color=tk["jev"], va="center", zorder=3)
        ax.text(660 + bar_w + 20, yy + 40, f"strict {strict}", fontproperties=fp("mono_med", 12),
                color=tk["muted"], va="center", zorder=3)
        ax.text(660 + bar_w + 20, yy + 70, f"abstain {abst}", fontproperties=fp("mono_med", 12),
                color=tk["muted"], va="center", zorder=3)
        ax.text(660 + bar_w + 20, yy + 95, "editable: NO", fontproperties=fp("semibold", 11),
                color=tk["bad"], va="center", zorder=3)

    ax.text(1080, 740, "Match numbers ≠ own the criterion", fontproperties=fp("semibold", 14),
            color=tk["bad"], ha="center", zorder=4)

    footer(ax, tk, "HARD_PROOF · O0-trained Train abstain stuck at 185 until retuned")
    return save(fig, tk, "05_label_cost", theme)


def main():
    drawers = [draw_00, draw_01, draw_02, draw_03, draw_04, draw_05]
    written = []
    for draw in drawers:
        for theme in ("light", "dark"):
            p = draw(theme)
            written.append(p)
            print("wrote", p, p.stat().st_size)
    print(f"done: {len(written)} PNGs → {OUT}")


if __name__ == "__main__":
    main()
