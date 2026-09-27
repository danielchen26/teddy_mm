#!/usr/bin/env python3
"""Pedagogical Care→Expectation→Problem→Solve infographics for TEDDY×ANM docs.

Style: site IBM Plex + teal tokens, anm-jev arm colours (ANM blue / Jev orange).
Numbers: only verified HARD_PROOF / SCOPE_REFINE / MISSING_MODALITY values — no invented metrics.

Usage:
  python scripts/make_story_infographics.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import to_rgb, to_hex

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "assets" / "infographics"
FONT_DIR = Path(__file__).resolve().parents[1] / "docs" / "assets" / "fonts"
for _cand in (FONT_DIR, Path("/workspace/teddy_fonts"), Path("/tmp/teddy_fonts")):
    if _cand.exists():
        FONT_DIR = _cand
        break

W, H, DPI = 1400, 820, 140

THEMES = {
    "light": dict(
        bg="#f7f6f2", panel="#ffffff", soft="#efeee8", ink="#1a1c1e", muted="#5a5f66",
        line="#d8d5cc", accent="#0d6e6e", accent_soft="#d8efef", warn="#8a4b12",
        warn_bg="#f5e6d4", good="#1f6b3a", good_bg="#dcefe3", bad="#8b2e2e",
        bad_bg="#f3e4e4", anm="#2a78d6", jev="#eb6834", gray="#8c959f",
    ),
    "dark": dict(
        bg="#121416", panel="#1b1e22", soft="#23272c", ink="#eceae4", muted="#a3a8b0",
        line="#2e333a", accent="#5ec4c4", accent_soft="#1a3333", warn="#e0a86a",
        warn_bg="#2f2418", good="#7dcea0", good_bg="#1a2e22", bad="#e08a8a",
        bad_bg="#2a1c1c", anm="#3987e5", jev="#d95926", gray="#6e7681",
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


def wrap(text, width=28):
    words = text.split()
    lines, cur = [], []
    for w in words:
        trial = " ".join(cur + [w])
        if len(trial) <= width:
            cur.append(w)
        else:
            if cur:
                lines.append(" ".join(cur))
            cur = [w]
    if cur:
        lines.append(" ".join(cur))
    return "\n".join(lines)


CASES = [
    {
        "slug": "00_global_frame",
        "eyebrow": "GLOBAL FRAME  ·  CITE / TEDDY FROZEN MAP",
        "title": "RNA → protein evidence. Editable decisions on top.",
        "subtitle": "NeurIPS 2021 BMMC CITE · site4/test n=16,750 · phase-1 full-ADT Pearson 0.610 (~0.61)",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "Trustworthy RNA→protein / cell-state readout",
                "body": "Same CITE cell: from this cell’s RNA, what is its protein panel and lineage-like state?",
                "chips": ["Biologists", "TEDDY FM scientists"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "One frozen map + one fixed score",
                "body": "Freeze TEDDY-G → z₅₁₂ → ADT preds. Once Pearson looks good enough, treat the readout as done forever.",
                "chips": ["Pearson ~0.61", "Answer every cell"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "The question is not fixed",
                "body": "Criteria shift, modalities drop, reviewers ask “why this call,” budgets force scope — a frozen head cannot edit the decision field.",
                "chips": ["Silent over-answer", "No declared LOO"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "Editable declared field on same δu",
                "body": "Typed δu events → ANM P_f / Q_f. Edit YAML observers — no TEDDY retrain, no claim Pearson > ~0.61.",
                "chips": ["Call / abstain", "Markers / scope"],
            },
        ],
        "footer": "Not clinical · not a Pearson win · holdout adt_true = verifier only · Train+Jev = log-loss stand-in",
    },
    {
        "slug": "01_edit_question",
        "eyebrow": "CASE 1  ·  EDIT THE QUESTION  ·  HARD_PROOF",
        "title": "Change the scientific question without retraining TEDDY",
        "subtitle": "Same frozen δu. O0 equal-panel → O1 tighter abstain → O2 key-marker rewrite",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "Lineage-style call from the panel",
                "body": "Is this cell B, T, or myeloid? One RNA→protein map, one scientific question — until collaborators change the criterion.",
                "chips": ["B / T / myeloid"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "One fixed decision rule forever",
                "body": "Maximize Pearson (~0.61), freeze one marker-weight rule, never touch the map again.",
                "chips": ["Fixed weights", "Same question weekly"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "O0→O2 rewrites expected_action",
                "body": "Disagree O0↔O1 = 0.0000 (tighter abstain only). O0↔O2 rewrites 1,794 / 14,379 cells (rate 0.1248). Frozen Pearson head cannot rewrite that field.",
                "chips": ["Rate 0.1248", "Needs label retune"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "YAML edit only — abstain tracks observer",
                "body": "ANM abstain 317 → 834 → 3,373 (O0→O1→O2). TEDDY alone still answers on 780 / 1,752 cells where O1 / O2 abstain. Train+Jev stuck at abstain 185 on every O*.",
                "chips": ["Declaration only", "Train stuck @185"],
            },
        ],
        "footer": "HARD_PROOF · site4/test n=16,750 · Q secondary · no retrain",
    },
    {
        "slug": "02_missing_modality",
        "eyebrow": "CASE 2  ·  MISSING MODALITY  ·  MISSING_MODALITY_ANM_DEMO",
        "title": "Abstain when RNA or ADT is weak — don’t silently over-answer",
        "subtitle": "Same site4 cells · masks rna_only / adt_only / joint · reverse diagnosis, not a Pearson contest",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "Can I trust the readout when a channel is offline?",
                "body": "Real CITE batches are messy. Biologists need honest silence when evidence is thin.",
                "chips": ["RNA weak", "ADT offline"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "Always answer from whatever is present",
                "body": "A “good” multimodal model keeps answering and keeps Pearson high under joint RNA+ADT. Silence looks like failure.",
                "chips": ["Always answer", "Joint Pearson"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "adt_only: TEDDY abstain = 0 while Q ≈ 0.63",
                "body": "Silent over-answer under weak ADT. Always-answering ≠ trustworthy. Panel Pearson is secondary; full-ADT phase-1 reference stays ~0.61.",
                "chips": ["TEDDY abstain 0", "Q ≈ 0.63"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "P_f / abstain track the mask",
                "body": "ANM O0: rna_only abstain 233 (P_f 0.9815); adt_only abstain 1,664 (P_f 0.8675); joint abstain 0 (P_f 1.0). O2 under adt_only abstains on all 12,563.",
                "chips": ["Abstain 1,664", "O2: all 12,563"],
            },
        ],
        "footer": "MISSING_MODALITY · hybrid export (phase-1 RNA + phase-2 ADT) · Train-retune plateaus Q≈0.68–0.70 @1k O2 labels",
    },
    {
        "slug": "03_attribution",
        "eyebrow": "CASE 3  ·  ATTRIBUTION  ·  HARD_PROOF",
        "title": "Which surface marker flipped this cell’s call?",
        "subtitle": "Closed-form leave-one-out on the declared field — auditable per cell",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "Scientific attribution on protein evidence",
                "body": "After a call: which marker actually flipped the decision? Not a residual plot — a yes/no flip with margin.",
                "chips": ["Why this call?"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "Saliency / SHAP / top expressed",
                "body": "Once Pearson is good enough, rank opaque coefficients. A “good” answer is a feature list — not a flip distance on an editable observer.",
                "chips": ["SHAP ranks", "Opaque coeffs"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "Need LOO on a declared decision field",
                "body": "Reviewer: “If we drop CD5, does the call flip — and how far was the margin?” Panel-mean / Train heads lack closed-form LOO flip-distance. Top-1 flips: CD5 3,882 · CD2 2,516 · CD36 2,394.",
                "chips": ["CD5 3,882", "No declared LOO"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "Full-n LOO + flip distances",
                "body": "n_attr=16,433 · top-1 flip rate 0.277 · flip-distance q50=0.75 · bootstrap B=200 (mode CD5) · permutation p≈0.0099. Matching abstain ≠ attribution.",
                "chips": ["n_attr 16,433", "p≈0.0099"],
            },
        ],
        "footer": "HARD_PROOF · declared O0 field · TEDDY alone / Train+Jev: no declared-field LOO in the proof",
    },
    {
        "slug": "04_scope_gate",
        "eyebrow": "CASE 4  ·  SCOPE GATE  ·  SCOPE_REFINE_PROOF",
        "title": "Spend limited follow-up budget on trustworthy cells",
        "subtitle": "soft_P = max(action_scores) gate + ANM LOO flip-sensitive scope",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "Which cells to confirm / label next?",
                "body": "Follow-up budget is limited. Conditional correctness under partial coverage — not only full-set average score.",
                "chips": ["Budget-aware", "Conditional Q"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "One global score at coverage = 1",
                "body": "Report one accuracy / Pearson on the full test set. Every cell equally worth answering.",
                "chips": ["Full-set Pearson", "cov = 1"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "Can confirm only ~40% of cells",
                "body": "Full-coverage Pearson hides barely-decided cells. Without a workability gate or attributed-hard scope, you cannot trade coverage for conditional correctness.",
                "chips": ["~40% confirmable", "Barely-decided hidden"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "Raise τ on soft_P · mark flip-sensitive",
                "body": "O0 soft_P: 0.9453 → peak 1.0000 (ΔQ +0.0547). O2: +0.0554. adt_only O0: +0.3206. Flip-sensitive Q 0.8099 vs random 0.9458 (ΔQ −0.1359). n=50 hard-scope logreg ΔQ ≈ +0.034.",
                "chips": ["ΔQ +0.0547", "Flip Q 0.8099"],
            },
        ],
        "footer": "SCOPE_REFINE_PROOF · TEDDY alone = cov=1 baseline · Train lacks declared soft_P field that produced the scope",
    },
    {
        "slug": "05_label_cost",
        "eyebrow": "CASE 5  ·  LABEL COST  ·  HARD_PROOF",
        "title": "Adopt a new criterion with zero new endpoint labels",
        "subtitle": "O0→O2 mid-project · numerical match ≠ criterion ownership",
        "panels": [
            {
                "k": "1 · Care about", "kind": "care",
                "headline": "New rubric without a labeling campaign",
                "body": "Scoring rubric changes O0→O2. Can biologists / FM teams adopt it without collecting a new labeled set?",
                "chips": ["O0 → O2", "Zero-label transfer"],
            },
            {
                "k": "2 · Expectation", "kind": "expect",
                "headline": "Collect labels · fit a small head",
                "body": "Cost = labeled cells. Success = matching accuracy / abstain numbers on the new observer.",
                "chips": ["Label campaign", "Calibrator fit"],
            },
            {
                "k": "3 · Where it breaks", "kind": "problem",
                "headline": "Match numbers, still no ownership",
                "body": "Even after ≈50–500 O2 labels to hit abstain≈0.20 / strict≈0.82, Train has no editable-field attribution and must re-spend on the next observer edit.",
                "chips": ["≈50–500 labels", "No editable field"],
            },
            {
                "k": "4 · TEDDY + ANM", "kind": "solve",
                "headline": "ANM already at the O2 operating point",
                "body": "ANM: 0 endpoint labels · abstain≈0.201 · strict≈0.821 · editable field. To match: grid≈50 · MLP≈200 · logistic≈500. O0-trained Train abstain stays 185 across O0/O1/O2 until retuned.",
                "chips": ["ANM: 0 labels", "Editable: Yes"],
            },
        ],
        "footer": "HARD_PROOF · Train can numerically match — never buys declared field or LOO flip distances",
    },
]


def panel_colors(tk, kind):
    if kind == "care":
        return tk["accent_soft"], tk["accent"], tk["accent"]
    if kind == "expect":
        return tk["soft"], tk["line"], tk["muted"]
    if kind == "problem":
        return tk["bad_bg"], mix(tk["bad"], tk["line"], 0.35), tk["bad"]
    if kind == "solve":
        return tk["good_bg"], mix(tk["good"], tk["line"], 0.35), tk["good"]
    return tk["panel"], tk["line"], tk["ink"]


def draw_case(case, theme="light"):
    tk = THEMES[theme]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=tk["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ax.set_facecolor(tk["bg"])

    box(ax, 28, 22, W - 56, 118, tk["panel"], tk["line"], lw=1.0, r=16, z=1)
    box(ax, 28, 22, 8, 118, tk["accent"], r=4, z=2)

    ax.text(56, 48, case["eyebrow"], fontproperties=fp("semibold", 11), color=tk["accent"], va="center", zorder=3)
    ax.text(56, 82, case["title"], fontproperties=fp("bold", 22), color=tk["ink"], va="center", zorder=3)
    ax.text(56, 112, case["subtitle"], fontproperties=fp("regular", 12.5), color=tk["muted"], va="center", zorder=3)

    n = 4
    gap = 16
    left = 28
    usable = W - 56 - gap * (n - 1)
    pw = usable / n
    top = 158
    ph = 560

    for i, p in enumerate(case["panels"]):
        x = left + i * (pw + gap)
        bg, ec, accent = panel_colors(tk, p["kind"])
        box(ax, x, top, pw, ph, bg, ec, lw=1.2, r=16, z=1)
        box(ax, x, top, pw, 6, accent, r=3, z=2)

        ax.text(x + 18, top + 36, p["k"].upper(), fontproperties=fp("semibold", 10.5), color=accent, va="center", zorder=3)
        ax.text(
            x + 18, top + 78, wrap(p["headline"], width=22 if pw < 340 else 26),
            fontproperties=fp("bold", 15.5), color=tk["ink"], va="top", linespacing=1.25, zorder=3,
        )
        ax.text(
            x + 18, top + 175, wrap(p["body"], width=28 if pw < 340 else 32),
            fontproperties=fp("regular", 12.2),
            color=tk["muted"] if p["kind"] == "care" else tk["ink"],
            va="top", linespacing=1.4, zorder=3,
        )

        cy = top + ph - 70
        cx = x + 18
        for chip in p.get("chips", []):
            tw = max(70, len(chip) * 7.2 + 22)
            if cx + tw > x + pw - 14:
                cy += 32
                cx = x + 18
            chip_bg = mix(accent, bg, 0.78) if p["kind"] != "expect" else tk["panel"]
            chip_ec = mix(accent, tk["line"], 0.5)
            box(ax, cx, cy, tw, 26, chip_bg, chip_ec, lw=0.8, r=13, z=2)
            ax.text(
                cx + tw / 2, cy + 13, chip, fontproperties=fp("medium", 10),
                color=tk["ink"], ha="center", va="center", zorder=3,
            )
            cx += tw + 8

        if i < n - 1:
            ax.annotate(
                "",
                xy=(x + pw + gap * 0.55, top + ph / 2),
                xytext=(x + pw + gap * 0.15, top + ph / 2),
                arrowprops=dict(arrowstyle="-|>", color=tk["line"], lw=1.4, mutation_scale=12),
                zorder=4,
            )

    ax.text(28, H - 28, case["footer"], fontproperties=fp("regular", 11), color=tk["muted"], va="center", zorder=3)
    ax.text(
        W - 28, H - 28, "TEDDY × ANM", fontproperties=fp("semibold", 11),
        color=tk["accent"], ha="right", va="center", zorder=3,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{case['slug']}_{theme}.png"
    fig.savefig(path, dpi=DPI, facecolor=tk["bg"], bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return path


def main():
    written = []
    for case in CASES:
        for theme in ("light", "dark"):
            p = draw_case(case, theme)
            written.append(p)
            print("wrote", p, p.stat().st_size)
    print(f"done: {len(written)} PNGs → {OUT}")


if __name__ == "__main__":
    main()
