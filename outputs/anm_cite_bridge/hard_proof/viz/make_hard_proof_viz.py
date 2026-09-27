#!/usr/bin/env python3
"""HARD_PROOF live comparison visuals for TEDDY×ANM bakeoff.

Style: import Canvas / Clip / themes from anm-jev/anim/style.py (same 1280×720
README GIF pipeline, tokens, chips, headers). Numbers: hard_proof_results.json
+ attr_compact.npz only — no invented rates.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[4]  # teddy_mm
PROOF = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
ANM_JEV = Path("/Users/tianchichen/Documents/GitHub/anm-jev")
sys.path.insert(0, str(ANM_JEV))

from anim.style import (  # noqa: E402
    THEMES, W, H, DPI, SCALE, TEXT_SCALE, Canvas, Clip, ease, ease_out, mix, pt, seg, write_gif,
)

JSON_PATH = PROOF / "hard_proof_results.json"
NPZ_PATH = PROOF / "attr_compact.npz"


def load_proof():
    with open(JSON_PATH) as f:
        return json.load(f)


# ── colours for three arms (anm-jev series + TEDDY as gray/ink accent) ──────
def arm_colors(tk):
    return {
        "teddy": tk["gray"],          # multimodal alone
        "anm": tk["anm"],             # declared field
        "jev": tk["jev"],             # logistic / Jev stand-in
        "third": tk["third"],
    }


# =============================================================================
# 1. Observer shift GIF — O0→O1→O2, three panels
# =============================================================================
class ObserverShift(Clip):
    name = "observer_shift_o0o1o2"
    duration = 10.0
    hold = 2.8
    fps = 20
    size = (W, H)

    def build(self):
        c, tk = self.c, self.tk
        d = load_proof()
        teddy = d["arms"]["teddy_alone"]
        anm = d["arms"]["anm"]
        train = d["arms"]["train_plus_anm"]
        flat = train["zero_shot_o0_model_abstain_does_not_track_observer"]
        silent = teddy["silent_over_answer"]
        n = d["n_site4"]

        self.crits = ["O0", "O1", "O2"]
        self.abstain = {
            "teddy": [teddy["criteria"][k]["abstain_count"] for k in self.crits],
            "anm": [anm["criteria"][k]["abstain_count"] for k in self.crits],
            "jev": [flat["abstain_counts"][k] for k in self.crits],
        }
        self.rates = {
            "teddy": [teddy["criteria"][k]["abstain_rate"] for k in self.crits],
            "anm": [anm["criteria"][k]["abstain_rate"] for k in self.crits],
            "jev": [flat["abstain_counts"][k] / n for k in self.crits],
        }
        self.Q = {
            "teddy": [teddy["criteria"][k]["Q_analogue"] for k in self.crits],
            "anm": [anm["criteria"][k]["Q_analogue"] for k in self.crits],
            "jev": [flat["Q"][k] for k in self.crits],
        }
        self.silent = silent
        self.n = n
        self.max_abs = max(max(self.abstain["teddy"]), max(self.abstain["anm"]), max(self.abstain["jev"]))

        col = arm_colors(tk)
        c.header(
            "HARD PROOF  ·  OBSERVER SHIFT  O0 → O1 → O2",
            "Same δu (site4/test). Who follows the declared observer?",
            f"n={n:,} cells  ·  abstain rises with declaration for ANM; TEDDY silently over-answers; "
            "Jev stand-in stays flat.",
        )

        # criterion strip
        c.box(40, 118, 1200, 44, r=10, fc=tk["panel"], ec=tk["edge"], lw=1)
        c.text(56, 146, "criterion", size=14, color=tk["muted"], weight="bold")
        self.crit_labels = []
        xs = [160, 420, 760]
        labels = [
            "O0  equal-panel",
            "O1  tighter abstain",
            "O2  key-marker rewrite",
        ]
        for i, (x, lab) in enumerate(zip(xs, labels)):
            t = c.text(x, 146, lab, size=17, color=tk["ink"], weight="bold" if i == 0 else "normal", alpha=1.0)
            self.crit_labels.append(t)
        self.arrow1 = c.text(370, 146, "→", size=20, color=tk["muted"], alpha=0.35)
        self.arrow2 = c.text(710, 146, "→", size=20, color=tk["muted"], alpha=0.35)
        disagree = d["disagreement_rates"]["O0_vs_O2"]["disagree_rate_among_both"]
        self.rewrite_chip = c.chip(
            1040, 127, f"O0↔O2 rewrite {disagree:.1%}",
            fc=mix(tk["anm"], tk["panel"], 0.82), color=tk["ink"], size=13,
        )
        for a in self.rewrite_chip[:2]:
            a.set_alpha(0)

        # three panels
        panels = [
            (40, col["teddy"], "TEDDY alone", "multimodal readout  ·  no declared field"),
            (440, col["anm"], "TEDDY + ANM", "declared field  ·  0 endpoint-label fit"),
            (840, col["jev"], "TEDDY + Jev stand-in", "O0-trained logistic  ·  one conf gate"),
        ]
        self.panels = []
        for x, color, name, desc in panels:
            c.box(x, 178, 380, 480, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
            c.system_header(x, 178, 380, name, desc, color)
            P = {"x": x, "col": color}
            P["big"] = c.text(x + 20, 290, "", size=48, weight="bold", color=tk["ink"])
            P["big_unit"] = c.text(x + 20, 318, "abstain count", size=14, color=tk["ink2"])
            P["rate"] = c.text(x + 360, 290, "", size=18, weight="bold", ha="right", color=tk["ink2"])
            P["q"] = c.text(x + 20, 360, "", size=16, color=tk["ink2"])
            # bar track for abstain fraction of max
            c.box(x + 20, 390, 340, 16, r=8, fc=tk["code_bg"], z=2)
            P["bar"] = c.box(x + 20, 390, 8, 16, r=8, fc=color, z=3)
            P["bar"].set_visible(False)
            P["note"] = c.text(x + 20, 440, "", size=14, color=tk["ink2"], linespacing=1.35)
            P["chip"] = None
            self.panels.append(P)

        # silent-over-answer annotation area on TEDDY panel (built as text)
        self.silent_txt = c.text(60, 520, "", size=13.5, color=tk["ink2"], linespacing=1.3)
        self.silent_chip_artists = []

        c.footer(
            f"HARD_PROOF.json  ·  abstain counts on site4/test n={n:,}  ·  "
            "Jev stand-in = zero_shot_o0_model_abstain_does_not_track_observer  ·  no retrain"
        )
        self.phase_chip = c.chip(1080, 26, "O0", fc=mix(tk["anm"], tk["panel"], 0.75), color=tk["ink"], size=14)
        for a in self.phase_chip[:2]:
            a.set_alpha(1)

    def _phase(self, t):
        # 0–2.4 O0, 2.4–5.4 O1, 5.4–8.8 O2, then hold story
        if t < 2.4:
            return 0, ease(seg(t, 0.4, 1.2))
        if t < 5.4:
            return 1, ease(seg(t, 2.4, 3.2))
        return 2, ease(seg(t, 5.4, 6.2))

    def update(self, t):
        tk = self.tk
        phase, fade = self._phase(t)
        # criterion highlight
        for i, lab in enumerate(self.crit_labels):
            if i == phase:
                lab.set_alpha(1.0)
                lab.set_fontweight("bold")
                lab.set_color(tk["ink"])
            elif i < phase:
                lab.set_alpha(0.45)
                lab.set_fontweight("normal")
                lab.set_color(tk["ink2"])
            else:
                lab.set_alpha(0.30)
                lab.set_fontweight("normal")
                lab.set_color(tk["muted"])
        self.arrow1.set_alpha(0.9 if phase >= 1 else 0.25)
        self.arrow2.set_alpha(0.9 if phase >= 2 else 0.25)
        for a in self.rewrite_chip[:2]:
            a.set_alpha(ease(seg(t, 5.6, 6.2)) if phase >= 2 else 0)

        # update phase chip text by rebuilding label string
        self.phase_chip[1].set_text(["O0", "O1", "O2"][phase])

        keys = ["teddy", "anm", "jev"]
        notes = [
            (
                "silent over-answer vs ANM declared abstain"
                if phase > 0 else
                "low abstain under O0 rule"
            ),
            (
                "abstain tracks declaration\n0 labels · edit YAML only"
                if phase > 0 else
                "declared field ready"
            ),
            (
                "flat abstain — does not\ntrack observer (O0 gate)"
            ),
        ]
        for i, (P, k) in enumerate(zip(self.panels, keys)):
            cnt = self.abstain[k][phase]
            rate = self.rates[k][phase]
            q = self.Q[k][phase]
            # interpolate from previous phase for motion
            if phase == 0:
                cnt_s, rate_s = cnt, rate
            else:
                prev_c = self.abstain[k][phase - 1]
                prev_r = self.rates[k][phase - 1]
                cnt_s = prev_c + (cnt - prev_c) * fade
                rate_s = prev_r + (rate - prev_r) * fade
            P["big"].set_text(f"{int(round(cnt_s)):,}")
            P["big"].set_alpha(0.35 + 0.65 * fade if phase > 0 else 1.0)
            P["rate"].set_text(f"{rate_s:.1%}")
            P["q"].set_text(f"Q_analogue  {q:.3f}")
            frac = cnt_s / max(self.max_abs, 1)
            P["bar"].set_visible(True)
            P["bar"].set_width(max(8, 340 * frac))
            P["note"].set_text(notes[i])

        # silent over-answer numbers (real JSON)
        if phase == 0:
            self.silent_txt.set_text("")
        elif phase == 1:
            s = self.silent["O0_over_answer_vs_O1_decl"]
            self.silent_txt.set_text(
                f"silent over-answer vs O1 declared:\n"
                f"{s['count']:,} cells  ({s['rate']:.2%})"
            )
            self.silent_txt.set_alpha(ease(seg(t, 3.0, 3.6)))
        else:
            s = self.silent["O0_over_answer_vs_O2_decl"]
            self.silent_txt.set_text(
                f"silent over-answer vs O2 declared:\n"
                f"{s['count']:,} cells  ({s['rate']:.2%})"
            )
            self.silent_txt.set_alpha(ease(seg(t, 6.0, 6.6)))


# =============================================================================
# 2. Label-cost curve — static + GIF
# =============================================================================
class LabelCost(Clip):
    name = "label_cost_curve"
    duration = 8.5
    hold = 2.5
    fps = 20
    size = (W, H)

    def build(self):
        c, tk = self.c, self.tk
        d = load_proof()
        tpa = d["arms"]["train_plus_anm"]
        targets = tpa["anm_o2_targets"]
        match = tpa["match_headline"]
        log_curve = tpa["curves"]["logistic_jev_standin"]
        mlp_curve = tpa["curves"]["mlp_small_head"]
        grid_curve = tpa["curves"]["teddy_thr_boost_grid"]

        def cal_rows(curve, key="abstain_calibrated_to_anm_o2"):
            rows = []
            for row in curve:
                if "skipped" in row:
                    continue
                cal = row.get(key) or row.get("abstain_calibrated_thr")
                if not cal:
                    continue
                rows.append((row["n_o2_labels"], cal["abstain_rate"], cal["Q_analogue"],
                             cal.get("accuracy_strict_labeled"), cal.get("Q_gap_vs_anm")))
            return rows

        self.log_rows = cal_rows(log_curve)
        self.mlp_rows = cal_rows(mlp_curve)
        self.grid_rows = cal_rows(grid_curve)
        self.targets = targets
        self.match = match

        col = arm_colors(tk)
        c.header(
            "HARD PROOF  ·  LABEL COST TO MATCH ANM O2",
            "Train+ANM needs N O2 labels. ANM needs 0.",
            "Target = ANM declared O2 operating point on the same eval "
            f"(abstain={targets['abstain_rate']:.4f}, Q={targets['Q_analogue']:.4f}).",
        )

        # left story cards
        c.box(40, 124, 380, 520, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
        c.system_header(40, 124, 380, "ANM declared O2", "edit YAML · 0 endpoint labels", col["anm"])
        c.text(58, 250, "0", size=72, weight="bold", color=tk["ink"])
        c.text(58, 290, "O2 labels consumed", size=16, color=tk["ink2"])
        c.text(58, 340, f"abstain  {targets['abstain_rate']:.4f}", size=18, weight="bold")
        c.text(58, 368, f"Q_analogue  {targets['Q_analogue']:.4f}", size=18, weight="bold")
        c.text(58, 396, f"strict  {targets['accuracy_strict_labeled']:.4f}", size=16, color=tk["ink2"])
        c.text(58, 450, "editable field · LOO attribution\npermission-free observer switch", size=14,
               color=tk["ink2"], linespacing=1.4)

        # match headlines
        y0 = 510
        for label, key, color in (
            ("grid thr/boost", "grid_calibrated", col["third"]),
            ("MLP head", "mlp_calibrated", col["teddy"]),
            ("logistic (Jev)", "logistic_calibrated", col["jev"]),
        ):
            m = match[key]
            c.box(58, y0, 12, 12, r=6, fc=color, z=4)
            c.text(80, y0 + 10, f"{label}: N={m['n_o2_labels']}  strict={m['accuracy_strict_labeled']:.3f}",
                   size=13.5, color=tk["ink"], va="center_baseline")
            y0 += 28

        # chart
        CH_X, CH_Y, CH_W, CH_H = 460, 160, 760, 420
        a = self.ax = c.chart(CH_X, CH_Y, CH_W, CH_H)
        xs = [r[0] for r in self.log_rows]
        a.set_xscale("log")
        a.set_xlim(40, 12000)
        a.set_ylim(0.92, 0.98)
        a.set_yticks([0.93, 0.94, 0.95, 0.96, 0.97, 0.98])
        a.axhline(targets["Q_analogue"], color=col["anm"], lw=2.0, ls="--", alpha=0.85, zorder=3)
        self.anm_hline_lab = c.text(CH_X + CH_W - 8, CH_Y + 18,
                                    f"ANM O2 Q={targets['Q_analogue']:.3f}  @ 0 labels",
                                    size=13.5, color=col["anm"], weight="bold", ha="right")

        self.log_line, = a.plot([], [], color=col["jev"], lw=2.6, zorder=4)
        self.log_pts = a.scatter([], [], s=70, color=col["jev"], edgecolor=tk["bg"], linewidth=1.8, zorder=5)
        self.mlp_line, = a.plot([], [], color=col["teddy"], lw=2.2, zorder=4)
        self.mlp_pts = a.scatter([], [], s=55, color=col["teddy"], edgecolor=tk["bg"], linewidth=1.5, zorder=5)
        self.grid_line, = a.plot([], [], color=col["third"], lw=2.0, zorder=4, ls=":")
        self.grid_pts = a.scatter([], [], s=55, color=col["third"], edgecolor=tk["bg"], linewidth=1.5, zorder=5, marker="D")

        # match markers (stars) — positions fixed, alpha animated
        self.match_marks = []
        for key, color, marker in (
            ("logistic_calibrated", col["jev"], "*"),
            ("mlp_calibrated", col["teddy"], "*"),
            ("grid_calibrated", col["third"], "*"),
        ):
            m = match[key]
            sc = a.scatter([m["n_o2_labels"]], [m["Q_analogue"]], s=220, marker=marker,
                           color=color, edgecolor=tk["bg"], linewidth=1.5, zorder=7, alpha=0)
            self.match_marks.append(sc)

        c.text(CH_X - 8, CH_Y + 10, "Q_analogue\n(calibrated abstain≈ANM O2)", size=13, color=tk["muted"], ha="right")
        c.text(CH_X + CH_W / 2, 650, "N O2 labels given to Train+ANM (log scale)", size=14, color=tk["muted"], ha="center")
        c.text(40, 690,
               "match = abstain+accuracy_strict vs ANM O2 declared  ·  "
               f"logistic N={match['logistic_calibrated']['n_o2_labels']}, "
               f"mlp N={match['mlp_calibrated']['n_o2_labels']}, "
               f"grid N={match['grid_calibrated']['n_o2_labels']}  ·  hard_proof_results.json",
               size=12.5, color=tk["muted"])

        self.n_reveal = c.text(1220, 200, "", size=22, weight="bold", ha="right", family="monospace")

    def update(self, t):
        # progressively reveal logistic then mlp then grid
        n_log = len(self.log_rows)
        p = ease(seg(t, 0.5, 5.5))
        k = max(1, int(round(1 + (n_log - 1) * p)))
        xs = [r[0] for r in self.log_rows[:k]]
        ys = [r[2] for r in self.log_rows[:k]]
        self.log_line.set_data(xs, ys)
        self.log_pts.set_offsets(np.c_[xs, ys])
        self.n_reveal.set_text(f"N ≤ {xs[-1]:,}")

        p2 = ease(seg(t, 3.0, 6.0))
        k2 = max(1, int(round(1 + (len(self.mlp_rows) - 1) * p2))) if self.mlp_rows else 0
        if k2:
            xs2 = [r[0] for r in self.mlp_rows[:k2]]
            ys2 = [r[2] for r in self.mlp_rows[:k2]]
            self.mlp_line.set_data(xs2, ys2)
            self.mlp_pts.set_offsets(np.c_[xs2, ys2])

        p3 = ease(seg(t, 4.0, 6.5))
        k3 = max(1, int(round(1 + (len(self.grid_rows) - 1) * p3))) if self.grid_rows else 0
        if k3:
            xs3 = [r[0] for r in self.grid_rows[:k3]]
            ys3 = [r[2] for r in self.grid_rows[:k3]]
            self.grid_line.set_data(xs3, ys3)
            self.grid_pts.set_offsets(np.c_[xs3, ys3])

        s = ease(seg(t, 6.2, 7.0))
        for sc in self.match_marks:
            sc.set_alpha(s)


def render_label_cost_static(theme="light"):
    """Final-frame static PNG of the label-cost story."""
    clip = LabelCost(theme)
    clip.update(clip.duration - 0.05)
    img = clip.c.frame_rgb()
    out = OUT / f"label_cost_curve_{theme}.png"
    img.save(out)
    plt.close(clip.c.fig)
    return out


# =============================================================================
# 3. Attribution — top-1 protein + flip-distance
# =============================================================================
def render_attribution(theme="light"):
    d = load_proof()
    attr = d["arms"]["anm"]["attribution"]
    z = np.load(NPZ_PATH, allow_pickle=True)
    top = z["top_protein"].astype(str)
    flipped = z["flipped"].astype(bool)
    flip_d = z["flip_distance"].astype(float)

    c = Canvas(theme)
    tk = c.tk
    col = arm_colors(tk)
    n_attr = int(attr["n_attr_cells"])
    dist = attr["top1_protein_distribution"]
    boot = attr["bootstrap_top1"]
    perm = attr["permutation_test"]
    fq = attr["flip_distance_quantiles"]

    c.header(
        "HARD PROOF  ·  ATTRIBUTION  (ANM O0 FIELD)",
        "Top-1 protein + flip-distance  ·  full-n LOO",
        f"n_attr={n_attr:,} / requested={attr['n_attr_requested']:,}  ·  "
        f"bootstrap B={boot['n_boot']}  ·  permutation n_perm={perm['n_perm']}  ·  "
        f"p≈{perm['p_value_one_sided']:.3f}",
    )

    # left: top-1 histogram
    c.box(40, 124, 600, 520, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
    c.system_header(40, 124, 600, "Top-1 protein", "closed-form field LOO  ·  mode CD5", col["anm"])

    proteins = sorted(dist.items(), key=lambda kv: -kv[1])
    names = [p for p, _ in proteins]
    counts = [v for _, v in proteins]
    CH_X, CH_Y, CH_W, CH_H = 70, 220, 540, 320
    a = c.chart(CH_X, CH_Y, CH_W, CH_H)
    ypos = np.arange(len(names))
    a.barh(ypos, counts, color=col["anm"], height=0.72, edgecolor="none", zorder=3)
    a.set_yticks(ypos)
    a.set_yticklabels(names)
    a.set_xlim(0, max(counts) * 1.12)
    a.invert_yaxis()
    a.grid(True, axis="x", color=tk["grid"], lw=1.0)
    a.spines["bottom"].set_visible(True)
    # annotate mode
    mode = boot["observed_mode_protein"]
    mode_frac = boot["observed_mode_fraction"]
    ci = boot["boot_mode_fraction_ci95"]
    c.text(58, 560,
           f"mode {mode}  frac={mode_frac:.4f}  boot CI95=[{ci[0]:.4f}, {ci[1]:.4f}]",
           size=13.5, color=tk["ink2"])
    c.text(58, 585,
           f"perm null_max_mean={perm['null_max_fraction_mean']:.4f}  "
           f"p={perm['p_value_one_sided']:.4f}  (n_perm={perm['n_perm']})",
           size=13.5, color=tk["ink2"])

    # right: flip distance
    c.box(660, 124, 580, 520, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
    c.system_header(660, 124, 580, "Flip distance", "among cells with a flip", col["jev"])

    with_flip = flip_d[flipped & np.isfinite(flip_d) & (flip_d >= 0)]  # exclude sentinel -1
    CH2_X, CH2_Y, CH2_W, CH2_H = 690, 220, 520, 280
    a2 = c.chart(CH2_X, CH2_Y, CH2_W, CH2_H)
    if len(with_flip):
        a2.hist(with_flip, bins=30, color=mix(col["jev"], tk["panel"], 0.15),
                edgecolor=col["jev"], linewidth=1.2, zorder=3)
    a2.axvline(fq["q25"], color=tk["ink2"], lw=1.4, ls="--", alpha=0.8)
    a2.axvline(fq["q50"], color=tk["ink"], lw=1.8, ls="-")
    a2.axvline(fq["q75"], color=tk["ink2"], lw=1.4, ls="--", alpha=0.8)
    a2.set_xlabel("")
    c.text(CH2_X + CH2_W / 2, 530, "flip distance (cells with flip)", size=14, color=tk["muted"], ha="center")

    flip_rate = attr["top1_flip_rate"]
    c.text(680, 560,
           f"flip rate={flip_rate:.4f}  ·  n_with_flip={fq['n_with_flip']:,}  ·  "
           f"n_no_flip={fq['n_no_flip']:,}",
           size=13.5, color=tk["ink2"])
    c.text(680, 585,
           f"quantiles  q25={fq['q25']:.3f}  q50={fq['q50']:.3f}  q75={fq['q75']:.3f}",
           size=13.5, color=tk["ink2"])
    # chip for p
    c.chip(1000, 140, f"perm p≈{perm['p_value_one_sided']:.3f}",
           fc=mix(col["anm"], tk["panel"], 0.78), color=tk["ink"], size=13)

    c.footer("attr_compact.npz + hard_proof_results.json attribution block  ·  ANM O0 field  ·  no retrain")
    img = c.frame_rgb()
    out = OUT / f"attribution_top1_flip_{theme}.png"
    img.save(out)
    plt.close(c.fig)
    return out


# =============================================================================
# 4. Summary dashboard
# =============================================================================
def render_dashboard(theme="light"):
    d = load_proof()
    teddy = d["arms"]["teddy_alone"]
    anm = d["arms"]["anm"]
    train = d["arms"]["train_plus_anm"]
    flat = train["zero_shot_o0_model_abstain_does_not_track_observer"]
    match = train["match_headline"]
    targets = train["anm_o2_targets"]
    attr = anm["attribution"]
    silent = teddy["silent_over_answer"]
    disagree = d["disagreement_rates"]["O0_vs_O2"]["disagree_rate_among_both"]
    n = d["n_site4"]

    # slightly taller canvas for dashboard
    c = Canvas(theme, w=W, h=H)
    tk = c.tk
    col = arm_colors(tk)

    c.header(
        "HARD PROOF  ·  TEDDY × ANM × JEV STAND-IN",
        "Declaration edits the observer. Labels buy numbers, not the field.",
        f"site4/test n={n:,}  ·  OOD n={d['n_ood']:,}  ·  n_attr={attr['n_attr_cells']:,}  ·  "
        f"perm p={attr['permutation_test']['p_value_one_sided']:.4f}  ·  {d['timestamp_local']}",
    )

    # row of three abstain sparklines (static mini bars O0/O1/O2)
    arms_spec = [
        (40, col["teddy"], "TEDDY alone",
         [teddy["criteria"][k]["abstain_count"] for k in ("O0", "O1", "O2")],
         f"silent OA vs O2: {silent['O0_over_answer_vs_O2_decl']['count']:,} "
         f"({silent['O0_over_answer_vs_O2_decl']['rate']:.1%})"),
        (440, col["anm"], "TEDDY + ANM",
         [anm["criteria"][k]["abstain_count"] for k in ("O0", "O1", "O2")],
         f"O0→O2 rewrite {disagree:.1%}  ·  0 labels"),
        (840, col["jev"], "Jev stand-in",
         [flat["abstain_counts"][k] for k in ("O0", "O1", "O2")],
         "flat 185 / 185 / 185 — no track"),
    ]
    max_a = max(max(a[3]) for a in arms_spec)
    for x, color, name, counts, note in arms_spec:
        c.box(x, 118, 380, 210, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
        c.box(x + 18, 138, 14, 14, r=7, fc=color, z=4)
        c.text(x + 42, 151, name, size=18, weight="bold", va="baseline")
        # mini bars
        for i, (lab, cnt) in enumerate(zip(("O0", "O1", "O2"), counts)):
            bx = x + 30 + i * 110
            h = 8 + 90 * (cnt / max_a)
            c.box(bx, 280 - h, 70, h, r=6, fc=color, z=3)
            c.text(bx + 35, 290, lab, size=13, color=tk["muted"], ha="center")
            c.text(bx + 35, 280 - h - 8, f"{cnt:,}", size=13, weight="bold", ha="center", va="baseline")
        c.text(x + 20, 310, note, size=12.5, color=tk["ink2"])

    # label cost strip
    c.box(40, 348, 600, 300, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
    c.system_header(40, 348, 600, "Label cost → ANM O2 point",
                    f"target abstain={targets['abstain_rate']:.4f}  Q={targets['Q_analogue']:.4f}",
                    col["anm"])
    c.text(60, 440, "ANM", size=16, weight="bold", color=col["anm"])
    c.text(200, 440, "0 labels", size=28, weight="bold")
    c.text(60, 490, "grid / MLP / logistic match at", size=14, color=tk["ink2"])
    c.text(60, 530,
           f"N = {match['grid_calibrated']['n_o2_labels']}  /  "
           f"{match['mlp_calibrated']['n_o2_labels']}  /  "
           f"{match['logistic_calibrated']['n_o2_labels']}",
           size=26, weight="bold")
    c.text(60, 575,
           f"strict @ match  "
           f"{match['grid_calibrated']['accuracy_strict_labeled']:.3f}  /  "
           f"{match['mlp_calibrated']['accuracy_strict_labeled']:.3f}  /  "
           f"{match['logistic_calibrated']['accuracy_strict_labeled']:.3f}",
           size=14, color=tk["ink2"])
    c.text(60, 610,
           "Train still lacks editable field + LOO flip-distance attribution.",
           size=13.5, color=tk["ink2"])

    # attribution strip
    c.box(660, 348, 580, 300, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
    c.system_header(660, 348, 580, "Attribution audit",
                    "full-n LOO + bootstrap + permutation", col["anm"])
    boot = attr["bootstrap_top1"]
    perm = attr["permutation_test"]
    fq = attr["flip_distance_quantiles"]
    c.text(680, 440, f"top-1 mode  {boot['observed_mode_protein']}", size=22, weight="bold")
    c.text(680, 475, f"frac={boot['observed_mode_fraction']:.4f}  "
                     f"CI95=[{boot['boot_mode_fraction_ci95'][0]:.3f}, "
                     f"{boot['boot_mode_fraction_ci95'][1]:.3f}]",
           size=14, color=tk["ink2"])
    c.text(680, 520, f"perm p = {perm['p_value_one_sided']:.4f}", size=26, weight="bold", color=col["anm"])
    c.text(680, 555,
           f"flip rate={attr['top1_flip_rate']:.3f}  ·  "
           f"q50={fq['q50']:.2f}  ·  n_flip={fq['n_with_flip']:,}",
           size=14, color=tk["ink2"])
    c.text(680, 590, "Accuracy is secondary. Not clinical.", size=14, color=tk["muted"], weight="bold")

    c.footer("Numbers from HARD_PROOF.md / hard_proof_results.json / attr_compact.npz  ·  no retrain of best.pt")
    img = c.frame_rgb()
    out = OUT / f"dashboard_hard_proof_{theme}.png"
    img.save(out)
    plt.close(c.fig)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    theme = "light"
    print("rendering observer_shift GIF…")
    ObserverShift(theme).render(str(OUT / f"observer_shift_o0o1o2_{theme}.gif"))
    print("rendering label_cost GIF + PNG…")
    LabelCost(theme).render(str(OUT / f"label_cost_curve_{theme}.gif"))
    render_label_cost_static(theme)
    print("rendering attribution PNG…")
    render_attribution(theme)
    print("rendering dashboard PNG…")
    render_dashboard(theme)
    # also dark theme statics for README dark mode (optional)
    print("rendering dark statics…")
    render_label_cost_static("dark")
    render_attribution("dark")
    render_dashboard("dark")
    ObserverShift("dark").render(str(OUT / "observer_shift_o0o1o2_dark.gif"))
    print("done →", OUT)
    for p in sorted(OUT.glob("*")):
        if p.name.startswith("_"):
            continue
        if p.is_file():
            print(f"  {p.name:48s} {p.stat().st_size/1024:8.1f} KB")


if __name__ == "__main__":
    main()
