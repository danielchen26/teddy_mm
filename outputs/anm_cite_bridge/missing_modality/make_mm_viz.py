#!/usr/bin/env python3
"""anm-jev-style GIF: three modality masks × ANM P_f / abstain vs TEDDY-alone."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
ANM_JEV = Path("/Users/tianchichen/Documents/GitHub/anm-jev")
sys.path.insert(0, str(ANM_JEV))

from anim.style import THEMES, W, H, Canvas, Clip, ease, mix, seg  # noqa: E402

JSON_PATH = OUT / "missing_modality_results.json"
MASKS = ["rna_only", "adt_only", "joint"]


def load():
    return json.loads(JSON_PATH.read_text())


class MissingModalityShift(Clip):
    name = "missing_modality_masks"
    duration = 8.0
    hold = 2.2
    fps = 16
    size = (W, H)

    def build(self):
        c, tk = self.c, self.tk
        d = load()
        self.n = d["n_base_cells"]
        self.source = d.get("source") or d.get("export_manifest", {}).get("source", "?")
        self.abstain = []
        self.teddy_abs = []
        self.pf = []
        self.q = []
        for m in MASKS:
            a0 = d["by_mask"][m]["anm"]["criteria"]["O0"]
            t0 = d["by_mask"][m]["teddy_alone"]["criteria"]["O0"]
            self.abstain.append(int(a0["abstain_count"]))
            self.teddy_abs.append(int(t0["abstain_count"]))
            self.pf.append(float(a0.get("mean_P_f") or 0.0))
            self.q.append(float(a0.get("Q_analogue") or 0.0))
        self.max_abs = max(max(self.abstain), 1)

        c.header(
            "MISSING MODALITY × ANM",
            "Same δu cells · modality_mask as source ablation",
            f"n={self.n:,} site4/test  ·  source={self.source}  ·  declared RNA reliability gate on adt_only",
        )

        notes = [
            "RNA present\nphase-1 δu (strong)",
            "RNA ABSENT\nreliability ×0.35",
            "joint fuse\nhybrid channels",
        ]
        self.panels = []
        xs = [40, 440, 840]
        for i, (x, m) in enumerate(zip(xs, MASKS)):
            c.box(x, 150, 380, 500, r=14, fc=tk["panel"], ec=tk["edge"], lw=1)
            c.system_header(x, 150, 380, m, notes[i].replace("\n", " · "), tk["anm"])
            P = {"x": x}
            P["pf_label"] = c.text(x + 20, 250, "ANM mean P_f", size=14, color=tk["muted"])
            c.box(x + 20, 270, 340, 18, r=8, fc=tk["code_bg"], z=2)
            P["pf_bar"] = c.box(x + 20, 270, max(8, int(340 * self.pf[i])), 18, r=8, fc=tk["anm"], z=3)
            P["pf_val"] = c.text(x + 20, 310, f"{self.pf[i]:.3f}", size=28, weight="bold", color=tk["ink"])
            P["abs_label"] = c.text(x + 20, 360, "abstain O0", size=14, color=tk["muted"])
            P["anm_abs"] = c.text(x + 20, 400, f"ANM  {self.abstain[i]}", size=22, weight="bold", color=tk["anm"])
            P["ted_abs"] = c.text(x + 20, 440, f"TEDDY  {self.teddy_abs[i]}", size=22, weight="bold", color=tk["gray"])
            P["q"] = c.text(x + 20, 490, f"Q≈{self.q[i]:.3f}", size=16, color=tk["ink2"])
            P["dim"] = c.box(x, 150, 380, 500, r=14, fc=tk["bg"], ec="none", alpha=0.0, z=20)
            self.panels.append(P)

        self.badge = c.chip(
            40, 670,
            "TEDDY-alone over-answers under adt_only · ANM P_f tracks missing RNA · no retrain",
            fc=mix(tk["anm"], tk["panel"], 0.85), color=tk["ink"], size=13,
        )
        c.footer(
            "Numbers from missing_modality_results.json  ·  do not claim win over Pearson 0.61  ·  not clinical"
        )
        self.phase_chip = c.chip(1120, 26, MASKS[0], fc=mix(tk["anm"], tk["panel"], 0.75), color=tk["ink"], size=13)

    def _phase(self, t):
        if t < 2.2:
            return 0, ease(seg(t, 0.3, 1.0))
        if t < 4.6:
            return 1, ease(seg(t, 2.2, 3.0))
        return 2, ease(seg(t, 4.6, 5.4))

    def update(self, t):
        tk = self.tk
        phase, _ = self._phase(t)
        # update phase chip text
        self.phase_chip[1].set_text(MASKS[phase])
        for i, P in enumerate(self.panels):
            # dim non-active panels slightly
            P["dim"].set_alpha(0.0 if i == phase else 0.35)


def render_dashboard(theme="light"):
    d = load()
    c = Canvas(theme)
    tk = c.tk
    n = d["n_base_cells"]
    source = d.get("source") or d.get("export_manifest", {}).get("source")
    c.header(
        "MISSING MODALITY DASHBOARD",
        "TEDDY alone vs TEDDY+ANM under rna_only / adt_only / joint",
        f"n={n:,}  ·  source={source}",
    )
    xs = [40, 440, 840]
    for x, m in zip(xs, MASKS):
        a0 = d["by_mask"][m]["anm"]["criteria"]["O0"]
        a2 = d["by_mask"][m]["anm"]["criteria"]["O2"]
        t0 = d["by_mask"][m]["teddy_alone"]["criteria"]["O0"]
        attr = d["by_mask"][m]["attribution"]
        c.box(x, 140, 380, 520, r=14, fc=tk["panel"], ec=tk["edge"])
        c.system_header(x, 140, 380, m, f"panel r≈{(d.get('export_manifest') or {}).get('panel_pearson_by_mask', {}).get(m)}", tk["anm"])
        c.text(x + 20, 250, f"ANM P_f  {a0.get('mean_P_f'):.3f}", size=20, weight="bold", color=tk["anm"])
        c.text(x + 20, 290, f"ANM O0 abstain  {a0['abstain_count']}", size=16, color=tk["ink"])
        c.text(x + 20, 320, f"ANM O2 abstain  {a2['abstain_count']}", size=16, color=tk["ink"])
        c.text(x + 20, 360, f"TEDDY O0 abstain  {t0['abstain_count']}", size=16, color=tk["gray"])
        c.text(x + 20, 410, f"ANM O0 Q  {a0.get('Q_analogue')}", size=15, color=tk["ink2"])
        c.text(x + 20, 440, f"TEDDY O0 Q  {t0.get('Q_analogue')}", size=15, color=tk["ink2"])
        c.text(
            x + 20, 490,
            f"LOO n_attr={attr.get('n_attr_cells')}  flip={attr.get('top1_flip_rate')}",
            size=13, color=tk["muted"],
        )
        silent = d["by_mask"][m]["teddy_alone"].get("silent_over_answer_O0_vs_O2", {})
        c.text(x + 20, 530, f"silent O0→O2: {silent.get('count')}", size=13, color=tk["muted"])
    c.footer("Do not claim win over phase-1 Pearson ~0.61. Not clinical.")
    img = c.frame_rgb()
    out = OUT / f"missing_modality_dashboard_{theme}.png"
    img.save(out)
    plt.close(c.fig)
    return out


def main():
    print("rendering missing_modality GIF (light)…")
    MissingModalityShift("light").render(str(OUT / "missing_modality_masks_light.gif"))
    print("rendering missing_modality GIF (dark)…")
    MissingModalityShift("dark").render(str(OUT / "missing_modality_masks_dark.gif"))
    print("rendering dashboards…")
    render_dashboard("light")
    render_dashboard("dark")
    print("done →", OUT)
    for p in sorted(OUT.glob("missing_modality*")):
        if p.is_file():
            print(f"  {p.name:48s} {p.stat().st_size/1024:8.1f} KB")


if __name__ == "__main__":
    main()
