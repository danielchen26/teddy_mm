#!/usr/bin/env python3
"""Write README.md and the figure of this folder from results/development.json, results/forecasts.json and, when it is
there, results/what_is_left.json. Nothing is fitted or scored here.

  python write_readme.py

Contrasts that scripts/read.py and scripts/grid_changed.py do not write are computed here with the function they use
(`scope.arms_family`), from the per-donor values of development.json and of results/grid_changed/*/scope.json: four
donors against one, 20,000 cells against 5,000, one donor of the same batch against one of the other; and, for D6,
the head's own output against a readout of the other scored donors. They are marked as such in the README and
written to results/readme_contrasts.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import common
import grid_changed
import plan
import write_final
from anm_diag import scope

RUN = Path(__file__).resolve().parents[1]
RES = RUN / "results"
D = json.loads((RES / "development.json").read_text())
FC = json.loads((RES / "forecasts.json").read_text())["forecast"]
LEFT = json.loads((RES / "what_is_left.json").read_text()) if (RES / "what_is_left.json").exists() else None
if LEFT and not (LEFT["seeds"] == list(plan.SEEDS) and LEFT["fit_sets"] == ["M1", "M2"]):
    LEFT = None
BOUND = D["numbers"]["bound_of_the_rule"]
GUARD = D["numbers"]["guard"]
SEEDS = D["seeds"]
G = D["gains_over_as_is"]
INK, MUTED, GRID = "#1b1f24", "#5d6672", "#e6e9ee"
HUE = {"claim": "#2f6fb5", "dose": "#c9692a", "control": "#7b5ea7"}        # checked with the palette validator: separable under colour-vision deficiency


def f(x, nd=3, sign=False):
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    if float(s) == 0:
        s = s.lstrip("+-")
    return s.replace("-", "−")


def iv(c, nd=3):
    return f"[{f(c[0], nd)}, {f(c[1], nd)}]"


def pct(x):
    return f"{round(100 * x)}%"


def U(arm):
    """u of an arm, donors by seeds, from the per-donor values the reading script wrote."""
    return np.array([[D["per_donor"][d][f"seed{s}"][arm]["u"] for s in SEEDS] for d in plan.DONORS])


def contrast(base, arm):
    """u of `base` minus u of `arm` over donors, as read.py computes a gain."""
    r = scope.arms_family({base: U(base), arm: U(arm)}, base=base)["arms"][arm]
    return {"gain": r["gain"], "ci": r["ci"]}


EXTRA = {"four_minus_one": contrast("one_donor_own_batch", "four_donors"), "more_cells": contrast("four_donors", "coverage_seen"),
         "own_minus_other_batch": contrast("one_donor_other_batch", "one_donor_own_batch")}
_other = np.array(G["one_donor_other_batch"]["per_unit"])
_b1 = np.array([plan.BATCH_OF[d].endswith("1") for d in plan.DONORS])
EXTRA["one_donor_other_batch_by_scored_batch"] = {"scored donors of batch 1 (fit donor of batch 2)": [float(_other[_b1].min()), float(_other[_b1].max())],
                                                  "scored donors of batch 2 (fit donor of batch 1)": [float(_other[~_b1].min()), float(_other[~_b1].max())]}


def left_gaps():
    """D6: u of the head's own output minus u of a readout of the other scored donors at the state, per head."""
    todo = grid_changed.runs(list(plan.SEEDS), ["M1", "M2"])
    donors = sorted({d for r in todo for d in r["scored"]})
    A = {h: {c: np.full((len(donors), len(plan.SEEDS)), np.nan) for c in ("frozen", "units")} for h in grid_changed.HEADS}
    for r in todo:
        for e in json.loads((RES / "grid_changed" / r["name"] / "scope.json").read_text())["results"]:
            T = e["targets"]["variation"]
            i, j = donors.index(e["event"]["unit"]), plan.SEEDS.index(r["seed"])
            A[r["head"]]["frozen"][i, j] = grid_changed.cell(T, "frozen", "system")
            A[r["head"]]["units"][i, j] = grid_changed.cell(T, "units", "Z")
    out = {}
    for h in grid_changed.HEADS:
        g = scope.arms_family({"head": A[h]["frozen"], "readout": A[h]["units"]}, base="head")["arms"]["readout"]
        out[h] = {"gain": g["gain"], "ci": g["ci"]}
    return out


if LEFT:
    EXTRA["D6_head_output_minus_readout_of_other_scored_donors"] = left_gaps()
HEADS = D["heads"]
by_arm = lambda arm: {n: t for n, t in HEADS.items() if t["arm"] == arm}
ROWS = [("h", "20,000 PBMC cells of four donors added"),
        ("coverage_seen", "the scored donor's batch in the fit", "claim"),
        ("coverage_not_seen", "the scored donor's batch not in the fit", "claim"),
        ("h", "5,000 PBMC cells added"),
        ("four_donors", "four donors", "dose"),
        ("one_donor_own_batch", "one donor of the scored donor's batch", "dose"),
        ("one_donor_other_batch", "one donor of the other batch", "dose"),
        ("h", "controls"),
        ("amount", "20,000 fewer bone-marrow cells, nothing added", "control"),
        ("function", "a wider head, twice the epochs", "control"),
        ("labels", "PBMC proteins shuffled within the donor", "control")]
COST = [("coverage", "20,000 PBMC cells", "claim"), ("four_donors", "5,000 cells, four donors", "dose"), ("one_donor", "5,000 cells, one donor", "dose"),
        ("amount", "fewer bone-marrow cells", "control"), ("function", "wider head, longer", "control"), ("labels", "proteins shuffled", "control")]


# ------------------------------------------------------------------ the figure
def clean(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#9aa1ab")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7, color="#9aa1ab", labelcolor=INK)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def head_text(ax, y, s):
    """A group heading in the column of the row labels."""
    ax.text(-0.018, y, s, fontsize=7.4, fontweight="bold", va="center", ha="right", color=INK, transform=ax.get_yaxis_transform())


def panel_gains(ax):
    lo, hi = -0.02, 0.09
    y, ticks, labels = 0.0, [], []
    ax.axvline(0, color="#9aa1ab", lw=0.9, zorder=1)
    ax.axvline(BOUND, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    for row in ROWS:
        if row[0] == "h":
            y += 0.25
            head_text(ax, y, row[1])
            y += 0.95
            continue
        k, name, cls = row
        g, c = G[k], HUE[cls]
        if k.startswith("coverage"):
            fc = FC["pbmc_seen" if k == "coverage_seen" else "pbmc_not_seen"]
            ax.fill_betweenx([y - 0.3, y + 0.3], fc["ci"][0], fc["ci"][1], color="#9aa1ab", alpha=0.22, lw=0, zorder=1)
            ax.plot([fc["gain"]] * 2, [y - 0.3, y + 0.3], color=MUTED, lw=1.1, zorder=2)
            if k == "coverage_seen":
                ax.text(fc["gain"], y - 0.42, "ceiling: a readout fitted\non the same donors", fontsize=6.2, color=MUTED, ha="center", va="bottom", linespacing=1.15)
        if g["ci"][1] < lo:
            ax.annotate("", xy=(lo + 0.001, y), xytext=(lo + 0.014, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.6, mutation_scale=9))
            ax.text(lo + 0.0165, y, "off the scale", fontsize=6.4, color=MUTED, va="center", ha="left", bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=6)
        else:
            ax.plot(g["ci"], [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=3)
            ax.scatter([g["gain"]], [y], s=34, color=c, edgecolor="white", linewidth=0.8, zorder=5)
            ax.scatter([p["gain"] for p in g["per_fit"]], [y + 0.3] * len(g["per_fit"]), s=9, facecolor="white", edgecolor=c, linewidth=0.8, zorder=4)
        ax.text(1.02, y, f"{f(g['gain'], 3, True)} {iv(g['ci'])}", fontsize=6.8, color=INK, va="center", ha="left", transform=ax.get_yaxis_transform())
        ticks.append(y)
        labels.append(name)
        y += 1.0
    top = -1.25
    ax.set_ylim(y - 0.35, top)
    ax.set_xlim(lo, hi)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks([-0.02, 0, 0.02, 0.04, 0.06, 0.08])
    ax.text(BOUND, top + 0.12, f"the rule's bound, {f(BOUND)}", fontsize=6.2, color=INK, ha="center", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    ax.text(1.02, top + 0.12, "mean [interval]", fontsize=6.2, color=MUTED, ha="left", va="top", transform=ax.get_yaxis_transform())
    ax.set_xlabel("fall of u from the head fitted as registered, same seed (higher is better)", fontsize=7.2, color=INK)
    clean(ax)


def panel_donors(ax):
    y, ticks, labels = 0.0, [], []
    a, s, n = U("as_is").mean(1), U("coverage_seen").mean(1), U("coverage_not_seen").mean(1)
    for b in ("hao2021_Batch1", "hao2021_Batch2"):
        for i, d in enumerate(plan.DONORS):
            if plan.BATCH_OF[d] != b:
                continue
            ax.plot([min(s[i], n[i]), a[i]], [y, y], color="#c9ced6", lw=1.0, zorder=1)
            ax.scatter([a[i]], [y], s=26, color=INK, zorder=3)
            ax.scatter([s[i]], [y - 0.15], s=24, color=HUE["claim"], edgecolor="white", linewidth=0.5, zorder=4)
            ax.scatter([n[i]], [y + 0.15], s=20, facecolor="white", edgecolor=HUE["claim"], linewidth=1.1, zorder=4)
            ticks.append(y)
            labels.append(f"{d.replace('hao2021_', '')} · batch {b[-1]}")
            y += 1.0
        y += 0.45
    ax.set_ylim(y - 0.9, -0.6)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(0.50, 0.70)
    ax.set_xticks([0.50, 0.55, 0.60, 0.65, 0.70])
    ax.set_xlabel("u of the donor (lower is better)", fontsize=7.2, color=INK)
    hs = [plt.Line2D([], [], marker="o", ls="", ms=5, mfc=INK, mec=INK, label="fitted as registered"),
          plt.Line2D([], [], marker="o", ls="", ms=5, mfc=HUE["claim"], mec=HUE["claim"], label="PBMC cells added, the donor's batch in the fit"),
          plt.Line2D([], [], marker="o", ls="", ms=4.6, mfc="white", mec=HUE["claim"], label="PBMC cells added, the donor's batch not in the fit")]
    ax.legend(handles=hs, fontsize=6.4, frameon=False, loc="upper left", bbox_to_anchor=(-0.42, -0.2), ncol=1, handletextpad=0.2, labelspacing=0.3)
    clean(ax)


def panel_cost(ax):
    ax.axvline(0, color="#9aa1ab", lw=0.9, zorder=1)
    ax.axvline(-GUARD, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    rng = np.random.default_rng(0)
    for y, (arm, name, cls) in enumerate(COST):
        v = np.array([t["gate6"]["change"] for t in by_arm(arm).values()])
        ax.scatter(v, y + rng.uniform(-0.16, 0.16, len(v)), s=13, color=HUE[cls], alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
        ax.plot([v.mean()] * 2, [y - 0.32, y + 0.32], color=INK, lw=1.4, zorder=4)
    ax.set_yticks(range(len(COST)))
    ax.set_yticklabels([name for _, name, _ in COST], fontsize=7)
    ax.set_ylim(len(COST) - 0.45, -1.15)
    ax.set_xlim(-0.024, 0.022)
    ax.set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
    ax.text(-GUARD + 0.0008, -0.95, f"the guard, −{GUARD:.2f}", fontsize=6.2, color=INK, ha="left", va="center")
    ax.set_xlabel("change of the bone-marrow score, site 4\n(a dot: one head; the bar: their mean)", fontsize=7.2, color=INK, linespacing=1.25)
    clean(ax)


def figure(path):
    fig = plt.figure(figsize=(11.2, 6.0), dpi=200)
    a = fig.add_axes([0.232, 0.135, 0.315, 0.775])
    b = fig.add_axes([0.790, 0.585, 0.195, 0.325])
    c = fig.add_axes([0.790, 0.135, 0.195, 0.215])
    panel_gains(a)
    panel_donors(b)
    panel_cost(c)
    for ax, x0, letter, title in ((a, 0.012, "a", "What each change of the training cells does on donors the head was not fitted on"),
                                  (b, 0.672, "b", "The same fall in each of the eight donors"), (c, 0.672, "c", "What it costs on the head's own task")):
        y1 = ax.get_position().y1
        fig.text(x0, y1 + 0.034, letter, fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")
        fig.text(x0 + 0.017, y1 + 0.036, title, fontsize=8.2, fontweight="bold", color=INK, va="bottom", ha="left")
    fig.text(0.012, 0.018, "Panel a. Filled mark and line: the mean over 8 donors and 3 seeds with its interval (1,000 resamples of donors, seeds resampled inside each arm). "
                           "Small open marks: each seed.", fontsize=6.4, color=MUTED, ha="left", va="bottom")
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"))
    plt.close(fig)


# ------------------------------------------------------------------ the README
def reads(k):
    g = G[k]
    return "above" if g["ci"][0] >= BOUND else ("below" if g["ci"][1] < BOUND else "crosses it")


def gain_row(k, name, cells):
    g = G[k]
    return (f"| {name} | {cells} | {f(D['u'][k]['mean'])} | {f(g['gain'], 3, True)} {iv(g['ci'])} | "
            + ", ".join(f(p["gain"], 3, True) for p in g["per_fit"]) + f" | {f(D['score_within_donor'][k]['mean'])} | {reads(k)} |")


def main():
    figure(RES / "development.svg")
    L = []
    a = L.append
    cs, cn = G["coverage_seen"], G["coverage_not_seen"]
    d1, d2, d4, d5 = D["D1_batch_seen"], D["D2_batch_not_seen"], D["D4_levels_or_relation"], D["D5_cost"]
    go = D["go_to_gate_5"]
    n_adm = sum(t["gate3"]["admitted"] for t in HEADS.values())
    fin = write_final.load() if write_final.have() else None
    if fin:
        write_final.figure(RES / "final.svg", fin[0], fin[1], fin[3])
        a("# Round 1: a head for donors it has not seen — development (16 November 2026) and the test (17 November 2026)\n")
        a("**One test, read once; everything else is exploratory.** The design is [`DESIGN.md`](DESIGN.md), frozen before any script of the round existed. "
          "The test (gate 5) comes first below; the development it rests on (gates 3, 4 and 6) follows unchanged.\n")
        a(write_final.headline(fin[0], fin[3]) + "\n")
        L.extend(write_final.lines(*fin))
        a("## Development (gates 3, 4 and 6), 16 November 2026\n")
        a("The forecasts ([`results/forecasts.json`](results/forecasts.json)) and the two scripts that read the heads ([`scripts/read.py`](scripts/read.py), "
          "[`scripts/grid_changed.py`](scripts/grid_changed.py)) were committed before any fitted head was scored on a PBMC cell. "
          "The eight PBMC donors below are the donors the forecast was read on, so **this part is development, not the test**.\n")
    else:
        a("# Round 1: a head for donors it has not seen — development (gates 3, 4 and 6) — 16 November 2026\n")
    if not fin:
        a("**Exploratory.** The design is [`DESIGN.md`](DESIGN.md), frozen before any script of the round existed. The forecasts "
          "([`results/forecasts.json`](results/forecasts.json)) and the two scripts that read the heads ([`scripts/read.py`](scripts/read.py), "
          "[`scripts/grid_changed.py`](scripts/grid_changed.py)) were committed before any fitted head was scored on a PBMC cell. "
          "The eight PBMC donors below are the donors the forecast was read on, so **this is development, not the test**. "
          "The final study F has not been opened and gate 5 has not been run.\n")
    a("![development](results/development.svg)\n")
    a(f"*u* is the share of the variation among a donor's cells that the head's output leaves unexplained: 1 minus the mean over "
      f"{D['proteins']} proteins of r² within the donor. Every donor is scored only by heads that were not fitted on it. "
      "A gain is the fall of *u* from the head fitted as registered (`as_is`) with the same seed, on the same cells.\n")
    a("## Outcome of development\n" if fin else "## Outcome\n")
    a(f"- **{n_adm} of {len(HEADS)} heads are admitted** (gate 3).")
    a(f"- **Adding PBMC cells to the head's training lowers *u* on donors it was not fitted on by {f(cs['gain'])} {iv(cs['ci'])}**, "
      f"from {f(D['u']['as_is']['mean'])} to {f(D['u']['coverage_seen']['mean'])}; in every seed "
      f"({', '.join(f(p['gain']) for p in cs['per_fit'])}). In the head's ordinary score, the mean correlation within a donor, "
      f"that is {f(D['score_within_donor']['as_is']['mean'])} to {f(D['score_within_donor']['coverage_seen']['mean'])}.")
    a(f"- **It does so as much when the donor's batch was not in the fit**: {f(cn['gain'])} {iv(cn['ci'])}; the difference between the two is "
      f"{f(d2['seen_minus_not_seen']['gain'], 3, True)} {iv(d2['seen_minus_not_seen']['ci'])}.")
    a(f"- **No control does**: 20,000 fewer bone-marrow cells {f(G['amount']['gain'], 3, True)} {iv(G['amount']['ci'])}; a wider head trained "
      f"twice as long {f(G['function']['gain'], 3, True)} {iv(G['function']['ci'])}; the PBMC proteins shuffled within the donor "
      f"{f(G['labels']['gain'], 3, True)} {iv(G['labels']['ci'])}.")
    a(f"- **The head's own task is not paid for it**: on bone marrow, site 4, the score of the heads with PBMC cells changes by "
      f"{f(d5['coverage']['site4_change'], 4, True)} on average; every head of the round stays within the guard of {GUARD:.2f}.")
    cond = [("`coverage` is admitted in every seed", go["coverage_admitted_in_every_seed"]), ("its gain with the batch not seen has its interval above 0",
            go["gain_with_batch_not_seen_above_0"]), ("the guard holds", go["guard_holds"])]
    a(f"- **The go rule of gate 4 {'holds' if go['go'] else 'does not hold'}**: " + "; ".join(f"{t} ({'yes' if ok else 'no'})" for t, ok in cond) + ". "
      + ("Gate 5 was run on that basis on 17 November (above).\n" if fin else
         "Gate 5, the one test on a study nobody has looked at, may be run once F is prepared and the owner of the repository approves opening it.\n"))
    a("With eight donors every interval here is approximate.\n")
    a("## The arms\n")
    a(f"| arm | PBMC cells in the fit | *u* | fall of *u* from `as_is` | by seed | ordinary score | against the rule's bound {f(BOUND)} |")
    a("|---|---|---|---|---|---|---|")
    a(f"| `as_is`: the registered recipe | none | {f(D['u']['as_is']['mean'])} | | | {f(D['score_within_donor']['as_is']['mean'])} | |")
    a(gain_row("coverage_seen", "`coverage`, the scored donor's batch in the fit (M1, M2)", "20,000 from 4 donors"))
    a(gain_row("coverage_not_seen", "`coverage`, the batch not in the fit (B1, B2)", "20,000 from 4 donors"))
    a(gain_row("four_donors", "`four_donors`", "5,000 from 4 donors"))
    a(gain_row("one_donor_own_batch", "`one_donor`, of the scored donor's batch", "5,000 from 1 donor"))
    a(gain_row("one_donor_other_batch", "`one_donor`, of the other batch", "5,000 from 1 donor"))
    a(gain_row("amount", "`amount`: 20,000 fewer bone-marrow cells", "none"))
    a(gain_row("function", "`function`: hidden width 1,024, 80 epochs", "none"))
    a(gain_row("labels", "`labels`: as `coverage`, proteins shuffled within the donor", "20,000 from 4 donors"))
    a("")
    a(f"The head in use (the phase-1 head, one seed) has *u* = {f(D['u']['head_in_use']['mean'])} on the same cells; `as_is` refits its recipe with three seeds "
      f"({', '.join(f(x) for x in D['u']['as_is']['by_seed'])}).\n")
    a("The rule of `donors_hao_coverage` asks the interval of `coverage` to lie above the bound and those of its controls `function` and `amount` below it. "
      "Read on these donors it would be met. These are the donors the forecast came from, so the rule is not decided here. "
      "`labels` is the control of a second hypothesis that gate 5 enters before F is opened.\n")
    a("## What development was read for\n")
    a(f"**D1. Does the arm deliver on donors not seen, and how much of the ceiling?** Yes: {f(cs['gain'])} {iv(cs['ci'])}, "
      f"{pct(d1['of_the_ceiling']['share'])} of the ceiling: what a ridge readout of the state, fitted on the same 20,000 PBMC cells, gains over the head in use "
      f"({f(d1['of_the_ceiling']['ceiling'])} {iv(FC['pbmc_seen']['ci'])}). In *u* itself the head with PBMC cells is at {f(D['u']['coverage_seen']['mean'])} and that readout at "
      f"{f(FC['pbmc_seen']['u'])}. One linear map fitted on bone-marrow and PBMC cells together, without a study "
      f"label, was forecast to gain {f(d1['one_linear_map_forecast'])}; the head gains more.\n")
    a(f"**D2. And when the batch was not seen?** The same: {f(cn['gain'])} {iv(cn['ci'])}, {pct(d2['of_the_ceiling']['share'])} of that ceiling "
      f"({f(d2['of_the_ceiling']['ceiling'])} {iv(FC['pbmc_not_seen']['ci'])}; *u* {f(D['u']['coverage_not_seen']['mean'])} against the readout's {f(FC['pbmc_not_seen']['u'])}; "
      f"one linear map: {f(d2['one_linear_map_forecast'])}). Paired by donor, "
      f"seen minus not seen is {f(d2['seen_minus_not_seen']['gain'], 3, True)} {iv(d2['seen_minus_not_seen']['ci'])}. "
      "With four donors and 20,000 cells in the fit, whether the scored donor's batch is among them makes no difference that these data resolve.\n")
    e1, e2, e3 = EXTRA["four_minus_one"], EXTRA["more_cells"], EXTRA["own_minus_other_batch"]
    eb = list(EXTRA["one_donor_other_batch_by_scored_batch"].values())
    a(f"**D3. Donors or cells?** Cells, as far as this resolves, with one exception. At 5,000 cells, four donors give {f(G['four_donors']['gain'])} "
      f"{iv(G['four_donors']['ci'])} and one donor of the scored donor's batch {f(G['one_donor_own_batch']['gain'])} {iv(G['one_donor_own_batch']['ci'])}: "
      f"{f(e1['gain'], 3, True)} {iv(e1['ci'])} between them, so more donors at the same number of cells add nothing resolved. "
      f"Four times the cells from the same four donors add {f(e2['gain'], 3, True)} {iv(e2['ci'])}. "
      f"The exception: one donor of the other batch gives {f(G['one_donor_other_batch']['gain'])} {iv(G['one_donor_other_batch']['ci'])}, "
      f"{f(e3['gain'])} {iv(e3['ci'])} less than one donor of the same batch; the loss is uneven (with the fit donor from batch 2, the donors of batch 1 gain "
      f"{f(eb[0][0])} to {f(eb[0][1])}; with it from batch 1, those of batch 2 gain {f(eb[1][0])} to {f(eb[1][1])}). "
      "So with a single donor the batch matters, and with four donors and 20,000 cells (D2) "
      "it does not. No arm has 20,000 cells of one donor, or 5,000 cells of four donors of the other batch, so whether donors or cells remove the batch's part is not "
      "separated. The three contrasts of this paragraph are not in `scripts/read.py`; they were computed for this README from the per-donor values with the same function.\n")
    a(f"**D4. Levels or relation?** The relation. With the PBMC proteins shuffled among the cells of each donor, the head's ordinary score on donors it was not fitted on is "
      f"{f(D['score_within_donor']['labels']['mean'])} and *u* is {f(D['u']['labels']['mean'])}: {f(d4['coverage_minus_labels']['gain'])} "
      f"{iv(d4['coverage_minus_labels']['ci'])} worse than `coverage`, and {f(-G['labels']['gain'])} worse than `as_is`. A head that is shown PBMC cells whose proteins do not "
      "follow their RNA stops following the RNA on PBMC cells of other donors too, and keeps its bone-marrow score (below). Read: the head tells the two tissues apart "
      "in the state and maps each in its own way. The gain of `coverage` comes from the relation between RNA and protein in the added cells, not from their protein levels.\n")
    a(f"**D5. What it costs.** Nothing that the guard sees. On bone marrow, site 4 (never fitted), the mean change of the score is {f(d5['coverage']['site4_change'], 4, True)} "
      f"for `coverage`; on the 22 proteins the PBMC study does not measure, {f(d5['coverage']['site4_proteins_the_second_study_lacks'], 4, True)}. "
      f"`amount` loses {f(-d5['amount']['site4_change'], 3)}. It has the same bone-marrow cells as `coverage` and 30% fewer cells in all, so 30% fewer updates in its 40 epochs; "
      f"`labels`, with the full number of cells and nothing to learn from the added ones, loses {f(-d5['labels']['site4_change'], 3)}. "
      "So most of the loss of `amount` is not from the bone-marrow cells taken away. The number of updates is the likely reason; no arm tests it. "
      f"The wider, longer `function` head gains {f(d5['function']['site4_change'], 3, True)} on site 4 "
      f"(validation donor {f(d5['function']['validation_change'], 3, True)}) and nothing on PBMC donors: more of the same map helps inside the study the head was trained on and not outside it. "
      f"The largest fall of a single protein in any `coverage` head is {d5['coverage']['site4_worst_protein']['name']}, "
      f"{f(d5['coverage']['site4_worst_protein']['change'], 3, True)}, in one head of twelve.\n")
    if LEFT:
        v = LEFT["variation"]
        st, fall, ins = v["steps_at_the_stage_the_head_reads"], v["fall_from_as_is_to_coverage"], v["added_inside_the_head_units_scope"]
        a("**D6. What is left.** The grid of `donors_hao` on the changed head and on `as_is` of the same seed, on the four donors the changed head was not fitted on "
          "(fit sets M1 and M2, three seeds; a run has four donors, so numbers compare inside a pair, not with `runs/scope_20261116`).\n")
        a("| *u* of `variation`, mean over 8 donors and 3 seeds | `as_is` | `coverage` | fall |")
        a("|---|---|---|---|")
        for key, name in (("frozen", "the head's own output"), ("train@Z", "a readout of the head's training cells, at the state"),
                          ("units@Z", "a readout of the other scored donors, at the state"), ("own@Z", "a readout of the donor's own cells, at the state"),
                          ("units@HEAD_LAT", "a readout of the other scored donors, inside the head"),
                          ("units@HEAD_OUT", "a readout of the other scored donors, at the head's output")):
            a(f"| {name} | {f(v['u']['as_is'][key])} | {f(v['u']['coverage'][key])} | {f(fall[key]['gain'], 3, True)} {iv(fall[key]['ci'])} |")
        a("")
        a("The training cells of `coverage` include its 20,000 PBMC cells, so the second row changes with the head. Rows three and four are readouts of the state, "
          "which is TEDDY's and the same for both heads.\n")
        s0, s1, gp = st["as_is"], st["coverage"], EXTRA["D6_head_output_minus_readout_of_other_scored_donors"]
        a(f"- **The gap between the head and what the other donors can teach is mostly closed.** A readout of the other scored donors at the state reads "
          f"{f(gp['as_is']['gain'])} {iv(gp['as_is']['ci'])} more than the head's own output before the change and {f(gp['coverage']['gain'])} {iv(gp['coverage']['ci'])} after it "
          "(computed for this README from the runs' files). That readout is the same for both heads: the state is TEDDY's and was not changed.")
        a(f"- **The changed head reads more than a linear readout of its own training cells.** Its own output leaves {f(-s1['frozen->train']['gain'])} "
          f"{iv([-s1['frozen->train']['ci'][1], -s1['frozen->train']['ci'][0]])} less unexplained than a ridge readout of those cells at the state does; for `as_is` the difference is "
          f"{f(-s0['frozen->train']['gain'])} {iv([-s0['frozen->train']['ci'][1], -s0['frozen->train']['ci'][0]])}.")
        a(f"- **What only the donor's own cells carry is unchanged**: {f(s0['units->own']['gain'])} {iv(s0['units->own']['ci'])}.")
        a(f"- Read at the head's output, a readout of the other scored donors gains {f(fall['units@HEAD_OUT']['gain'], 3, True)} {iv(fall['units@HEAD_OUT']['ci'])} from the change, "
          f"and inside the changed head it reads {f(-ins['coverage']['gain'])} more at the output than at the state.")
        a(f"- With the donor's own cells in the readout *u* is still {f(v['u']['coverage']['own@Z'])}. That part is not reached by any ridge readout of the state; "
          "how much of it is measurement noise of the proteins, these data do not say.")
        a(f"- At the tool's margin of 0.10 every donor reads \"unreadable\" in both heads ({sum(v['readings']['as_is'].values())} of {sum(v['readings']['as_is'].values())} "
          "donor-by-seed readings each). That margin does not fit this target; the reading is in the gaps above. "
          "The `shift` target is written to the file for completeness and is not read: with three control donors it has almost no resolution.\n")
        a("For a next round this leaves, on these donors: a small remainder between the head and the other donors, and the part specific to a donor. Round 2 is not designed here; the test comes first.\n")
    else:
        a("**D6. What is left.** Not read yet: `scripts/grid_changed.py` writes `results/what_is_left.json`.\n")
    a("## Against the forecasts\n")
    a("| | forecast, written before any fit | measured |")
    a("|---|---|---|")
    a(f"| `coverage`, hypothesis `donors_hao_coverage` (a ceiling) | 0.068 to 0.082 | {f(cs['gain'])} {iv(cs['ci'])}: **below the forecast**, above the rule's bound |")
    a(f"| `function`, the same hypothesis | −0.004 to 0.001 | {f(G['function']['gain'], 3, True)} {iv(G['function']['ci'])} |")
    a(f"| `amount`, the same hypothesis (a chosen band, not read from the grid) | −0.05 to 0.05 | {f(G['amount']['gain'], 3, True)} {iv(G['amount']['ci'])} |")
    a(f"| ceiling with the batch seen (gate 0c) | {f(FC['pbmc_seen']['gain'])} {iv(FC['pbmc_seen']['ci'])} | {pct(d1['of_the_ceiling']['share'])} of it delivered |")
    a(f"| ceiling with the batch not seen (gate 0c) | {f(FC['pbmc_not_seen']['gain'])} {iv(FC['pbmc_not_seen']['ci'])} | {pct(d2['of_the_ceiling']['share'])} of it delivered |")
    a(f"| one linear map on both studies, batch seen / not seen (gate 0c) | {f(FC['pooled_seen']['gain'])} / {f(FC['pooled_not_seen']['gain'])} | the head delivers "
      f"{f(cs['gain'])} / {f(cn['gain'])} |")
    a("")
    a("The forecast for `coverage` was read from a readout of the other seven donors (up to 30,000 cells); the arm has four donors and 20,000 cells. "
      f"Of the ceiling of its own configuration, about a {'quarter' if d1['of_the_ceiling']['share'] < 0.8 else 'fifth'} is not delivered with the batch seen. "
      f"The ceilings are gains over the head in use (*u* = {f(FC['head']['u'])}); the arms' gains are over `as_is` of the same seed (mean {f(D['u']['as_is']['mean'])}).\n")
    a("## Gate 3: admission\n")
    a("| arm | heads | admitted | 1. validation score, last 5 epochs minus the 5 before (not below −0.005) | 2. own objective | 4. validation donor against `as_is` (not below −0.02) |")
    a("|---|---|---|---|---|---|")
    for arm in ("as_is", "coverage", "four_donors", "one_donor", "amount", "function", "labels"):
        hs = by_arm(arm)
        ran = [t["gate3"]["ran"]["last5"] - t["gate3"]["ran"]["before5"] for t in hs.values()]
        gd = [t["gate3"]["guard_validation"]["change"] for t in hs.values()]
        if arm == "coverage":
            own = [t["gate3"]["own_objective"]["gain"] for t in hs.values()]
            o = f"score on its own PBMC fit cells above `as_is` by {f(min(own))} to {f(max(own))} (at least {D['numbers']['own_objective']})"
        elif arm == "function":
            o = "training loss below `as_is` of the seed: " + ", ".join(f"{f(t['gate3']['own_objective']['train_loss'])} against {f(t['gate3']['own_objective']['as_is'])}"
                                                                         for t in hs.values())
        else:
            o = "—"
        a(f"| `{arm}` | {len(hs)} | {sum(t['gate3']['admitted'] for t in hs.values())} | {f(min(ran), 4, True)} to {f(max(ran), 4, True)} | {o} | "
          + ("—" if arm == "as_is" else f"{f(min(gd), 4, True)} to {f(max(gd), 4, True)}") + " |")
    a("")
    a("Criterion 3 (no scored cell in a fit, by donor and barcode; bone-marrow fit cells are training cells) holds for every head. No cell of F is in any file a fit reads: F has not been downloaded.\n")
    closest = min(HEADS.items(), key=lambda kv: kv[1]["gate3"]["ran"]["last5"] - kv[1]["gate3"]["ran"]["before5"])
    a(f"The closest call of criterion 1 is `{closest[0]}`: {f(closest[1]['gate3']['ran']['last5'] - closest[1]['gate3']['ran']['before5'], 4, True)} against a limit of −0.005.\n")
    a("## Gate 6: the guard on bone marrow, site 4\n")
    a("| arm | heads | mean change of the score | range over heads | the 22 proteins the PBMC study lacks | largest fall of one protein in one head |")
    a("|---|---|---|---|---|---|")
    for arm in ("coverage", "four_donors", "one_donor", "amount", "function", "labels"):
        hs = by_arm(arm)
        ch = [t["gate6"]["change"] for t in hs.values()]
        w = d5[arm]["site4_worst_protein"]
        a(f"| `{arm}` | {len(hs)} | {f(np.mean(ch), 4, True)} | {f(min(ch), 4, True)} to {f(max(ch), 4, True)} | "
          f"{f(d5[arm]['site4_proteins_the_second_study_lacks'], 4, True)} | {w['name']} {f(w['change'], 3, True)} |")
    a("")
    a(f"The score is the mean over 134 proteins of the correlation between measured and predicted on the site-4 cells; the change is against `as_is` of the same seed. "
      f"The guard is {GUARD:.2f}; every head passes.\n")
    a("## Each donor\n")
    a("| donor | batch | *u*, `as_is` | *u*, `coverage`, batch in the fit | *u*, `coverage`, batch not in the fit | fall, batch in the fit | fall, batch not in the fit |")
    a("|---|---|---|---|---|---|---|")
    ua, us, un = U("as_is").mean(1), U("coverage_seen").mean(1), U("coverage_not_seen").mean(1)
    for i, dn in enumerate(plan.DONORS):
        a(f"| {dn.replace('hao2021_', '')} | {plan.BATCH_OF[dn][-1]} | {f(ua[i])} | {f(us[i])} | {f(un[i])} | {f(ua[i] - us[i])} | {f(ua[i] - un[i])} |")
    a("")
    a("Means over the three seeds. The fall is positive in every donor.\n")
    a("## How the fits were run\n")
    hours = sum(t["minutes"] for t in HEADS.values()) / 60
    a(f"- 45 heads, {hours:.1f} hours of fitting on 2 CPU threads (torch {HEADS['s0_as_is']['machine']['torch']}), in seven jobs of at most two hours. "
      "Four times a job's time ended during a head and twice the machine restarted; each time the head went on from the state kept after its last finished epoch. "
      "The trainer's test shows a resumed fit gives the same weights as an uninterrupted one.")
    a("- The design asks that all arms of a seed are fitted on one machine. Every head records the same device, thread count and library version. "
      "Whether the hardware under it stayed the same across the two restarts is not known.")
    a("- Nothing in `scripts/train_arm.py`, `scripts/plan.py` or `scripts/run_all.py` changed after the first fit started.\n")
    a("## What development does not show\n" if fin else "## What this does not show\n")
    for s in ("It is not the test. The eight donors are the donors the forecast was read on, and the rule is decided on F only.",
              "Eight donors of one study, all healthy, two batches. Nothing here is about disease, another laboratory or another assay.",
              "TEDDY's pretraining data probably held the RNA of this study. \"Not seen\" is a statement about the head, the same for every arm.",
              "The gain is in the variation among a donor's cells. A donor's shift in level against other donors is not read here.",
              "The 112 proteins both panels share; one head recipe; targets of the two studies on their own scales.",
              "A quarter of the ceiling is not delivered, and the ceiling is itself a ridge readout: what a head could reach is not bounded by it.",
              "D3 does not separate donors from cells for a batch that is not in the fit."):
        a(f"- {s}")
    a("")
    a("## Next\n")
    if fin:
        a("The test has been read and is not read again. A next round is not designed here. For what the admitted seeds read, the design names two candidates: a second study "
          "in the fit, or scales of the loss that are specific to a study. For the outcome the rule states, it names the reserve donors, read once; that is left to the owner "
          f"of the repository (the {len(fin[2].get('reserve', []))} donors that were not drawn stay unopened). "
          + ("What the grid on the final study leaves open is under \"What is left on the final study\" above.\n" if fin[4] else
             "The grid on the final study (step 6 of gate 5) was running when this was written; its reading is added when it ends.\n"))
    else:
        a("Gate 5, once: the final study is prepared (criteria checked without reading a protein value, donors listed, cells embedded), the recipe is frozen by a commit, "
          "`coverage` and `labels` are refitted on all eight donors, the second hypothesis and the sealed declaration are written, and then the heads are run on F. "
          "Opening F needs the owner's approval.\n")
    a("## Files\n")
    a("- [`DESIGN.md`](DESIGN.md)")
    a("- [`scripts/plan.py`](scripts/plan.py), [`scripts/train_arm.py`](scripts/train_arm.py), [`scripts/run_all.py`](scripts/run_all.py): the heads and their fits")
    a("- [`scripts/forecasts.py`](scripts/forecasts.py), [`results/forecasts.json`](results/forecasts.json): gate 0c")
    a("- [`scripts/read.py`](scripts/read.py) → [`results/development.json`](results/development.json) (every head, donor and seed) and "
      "[`results/ledger.jsonl`](results/ledger.jsonl) (one intervention event per head)")
    a("- [`scripts/grid_changed.py`](scripts/grid_changed.py) → [`results/what_is_left.json`](results/what_is_left.json) and `results/grid_changed/<run>/`")
    a("- [`scripts/write_readme.py`](scripts/write_readme.py): this file, the figure and [`results/readme_contrasts.json`](results/readme_contrasts.json)")
    if fin:
        a("- [`GATE5.md`](GATE5.md): the record of gate 5, written in parts, each committed before the step it governs")
        a("- [`scripts/final_study.py`](scripts/final_study.py) → [`results/final_selection.json`](results/final_selection.json): the criteria, the donors, the cells and the pack of the final study")
        a("- [`scripts/declare_final.py`](scripts/declare_final.py) → [`results/gate5_declaration.json`](results/gate5_declaration.json) and [`results/ledger_final.jsonl`](results/ledger_final.jsonl)")
        a("- [`scripts/read_final.py`](scripts/read_final.py) → [`results/final.json`](results/final.json) (every head, donor and seed)")
        a("- [`scripts/grid_final.py`](scripts/grid_final.py) → [`results/what_is_left_final.json`](results/what_is_left_final.json) and `results/grid_final/<run>/`")
        a("- [`scripts/write_final.py`](scripts/write_final.py): the part of this file about the test and its figure")
        a("- [`results/gate5_heads.json`](results/gate5_heads.json): the sha256 of the heads of development that the test reuses")
    a(f"- The heads' weights are not in the repository ({48 if fin else 42} files of 5.6 MB and three of 15 MB).")
    (RUN / "README.md").write_text("\n".join(L) + "\n")
    print("wrote README.md and the figure; D6", "included" if LEFT else "not yet there")
    (RES / "readme_contrasts.json").write_text(json.dumps({"what": "contrasts computed by scripts/write_readme.py with scope.arms_family; not in read.py or grid_changed.py",
                                                         **EXTRA}, indent=1))
    print(json.dumps(EXTRA)[:900])


if __name__ == "__main__":
    main()
