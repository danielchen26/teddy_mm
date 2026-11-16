#!/usr/bin/env python3
"""Write README.md and the figure of this folder from results/gate0.json and results/public_data_census.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RES = HERE / "results"
G0 = json.loads((RES / "gate0.json").read_text())
CEN = json.loads((RES / "public_data_census.json").read_text())
THR = 0.02
MARGIN = 0.10


def f(x, nd=3):
    s = f"{x:.{nd}f}"
    if float(s) == 0:
        s = s.lstrip("-")
    return s.replace("-", "−")


def iv(ci, nd=3):
    return f"[{f(ci[0], nd)}, {f(ci[1], nd)}]"


def read(t):
    lo, hi = t["ci"]
    return "resolved" if lo > THR else ("absent" if hi < THR else "unresolved")


def label(g):
    if g["dataset"] == "bm":
        return f"bone marrow, donor {g['donor']}, {g['site'].replace('site', 'site ')}"
    return f"PBMC, donor {g['donor'].replace('hao2021_', '')}"


def kind(g):
    if g["dataset"] == "ext":
        return "external"
    if g["seen_by_head"]:
        return "trained"
    if g["split"] == "val":
        return "validation"
    return "site4_seen_donor" if g["donor"] == "15078" else "site4"


GROUPS = G0["groups"]
UNSEEN = [g for g in GROUPS if not g["seen_by_head"]]
TRAINED = [g for g in GROUPS if g["seen_by_head"]]
STUDY = {"bm": "bone marrow (4 groups the head was not trained on)", "ext": "external PBMC (8 donors)"}
NOTE = {"validation": "donor unseen, site seen; used to choose the checkpoint", "site4": "donor and site unseen",
        "site4_seen_donor": "site unseen; the donor is in training at other sites", "external": "another study and tissue"}
TERMS = [("readout_minus_head", "the head leaves transferable signal undelivered", "readout of the head's training cells − head"),
         ("same_study_donors", "more donors of the same kind are needed", "other donors of the study − the head's training cells"),
         ("own_cells", "part is specific to the donor", "own cells − other donors of the study"),
         ("raw_minus_state_own", "part is lost before the head", "RNA components − state, own cells")]


def sc(g, k):
    return g["scores"][k]["mean"] if k in g["scores"] else None


def cell(g, k):
    return f(sc(g, k)) if sc(g, k) is not None else "—"


def term_cell(g, k):
    if k not in g["terms"]:
        return "—"
    t = g["terms"][k]
    return f"{f(t['value'])} {iv(t['ci'])} {read(t)}"


def study_reading(ds, k):
    gs = [g for g in UNSEEN if g["dataset"] == ds and k in g["terms"]]
    r = [read(g["terms"][k]) for g in gs]
    n = len(r)
    if r and all(x == "resolved" for x in r):
        verdict = "**holds**"
    elif r and all(x == "absent" for x in r):
        verdict = "does not hold (absent in every donor)"
    else:
        verdict = "not stated"
    vals = [g["terms"][k]["value"] for g in gs]
    return f"{verdict}: resolved in {r.count('resolved')} of {n}, absent in {r.count('absent')}, range {f(min(vals))} to {f(max(vals))}"


def figure(path: Path):
    order = [("trained", "trained on (the head is in-sample)"), ("validation", "validation donor"), ("site4", "site 4, not trained on"),
             ("site4_seen_donor", "site 4, donor seen at other sites"), ("external", "external PBMC, not trained on")]
    rows = []
    for k, title in order:
        gs = sorted([g for g in GROUPS if kind(g) == k], key=lambda g: (g["site"], g["donor"]))
        if gs:
            rows.append((title, gs))
    n = sum(len(gs) for _, gs in rows) + len(rows)
    fig, ax = plt.subplots(figsize=(7.6, 0.27 * n + 1.5), dpi=200)
    series = [("A", "the frozen head", "o", "#1b1f24", True), ("Z_T", "ridge readout of the head's training cells", "^", "#2f6fb5", True),
              ("Z_L", "ridge readout of the study's other donors", "s", "#1a8f7a", True),
              ("Z_C", "ridge readout of the donor's own cells", "D", "#c9692a", True),
              ("RAW_C", "the same from 512 RNA components", "D", "#8a8f98", False)]
    y, ticks, labels = 0, [], []
    for title, gs in rows:
        ax.text(0.0, y, title, fontsize=7.6, fontweight="bold", va="center", ha="left", color="#1b1f24",
                transform=ax.get_yaxis_transform())
        y += 1
        for g in gs:
            vals = [sc(g, k) for k, *_ in series if sc(g, k) is not None]
            ax.plot([min(vals), max(vals)], [y, y], color="#c9ced6", lw=1.0, zorder=1)
            for k, _, m, c, filled in series:
                v = sc(g, k)
                if v is None:
                    continue
                ax.scatter([v], [y], marker=m, s=26 if m != "D" else 20, facecolor=c if filled else "white", edgecolor=c,
                           linewidth=0.9, zorder=3)
            ticks.append(y)
            labels.append(f"{g['donor'].replace('hao2021_', '')}" + (f" · {g['site'].replace('site', 'site ')}" if g["site"] else ""))
            y += 1
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(y - 0.4, -0.8)
    allv = [sc(g, k) for g in GROUPS for k, *_ in series if sc(g, k) is not None]
    ax.set_xlim(max(0.0, round(min(allv) - 0.06, 1)), min(1.0, round(max(allv) + 0.06, 1)))
    ax.set_xlabel("mean per-protein Pearson correlation within the group (112 proteins)", fontsize=7.6)
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(axis="x", color="#e6e9ee", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Line2D([], [], marker=m, linestyle="", markersize=5, markerfacecolor=c if filled else "white",
                          markeredgecolor=c, label=lab) for _, lab, m, c, filled in series]
    ax.legend(handles=handles, fontsize=6.8, frameon=False, loc="upper center", bbox_to_anchor=(0.42, -0.09), ncol=2,
              handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def main():
    figure(RES / "gate0_scores.svg")
    L = []
    a = L.append
    a("# Generalization to unseen donors: where is the signal the head does not deliver? — 12 November 2026\n")
    a("Exploratory. The design and the readings are in [`DESIGN.md`](DESIGN.md), committed with "
      "[`scripts/gate0.py`](scripts/gate0.py) before the script was run on data. No model was trained. Every donor here has "
      "been looked at before; nothing is final data.\n")
    a("The first run was killed for lack of memory during the fits, before any result was written. The script was changed to "
      "hold less in memory (fits, folds and seeds unchanged; identical output on synthetic data) and run again; the results "
      "below are from that run.\n")
    a("![scores by group](results/gate0_scores.svg)\n")
    a("Each row is one donor (bone marrow: one donor at one site). The score is the repository's: the Pearson correlation "
      "between predicted and measured protein within the group, averaged over the 112 proteins measured in both studies. "
      "Readouts are ridge regressions on TEDDY's pooled state, scored on cells they were not fitted on.\n")
    a("## Outcome by the readings written in advance\n")
    a("A term is resolved in a donor if its interval lies above 0.02, absent if it lies below. A reading is stated for a "
      "study only if it holds in every donor of that study the head was not trained on.\n")
    a("| reading | term | " + " | ".join(STUDY.values()) + " |")
    a("|---|---|---|---|")
    for k, name, how in TERMS:
        a(f"| {name} | {how} | " + " | ".join(study_reading(ds, k) for ds in STUDY) + " |")
    a("")
    # ---- what the numbers show (computed)
    def rg(vals, nd=3):
        return f"{f(min(vals), nd)} to {f(max(vals), nd)}"

    def pct(vals):
        return f"{round(100 * min(vals))}% to {round(100 * max(vals))}%"

    ub, ue = [g for g in UNSEEN if g["dataset"] == "bm"], [g for g in UNSEEN if g["dataset"] == "ext"]
    tv = lambda gs, k: [g["terms"][k]["value"] for g in gs if k in g["terms"]]
    rd = lambda gs, k: [read(g["terms"][k]) for g in gs if k in g["terms"]]
    g15 = next(g for g in GROUPS if g["id"] == "bm:15078:site4")
    a("## What the numbers show\n")
    a(f"- **The head already delivers most of what is readable.** On the 12 groups it was not trained on, its score is "
      f"{pct([sc(g, 'A') / sc(g, 'Z_C') for g in UNSEEN])} of the score of a readout fitted on the donor's own cells "
      f"(bone marrow {pct([sc(g, 'A') / sc(g, 'Z_C') for g in ub])}, external {pct([sc(g, 'A') / sc(g, 'Z_C') for g in ue])}). "
      f"What separates the groups is mostly that ceiling: {rg([sc(g, 'Z_C') for g in GROUPS], 2)} over the 20 groups, by site. "
      "That explains the table that opened the design.")
    r1 = rd(UNSEEN, "readout_minus_head")
    a(f"- **No sign that the head leaves transferable signal undelivered.** A readout fitted on the head's own training cells "
      f"differs from the head by {rg(tv(UNSEEN, 'readout_minus_head'))}: absent in {r1.count('absent')} of 12, unresolved in "
      f"{r1.count('unresolved')}, resolved in {r1.count('resolved')}. On the same cells a kernel ridge does not beat the head by more than 0.03.")
    a(f"- **External donors: the study is missing, not the donor.** With the other seven donors of the study in the fit the "
      f"score rises by {rg(tv(ue, 'same_study_donors'))} (resolved in 8 of 8), and {rg(tv(ue, 'own_cells'))} then remains to "
      "the donor's own cells. Donors of one study teach each other almost everything. Tissue, laboratory and antibody panel "
      "all differ between the two studies, and these data cannot separate them.")
    a(f"- **Bone marrow: the site and the sample.** Leaving a group's site out of the fit costs "
      f"{rg([sc(g, 'Z_T') - sc(g, 'Z_S') for g in TRAINED] + [g['terms']['unseen_site_cost']['value'] for g in ub if g['split'] == 'val'], 2)}"
      f". The donor's own cells add {rg(tv(ub, 'own_cells'))} over the other donors (resolved in 4 of 4). A bone-marrow group is "
      "one donor at one site, one sample. For donor 15078 at site 4 the three can be told apart: adding that donor's cells from "
      f"the other sites adds {f(g15['terms']['unseen_site_cost']['value'])}, adding the other donors of site 4 adds "
      f"{f(sc(g15, 'Z_L') - sc(g15, 'Z_S'))}, and the donor's own site-4 cells add {f(sc(g15, 'Z_C') - sc(g15, 'Z_L'))} more. "
      "That points to the sample, not the person. The external donors share runs, and there the term is absent.")
    a(f"- **Nothing is lost before the head, against 512 RNA components.** Fitted on a donor's own cells, the state scores "
      f"higher than the RNA components in {sum(sc(g, 'RAW_C') < sc(g, 'Z_C') for g in GROUPS)} of 20 groups, by "
      f"{rg([sc(g, 'Z_C') - sc(g, 'RAW_C') for g in GROUPS])}. In transfer from the head's training cells the two are within "
      f"{f(max(abs(x) for x in tv(UNSEEN, 'raw_minus_state_T')))} of each other.")
    uh = lambda gs, k: [g["shift"]["by"][k]["u"] for g in gs]
    a(f"- **A donor's shift in protein level (external donors).** The head leaves {pct(uh(ue, 'A'))} of it unexplained. The "
      f"state, read by a readout of the other donors, leaves {pct(uh(ue, 'Z_L'))} (below the 0.10 margin in "
      f"{sum(g['shift']['by']['Z_L']['u_ci'][1] < MARGIN for g in ue)} of 8); a readout of the head's training cells leaves "
      f"{pct(uh(ue, 'Z_T'))}. The shift is in the state and the head does not carry it, the pattern of the PD-L1 audits, and "
      f"here it takes data of the same study to read it. In bone marrow the state leaves {pct(uh(ub, 'Z_L'))} unexplained: "
      "part of a sample's shift is not in the RNA.")
    a("")
    a("## What follows for a loop aimed at donors and diseases\n")
    for line in (
            "The unit to hold out is the study or the batch. Holding out donors inside a study the head has seen would test little.",
            "The change to test is coverage: cells of another study in the head's training. Its control is the same number of "
            "cells from studies already covered, so that coverage is separated from amount.",
            "The size to expect on a new study of a kind already covered is the same-study term above, 0.07 to 0.10 here, "
            "measured with a ridge readout. Whether the head gains as much is not shown; that needs a training run.",
            "Repairing the head at a located step, as in the knockout loop, is not what these numbers call for.",
            "Disease is untested. Every donor here is healthy."):
        a(f"- {line}")
    a("")
    a("## Donors the head was not trained on\n")
    a("| donor | cells | the head | readout of the head's training cells | readout of the study's other donors | readout of "
      "the donor's own cells | RNA components, own cells |")
    a("|---|---|---|---|---|---|---|")
    for g in UNSEEN:
        a(f"| {label(g)} ({NOTE[kind(g)]}) | {g['n']:,} | {cell(g, 'A')} | {cell(g, 'Z_T')} | {cell(g, 'Z_L')} | {cell(g, 'Z_C')} | "
          f"{cell(g, 'RAW_C')} |")
    a("")
    a("The differences, with intervals from 300 resamples of the donor's cells:\n")
    a("| donor | readout − head | other donors of the study − head's training cells | own cells − other donors | RNA components "
      "− state (own cells) | cost of an unseen site |")
    a("|---|---|---|---|---|---|")
    for g in UNSEEN:
        a(f"| {label(g)} | " + " | ".join(term_cell(g, k) for k in ("readout_minus_head", "same_study_donors", "own_cells",
                                                                     "raw_minus_state_own", "unseen_site_cost")) + " |")
    a("")
    a("## Groups the head was trained on, for reference\n")
    a("The head is in-sample here, so its score is not comparable with the readouts, which never contain the donor they are "
      "scored on. \"Without the site\" is the readout of the head's training cells fitted without the group's donor and "
      "without its site.\n")
    a("| group | cells | the head (in-sample) | readout, donor left out | readout, donor and site left out | readout of the "
      "study's other donors | own cells |")
    a("|---|---|---|---|---|---|---|")
    for g in TRAINED:
        a(f"| {label(g)} | {g['n']:,} | {cell(g, 'A')} | {cell(g, 'Z_T')} | {cell(g, 'Z_S')} | {cell(g, 'Z_L')} | {cell(g, 'Z_C')} |")
    a("")
    a("## The donor's shift\n")
    a("For each donor, the shift of every protein's mean against the other donors of the same study, in units of their "
      "spread, and the share of it that a prediction leaves unexplained (*u*, over proteins). A shift is carried if the upper "
      f"bound of *u* is below {MARGIN:.2f}. A donor's shift mixes biology with staining and batch.\n")
    a("| donor | size of the shift (RMS, in SD) | the head | state, readout of the study's other donors | RNA components, "
      "readout of the study's other donors | state, readout of the head's training cells |")
    a("|---|---|---|---|---|---|")

    def ucell(g, k):
        b = g["shift"]["by"].get(k)
        if not b:
            return "—"
        tag = "carried" if b["u_ci"][1] < MARGIN else ("not carried" if b["u_ci"][0] > MARGIN else "unresolved")
        return f"{f(b['u'], 2)} {iv(b['u_ci'], 2)} {tag}"

    for g in GROUPS:
        a(f"| {label(g)}{'' if not g['seen_by_head'] else ' (trained on)'} | {f(g['shift']['size_rms_sd'], 2)} | {ucell(g, 'A')} | "
          f"{ucell(g, 'Z_L')} | {ucell(g, 'RAW_L')} | {ucell(g, 'Z_T')} |")
    a("")
    a("## Readout classes and penalties\n")
    cls = {}
    for g in GROUPS:
        for k, v in g["scores"].items():
            if "chosen" in v:
                cls.setdefault(k, []).append(v["chosen"])
    a("Chosen class by fit (of 20 groups): " + "; ".join(f"`{k}` ridge {v.count('ridge')}, with random features {v.count('rff')}"
                                                         for k, v in sorted(cls.items())) + ".\n")
    a("## Public data that could carry a loop aimed at donors and diseases\n")
    a(CEN["what"].capitalize() + ". TEDDY's pretraining data: " + CEN["teddy_pretraining"]["source"] + ".\n")
    a("| data | tissue | donors | disease | proteins | overlap with the head's 134 | access | what it could serve | status |")
    a("|---|---|---|---|---|---|---|---|---|")
    for d in CEN["datasets"]:
        a(f"| {d['name']} ({d.get('accession', '')}) | {d.get('tissue', '')} | {d.get('donors', '')} | {d.get('disease', '')} | "
          f"{d.get('proteins', '')} | {d.get('overlap_with_head', '')} | {d.get('access', '')} | {d.get('role', '')} | "
          f"{d.get('status', '')} |")
    a("")
    for r in CEN["reading"]:
        a(f"- {r}")
    a("")
    a("Sources are listed in [`results/public_data_census.json`](results/public_data_census.json).\n")
    a("## Limits\n")
    for s in ("Twelve donors the head was not trained on, all healthy, none of them new to this repository. The intervals "
              "describe each donor's cells; they say nothing about donors in general and nothing about disease.",
              "The readouts are ridge regressions, not the head's architecture. Where a readout scores below the head, that "
              "does not show that nothing more is available.",
              "The RNA is represented by 512 principal components. Signal outside them is not counted as readable from the RNA.",
              "The score ignores level and scale within a group. The donor's shift looks at level, within a study only, and it "
              "mixes biology with staining.",
              "The external study's RNA may have been in TEDDY's pretraining data.",
              "One head, one training seed."):
        a(f"- {s}")
    a("")
    a("## Files\n")
    a("- [`DESIGN.md`](DESIGN.md), [`scripts/gate0.py`](scripts/gate0.py), [`scripts/write_readme.py`](scripts/write_readme.py)")
    a("- [`results/gate0.json`](results/gate0.json): every score, term and interval; "
      "[`results/gate0_per_protein.json`](results/gate0_per_protein.json): the same by protein")
    a("- [`results/public_data_census.json`](results/public_data_census.json)")
    (HERE / "README.md").write_text("\n".join(L) + "\n")
    print("wrote README.md and the figure;", len(UNSEEN), "unseen groups,", len(TRAINED), "trained")


if __name__ == "__main__":
    sys.exit(main())
