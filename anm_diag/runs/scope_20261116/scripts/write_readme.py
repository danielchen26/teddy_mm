#!/usr/bin/env python3
"""Write README.md and the figure of this folder from results/summary.json and the six scope.json files."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RUN = Path(__file__).resolve().parents[1]
RES = RUN / "results"
SUM = json.loads((RES / "summary.json").read_text())
R = {p.parent.name: json.loads(p.read_text()) for p in RES.glob("*/scope.json")}
SCOPES = ["frozen", "train", "other", "units", "own"]
NAMES = {"frozen": "the head\nas it is", "train": "its training\ncells", "other": "another\ncell type", "units": "other units\nof the study", "own": "the unit's\nown cells"}


def f(x, nd=2):
    s = f"{x:.{nd}f}"
    if float(s) == 0:
        s = s.lstrip("-")
    return s.replace("-", "−")


def iv(c, nd=2):
    return f"[{f(c[0], nd)}, {f(c[1], nd)}]"


def yes(b):
    return "holds" if b else "**does not hold**"


def uc(x, nd=2):
    return f"{f(x['u'], nd)} {iv(x['ci'], nd)} {x['verdict']}"


# ------------------------------------------------------------------ the figure: |u| by stage and scope
def panel(ax, title, stages, get, note):
    M = np.full((len(stages), len(SCOPES)), np.nan)
    for i, st in enumerate(stages):
        for j, s in enumerate(SCOPES):
            v = get(st, s)
            if v is not None:
                M[i, j] = abs(v)
    ax.imshow(np.clip(M, 0, 1), cmap="YlOrBr", vmin=0, vmax=1, aspect="auto")
    for i in range(len(stages)):
        for j in range(len(SCOPES)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f(M[i, j]), ha="center", va="center", fontsize=7.5, color="#1b1f24" if M[i, j] < 0.6 else "white")
            else:
                ax.text(j, i, "·", ha="center", va="center", fontsize=8, color="#9aa1ab")
    ax.set_xticks(range(len(SCOPES)))
    ax.set_xticklabels([NAMES[s] for s in SCOPES], fontsize=6.3)
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels(stages, fontsize=7)
    ax.tick_params(length=0)
    ax.xaxis.tick_top()
    ax.set_title(title, fontsize=8, fontweight="bold", loc="left", pad=30)
    ax.text(0, -0.08, note, transform=ax.transAxes, fontsize=6.5, color="#5d6672", va="top")
    for s_ in ax.spines.values():
        s_.set_visible(False)


def figure(path):
    fig, axes = plt.subplots(2, 3, figsize=(10.4, 5.6), dpi=200)

    def pooled(name, target="shift"):
        T = R[name]["results"][0]["targets"][target]["stages"]
        return lambda st, s: (T.get("system", {}).get("frozen", {}).get("u") if (s == "frozen" and st == "output") else
                              (T.get(st, {}).get(s, {}).get("u") if s != "frozen" and st != "output" else None))

    def fam(name, target):
        T = SUM["map"][name][target]["stages"]
        return lambda st, s: (T.get("system", {}).get("frozen", {}).get("mean") if (s == "frozen" and st == "output") else
                              (T.get(st, {}).get(s, {}).get("mean") if s != "frozen" and st != "output" else None))

    ko = ["Z", "HEAD_LAT", "HEAD_OUT", "output"]
    dn = ["RNA512", "Z", "HEAD_LAT", "HEAD_OUT", "output"]
    panel(axes[0, 0], "Knockouts, melanoma: PD-L1 response", ko, pooled("ko_melanoma_pooled"), "five genes pooled")
    panel(axes[0, 1], "Knockouts, THP-1: PD-L1 response", ko, pooled("ko_thp1_pooled"), "five genes pooled")
    panel(axes[0, 2], "External donors: shift in 112 proteins", dn, fam("donors_hao", "shift"), "mean over 8 donors")
    panel(axes[1, 0], "External donors: variation among cells", dn, fam("donors_hao", "variation"), "mean over 8 donors; 1 − mean r²")
    panel(axes[1, 1], "Bone marrow, groups not trained on: shift", dn, lambda st, s: _fam_unseen("shift", st, s), "mean over 4 groups")
    panel(axes[1, 2], "Bone marrow, groups not trained on: variation", dn, lambda st, s: _fam_unseen("variation", st, s), "mean over 4 groups")
    fig.suptitle("Unexplained share |u| by stage (rows) and by the cells the readout was learned from (columns)", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=2.2, w_pad=1.6)
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def _fam_unseen(target, st, s):
    T = SUM["map"]["donors_bm_unseen"][target]["stages"]
    if s == "frozen":
        return T.get("system", {}).get("frozen", {}).get("mean") if st == "output" else None
    return None if st == "output" else T.get(st, {}).get(s, {}).get("mean")


# ------------------------------------------------------------------ the text
def rng_(vals, nd=2):
    return f"{f(min(vals), nd)} to {f(max(vals), nd)}"


def main():
    figure(RES / "scope_map.svg")
    G1, G2, G3, GA = SUM["G1"], SUM["G2"], SUM["G3"], SUM.get("G3_added", {})
    CTS = [k for k, v in G3.items() if isinstance(v, dict)]
    kn = G2["knockouts"]
    fr = SUM["ledger"]["forecast_record"]
    L = []
    a = L.append
    a("# Two finished pieces of work re-read with the version 2 tool — 16 November 2026\n")
    a("Exploratory. The design and the three criteria are in [`DESIGN.md`](DESIGN.md), committed before any script here was run "
      "on data. No model was trained and nothing new was audited: the cells, the head and the measurements are those of earlier "
      "runs. The tool is `anm_diag scope` with the hypotheses and tests of [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md), version 2.\n")
    a("![the map](results/scope_map.svg)\n")
    a("Each cell is |*u*|: the share of a response that a prediction leaves unexplained (0: all of it is carried; 1: none). Rows "
      "are stages of the chain; `output` is the head's own prediction. Columns are the cells the readout was learned from.\n")
    a("## The three criteria\n")
    a("| | criterion | result |")
    a("|---|---|---|")
    a(f"| G1 | one code and one file format for a perturbation and for a contrast nobody assigned | {yes(G1['holds'])}: six SPEC files "
      f"({', '.join(sorted(set(G1['kinds'].values())))}) through one command; none of the functions that compute uses the kind of event |")
    a(f"| G2 | what is already known comes back | {yes(G2['holds'])}: knockouts {yes(all(k['holds'] for k in kn.values()))}; external donors "
      f"{yes(G2['donors']['holds'])} |")
    by = ", ".join(f"{ct} reads '{G3[ct]['reading_of_the_grid']}'" for ct in CTS)
    a(f"| G3 | the grid agrees with interventions already made | **in part**. The rule of the design is met in both cell types (test "
      f"{', '.join(G3[ct]['outcome'] for ct in CTS)}). The reading of the grid: {by}; where it is unresolved, the gap between the scopes "
      f"is resolved and says coverage. The forecast intervals did not hold up: {fr['inside']} of {fr['n']} measured values fell inside "
      f"(see below) |")
    a("")
    a("## What the re-reading says about the tool\n")
    cov = [GA[ct]["share_of_forecast_delivered"] for ct in GA]
    for line in (
            "**The scope of version 1 gives back version 1.** A different estimator on the same cells puts the loss at the same "
            "stage with almost the same numbers (first table below).",
            "**The new axis says something version 1 could not, and it agrees with what the repairs showed.** A readout learned from "
            "the head's training cells reads almost none of the PD-L1 response from the state; one learned from the other cell type "
            "reads little; one learned from the other genes of the same cell type reads it. The repairs that worked were the ones "
            "fitted on cells of the same cell type, and the one fitted on the other cell type did not transfer.",
            f"**A forecast from the grid is a ceiling, not a prediction.** The repairs with cells of this kind delivered "
            f"{' and '.join(str(round(100 * c)) + '%' for c in cov)} of what the grid said a readout could carry at the head's own "
            f"output. `scope.forecast` now says so, and a rule has to ask for less than the ceiling.",
            "**The design's forecast was written for the wrong quantity.** It set |u| at the head's 134 outputs (units scope) against "
            "|u| at the state in another scope, which goes negative for an arm that cannot help. A forecast and its test have to be "
            "about the same quantity. The added check below uses the head's own output on both sides.",
            "**The margin of 0.10 does not fit every target.** For a donor's shift over 112 proteins and for the variation among a "
            "donor's cells, no scope comes within 0.10, not even the donor's own cells, so every donor event is read `unreadable`. "
            "The answer is in the gaps between scopes. A margin has to be declared per target, on its own scale.",
            "**The same picture for an assigned and an unassigned event.** Between groups (a knockout, an external donor) the state "
            "carries more than the head delivers and it takes cells of the same kind to read it. Among the cells of a donor the "
            "head's stages lose nothing; what is missing there is coverage alone."):
        a(f"- {line}")
    a("")
    a("## Knockouts: the scope of version 1 gives back version 1\n")
    a("Pooled event, PD-L1, `units` scope (a readout that never saw the gene it is scored on). The sealed values of version 1 "
      "are beside the new ones.\n")
    a("| cell type | stage | version 2, `units` | version 1, sealed |")
    a("|---|---|---|---|")
    for ct, k in kn.items():
        for st in ("Z", "HEAD_LAT", "HEAD_OUT"):
            a(f"| {ct} | {st} | {uc(k['units'][st])} | {uc(k['version_1_sealed'][st])} |")
    a("")
    for ct, k in kn.items():
        a(f"- {ct}: the head adds {f(k['added_Z_to_HEAD_OUT']['value'])} {iv(k['added_Z_to_HEAD_OUT']['ci'])} from `Z` to `HEAD_OUT` "
          f"(interval above 0.10: {yes(k['added_loss_above_margin'])}); `HEAD_OUT` rejected: {yes(k['head_out_rejected'])}; "
          f"`Z` not rejected: {yes(k['z_not_rejected'])}. Treated side resampled by {k['resampled']['treated']}, so the intervals are "
          f"approximate.")
    a("")
    a("## Knockouts: what is missing\n")
    a("The scopes at `Z`, the stage the head reads, and the head's own PD-L1 output:\n")
    a("| cell type | the head as it is | readout of its training cells | of the other cell type | of the other genes | with the gene's own cells | reading |")
    a("|---|---|---|---|---|---|---|")
    for ct in CTS:
        v = G3[ct]
        s = v["scopes_at_Z"]
        a(f"| {ct} | {uc(v['frozen'])} | {uc(s['train'])} | {uc(s['other'])} | {uc(s['units'])} | {uc(s['own'])} | **{v['reading_of_the_grid']}**"
          + (f"; resolved gaps: {', '.join(v['resolved_gaps'])}" if v["resolved_gaps"] else "") + " |")
    a("")
    for name, ct in (("ko_melanoma_each", "melanoma"), ("ko_thp1_each", "THP-1")):
        m = SUM["map"][name]["shift"]
        zz = m["stages"]["Z"]
        a(f"- {ct}, gene by gene ({m['n_events']} events): mean |u| at `Z` is {f(zz['train']['mean'])} for a readout of the training cells, "
          f"{f(zz['other']['mean'])} for the other cell type, {f(zz['units']['mean'])} for the other genes, {f(zz['own']['mean'])} with the "
          f"gene's own cells. The gap from training cells to other genes is resolved in {m['gaps']['train->units']['resolved_in']} of "
          f"{m['n_events']}. Read as: " + ", ".join(f"{k} {v}" for k, v in m["missing"].items()) + ".")
    a("")
    a("## Knockouts: the grid against the repairs already made\n")
    a("A hypothesis was written from each grid and read against the recorded audits of the repaired heads. It is retrospective: "
      "the outcomes were known. The rule was fixed in the design: a fall above 0.20 for the repair with cells of this kind, below "
      "0.20 for the repair with the other cell type.\n")
    a("**As designed**: the fall of |u| at `HEAD_OUT` in the audit of the changed head (version 1 statistic).\n")
    a("| cell type | arm | forecast from the grid | measured | inside the forecast | test |")
    a("|---|---|---|---|---|---|")
    for ct in CTS:
        v = G3[ct]
        for arm in ("coverage", "other study", "function"):
            fc, m = v["forecast"][arm], v["arms"].get(arm)
            a(f"| {ct} | {arm} | {f(fc['point'])} {iv(fc['ci'])} | " + (f"{f(m['gain'])} {iv(m['ci'])}" if m else "no head was retrained on bone marrow alone")
              + f" | {('yes' if m['inside_forecast'] else 'no') if m else '—'} | {v['outcome'] if arm == 'coverage' else ''} |")
    a("")
    a("The forecast for an arm that cannot help comes out negative: the design subtracted quantities of two different scopes.\n")
    if GA:
        a("**Added after the design** ([`scripts/frozen_arms.py`](scripts/frozen_arms.py)): the same hypothesis on the head's own PD-L1 "
          "output, the quantity `scope.forecast` is about, measured from the saved outputs of the repaired heads.\n")
        a("| cell type | the head's own output as trained | arm | ceiling from the grid | measured fall | test |")
        a("|---|---|---|---|---|---|")
        for ct, v in GA.items():
            for arm in ("coverage", "other study", "function"):
                fc, m = v["forecast"][arm], v["arms"].get(arm)
                a(f"| {ct} | {f(v['as_trained']['u'])} | {arm} | {f(fc['point'])} {iv(fc['ci'])} | "
                  + (f"{f(m['gain'])} {iv(m['ci'])}" if m else "no head was retrained on bone marrow alone") + f" | {v['outcome'] if arm == 'coverage' else ''} |")
        a("")
        a("The repair made elsewhere in the head (its decoder, same cells) lowered the head's own output's |u| by "
          + " and ".join(f"{f(v['elsewhere']['fall'])} ({ct})" for ct, v in GA.items() if v.get("elsewhere"))
          + ". Version 1 found that it removed little of the loss at the head's 134 outputs. The two are different questions: whether "
            "the head's answer is right, and whether the information is in its outputs.\n")
    a(f"Forecast record over both tables: {fr['inside']} of {fr['n']} measured values inside their interval; mean absolute error "
      f"{f(fr['mean_abs_error'])}. Every miss is on the side of the forecast being too high.\n")
    a("## Donors: the same tool on a contrast nobody assigned\n")
    d = G2["donors"]["donors_hao"]
    a(f"External PBMC, 8 donors, variation among a donor's cells (the head's ordinary score), at `Z`. Criteria of the design: the "
      f"gap from the head to a readout of its training cells is smaller than the gap from there to the other donors "
      f"({d['a_head_not_the_bottleneck']} of 8); that gap is positive ({d['b_coverage_positive']} of 8); the gap from the other donors "
      f"to the donor's own cells is smaller ({d['c_own_smaller']} of 8).\n")
    a("| donor | head → its training cells | training cells → other donors | other donors → own cells |")
    a("|---|---|---|---|")
    for r in d["rows"]:
        a(f"| {r['unit'].replace('hao2021_', '')} | {f(r['frozen->train'], 3)} | {f(r['train->units'], 3)} | {f(r['units->own'], 3)} |")
    a("")
    mv, ms = SUM["map"]["donors_hao"]["variation"], SUM["map"]["donors_hao"]["shift"]
    a(f"Over the 8 donors (intervals from resampling donors, approximate): head → training cells "
      f"{f(mv['gaps']['frozen->train']['mean'], 3)} {iv(mv['gaps']['frozen->train']['ci'], 3)}; training cells → other donors "
      f"{f(mv['gaps']['train->units']['mean'], 3)} {iv(mv['gaps']['train->units']['ci'], 3)}; other donors → own cells "
      f"{f(mv['gaps']['units->own']['mean'], 3)} {iv(mv['gaps']['units->own']['ci'], 3)}. These are differences of *u* = 1 − mean *r*², "
      f"not of the correlations reported in `runs/generalization_20261112`; the order is the same. None reaches the step margin of "
      f"0.10, and no scope reads the variation to within 0.10 (the donor's own cells leave {f(mv['stages']['Z']['own']['mean'])}).\n")
    bm = G2["donors"]["donors_bm"]
    a(f"Bone marrow, 12 donor × site groups, the same three comparisons: {bm['a_head_not_the_bottleneck']}, {bm['b_coverage_positive']} "
      f"and {bm['c_own_smaller']} of 12. There the head was trained on the study, the other donors add little, and a group's own "
      f"cells add more than the other donors do.\n")
    a("A donor's shift in protein level (112 proteins), external donors, mean |u|:\n")
    a("| stage | the head as it is | readout of its training cells | of the other donors | with the donor's own cells |")
    a("|---|---|---|---|---|")
    for st in ("RNA512", "Z", "HEAD_LAT", "HEAD_OUT"):
        x = ms["stages"][st]
        a(f"| {st} | {f(ms['stages']['system']['frozen']['mean']) if st == 'HEAD_OUT' else ''} | {f(x['train']['mean'])} | {f(x['units']['mean'])} | "
          f"{f(x['own']['mean'])} |")
    a("")
    hs = SUM["head_stages_for_donors"]
    a(f"Read by the other donors, a donor's shift is carried less at the head's outputs than at the state: |u| rises by "
      f"{f(hs['donors_hao:shift']['mean'])} ({rng_([hs['donors_hao:shift']['min'], hs['donors_hao:shift']['max']])}) from `Z` to "
      f"`HEAD_OUT`, resolved in {hs['donors_hao:shift']['readings']['loss']} of 8 and unresolved in "
      f"{hs['donors_hao:shift']['readings']['unresolved']}. The variation among a donor's cells is not: "
      f"{hs['donors_hao:variation']['readings']['no added loss']} of 8 add no loss. This was not a criterion; it is reported.\n")
    a("## The correction to `runs/generalization_20261112`\n")
    ext = [c for c in SUM["correction"] if c["study"] == "donors_hao"]
    a("That run reported how much of a donor's shift in protein level is left unexplained, as a share of squared differences and, "
      "for the external donors, without a slope between the head's scale and the measurement's. On the scale of *u* (the square "
      "root of the old values) and with the slope (version 2):\n")
    a("| donor | the head: old | version 2 | state, other donors: old | version 2 | state, training cells: old | version 2 | state, own cells: version 2 |")
    a("|---|---|---|---|---|---|---|---|")
    for c in ext:
        o, n = c["old_on_the_scale_of_u"], c["new_u"]
        a(f"| {c['unit'].replace('hao2021_', '')} | {f(o['head'])} | {f(n['head']['u'])} {iv(n['head']['ci'])} | {f(o['state_other_donors'])} | "
          f"{f(n['state_other_donors']['u'])} {iv(n['state_other_donors']['ci'])} | {f(o['state_training_cells'])} | "
          f"{f(n['state_training_cells']['u'])} {iv(n['state_training_cells']['ci'])} | {f(n['state_own_cells']['u'])} {iv(n['state_own_cells']['ci'])} |")
    a("")
    a(f"- The slope changes little: the head leaves {rng_([c['new_u']['head']['u'] for c in ext])} of a donor's shift unexplained "
      f"(old, on this scale: {rng_([c['old_on_the_scale_of_u']['head'] for c in ext])}).")
    a(f"- The scale changes the sentence. The old run said the state, read by the other donors, leaves \"8% to 29%\" and is below the "
      f"margin in 1 of 8. On the scale the margin is defined for it leaves {rng_([c['new_u']['state_other_donors']['u'] for c in ext])}, "
      f"below the margin in {sum(c['new_u']['state_other_donors']['verdict'] == 'support' for c in ext)} of 8. With the donor's own cells "
      f"it leaves {rng_([c['new_u']['state_own_cells']['u'] for c in ext])}.")
    a("- What stands: the state carries much more of a donor's shift than the head delivers, and it takes cells of the same study to "
      "read it. What does not stand: that the shift \"is in the state\". About half of it is, by these readouts.")
    a("- The conclusions of that run that rest on correlation within a donor are not affected.\n")
    a("## An open hypothesis for the donors\n")
    dh = SUM["donor_hypothesis"]
    a(f"Written from the external-donor grid, not retrospective, not tested: **coverage**. Ceiling for the fall of *u* (variation) at "
      f"the head's own output on donors it was not trained on: {f(dh['forecast']['coverage']['mean'], 3)} "
      f"{iv(dh['forecast']['coverage']['ci'], 3)} if PBMC cells of other donors are added to its training; "
      f"{f(dh['forecast']['function']['mean'], 3)} {iv(dh['forecast']['function']['ci'], 3)} for another map on the same bone-marrow "
      f"cells. Rule, fixed in the design: confirmed if the first lies above {f(dh['rule']['at_least'], 3)} and the other two arms below "
      f"it. If the repairs of the knockout loop are a guide, about two thirds of the ceiling is what to expect.\n")
    a("## The ledger\n")
    lg = SUM["ledger"]
    a(f"[`results/ledger.jsonl`](results/ledger.jsonl): {lg['events']} typed events ("
      + ", ".join(f"{v} {k}" for k, v in sorted(lg["by_type"].items())) + f"). Opposite verdicts at the same readout, stage, scope and "
      f"target: {lg['opposite_verdicts']}, all between different events or studies; contradictions: {lg['contradictions']}. Ledger rules: "
      + "; ".join(sorted({f"{c['code']} ({c['severity']})" for c in lg["checks"]})) + ": the four tests are marked retrospective, and the donor "
        "hypothesis is open.\n")
    a("## Limits\n")
    for s_ in ("A re-reading of finished work. The hypotheses for the knockouts were written after their outcomes were known; "
               "agreement is a check of consistency, not a confirmation.",
               "In the table as designed, the measured falls use version 1's estimator and the forecasts version 2's.",
               "The `other study` arm's interval in that table takes its baseline at a point value.",
               "Five genes of one pathway (intervals over five clusters are approximate); eight external donors of one study, all "
               "healthy; one head, one training seed.",
               "A donor contrast is confounded with site, sample and run in bone marrow.",
               "Readouts are ridge regressions, and boosting for the pooled knockout events. Donor fits use at most 30,000 cells.",
               "The first launch of the runs was stopped before it wrote anything, to add the stage-to-stage loss to the output; a "
               "second stopped at a programming error in the summary of a family, before writing. Neither had been read."):
        a(f"- {s_}")
    a("")
    a("## Files\n")
    a("- [`DESIGN.md`](DESIGN.md); [`scripts/build.py`](scripts/build.py) (arrays and SPEC files), [`scripts/run.py`](scripts/run.py) "
      "(the six runs), [`scripts/frozen_arms.py`](scripts/frozen_arms.py) (added), [`scripts/read.py`](scripts/read.py) (criteria, "
      "hypotheses, ledger), [`scripts/write_readme.py`](scripts/write_readme.py)")
    a("- `results/<name>/`: `SPEC.yaml`, `scope.json`, `scope_events.jsonl` for each of the six runs (the resamples are not committed)")
    a("- [`results/summary.json`](results/summary.json), [`results/ledger.jsonl`](results/ledger.jsonl), "
      "[`results/frozen_arms.json`](results/frozen_arms.json)")
    (RUN / "README.md").write_text("\n".join(L) + "\n")
    print("wrote README.md and the figure")


if __name__ == "__main__":
    main()
