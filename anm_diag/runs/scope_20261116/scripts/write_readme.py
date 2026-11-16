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
def main():
    figure(RES / "scope_map.svg")
    G1, G2, G3 = SUM["G1"], SUM["G2"], SUM["G3"]
    L = []
    a = L.append
    a("# Two finished pieces of work re-read with the version 2 tool — 16 November 2026\n")
    a("Exploratory. The design and the three criteria are in [`DESIGN.md`](DESIGN.md), committed before any script here was run "
      "on data. No model was trained and nothing new was audited: the cells, the head and the measurements are those of earlier "
      "runs. The tool is `anm_diag scope` with the hypotheses and tests of [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md), version 2.\n")
    a("![the map](results/scope_map.svg)\n")
    a("Each cell is |*u*|: the share of the response a prediction leaves unexplained (0: all of it is carried; 1: none). Rows are "
      "stages of the chain; `output` is the head's own prediction. Columns are the cells the readout was learned from.\n")
    a("## The three criteria\n")
    a("| | criterion | result |")
    a("|---|---|---|")
    a(f"| G1 | one code and one file format for a perturbation and for a contrast nobody assigned | {yes(G1['holds'])}: six SPEC files "
      f"({', '.join(sorted(set(G1['kinds'].values())))}) through one command; none of the functions that compute reads the kind of event |")
    kn = G2["knockouts"]
    a(f"| G2 | what is already known comes back | {yes(G2['holds'])}: knockouts {yes(all(k['holds'] for k in kn.values()))}; donors "
      f"{yes(G2['donors']['holds'])} (see below) |")
    a(f"| G3 | the grid agrees with interventions already made | {yes(G3['holds'])} by the verdicts"
      + ("" if G3["holds"] else f"; {yes(G3['holds_by_gap'])} if the resolved gap between scopes is taken where no verdict is reached") + " (see below) |")
    a("")
    a("## Knockouts: the scope of version 1 gives back version 1\n")
    a("Pooled event, PD-L1, `units` scope (a readout that never saw the gene it is scored on). Version 1 used a different estimator "
      "on the same cells; the sealed values are beside the new ones.\n")
    a("| cell type | stage | version 2, `units` | version 1, sealed |")
    a("|---|---|---|---|")
    for ct, k in kn.items():
        for st in ("Z", "HEAD_LAT", "HEAD_OUT"):
            a(f"| {ct} | {st} | {uc(k['units'][st])} | {uc(k['version_1_sealed'][st])} |")
    a("")
    for ct, k in kn.items():
        a(f"- {ct}: the head adds {f(k['added_Z_to_HEAD_OUT']['value'])} {iv(k['added_Z_to_HEAD_OUT']['ci'])} from `Z` to `HEAD_OUT` "
          f"(criterion: interval above 0.10: {yes(k['added_loss_above_margin'])}); `HEAD_OUT` rejected: {yes(k['head_out_rejected'])}; "
          f"`Z` not rejected: {yes(k['z_not_rejected'])}. Treated side resampled by {k['resampled']['treated']}.")
    a("")
    a("## Knockouts: what is missing, and the repairs that were already made\n")
    a("The scopes at `Z`, the stage the head reads:\n")
    a("| cell type | the head as it is | readout of the head's training cells | of the other cell type | of the other genes | with the gene's own cells | reading |")
    a("|---|---|---|---|---|---|---|")
    for ct, v in G3.items():
        if not isinstance(v, dict):
            continue
        s = v["scopes_at_Z"]
        a(f"| {ct} | {uc(v['frozen'])} | {uc(s['train'])} | {uc(s['other'])} | {uc(s['units'])} | {uc(s['own'])} | **{v['reading_of_the_grid']}**"
          + (f"; resolved gaps: {', '.join(v['resolved_gaps'])}" if v["resolved_gaps"] else "") + " |")
    a("")
    a("A hypothesis was written from each grid and read against the recorded audits of the repaired heads. It is retrospective: "
      "the outcomes were known. The rule (fall above 0.20 for the repair with cells of this kind, below 0.20 for the repair with "
      "the other cell type) was fixed in the design.\n")
    a("| cell type | arm | forecast from the grid | measured fall of \\|u\\| at `HEAD_OUT` | inside the forecast | outcome |")
    a("|---|---|---|---|---|---|")
    for ct, v in G3.items():
        if not isinstance(v, dict):
            continue
        for arm in ("coverage", "other study", "function"):
            fc = v["forecast"][arm]
            m = v["arms"].get(arm)
            a(f"| {ct} | {arm} | {f(fc['point'])} {iv(fc['ci'])} | " + (f"{f(m['gain'])} {iv(m['ci'])}" if m else "no head was retrained on bone marrow alone")
              + f" | {('yes' if m['inside_forecast'] else 'no') if m else '—'} | {v['outcome'] if arm == 'coverage' else ''} |")
    a("")
    fr = SUM["ledger"]["forecast_record"]
    a(f"Forecast record: {fr['inside']} of {fr['n']} measured falls inside their forecast interval; mean absolute error {f(fr['mean_abs_error'])}.\n")
    a("## Donors: the same tool on a contrast nobody assigned\n")
    for name, title in (("donors_hao", "External PBMC, 8 donors"), ("donors_bm", "Bone marrow, 12 donor × site groups")):
        d = G2["donors"][name]
        a(f"**{title}**, variation target (the head's ordinary score) at `Z`: the gap from the head to a readout of its training "
          f"cells is smaller than the gap from there to the other donors in {d['a_head_not_the_bottleneck']} of {d['n']}; the gap to the "
          f"other donors is positive in {d['b_coverage_positive']} of {d['n']}; the gap from the other donors to the donor's own cells "
          f"is smaller in {d['c_own_smaller']} of {d['n']}.\n")
    a("| donor | head → training cells | training cells → other donors | other donors → own cells |")
    a("|---|---|---|---|")
    for r in G2["donors"]["donors_hao"]["rows"]:
        a(f"| {r['unit'].replace('hao2021_', '')} | {f(r['frozen->train'], 3)} | {f(r['train->units'], 3)} | {f(r['units->own'], 3)} |")
    a("")
    a("Gaps are differences of *u* = 1 − mean *r*² within the donor, so they are not on the scale of the correlation differences "
      "reported in `runs/generalization_20261112`.\n")
    for name, title in (("donors_hao", "external donors"), ("donors_bm_unseen", "bone-marrow groups not trained on")):
        for tg in ("variation", "shift"):
            m = SUM["map"][name][tg]
            a(f"- {title}, {tg}: read as " + ", ".join(f"{k} in {v}" for k, v in sorted(m["missing"].items(), key=lambda kv: -kv[1])) + f" of {m['n_events']}.")
    a("")
    a("## The correction to `runs/generalization_20261112`\n")
    a("That run reported how much of a donor's shift in protein level is left unexplained, as a share of squared differences and, "
      "for the external donors, without a slope between the head's scale and the measurement's. The values on the scale of *u* "
      "(the square root of the old ones) and the values of version 2 (with the slope):\n")
    a("| donor | the head: old, on the scale of u | the head: version 2 | state, other donors: old | version 2 | state, training cells: old | version 2 |")
    a("|---|---|---|---|---|---|---|")
    for c in SUM["correction"]:
        if c["study"] != "donors_hao":
            continue
        o, n = c["old_on_the_scale_of_u"], c["new_u"]
        a(f"| {c['unit'].replace('hao2021_', '')} | {f(o['head'])} | {uc(n['head'])} | {f(o['state_other_donors'])} | {uc(n['state_other_donors'])} | "
          f"{f(o['state_training_cells'])} | {uc(n['state_training_cells'])} |")
    a("")
    ext = [c for c in SUM["correction"] if c["study"] == "donors_hao"]
    rg = lambda k, src: f"{f(min(c[src][k] if src == 'old_on_the_scale_of_u' else c[src][k]['u'] for c in ext))} to " \
                        f"{f(max(c[src][k] if src == 'old_on_the_scale_of_u' else c[src][k]['u'] for c in ext))}"
    a(f"On the external donors the head leaves {rg('head', 'new_u')} of a donor's shift unexplained (the old run said 64% to 86% of "
      f"squared differences, {rg('head', 'old_on_the_scale_of_u')} on this scale). The state, read by the other donors, leaves "
      f"{rg('state_other_donors', 'new_u')}; by a readout of the head's training cells, {rg('state_training_cells', 'new_u')}. "
      f"Supported at the state by the other donors in {sum(c['new_u']['state_other_donors']['verdict'] == 'support' for c in ext)} of 8.\n")
    a("## An open hypothesis for the donors\n")
    dh = SUM["donor_hypothesis"]
    a(f"Written from the external-donor grid, not retrospective, and not tested: **coverage**. Forecast for the fall of *u* "
      f"(variation) at the head's own output on donors it was not trained on: {f(dh['forecast']['coverage']['mean'], 3)} "
      f"{iv(dh['forecast']['coverage']['ci'], 3)} if PBMC cells of other donors are added to its training; "
      f"{f(dh['forecast']['function']['mean'], 3)} {iv(dh['forecast']['function']['ci'], 3)} for another map on the same bone-marrow "
      f"cells. Rule: confirmed if the first lies above {f(dh['rule']['at_least'], 3)} and the others below it. The intervals come "
      f"from resampling 8 donors and are approximate.\n")
    a("## The ledger\n")
    lg = SUM["ledger"]
    a(f"[`results/ledger.jsonl`](results/ledger.jsonl): {lg['events']} typed events ("
      + ", ".join(f"{v} {k}" for k, v in sorted(lg["by_type"].items())) + f"). Opposite verdicts at the same readout, stage, scope and target: "
      f"{lg['conflicts']}. Ledger rules: " + ("; ".join(f"{c['code']} ({c['severity']}): {c['message']}" for c in lg["checks"]) or "nothing to report") + ".\n")
    a("## Limits\n")
    for s in ("A re-reading of finished work. The hypotheses for the knockouts were written after their outcomes were known; "
              "agreement is a check of consistency, not a confirmation.",
              "The measured falls use version 1's estimator and the forecasts use version 2's. They are the same quantity estimated "
              "in two ways.",
              "The `other study` arm's interval takes its baseline at a point value.",
              "Five genes of one pathway; eight external donors of one study, all healthy; one head, one training seed.",
              "A donor contrast is confounded with site, sample and run in bone marrow.",
              "Forecasts come from ridge and boosting readouts. A trained head can do better or worse."):
        a(f"- {s}")
    a("")
    a("## Files\n")
    a("- [`DESIGN.md`](DESIGN.md); [`scripts/build.py`](scripts/build.py) (arrays and SPEC files), [`scripts/run.py`](scripts/run.py) "
      "(the six runs), [`scripts/read.py`](scripts/read.py) (criteria, hypotheses, ledger), [`scripts/write_readme.py`](scripts/write_readme.py)")
    a("- `results/<name>/`: `SPEC.yaml`, `scope.json`, `scope_events.jsonl` for each of the six runs (the resamples are not committed)")
    a("- [`results/summary.json`](results/summary.json), [`results/ledger.jsonl`](results/ledger.jsonl)")
    (RUN / "README.md").write_text("\n".join(L) + "\n")
    print("wrote README.md and the figure")


if __name__ == "__main__":
    main()
