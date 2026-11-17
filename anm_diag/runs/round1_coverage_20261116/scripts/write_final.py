#!/usr/bin/env python3
"""The part of README.md about gate 5, and its figure, from results/final.json, results/gate5_declaration.json,
results/final_selection.json and, when it is there, results/what_is_left_final.json. Nothing is fitted or scored here.
Used by write_readme.py; run alone it writes the figure and prints the section.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RUN = Path(__file__).resolve().parents[1]
RES = RUN / "results"
INK, MUTED, GRID = "#1b1f24", "#5d6672", "#e6e9ee"
HUE = {"claim": "#2f6fb5", "dose": "#c9692a", "control": "#7b5ea7"}
ARMS = [("coverage", "20,000 PBMC cells of the 8 development donors added", "claim", "coverage_not_seen"),
        ("amount", "control: 20,000 fewer bone-marrow cells, nothing added", "control", "amount"),
        ("function", "control: a wider head, twice the epochs", "control", "function"),
        ("labels", "control: the PBMC proteins shuffled within the donor", "control", "labels")]


def f(x, nd=3, sign=False):
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    if float(s) == 0:
        s = s.lstrip("+-")
    return s.replace("-", "−")


def iv(c, nd=3):
    return f"[{f(c[0], nd)}, {f(c[1], nd)}]"


def have():
    return (RES / "final.json").exists()


def load():
    F = json.loads((RES / "final.json").read_text())
    decl = json.loads((RES / "gate5_declaration.json").read_text())
    sel = json.loads((RES / "final_selection.json").read_text())
    dev = json.loads((RES / "development.json").read_text())
    left = json.loads((RES / "what_is_left_final.json").read_text()) if (RES / "what_is_left_final.json").exists() else None
    return F, decl, sel, dev, left


def clean(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#9aa1ab")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7, color="#9aa1ab", labelcolor=INK)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def figure(path, F, decl, dev):
    G, bound = F["gains_over_as_is"], decl["hypotheses"][0]["rule"]["at_least"]
    fig = plt.figure(figsize=(11.2, 4.5), dpi=200)
    a = fig.add_axes([0.262, 0.2, 0.3, 0.63])
    b = fig.add_axes([0.745, 0.2, 0.24, 0.63])
    # ---- a: the arms on the final study
    lo, hi = -0.02, 0.09
    a.axvline(0, color="#9aa1ab", lw=0.9, zorder=1)
    a.axvline(bound, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    for y, (k, name, cls, dk) in enumerate(ARMS):
        g, c, d = G[k], HUE[cls], dev["gains_over_as_is"][dk]["gain"]
        if lo <= d <= hi:
            a.scatter([d], [y - 0.3], marker="D", s=20, facecolor="white", edgecolor=MUTED, linewidth=0.9, zorder=3)
        if g["ci"][1] < lo:
            a.annotate("", xy=(lo + 0.001, y), xytext=(lo + 0.014, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.6, mutation_scale=9))
            a.text(lo + 0.0165, y, "off the scale", fontsize=6.4, color=MUTED, va="center", ha="left", bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=6)
        else:
            a.plot(np.clip(g["ci"], lo, hi), [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=3)
            a.scatter([g["gain"]], [y], s=34, color=c, edgecolor="white", linewidth=0.8, zorder=5)
            a.scatter([p["gain"] for p in g["per_fit"]], [y + 0.3] * len(g["per_fit"]), s=9, facecolor="white", edgecolor=c, linewidth=0.8, zorder=4)
        a.text(1.02, y, f"{f(g['gain'], 3, True)} {iv(g['ci'])}", fontsize=6.8, color=INK, va="center", ha="left", transform=a.get_yaxis_transform())
    a.set_yticks(range(len(ARMS)))
    a.set_yticklabels([n for _, n, _, _ in ARMS], fontsize=7)
    a.set_ylim(len(ARMS) - 0.4, -1.0)
    a.set_xlim(lo, hi)
    a.set_xticks([-0.02, 0, 0.02, 0.04, 0.06, 0.08])
    a.text(bound, -0.88, f"the rule's bound, {f(bound)}", fontsize=6.2, color=INK, ha="center", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    a.text(1.02, -0.88, "mean [interval]", fontsize=6.2, color=MUTED, ha="left", va="top", transform=a.get_yaxis_transform())
    a.set_xlabel("fall of u from the head fitted as registered, same seed (higher is better)", fontsize=7.2, color=INK)
    clean(a)
    # ---- b: the claim's arm, donor by donor
    st = F["status_of_donor"]
    donors = list(st)                                         # the order of the declaration, which is the order of per_unit
    gain = dict(zip(donors, G["coverage"]["per_unit"]))
    rng = np.random.default_rng(0)
    b.axvline(0, color="#9aa1ab", lw=0.9, zorder=1)
    b.axvline(bound, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    for y, (grp, key) in enumerate((("healthy", "healthy"), ("patient", "patients"))):
        v = np.array([gain[d] for d in donors if st[d] == grp])
        b.scatter(v, y + rng.uniform(-0.2, 0.2, len(v)), s=14, color=HUE["claim"], alpha=0.75, edgecolor="white", linewidth=0.4, zorder=3)
        hp = F["healthy_and_patients"]["coverage"][key]
        b.plot(hp["ci"], [y - 0.36, y - 0.36], color=INK, lw=1.6, solid_capstyle="round", zorder=4)
        b.scatter([hp["gain"]], [y - 0.36], s=22, color=INK, zorder=5)
    b.set_yticks([0, 1])
    b.set_yticklabels([f"healthy donors ({F['healthy']})", f"patients ({F['patients']})"], fontsize=7)
    b.set_ylim(1.55, -0.85)
    allv = np.array(G["coverage"]["per_unit"])
    b.set_xlim(min(-0.02, float(allv.min()) - 0.01), max(0.09, float(allv.max()) + 0.01))
    b.set_xlabel("fall of u with PBMC cells added\n(blue: a donor; black: the group's mean and interval)", fontsize=7.0, color=INK, linespacing=1.25)
    clean(b)
    for ax, x0, letter, title in ((a, 0.012, "a", "The arms on the final study"), (b, 0.63, "b", "The claim's arm in healthy donors and in patients")):
        y1 = ax.get_position().y1
        fig.text(x0, y1 + 0.05, letter, fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")
        fig.text(x0 + 0.017, y1 + 0.053, title, fontsize=8.2, fontweight="bold", color=INK, va="bottom", ha="left")
    fig.text(0.012, 0.02, f"Panel a. Filled mark and line: the mean over {F['donors']} donors and 3 seeds with its interval (1,000 resamples of donors, seeds resampled inside each arm). "
                          "Small open marks: each seed. Open diamond: the same arm in development.", fontsize=6.4, color=MUTED, ha="left", va="bottom")
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"))
    plt.close(fig)


def row_of_design(F):
    """The row of DESIGN.md's table of outcomes that the reading falls in, stated from the numbers alone."""
    o, hp = F["outcome"]["donors_hao_coverage"], F["healthy_and_patients"]["coverage"]["patients_minus_healthy"]
    t = F["tests"]["donors_hao_coverage"]["together"]
    if not F["guard_holds"]:
        return "guard broken"
    if o == "confirmed":
        if hp["ci"][1] < 0:
            return "confirmed among healthy, smaller among patients"
        if hp["ci"][0] > 0:
            return "confirmed; larger among patients than among healthy donors (a row the design did not list)"
        return "confirmed, guard holds, patients as healthy (the difference between the two groups is not resolved)"
    if o == "refuted":
        return "a control arm above its bound" if "control arm gains more" in t["why"] else "own arm below the bound, development had delivered"
    return "unresolved"


def lines(F, decl, sel, dev, left):
    L = []
    a = L.append
    G, T, hyp = F["gains_over_as_is"], F["tests"], {h["id"]: h for h in decl["hypotheses"]}
    bound = hyp["donors_hao_coverage"]["rule"]["at_least"]
    c = G["coverage"]
    a("## The test (gate 5), 17 November 2026\n")
    a(f"Read once, on a study nobody in this repository had looked at: {sel['name']}. {F['donors']} donors chosen by the rule of gate 0b "
      f"({F['healthy']} healthy, {F['patients']} patients), 1,000 cells each, {F['proteins']} proteins. Every step is recorded in [`GATE5.md`](GATE5.md), "
      "each part committed before the step it governs; the declaration ([`results/gate5_declaration.json`](results/gate5_declaration.json)) names every "
      "input by hash, and [`scripts/read_final.py`](scripts/read_final.py) refuses any other.\n")
    a("![the test](results/final.svg)\n")
    a(f"*u* is the share of the variation among a donor's cells that the head's output leaves unexplained: 1 minus the mean over the {F['proteins']} proteins of r² "
      "within the donor. A gain is the fall of *u* from the head fitted as registered (`as_is`) with the same seed, on the same cells. No head was fitted on a cell of this study.\n")
    a("### Outcome\n")
    for hid, own_txt in (("donors_hao_coverage", "`coverage` above the bound; `function` and `amount` below it"), ("coverage_not_levels", "`coverage` above the bound; `labels` below it")):
        t = T[hid]
        a(f"- **`{hid}`: {F['outcome'][hid]}.** The rule: {own_txt} ({f(bound)}). Over the three seeds together: {t['together']['outcome']} ({t['together']['why'].replace(str(bound), f(bound))}). "
          f"Each seed alone: {', '.join(x['outcome'] for x in t['by_seed'])}.")
    a(f"- **The guard {'holds' if F['guard_holds'] else 'does not hold'}**: on bone marrow, site 4, " + ("every head of the test is within 0.02 of `as_is` of its seed."
                                                                                                         if F["guard_holds"] else "a head of the test falls more than 0.02 below `as_is` of its seed."))
    adm = F["admission"]
    bad = [f"{k} ({s})" for k, v in adm.items() for s, ok in v.items() if not ok]
    a("- **Admission**: " + ("all fifteen heads are admitted." if not bad else "not admitted: " + ", ".join(bad) + "."))
    a(f"- In the table of outcomes of [`DESIGN.md`](DESIGN.md) this is the row: **{row_of_design(F)}**.\n")
    a(f"With PBMC cells of the development study in its training, the head's *u* on the donors of the final study goes from {f(F['u']['as_is']['mean'])} to "
      f"{f(F['u']['coverage']['mean'])}: a fall of {f(c['gain'])} {iv(c['ci'])} ({', '.join(f(p['gain']) for p in c['per_fit'])} by seed). "
      f"In the head's ordinary score, the mean correlation within a donor: {f(F['score_within_donor']['as_is']['mean'])} to {f(F['score_within_donor']['coverage']['mean'])}. "
      f"In development the same arm gained {f(dev['gains_over_as_is']['coverage_not_seen']['gain'])} on donors of the study it was fitted on.\n")
    a("### The arms on the final study\n")
    a(f"| arm | *u* | fall of *u* from `as_is` | by seed | read against {f(bound)}, by seed | forecast | ordinary score |")
    a("|---|---|---|---|---|---|---|")
    a(f"| `as_is` | {f(F['u']['as_is']['mean'])} | | | | | {f(F['score_within_donor']['as_is']['mean'])} |")
    fc = {**hyp["donors_hao_coverage"]["forecast"], "labels": hyp["coverage_not_levels"]["forecast"]["labels"]}
    for k, name, _, _ in ARMS:
        g = G[k]
        rd = ["above" if p["ci"][0] >= bound else ("below" if p["ci"][1] < bound else "crosses") for p in g["per_fit"]]
        fct = f"{f(fc[k][0])} to {f(fc[k][1])}" + (f"; from development {f(hyp['coverage_not_levels']['forecast']['coverage'][0])} to {f(hyp['coverage_not_levels']['forecast']['coverage'][1])}" if k == "coverage" else "")
        a(f"| `{k}`: {name.replace('control: ', '')} | {f(F['u'][k]['mean'])} | {f(g['gain'], 3, True)} {iv(g['ci'])} | {', '.join(f(p['gain'], 3, True) for p in g['per_fit'])} | "
          f"{', '.join(rd)} | {fct} | {f(F['score_within_donor'][k]['mean'])} |")
    a("")
    a(f"The head in use has *u* = {f(F['u']['head_in_use']['mean'])} on these cells. With {F['donors']} donors the intervals over donors are "
      + ("approximate.\n" if F["approximate"] else "not marked approximate (10 donors or more).\n"))
    a("### Healthy donors and patients\n")
    a("| arm | healthy donors | patients | patients minus healthy |")
    a("|---|---|---|---|")
    for k, _, _, _ in ARMS:
        h = F["healthy_and_patients"][k]
        a(f"| `{k}` | {f(h['healthy']['gain'], 3, True)} {iv(h['healthy']['ci'])} | {f(h['patients']['gain'], 3, True)} {iv(h['patients']['ci'])} | "
          f"{f(h['patients_minus_healthy']['gain'], 3, True)} {iv(h['patients_minus_healthy']['ci'])} |")
    a("")
    a("Donors are resampled inside each group. This reading was declared and does not change the outcome. Healthy donors and patients differ in more than the disease "
      "(centre, age, time of sampling); the difference is not an effect of the disease alone.\n")
    e = F["exact_name_proteins_only"]
    a("### Without the proteins matched under another name\n")
    a(f"Of the {F['proteins']} proteins, {e['proteins']} carry the same name in the study and in the head's panel; the others were matched as the same target under another name "
      f"(GATE5.md, part 2). On those {e['proteins']} alone: `coverage` {f(e['gains_over_as_is']['coverage']['gain'], 3, True)} {iv(e['gains_over_as_is']['coverage']['ci'])}, "
      f"`amount` {f(e['gains_over_as_is']['amount']['gain'], 3, True)} {iv(e['gains_over_as_is']['amount']['ci'])}, "
      f"`function` {f(e['gains_over_as_is']['function']['gain'], 3, True)} {iv(e['gains_over_as_is']['function']['ci'])}, "
      f"`labels` {f(e['gains_over_as_is']['labels']['gain'], 3, True)} {iv(e['gains_over_as_is']['labels']['ci'])}. Declared before the reading; not part of the outcome.\n")
    fr = F["forecasts"]
    a("### How the forecasts held\n")
    a(f"{fr['inside']} of {fr['n']} measured gains fall inside their forecast.\n")
    a("| hypothesis | arm | forecast | measured | inside |")
    a("|---|---|---|---|---|")
    for r in fr["rows"]:
        a(f"| `{r['hypothesis']}` | `{r['arm']}` | {f(r['forecast'][0])} to {f(r['forecast'][1])} | {f(r['gain'], 3, True)} | {'yes' if r['inside_forecast'] else 'no'} |")
    a("")
    new = {n: t for n, t in F["heads"].items() if t["fitted"] == "gate 5"}
    a("### The two arms fitted for the test\n")
    a("| head | checkpoint epoch | validation score, last 5 epochs minus the 5 before | own objective | validation donor against `as_is` | bone marrow, site 4, against `as_is` | admitted |")
    a("|---|---|---|---|---|---|---|")
    for n, t in sorted(new.items()):
        g3, g6 = t["gate3"], t["gate6"]
        own = f"{f(g3['own_objective']['gain'])} above `as_is` on its PBMC fit cells" if "own_objective" in g3 else "—"
        a(f"| `{n}` | {t['checkpoint_epoch']} | {f(g3['ran']['last5'] - g3['ran']['before5'], 4, True)} | {own} | {f(g3['guard_validation']['change'], 4, True)} | "
          f"{f(g6['change'], 4, True)} | {'yes' if g3['admitted'] else 'no'} |")
    a("")
    a("No cell of the final study is in a file a fit reads: the fits read the bone-marrow pack and the pack of the development study, and the reading script checks "
      "that neither is a file of the final study.\n")
    if left:
        v = left["variation"]
        ce, lf = v["ceiling_as_is"], v["left_after_the_change"]
        a("### What is left on the final study (the grid, read after the test)\n")
        a(f"Each donor against the other donors of the final study, for `as_is` and `coverage` of seed {left['seeds'][0]} "
          "([`results/what_is_left_final.json`](results/what_is_left_final.json)). This is for a next round and changes no outcome.\n")
        a("| *u* of `variation`, mean over the donors | `as_is` | `coverage` |")
        a("|---|---|---|")
        for key, name in (("frozen", "the head's own output"), ("train@Z", "a readout of the head's training cells, at the state"),
                          ("units@RNA512", "a readout of the other donors of the final study, on 512 RNA components"),
                          ("units@Z", "a readout of the other donors of the final study, at the state"), ("own@Z", "with the donor's own cells, at the state")):
            a(f"| {name} | {f(v['u']['as_is'][key])} | {f(v['u']['coverage'][key])} |")
        a("")
        share = c["gain"] / ce["all"]["gain"] if ce["all"]["gain"] > 0 else float("nan")
        a(f"- **The ceiling of the final study**: a readout of its other donors reads {f(ce['all']['gain'])} {iv(ce['all']['ci'])} more than `as_is` "
          f"(healthy donors {f(ce['healthy']['gain'])}, patients {f(ce['patients']['gain'])}). The fall measured in the test, {f(c['gain'])}, is {round(100 * share)}% of it.")
        a(f"- **Left after the change**: {f(lf['all']['gain'])} {iv(lf['all']['ci'])} between the changed head and that readout "
          f"(healthy donors {f(lf['healthy']['gain'])}, patients {f(lf['patients']['gain'])}).")
        sr = v["state_against_rna"]["as_is"]
        a(f"- **The state against the RNA**: at the state the readout of the other donors leaves {f(-sr['gain'], 3, True)} {iv([-sr['ci'][1], -sr['ci'][0]])} "
          "more unexplained than on 512 RNA components (negative: the state carries more).")
        oc = v["own_cells_after_other_donors"]
        a(f"- **The donor's own cells** add {f(oc['gain'])} {iv(oc['ci'])} after the other donors.")
        if left.get("disease"):
            d = left["disease"]
            parts = []
            for h in ("as_is", "coverage"):
                if h in d and "shift" in d[h]:
                    u = d[h]["shift"]["u"]
                    parts.append(f"`{h}`: the head's own output leaves {f(u['frozen'], 2)} of it unexplained, a readout of its training cells {f(u.get('train@Z', float('nan')), 2)}, "
                                 f"a readout crossed over donors {f(u.get('units@Z', float('nan')), 2)}")
            if parts:
                a("- **The disease as an event** (patients against healthy donors, the shift in protein level): " + "; ".join(parts) + ". The two groups differ in centre and "
                  "sampling as well, and the reading is of one seed.")
        a("")
    a("### What the test does not show\n")
    for t in (f"One study. From the development study to this one the laboratory, the chemistry of the assay (3′ to 5′), the antibody panel, the processing and the disease "
              "change together; the test does not say which of them the gain survives, only that it survives all of them at once.",
              f"The rule is about the gain against a bound. It is not a statement that the head is good on this study: *u* stays at {f(F['u']['coverage']['mean'])} with the change.",
              f"{F['proteins']} of the head's proteins, {F['proteins'] - F['exact_name_proteins_only']['proteins']} of them matched under another name of the same target; clones were not compared.",
              f"{F['donors']} donors taken by a rule fixed before the file was opened, from donors with one status, one centre and at least 1,000 cells in their first sample. "
              f"They are not a random sample of the study; the {len(sel.get('reserve', []))} others stay unopened.",
              "Healthy donors against patients is not the effect of the disease: the two groups differ in centre and sampling.",
              "TEDDY's pretraining data may have held the RNA of this study (not checked). \"Not seen\" is a statement about the head, the same for every arm.",
              "The gain is in the variation among a donor's cells. A donor's shift in level against other donors is not part of the outcome.",
              "The grid on the final study is of one seed, was read after the test and decides nothing."):
        a(f"- {t}")
    a("")
    a("Departures from what was written before each step are listed in [`GATE5.md`](GATE5.md), parts 3 and 4.\n")
    return L


def main():
    F, decl, sel, dev, left = load()
    figure(RES / "final.svg", F, decl, dev)
    print("\n".join(lines(F, decl, sel, dev, left)))


if __name__ == "__main__":
    main()
