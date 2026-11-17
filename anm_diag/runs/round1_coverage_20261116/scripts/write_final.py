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
    a = fig.add_axes([0.262, 0.22, 0.3, 0.61])
    b = fig.add_axes([0.745, 0.22, 0.24, 0.61])
    # ---- a: the arms on the final study
    dv = [dev["gains_over_as_is"][dk]["gain"] for _, _, _, dk in ARMS]
    lo = np.floor((min(min(G[k]["ci"][0] for k, _, _, _ in ARMS), min(p["gain"] for k, _, _, _ in ARMS for p in G[k]["per_fit"])) - 0.004) / 0.02) * 0.02
    hi = np.ceil((max(bound, max(G[k]["ci"][1] for k, _, _, _ in ARMS), max(d for d in dv if d < 0.2)) + 0.004) / 0.02) * 0.02
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
    a.set_xticks([round(x, 2) for x in np.arange(lo, hi + 1e-9, 0.02)])
    a.text(bound, -0.88, f"the rule's bound, {f(bound)}", fontsize=6.2, color=INK, ha="center", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    a.text(1.02, -0.88, "mean [interval]", fontsize=6.2, color=MUTED, ha="left", va="top", transform=a.get_yaxis_transform())
    a.set_xlabel("fall of u from the head fitted as registered, same seed (above 0: better than as_is)", fontsize=7.2, color=INK)
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
    b.set_xlim(min(-0.02, float(allv.min()) - 0.008), max(bound + 0.012, float(allv.max()) + 0.008))
    b.set_xlabel("fall of u with PBMC cells added\n(blue: a donor; black: the group's mean and interval)", fontsize=7.0, color=INK, linespacing=1.25)
    clean(b)
    for ax, x0, letter, title in ((a, 0.012, "a", "The arms on the final study"), (b, 0.63, "b", "The claim's arm in healthy donors and in patients")):
        y1 = ax.get_position().y1
        fig.text(x0, y1 + 0.05, letter, fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")
        fig.text(x0 + 0.017, y1 + 0.053, title, fontsize=8.2, fontweight="bold", color=INK, va="bottom", ha="left")
    fig.text(0.012, 0.02, f"Panel a. Filled mark and line: the mean over {F['donors']} donors and 3 seeds with its interval (1,000 resamples of donors, seeds resampled inside each arm). "
                          "Small open marks: each seed.\nOpen diamond: the same arm in development" + ("" if all(lo <= d <= hi for d in dv) else " (not drawn where it is off the scale)")
                          + ". Dashed line: the bound of the rule. Panel b. Donors are resampled inside each group.",
             fontsize=6.4, color=MUTED, ha="left", va="bottom", linespacing=1.4)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"))
    plt.close(fig)


def reading_of_seeds(F, hid):
    """What the rule reads in the seeds whose arms are all admitted, when the seeds together are not testable."""
    t = F["tests"][hid]
    return sorted({x["outcome"] for x in t["by_seed"] if x["outcome"] != "not testable"})


def row_of_design(F):
    """The row of DESIGN.md's table of outcomes, stated from the numbers alone. The table has no row for a head that
    is not admitted, so the row is given as the rule states it and, beside it, as the admitted seeds read."""
    o, hp = F["outcome"]["donors_hao_coverage"], F["healthy_and_patients"]["coverage"]["patients_minus_healthy"]
    t = F["tests"]["donors_hao_coverage"]["together"]
    if o == "confirmed":
        if hp["ci"][1] < 0:
            return "confirmed among healthy, smaller among patients"
        if hp["ci"][0] > 0:
            return "confirmed; larger among patients than among healthy donors (a row the design did not list)"
        return "confirmed, guard holds, patients as healthy (the difference between the two groups is not resolved)"
    if o == "refuted":
        return "a control arm above its bound" if "control arm gains more" in t["why"] else "own arm below the bound, development had delivered"
    return "unresolved"


def guard_text(F):
    bad = sorted((n, t["gate6"]["change"], t["arm"]) for n, t in F["heads"].items() if not t["gate6"]["pass"])
    cov = [t["gate6"]["change"] for n, t in sorted(F["heads"].items()) if t["arm"] == "coverage"]
    if not bad:
        return "- **The guard holds**: on bone marrow, site 4, every head of the test is within 0.02 of `as_is` of its seed."
    own = [b for b in bad if b[2] == "coverage"]
    s = ("- **The guard, stated over every head of the test, does not hold**: on bone marrow, site 4, "
         + "; ".join(f"`{n}` ({'the claim’s arm' if a == 'coverage' else 'a control'}) scores {f(-c)} below `as_is` of its seed" for n, c, a in bad) + ", against a limit of 0.02. ")
    s += ("The three `coverage` heads are within it (" + ", ".join(f(c, 3, True) for c in cov) + ")." if not own else "")
    return s


def post_hoc_contrasts(F):
    """Contrasts read after the test from the per-donor values of results/final.json, with the function the reading
    script uses (`scope.arms_family`). Written to results/final_contrasts.json; they change no outcome."""
    import sys
    sys.path.insert(0, str(RUN.parents[1]))
    from anm_diag import scope
    donors = list(F["status_of_donor"])
    seeds = sorted(F["per_donor"][donors[0]])
    U = {k: np.array([[F["per_donor"][d][s][k]["u"] for s in seeds] for d in donors]) for k in ("as_is", "coverage", "amount", "labels", "function")}
    g = lambda base, arm: {x: scope.arms_family({"base": U[base], "arm": U[arm]}, base="base", bootstrap=F["decision"]["bootstrap"], seed=F["decision"]["seed"],
                                                alpha=F["decision"]["alpha"])["arms"]["arm"][x] for x in ("gain", "ci")}
    out = {"what": "contrasts computed after the test by scripts/write_final.py with scope.arms_family, from the per-donor values of final.json; not in read_final.py",
           "coverage_against_amount": g("amount", "coverage"), "labels_against_coverage": g("coverage", "labels")}
    (RES / "final_contrasts.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def headline(F, dev, left=None):
    """One sentence for the top of the README, from the numbers."""
    c, d = F["gains_over_as_is"]["coverage"], dev["gains_over_as_is"]["coverage_not_seen"]
    o = F["outcome"]["donors_hao_coverage"]
    if o == "unresolved" and reading_of_seeds(F, "donors_hao_coverage") == ["refuted"]:
        s = (f"**The test did not bear the round out.** The gain read in development ({f(d['gain'], 3, True)}) is not there on the study of the test "
             f"({f(c['gain'], 3, True)} {iv(c['ci'])}). By the rule written first the outcome is *unresolved*, because one of the three heads of the claim is not admitted; "
             "the two seeds that can be read are each *refuted*.")
        if left:
            ce = left["variation"]["ceiling_as_is"]["all"]
            s += (f" **The grid read afterwards shows that the bound could not have been met on this study**: a readout fitted on its own donors reads only {f(ce['gain'])} {iv(ce['ci'])} "
                  f"more than the head, against {f(dev['D1_batch_seen']['of_the_ceiling']['ceiling'])} on the development study. On this study the head is not short of coverage.")
        return s
    return f"**The test: {o}.** The gain of the claim's arm is {f(c['gain'], 3, True)} {iv(c['ci'])}; in development it was {f(d['gain'], 3, True)}."


def lines(F, decl, sel, dev, left):
    L = []
    a = L.append
    G, T, hyp = F["gains_over_as_is"], F["tests"], {h["id"]: h for h in decl["hypotheses"]}
    bound = hyp["donors_hao_coverage"]["rule"]["at_least"]
    c = G["coverage"]
    dv = dev["gains_over_as_is"]["coverage_not_seen"]
    pu = np.array(c["per_unit"])
    adm = F["admission"]
    bad = [(k, s) for k, v in adm.items() for s, ok in v.items() if not ok]
    a("## The test (gate 5), 17 November 2026\n")
    a(f"Read once, on a study nobody in this repository had looked at: {sel['name']}. {F['donors']} donors chosen by the rule of gate 0b "
      f"({F['healthy']} healthy, {F['patients']} patients), 1,000 cells each, {F['proteins']} proteins. Every step is recorded in [`GATE5.md`](GATE5.md), "
      "each part committed before the step it governs; the declaration ([`results/gate5_declaration.json`](results/gate5_declaration.json)) names every "
      "input by hash, and [`scripts/read_final.py`](scripts/read_final.py) refuses any other.\n")
    a("![the test](results/final.svg)\n")
    a(f"*u* is the share of the variation among a donor's cells that the head's output leaves unexplained: 1 minus the mean over the {F['proteins']} proteins of r² "
      "within the donor. A gain is the fall of *u* from the head fitted as registered (`as_is`) with the same seed, on the same cells. No head was fitted on a cell of this study.\n")
    a("### Outcome\n")
    worse = c["ci"][1] < 0
    a(f"- **The gain of development is not there on this study.** With the PBMC cells of the development study in its training, the head's *u* on the {F['donors']} donors goes from "
      f"{f(F['u']['as_is']['mean'])} to {f(F['u']['coverage']['mean'])}: a gain of {f(c['gain'], 3, True)} {iv(c['ci'])}"
      + (", below 0" if worse else "") + f"; by seed {', '.join(f(p['gain'], 3, True) for p in c['per_fit'])}. The bound of the rule is {f(bound)}; in development the same arm gained "
      f"{f(dv['gain'], 3, True)} {iv(dv['ci'])} on donors of the study it was fitted on. The head with PBMC cells is the worse of the two in {int((pu < 0).sum())} of the {len(pu)} donors. "
      f"In the head's ordinary score, the mean correlation within a donor: {f(F['score_within_donor']['as_is']['mean'])} against {f(F['score_within_donor']['coverage']['mean'])}.")
    for hid, own_txt in (("donors_hao_coverage", "`coverage` at or above the bound; `function` and `amount` below it"), ("coverage_not_levels", "`coverage` at or above the bound; `labels` below it")):
        t = T[hid]
        a(f"- **`{hid}`: {F['outcome'][hid]}, by the rule written first.** The rule: {own_txt} ({f(bound)}); the outcome over the three seeds together counts only if each seed alone reads the same. "
          f"Over the three seeds together: {t['together']['outcome']} ("
          + ("an arm of the rule is not admitted in one seed, below" if t["together"]["outcome"] == "not testable" else t["together"]["why"].replace(str(bound), f(bound)))
          + f"). Each seed alone: {', '.join(x['outcome'] for x in t['by_seed'])}.")
    if bad:
        for k, sd in bad:
            n = f"s{sd.replace('seed', '')}_{ {'coverage': 'coverage_all', 'labels': 'labels_all'}.get(k, k)}"
            g3 = F["heads"][n]["gate3"]
            why = [x for x in ("ran", "own_objective", "no_scored_cell_in_fit", "guard_validation") if x in g3 and not g3[x]["pass"]]
            txt = {"ran": f"its validation score falls at the end of training (last 5 epochs {f(g3['ran']['last5'])}, the 5 before {f(g3['ran']['before5'])}; the limit is a fall of 0.005)",
                   "own_objective": "it does not reach its own objective", "no_scored_cell_in_fit": "a scored cell is in its fit",
                   "guard_validation": "its validation donor falls more than 0.02 below `as_is`"}
            hist = RES / "gate5_fit_histories.json"
            more = ""
            if "ran" in why and hist.exists():
                h = json.loads(hist.read_text())["heads"][n]
                v = h["validation_score"]
                more = (f" The fall is in the last epoch alone: {f(v[-1])} after {f(np.mean(v[-6:-1]))} over the five epochs before it "
                        f"([`results/gate5_fit_histories.json`](results/gate5_fit_histories.json)). The head that is scored is the checkpoint of epoch {h['checkpoint_epoch']}. "
                        "The rule of admission is applied as written.")
            a(f"- **Admission: `{n}` is not admitted** (gate 3): " + "; ".join(txt[x] for x in why) + "." + more + " The other fourteen heads are admitted. A head that is not admitted was not tested, "
              "so the three seeds together are not testable and the round's outcome is not stated as refuted.")
    else:
        a("- **Admission**: all fifteen heads are admitted.")
    a(guard_text(F))
    rs = reading_of_seeds(F, "donors_hao_coverage")
    row = row_of_design(F)
    if row == "unresolved" and rs == ["refuted"]:
        a("- **In the table of outcomes of [`DESIGN.md`](DESIGN.md)** the rule gives the row *unresolved*. The table wrote that row for intervals that do not decide; here the intervals decide in every "
          "seed that can be read, and they read the row *own arm below the bound, development had delivered*. The table gives that row the meaning that what was learned belongs to the study. "
          + ("The grid, read after the test, does not support putting it that way: on this study there was next to nothing of that kind to gain (below). " if left else "")
          + "The step the table names after *unresolved*, the reserve donors read once, has not been taken: more donors do not change the admission of a head, and opening them is the owner's decision.\n")
    else:
        a(f"- In the table of outcomes of [`DESIGN.md`](DESIGN.md) this is the row: **{row}**.\n")
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
    a(f"The head in use has *u* = {f(F['u']['head_in_use']['mean'])} on these cells. Every head leaves about four fifths of the variation unexplained here, against "
      f"{f(dev['u']['as_is']['mean'])} for `as_is` on the development study. With {F['donors']} donors the intervals over donors are "
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
      f"(GATE5.md, part 2). On those {e['proteins']} alone *u* of `as_is` is {f(e['u']['as_is'])}, and the gains are: `coverage` {f(e['gains_over_as_is']['coverage']['gain'], 3, True)} {iv(e['gains_over_as_is']['coverage']['ci'])}, "
      f"`amount` {f(e['gains_over_as_is']['amount']['gain'], 3, True)} {iv(e['gains_over_as_is']['amount']['ci'])}, "
      f"`function` {f(e['gains_over_as_is']['function']['gain'], 3, True)} {iv(e['gains_over_as_is']['function']['ci'])}, "
      f"`labels` {f(e['gains_over_as_is']['labels']['gain'], 3, True)} {iv(e['gains_over_as_is']['labels']['ci'])}. Declared before the reading; not part of the outcome. "
      "The reading is the same without the matched names.\n")
    fr = F["forecasts"]
    a("### How the forecasts held\n")
    a(f"{fr['inside']} of {fr['n']} measured gains {'falls' if fr['inside'] == 1 else 'fall'} inside {'its' if fr['inside'] == 1 else 'their'} forecast.\n")
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
          f"{f(g6['change'], 4, True)} | {'yes' if g3['admitted'] else '**no**'} |")
    a("")
    a("No cell of the final study is in a file a fit reads: the fits read the bone-marrow pack and the pack of the development study, and the reading script checks "
      "that neither is a file of the final study.\n")
    if left:
        v = left["variation"]
        ce, lf = v["ceiling_as_is"], v["left_after_the_change"]
        a("### What is left on the final study (the grid, read after the test)\n")
        a(f"For each of the {left['donors']} donors, *u* under readouts that were not fitted on the donor, for `as_is` and `coverage` of seed {left['seeds'][0]} "
          "([`results/what_is_left_final.json`](results/what_is_left_final.json); what was read and what was not is in [`GATE5.md`](GATE5.md), part 5). "
          "This is for a next round and changes no outcome. The readouts fitted on cells of the final study do not depend on the head.\n")
        a("| *u* of `variation`, mean over the donors | `as_is` | `coverage` |")
        a("|---|---|---|")
        for key, name in (("frozen", "the head's own output"), ("train@Z", "a readout of the head's training cells, at the state"),
                          ("units@RNA512", "a readout of the other donors of the final study, on 512 RNA components"),
                          ("units@Z", "a readout of the other donors of the final study, at the state"), ("own@Z", "with the donor's own cells, at the state")):
            a(f"| {name} | {f(v['u']['as_is'][key])} | {f(v['u']['coverage'][key])} |")
        a("")
        dc = dev["D1_batch_seen"]["of_the_ceiling"]["ceiling"]
        a(f"- **On this study the head is not short of coverage.** A readout fitted on the other donors of the final study leaves {f(v['u']['as_is']['units@Z'])} unexplained and "
          f"`as_is` {f(v['u']['as_is']['frozen'])}: {f(ce['all']['gain'])} {iv(ce['all']['ci'])} between them (healthy donors {f(ce['healthy']['gain'])}, patients {f(ce['patients']['gain'])}), "
          f"against {f(dc)} on the development study. A readout of the head's own training cells, which are bone marrow, reads as much ({f(v['u']['as_is']['train@Z'])}).")
        own_gain = v["u"]["as_is"]["frozen"] - v["u"]["as_is"]["own@Z"]
        a(f"- **The bound was out of reach here.** The bound of the rule, {f(bound)}, was set from the development study. On the final study no readout of the state leaves {f(bound)} "
          f"less unexplained than `as_is`: the readout of the other donors gains {f(ce['all']['gain'])}" + (f" (its interval ends at {f(ce['all']['ci'][1])})" if ce["all"]["ci"][1] < bound else "")
          + f", and {f(own_gain)} with the donor's own cells in the fit. The test could show that the gain was absent. It could not have shown it present, and so it does not say "
          "whether coverage is what donors of other PBMC studies lack.")
        gf = RES / "grid_final" / f"s{left['seeds'][0]}_donors_variation.json"
        dead = ""
        if gf.exists():
            pd_ = json.loads(gf.read_text())["heads"]["as_is"]["donors"]
            n_dead = sum(x["frozen"]["u"] > 0.99 and x["own@Z"]["u"] > 0.99 for x in pd_.values())
            dead = f" In {n_dead} of the {len(pd_)} donors nothing is read at all (*u* above 0.99 under every readout)." if n_dead else ""
        a(f"- **Most of the variation is beyond every readout**: {f(v['u']['as_is']['own@Z'])} is left even with the donor's own cells in the fit." + dead
          + " Whether that is noise of the protein measurements of this study or signal that a linear readout of the state does not reach, these data do not say.")
        a(f"- **Left after the change**: {f(lf['all']['gain'])} {iv(lf['all']['ci'])} between the head with PBMC cells and the readout of the other donors "
          f"(healthy donors {f(lf['healthy']['gain'])}, patients {f(lf['patients']['gain'])}).")
        tr = v["fall_from_as_is_to_coverage"]["train@Z"]
        a(f"- **Adding the PBMC cells to a linear readout changes nothing here**: a readout of the head's training cells leaves {f(v['u']['as_is']['train@Z'])} for `as_is` (bone marrow) and "
          f"{f(v['u']['coverage']['train@Z'])} for `coverage` (bone marrow and the PBMC cells of the development study), a difference of {f(tr['gain'], 3, True)} {iv(tr['ci'])}.")
        hr = v["head_against_readout_of_its_training_cells"]
        a(f"- **The heads against a linear readout of their own training cells**: `as_is` {f(hr['as_is']['gain'], 3, True)} {iv(hr['as_is']['ci'])}, "
          f"`coverage` {f(hr['coverage']['gain'], 3, True)} {iv(hr['coverage']['ci'])} (negative: the head is the worse of the two). On the development study the head with PBMC cells was the better of the two.")
        if gf.exists():
            G_ = json.loads(gf.read_text())["heads"]
            grp = {}
            for d_ in decl["donors"]:
                grp.setdefault((d_["status"], d_["centre"]), []).append(d_["donor"])
            a("")
            a("| donors | n | `as_is`, own output | `coverage`, own output | readout of the other donors, at the state | with the donor's own cells |")
            a("|---|---|---|---|---|---|")
            for (st_, ce_), ds in sorted(grp.items()):
                m = lambda h, k: float(np.mean([G_[h]["donors"][x][k]["u"] for x in ds]))
                a(f"| {st_}, centre {ce_} | {len(ds)} | {f(m('as_is', 'frozen'))} | {f(m('coverage', 'frozen'))} | {f(m('as_is', 'units@Z'))} | {f(m('as_is', 'own@Z'))} |")
            a("")
            a("By status and centre (seed 0, means over donors; no interval is given for these groups). The loss of `coverage` sits where the table shows it; the grid does not say why.\n")
        sr = v["state_against_rna"]["as_is"]
        a(f"- **The state against the RNA**: at the state the readout of the other donors leaves {f(-sr['gain'], 3, True)} {iv([-sr['ci'][1], -sr['ci'][0]])} "
          "more unexplained than on 512 RNA components (negative: the state carries more).")
        oc = v["own_cells_after_other_donors"]
        a(f"- **The donor's own cells** add {f(oc['gain'])} {iv(oc['ci'])} after the other donors.")
        a(f"- **Not read on the final study**: " + "; ".join(left["not_read"]) + ".")
        a("")
        a("This grid was read after the test, on one seed, and is exploratory. It was not available when the bound was set: the final study was not to be opened before the test.\n")
    chk = RES / "final_pack_check.json"
    if chk.exists():
        k = json.loads(chk.read_text())
        fs, ds = k["final_study"], k["development_study"]
        both = [x for x in fs["pairs"] if x in ds["pairs"]]
        show = [x for x in ("CD3", "CD8", "CD20", "HLA-DR", "CD19", "CD16") if x in fs["pairs"]]
        a("### Checks made after the reading (post hoc)\n")
        a(f"- **Is the pack of the final study sound?** For {len(fs['pairs'])} proteins whose coding gene is plain, the correlation within a donor between the protein and the RNA of its own gene "
          f"([`results/final_pack_check.json`](results/final_pack_check.json), `scripts/check_final_pack.py`): " + ", ".join(f"{x} {f(fs['pairs'][x], 2)}" for x in show)
          + f". Rows and names are matched. Over the {len(both)} pairs measured in both studies the mean is {f(np.mean([fs['pairs'][x] for x in both]), 2)} on the final study and "
          f"{f(np.mean([ds['pairs'][x] for x in both]), 2)} on the development study, with the same code: protein and RNA agree less, pair by pair, in the final study.")
        a("- **Is it the matched names?** No: on the proteins with the same name in both panels the gain of `coverage` is the same (above).")
        pc = post_hoc_contrasts(F)
        ca, lc = pc["coverage_against_amount"], pc["labels_against_coverage"]
        a(f"- **How much of the loss is the bone-marrow cells that `coverage` gives up?** `coverage` is fitted with 20,000 fewer bone-marrow cells, as `amount` is. Against `amount`, "
          f"`coverage` gains {f(ca['gain'], 3, True)} {iv(ca['ci'])}: " + ("the added PBMC cells neither help nor cost beyond that. " if ca["ci"][0] < 0 < ca["ci"][1] else
                                                                           ("the added PBMC cells cost beyond that. " if ca["ci"][1] < 0 else "the added PBMC cells help beyond that. "))
          + f"Against `coverage`, shuffling the proteins of the added cells costs {f(-lc['gain'])} {iv([-lc['ci'][1], -lc['ci'][0]])} "
          "([`results/final_contrasts.json`](results/final_contrasts.json)).")
        a("- **Was the fit of the head that is not admitted disturbed?** It finished before any other process ran on the machine; see GATE5.md, part 4.\n")
    a("### What the test does not show\n")
    for t in ("One study. From the development study to this one the laboratory, the chemistry of the assay (3′ to 5′), the antibody panel, the processing and the disease "
              "change together; nothing here separates them.",
              ("It could not have borne the round out on this study: the grid above finds almost no shortfall of coverage there. It is a test of whether the gain is there, "
               "not of whether coverage is what donors of other PBMC studies lack." if left else
               "It does not show that added cells cannot help on this study: no cell of this study was added. What cells of the study itself can read is for the grid on the "
               "final study, which was not finished when this was written."),
              f"By the rule the outcome is unresolved, not refuted: one of the three `coverage` heads is not admitted. A further seed is not in the design and was not fitted.",
              f"{F['proteins']} of the head's proteins, {F['proteins'] - F['exact_name_proteins_only']['proteins']} of them matched under another name of the same target; clones were not compared.",
              f"{F['donors']} donors taken by a rule fixed before the file was opened, from donors with one status, one centre and at least 1,000 cells in their first sample. "
              f"They are not a random sample of the study; the {len(sel.get('reserve', []))} others stay unopened.",
              "Healthy donors against patients is not the effect of the disease: the two groups differ in centre and sampling.",
              "TEDDY's pretraining data may have held the RNA of this study (not checked). \"Not seen\" is a statement about the head, the same for every arm.",
              "The gain is in the variation among a donor's cells. A donor's shift in level against other donors is not part of the outcome.",
              "The grid on the final study is of one seed, is read after the test and decides nothing."):
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
