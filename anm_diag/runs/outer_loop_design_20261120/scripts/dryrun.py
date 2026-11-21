#!/usr/bin/env python3
"""A dry run of the outer loop's routing on the finished round of 20 November (runs/teg_loop_20261120).

  python dryrun.py --round ../../teg_loop_20261120/results --out ../results

No cell is read and nothing is fitted: the script reads the three result files of the round and applies the rules of
OUTER_LOOP.md to them. It writes

  open_items.json   the status of every target, every open item with its reading, lever and actions, and the counts
  open_items.md     the same as tables

The routing table (ROUTES) is the one printed in OUTER_LOOP.md; `--routes` prints it as markdown.
"""
import argparse
import json
import math
from pathlib import Path

MARGIN, STEP, MIN_SNR, MIN_UNITS = 0.10, 0.10, 3.0, 10

LEVERS = {"M": "the model (inner loop)", "E": "events: which units of the pool are read next, and in which role",
          "V": "views: the readers of a stage, the input, the readouts", "S": "schema: the fields of an event and what is pooled",
          "C": "claim: the targets, the grade and the tolerance asked for", "P": "protocol: the gates and the rules"}

# reading -> (lever, what the reading means, admissible actions in order of cost, what not to do)
ROUTES = {
    "invalid": ("V", "the reading cannot be used (Section 6)", ["repair the reading before any action, as Section 6 says for the check that failed"], "act on the reading"),
    "unresolved": ("E", "the interval does not decide between support and rejection",
                   ["if the units of the pool differ among themselves, the split comes first", "more independent units of the same type from the pool",
                    "if the pool cannot supply them, report the estimate as it is: at the margin"],
                   "change the model or the views"),
    "small response": ("E", "the response is too small against its interval to be followed",
                       ["more cells per event, or more events that move this readout", "otherwise the readout is not asked for this type"], "make it a target"),
    "unreadable": ("V", "a resolved loss at the input view, with the event's own cells, in a view that is not dominated",
                   ["a richer reader of the same input", "another input already measured on the same cells, never the readout itself", "a new measurement", "a stated limit of the claim"],
                   "train"),
    "upstream": ("M", "the input view reads it and the stage the head reads does not",
                 ["change what comes before the stage: an inner loop of another depth", "if that part cannot be changed, the item stays open under its owner"], "change the head"),
    "no rows": ("E", "cells of this kind are needed and the table has no admitted event of the type", ["events of the type from the pool"], "carry rows over from another type"),
    "pooled, and differs": ("S", "a distinction the schema pools is visible in the readings",
                            ["split along a field from the record, with rival fields, tested on new units"], "repair across the pool"),
    "pooled, may be the same": ("S", "a distinction the schema keeps looks invisible to the readout", ["merge after bounded support on new units, twice"], "merge on development readings"),
    "mistyped": ("S", "cells get a type from a reader of their state, and some get the wrong one", ["take the type from the record"], "tune the reader"),
    "self event": ("C", "the event acts on the readout's own gene: an RNA-only system may have nothing to read (rule S3)",
                   ["a stated limit for this type", "or the assignment as an input, which is another system: the owner's decision"], "use it as a row"),
    "rule": ("P", "a gate or a test could not be evaluated as written", ["amend the protocol before the next iteration"], "re-grade the finished run"),
    "selection": ("P", "an event was selected by a measurement on the cells that were then audited", ["a selection split"], "keep selecting on audited cells"),
    "no final data": ("E", "no unit of the claim is unread", ["a census of the pool; units set aside as final before anything else is read"], "call a development reading a confirmation"),
    "confounded": ("E", "two fields vary together in the units that were read", ["units of the pool in which they vary apart", "otherwise the claim names their product"],
                   "attribute the difference to one of them"),
    "guard conflict": ("C", "every change that removes the loss breaks the guard on the system's own task",
                       ["a change that holds the guard: the inner loop again", "if none does, the owner decides between the target and the guard"], "widen the guard after the result"),
}

# What the pool holds for each type, from the records of earlier runs (design labels and counts only; no response is read here).
POOL = {
    "melanoma|IFN-gamma": {"with_embeddings": "206 knockout targets. On the readout's own scale the other 200 genes change PD-L1 by 0.35 at most, median 0.04 "
                                               "(runs/teg_repair_test_20261119), where the five pathway genes change it by 0.69 to 0.96: none is an event of the same kind",
                           "more_units_of_the_type": {"PDL1": 0}},
    "melanoma|co-culture": {"with_embeddings": "5 knockout targets (6,000 cells)", "more_units_of_the_type": {"PDL1": None, "MHC1": None, "HLAE": None},
                            "note": "the rest of the co-culture screen has no embeddings yet"},
    "THP-1|stimulated": {"with_embeddings": "5 knockout targets (7,478 cells)", "more_units_of_the_type": {"PDL1": None},
                         "note": "the rest of the screen has no embeddings yet"},
}
TYPE = {"mel_ifn": ("melanoma", "IFN-gamma"), "mel_co": ("melanoma", "co-culture"), "thp1": ("THP-1", "stimulated")}
GROUPS_OTHER = {"mel_ifn": "mel_co", "mel_co": "mel_ifn", "thp1": "mel_ifn"}      # the other group of each group in the round (build.py)


def needed_units(u, ci, n, bound=MARGIN):
    """Units needed for the interval to leave the bound on the side of the estimate, if the estimate holds and the
    half-width falls as 1/sqrt(n). A planning number, not a promise."""
    q = abs(u)
    lo, hi = sorted((abs(ci[0]) if ci[0] * ci[1] > 0 else 0.0, max(abs(ci[0]), abs(ci[1]))))
    if q < bound:
        half, gap = hi - q, bound - q
    else:
        half, gap = q - lo, q - bound
    if half <= 0 or gap >= half:
        return n
    if gap <= 0:
        return None
    return int(math.ceil(n * (half / gap) ** 2))


NEEDS = {"confounded": "units of the pool", "guard conflict": "the owner's decision", "invalid": "no new data", "rule": "no new data", "selection": "no new data", "mistyped": "no new data", "unreadable": "no new data first; then the owner's decision",
         "self event": "the owner's decision", "small response": "units of the pool, or the owner's decision", "no final data": "units of the pool",
         "pooled, and differs": "units of the pool", "pooled, may be the same": "units of the pool", "no rows": "units of the pool", "upstream": "a change of the encoder"}


def item(items, target, reading, quantity, value=None, ci=None, n=None, evidence=None, note=None, event="pool", needs=None, question=None, standing=None):
    """One open item. Observations that bear on the same question are kept under one item."""
    lever, meaning, actions, dont = ROUTES[reading]
    r, g = target
    obs = {"readout": r, "group": g, "event": event, "quantity": quantity, "value": value, "ci": ci, "evidence": evidence, "note": note}
    if question:
        for it in items:
            if it.get("question") == question and it["reading"] == reading:
                if not any(x["group"] == g and x["value"] == value for x in it["observations"]):
                    it["observations"].append(obs)
                return it
    ct, st = TYPE.get(g, (None, None))
    have = (POOL.get(f"{ct}|{st}", {}).get("more_units_of_the_type", {}) or {}).get(r) if g else None
    it = {"id": f"O{len(items) + 1:02d}", "question": question, "readout": r, "group": g, "cell_type": ct, "stratum": st, "event": event, "reading": reading, "lever": lever,
          "quantity": quantity, "value": value, "ci": ci, "resampling_units": n, "units_needed": needs, "more_units_in_pool": have,
          "next_step_needs": NEEDS.get(reading), "standing": standing, "actions": actions, "not": dont, "evidence": evidence, "note": note, "observations": [obs]}
    if reading == "unresolved" and needs is not None and n is not None:
        short = needs - n
        if have is not None and short > have:
            it["forecast"], it["next_step_needs"] = "at the margin: the pool has no further unit of the kind", "at the margin"
        elif short > 3 * n:
            it["forecast"], it["next_step_needs"] = "at the margin: more than four times the units would be needed", "at the margin"
        elif have is None:
            it["forecast"], it["next_step_needs"] = "within reach if the rest of the screen holds units of the type: a census first", "units of the pool"
        else:
            it["forecast"], it["next_step_needs"] = "within reach of the pool", "units of the pool"
    elif reading == "unresolved" and needs is None:
        it["forecast"], it["next_step_needs"] = "at the margin: the estimate sits on the bound", "at the margin"
    items.append(it)
    return it


def read_round(R, files):
    targets, items = [], []
    for fname, tag in files:
        D = json.loads((R / fname).read_text())
        for r, X in D.items():
            for g, e in X["groups"].items():
                adm, grid = e["admission"], e["grid"]
                ev = f"{fname}:{r}.groups.{g}"
                Z = grid["scopes"]["units"]
                n_ev = len(e["events"])
                if adm["diagnostic"] != "ready":
                    snr = grid["response"]["snr"]
                    n_ok = sum(x["admission"]["status"] == "admitted" for x in e["events"])
                    targets.append({"readout": r, "group": g, "part": tag, "status": "not a target", "why": adm["reasons"]})
                    item(items, (r, g), "small response", "the change of the readout against the half-width of its interval, pooled", snr, None, n_ev, ev + ".admission",
                         f"pooled {snr:.1f} (wanted {MIN_SNR:g}); events that pass on their own: {n_ok} of {n_ev}",
                         needs=None if snr >= MIN_SNR else int(math.ceil(n_ev * (MIN_SNR / snr) ** 2)))
                    continue
                rows = e["arms"]["rows"]
                beyond = rows["beyond_units_at_Z"]
                grade = e["closure"]["grade"]
                status = {"relative": "the head's part closed (relative)", "absolute": "closed (absolute)"}.get(grade, "not graded")
                targets.append({"readout": r, "group": g, "part": tag, "status": status, "attainable": adm["attainable"], "own_output": {"u": rows["u"], "ci": rows["ci"], "verdict": rows["verdict"]},
                                "added_beyond_Z": {"added": beyond["added"], "ci": beyond["ci"], "reading": beyond["reading"]}})
                # the rules that could not be evaluated
                if not e["gate3"]["objective_reached"] and beyond["reading"] == "no added loss":
                    item(items, (r, g), "rule", "gate 3: the change reaches its objective on the cells it was fitted on", e["gate3"]["on_the_cells_it_was_fitted_on"], None, None, ev + ".gate3",
                         "the check is absolute; the reference stage itself leaves " + f"{Z['u']:.2f}, so no change of the head can pass it. State it against the reference stage.")
                if e["test"]["outcome"] == "not testable":
                    item(items, (r, g), "rule", "the test of the hypothesis", None, None, None, ev + ".test", e["test"]["why"] + ": the other group has no admitted event, and the protocol names no other control")
                # the question of an absolute pass
                if adm["attainable"] == "absolute" and rows["verdict"] != "support":
                    item(items, (r, g), "unresolved", "u of the repaired output, against the margin", rows["u"], rows["ci"], n_ev, ev + ".arms.rows",
                         "the reference stage is supported: nothing is left before the head", needs=needed_units(rows["u"], rows["ci"], n_ev))
                # what is left before the head
                iv = (adm.get("input_view") or {}).get("pool")
                dominated = bool(iv and iv["added_from_RNA_to_Z_own_scope"]["ci"][1] < 0)
                if Z["verdict"] == "unresolved":
                    item(items, (r, g), "unresolved", "u at Z, read from the other genes of the group: is anything lost before the head", Z["u"], Z["ci"], n_ev, ev + ".grid.scopes.units",
                         needs=needed_units(Z["u"], Z["ci"], n_ev))
                elif Z["verdict"] == "reject":
                    item(items, (r, g), "unreadable", "u at Z and at the input view, with the event's own cells", grid["scopes"]["own"]["u"], grid["scopes"]["own"]["ci"], n_ev, ev + ".grid.scopes",
                         "one reader class (ridge) at both stages; the input view was read after the run (after_replication.json): a richer reader comes before any conclusion about RNA")
                if dominated:
                    a = iv["added_from_RNA_to_Z_own_scope"]
                    item(items, (r, g), "invalid", "from the input view to Z, own scope: a later stage reads better than an earlier one", a["added"], a["ci"], n_ev, ev + ".admission.input_view.pool",
                         f"the RNA view (512 components, ridge) leaves {iv['own']['u']:.2f} and Z, which is computed from the same RNA, leaves {grid['scopes']['own']['u']:.2f}: "
                         "the RNA view is the weaker reader, and its rejection says nothing about RNA")
                # events inside the pool that read differently from it
                for s, x in rows["events"].items():
                    if x["verdict"] == "reject" and x["ci"][0] > rows["ci"][1]:
                        item(items, (r, g), "pooled, and differs", "u of the repaired output for one event, against the pool's interval", x["u"], x["ci"], None, ev + f".arms.rows.events.{s}",
                             f"the pool reads {rows['u']:.2f} [{rows['ci'][0]:.2f}, {rows['ci'][1]:.2f}]; this event's interval lies above it" +
                             ("; for this event the loss is already at Z" if x["beyond_units_at_Z"]["reading"] == "no added loss" else ""), event=s,
                             question=f"one event inside the pool: {s}, {' '.join(TYPE[g])}")
                # what carries across types
                o, ro = e["arms"].get("other"), e["arms"].get("rows_other")
                if o and o["beyond_units_at_Z"]["reading"] == "loss":
                    b = o["beyond_units_at_Z"]
                    pair = " and ".join(sorted({" ".join(TYPE[g]), " ".join(TYPE[GROUPS_OTHER[g]])}))
                    item(items, (r, g), "pooled, and differs", "added beyond Z by the repair of the other group, applied as it is", b["added"], b["ci"], n_ev, ev + ".arms.other",
                         "the whole repair does not carry across the two types: the field that separates them stays in the schema", question=f"the whole repair across: {pair}",
                         standing="the schema already keeps them apart; seen on development cells, to be confirmed on new units")
                if ro:
                    b = ro["beyond_units_at_Z"]
                    if b["reading"] == "no added loss":
                        pair = " and ".join(sorted({" ".join(TYPE[g]), " ".join(TYPE[GROUPS_OTHER[g]])}))
                        item(items, (r, g), "pooled, may be the same", "added beyond Z with this type's control cells and the other type's rows", b["added"], b["ci"], n_ev, ev + ".arms.rows_other",
                             "a candidate for merging the rows of the two types; a merge needs bounded support on new units, twice", question=f"the rows across: {pair}")
                    elif b["reading"] == "unresolved":
                        item(items, (r, g), "unresolved", "added beyond Z with this type's control cells and the other type's rows: can the rows be merged", b["added"], b["ci"], n_ev,
                             ev + ".arms.rows_other", needs=needed_units(b["added"], b["ci"], n_ev, STEP),
                             question="the rows across: " + " and ".join(sorted({" ".join(TYPE[g]), " ".join(TYPE[GROUPS_OTHER[g]])})))
            if "probe" in X:
                p = X["probe"]
                item(items, (r, p["group"]), "self event", "u at Z for the knockout of the readout's own gene", p["grid"]["scopes"]["own"]["u"], p["grid"]["scopes"]["own"]["ci"], None,
                     f"{fname}:{r}.probe", "not admitted (A2); no reader of Z reads it", event=p["site"])
            if "gate" in X:
                worst = min((q for q in X["gate"] if "own_on_knockout_cells" in q), key=lambda q: q["own_on_knockout_cells"])
                if worst["own_on_knockout_cells"] < 0.99:
                    item(items, ("all", worst["group"]), "mistyped", "share of knockout cells the gate gives their own type", worst["own_on_knockout_cells"], None, None, f"{fname}:{r}.gate",
                         "the type is guessed from the state; the record has it", question="how a cell gets its type")
    return targets, items


def static_items(items):
    """Open items that follow from how the round was made, not from one of its numbers."""
    item(items, ("all", None), "selection", "A1 (the response against its interval) was read on the cells that were then audited", None, None, None, "runs/teg_loop_20261120/scripts/round.py: admit",
         "negligible where the response is ten times its half-width (PD-L1); material near the threshold of three, where no event was admitted")
    item(items, ("all", None), "rule", "a group ready as a pool with no admitted event", None, None, None, "runs/teg_loop_20261120/README.md",
         "resolved before the replication as 'not a target'; the protocol does not name the case")
    item(items, ("all", None), "no final data", "units of the knockout claim that no run has read", 0, None, None, "runs/teg_loop_20261120/DESIGN.md",
         "every cell with an embedding has been read; the rest of the two screens has no embeddings yet")


def md(targets, items, counts):
    f = lambda v: "" if v is None else (f"{v:.2f}".replace("-", "−") if isinstance(v, float) else str(v))
    ci = lambda c: "" if not c else f"[{f(c[0])}, {f(c[1])}]"
    o = ["## Targets at the end of the round", "", "| readout | group | part | status | own output after the repair | added beyond Z |", "|---|---|---|---|---|---|"]
    for t in targets:
        oo, b = t.get("own_output"), t.get("added_beyond_Z")
        o.append(f"| {t['readout']} | {t['group']} | {t['part']} | {t['status']} | {(f(oo['u']) + ' ' + ci(oo['ci']) + ' ' + oo['verdict']) if oo else ''} | {(f(b['added']) + ' ' + ci(b['ci']) + ' ' + b['reading']) if b else ''} |")
    o += ["", "## Open items", "", "| id | readout | group | event | reading | lever | quantity | value | units now → needed | the next step needs |", "|---|---|---|---|---|---|---|---|---|---|"]
    for i in items:
        need = "" if (i["resampling_units"] is None or i["reading"] not in ("unresolved", "small response")) else (
            f"{i['resampling_units']} → {'not reachable' if (i['reading'] == 'unresolved' and i['units_needed'] is None) else (f(i['units_needed']) or 'see note')}")
        val = "; ".join(f"{f(x['value'])} {ci(x['ci'])}".strip() + (f" ({x['group']})" if len(i["observations"]) > 1 else "") for x in i["observations"])
        o.append(f"| {i['id']} | {i['readout']} | {i['group'] or ''} | {i['event']} | {i['reading']} | {i['lever']} | {i['question'] or i['quantity']} | {val} | {need} | "
                 f"{i.get('forecast') or i['next_step_needs']} |")
    o += ["", "## Count by lever", "", "| lever | what it is | open items |", "|---|---|---|"]
    o += [f"| {k} | {LEVERS[k]} | {counts[k]} |" for k in LEVERS]
    return "\n".join(o) + "\n"


def routes_md():
    o = ["| reading | what it means | lever | actions, cheapest first | what not to do |", "|---|---|---|---|---|"]
    for k, (lever, meaning, actions, dont) in ROUTES.items():
        o.append(f"| {k} | {meaning} | {lever} | " + (actions[0] if len(actions) == 1 else "; ".join(f"{n + 1}. {a}" for n, a in enumerate(actions))) + f" | {dont} |")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--round", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--routes", action="store_true")
    a = ap.parse_args()
    if a.routes:
        print(routes_md())
        doc = Path(__file__).resolve().parents[3] / "OUTER_LOOP.md"
        missing = [ln for ln in routes_md().split("\n")[2:] if ln not in doc.read_text()]
        assert not missing, missing          # the table in OUTER_LOOP.md is this table
        return
    targets, items = read_round(a.round, [("round_development.json", "development"), ("round_replication.json", "replication")])
    static_items(items)
    counts = {k: sum(i["lever"] == k for i in items) for k in LEVERS}
    # the checks of the dry run: every target has one status, every item one reading that the table routes
    assert len({(t["readout"], t["group"]) for t in targets}) == len(targets)
    assert all(i["reading"] in ROUTES and i["lever"] in LEVERS and i["actions"] and i["next_step_needs"] for i in items)
    for n, i in enumerate(items):
        i["id"] = f"O{n + 1:02d}"
    for i in items:
        i["approximate"] = bool(i["resampling_units"] is not None and i["resampling_units"] < MIN_UNITS)
        if i["reading"] == "unresolved":
            i["after"] = [j["id"] for j in items if j["reading"] == "pooled, and differs" and j["group"] == i["group"] and j["readout"] == i["readout"] and j["event"] != "pool"]
    out = {"source": "runs/teg_loop_20261120/results", "margin": MARGIN, "step_margin": STEP, "levers": LEVERS,
           "routes": {k: {"lever": v[0], "meaning": v[1], "actions": v[2], "not": v[3]} for k, v in ROUTES.items()},
           "pool": POOL, "targets": targets, "items": items, "count_by_lever": counts,
           "count_by_reading": {k: sum(i["reading"] == k for i in items) for k in ROUTES},
           "count_by_what_the_next_step_needs": {k: sum(i["next_step_needs"] == k for i in items) for k in sorted({i["next_step_needs"] for i in items})},
           "checks": {"targets": len(targets), "every_target_has_one_status": True, "items": len(items), "every_item_is_routed": True,
                      "readings_never_met_in_this_round": [k for k in ROUTES if not any(i["reading"] == k for i in items)],
                      "items_whose_interval_is_approximate": [i["id"] for i in items if i["approximate"]],
                      "note": f"an interval over fewer than {MIN_UNITS} events is approximate (protocol v2); every pooled interval of this round rests on five"}}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "open_items.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    (a.out / "open_items.md").write_text(md(targets, items, counts))
    print(md(targets, items, counts))
    print(out["checks"])


if __name__ == "__main__":
    main()
