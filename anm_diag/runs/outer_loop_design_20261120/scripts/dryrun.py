#!/usr/bin/env python3
"""A dry run of the outer loop's routing on the finished round of 20 November (runs/teg_loop_20261120).

  python dryrun.py --round ../../teg_loop_20261120/results --out ../results
  python dryrun.py --routes        # the routing table; fails if OUTER_LOOP.md prints another one

No cell is read and nothing is fitted: the script turns the three result files of the round into the plain target
format of `anm_diag.outer`, and `outer.open_items` applies the rules of OUTER_LOOP.md. It writes

  open_items.json   the status of every target, every open item with its reading, lever and actions, and the counts
  open_items.md     the same as tables

If results/pool_screen.json exists (scripts/pool_screen.py), it is used for what the pool at hand can supply.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from anm_diag import outer  # noqa: E402

TYPE = {"mel_ifn": ("melanoma", "IFN-gamma"), "mel_co": ("melanoma", "co-culture"), "thp1": ("THP-1", "stimulated")}
OTHER = {"mel_ifn": "mel_co", "mel_co": "mel_ifn", "thp1": "mel_ifn"}      # the other group of each group in the round (build.py)
pair = lambda g: " and ".join(sorted({" ".join(TYPE[g]), " ".join(TYPE[OTHER[g]])}))
sub = lambda v, *keys: {k: v[k] for k in keys}


def targets_of_round(R, screen=None):
    """The round's result files in the target format of anm_diag.outer."""
    after = json.loads((R / "after_replication.json").read_text())
    targets, world = [], {"final_units": 0, "final_note": "every cell with an embedding has been read; the rest of the two screens has no embeddings yet",
                           "selection_on_audited_cells": "A1 (the response against its interval) was read on the cells that were then audited; negligible where the response is ten times "
                                                         "its half-width (PD-L1), material near the threshold of three, where no event was admitted",
                           "undeclared_cases": ["a group ready as a pool with no admitted event"]}
    typed = []
    for fname, part in (("round_development.json", "development"), ("round_replication.json", "replication")):
        D = json.loads((R / fname).read_text())
        for r, X in D.items():
            typed += [(q["own_on_knockout_cells"], q["group"]) for q in X.get("gate", []) if "own_on_knockout_cells" in q]
            for g, e in X["groups"].items():
                adm, grid, ev = e["admission"], e["grid"], f"{fname}:{r}.groups.{g}"
                now = (screen or {}).get("readouts", {}).get(r, {}).get("admissible_not_in_the_round") if g == "mel_ifn" else None
                t = {"readout": r, "group": g, "part": part, "type": dict(zip(("cell_type", "stratum"), TYPE[g])), "events_n": len(e["events"]), "evidence": ev + ".admission",
                     "admitted": adm["diagnostic"] == "ready", "response": {"snr": grid["response"]["snr"]},
                     "pool": {"more_units_of_the_type": None if now is None else len(now), "admissible_now": now}}
                if not t["admitted"]:
                    n_ok = sum(x["admission"]["status"] == "admitted" for x in e["events"])
                    t["note"] = f"pooled {grid['response']['snr']:.1f} (wanted {outer.DEFAULTS['min_snr']:g}); events that pass on their own: {n_ok} of {len(e['events'])}"
                    targets.append(t)
                    continue
                rows = e["arms"]["rows"]
                t["attainable"] = adm["attainable"]
                t["reference"] = {"units": grid["scopes"]["units"], "own": grid["scopes"]["own"], "evidence": ev + ".grid.scopes"}
                t["model"] = {"own": sub(rows, "u", "ci", "verdict"), "added": rows["beyond_units_at_Z"], "reading": None, "rows_of_type": len(e["events"]), "evidence": ev + ".arms.rows",
                              "events": {s: {**sub(x, "u", "ci", "verdict"), "added": x["beyond_units_at_Z"]} for s, x in rows["events"].items()}}
                iv = (adm.get("input_view") or {}).get("pool") or (after.get(f"{g}_{r}") or {}).get("pool")
                if iv:
                    own = iv["own"] if "own" in iv else iv["from_RNA"]["own"]
                    z_own = grid["scopes"]["own"]["u"]
                    t["input_view"] = {"own": own, "gap_to_later": iv["added_from_RNA_to_Z_own_scope"], "classes": ["ridge"], "nonlinear": False, "evidence": ev + ".admission.input_view",
                                       "note": f"the RNA view (512 components, ridge) leaves {own['u']:.2f} and Z, which is computed from the same RNA, leaves {z_own:.2f}"}
                why = []
                if not e["gate3"]["objective_reached"] and rows["beyond_units_at_Z"]["reading"] == "no added loss":
                    why.append({"question": "gate 3 is absolute: the change reaches its objective on the cells it was fitted on", "value": e["gate3"]["on_the_cells_it_was_fitted_on"], "evidence": ev + ".gate3",
                                "note": f"the reference stage itself leaves {grid['scopes']['units']['u']:.2f}, so no change of the head can pass it; state it against the reference stage"})
                if e["test"]["outcome"] == "not testable":
                    why.append({"question": "the hypothesis had no control arm", "evidence": ev + ".test", "note": e["test"]["why"] + ": the other group has no admitted event"})
                t["gates"] = {"evaluable": not e["closure"]["grade"].startswith("not tested"), "why": why, "guard_kept": e["gate6"]["guard_on_test_cells"]["kept"]}
                t["arms"] = {}
                for a in ("other", "rows_other"):
                    if a in e["arms"]:
                        t["arms"][a] = {**e["arms"][a]["beyond_units_at_Z"], "pair": pair(g), "evidence": ev + ".arms." + a}
                if "probe" in X and X["probe"]["group"] == g:
                    p = X["probe"]
                    t["self_events"] = [{"site": p["site"], **sub(p["grid"]["scopes"]["own"], "u", "ci"), "evidence": f"{fname}:{r}.probe", "note": "not admitted (A2); no reader of Z reads it"}]
                targets.append(t)
    world["typed_right"], world["typed_group"] = min(typed)
    world["typed_evidence"] = "round_development.json:PDL1.gate"
    return targets, world


def md(targets, items, c):
    f = lambda v: "" if v is None else (f"{v:.2f}".replace("-", "−") if isinstance(v, float) else str(v))
    ci = lambda x: "" if not x else f"[{f(x[0])}, {f(x[1])}]"
    o = ["## Targets at the end of the round", "", "| readout | group | part | status | own output after the repair | added beyond Z |", "|---|---|---|---|---|---|"]
    for t in targets:
        m = t.get("model")
        o.append(f"| {t['readout']} | {t['group']} | {t['part']} | {outer.target_status(t)} | {(f(m['own']['u']) + ' ' + ci(m['own']['ci']) + ' ' + m['own']['verdict']) if m else ''} | "
                 f"{(f(m['added']['added']) + ' ' + ci(m['added']['ci']) + ' ' + m['added']['reading']) if m else ''} |")
    o += ["", "## Open items", "", "| id | readout | group | event | reading | lever | question | value | units now → needed | after | the next step needs |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for i in items:
        need = "" if (i["resampling_units"] is None or i["reading"] not in ("unresolved", "small response")) else (
            f"{i['resampling_units']} → {'not reachable' if (i['reading'] == 'unresolved' and i['units_needed'] is None) else (f(i['units_needed']) or 'see note')}")
        val = "; ".join(f"{f(x['value'])} {ci(x['ci'])}".strip() + (f" ({x['group']})" if len(i["observations"]) > 1 else "") for x in i["observations"])
        o.append(f"| {i['id']} | {i['readout']} | {i['group'] or ''} | {i['event']} | {i['reading']} | {i['lever']} | {i['question']} | {val} | {need} | {' '.join(i['after'])} | "
                 f"{i.get('forecast') or i['next_step_needs']} |")
    o += ["", "## Counts", "", "| lever | what it is | open items |", "|---|---|---|"]
    o += [f"| {k} | {outer.LEVERS[k]} | {c['by_lever'][k]} |" for k in outer.LEVERS]
    o += ["", "| what the next step needs | open items |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c["by_what_the_next_step_needs"].items()]
    return "\n".join(o) + "\n"


def routes_md():
    o = ["| reading | what it means | lever | actions, cheapest first | what not to do |", "|---|---|---|---|---|"]
    for k, (lever, meaning, actions, dont) in outer.ROUTES.items():
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
        doc = HERE.parents[2] / "OUTER_LOOP.md"
        missing = [ln for ln in routes_md().split("\n")[2:] if ln not in doc.read_text()]
        assert not missing, missing          # the table in OUTER_LOOP.md is this table
        return
    sc = a.out / "pool_screen.json"
    screen = json.loads(sc.read_text()) if sc.exists() else None
    targets, world = targets_of_round(a.round, screen)
    items = outer.open_items(targets, world)
    c = outer.counts(items)
    findings = outer.outer_check(targets, items, roles={"final": []})
    out = {"source": "runs/teg_loop_20261120/results", "defaults": outer.DEFAULTS, "levers": outer.LEVERS,
           "routes": {k: {"lever": v[0], "meaning": v[1], "actions": v[2], "not": v[3]} for k, v in outer.ROUTES.items()},
           "pool_screen": None if screen is None else {r: sub(v, "pass_A1", "admissible", "admissible_not_in_the_round", "largest_snr_not_admissible_by_A1") for r, v in screen["readouts"].items()},
           "targets": [{**sub(t, "readout", "group", "part"), "status": outer.target_status(t), "attainable": t.get("attainable"),
                        "own_output": (t.get("model") or {}).get("own"), "added_beyond_Z": (t.get("model") or {}).get("added")} for t in targets],
           "items": items, "count_by_lever": c["by_lever"], "count_by_reading": c["by_reading"], "count_by_what_the_next_step_needs": c["by_what_the_next_step_needs"],
           "coverage": outer.coverage(targets, items), "ledger_check": findings,
           "checks": {"targets": len(targets), "items": len(items), "readings_never_met_in_this_round": [k for k in outer.ROUTES if not any(i["reading"] == k for i in items)],
                      "items_whose_interval_is_approximate": [i["id"] for i in items if i["approximate"]],
                      "note": f"an interval over fewer than {outer.DEFAULTS['min_units']} events is approximate (protocol v2); every pooled interval of this round rests on five"}}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "open_items.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    (a.out / "open_items.md").write_text(md(targets, items, c))
    print(md(targets, items, c))
    print(out["checks"], out["coverage"], findings)


if __name__ == "__main__":
    main()
