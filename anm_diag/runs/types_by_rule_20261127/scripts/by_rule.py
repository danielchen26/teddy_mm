#!/usr/bin/env python3
"""Types by rule: a law that is told of no type flags, proposes and tests its types. Nothing never read is opened.

  python by_rule.py --part twin --library LIBRARY.npz --cells batch4_cells.tsv --states B2 B3 B4 --packs PACK_A PACK_B PACK_4
  python by_rule.py --part real (the same)

Until now the one type of the system (the interferon pathway: five genes) was declared by the author. Here a type is
found as an edge is (`events.propose_types`, `account.test_types`):

  flag     the law of a library, told of no type, flags an event that moves the stage where a readout of control cells
           misreads what the event does
  propose  a flagged event, with every event that moves the stage the same way, is proposed as one type
  test     rows that did not propose it (another condition) admit it if, on the pairs it teaches, the records of its
           other events bring the forecast closer more often than chance and the error falls; else they refuse it

twin -> results/twin.json. A planted system with the real measurement (real control cells; made-up events), where the
answer is known: two types; one group of units whose cells stand to one side of the state as control cells do (the
readout reads them right) and that each do something else (to be proposed, since they are flagged and move the stage
the same way, and refused); units that each move the stage their own way (no type); own-gene units and plain units.
real -> results/real.json. The library under interferon, told of no type: what it proposes, and what the rows of
the co-culture and of the cells with no signal say; then the other way round. Every row was read before
(retrospective): the author's declaration came from the same knowledge, so this shows that the rule can replace the
hand, and it is no test on data never read.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parent / "batch4_frozen_20261127" / "scripts"))
import _b4 as B  # noqa: E402
from _b4 import V, V6  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

CELLS, SPREAD, HIDDEN = 200, {"library": 0.02, "stratum": 0.05}, 0.04
OWN = ["CD119", "PDL1", "CD29", "CD44", "CD47", "CD49f"]
# what an event of a type does to proteins, in spreads of the readout, when it moves the state as far as a real pathway knockout
EFFECT = {"A": {"MHC1": -0.8, "PDL1": -0.7, "CD119": 0.3}, "B": {"CD47": -0.6, "CD58": 0.5, "CD44": -0.4}}
UNITS = {"A": [f"A{i}" for i in range(1, 6)], "B": [f"B{i}" for i in range(1, 5)], "C": [f"C{i}" for i in range(1, 4)], "M": [f"M{i}" for i in range(1, 4)],
         "O": [f"O{i:02d}" for i in range(1, 7)], "N": [f"N{i:02d}" for i in range(1, 21)]}
ALL = [u for v in UNITS.values() for u in v]
SETTLED = 8                                                     # the worlds looked at while the rule of the test was settled


def condition(rng, Z, L, panel, plan, name, n_control):
    """The cells of one condition: control cells as they are, and every unit with what is planted for it."""
    col = {r: panel.index(ab) for r, (ab, _g) in V.MELANOMA.items()}
    iso = [panel.index(k) for k in V.ISOTYPES]
    order = rng.permutation(len(Z))
    ctl, rest = order[:n_control], order[n_control:]
    clr_c = L[ctl] - L[ctl].mean(1, keepdims=True)
    sd = {r: float(clr_c[:, c].std()) for r, c in col.items()}
    excess = {r: float((L[ctl][:, c] - L[ctl][:, iso].mean(1)).mean()) for r, c in col.items()}
    xs, Ls, un = [Z[ctl]], [L[ctl]], ["control"] * n_control
    side = rest[np.argsort(-(Z[rest] @ plan["C1"]["side"]))[:CELLS * len(UNITS["C"])]]     # the cells that stand furthest to one side, as they are
    side = dict(zip(UNITS["C"], np.array_split(rng.permutation(side), len(UNITS["C"]))))
    for u in ALL:
        idx = side[u] if u in side else rng.choice(rest, CELLS, replace=False)
        x, y = Z[idx].copy(), L[idx].copy()
        p = plan[u]
        for r in V.MELANOMA:
            d = rng.normal(0, SPREAD[name]) * sd[r]
            if p.get("own") == r:
                d = -p["share"] * excess[r]
            elif (name, r) in p["hidden"]:
                d = p["hidden"][(name, r)] * sd[r]
            d += p.get("does", {}).get(r, 0.0) * sd[r]
            y[:, col[r]] += d
        if "moves" in p:
            x = V.unit_length(x + p["moves"])
        xs.append(x); Ls.append(y); un += [u] * CELLS
    x, Lall, un = np.vstack(xs), np.vstack(Ls), np.array(un)
    clr = Lall - Lall.mean(1, keepdims=True)
    bg, C = Lall[:, iso].mean(1), un == "control"
    return x, {r: clr[:, c] for r, c in col.items()}, un, {r: V.scale(Lall[C][:, c], bg[C], clr[C][:, c]) for r, c in col.items()}


def world(seed, pools, panel, dirs, amp0):
    rng = np.random.default_rng([20261127, seed])
    plan = {}
    for u in ALL:
        k, p = u[0], {"hidden": {}}
        amp = float(rng.uniform(0.6, 1.4))                           # how far the unit moves the state, as a share of a real pathway knockout
        if k in "AB":
            p.update(moves=amp * amp0 * dirs[k], does={r: amp * v for r, v in EFFECT[k].items()})
        elif k == "C":                                              # cells that stand to one side as they are, and every unit does something else
            rs = rng.choice(list(V.MELANOMA), 3, replace=False)
            p.update(side=dirs["C"], does={str(r): float(rng.choice([-1, 1]) * rng.uniform(0.4, 0.9)) for r in rs})
        elif k == "M":                                              # a direction of its own, and nothing else
            p.update(moves=amp * amp0 * dirs[u])
        elif k == "O":
            p.update(own=OWN[int(u[1:]) - 1], share=float(rng.uniform(0.2, 0.8)))
        for name in ("library", "stratum"):
            for r in V.MELANOMA:
                if p.get("own") != r and r not in p.get("does", {}) and rng.random() < HIDDEN:
                    p["hidden"][(name, r)] = float(rng.choice([-1, 1]) * rng.uniform(0.3, 1.5))
        plan[u] = p
    table = {(u, p["own"]): "own" for u, p in plan.items() if "own" in p}
    xl, Yl, ul, scl = condition(rng, *pools["library"], panel, plan, "library", 1500)
    xs, Ys, us, scs = condition(rng, *pools["stratum"], panel, plan, "stratum", 2000)
    lib = V6.library_law(xl, Yl, ul, table, scl, types={})          # told of no type
    flags = lib["law"].flags()
    proposed = EV.propose_types(flags, lib["law"].shifts(), moving=lib["law"].moving())
    laws, known, _first = V6.stratum_laws({"all": {"x": xs, "Y": Ys, "unit": us, "table": table, "scales": scs}}, lib, types={})
    _law, ed, sc = laws["all"]
    made_up = {"type_id": "+".join(UNITS["M"]), "kind": "type", "members": UNITS["M"], "basis": "made up", "status": "proposed", "stratum": None, "evidence": None, "provenance": "made up: events that each move the stage their own way"}
    tested = AC.test_types(proposed + [made_up], xs, Ys, us, edges=ed, scale=sc, known=V6.values(known), fallback=lib["law"].sigma0)
    movers = sorted({f["event"] for f in flags if f["flag"] == "a direction the law does not cover"})
    return {"flagged_as_a_direction_the_law_does_not_cover": movers, "flagged_as_a_missing_edge": len([f for f in flags if f["flag"] == "a missing edge"]),
            "types": [{"members": t["members"], "basis": t["basis"], "status": t["status"], "cosine": (t.get("evidence") or {}).get("cosine_of_the_shifts"),
                       **{k: t["test"].get(k) for k in ("events_that_move_the_stage", "pairs", "pairs_it_teaches", "events_with_a_pair_it_teaches", "on_the_pairs_it_teaches", "on_all_its_pairs", "why")}} for t in tested]}


def say(tt):
    on, al = tt["on_the_pairs_it_teaches"], tt["on_all_its_pairs"]
    if not tt["pairs_it_teaches"]:
        return f"teaches 0 of {tt['pairs']} pairs"
    return (f"teaches {tt['pairs_it_teaches']} of {tt['pairs']} pairs of {len(tt['events_with_a_pair_it_teaches'])} events: error {on['error_without_the_type']:.3f} -> {on['error_with_the_type']:.3f}, removed {on['share_of_the_squared_error_removed']:.2f}, "
            f"z {on['z']:.1f}, closer {on['closer']} (by chance {on['by_chance']:.3f}) | all pairs: {al['error_without_the_type']:.3f} -> {al['error_with_the_type']:.3f}, removed {al['share_of_the_squared_error_removed']:.2f}")


def part_twin(a):
    packs = B.Packs(a.packs)
    panel = packs.panel[B.STUDY]
    logs = lambda ids: np.array([packs.log[B.STUDY][i] for i in ids])
    d = np.load(a.library)
    c0 = d["A"] == 0
    cells, st = B.table(a.cells), B.States(a.states)
    c = B.stratum("mel_co", cells, st)
    ctl = c["kind"] == "control"
    pools = {"library": (V.unit_length(d["Z"][c0].astype(np.float64)), logs(d["ids"].astype(str)[c0])), "stratum": (c["Z"][ctl], logs(c["ids"][ctl]))}
    # real directions of the state: how the cells of real units moved in the co-culture
    shift = {g: c["Z"][c["unit"] == g].mean(0) - c["Z"][ctl].mean(0) for g in sorted(set(c["unit"].tolist()) - {"control"}) if (c["unit"] == g).sum() >= 100}
    path = np.array([shift[g] for g in V.PATHWAY if g in shift])
    amp0 = float(np.median(np.linalg.norm(path, axis=1)))
    basis = [path.mean(0) / np.linalg.norm(path.mean(0))]
    rest = {g: v for g, v in shift.items() if g not in V.PATHWAY}
    while len(basis) < 6:                                           # the unit that moves the state furthest in a direction not taken yet
        left = {g: v - sum((v @ b) * b for b in basis) for g, v in rest.items()}
        g = max(left, key=lambda k: np.linalg.norm(left[k]))
        basis.append(left[g] / np.linalg.norm(left[g]))
        rest.pop(g)
    dirs = dict(zip(["A", "B", "C"] + UNITS["M"], basis))
    print(f"pools: library {len(pools['library'][0])} control cells, stratum {len(pools['stratum'][0])}; a pathway knockout moves the state by {amp0:.3f}", flush=True)
    res = {"planted": {"units": UNITS, "effect_of_a_type": EFFECT, "cells_per_unit": CELLS, "spread_with_no_edge": SPREAD, "the_state_moves_by": amp0, "min_cosine": EV.MIN_COSINE, "a_pair_is_taught_from": AC.TAUGHT}, "worlds": []}
    want = {k: UNITS[k] for k in "ABC"}
    for seed in range(a.worlds):
        w = world(seed, pools, panel, dirs, amp0)
        got = {tuple(t["members"]): t for t in w["types"]}
        w["as_planted"] = {k: (got[tuple(v)]["status"] if tuple(v) in got else "not proposed") for k, v in want.items()}
        w["made_up"] = got[tuple(UNITS["M"])]["status"]
        w["other_types_proposed"] = [list(m) for m, t in got.items() if list(m) not in list(want.values()) and t["basis"] != "made up"]
        res["worlds"].append(w)
        print(f"world {seed}: flagged {w['flagged_as_a_direction_the_law_does_not_cover']} | A {w['as_planted']['A']}, B {w['as_planted']['B']}, C {w['as_planted']['C']}, made up {w['made_up']} | others {w['other_types_proposed']}", flush=True)
        for t in w["types"]:
            print(f"    {''.join(sorted({m[0] for m in t['members']}))} {t['status']:8s} {say(t)}", flush=True)
    def summary(ws):
        count = lambda f: sum(1 for w in ws if f(w))
        return {"worlds": len(ws), "type_A_proposed_as_planted": count(lambda w: w["as_planted"]["A"] != "not proposed"), "type_A_admitted": count(lambda w: w["as_planted"]["A"] == "admitted"),
                "type_B_proposed_as_planted": count(lambda w: w["as_planted"]["B"] != "not proposed"), "type_B_admitted": count(lambda w: w["as_planted"]["B"] == "admitted"),
                "group_C_proposed": count(lambda w: w["as_planted"]["C"] != "not proposed"), "group_C_refused": count(lambda w: w["as_planted"]["C"] == "refused"),
                "made_up_type_refused": count(lambda w: w["made_up"] == "refused"), "worlds_with_another_type_proposed": count(lambda w: w["other_types_proposed"]),
                "worlds_with_a_wrong_type_admitted": count(lambda w: any(t["status"] == "admitted" and t["members"] not in (UNITS["A"], UNITS["B"]) for t in w["types"]))}
    # the rule of the test was settled while looking at worlds 0 to 7; the worlds after them were run once, with the rule as it is
    res["summary"] = {"worlds_the_rule_was_settled_on": summary(res["worlds"][:SETTLED]), "worlds_run_once_after": summary(res["worlds"][SETTLED:])}
    (RUN / "results" / "twin.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    for k, v in res["summary"].items():
        print(k, v)


def part_real(a):
    packs = B.Packs(a.packs)
    cells, st = B.table(a.cells), B.States(a.states)
    d = np.load(a.library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))
    S = {g: B.stratum(g, cells, st) for g in B.STRATA}
    for c in S.values():
        c["Y"], c["log"], c["bg"] = packs.readouts(B.STUDY, c["ids"])
    YL, logL, bgL = packs.readouts(B.STUDY, lib_ids)
    ZL = V.unit_length(d["Z"].astype(np.float64))
    scale = {r: {"library": V.scale(logL[r][un == "control"], bgL[un == "control"], YL[r][un == "control"])} for r in B.READOUTS}
    for g, c in S.items():
        C = c["kind"] == "control"
        for r in B.READOUTS:
            scale[r][g] = V.scale(c["log"][r][C], c["bg"][C], c["Y"][r][C])
    sc_of = lambda g: {r: scale[r][g] for r in B.READOUTS}
    sites = sorted((set(un.tolist()) | {u for c in S.values() for u in c["unit"].tolist()}) - {"control"})
    declared, tab = V.declared(sites)
    on_record, _w = B.edges_on_record()
    b4 = json.loads((RUN.parent / "batch4_frozen_20261127" / "results" / "reading.json").read_text())["edges"]["tested"]
    admitted = [e for e in on_record + b4 if e["status"] == "admitted"]
    tab_of = {g: {**tab, **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}} for g in S}
    hand = {"type_id": "declared by hand", "kind": "type", "members": sorted(V.PATHWAY), "basis": "declaration", "status": "proposed", "stratum": None, "evidence": None, "provenance": "run 40"}
    res = {"declared_by_hand": sorted(V.PATHWAY), "min_cosine": EV.MIN_COSINE, "a_pair_is_taught_from": AC.TAUGHT, "from": {}}

    def brief(t):
        return {"type": t["type_id"], "members": t["members"], "status": t["status"], "cosine": (t.get("evidence") or {}).get("cosine_of_the_shifts"), "seeded_by": (t.get("evidence") or {}).get("seeded_by"), "joined_by_direction": (t.get("evidence") or {}).get("joined_by_the_direction_of_their_shift"),
                **{k: t.get("test", {}).get(k) for k in ("events_that_move_the_stage", "pairs", "pairs_it_teaches", "events_with_a_pair_it_teaches", "on_the_pairs_it_teaches", "on_all_its_pairs", "by_event", "why")}}

    # --- the library under interferon, told of no type
    lib = V6.library_law(ZL, YL, un, tab, sc_of("library"), types={})
    flags = lib["law"].flags()
    proposed = EV.propose_types(flags, lib["law"].shifts(), moving=lib["law"].moving())
    movers = sorted({f["event"] for f in flags if f["flag"] == "a direction the law does not cover"})
    sh = lib["law"].shifts()
    cos = lambda p, q: float(sh[p] @ sh[q] / np.linalg.norm(sh[p]) / np.linalg.norm(sh[q]))
    mv = lib["law"].moving()
    out = {"events": len(lib["law"].events), "flags": len(flags), "events_flagged_as_a_direction_the_law_does_not_cover": movers, "events_that_move_the_stage": mv,
           "cosines_of_their_shifts": {f"{p}|{q}": cos(p, q) for i, p in enumerate(mv) for q in mv[i + 1:]},
           "proposed": [brief(t) for t in proposed], "the_same_events_as_declared_by_hand": any(sorted(t["members"]) == sorted(V.PATHWAY) for t in proposed), "asked_of": {}}
    print("library: flagged as a direction", movers, "| proposed", [t["members"] for t in proposed], flush=True)
    known_of = {}
    for g, c in S.items():
        laws, known, _first = V6.stratum_laws({"all": {"x": c["Z"], "Y": c["Y"], "unit": c["unit"], "table": tab_of[g], "scales": sc_of(g)}}, lib, types={})
        law, ed, sc = laws["all"]
        known_of[g] = (law, ed, sc, known)
        tested = AC.test_types(proposed + [hand], c["Z"], c["Y"], c["unit"], edges=ed, scale=sc, known=V6.values(known), fallback=lib["law"].sigma0)
        out["asked_of"][g] = [brief(t) for t in tested]
        for t in tested:
            tt = t["test"]
            print(f"  asked of {g}: {t['type_id']}: {t['status']} | " + (tt["why"] if "on_the_pairs_it_teaches" not in tt else say(tt)), flush=True)
    res["from"]["library"] = out
    # --- the other way: each stratum, told of no type; asked of the library
    for g, c in S.items():
        law, ed, sc, known = known_of[g]
        flags = law.flags()
        proposed = EV.propose_types(flags, law.shifts(), moving=law.moving())
        movers = sorted({f["event"] for f in flags if f["flag"] == "a direction the law does not cover"})
        tested = AC.test_types(proposed, ZL, YL, un, edges=lib["edges"], scale=lib["scales"], known=V6.values(lib["known"]))
        res["from"][g] = {"events": len(law.events), "flags": len(flags), "events_flagged_as_a_direction_the_law_does_not_cover": movers, "events_that_move_the_stage": law.moving(), "proposed": [brief(t) for t in proposed],
                          "the_same_events_as_declared_by_hand": any(sorted(t["members"]) == sorted(V.PATHWAY) for t in proposed), "asked_of": {"library": [brief(t) for t in tested]}}
        print(f"{g}: flagged as a direction {movers} | proposed {[t['members'] for t in proposed]} | on the library's rows: {[(t['type_id'], t['status'], say(t['test']) if 'on_the_pairs_it_teaches' in t['test'] else t['test']['why']) for t in tested]}", flush=True)
    (RUN / "results" / "real.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["twin", "real"], required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--worlds", type=int, default=24)
    a = ap.parse_args()
    {"twin": part_twin, "real": part_real}[a.part](a)


if __name__ == "__main__":
    main()
