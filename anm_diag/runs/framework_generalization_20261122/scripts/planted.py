#!/usr/bin/env python3
"""The general cycle on planted systems: is every loss put where it was planted, at three sizes and over seeds?

  python planted.py --seeds 12 --part 0 --of 2 --out ../results/planted_0.json      (one part of the runs)
  python planted.py --merge ../results/planted_0.json ../results/planted_1.json     (writes results/planted.json)

A planted system (anm_diag/toy.py) has two modalities and two types of cell; `sound` is the same system with nothing
wrong. For every planted fact a run is scored, by the rule written here before any run:

  right      the loop says what was planted
  undecided  the loop says less: no item, an interval that does not decide, a reading it refuses to use, a weaker status
  wrong      the loop says something else that is definite

The property asked for is that `wrong` stays at zero as the data shrink: with fewer cells the loop may say less, and it
may not say something false.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from anm_diag import toy  # noqa: E402

SIZES = {"large": dict(n_control=4000, n_event=800), "medium": dict(n_control=1500, n_event=300), "small": dict(n_control=500, n_event=100)}
DEFINITE = ("upstream", "unreadable", "pooled, and differs", "pooled, may be the same", "self event", "mistyped", "no rows", "confounded", "guard conflict")
WEAKER = ("the model's part closed (relative)", "undecided in the model", "undecided before the model", "not graded")
FACTS = ["head of type one: closed by its repair", "enc: upstream", "modb: unreadable, and the other input reads it", "post: unreadable, and no input reads it", "weak: refused at selection",
         "the repair of one type does not carry to the other", "odd: differs inside its pool", "odd: upstream", "no other definite item",
         "early: upstream", "early: lost in the first step", "enc: lost in the second step", "odd: lost in the second step"]


def score(r):
    it = r["items"]
    of = lambda g, e: [i for i in it if i["group"] == g and i["event"] == e]
    rd = lambda g, e: {i["reading"] for i in of(g, e)}
    out = {}
    st = r["types"]["one"]["status"]
    out[FACTS[0]] = "right" if st == "closed (absolute)" else ("undecided" if st in WEAKER else "wrong")
    e = rd("one", "enc")
    out[FACTS[1]] = "right" if "upstream" in e else ("wrong" if "unreadable" in e else "undecided")
    for k, name, wants_other in ((2, "modb", True), (3, "post", False)):
        x = [i for i in of("one", name) if i["reading"] == "unreadable"]
        if "upstream" in rd("one", name):
            out[FACTS[k]] = "wrong"
        elif not x:
            out[FACTS[k]] = "undecided"
        else:
            named = "view_b" in (x[0].get("forecast") or "")
            out[FACTS[k]] = "right" if named == wants_other else ("undecided" if wants_other else "wrong")
    out[FACTS[4]] = "wrong" if "weak" in r["admitted"]["one"] else "right"
    across = [i for i in it if i["question"].startswith("the whole repair across")]
    merge = [i for i in it if i["reading"] == "pooled, may be the same"]
    out[FACTS[5]] = "wrong" if merge else ("right" if across else "undecided")
    o = rd("two", "odd")
    out[FACTS[6]] = "right" if "pooled, and differs" in o else "undecided"
    out[FACTS[7]] = "right" if "upstream" in o else ("wrong" if "unreadable" in o else "undecided")
    e = rd("one", "early")
    out[FACTS[9]] = "right" if "upstream" in e else ("wrong" if "unreadable" in e else "undecided")
    for k, (g, name) in ((10, ("one", "early")), (11, ("one", "enc")), (12, ("two", "odd"))):          # between which two stages
        x = [i for i in of(g, name) if i["reading"] == "upstream"]
        a, b = toy.LOST[g, name]
        said = (x[0].get("note") or "") if x else ""
        out[FACTS[k]] = "undecided" if "between" not in said else ("right" if said == f"a step is lost between {a} and {b}" else "wrong")
    planted = {("one", "early", "upstream"), ("one", "enc", "upstream"), ("one", "modb", "unreadable"), ("one", "post", "unreadable"), ("two", "odd", "upstream"), ("two", "odd", "pooled, and differs"),
               ("one", "pool", "pooled, and differs"), ("two", "pool", "pooled, and differs")}
    extra = [(i["group"], i["event"], i["reading"], i["value"], i["ci"]) for i in it if i["reading"] in DEFINITE and (i["group"], i["event"], i["reading"]) not in planted
             and not (i["event"] == "pool" and i["question"].startswith("the whole repair across"))]
    out[FACTS[8]] = "wrong" if extra else "right"
    return out, {"statuses": {k: T["status"] for k, T in r["types"].items()}, "items": [(i["group"], i["event"], i["reading"]) for i in it], "extra_definite": extra,
                 "not_admitted": {k: sorted(set(toy.events_of(k)) - set(v)) for k, v in r["admitted"].items()}}


def score_sound(r):
    definite = [(i["group"], i["event"], i["reading"], i["value"], i["ci"]) for i in r["items"] if i["reading"] in DEFINITE]
    st = [T["status"] for T in r["types"].values()]
    rp = lambda x: x == "read" if isinstance(x, str) else "read" in (x or {}).values()
    repaired = [k for k, T in r["types"].items() if rp(T["target"].get("repair"))]
    s = "wrong" if definite or repaired or any(x in ("open in the model", "open before the model") for x in st) else ("right" if all(x == "closed (absolute)" for x in st) and not r["items"] else "undecided")
    return {"a system with nothing wrong: nothing is open": s}, {"statuses": st, "items": [(i["group"], i["event"], i["reading"]) for i in r["items"]], "extra_definite": definite, "repair_read_for": repaired}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--merge", type=Path, nargs="*")
    a = ap.parse_args()
    if a.merge:
        runs = sorted((x for f in a.merge for x in json.loads(f.read_text())), key=lambda x: (list(SIZES).index(x["size"]), x["seed"], x["sound"]))
        table = {}
        for x in runs:
            for fact, s in x["score"].items():
                c = table.setdefault(x["size"], {}).setdefault(fact, {"right": 0, "undecided": 0, "wrong": 0})
                c[s] += 1
        total = {s: sum(c[s] for T in table.values() for c in T.values()) for s in ("right", "undecided", "wrong")}
        out = HERE.parent / "results" / "planted.json"
        out.write_text(json.dumps({"sizes": SIZES, "seeds": sorted({x["seed"] for x in runs}), "table": table, "total": total, "runs": runs}, indent=1) + "\n")
        facts = FACTS + ["a system with nothing wrong: nothing is open"]
        print("| planted fact | " + " | ".join(f"{k}: right / undecided / wrong" for k in table) + " |\n|---|" + "---|" * len(table))
        for f in facts:
            print(f"| {f} | " + " | ".join(f"{table[k][f]['right']} / {table[k][f]['undecided']} / {table[k][f]['wrong']}" for k in table) + " |")
        print(total)
        return
    jobs = [(size, seed, sound) for size in SIZES for seed in range(a.seeds) for sound in (False, True)]
    res = []
    for k, (size, seed, sound) in enumerate(jobs):
        if k % a.of != a.part:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            r = toy.run(tmp, seed=seed, sound=sound, **SIZES[size])
        s, detail = (score_sound if sound else score)(r)
        res.append({"size": size, "seed": seed, "sound": sound, "score": s, **detail})
        print(size, seed, "sound" if sound else "planted", s, flush=True)
        a.out.write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
