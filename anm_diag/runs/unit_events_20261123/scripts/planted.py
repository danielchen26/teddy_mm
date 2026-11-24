#!/usr/bin/env python3
"""Events without control cells on planted systems: is every loss put where it was planted, at three sizes and over seeds?

  python planted.py --part 0 --of 2 --out ../results/planted_0.json      (one part of the runs)
  python planted.py --merge ../results/planted_0.json ../results/planted_1.json     (writes results/planted.json)

Two planted systems (anm_diag/toy.py): a study of units in which no cell is a control cell (`run_units`), and a disease,
patients against healthy donors (`run_groups`). `sound` is the same system with nothing wrong. The rule is in DESIGN.md,
Section 5, written before these runs:

  right      the loop says what was planted
  undecided  the loop says less: no item, an interval that does not decide, a weaker status, a unit that is not admitted
  wrong      the loop says something else that is definite
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from anm_diag import toy  # noqa: E402

SIZES = {"large": (dict(n_ordinary=72, n_cell=300), dict(n_healthy=12, n_patients=8, n_cell=400)),
         "medium": (dict(n_ordinary=40, n_cell=200), dict(n_healthy=10, n_patients=6, n_cell=200)),
         "small": (dict(n_ordinary=20, n_cell=100), dict(n_healthy=8, n_patients=5, n_cell=80))}
SEEDS = list(range(1, 13))                                       # seed 0 was used while the code was written
DEFINITE = ("upstream", "unreadable", "pooled, and differs", "pooled, may be the same", "self event", "mistyped", "no rows", "confounded", "guard conflict")
WEAKER = ("the model's part closed (relative)", "undecided in the model", "undecided before the model", "not graded")
CLOSED = ("relative", "absolute")
HEAD = ["h1", "h2", "h3"]
FACTS = [f"{h}: a loss in the head, closed by the repair" for h in HEAD] + ["early: upstream", "early: lost in the first step", "enc: upstream", "enc: lost in the second step",
         "modb: unreadable", "post: unreadable", "solo: a loss in the head, nothing placed before the model", "solo: the rows of the other units do not carry", "weak: refused at selection", "no other definite item"]
SOUND, DISEASE, DISEASE_SOUND = "units, nothing wrong: nothing is open and no repair is read", "disease: closed by the repair", "disease, nothing wrong: nothing is open"


def score_units(r):
    it = r["items"]
    rd = lambda e: {i["reading"] for i in it if i["event"] == e}
    adm, reads = r["admitted"], r.get("reading") or {}
    grade = lambda s: reads[s]["closure"]["grade"]
    out = {}
    for k, h in enumerate(HEAD):                                 # a loss in the head, closed by the repair
        out[FACTS[k]] = "undecided" if h not in adm else ("wrong" if rd(h) & {"upstream", "unreadable"} else ("right" if grade(h) in CLOSED else "undecided"))
    out["solo: a loss in the head, nothing placed before the model"] = ("undecided" if "solo" not in adm else ("wrong" if rd("solo") & {"upstream", "unreadable"} else
                                                                        ("right" if reads["solo"]["A3"]["holds"] else "undecided")))
    for name, want, other in (("early", "upstream", "unreadable"), ("enc", "upstream", "unreadable"), ("modb", "unreadable", "upstream"), ("post", "unreadable", "upstream")):
        out[f"{name}: {want}"] = "right" if want in rd(name) else ("wrong" if other in rd(name) else "undecided")
    for name, fact in (("early", "early: lost in the first step"), ("enc", "enc: lost in the second step")):
        x = [i for i in it if i["event"] == name and i["reading"] == "upstream"]
        a, b = toy.LOST_UNITS[name]
        said = (x[0].get("note") or "") if x else ""
        out[fact] = "undecided" if "between" not in said else ("right" if said == f"a step is lost between {a} and {b}" else "wrong")
    out["solo: the rows of the other units do not carry"] = "right" if "pooled, and differs" in rd("solo") else "undecided"
    out["weak: refused at selection"] = "wrong" if "weak" in adm else "right"
    planted = {(k, v) for k, v in toy.TRUTH_UNITS.items()}
    extra = [(i["event"], i["reading"], i["value"], i["ci"]) for i in it if i["reading"] in DEFINITE and (i["event"], i["reading"]) not in planted]
    extra += [(s, "admitted") for s in adm if s not in toy.UNITS]                              # an ordinary unit that passes admission
    out["no other definite item"] = "wrong" if extra else "right"
    return out, {"status": r["status"], "admitted": adm, "grades": {s: grade(s) for s in reads}, "items": [(i["event"], i["reading"]) for i in it], "extra_definite": extra}


def score_sound(r, fact):
    definite = [(i["event"], i["reading"], i["value"], i["ci"]) for i in r["items"] if i["reading"] in DEFINITE]
    rp = r["target"].get("repair")
    repaired = rp == "read" if isinstance(rp, str) else "read" in (rp or {}).values()
    s = "wrong" if definite or repaired or r["status"] in ("open in the model", "open before the model") else ("right" if r["status"] == "closed (absolute)" and not r["items"] else "undecided")
    return {fact: s}, {"status": r["status"], "admitted": r["admitted"], "items": [(i["event"], i["reading"]) for i in r["items"]], "extra_definite": definite, "repair_read": bool(repaired)}


def score_disease(r):
    st = r["status"]
    definite = [(i["event"], i["reading"]) for i in r["items"] if i["reading"] in DEFINITE]
    s = "wrong" if definite or st in ("open in the model", "open before the model") else ("right" if st == "closed (absolute)" else "undecided")
    return {DISEASE: s}, {"status": st, "admitted": r["admitted"], "read_as": r.get("read_as"), "items": [(i["event"], i["reading"]) for i in r["items"]], "extra_definite": definite}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--merge", type=Path, nargs="*")
    a = ap.parse_args()
    if a.merge:
        order = ["units", "units sound", "disease", "disease sound"]
        runs = sorted((x for f in a.merge for x in json.loads(f.read_text())), key=lambda x: (list(SIZES).index(x["size"]), x["seed"], order.index(x["system"])))
        table = {}
        for x in runs:
            for fact, s in x["score"].items():
                c = table.setdefault(x["size"], {}).setdefault(fact, {"right": 0, "undecided": 0, "wrong": 0})
                c[s] += 1
        total = {k: {s: sum(c[s] for c in T.values()) for s in ("right", "undecided", "wrong")} for k, T in table.items()}
        out = HERE.parent / "results" / "planted.json"
        out.write_text(json.dumps({"sizes": {k: {"units": v[0], "disease": v[1]} for k, v in SIZES.items()}, "seeds": sorted({x["seed"] for x in runs}), "table": table, "total": total, "runs": runs}, indent=1) + "\n")
        facts = FACTS + [SOUND, DISEASE, DISEASE_SOUND]
        print("| planted fact | " + " | ".join(f"{k}: right / undecided / wrong" for k in table) + " |\n|---|" + "---|" * len(table))
        for f in facts:
            print(f"| {f} | " + " | ".join(f"{table[k][f]['right']} / {table[k][f]['undecided']} / {table[k][f]['wrong']}" for k in table) + " |")
        print(total)
        return
    jobs = [(size, seed, system) for size in SIZES for seed in SEEDS for system in ("units", "units sound", "disease", "disease sound")]
    res = []
    for k, (size, seed, system) in enumerate(jobs):
        if k % a.of != a.part:
            continue
        su, sg = SIZES[size]
        with tempfile.TemporaryDirectory() as tmp:
            if system.startswith("units"):
                r = toy.run_units(tmp, seed=seed, sound=system.endswith("sound"), **su)
                s, detail = score_sound(r, SOUND) if system.endswith("sound") else score_units(r)
            else:
                r = toy.run_groups(tmp, seed=seed, sound=system.endswith("sound"), **sg)
                s, detail = score_sound(r, DISEASE_SOUND) if system.endswith("sound") else score_disease(r)
        res.append({"size": size, "seed": seed, "system": system, "score": s, **detail})
        print(size, seed, system, s, flush=True)
        a.out.write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
