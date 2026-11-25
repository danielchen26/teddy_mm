#!/usr/bin/env python3
"""Where the changes read so far are, by part: the accounts of run 41 and of this run in one table. results/parts.json.

  python census.py

No reading is made here. Run 41: the five earlier targets (IRF1 and the four pathway events in THP-1; the five pathway
events on MHC-I and PD-L1 in co-culture and on PD-L1 under interferon). This run: the pairs of batch 1 and of batch 2
with a change of at least four errors. For every account: the four parts, which of them are open (interval above a
step), and for a part not in the input what the recorded events say.
"""
import json

from _b2 import RUN, RUNS

PARTS = ("not_in_the_input", "lost_before_the_state", "held_not_read", "read")


def row(source, stratum, event, readout, a, own=False):
    ch = a["chain"]
    return {"source": source, "stratum": stratum, "event": event, "readout": readout, "own_gene": own, "change": ch["change"],
            "parts": {k: ch["parts"][k]["share"] for k in PARTS}, "intervals": {k: ch["parts"][k]["ci"] for k in PARTS},
            "open": [p["part"] for p in ch["open"]], "unresolved": [p["part"] for p in ch["unresolved"]],
            "across_events": (ch.get("across_events") or {}).get("reading"), "says": ch["says"],
            "carried": {k: v["shares"]["carried"]["share"] for k, v in a["stages"].items()}}


def main():
    rows = []
    t = json.loads((RUNS / "account_20261125" / "results" / "thp1.json").read_text())
    rows += [row("run 41, development", "thp1", e, "PDL1", a) for e, a in t["events"].items()]
    m = json.loads((RUNS / "account_20261125" / "results" / "melanoma.json").read_text())
    for tg, v in m["targets"].items():
        rows += [row("run 41, development", v["stratum"], e, v["readout"], a) for e, a in v["events"].items()]
    R = json.loads((RUN / "results" / "reading.json").read_text())
    for g, s in R["strata"].items():
        for key, a in s["accounts"].items():
            e, r = key.split("|")
            rows.append(row("batch 2" if a["of_the_batch"] else "batch 1", g, e, r, a, own=a["own_gene"]))
    count = lambda sel: {"accounts": len(sel), **{f"open: {k}": sum(k in r["open"] for r in sel) for k in PARTS[:3]}, "nothing_open": sum(not r["open"] for r in sel),
                         "not_in_the_input_and_not_forecast_by_other_events": sum("not_in_the_input" in r["open"] and r["across_events"] == "the other events do not forecast it" for r in sel)}
    out = {"accounts": rows, "all": count(rows), "by_source": {s: count([r for r in rows if r["source"] == s]) for s in sorted({r["source"] for r in rows})},
           "own_gene": count([r for r in rows if r["own_gene"]]), "not_own_gene": count([r for r in rows if not r["own_gene"]])}
    (RUN / "results" / "parts.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in ("all", "by_source", "own_gene", "not_own_gene")}, indent=1))
    for r in rows:
        print(f"{r['source']:20s} {r['stratum']:7s} {r['event']:9s} {r['readout']:5s} " + " ".join(f"{k.split('_')[0]} {r['parts'][k]:5.0%}" for k in PARTS) + f" | open: {', '.join(r['open']) or '-'}")


if __name__ == "__main__":
    main()
