"""Verdict ledger: every run's verdict per readout and stage as one table, plus the opposite pairs."""
import json, os
from _paths import RES, RUNS, SPLIT
from anm_diag import events as ev
R = str(RUNS); O = str(SPLIT)
FR = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]; PA = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]; BA = ["rep1", "rep3", "rep4"]
runs = [("frangieh2021_headfix", f"{R}/frangieh2021_headfix/anm_record.json", "pool", "mel", "cf"), ("papalexi2021_headfix", f"{R}/papalexi2021_headfix/anm_record.json", "pool", "thp", "cf")]
runs += [(f"frangieh2021_target_{g}", f"{O}/frangieh2021_target_{g}/anm_record.json", g, "mel", "ex") for g in FR]
runs += [(f"papalexi2021_target_{g}", f"{O}/papalexi2021_target_{g}/anm_record.json", g, "thp", "ex") for g in PA]
runs += [(f"papalexi2021_batch_{b}", f"{O}/papalexi2021_batch_{b}/anm_record.json", "pool", f"thp:{b}", "ex") for b in BA]
evs, rows, used, missing = [], [], [], []
for run, path, site, ctx, grade in runs:
    if not os.path.exists(path): missing.append(run); continue
    rec = json.load(open(path)); used.append(run)
    for e in ev.verdict_events(rec, run, site=site, context=ctx):
        evs.append(e); rows.append(dict(run=run, grade=grade, ctx=ctx, src=site, ro=e["readout"], st=e["site"], v=e["verdict"], u=round(e["amplitude"], 3), ci=[round(x, 3) for x in e["ci"]], cls=e["channel"]))
conf = ev.ledger_conflicts(evs)
json.dump(dict(rows=rows, conflicts=conf, runs=used, missing=missing), open(RES / "ledger.json", "w"), ensure_ascii=False)
print("runs", len(used), "rows", len(rows), "conflicts", len(conf), "missing", missing)
for c in conf: print("  ", c["readout"], c["stage"], c["runs"], c["verdicts"], [round(u, 3) for u in c["u"]], c["differ_in"])
