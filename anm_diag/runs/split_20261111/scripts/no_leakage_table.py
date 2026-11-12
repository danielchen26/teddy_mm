"""No-leakage entries of every declared source: the two TEDDY data sets, the two KS-2 declarations, the forward run."""
import collections
import json

from _paths import RES, ROOT
from anm_diag import events as ev
from anm_diag.spec import load_declaration

out = {"sufficiency": {}, "forward": {}}
for f in ("teddy_mm_frangieh2021_headfix.sealed.yaml", "teddy_mm_papalexi2021_headfix.sealed.yaml", "ks2_teddy_ribo.yaml", "ks2_teddy_mitoribo.yaml"):
    d = load_declaration(str(ROOT / "declarations" / f))
    evs = [e for e in ev.source_events(d) if e["role"] == "treated"]
    out["sufficiency"][f] = dict(
        n_events=len(evs), status=dict(collections.Counter(e["no_leakage"]["status"] for e in evs)),
        routes={r: sorted({e["no_leakage"]["routes"][r]["status"] for e in evs}) for r in ev.LEAK_ROUTES},
        reasons={r: sorted({e["no_leakage"]["routes"][r]["reason"] for e in evs}) for r in ev.LEAK_ROUTES},
        s4=[p for p in ev.check_source(d) if p["code"] == "S4"])
fwd = json.load(open(ROOT / "runs" / "forward_typed_events" / "forward_record.json"))
for ds in fwd["datasets"]:
    evs = ev.forward_events(fwd, ds)
    out["forward"][ds] = dict(n_events=len(evs), status=dict(collections.Counter(e["no_leakage"]["status"] for e in evs)),
                              flagged_sites=sorted({e["site"] for e in evs if e["no_leakage"]["status"] == "flag"}))
json.dump(out, open(RES / "no_leakage.json", "w"), indent=1)
print({k: v["status"] for k, v in out["sufficiency"].items()}, {k: v["status"] for k, v in out["forward"].items()})
