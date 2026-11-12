"""Command line.

  python -m anm_diag check DECLARATION
  python -m anm_diag run DECLARATION --adapter ks2_teddy|perturb_cite --paths PATHS.yaml --out DIR [--no-exact]
  python -m anm_diag export-tokens DECLARATION --adapter perturb_cite --paths PATHS.yaml --out tokens.npz
  python -m anm_diag inspect --paths PATHS.yaml [--show col1,col2]       metadata only, before sealing
  python -m anm_diag seal DECLARATION --paths PATHS.yaml --out SEALED.yaml
"""
import argparse
import json
import sys

from .checks import check, errors
from .spec import load_declaration


def _adapter(name, paths):
    if name == "ks2_teddy":
        from .adapters.ks2_teddy import KS2TeddyAdapter
        return KS2TeddyAdapter(paths)
    if name == "perturb_cite":
        from .adapters.perturb_cite import PerturbCiteAdapter
        return PerturbCiteAdapter(paths)
    raise SystemExit(f"unknown adapter {name}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="anm_diag")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check"); c.add_argument("declaration")
    r = sub.add_parser("run"); r.add_argument("declaration"); r.add_argument("--adapter", required=True)
    r.add_argument("--paths", required=True); r.add_argument("--out", required=True)
    r.add_argument("--no-exact", action="store_true")
    e = sub.add_parser("export-tokens"); e.add_argument("declaration"); e.add_argument("--adapter", required=True)
    e.add_argument("--paths", required=True); e.add_argument("--out", required=True)
    i = sub.add_parser("inspect"); i.add_argument("--paths", required=True); i.add_argument("--show", default="")
    s = sub.add_parser("seal"); s.add_argument("declaration"); s.add_argument("--paths", required=True)
    s.add_argument("--out", required=True)
    v = sub.add_parser("events", help="write a declaration's source as typed events and check its closure")
    v.add_argument("declaration"); v.add_argument("--record"); v.add_argument("--forward-record")
    v.add_argument("--dataset"); v.add_argument("--out")
    a = ap.parse_args(argv)
    if a.cmd == "inspect":
        from .inspect_meta import inspect
        print(inspect(a.paths, [c for c in a.show.split(",") if c]))
        return 0
    decl = load_declaration(a.declaration)
    if a.cmd == "events":
        from . import events as ev
        rec = json.load(open(a.record)) if a.record else None
        evs = ev.source_events(decl, rec)
        report = {"closure": ev.closure(evs), "problems": ev.check_source(decl)}
        if a.forward_record and a.dataset:
            fwd = json.load(open(a.forward_record))
            evs += ev.forward_events(fwd, a.dataset)
            sites = [s for s in ([decl["source"]["treated"]] if isinstance(decl["source"]["treated"], str)
                                 else decl["source"]["treated"])]
            report["admission"] = ev.admission(fwd, a.dataset, sites)
            if rec:
                report["join"] = {rid: ev.join(fwd, rec, a.dataset, rid) for rid in rec["readouts"]}
        if a.out:
            ev.write_jsonl(evs, a.out)
        report["n_events"] = len(evs)
        print(json.dumps(report, indent=1, ensure_ascii=False))
        return 0
    if a.cmd == "seal":
        probs = check(decl)
        if errors(probs):
            print(json.dumps(probs, indent=1))
            return 1
        from .seal import seal
        d = seal(decl, a.paths, a.out)
        print(f"sealed {d['name']}: {len(d['data_hashes'])} files, code {d['code_hash'][:12]}, at {d['sealed_at']} -> {a.out}")
        return 0
    if a.cmd == "check":
        probs = check(decl)
        print(json.dumps(probs, indent=1) if probs else "no problems")
        return 1 if errors(probs) else 0
    ad = _adapter(a.adapter, a.paths)
    if a.cmd == "export-tokens":
        print(ad.export_tokens(decl, a.out))
        return 0
    if a.no_exact:
        ad.exact = lambda *args, **kw: None
    from .runner import run
    rec = run(decl, ad, a.out)
    for rid, rr in rec["readouts"].items():
        print(f"{rid}: {rr['root_cause']['cause']} — {rr['root_cause'].get('detail', '')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
