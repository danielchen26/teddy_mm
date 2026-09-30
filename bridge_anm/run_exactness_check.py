#!/usr/bin/env python3
"""ANM vs the fixed TEDDY-alone rule, cell by cell, on a bridge export (site4 cells).

For each question (O0 soft, O1 strict, O2 B/T-priority): how many cells get exactly the same call
(or the same no-call) from ANM and from the rule, no-call counts, accuracy of calls made and strict
accuracy at matched coverage (both arms abstain on the same number of cells, lowest top score
first), and the annotation-graded view (in-scope B/T/myeloid accuracy, share of out-of-scope cells
that still receive a call). Same computation as the v2 run's
outputs/anm_cite_bridge_v2/hard_proof/exactness_matched_coverage.py, with --workers.

usage: run_exactness_check.py BRIDGE_DIR OUT_JSON [--schema PATH] [--workers N]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_hard_proof as H  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bridge", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--schema", type=Path,
                    default=Path(__file__).resolve().parent / "schemas/cite_lineage_finite_field_v0.json")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    schema = json.loads(a.schema.read_text())
    bridge = a.bridge
    man = json.loads((bridge / "export_manifest.json").read_text())
    cells = H.load_cells(bridge / "cite_cells_meta.jsonl", man.get("z_rna_export"), site4_only=True)
    ids = sorted(cells)
    ev = bridge / "cite_typed_events.jsonl"
    p95 = H.p95_from_events(ev, ids)
    by = H.events_by_cell(ev, ids)
    res = {"bridge_dir": str(bridge), "schema": str(a.schema)}
    with H.Runner(a.workers, schema) as R:
        for c in H.CRIT_IDS:
            crit = H.CRITERIA[c]
            t_pred, t_top, lab = [], [], []
            for cid in ids:
                sc = H.lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
                t_pred.append(H.decide_argmax(sc, crit["readout_threshold"]))
                t_top.append(max(sc.values()))
                lab.append(H.expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit))
            outs = R.map(H._t_anm, ((c, H.make_instance(cid, cells[cid], by[cid], crit, p95)) for cid in ids))
            a_pred = [o["recommended_action"] for o in outs]
            a_top = [max(o["action_scores"].values()) if o["action_scores"] else -1e9 for o in outs]
            n = len(ids)
            same = sum(1 for x, y in zip(t_pred, a_pred) if x == y)
            both_called_diff = sum(1 for x, y in zip(t_pred, a_pred) if x is not None and y is not None and x != y)
            only_t = sum(1 for x, y in zip(t_pred, a_pred) if x is not None and y is None)
            only_a = sum(1 for x, y in zip(t_pred, a_pred) if x is None and y is not None)
            kt, ka = t_pred.count(None), a_pred.count(None)
            k = max(kt, ka)

            def at_k(pred, top, k):
                order = sorted(range(n), key=lambda i: top[i])
                drop = set(order[:k])
                p = [None if i in drop else pred[i] for i in range(n)]
                d = [i for i in range(n) if p[i] is not None and lab[i] is not None]
                cor = sum(1 for i in d if p[i] == lab[i])
                nl = sum(1 for y in lab if y is not None)
                return {"abstain": sum(1 for x in p if x is None), "acc_called": cor / max(len(d), 1), "acc_strict": cor / nl}

            res[c] = {"n": n, "exact_same_call": same, "exact_rate": same / n,
                      "both_called_disagree": both_called_diff, "rule_only_called": only_t, "anm_only_called": only_a,
                      "abstain_rule": kt, "abstain_anm": ka, "matched_k": k,
                      "rule_at_k": at_k(t_pred, t_top, k), "anm_at_k": at_k(a_pred, a_top, k)}
            ct = [cells[i].get("cell_type") for i in ids]
            for nm, pp in (("rule", t_pred), ("anm", a_pred)):
                g = H.annotation_graded(pp, ct)
                res[c]["ann_" + nm] = {"in_n": g["in_scope"]["n"], "in_called": g["in_scope"]["n_called"],
                                       "in_acc_strict": g["in_scope"]["accuracy_strict"],
                                       "in_acc_called": g["in_scope"]["accuracy_among_called"],
                                       "out_n": g["out_of_scope"]["n"], "out_called": g["out_of_scope"]["n_called"],
                                       "out_share": g["out_of_scope"]["call_share"]}
            print(c, json.dumps(res[c]), flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
