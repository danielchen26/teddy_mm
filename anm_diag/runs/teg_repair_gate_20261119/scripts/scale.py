"""After the gate was read (post hoc): the scale on which the head's own PD-L1 output is read.

The measured protein is on a log scale. The reading of the own output introduced in runs/scope_20261116 and used in
runs/teg_training_20261119 takes the head's predicted mean, the exponential of its output, and compares it with the
measurement through one slope fitted among the control cells. The statistic is unchanged by a linear change of the
prediction, not by an exponential. This script reads the heads already made on both scales, on the same cells and
resamples. No head is fitted.

  ANM_FEATURES=... ANM_TEG_OUT=<the heads of runs/teg_training_20261119> python scale.py     (writes results/scale.json)
"""
import json

import numpy as np

import gate as G
from gate import CELLS, TARGETS, forward_np, interval, load, own_output_u, reading, weights
from _common import ARMS, OUT

ALL = ["as_trained"] + ARMS


def main():
    P0, names, j = weights()
    res = {"note": "post hoc; the heads of runs/teg_training_20261119 read on two scales; no head is fitted", "data_sets": {}}
    for ds, genes in TARGETS.items():
        d = load(ds); A, tg, y = d["A"], d["tg"], d["y"]
        use = (A == 1) | d["audit_ctl"]; Cc = d["audit_ctl"]; T = A == 1
        entry = {"measured": {"controls_mean": float(y[Cc].mean()), "controls_sd": float(y[Cc].std()), "knockout_mean": float(y[T].mean()), "min": float(y.min()), "max": float(y.max())}, "arms": {}}
        for a in ALL:
            out = {g: forward_np(P0 if a == "as_trained" else {k: v for k, v in np.load(OUT / "heads" / f"{ds}__{g}__{a}.npz").items()}, d["x"])[1][:, j] for g in genes}
            pooled = np.zeros(len(A))
            for k, g in enumerate(genes):
                own = ((A == 1) & (tg == g)) | (d["ctl_group"] == k)
                pooled[own] = out[g][own]
            row = {}
            for scale, f in (("mean", np.exp), ("log", lambda v: v)):
                p = f(pooled)
                u, slope = own_output_u(y[use], p[use], A[use], tg[use])
                iv = interval(u[1:])
                valid = bool(slope > 0 and np.sign(p[T].mean() - p[Cc].mean()) == np.sign(y[T].mean() - y[Cc].mean()))
                row[scale] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope": slope, "valid": valid,
                              "correlation_among_audit_controls": float(np.corrcoef(p[Cc], y[Cc])[0, 1]),
                              "events": {g: G.u_point(y[((A == 1) & (tg == g)) | Cc], f(out[g])[((A == 1) & (tg == g)) | Cc], A[((A == 1) & (tg == g)) | Cc])[0] for g in genes}}
            row["output_sd_among_audit_controls_log_scale"] = float(pooled[Cc].std())
            row["output_shift_log_scale"] = float(pooled[T].mean() - pooled[Cc].mean())
            entry["arms"][a] = row
        res["data_sets"][ds] = entry
        print(f"\n{CELLS[ds]}: measured protein, controls {entry['measured']['controls_mean']:.2f} ± {entry['measured']['controls_sd']:.2f}, knockouts {entry['measured']['knockout_mean']:.2f}")
        print(f"{'arm':11s} {'read as the predicted mean (as before)':40s} {'read on the log scale':40s} {'sd, shift of the log output':>28s}")
        for a in ALL:
            r = entry["arms"][a]
            s = lambda q: f"{q['u']:6.3f} [{q['ci'][0]:6.3f}, {q['ci'][1]:6.3f}] {q['reading'][:5]} {'ok' if q['valid'] else 'NO'} r {q['correlation_among_audit_controls']:5.2f}"
            print(f"{a:11s} {s(r['mean']):40s} {s(r['log']):40s} {r['output_sd_among_audit_controls_log_scale']:10.3f} {r['output_shift_log_scale']:8.3f}")
    (G.RES / "scale.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
