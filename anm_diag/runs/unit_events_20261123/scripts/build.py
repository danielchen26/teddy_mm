#!/usr/bin/env python3
"""The round on eight donors: admission on the selection third, the arrays of the audit, and its SPEC. No fit is made here.

  python build.py --scope-work SCOPE_WORK --pack PACK --work WORK

SCOPE_WORK: the folder runs/scope_20261116/scripts/build.py wrote (hao_donors.npz: the state, 512 components of RNA, the 112
proteins and the head's own outputs for the 40,000 cells; bm_train.npz: the head's training cells). PACK: the study's pack
(cite_arrays.npz), read for the cell names and the batch only. Writes results/selection.json and, in WORK, `donors.npz` (the
cells outside the selection third), `train.npz` and `donors.yaml`.

Admission (A1) is `anm_diag.cycle.select(control=None)`: every donor against the other donors, for every protein, on the
cells whose SHA-256 of '20261123-units|<barcode>' is divisible by 3. Those cells are written to no array.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
from anm_diag import cycle  # noqa: E402

SEED, SEED_PERM = "20261123-units", 20261123
ALPHAS = [float(a) for a in np.logspace(-1, 5, 13)]
DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}
STAGES = ["RNA512", "Z", "PERM"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope-work", type=Path, required=True)
    ap.add_argument("--pack", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    a = ap.parse_args()
    out_file = RUN / "results" / "selection.json"
    if out_file.exists():
        raise SystemExit("selection.json exists: admission was made once and is not made again")
    a.work.mkdir(parents=True, exist_ok=True)
    H, P = np.load(a.scope_work / "hao_donors.npz"), np.load(a.pack, allow_pickle=False)
    names = (a.scope_work / "proteins.txt").read_text().split()
    unit, ids, batch = H["unit"].astype(str), P["barcodes"].astype(str), P["sites"].astype(str)
    assert np.array_equal(unit, P["donors"].astype(str)) and len(set(ids)) == len(ids), "the arrays and the pack are not of the same cells"
    Y = H["Y"].astype(np.float64)
    assert Y.shape[1] == len(names) == H["frozen"].shape[1]
    res = {"rule": f"a cell is in the selection third if sha256('{SEED}|<barcode>') is divisible by 3; a donor is screened against the other donors, which count as units",
           "min_snr": cycle.MIN_SNR, "donors": sorted(set(unit)), "batch_of_donor": {u: sorted(set(batch[unit == u])) for u in sorted(set(unit))}, "proteins": {}}
    share = None
    for j, name in enumerate(names):
        sel = cycle.select(ids, unit, Y[:, j], SEED, control=None)
        share = sel["share"]
        res["proteins"][name] = {"index": j, "admitted": sel["admitted"], "table": sel["table"]}
    keep = ~share
    res["cells"] = {"all": int(len(ids)), "selection_third": int(share.sum()), "read": int(keep.sum()), "read_by_donor": {u: int((keep & (unit == u)).sum()) for u in sorted(set(unit))}}
    res["admitted_pairs"] = int(sum(len(x["admitted"]) for x in res["proteins"].values()))
    res["proteins_with_an_admitted_donor"] = int(sum(bool(x["admitted"]) for x in res["proteins"].values()))
    rna = H["RNA512"][keep].astype(np.float32)
    perm = np.random.default_rng(SEED_PERM).permutation(len(rna))
    np.savez(a.work / "donors.npz", ids=ids[keep], unit=unit[keep], batch=batch[keep], Z=H["Z"][keep].astype(np.float32), RNA512=rna, PERM=rna[perm], Y=H["Y"][keep].astype(np.float32),
             frozen=H["frozen"][keep].astype(np.float32))
    T = np.load(a.scope_work / "bm_train.npz")                  # the head's training cells, for the scope `train`; its permuted view is its own
    tr = T["RNA512"].astype(np.float32)
    np.savez(a.work / "train.npz", Z=T["Z"].astype(np.float32), RNA512=tr, PERM=tr[np.random.default_rng(SEED_PERM + 1).permutation(len(tr))], Y=T["Y"].astype(np.float32))
    spec = dict(name="donors", kind="donor", study="hao2021", readout_name="112 proteins", events="each", each_readout=True, reads="Z", input_view="RNA512", max_fit=30000,
                decision=DECISION, label_basis="design", confounded_with=["batch"], keep_draws=True, readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS, "alpha_per_target": True}},
                keys=dict(readouts="Y", stages={k: k for k in STAGES}, unit="unit", frozen={"array": "frozen"}), arrays={"audit": "donors.npz", "train": "train.npz"})
    (a.work / "donors.yaml").write_text(yaml.safe_dump(spec, sort_keys=False))
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("cells", res["cells"])
    print("admitted pairs", res["admitted_pairs"], "| proteins with an admitted donor", res["proteins_with_an_admitted_donor"], "of", len(names))
    by = {u: sum(u in x["admitted"] for x in res["proteins"].values()) for u in res["donors"]}
    print("admitted proteins by donor", by)


if __name__ == "__main__":
    main()
