"""Where the scripts in this folder find things. Nothing here is machine-specific:

  ROOT          the anm_diag project folder
  RUN           this run folder (one sub-folder per audit, each with declaration.yaml and anm_record.json)
  LOOP, GUARD   runs/loop_20261112 and runs/loop_guard_20261112, the earlier rounds: declarations, audit controls, results
  SPLIT         runs/split_20261111, whose adapter is reused
  FEATURES      exported cell features and head_weights.npz, as for runs/split_20261111 (environment variable ANM_FEATURES)
  OUT_FEATURES  where heads and the features of the repaired heads are written (ANM_REPAIR_FEATURES, default: features/ here, not in git)
  BM_FIT        bone-marrow cells for the holding term and for choosing its weight (ANM_BM_FIT, not in git): arrays z_train,
                z_val, adt_val, adt_names; cells of the head's own training and validation splits
  BM_TEST       bone-marrow test cells for the guard (ANM_BM_TEST, not in git): arrays z, adt, adt_names
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
ROOT = RUN.parents[1]
LOOP = ROOT / "runs" / "loop_20261112"
GUARD = ROOT / "runs" / "loop_guard_20261112"
SPLIT = ROOT / "runs" / "split_20261111"
RES = RUN / "results"
FEATURES = Path(os.environ.get("ANM_FEATURES", SPLIT / "features"))
OUT_FEATURES = Path(os.environ.get("ANM_REPAIR_FEATURES", RUN / "features"))
BM_FIT = os.environ.get("ANM_BM_FIT")
BM_TEST = os.environ.get("ANM_BM_TEST")
ARM = "joint"
SEEDS = (0, 1, 2)                                   # splits of the control cells; 0 is the split of the earlier rounds
TARGETS = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}
GENES = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1"]      # held out one at a time, in both data sets at once
TRAINED = {"frangieh2021": {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}, "papalexi2021": {"PDL1": "CD274"}}   # readout -> head output


def audits():
    """Every audit of the design: name -> (kind, data set, held-out target or None, seed)."""
    out = {}
    for s in SEEDS:
        for ds, genes in TARGETS.items():
            out[f"pooled_{ds}_{ARM}_s{s}"] = ("pooled", ds, None, s)
            if s == 0:
                for g in genes:
                    out[f"target_{ds}_{g}_{ARM}_s0"] = ("target", ds, g, 0)
    return out


sys.path.insert(0, str(ROOT))
sys.path.append(str(SPLIT / "scripts"))
