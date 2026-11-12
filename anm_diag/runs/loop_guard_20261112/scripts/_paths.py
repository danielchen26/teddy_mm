"""Where the scripts in this folder find things. Nothing here is machine-specific:

  ROOT          the anm_diag project folder
  RUN           this run folder (one sub-folder per audit, each with declaration.yaml and anm_record.json)
  LOOP          runs/loop_20261112, the first closed-loop run: its declarations, audit controls and results are reused
  SPLIT         runs/split_20261111, whose adapter is reused
  FEATURES      exported cell features and head_weights.npz, as for runs/split_20261111 (environment variable ANM_FEATURES)
  OUT_FEATURES  where heads and the features of the guarded heads are written (ANM_GUARD_FEATURES, default: features/ here, not in git)
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
SPLIT = ROOT / "runs" / "split_20261111"
RES = RUN / "results"
FEATURES = Path(os.environ.get("ANM_FEATURES", SPLIT / "features"))
OUT_FEATURES = Path(os.environ.get("ANM_GUARD_FEATURES", RUN / "features"))
BM_FIT = os.environ.get("ANM_BM_FIT")
BM_TEST = os.environ.get("ANM_BM_TEST")
ARM = "mlp_guarded"
TARGETS = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}


def audits():
    """Every audit of the design: name -> (kind, data set, held-out target or None)."""
    out = {}
    for ds, genes in TARGETS.items():
        out[f"pooled_{ds}_{ARM}"] = ("pooled", ds, None)
        for g in genes:
            out[f"target_{ds}_{g}_{ARM}"] = ("target", ds, g)
    return out


sys.path.insert(0, str(ROOT))
sys.path.append(str(SPLIT / "scripts"))
