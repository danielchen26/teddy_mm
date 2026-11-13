"""Where the scripts in this folder find things. Nothing here is machine-specific:

  ROOT          the anm_diag project folder
  RUN           this run folder (one sub-folder per audit, each with declaration.yaml and anm_record.json)
  SPLIT         runs/split_20261111, whose declarations and adapter are reused
  FEATURES      exported cell features and head_weights.npz, as for runs/split_20261111 (environment variable ANM_FEATURES)
  OUT_FEATURES  where the features of the retrained heads are written (ANM_LOOP_FEATURES, default: features/ here, not in git)
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
ROOT = RUN.parents[1]
SPLIT = ROOT / "runs" / "split_20261111"
RES = RUN / "results"
FEATURES = Path(os.environ.get("ANM_FEATURES", SPLIT / "features"))
OUT_FEATURES = Path(os.environ.get("ANM_LOOP_FEATURES", RUN / "features"))
ARMS = ("as_trained", "mlp_only", "decoder_only")
TARGETS = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}
OTHER = {"frangieh2021": "papalexi2021", "papalexi2021": "frangieh2021"}


def audits():
    """Every audit of the design: name -> (kind, data set, arm, held-out target or None)."""
    out = {}
    for ds, genes in TARGETS.items():
        for arm in ARMS:
            out[f"pooled_{ds}_{arm}"] = ("pooled", ds, arm, None)
            for g in genes:
                out[f"target_{ds}_{g}_{arm}"] = ("target", ds, arm, g)
        for arm in ARMS[1:]:
            out[f"transfer_{ds}_{arm}"] = ("transfer", ds, arm, None)
    return out


sys.path.insert(0, str(ROOT))
sys.path.append(str(SPLIT / "scripts"))
