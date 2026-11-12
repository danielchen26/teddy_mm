"""Where the scripts in this folder find things. Nothing here is machine-specific:

  ROOT          the anm_diag project folder
  RUN           this run folder (one sub-folder per run, each with declaration.yaml and anm_record.json)
  SPLIT         runs/split_20261111, whose declarations and adapter are reused
  FEATURES      exported cell features and head_weights.npz, as for runs/split_20261111 (environment variable ANM_FEATURES)
  OUT_FEATURES  where the features of the three heads on the audit half are written (ANM_EXPOSURE_FEATURES,
                default: features/ in this run folder, not kept in git)
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
OUT_FEATURES = Path(os.environ.get("ANM_EXPOSURE_FEATURES", RUN / "features"))
ARMS = ("as_trained", "without_IRF1", "with_IRF1")
TARGETS = ("IRF1", "JAK2")
sys.path.insert(0, str(ROOT))
sys.path.append(str(SPLIT / "scripts"))
