"""Where the scripts in this folder find things. Nothing here is machine-specific:

  ROOT      the anm_diag project folder (holds the anm_diag package, declarations/, runs/)
  SPLIT     this run folder (one sub-folder per run, each with declaration.yaml and anm_record.json)
  RES       results/ in this run folder
  FEATURES  exported cell features, one sub-folder per data set (frangieh2021/, papalexi2021/) plus head_weights.npz;
            set the environment variable ANM_FEATURES, default: features/ in this run folder (not kept in git)
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPLIT = HERE.parent
ROOT = SPLIT.parents[1]
RUNS = ROOT / "runs"
RES = SPLIT / "results"
FEATURES = Path(os.environ.get("ANM_FEATURES", SPLIT / "features"))
DATASET_DIR = {"frangieh2021": "frangieh2021", "papalexi2021": "papalexi2021"}
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))


def record(name):
    import json
    p = SPLIT / name / "anm_record.json"
    return json.load(open(p)) if p.exists() else None
