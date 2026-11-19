"""Everything of this run, in order (DESIGN.md): the check, the fits of every arm in both data sets, the reading.

  ANM_FEATURES=... ANM_BM_FIT=... ANM_BM_TEST=... python run_all.py

About a minute per fit on a laptop CPU, 90 fits; no GPU is needed. The reading is made once, at the end.
"""
import subprocess
import sys

from _common import ARMS, TARGETS

py = [sys.executable]
subprocess.run(py + ["fit.py", "--check"], check=True)
for ds in TARGETS:
    subprocess.run(py + ["fit.py", ds] + ARMS, check=True)
subprocess.run(py + ["read.py"], check=True)
