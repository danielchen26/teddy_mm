"""anm_diag: a generic ANM response-sufficiency diagnosis (arXiv 2610.14549, §2.2-§2.3, Table 2, §8.3).

A study is a *declaration* (YAML or JSON) plus an *adapter* that turns data into arrays.
The declaration names the system, the measurements, the source perturbation, the readouts and the
chain of candidate retained states. `check()` refuses designs that cannot answer the ANM question
(no controlled source, a readout taken from the model's own output, approximations or analyst
compression on the computation path, input masking, partition statistics, a readout that feeds a
candidate). `run()` then computes the forward and inverse tests, the three-way verdicts and the
root-cause reading, and writes a JSON record and a Markdown report.
"""
__version__ = "0.1.0"

from .spec import load_declaration, declaration_hash  # noqa: E402
from .checks import check, CheckError  # noqa: E402
from .runner import run  # noqa: E402

__all__ = ["load_declaration", "declaration_hash", "check", "CheckError", "run"]
