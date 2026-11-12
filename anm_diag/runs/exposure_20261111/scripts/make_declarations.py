"""Write the six declarations: the per-target declarations of runs/split_20261111 for IRF1 and JAK2, restricted to the
audit half, at Z, head latent and head output, readout PD-L1, once per head."""
import copy

import yaml

from _paths import ARMS, RUN, SPLIT, TARGETS

HEAD = {"as_trained": "the head as trained",
        "without_IRF1": "the head fine-tuned on the training half without IRF1 knockout cells",
        "with_IRF1": "the head fine-tuned on the training half with IRF1 knockout cells"}
for g in TARGETS:
    base = yaml.safe_load(open(SPLIT / f"papalexi2021_target_{g}" / "declaration.yaml"))
    for arm in ARMS:
        d = copy.deepcopy(base)
        d["name"] = f"exposure_{arm}_{g}"
        d["nodes"] = [n for n in d["nodes"] if n["id"] in ("Z", "HEAD_LAT", "HEAD_OUT")]
        for n in d["nodes"]:
            if n["id"] == "Z":
                n["from"] = "rna_counts"
        d["chain"] = ["Z", "HEAD_LAT", "HEAD_OUT"]
        d["readouts"] = [r for r in d["readouts"] if r["id"] == "PDL1"]
        d["system"]["notes"] = (f"Exposure test, see DESIGN.md. Head: {HEAD[arm]}. Cells: the audit half of split seed 0 "
                                "(results/audit_half_seed0.json); no head was fitted on any of them. " + d["system"]["notes"])
        (RUN / f"{arm}_{g}").mkdir(exist_ok=True)
        yaml.safe_dump(d, open(RUN / f"{arm}_{g}" / "declaration.yaml", "w"), allow_unicode=True, sort_keys=False)
        print("wrote", f"{arm}_{g}")
