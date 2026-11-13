"""Write one declaration per audit of the design, from the declarations of runs/split_20261111: the quick pooled declaration
(all readouts; Z, head latent, head output; resampling by target) for the pooled and transfer audits, and the per-target
declaration reduced to the same three stages and PD-L1 for the single-head audits."""
import copy

import yaml

from _paths import OTHER, RUN, SPLIT, audits

HEAD = {"as_trained": "the head as trained", "mlp_only": "the head with only its latent network (MLP) retrained", "decoder_only": "the head with only its decoder retrained"}
for name, (kind, ds, arm, g) in audits().items():
    d = yaml.safe_load(open(SPLIT / (f"{ds}_target_{g}" if kind == "target" else f"{ds}_pooled_quick") / "declaration.yaml"))
    d = copy.deepcopy(d); d["name"] = "loop_" + name
    if kind == "target":
        d["nodes"] = [n for n in d["nodes"] if n["id"] in ("Z", "HEAD_LAT", "HEAD_OUT")]
        for n in d["nodes"]:
            if n["id"] == "Z":
                n["from"] = "rna_counts"
        d["chain"] = ["Z", "HEAD_LAT", "HEAD_OUT"]
        d["readouts"] = [r for r in d["readouts"] if r["id"] == "PDL1"]
    cells = {"pooled": "all knockout cells and the audit half of the control cells; each cell is passed through the head fitted with its target held out "
                       "(control cells: one of the five heads, assigned at random)",
             "target": f"the knockout cells of {g}, held out when this head was fitted, and the audit half of the control cells",
             "transfer": f"all cells; the head was fitted on {OTHER[ds]} only"}[kind]
    d["system"]["notes"] = f"Closed-loop test, see DESIGN.md. Head: {HEAD[arm]}. Cells: {cells}. No head was fitted on any audited cell. " + d["system"]["notes"]
    (RUN / name).mkdir(exist_ok=True)
    yaml.safe_dump(d, open(RUN / name / "declaration.yaml", "w"), allow_unicode=True, sort_keys=False)
print("wrote", len(audits()), "declarations")
