"""Write the 12 declarations: those of runs/loop_20261112 for the head with its MLP repaired, with the name and the note
on the head changed. Cells, stages, readouts, folds, bootstrap, margins and seed are the same."""
import yaml

from _paths import ARM, LOOP, RUN, audits

OLD = "Head: the head with only its latent network (MLP) retrained."
NEW = "Head: the head with only its latent network (MLP) retrained, with a term in the loss that holds its outputs on bone-marrow training cells (guarded repair, see DESIGN.md of this folder)."
for name in audits():
    d = yaml.safe_load(open(LOOP / name.replace(ARM, "mlp_only") / "declaration.yaml"))
    d["name"] = "loopguard_" + name
    assert d["system"]["notes"].count(OLD) == 1
    d["system"]["notes"] = d["system"]["notes"].replace(OLD, NEW).replace("Closed-loop test, see DESIGN.md.", "Closed-loop test with a guard.")
    (RUN / name).mkdir(exist_ok=True)
    yaml.safe_dump(d, open(RUN / name / "declaration.yaml", "w"), allow_unicode=True, sort_keys=False)
print("wrote", len(audits()), "declarations")
