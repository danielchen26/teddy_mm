"""Write the declarations: those of runs/loop_20261112 for the head with its MLP repaired, with the name and the note on
the head changed. Cells, stages, readouts, folds, bootstrap, margins and seed are the same; with split seeds 1 and 2 the
audited half of the control cells is another half."""
import yaml

from _paths import ARM, LOOP, RUN, audits

OLD = "Head: the head with only its latent network (MLP) retrained."
NEW = ("Head: the head with only its latent network (MLP) retrained, in one fit for both cell types and for every measured protein the head predicts, "
       "with a term in the loss that holds its outputs on bone-marrow training cells (round 3, see DESIGN.md of this folder). Split of the control cells: seed {seed}.")
for name, (kind, ds, g, seed) in audits().items():
    src = f"pooled_{ds}_mlp_only" if kind == "pooled" else f"target_{ds}_{g}_mlp_only"
    d = yaml.safe_load(open(LOOP / src / "declaration.yaml"))
    d["name"] = "repair_" + name
    assert d["system"]["notes"].count(OLD) == 1
    d["system"]["notes"] = d["system"]["notes"].replace(OLD, NEW.format(seed=seed)).replace("Closed-loop test, see DESIGN.md.", "Closed-loop test, round 3.")
    (RUN / name).mkdir(exist_ok=True)
    yaml.safe_dump(d, open(RUN / name / "declaration.yaml", "w"), allow_unicode=True, sort_keys=False)
print("wrote", len(audits()), "declarations")
