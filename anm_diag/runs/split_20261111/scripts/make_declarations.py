"""Write the declarations of this folder from the sealed head-fix declarations: one per target, one per batch, and a
quick pooled rerun at Z, head latent and head output. Each is the sealed design with only the treated class, the
strata or the resampling unit changed; the seal fields are dropped and the status is exploratory."""
import copy

import yaml

from _paths import ROOT, SPLIT

SEAL = ("status", "sealed_at", "data_hashes", "code_hash", "seal", "sealed", "seal_digest")


def base(name):
    d = yaml.safe_load(open(ROOT / "declarations" / f"teddy_mm_{name}_headfix.sealed.yaml"))
    for k in SEAL:
        d.pop(k, None)
    d["status"] = "exploratory"
    return d


def write(d, name, treated, cluster, strata=None, quick=False, drop=()):
    x = copy.deepcopy(d)
    x["name"], x["source"]["treated"], x["source"]["cluster"] = name, treated, cluster
    if strata is not None:
        x["source"]["stratum"] = strata
    if quick:
        x["nodes"] = [n for n in x["nodes"] if n["id"] in ("Z", "HEAD_LAT", "HEAD_OUT")]
        for n in x["nodes"]:
            if n["id"] == "Z":
                n["from"] = "rna_counts"
        x["chain"] = ["Z", "HEAD_LAT", "HEAD_OUT"]
    else:
        x["nodes"] = [n for n in x["nodes"] if n["id"] not in drop]
    (SPLIT / name).mkdir(exist_ok=True)
    yaml.safe_dump(x, open(SPLIT / name / "declaration.yaml", "w"), allow_unicode=True, sort_keys=False)


fr = base("frangieh2021")
G = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
write(fr, "frangieh2021_pooled_full", G, "target")              # every node: the selection the features are exported for
write(fr, "frangieh2021_pooled_quick", G, "target", quick=True)
write(fr, "frangieh2021_pooled_quick_byguide", G, "guide", quick=True)
for g in G:
    write(fr, f"frangieh2021_target_{g}", [g], "guide")
pa = base("papalexi2021")
G = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]
write(pa, "papalexi2021_pooled_full", G, "target")
write(pa, "papalexi2021_pooled_quick", G, "target", quick=True)
for g in G:
    write(pa, f"papalexi2021_target_{g}", [g], "guide", drop=("SET",))
for s in ("rep1-tx", "rep3-tx", "rep4-tx"):                 # rep2-tx has 73 cells
    write(pa, f"papalexi2021_batch_{s.split('-')[0]}", G, "target", strata=[s], drop=("SET",))
