"""The heads of the round, as DESIGN.md lists them. Nothing here reads data."""
from __future__ import annotations

DONORS = [f"hao2021_P{i}" for i in range(1, 9)]
BATCH_OF = {d: ("hao2021_Batch1" if i < 4 else "hao2021_Batch2") for i, d in enumerate(DONORS)}
FIT_SETS = {"M1": [1, 2, 5, 6], "M2": [3, 4, 7, 8], "B1": [1, 2, 3, 4], "B2": [5, 6, 7, 8]}
BATCH_SEEN = {"M1": True, "M2": True, "B1": False, "B2": False}
SEEDS = (0, 1, 2)
N_ADDED, N_DOSE = 20000, 5000
ONE_DONOR = {"M1": [1, 5], "M2": [3, 7]}          # the first fit donor of each batch


def fit_donors(fit_set):
    return [f"hao2021_P{i}" for i in FIT_SETS[fit_set]]


def scored_donors(fit_set):
    return [d for d in DONORS if d not in fit_donors(fit_set)]


def heads(seeds=SEEDS, n_train=67405, final=False):
    """Development heads in the order they are fitted: seed by seed, the long `function` fit last in its seed.
    final=True: the two heads gate 5 adds per seed (all eight donors, 2,500 cells each)."""
    out = []
    for s in seeds:
        def h(name, **kw):
            out.append({"name": f"s{s}_{name}", "seed": s, "arm": name.split("_M")[0].split("_B")[0].split("_all")[0], "n_bm": None,
                        "pbmc": {}, "shuffle_labels": False, "hidden": 512, "epochs": 40, "fit_set": None, **kw})
        if final:
            every = {d: N_ADDED // len(DONORS) for d in DONORS}
            h("coverage_all", arm="coverage", n_bm=n_train - N_ADDED, pbmc=every, fit_set="all")
            h("labels_all", arm="labels", n_bm=n_train - N_ADDED, pbmc=every, shuffle_labels=True, fit_set="all")
            continue
        h("as_is", arm="as_is")
        h("amount", arm="amount", n_bm=n_train - N_ADDED)
        for fs in FIT_SETS:
            h(f"coverage_{fs}", arm="coverage", n_bm=n_train - N_ADDED, pbmc={d: N_ADDED // 4 for d in fit_donors(fs)}, fit_set=fs)
        for fs in ("M1", "M2"):
            h(f"labels_{fs}", arm="labels", n_bm=n_train - N_ADDED, pbmc={d: N_ADDED // 4 for d in fit_donors(fs)}, shuffle_labels=True,
              fit_set=fs)
        for fs in ("M1", "M2"):
            for i in ONE_DONOR[fs]:
                h(f"one_donor_{fs}_P{i}", arm="one_donor", n_bm=n_train - N_DOSE, pbmc={f"hao2021_P{i}": N_DOSE}, fit_set=fs)
        for fs in ("M1", "M2"):
            h(f"four_donors_{fs}", arm="four_donors", n_bm=n_train - N_DOSE, pbmc={d: N_DOSE // 4 for d in fit_donors(fs)}, fit_set=fs)
        h("function", arm="function", hidden=1024, epochs=80)
    return out


def head_by_name(name, n_train=67405):
    for final in (False, True):
        for c in heads(n_train=n_train, final=final):
            if c["name"] == name:
                return c
    raise SystemExit(f"no head named {name}")
