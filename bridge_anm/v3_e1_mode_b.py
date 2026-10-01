#!/usr/bin/env python3
"""E1, Mode B rerun on the v3 key: Experiments 1 (change the question), 3 (why this call), 4 (which calls to
trust), 5 (label cost) and C2 (nested readouts), exactly as registered.

Registration: registration/registration_v3.json (experiments.E1), loaded with
v3_amend.load_registration_amended() (amendment A1 applied), and the E1 addendum
registration/addenda/E1.json (bridge_anm/v3_e1_addendum.py: sensitivity-panel bars and classifier C,
Q3 operating-point bars, declared implementation choices; train/val only).

Arms (same evidence for every method: the phase-1 head's predicted panel / training q95, clipped):
  rule        TEDDY + fixed rule (re-coded question): equal-weight panel mean, argmax, val bar
  anm         TEDDY + ANM: finite_graph_scalar engine (ANM_ROOT), one support event per panel protein at t = 0;
              by construction its calls equal the rule's (E1.1a checks it cell by cell)
  closure     ANM closure readout (readout_coordinates over each action's event sites), C2 only
  classifier  logistic regression on the same evidence, trained on training-site labels only
  margin      S_top1 - S_top2 (a confidence score)

Endpoints (registration experiments.E1 + A1): E1.1a-d, E1.3a-c, E1.4a-b, E1.5, E1.C2; key validity is reported
before any method result. Matched coverage with the A1.1 tie-break; two-stage donor -> cell bootstrap (B = 2000,
seeds.bootstrap) with every endpoint recomputed per replicate; registered win / loss / equivalent rules.

Stages (each resumable; caches carry a fingerprint of registration, amendment, addendum, code and inputs):
  prepare  evidence, keys, key validity (written first), classifiers and label-cost classifiers (train only)
  anm      ANM engine per evaluated cell (Q1 / Q2 / Q3, closure, leave-one-out of the called class)
  boot     bootstrap replicates in chunks of --chunk replicates (pooled primary, per donor, secondary)
  report   E1_results.json + REPORT.md
  all      prepare, anm, boot, report

A site4 run is refused unless registration_v3.json, amendment_A1.json and addenda/E1.json are committed with
matching hashes and the addendum's computed core reproduces from train/val. --smoke evaluates val donor 18303
cells only (two label-free pseudo-donors), with small B; it never reads a site4 row's protein, label or embedding.

Full run (repo root):
  PY=<venv>/bin/python
  ANM_ROOT=<anm_v2_fix> $PY bridge_anm/v3_e1_mode_b.py --stage all \
      --out-dir /Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E1 --threads 4
"""
from __future__ import annotations

import os
import sys


def _early_threads(argv) -> str:
    for i, a in enumerate(argv):
        if a == "--threads" and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith("--threads="):
            return a.split("=", 1)[1]
    return "4"


_T = _early_threads(sys.argv)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = _T

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e1 as e1  # noqa: E402
from lib import v3_key as vk  # noqa: E402

VERSION = "v3_e1_mode_b 1.0"
LIN = np.asarray(vk.LINEAGES)
K1 = ("primary", "annotation_only", "primary_per_batch", "primary_no_gdT158")   # Q1 / Q2 keys
K3 = ("q3", "q3_annotation_only")                                               # Q3 keys
KEY_NAMES = K1 + K3 + ("annotation", "gated", "gated_per_batch")
LC_N = [0, 25, 50, 100, 200, 500, 1000, 2000, 5000, "all"]
SCHEMA = ROOT / "bridge_anm" / "schemas" / "cite_lineage_finite_field_v0.json"
DEFAULT_OUT = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E1")
SMOKE_OUT = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E1_smoke")
SMOKE = {"n_primary": 1600, "n_secondary": 600, "seed": 0, "n_boot": 40, "n_perm": 40,
         "lc_n": [0, 25, 200, "all"], "lc_draws": 2, "mlp_seeds": 1}

_LOG: Path | None = None


def say(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG is not None:
        with _LOG.open("a") as f:
            f.write(line + "\n")


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def write_json(path: Path, obj) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(jsonable(obj), indent=1) + "\n")
    tmp.replace(path)


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez(tmp, **arrays)
    tmp.replace(path)


def ratio(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(b > 0, a / np.where(b > 0, b, 1.0), np.nan)
    return float(r) if r.ndim == 0 else r


# ============================================================================ args, guard, fingerprint
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", default="all", choices=["all", "prepare", "anm", "boot", "report"])
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration", type=Path, default=ROOT / "registration/registration_v3.json")
    p.add_argument("--amendment", type=Path, default=ROOT / "registration/amendment_A1.json")
    p.add_argument("--addendum", type=Path, default=ROOT / "registration/addenda/E1.json")
    p.add_argument("--anm-root", type=Path, default=None, help="default: $ANM_ROOT")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--smoke", action="store_true", help="val cells only, small B (workflow smoke test)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--chunk", type=int, default=100, help="bootstrap replicates per checkpoint file")
    a = p.parse_args(argv)
    a.out_dir = a.out_dir or (SMOKE_OUT if a.smoke else DEFAULT_OUT)
    return a


def _git(*args) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    except OSError:
        return None


def _committed(path: Path) -> bool:
    r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
    r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
    return r1.returncode == 0 and r2.returncode == 0


def registration_state(a) -> dict:
    """Hashes and commit state of everything E1 reads; refuses a site4 run unless all is frozen."""
    reg_sha = vk.sha256_file(a.registration)
    am_sha = vk.sha256_file(a.amendment)
    st = {"registration_sha256": reg_sha, "amendment_A1_sha256": am_sha,
          "registration_hash_file_matches": a.registration.with_name(a.registration.name + ".sha256").read_text().split()[0] == reg_sha,
          "amendment_hash_file_matches": a.amendment.with_name(a.amendment.name + ".sha256").read_text().split()[0] == am_sha,
          "registration_committed": _committed(a.registration), "amendment_committed": _committed(a.amendment)}
    if a.addendum.exists():
        add_sha = vk.sha256_file(a.addendum)
        hashes = a.addendum.parent / "HASHES.txt"
        doc = json.loads(a.addendum.read_text())
        st.update({"addendum_sha256": add_sha,
                   "addendum_hash_recorded": hashes.exists() and f"{add_sha}  E1.json" in hashes.read_text().splitlines(),
                   "addendum_committed": _committed(a.addendum),
                   "addendum_amends_these_files": doc.get("registration_sha256") == reg_sha and doc.get("amendment_A1_sha256") == am_sha,
                   "addendum_core_sha256": doc["provenance"]["core_sha256"]})
    else:
        st.update({"addendum_sha256": None, "addendum_hash_recorded": False, "addendum_committed": False,
                   "addendum_amends_these_files": False})
    code = [Path(__file__).resolve(), ROOT / "bridge_anm/lib/v3_e1.py", ROOT / "bridge_anm/v3_e1_addendum.py"]
    st["code_sha256"] = {str(c.relative_to(ROOT)): vk.sha256_file(c) for c in code}
    st["code_committed_clean"] = all(_committed(c) for c in code)
    st["git_head"] = _git("rev-parse", "HEAD")
    return st


def guard_site4(a, st: dict) -> None:
    need = {"registration hash file matches": st["registration_hash_file_matches"],
            "registration committed": st["registration_committed"],
            "amendment A1 hash file matches": st["amendment_hash_file_matches"],
            "amendment A1 committed": st["amendment_committed"],
            "addendum E1.json present and recorded in addenda/HASHES.txt": st["addendum_hash_recorded"],
            "addendum committed": st["addendum_committed"],
            "addendum amends the registration and A1 on disk": st["addendum_amends_these_files"]}
    bad = [k for k, ok in need.items() if not ok]
    if bad:
        raise SystemExit("site4 run refused (register and commit before any site4 evaluation): " + "; ".join(bad))


def check_addendum_reproduces(a, st: dict) -> dict:
    """Rebuild the addendum's computed core (train/val only) and compare with the committed one."""
    import v3_e1_addendum as ad

    ns = argparse.Namespace(processed=a.processed, z=a.z, ckpt=a.ckpt, registration=a.registration,
                            amendment=a.amendment)
    core = ad.build_core(ns)
    sha = ad.sha256_bytes(ad.content_bytes(core))
    ok = sha == st.get("addendum_core_sha256")
    if not ok and not a.smoke:
        raise SystemExit(f"addendum core does not reproduce: {sha} != {st.get('addendum_core_sha256')}")
    return {"recomputed_core_sha256": sha, "matches": ok}


def file_fingerprint(a, st: dict, cfg: dict) -> str:
    h = hashlib.sha256()
    for k in ("registration_sha256", "amendment_A1_sha256", "addendum_sha256"):
        h.update(str(st.get(k)).encode())
    for k, v in sorted(st["code_sha256"].items()):
        h.update(f"{k}={v}".encode())
    for p in (a.processed / "cite_arrays.npz", a.z, a.ckpt):
        s = p.stat()
        h.update(f"{p}:{s.st_size}:{int(s.st_mtime)}".encode())
    h.update(json.dumps(cfg, sort_keys=True).encode())
    return h.hexdigest()[:20]


# ============================================================================ config
def run_config(a, reg) -> dict:
    seeds = reg["seeds"]
    if a.smoke:
        return {"smoke": True, "n_boot": SMOKE["n_boot"], "n_perm": SMOKE["n_perm"], "lc_n": SMOKE["lc_n"],
                "lc_draws": list(seeds["label_cost_draws"])[: SMOKE["lc_draws"]],
                "mlp_seeds": list(seeds["classifier"])[: SMOKE["mlp_seeds"]], "boot_seed": int(seeds["bootstrap"]),
                "perm_seed": int(seeds["global"]), "smoke_cells": SMOKE}
    return {"smoke": False, "n_boot": 2000, "n_perm": 1000, "lc_n": LC_N, "lc_draws": list(seeds["label_cost_draws"]),
            "mlp_seeds": list(seeds["classifier"]), "boot_seed": int(seeds["bootstrap"]), "perm_seed": int(seeds["global"])}


# ============================================================================ data
def eval_splits(a, reg, meta) -> dict:
    """Evaluated splits: {name: {"idx": global rows, "group": bootstrap group per cell, "groups": order}}."""
    split, sites, donors = meta["split"], meta["sites"], meta["donors"]
    if a.smoke:
        va_idx = vk.split_indices(split, sites, donors, reg, "val")
        rng = np.random.default_rng(SMOKE["seed"])
        pick = np.sort(rng.choice(va_idx, size=SMOKE["n_primary"] + SMOKE["n_secondary"], replace=False))
        perm = rng.permutation(pick.size)
        prim = np.sort(pick[perm[: SMOKE["n_primary"]]])
        sec = np.sort(pick[perm[SMOKE["n_primary"]:]])
        grp = np.where(np.random.default_rng(SMOKE["seed"] + 1).random(prim.size) < 0.5, "valA", "valB")
        return {"primary": {"name": "smoke_val_primary", "idx": prim, "group": grp, "groups": ["valA", "valB"]},
                "secondary": {"name": "smoke_val_secondary", "idx": sec, "group": np.full(sec.size, "valC"),
                              "groups": ["valC"]}}
    out = {}
    for role, name in (("primary", "test_primary"), ("secondary", "test_secondary")):
        idx = vk.split_indices(split, sites, donors, reg, name)
        d = donors[idx].astype(str)
        out[role] = {"name": name, "idx": idx, "group": d, "groups": sorted(set(d))}
    return out


def load_meta(a) -> dict:
    npz = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
    return {"split": npz["split"].astype(str), "sites": npz["sites"].astype(str), "donors": npz["donors"].astype(str),
            "adt_names": [str(x) for x in npz["adt_names"]], "npz": npz}


def q3_annotation_key(annot: np.ndarray, ct: np.ndarray, reg) -> np.ndarray:
    spec = reg["questions"]["Q3"]["key"]
    out = np.asarray(annot).astype(object).copy()
    my = out == "myeloid"
    cl = np.isin(np.asarray(ct).astype(str), spec["classical_types"])
    out[my & cl] = "myeloid"
    out[my & ~cl] = "OUT"
    return out.astype(str)


def key_validity(keys: dict, ct: np.ndarray, group: np.ndarray, groups: list, reg) -> dict:
    rep = {"all": vk.agreement_report(keys["annotation"], keys["gated"]), "per_group": {}, "per_type": {}}
    for g in groups:
        m = group == g
        rep["per_group"][g] = vk.agreement_report(keys["annotation"][m], keys["gated"][m])
    rep["per_batch_gate_vs_annotation"] = vk.agreement_report(keys["annotation"], keys["gated_per_batch"])
    for t in sorted(set(ct)):
        m = ct == t
        rep["per_type"][t] = {"class": reg["annotation_map"][t], "n": int(m.sum()),
                              "gated": {c: int(np.sum(keys["gated"][m] == c)) for c in vk.CLASSES},
                              "in_primary_key": int(np.sum(keys["primary"][m] != vk.UNSCORED))}
    rep["key_counts"] = {k: {c: int(np.sum(keys[k] == c)) for c in list(vk.CLASSES) + [vk.UNSCORED]} for k in K1 + K3}
    kap = rep["all"]["kappa5"]
    rep["kappa_flag"] = {"threshold": 0.85, "kappa5": kap, "below": bool(kap is not None and kap < 0.85),
                         "consequence": "every result is reported with equal prominence on the annotation-only key "
                                        "(the primary key is not switched)"}
    return rep


# ============================================================================ stage: prepare
def stage_prepare(a, reg, cfg, fp, splits, meta) -> None:
    cache = a.out_dir / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    done = [cache / f"prepare_{r}.npz" for r in splits]
    if all(p.exists() and str(np.load(p)["fingerprint"]) == fp for p in done) and (cache / "prepare_meta.json").exists():
        say("prepare: cached")
        return
    t0 = time.time()
    amend = reg["amendment_A1"]
    addendum = json.loads(a.addendum.read_text())["computed"] if a.addendum.exists() else None
    if addendum is None:
        raise SystemExit(f"E1 addendum missing: {a.addendum} (run bridge_anm/v3_e1_addendum.py first)")
    npz = meta["npz"]
    names = meta["adt_names"]
    j = {n: i for i, n in enumerate(names)}
    split, sites, donors = meta["split"], meta["sites"], meta["donors"]
    tr_idx = vk.split_indices(split, sites, donors, reg, "train")
    ct_all = npz["cell_types"].astype(str)
    adt_all = np.asarray(npz["adt"], dtype=np.float32)
    zmm = np.load(a.z, mmap_mode="r")

    # ---- keys and key validity first (verifier side; before any method result)
    keys_ev, kv = {}, {}
    for role, sp in splits.items():
        idx = sp["idx"]
        k = vk.build_keys(ct_all[idx], adt_all[idx], names, reg, sites=sites[idx], donors=donors[idx])
        k["annotation_only"] = k["annotation"]
        k["q3_annotation_only"] = q3_annotation_key(k["annotation"], ct_all[idx], reg)
        keys_ev[role] = k
        kv[role] = {"split": sp["name"], **key_validity(k, ct_all[idx], sp["group"], sp["groups"], reg)}
        say(f"key validity {sp['name']}: n {idx.size}, agreement {kv[role]['all']['agreement']:.4f}, "
            f"kappa {kv[role]['all']['kappa5']:.4f}, flag {kv[role]['kappa_flag']['below']}")
    write_json(a.out_dir / "key_validity.json", {"registration_sha256": vk.sha256_file(a.registration),
                                                 "smoke": cfg["smoke"], "splits": kv})

    # ---- evidence (same for every method)
    def evid(idx):
        return vk.evidence(vk.head_predict(np.asarray(zmm[idx], dtype=np.float32), a.ckpt, device="cpu", size_factor=1.0), reg)

    v_tr = evid(tr_idx)
    tr_keys = vk.build_keys(ct_all[tr_idx], adt_all[tr_idx], names, reg)
    prim_tr, q3_tr = tr_keys["primary"], tr_keys["q3"]
    say(f"training cells {tr_idx.size}: primary-key labelled {int((prim_tr != vk.UNSCORED).sum())}, "
        f"q3 labelled {int((q3_tr != vk.UNSCORED).sum())}")

    # ---- classifiers (training-site labels only)
    a1c = reg["amendment_A1"]["computed"]["classifiers"]
    f12 = list(reg["classifier"]["primary"]["feature_proteins"])
    f14 = list(reg["classifier"]["q3"]["feature_proteins"])
    fanch = list(a1c["q1labels_anchors"]["features"])
    specs = {
        "primary": (f12, prim_tr, float(reg["classifier"]["primary"]["C"])),
        "q1labels_q3features": (list(a1c["q1labels_q3features"]["features"]), prim_tr, float(a1c["q1labels_q3features"]["C"])),
        "q3labels": (f14, q3_tr, float(reg["classifier"]["q3"]["C"])),
        "q1labels_anchors": (fanch, prim_tr, float(a1c["q1labels_anchors"]["C"])),
        "q3labels_anchors": (list(a1c["q3labels_anchors"]["features"]), q3_tr, float(a1c["q3labels_anchors"]["C"])),
    }
    for pn, ps in addendum["sensitivity_panels"].items():
        specs[f"sens_{pn}"] = (list(ps["classifier"]["features"]), prim_tr, float(ps["classifier"]["C"]))
    models, clf_meta = {}, {}
    for nm, (feats, y, C) in specs.items():
        m = y != vk.UNSCORED
        models[nm] = (e1.fit_logreg(v_tr[m][:, [j[p] for p in feats]], y[m], C), feats)
        clf_meta[nm] = {"features": feats, "C": C, "n_train_labelled": int(m.sum()),
                        "classes": [str(c) for c in models[nm][0].classes_]}
    if clf_meta["primary"]["n_train_labelled"] != int(reg["classifier"]["primary"]["n_train_labelled"]):
        raise SystemExit("primary classifier training set differs from the registration")
    say(f"classifiers fitted: {list(models)}")

    from sklearn.neural_network import MLPClassifier

    m12 = prim_tr != vk.UNSCORED
    X12 = v_tr[:, [j[p] for p in f12]]
    mlps = []
    for s in cfg["mlp_seeds"]:
        mlp = MLPClassifier(hidden_layer_sizes=(64,), alpha=1e-4, early_stopping=True, validation_fraction=0.1,
                            max_iter=500, random_state=int(s)).fit(X12[m12], prim_tr[m12])
        mlps.append(mlp)
        say(f"MLP seed {s}: {mlp.n_iter_} iterations")

    # ---- label-cost classifiers (E1.5): Q3 curve (14 features), anchors curve, Q1 reference curve
    pool_q3 = np.where(q3_tr != vk.UNSCORED)[0]   # positions; monotone in global cell id
    pool_q1 = np.where(prim_tr != vk.UNSCORED)[0]
    curves = {
        "q3": {"features": f14, "labels": q3_tr, "pool": pool_q3, "C": float(reg["classifier"]["q3"]["C"]),
               "n0": "q1labels_q3features", "nall": "q3labels", "ns": list(cfg["lc_n"])},
        "anch": {"features": fanch, "labels": q3_tr, "pool": pool_q3, "C": float(a1c["q3labels_anchors"]["C"]),
                 "n0": "q1labels_anchors", "nall": "q3labels_anchors", "ns": list(cfg["lc_n"])},
        "q1ref": {"features": f12, "labels": prim_tr, "pool": pool_q1, "C": float(reg["classifier"]["primary"]["C"]),
                  "n0": None, "nall": "primary", "ns": [n for n in cfg["lc_n"] if n != 0]},
    }
    lc_models: dict = {}
    lc_meta: dict = {}
    for cn, c in curves.items():
        Xc = v_tr[:, [j[p] for p in c["features"]]]
        lc_meta[cn] = {"features": c["features"], "C": c["C"], "ns": [str(n) for n in c["ns"]],
                       "draw_seeds": cfg["lc_draws"], "pool_size": int(c["pool"].size), "draws": {}}
        for n in c["ns"]:
            for d in cfg["lc_draws"]:
                if n == 0:
                    mdl = models[c["n0"]][0]
                    cls_present = clf_meta[c["n0"]]["classes"]
                elif n == "all":
                    mdl = models[c["nall"]][0]
                    cls_present = clf_meta[c["nall"]]["classes"]
                else:
                    pos = e1.nested_draw(c["pool"], n, d)
                    mdl = e1.fit_logreg(Xc[pos], c["labels"][pos], c["C"])
                    cls_present = [str(x) for x in mdl.classes_]
                lc_models[(cn, str(n), int(d))] = mdl
                lc_meta[cn]["draws"][f"n{n}_d{d}"] = {"classes_in_draw": cls_present}
        say(f"label-cost curve {cn}: {len(c['ns'])} n x {len(cfg['lc_draws'])} draws")

    # ---- per evaluated split: scores and classifier outputs
    pan_q1 = vk.question_panel(reg, "Q1")
    for role, sp in splits.items():
        idx = sp["idx"]
        v = evid(idx)
        S1 = vk.class_scores(v, names, reg, "Q1")
        S3 = vk.class_scores(v, names, reg, "Q3")
        arr = {"fingerprint": np.asarray(fp), "ids": idx.astype(np.int64), "group": sp["group"].astype(str),
               "ct": ct_all[idx], "donor": donors[idx], "v": v, "S1": S1, "S3": S3,
               "closure_recoded": vk.anm_action_scores(v, names, reg, "Q1", readout="closure")}
        for kn in KEY_NAMES:
            arr[f"key_{kn}"] = keys_ev[role][kn].astype(str)
        for pn in addendum["sensitivity_panels"]:
            arr[f"S_sens_{pn}"] = vk.class_scores(v, names, reg, "Q1", panel=pn)
        for nm, (mdl, feats) in models.items():
            arr[f"P5_{nm}"] = e1.proba5(mdl, v[:, [j[p] for p in feats]])
        mc, mk = [], []
        for mlp in mlps:
            c_, k_ = e1.lineage_conf_call(e1.proba5(mlp, v[:, [j[p] for p in f12]]))
            mc.append(c_)
            mk.append(np.asarray([vk.LINEAGES.index(x) for x in k_], dtype=np.int8))
        arr["mlp_conf"], arr["mlp_call"] = np.asarray(mc), np.asarray(mk)
        for cn, c in curves.items():
            Xe = v[:, [j[p] for p in c["features"]]]
            conf = np.zeros((len(c["ns"]), len(cfg["lc_draws"]), idx.size))
            call = np.zeros(conf.shape, dtype=np.int8)
            for ni, n in enumerate(c["ns"]):
                for di, d in enumerate(cfg["lc_draws"]):
                    c_, k_ = e1.lineage_conf_call(e1.proba5(lc_models[(cn, str(n), int(d))], Xe))
                    conf[ni, di], call[ni, di] = c_, [vk.LINEAGES.index(x) for x in k_]
            arr[f"lc_{cn}_conf"], arr[f"lc_{cn}_call"] = conf, call
        atomic_savez(cache / f"prepare_{role}.npz", **arr)
        say(f"prepare {sp['name']}: {idx.size} cells cached")
    write_json(cache / "prepare_meta.json", {"fingerprint": fp, "classifiers": clf_meta, "label_cost": lc_meta,
                                             "mlp_seeds": cfg["mlp_seeds"], "mlp_iterations": [int(m.n_iter_) for m in mlps],
                                             "panel_q1": pan_q1, "seconds": time.time() - t0})
    say(f"prepare done in {time.time() - t0:.0f}s")


# ============================================================================ stage: anm
def stage_anm(a, reg, cfg, fp, splits, chunk_cells: int = 2000) -> None:
    cache = a.out_dir / "cache"
    bridge = e1.AnmBridge(reg, SCHEMA, a.anm_root)
    names = reg["evidence"]["teddy_head"]["adt_names"]
    j = {n: i for i, n in enumerate(names)}
    p1, p3 = vk.question_panel(reg, "Q1"), vk.question_panel(reg, "Q3")
    for role, sp in splits.items():
        out = cache / f"anm_{role}.npz"
        if out.exists() and str(np.load(out)["fingerprint"]) == fp:
            say(f"anm {sp['name']}: cached")
            continue
        P = np.load(cache / f"prepare_{role}.npz")
        v = P["v"]
        N = v.shape[0]
        parts = []
        for c0 in range(0, N, chunk_cells):
            cf = cache / f"anm_{role}_part{c0 // chunk_cells:04d}.npz"
            if cf.exists() and str(np.load(cf)["fingerprint"]) == fp:
                parts.append(cf)
                continue
            t0 = time.time()
            rows = range(c0, min(N, c0 + chunk_cells))
            n = len(rows)
            R = {"q1_scores": np.zeros((n, 4)), "q1_call": np.full(n, -1, np.int8), "q2_call": np.full(n, -1, np.int8),
                 "q3_scores": np.zeros((n, 4)), "q3_call": np.full(n, -1, np.int8), "cl_scores": np.zeros((n, 4)),
                 "cl_q1_call": np.full(n, -1, np.int8), "cl_q2_call": np.full(n, -1, np.int8),
                 "n_rejected": np.zeros(n, np.int16), "loo_after": np.full((n, 3), np.nan),
                 "loo_q1_call": np.full((n, 3), -2, np.int8), "loo_q2_call": np.full((n, 3), -2, np.int8)}
            for r, i in enumerate(rows):
                vals1 = {k: [v[i, j[p]] for p in p1[k]] for k in vk.LINEAGES}
                o = bridge.run(p1, vals1, readouts=("Q1", "Q2"), closure_readouts=("closure_Q1", "closure_Q2"))
                R["q1_scores"][r], R["q1_call"][r] = o["Q1"]
                R["q2_call"][r] = o["Q2"][1]
                R["cl_scores"][r], R["cl_q1_call"][r] = o["closure_Q1"]
                R["cl_q2_call"][r] = o["closure_Q2"][1]
                o3 = bridge.run(p3, {k: [v[i, j[p]] for p in p3[k]] for k in vk.LINEAGES}, readouts=("Q3",))
                R["q3_scores"][r], R["q3_call"][r] = o3["Q3"]
                R["n_rejected"][r] = o["n_rejected"] + o3["n_rejected"]
                c = int(o["Q1"][1])
                if c >= 0:  # leave-one-out over the called class's events (Q2 calls are Q1 calls of the same class)
                    k = vk.LINEAGES[c]
                    for pi, p in enumerate(p1[k]):
                        od = bridge.run(p1, vals1, readouts=("Q1", "Q2"), drop=(k, p))
                        R["loo_after"][r, pi] = od["Q1"][0][c]
                        R["loo_q1_call"][r, pi] = od["Q1"][1]
                        R["loo_q2_call"][r, pi] = od["Q2"][1]
            atomic_savez(cf, fingerprint=np.asarray(fp), **R)
            parts.append(cf)
            say(f"anm {sp['name']}: cells {c0}-{c0 + n - 1} of {N} in {time.time() - t0:.1f}s")
        merged = {k: np.concatenate([np.load(p)[k] for p in parts]) for k in np.load(parts[0]).files if k != "fingerprint"}
        atomic_savez(out, fingerprint=np.asarray(fp), anm_commit=np.asarray(str(bridge.anm_commit())), **merged)
        for p in parts:
            p.unlink()
        say(f"anm {sp['name']}: done ({N} cells)")


# ============================================================================ evaluation model
class Method:
    def __init__(self, name, conf, call_idx, ids, amend, acc_covs=(), acc_keys=(), c2_covs=(), c2_keys=()):
        self.name = name
        self.order = e1.ranked_order(conf, ids, amend)
        self.call = LIN[np.asarray(call_idx, dtype=np.int64)]
        self.acc_covs, self.acc_keys = list(acc_covs), list(acc_keys)
        self.c2_covs, self.c2_keys = list(c2_covs), list(c2_keys)
        self.covs = list(dict.fromkeys(self.acc_covs + self.c2_covs))
        self.arr: dict = {}

    def bind(self, keys: dict):
        for kn in set(self.acc_keys) | set(self.c2_keys):
            key = keys[kn][self.order]
            call = self.call[self.order]
            self.arr[kn] = {"scored": (key != vk.UNSCORED).astype(np.float64),
                            "correct": (call == key).astype(np.float64),
                            "out": (key == "OUT").astype(np.float64),
                            "inscope": np.isin(key, vk.LINEAGES).astype(np.float64)}


class SplitEval:
    """Every bootstrapped quantity of E1 for one evaluated split, as a function of a count vector w."""

    def __init__(self, P, A, reg, cfg, addendum):
        self.reg, self.cfg = reg, cfg
        amend = reg["amendment_A1"]
        self.grid = list(reg["experiments"]["common"]["matched_coverage"]["grid"])
        self.gn = [f"g{c:.2f}" for c in self.grid]
        ids = P["ids"]
        self.N = ids.size
        self.keys = {kn: P[f"key_{kn}"].astype(str) for kn in KEY_NAMES}
        S1, S3 = P["S1"], P["S3"]
        self.called1 = (S1.max(axis=1) >= float(reg["questions"]["Q1"]["bar"])).astype(np.float64)
        self.called3 = (S3.max(axis=1) >= float(reg["questions"]["Q3"]["bar"])).astype(np.float64)
        rc1 = S1.argmax(axis=1)
        mg = vk.margin(S1)
        M = []
        add = M.append
        allg = self.gn
        c2c = ["c1", "g0.80", "g0.70"]
        add(Method("top", S1.max(axis=1), rc1, ids, amend, allg + ["c1"], K1, c2c, K1))
        add(Method("margin", mg, rc1, ids, amend, allg, K1))
        add(Method("entropy", e1.entropy_confidence(S1), rc1, ids, amend, allg, K1))
        cconf, ccall = e1.lineage_conf_call(P["P5_primary"])
        cidx = np.asarray([vk.LINEAGES.index(x) for x in ccall])
        add(Method("classifier", cconf, cidx, ids, amend, allg, K1, c2c, K1))
        add(Method("classifier_rule_calls", cconf, rc1, ids, amend, allg, ("primary",)))
        add(Method("closure", A["cl_scores"].max(axis=1), A["cl_scores"].argmax(axis=1), ids, amend,
                   ["c1"], K1, c2c, K1))
        for si in range(P["mlp_conf"].shape[0]):
            add(Method(f"mlp_s{si}", P["mlp_conf"][si], P["mlp_call"][si], ids, amend, allg, ("primary",)))
        for pn in addendum["sensitivity_panels"]:
            Sp = P[f"S_sens_{pn}"]
            add(Method(f"top@{pn}", Sp.max(axis=1), Sp.argmax(axis=1), ids, amend, allg, ("primary",)))
            add(Method(f"margin@{pn}", vk.margin(Sp), Sp.argmax(axis=1), ids, amend, allg, ("primary",)))
            add(Method(f"entropy@{pn}", e1.entropy_confidence(Sp), Sp.argmax(axis=1), ids, amend, allg, ("primary",)))
            pc, pk = e1.lineage_conf_call(P[f"P5_sens_{pn}"])
            add(Method(f"classifier@{pn}", pc, [vk.LINEAGES.index(x) for x in pk], ids, amend, allg, ("primary",)))
        add(Method("rule_q3", S3.max(axis=1), S3.argmax(axis=1), ids, amend, ["c3"], K3))
        for nm in ("q1labels_q3features", "q3labels", "q1labels_anchors", "q3labels_anchors"):
            pc, pk = e1.lineage_conf_call(P[f"P5_{nm}"])
            add(Method(f"clf_{nm}", pc, [vk.LINEAGES.index(x) for x in pk], ids, amend, ["c3"], K3))
        self.lc = {}
        for cn, cov, keys in (("q3", "c3", K3), ("anch", "c3", K3), ("q1ref", "c1", ("primary",))):
            conf, call = P[f"lc_{cn}_conf"], P[f"lc_{cn}_call"]
            ns = [n for n in self.cfg["lc_n"] if not (cn == "q1ref" and n == 0)]
            self.lc[cn] = {"ns": ns, "cov": cov, "keys": keys}
            for ni, n in enumerate(ns):
                for di in range(conf.shape[1]):
                    add(Method(f"lc_{cn}_n{n}_d{di}", conf[ni, di], call[ni, di], ids, amend, [cov], keys))
        for m in M:
            m.bind(self.keys)
        self.methods = M
        self.n_out_keys = {kn: (self.keys[kn] == "OUT").astype(np.float64) for kn in K1}
        # E1.3c: flip-sensitive calls (exact closed form, A1.3) vs as many lowest-margin calls
        fld = reg["anm"]["field_representation"]
        names = reg["evidence"]["teddy_head"]["adt_names"]
        jj = {n: i for i, n in enumerate(names)}
        pan = vk.question_panel(reg, "Q1")
        vcal = np.zeros((self.N, 3))
        for k, lin in enumerate(vk.LINEAGES):
            r = rc1 == k
            vcal[r] = P["v"][r][:, [jj[p] for p in pan[lin]]]
        self.e13 = {}
        rule_call = LIN[rc1]
        for q in ("Q1", "Q2"):
            bar = float(reg["questions"][q]["bar"])
            called = S1.max(axis=1) >= bar
            flip = np.zeros(self.N, dtype=bool)
            flip[called] = va.loo_flip_exact(S1[called], vcal[called], rc1[called], bar, fld)
            low = e1.select_lowest(mg, called, ids, amend)
            ent = {"flip": flip.astype(np.float64), "low": low, "keys": {}}
            for kn in K1:
                key = self.keys[kn]
                sc = (key != vk.UNSCORED)
                wr = sc & (rule_call != key)
                ent["keys"][kn] = (sc.astype(np.float64), wr.astype(np.float64))
            self.e13[q] = ent
        self.vcal, self.rc1 = vcal, rc1

    def stats(self, w: np.ndarray) -> dict:
        w = np.asarray(w, dtype=np.float64)
        Nb = w.sum()
        cov = dict(zip(self.gn, self.grid))
        cov["c1"] = float((w * self.called1).sum() / Nb)
        cov["c3"] = float((w * self.called3).sum() / Nb)
        out = {"N": Nb, "cov|c1": cov["c1"], "cov|c3": cov["c3"]}
        n_out_tot = {kn: float((w * a).sum()) for kn, a in self.n_out_keys.items()}
        for m in self.methods:
            ks = [e1.coverage_k(cov[c], Nb) for c in m.covs]
            s = e1.select_ordered(w[m.order], ks)
            ci = {c: i for i, c in enumerate(m.covs)}
            for kn in m.acc_keys:
                A = m.arr[kn]
                nsc, nco = s @ A["scored"], s @ A["correct"]
                for c in m.acc_covs:
                    out[f"{m.name}|{c}|{kn}|acc"] = ratio(nco[ci[c]], nsc[ci[c]])
            for kn in m.c2_keys:
                A = m.arr[kn]
                no, ni, nco = s @ A["out"], s @ A["inscope"], s @ A["correct"]
                for c in m.c2_covs:
                    out[f"{m.name}|{c}|{kn}|outdecl"] = 1.0 - ratio(no[ci[c]], n_out_tot[kn])
                    out[f"{m.name}|{c}|{kn}|inscope"] = ratio(nco[ci[c]], ni[ci[c]])
        for q, ent in self.e13.items():
            wf = w * ent["flip"]
            mflip = float(wf.sum())
            low = ent["low"]
            sl = e1.select_ordered(w[low], [mflip])[0]
            out[f"e13c|{q}|n_flip"] = mflip
            for kn, (sc, wr) in ent["keys"].items():
                ef = ratio((wf * wr).sum(), (wf * sc).sum())
                el = ratio((sl * wr[low]).sum(), (sl * sc[low]).sum())
                out[f"e13c|{q}|{kn}|err_flip"] = ef
                out[f"e13c|{q}|{kn}|err_low"] = el
        return out


# ============================================================================ stage: boot
def schemes(sp: dict, role: str) -> dict:
    """Bootstrap schemes for a split: pooled two-stage (primary only) and one-stage per group."""
    grp = sp["group"]
    groups = [np.where(grp == g)[0] for g in sp["groups"]]
    out = {}
    if role == "primary" and len(groups) > 1:
        out["pooled"] = (groups, True)
    for g, gi in zip(sp["groups"], groups):
        out[f"group_{g}"] = ([gi], False)
    return out


def stage_boot(a, reg, cfg, fp, splits, addendum) -> None:
    bdir = a.out_dir / "boot"
    bdir.mkdir(parents=True, exist_ok=True)
    cache = a.out_dir / "cache"
    for role, sp in splits.items():
        P = np.load(cache / f"prepare_{role}.npz")
        A = np.load(cache / f"anm_{role}.npz")
        sp = {**sp, "group": P["group"].astype(str)}
        ev = SplitEval(P, A, reg, cfg, addendum)
        point = ev.stats(np.ones(ev.N))
        names = sorted(point)
        (bdir / f"{role}_names.json").write_text(json.dumps(names))
        for sname, (groups, two) in schemes(sp, role).items():
            final = bdir / f"{role}_{sname}.npy"
            if final.exists() and (bdir / f"{role}_{sname}.fp").exists() and (bdir / f"{role}_{sname}.fp").read_text() == fp:
                say(f"boot {sp['name']} {sname}: cached")
                continue
            B, ch = int(cfg["n_boot"]), int(a.chunk)
            gen = e1.replicate_counts(groups, ev.N, B, cfg["boot_seed"], two_stage=two)
            t0 = time.time()
            done_rows = []
            for c0 in range(0, B, ch):
                cf = bdir / f"{role}_{sname}_chunk{c0 // ch:04d}.npy"
                fpf = cf.with_suffix(".fp")
                ws = [next(gen) for _ in range(min(ch, B - c0))]   # always drawn: the sequence is the same on resume
                if cf.exists() and fpf.exists() and fpf.read_text() == fp:
                    done_rows.append(np.load(cf))
                    continue
                M = np.asarray([[st[k] for k in names] for st in (ev.stats(w) for w in ws)])
                np.save(cf, M)
                fpf.write_text(fp)
                done_rows.append(M)
                el = time.time() - t0
                n_done = c0 + len(ws)
                say(f"boot {sp['name']} {sname}: {n_done}/{B} replicates, {el:.0f}s elapsed")
            np.save(final, np.concatenate(done_rows, axis=0))
            (bdir / f"{role}_{sname}.fp").write_text(fp)
            for f in bdir.glob(f"{role}_{sname}_chunk*"):
                f.unlink()
            say(f"boot {sp['name']} {sname}: done")


# ============================================================================ stage: report
class Stats:
    """Point / per-group / bootstrap views of one split's statistics."""

    def __init__(self, point: dict, groups: dict, boots: dict, names: list):
        self.point, self.groups, self.names = point, groups, names
        self.boots = {s: {k: M[:, i] for i, k in enumerate(names)} for s, M in boots.items()}


def aurc_of(get, method: str, key: str, gn: list, grid: list):
    """AURC (registered: trapezoid over the grid / 0.45) for a point (scalars) or all replicates (arrays)."""
    cols = np.stack([np.asarray(get(f"{method}|{g}|{key}|acc"), dtype=np.float64) for g in gn])
    dg = np.abs(np.diff(np.asarray(grid, dtype=np.float64))).reshape((-1,) + (1,) * (cols.ndim - 1))
    return np.sum(0.5 * (cols[1:] + cols[:-1]) * dg, axis=0) / (max(grid) - min(grid))


def comparison(st: Stats, fn, margin: float, primary_scheme: str | None, group_names: list) -> dict:
    """Registered comparison: point, bootstrap CI (pooled two-stage, or one-stage for a one-donor split),
    per-donor points and one-stage per-donor CIs, verdict."""
    point = float(fn(lambda k: st.point[k]))
    per = {g: float(fn(lambda k, g=g: st.groups[g][k])) for g in group_names}
    scheme = primary_scheme or f"group_{group_names[0]}"
    boot = np.asarray(fn(lambda k: st.boots[scheme][k]), dtype=np.float64)
    out = e1.verdict(point, boot, per, margin)
    out["per_donor_ci95"] = {}
    for g in group_names:
        bg = st.boots.get(f"group_{g}")
        if bg is not None:
            lo, hi, _ = e1.percentile_ci(np.asarray(fn(lambda k, bg=bg: bg[k]), dtype=np.float64))
            out["per_donor_ci95"][g] = [lo, hi]
    out["bootstrap"] = scheme
    return out


def value_block(st: Stats, fn, primary_scheme: str | None, group_names: list) -> dict:
    point = fn(lambda k: st.point[k])
    scheme = primary_scheme or f"group_{group_names[0]}"
    lo, hi, n = e1.percentile_ci(np.asarray(fn(lambda k: st.boots[scheme][k]), dtype=np.float64))
    return {"point": point, "ci95": [lo, hi], "per_donor": {g: fn(lambda k, g=g: st.groups[g][k]) for g in group_names}}


def call_table(calls: np.ndarray, ct: np.ndarray) -> dict:
    out = {}
    for t in sorted(set(ct)):
        m = ct == t
        out[t] = {"n": int(m.sum()), **{c: int(np.sum(calls[m] == c)) for c in vk.LINEAGES},
                  "no_call": int(np.sum(calls[m] == vk.NO_CALL))}
    return out


def operating_point(calls: np.ndarray, keys: dict, key_names, ct, group, groups) -> dict:
    out = {"coverage": float(np.mean(calls != vk.NO_CALL))}
    for kn in key_names:
        key = keys[kn]
        r = vk.selective_accuracy(calls, key)
        r["per_class_recall"] = {c: (float(np.mean(calls[key == c] == (c if c != "OUT" else vk.NO_CALL)))
                                     if np.any(key == c) else None) for c in vk.CLASSES}
        uns = key == vk.UNSCORED
        r["calls_on_unscored"] = {"n_unscored": int(uns.sum()), **{c: int(np.sum(calls[uns] == c)) for c in vk.LINEAGES}}
        if kn in ("primary", "q3"):
            r["call_table_unscored_per_type"] = call_table(calls[uns], ct[uns])
        r["per_donor"] = {g: vk.selective_accuracy(calls[group == g], key[group == g]) for g in groups}
        out[kn] = r
    out["call_table_per_type"] = call_table(calls, ct)
    return out


def e13_block(P, A, reg, ev: SplitEval, group, groups, cfg, question: str) -> dict:
    """E1.3a (engine vs closed form, sets and flips) and E1.3b (marker shares vs permutation null)."""
    S1 = P["S1"]
    bar = float(reg["questions"][question]["bar"])
    called = np.where(S1.max(axis=1) >= bar)[0]
    ci = ev.rc1[called]
    vcal = ev.vcal[called]
    fld = reg["anm"]["field_representation"]
    pan = vk.question_panel(reg, "Q1")
    base = A["q1_scores"][called, ci]
    dec = base[:, None] - A["loo_after"][called]
    eng_sets = e1.engine_top_sets(dec)
    cf_sets = e1.closed_form_top_sets(vcal)
    loo_call = A["loo_q1_call" if question == "Q1" else "loo_q2_call"][called]
    first_top = np.argmax(eng_sets, axis=1)
    eng_new = loo_call[np.arange(called.size), first_top]
    eng_flip = eng_new != ci
    # consistency among tied top events (they must give the same call)
    tied_inconsistent = int(np.sum([len(set(loo_call[i][eng_sets[i]].tolist())) > 1 for i in range(called.size)]))
    cf_flip = va.loo_flip_exact(S1[called], vcal, ci, bar, fld)
    out = {"question": question, "n_calls": int(called.size),
           "E1.3a": {"top_set_agreement": float(np.mean(np.all(eng_sets == cf_sets, axis=1))) if called.size else None,
                     "flip_agreement": float(np.mean(eng_flip == cf_flip)) if called.size else None,
                     "share_calls_multi_marker_set": float(np.mean(cf_sets.sum(axis=1) > 1)) if called.size else None,
                     "n_flip_closed_form": int(cf_flip.sum()), "n_flip_engine": int(eng_flip.sum()),
                     "tied_top_events_with_different_calls": tied_inconsistent,
                     "registered_form_flips (zeroing, superseded by A1.3)": int(va.loo_flip_registered_form(S1[called], vcal, ci, bar).sum()),
                     "expectation": "1.000; below 0.99 is reported as a discrepancy"}}
    a13 = out["E1.3a"]
    a13["discrepancy"] = bool(called.size and (a13["top_set_agreement"] < 0.99 or a13["flip_agreement"] < 0.99))
    # E1.3b shares with fractional credit, per donor, vs the within-class permutation null (engine sets)
    shares = {}
    grp = group[called]
    for g in list(groups) + ["all"]:
        sel = np.ones(called.size, bool) if g == "all" else grp == g
        blk = {}
        for subset, sm in (("all_calls", sel), ("single_marker_calls", sel & (eng_sets.sum(axis=1) == 1))):
            if not sm.any():
                blk[subset] = None
                continue
            obs = e1.marker_shares(eng_sets[sm], ci[sm])
            null = e1.marker_share_null(eng_sets[sm], ci[sm], cfg["n_perm"], cfg["perm_seed"])
            tab = {}
            for k, lin in enumerate(vk.LINEAGES):
                if obs["n_calls_per_class"][k] == 0:
                    continue
                for pi, p in enumerate(pan[lin]):
                    nw = null["within_class"][:, k, pi]
                    lo, hi, _ = e1.percentile_ci(nw)
                    tab[f"{lin}:{p}"] = {"share_all_calls": obs["all"][k, pi], "share_within_class": obs["within_class"][k, pi],
                                         "null_mean_within_class": float(np.nanmean(nw)), "null_ci95_within_class": [lo, hi],
                                         "p_two_sided": e1.perm_p_two_sided(obs["within_class"][k, pi], nw)}
            blk[subset] = {"n_calls": int(sm.sum()), "n_calls_per_class": dict(zip(vk.LINEAGES, obs["n_calls_per_class"].tolist())),
                           "markers": tab}
        shares[g] = blk
    out["E1.3b"] = {"null": f"within-class permutation of marker identities per cell, {cfg['n_perm']} permutations, "
                            f"seed {cfg['perm_seed']} (seeds.global)", "credit": "fractional (1/m per marker of a set of m)",
                    "sets_from": "ANM engine leave-one-out", "per_donor": shares}
    if not np.all(eng_sets == cf_sets):
        out["E1.3b"]["closed_form_sets_differ"] = True
    return out


def stage_report(a, reg, cfg, fp, splits, addendum, st_reg: dict) -> dict:
    cache, bdir = a.out_dir / "cache", a.out_dir / "boot"
    grid = list(reg["experiments"]["common"]["matched_coverage"]["grid"])
    gn = [f"g{c:.2f}" for c in grid]
    E = reg["experiments"]["E1"]
    kv = json.loads((a.out_dir / "key_validity.json").read_text())["splits"]
    pmeta = json.loads((cache / "prepare_meta.json").read_text())
    res: dict = {"experiment": "E1", "builder": VERSION, "smoke": cfg["smoke"],
                 "registration_sha256": st_reg["registration_sha256"], "amendment_A1_sha256": st_reg["amendment_A1_sha256"],
                 "addendum_E1_sha256": st_reg.get("addendum_sha256"), "state": st_reg, "config": cfg,
                 "key_validity": kv, "splits": {}}
    if cfg["smoke"]:
        res["SMOKE"] = "val donor 18303 cells with label-free pseudo-donors and small B; not a result"
    for role, sp in splits.items():
        P = np.load(cache / f"prepare_{role}.npz")
        A = np.load(cache / f"anm_{role}.npz")
        group = P["group"].astype(str)
        groups = sp["groups"]
        ev = SplitEval(P, A, reg, cfg, addendum)
        names = json.loads((bdir / f"{role}_names.json").read_text())
        point = ev.stats(np.ones(ev.N))
        gpts = {g: ev.stats((group == g).astype(np.float64)) for g in groups}
        boots = {s: np.load(bdir / f"{role}_{s}.npy") for s in schemes({**sp, "group": group}, role)}
        st = Stats(point, gpts, boots, names)
        ps = "pooled" if "pooled" in boots else None
        keys = ev.keys
        ct = P["ct"].astype(str)
        S1, S3 = P["S1"], P["S3"]
        R: dict = {"split": sp["name"], "n_cells": int(ev.N), "groups": {g: int(np.sum(group == g)) for g in groups},
                   "anm_commit": str(A["anm_commit"])}

        # ---------------- E1.1a: ANM engine vs re-coded rule, cell by cell
        rq = {q: vk.rule_calls(S1 if q != "Q3" else S3, reg["questions"][q]["bar"]) for q in ("Q1", "Q2", "Q3")}
        eng = {q: np.where(A[f"{q.lower()}_call"] >= 0, LIN[np.maximum(A[f"{q.lower()}_call"], 0)], vk.NO_CALL)
               for q in ("Q1", "Q2", "Q3")}
        G = reg["anm"]["field_gain"]
        clr = P["closure_recoded"]
        cl_rule = {q: vk.rule_calls(clr, reg["anm"][f"closure_bar_{q}"]) for q in ("Q1", "Q2")}
        cl_eng = {q: np.where(A[f"cl_{q.lower()}_call"] >= 0, LIN[np.maximum(A[f"cl_{q.lower()}_call"], 0)], vk.NO_CALL)
                  for q in ("Q1", "Q2")}
        e11a = {"registered_value": 0,
                "mismatches": {q: int(np.sum(eng[q] != rq[q])) for q in rq},
                "closure_mismatches": {q: int(np.sum(cl_eng[q] != cl_rule[q])) for q in cl_rule},
                "max_rel_diff_engine_vs_Gn_times_S": {
                    "Q1": float(np.max(np.abs(A["q1_scores"] - vk.field_gain(3, reg["anm"]["field_representation"]) * 3 * S1)
                                       / np.maximum(1e-12, np.abs(A["q1_scores"])))),
                    "Q3": float(np.max(np.abs(A["q3_scores"] - vk.field_gain(1, reg["anm"]["field_representation"]) * S3)
                                       / np.maximum(1e-12, np.abs(A["q3_scores"]))))},
                "max_abs_diff_engine_vs_recoded_closure": float(np.max(np.abs(A["cl_scores"] - clr))),
                "engine_events_rejected": int(A["n_rejected"].sum()),
                "registered_field_gain": G,
                "top_score_order_identical": bool(np.array_equal(
                    e1.ranked_order(S1.max(axis=1), P["ids"], reg["amendment_A1"]),
                    e1.ranked_order(A["q1_scores"].max(axis=1), P["ids"], reg["amendment_A1"]))),
                "closure_order_identical_engine_vs_recoded": bool(np.array_equal(
                    e1.ranked_order(clr.max(axis=1), P["ids"], reg["amendment_A1"]),
                    e1.ranked_order(A["cl_scores"].max(axis=1), P["ids"], reg["amendment_A1"]))),
                "tied_top_classes_share": {"Q1": float(np.mean(e1.tied_top_classes(S1))),
                                           "Q3": float(np.mean(e1.tied_top_classes(S3)))}}
        e11a["passed"] = bool(sum(e11a["mismatches"].values()) == 0 and sum(e11a["closure_mismatches"].values()) == 0
                              and e11a["engine_events_rejected"] == 0)
        R["E1.1a"] = e11a

        # ---------------- E1.1b: calls that change Q1 -> Q2 and Q1 -> Q3, per annotated type
        e11b = {"overall": {"Q1_to_Q2": float(np.mean(rq["Q1"] != rq["Q2"])), "Q1_to_Q3": float(np.mean(rq["Q1"] != rq["Q3"]))},
                "per_type": {}, "transitions_Q1_Q3": {}}
        for t in sorted(set(ct)):
            m = ct == t
            e11b["per_type"][t] = {"n": int(m.sum()), "Q1_to_Q2": float(np.mean(rq["Q1"][m] != rq["Q2"][m])),
                                   "Q1_to_Q3": float(np.mean(rq["Q1"][m] != rq["Q3"][m]))}
        for c1 in list(vk.LINEAGES) + [vk.NO_CALL]:
            e11b["transitions_Q1_Q3"][c1 or "no_call"] = {(c3 or "no_call"): int(np.sum((rq["Q1"] == c1) & (rq["Q3"] == c3)))
                                                         for c3 in list(vk.LINEAGES) + [vk.NO_CALL]}
        R["E1.1b"] = e11b

        # ---------------- E1.1c / E1.1d: Q3 at the rule's realised Q3 coverage
        e1cd = {"coverage_rule_Q3": value_block(st, lambda g: g("cov|c3"), ps, groups), "margin": E["exp1_change_the_question"]["margin"]}
        mg1 = float(E["exp1_change_the_question"]["margin"])
        for kn in K3:
            blk = {"selective_accuracy": {}}
            for arm in ("rule_q3", "clf_q1labels_q3features", "clf_q3labels", "clf_q1labels_anchors", "clf_q3labels_anchors"):
                blk["selective_accuracy"][arm] = value_block(st, lambda g, arm=arm, kn=kn: g(f"{arm}|c3|{kn}|acc"), ps, groups)
            for ep, arm in (("E1.1c", "clf_q1labels_q3features"), ("E1.1c_secondary_anchors", "clf_q1labels_anchors"),
                            ("E1.1d", "clf_q3labels"), ("E1.1d_secondary_anchors", "clf_q3labels_anchors")):
                blk[ep] = comparison(st, lambda g, arm=arm, kn=kn: np.asarray(g(f"rule_q3|c3|{kn}|acc")) - np.asarray(g(f"{arm}|c3|{kn}|acc")),
                                     mg1, ps, groups)
                blk[ep]["difference"] = f"rule minus {arm} (win = the rule is better)"
            e1cd[kn] = blk
        R["E1.1c_d"] = e1cd

        # ---------------- E1.3
        R["E1.3"] = {"Q1": e13_block(P, A, reg, ev, group, groups, cfg, "Q1"),
                     "Q2_secondary": e13_block(P, A, reg, ev, group, groups, cfg, "Q2")}
        mg3 = float(E["exp3_why_this_call"]["margin"])
        for q, tag in (("Q1", "Q1"), ("Q2", "Q2_secondary")):
            blk = {"n_flip_sensitive": value_block(st, lambda g, q=q: g(f"e13c|{q}|n_flip"), ps, groups)}
            for kn in K1:
                blk[kn] = {"err_flip": value_block(st, lambda g, q=q, kn=kn: g(f"e13c|{q}|{kn}|err_flip"), ps, groups),
                           "err_lowest_margin": value_block(st, lambda g, q=q, kn=kn: g(f"e13c|{q}|{kn}|err_low"), ps, groups),
                           "difference": comparison(st, lambda g, q=q, kn=kn: np.asarray(g(f"e13c|{q}|{kn}|err_flip"))
                                                    - np.asarray(g(f"e13c|{q}|{kn}|err_low")), mg3, ps, groups)}
            R["E1.3"][tag]["E1.3c"] = blk

        # ---------------- E1.4: which calls to trust (Q1 = Q2, computed once)
        mg4 = float(E["exp4_which_calls_to_trust"]["margin"])
        e14: dict = {"note": E["exp4_which_calls_to_trust"]["interpretation"], "margin": mg4, "keys": {}}
        for kn in K1:
            blk = {"AURC": {}, "acc@0.90": {}, "acc@0.70": {}}
            for arm in ("top", "margin", "entropy", "classifier", "classifier_rule_calls"):
                if kn != "primary" and arm == "classifier_rule_calls":
                    continue
                blk["AURC"][arm] = value_block(st, lambda g, arm=arm, kn=kn: aurc_of(g, arm, kn, gn, grid), ps, groups)
                blk["acc@0.90"][arm] = value_block(st, lambda g, arm=arm, kn=kn: g(f"{arm}|g0.90|{kn}|acc"), ps, groups)
                blk["acc@0.70"][arm] = value_block(st, lambda g, arm=arm, kn=kn: g(f"{arm}|g0.70|{kn}|acc"), ps, groups)
            for other in ("margin", "classifier", "entropy"):
                blk[f"top_vs_{other}"] = {
                    "E1.4a_AURC": comparison(st, lambda g, o=other, kn=kn: aurc_of(g, "top", kn, gn, grid) - aurc_of(g, o, kn, gn, grid), mg4, ps, groups),
                    "E1.4b_acc@0.90": comparison(st, lambda g, o=other, kn=kn: np.asarray(g(f"top|g0.90|{kn}|acc")) - np.asarray(g(f"{o}|g0.90|{kn}|acc")), mg4, ps, groups),
                    "E1.4b_acc@0.70": comparison(st, lambda g, o=other, kn=kn: np.asarray(g(f"top|g0.70|{kn}|acc")) - np.asarray(g(f"{o}|g0.70|{kn}|acc")), mg4, ps, groups)}
            if kn == "primary":
                nm = sum(1 for m in ev.methods if m.name.startswith("mlp_s"))
                blk["AURC"]["mlp_mean_over_seeds"] = value_block(
                    st, lambda g: np.mean([aurc_of(g, f"mlp_s{s}", "primary", gn, grid) for s in range(nm)], axis=0), ps, groups)
                for cc_ in ("0.90", "0.70"):
                    blk[f"acc@{cc_}"]["mlp_mean_over_seeds"] = value_block(
                        st, lambda g, cc_=cc_: np.mean([np.asarray(g(f"mlp_s{s}|g{cc_}|primary|acc")) for s in range(nm)], axis=0),
                        ps, groups)
                blk["acc_grid"] = {arm: {c: st.point.get(f"{arm}|{c}|primary|acc") for c in gn}
                                   for arm in ("top", "margin", "entropy", "classifier")}
            e14["keys"][kn] = blk
        sens = {}
        for pn in addendum["sensitivity_panels"]:
            sens[pn] = {"AURC": {arm: value_block(st, lambda g, arm=arm, pn=pn: aurc_of(g, f"{arm}@{pn}", "primary", gn, grid), ps, groups)
                                 for arm in ("top", "margin", "entropy", "classifier")},
                        "top_vs_margin_AURC": comparison(st, lambda g, pn=pn: aurc_of(g, f"top@{pn}", "primary", gn, grid)
                                                         - aurc_of(g, f"margin@{pn}", "primary", gn, grid), mg4, ps, groups)}
        e14["sensitivity_panels"] = sens
        e14["verdict_primary"] = e14["keys"]["primary"]["top_vs_margin"]["E1.4a_AURC"]["verdict"]
        R["E1.4"] = e14

        # ---------------- E1.5: label cost
        e15: dict = {"tolerance": 0.005}
        for cn, kn_list in (("q3", K3), ("anch", K3), ("q1ref", ("primary",))):
            lc = ev.lc[cn]
            ndraw = len(cfg["lc_draws"])
            rule_arm, cov = ("rule_q3", "c3") if cn != "q1ref" else ("top", "c1")
            cblk = {}
            for kn in kn_list:
                rows = {}

                def med(g, n, kn=kn, cn=cn, cov=cov):
                    return np.median(np.stack([np.asarray(g(f"lc_{cn}_n{n}_d{d}|{cov}|{kn}|acc"), dtype=np.float64)
                                               for d in range(ndraw)]), axis=0)

                for n in lc["ns"]:
                    rows[str(n)] = {"median_acc": value_block(st, lambda g, n=n: med(g, n), ps, groups),
                                    "gap_vs_rule": value_block(st, lambda g, n=n, kn=kn: med(g, n) - np.asarray(g(f"{rule_arm}|{cov}|{kn}|acc")), ps, groups),
                                    "per_draw_point": [st.point[f"lc_{cn}_n{n}_d{d}|{cov}|{kn}|acc"] for d in range(ndraw)]}

                def nstar(g, kn=kn):
                    rule = np.asarray(g(f"{rule_arm}|{cov}|{kn}|acc"), dtype=np.float64)
                    res_ = np.full(rule.shape, np.nan)
                    for i_n, n in enumerate(lc["ns"]):
                        ok = (med(g, n) >= rule - 0.005) & np.isnan(res_)
                        res_ = np.where(ok, i_n, res_)
                    return res_

                ns_point = nstar(lambda k: st.point[k])
                bdist = np.asarray(nstar(lambda k: st.boots[ps or f"group_{groups[0]}"][k]))
                cblk[kn] = {"rule_acc": value_block(st, lambda g, kn=kn: g(f"{rule_arm}|{cov}|{kn}|acc"), ps, groups),
                            "curve": rows,
                            "n_star": None if np.isnan(ns_point) else str(lc["ns"][int(ns_point)]),
                            "n_star_per_donor": {g: (None if np.isnan(v_) else str(lc["ns"][int(v_)]))
                                                 for g in groups for v_ in [nstar(lambda k, g=g: st.groups[g][k])]},
                            "n_star_bootstrap_distribution": {str(lc["ns"][i]): float(np.mean(bdist == i)) for i in range(len(lc["ns"]))}
                            | {"not_reached": float(np.mean(np.isnan(bdist)))}}
            e15[cn] = cblk
        e15["note"] = "a measurement, not a contest; the rule (and ANM, identical to it) needs 0 labels"
        R["E1.5"] = e15

        # ---------------- C2: nested readouts
        mg2 = float(E["c2_nested_readouts"]["margin"])
        c2: dict = {"c_star": value_block(st, lambda g: g("cov|c1"), ps, groups), "margin": mg2, "guard": -0.005, "keys": {}}
        for kn in K1:
            blk = {}
            for c in ("c1", "g0.80", "g0.70"):
                cb = {arm: {"out_decline": value_block(st, lambda g, arm=arm, c=c, kn=kn: g(f"{arm}|{c}|{kn}|outdecl"), ps, groups),
                            "inscope_acc": value_block(st, lambda g, arm=arm, c=c, kn=kn: g(f"{arm}|{c}|{kn}|inscope"), ps, groups)}
                      for arm in ("closure", "top", "classifier")}
                cb["closure_minus_rule_out_decline"] = comparison(
                    st, lambda g, c=c, kn=kn: np.asarray(g(f"closure|{c}|{kn}|outdecl")) - np.asarray(g(f"top|{c}|{kn}|outdecl")), mg2, ps, groups)
                cb["closure_minus_rule_inscope_acc"] = value_block(
                    st, lambda g, c=c, kn=kn: np.asarray(g(f"closure|{c}|{kn}|inscope")) - np.asarray(g(f"top|{c}|{kn}|inscope")), ps, groups)
                cb["classifier_minus_rule_out_decline_secondary"] = value_block(
                    st, lambda g, c=c, kn=kn: np.asarray(g(f"classifier|{c}|{kn}|outdecl")) - np.asarray(g(f"top|{c}|{kn}|outdecl")), ps, groups)
                guard_ok = cb["closure_minus_rule_inscope_acc"]["point"] >= -0.005
                cb["guard_inscope_ok"] = bool(guard_ok)
                vo = cb["closure_minus_rule_out_decline"]["verdict"]
                if vo == "win":
                    cb["verdict"] = "win" if guard_ok else "not a win (in-scope guard failed)"
                else:
                    cb["verdict"] = vo + ("" if guard_ok else " (in-scope guard also failed)")
                blk[{"c1": "c_star", "g0.80": "0.80", "g0.70": "0.70"}[c]] = cb
            c2["keys"][kn] = blk
        # nested pair at the deployed bars (rule and engine closure)
        nest = {}
        for arm, cq1, cq2 in (("rule", rq["Q1"], rq["Q2"]), ("anm_closure", cl_eng["Q1"], cl_eng["Q2"])):
            k = keys["primary"]
            nest[arm] = {"Q2_calls_nested_in_Q1_same_class": bool(np.all((cq2 == vk.NO_CALL) | (cq2 == cq1))),
                         **{q: {"P_f": float(np.mean(cc != vk.NO_CALL)),
                                "Q_f_among_scored": float(np.mean((cc == k)[k != vk.UNSCORED]))} for q, cc in (("Q1", cq1), ("Q2", cq2))}}
        c2["nested_pair_at_bars"] = nest
        c2["verdict_primary"] = c2["keys"]["primary"]["c_star"]["verdict"]
        R["C2"] = c2

        # ---------------- deployed operating points
        ops = {}
        for q in ("Q1", "Q2"):
            ops[f"rule_{q}"] = operating_point(rq[q], keys, K1, ct, group, groups)
            ops[f"anm_{q}"] = {"identical_to_rule": bool(np.array_equal(eng[q], rq[q]))}
            ops[f"anm_closure_{q}"] = operating_point(cl_eng[q], keys, K1, ct, group, groups)
            cc, ck = e1.lineage_conf_call(P["P5_primary"])
            ops[f"classifier_{q}"] = operating_point(np.where(cc >= reg["classifier"]["primary"]["bars"][q], ck, vk.NO_CALL),
                                                     keys, K1, ct, group, groups)
            for pn, ps_ in addendum["sensitivity_panels"].items():
                Sp = P[f"S_sens_{pn}"]
                ops[f"rule_{q}@{pn}"] = operating_point(vk.rule_calls(Sp, ps_["rule_bars"][q]), keys, ("primary",), ct, group, groups)
                pc, pk = e1.lineage_conf_call(P[f"P5_sens_{pn}"])
                ops[f"classifier_{q}@{pn}"] = operating_point(np.where(pc >= ps_["classifier"]["bars"][q], pk, vk.NO_CALL),
                                                              keys, ("primary",), ct, group, groups)
        ops["rule_Q3"] = operating_point(rq["Q3"], keys, K3, ct, group, groups)
        ops["anm_Q3"] = {"identical_to_rule": bool(np.array_equal(eng["Q3"], rq["Q3"]))}
        q3bars = addendum["q3_operating_point_bars"]["classifiers"]
        for nm, bn in (("q3labels", "registered_q3_classifier"), ("q1labels_q3features", "q1labels_q3features"),
                       ("q1labels_anchors", "q1labels_anchors"), ("q3labels_anchors", "q3labels_anchors")):
            pc, pk = e1.lineage_conf_call(P[f"P5_{nm}"])
            ops[f"classifier_{nm}_Q3"] = operating_point(np.where(pc >= q3bars[bn]["bar_Q3"], pk, vk.NO_CALL), keys, K3, ct, group, groups)
        R["operating_points"] = ops
        # internal consistency: matched-coverage selection at the rule's own realised coverage = its deployed calls
        cons = {"rule_Q1_at_c1_equals_operating_point": abs(st.point["top|c1|primary|acc"]
                                                             - ops["rule_Q1"]["primary"]["accuracy_of_calls"]) < 1e-12,
                "rule_Q3_at_c3_equals_operating_point": abs(st.point["rule_q3|c3|q3|acc"]
                                                             - ops["rule_Q3"]["q3"]["accuracy_of_calls"]) < 1e-12,
                "out_decline_at_c1_equals_operating_point": abs(st.point["top|c1|primary|outdecl"]
                                                                 - ops["rule_Q1"]["primary"]["out_decline_rate"]) < 1e-12}
        cons["passed"] = all(cons.values())
        R["consistency_checks"] = cons
        if not cons["passed"]:
            say(f"WARNING consistency checks failed on {sp['name']}: {cons}")
        res["splits"][role] = R
        say(f"report: {sp['name']} endpoints computed")

    # ---------------- summary of registered verdicts (primary split, primary key)
    Rp = res["splits"]["primary"]
    interp = Rp["E1.1a"]["passed"] and res["splits"].get("secondary", {}).get("E1.1a", {}).get("passed", True)
    v14 = Rp["E1.4"]["verdict_primary"]
    vc2 = Rp["C2"]["verdict_primary"]
    summ = {
        "E1.1a_bridge_equals_rule": Rp["E1.1a"]["passed"],
        "E1.1c_rule_vs_Q1label_classifier": Rp["E1.1c_d"]["q3"]["E1.1c"]["verdict"],
        "E1.1d_rule_vs_Q3label_classifier": Rp["E1.1c_d"]["q3"]["E1.1d"]["verdict"],
        "E1.3a_discrepancy": Rp["E1.3"]["Q1"]["E1.3a"]["discrepancy"],
        "E1.3c_flip_vs_lowest_margin": Rp["E1.3"]["Q1"]["E1.3c"]["primary"]["difference"]["verdict"],
        "E1.4a_top_score_vs_margin_AURC": v14,
        "E1.5_n_star_q3": Rp["E1.5"]["q3"]["q3"]["n_star"],
        "E1.C2_closure_vs_mean_rule_at_c_star": vc2,
        "key_validity_flag_primary": kv["primary"]["kappa_flag"]["below"],
    }
    if not interp:
        summ["falsification"] = ("not interpretable: E1.1a found a mismatch between ANM's engine and its declared rule "
                                 "(registered: stop and fix before any interpretation)")
    elif v14 != "win" and vc2 != "win":
        summ["falsification"] = ("falsified for v3 (A1 wording): neither E1.4 nor C2 is a win, so the readout forms ANM "
                                 "provides (top action score, closure readout) add no decision value over TEDDY's margin "
                                 "and the mean rule. Every ANM arm equals a re-coded rule cell by cell (E1.1a).")
    else:
        won = [n for n, v in (("E1.4 (top-score readout)", v14), ("C2 (closure readout)", vc2)) if v == "win"]
        summ["falsification"] = ("not falsified: " + ", ".join(won) + " is a win; credited to the readout form, which a "
                                 "written rule implements (A1), not to an ANM-specific mechanism")
    summ["E1.3c_reading"] = ("attribution flags errors beyond the margin" if summ["E1.3c_flip_vs_lowest_margin"] == "win"
                             else "attribution adds no flagging value beyond TEDDY's margin")
    res["summary"] = summ
    res["implementation_choices"] = json.loads(a.addendum.read_text())["implementation_choices"]
    return res


# ============================================================================ REPORT.md
def _f(x, n=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{n}f}"


def _ci(c, n=3):
    return f"[{_f(c[0], n)}, {_f(c[1], n)}]"


def _pd(d, n=3):
    return ", ".join(f"{g} {_f(v, n)}" for g, v in d.items())


def write_report(path: Path, res: dict) -> None:
    L = []
    title = "E1 (v3): Mode B rerun on the v3 key" + (" - SMOKE TEST, not a result" if res["smoke"] else "")
    L += [f"# {title}", "",
          f"Registration `{res['registration_sha256'][:16]}`, amendment A1 `{res['amendment_A1_sha256'][:16]}`, "
          f"E1 addendum `{(res['addendum_E1_sha256'] or 'none')[:16]}`; builder {res['builder']}; git `{(res['state'].get('git_head') or '')[:10]}`.", ""]
    if res["smoke"]:
        L += [f"**Smoke run:** {res['SMOKE']}.", ""]
    # key validity first
    L += ["## Key validity (reported before any method result)", "",
          "| split | cells | agreement | kappa (5 classes) | kappa < 0.85? | per-batch gate kappa |", "|---|---:|---:|---:|---|---:|"]
    for role, kv in res["key_validity"].items():
        L.append(f"| {kv['split']} | {kv['all']['n']:,} | {_f(kv['all']['agreement'], 4)} | {_f(kv['all']['kappa5'], 4)} | "
                 f"{'yes: annotation-only rows have equal prominence' if kv['kappa_flag']['below'] else 'no'} | "
                 f"{_f(kv['per_batch_gate_vs_annotation']['kappa5'], 4)} |")
        for g, r in kv["per_group"].items():
            L.append(f"| &nbsp; {g} | {r['n']:,} | {_f(r['agreement'], 4)} | {_f(r['kappa5'], 4)} | | |")
    L += ["", "Per-type gate tables, key counts and per-class precision/recall: `key_validity.json`.", ""]
    s = res["summary"]
    L += ["## Registered verdicts (primary split, primary key)", "", "| endpoint | result |", "|---|---|"]
    for k, v in s.items():
        L.append(f"| {k} | {v} |")
    L.append("")
    for role, R in res["splits"].items():
        L += [f"## {R['split']} ({R['n_cells']:,} cells; " + ", ".join(f"{g}: {n:,}" for g, n in R["groups"].items()) + ")", ""]
        if role == "secondary":
            L += ["Secondary split: reported with the same rules; never changes a primary verdict.", ""]
        a = R["E1.1a"]
        L += ["### E1.1a ANM engine vs re-coded rule", "",
              f"Call mismatches (registered 0): Q1 {a['mismatches']['Q1']}, Q2 {a['mismatches']['Q2']}, Q3 {a['mismatches']['Q3']}; "
              f"closure Q1 {a['closure_mismatches']['Q1']}, Q2 {a['closure_mismatches']['Q2']}. Rejected events {a['engine_events_rejected']}. "
              f"Max relative difference engine vs G(n)·n·S: Q1 {a['max_rel_diff_engine_vs_Gn_times_S']['Q1']:.2e}, "
              f"Q3 {a['max_rel_diff_engine_vs_Gn_times_S']['Q3']:.2e}; closure {a['max_abs_diff_engine_vs_recoded_closure']:.2e}. "
              f"Top-score order identical: {a['top_score_order_identical']}. **{'pass' if a['passed'] else 'FAIL: interpretation stops'}**. "
              f"Share of cells with tied top classes: Q1 {_f(a['tied_top_classes_share']['Q1'])}, Q3 {_f(a['tied_top_classes_share']['Q3'])}.", ""]
        b = R["E1.1b"]
        L += [f"### E1.1b calls that change (rule = ANM): Q1→Q2 {_f(b['overall']['Q1_to_Q2'])}, Q1→Q3 {_f(b['overall']['Q1_to_Q3'])} "
              "(per annotated type in the JSON)", ""]
        cd = R["E1.1c_d"]
        L += ["### E1.1c / E1.1d Q3 selective accuracy at the rule's realised Q3 coverage "
              f"({_f(cd['coverage_rule_Q3']['point'])}); margin {cd['margin']}", "",
              "| key | arm | selective accuracy [95% CI] | rule minus arm [95% CI] | per donor | verdict |", "|---|---|---|---|---|---|"]
        for kn in K3:
            blk = cd[kn]
            L.append(f"| {kn} | rule (Q3 anchors) | {_f(blk['selective_accuracy']['rule_q3']['point'])} {_ci(blk['selective_accuracy']['rule_q3']['ci95'])} | | | |")
            for ep, arm in (("E1.1c", "clf_q1labels_q3features"), ("E1.1c_secondary_anchors", "clf_q1labels_anchors"),
                            ("E1.1d", "clf_q3labels"), ("E1.1d_secondary_anchors", "clf_q3labels_anchors")):
                c = blk[ep]
                L.append(f"| {kn} | {ep}: {arm} | {_f(blk['selective_accuracy'][arm]['point'])} {_ci(blk['selective_accuracy'][arm]['ci95'])} | "
                         f"{_f(c['point'])} {_ci(c['ci95'])} | {_pd(c['per_donor'])} | {c['verdict']} |")
        L += ["", "'Zero labels to change the question' holds for any written rule (the rule and ANM alike); it is not an ANM property.", ""]
        L += ["### E1.3 why this call (Q1 primary, Q2 secondary)", ""]
        for tag in ("Q1", "Q2_secondary"):
            e = R["E1.3"][tag]
            a3 = e["E1.3a"]
            L.append(f"- **{tag}** ({e['n_calls']:,} calls). E1.3a: top-1 set agreement {_f(a3['top_set_agreement'], 4)}, flip agreement "
                     f"{_f(a3['flip_agreement'], 4)} (expected 1.000{'; DISCREPANCY' if a3['discrepancy'] else ''}); calls with a multi-marker "
                     f"set {_f(a3['share_calls_multi_marker_set'])}; flip-sensitive calls {a3['n_flip_closed_form']:,} "
                     f"(engine {a3['n_flip_engine']:,}).")
            for kn in ("primary", "annotation_only"):
                c = e["E1.3c"][kn]
                L.append(f"  E1.3c ({kn} key): error rate of flip-sensitive calls {_f(c['err_flip']['point'])} {_ci(c['err_flip']['ci95'])} "
                         f"vs as many lowest-margin calls {_f(c['err_lowest_margin']['point'])} {_ci(c['err_lowest_margin']['ci95'])}; "
                         f"difference {_f(c['difference']['point'])} {_ci(c['difference']['ci95'])}, per donor "
                         f"{_pd(c['difference']['per_donor'])}; margin {c['difference']['margin']}: **{c['difference']['verdict']}**.")
        L.append("")
        L += ["E1.3b marker shares (within the called class; fractional credit; null mean 1/3 by construction):", "",
              "| donor | marker | share within class | null 95% interval | p |", "|---|---|---:|---|---:|"]
        for g, blk in R["E1.3"]["Q1"]["E1.3b"]["per_donor"].items():
            if blk.get("all_calls") is None:
                continue
            for mk, r in blk["all_calls"]["markers"].items():
                L.append(f"| {g} | {mk} | {_f(r['share_within_class'])} | {_ci(r['null_ci95_within_class'])} | {_f(r['p_two_sided'], 4)} |")
        L.append("")
        e14 = R["E1.4"]
        L += [f"### E1.4 which calls to trust (Q1 = Q2; margin {e14['margin']})", "", f"{e14['note']}.", "",
              "| key | score | AURC [95% CI] | acc@0.90 | acc@0.70 |", "|---|---|---|---:|---:|"]
        for kn in K1:
            blk = e14["keys"][kn]
            for arm, r in blk["AURC"].items():
                a90 = blk["acc@0.90"].get(arm, {}).get("point")
                a70 = blk["acc@0.70"].get(arm, {}).get("point")
                L.append(f"| {kn} | {arm} | {_f(r['point'])} {_ci(r['ci95'])} | {_f(a90)} | {_f(a70)} |")
        L += ["", "| key | comparison | difference [95% CI] | per donor | verdict |", "|---|---|---|---|---|"]
        for kn in K1:
            for other in ("margin", "classifier", "entropy"):
                for ep, c in e14["keys"][kn][f"top_vs_{other}"].items():
                    if kn != "primary" and ep != "E1.4a_AURC":
                        continue
                    L.append(f"| {kn} | top score vs {other}, {ep} | {_f(c['point'])} {_ci(c['ci95'])} | {_pd(c['per_donor'])} | {c['verdict']} |")
        L += ["", "Sensitivity panels (primary key; never decide a verdict):", "",
              "| panel | AURC top | AURC margin | AURC classifier | top minus margin [95% CI] | verdict |", "|---|---:|---:|---:|---|---|"]
        for pn, blk in e14["sensitivity_panels"].items():
            c = blk["top_vs_margin_AURC"]
            L.append(f"| {pn} | {_f(blk['AURC']['top']['point'])} | {_f(blk['AURC']['margin']['point'])} | "
                     f"{_f(blk['AURC']['classifier']['point'])} | {_f(c['point'])} {_ci(c['ci95'])} | {c['verdict']} |")
        L.append("")
        e15 = R["E1.5"]
        L += ["### E1.5 label cost (a measurement, not a contest)", ""]
        for cn, lab in (("q3", "Q3 classifier (14 features)"), ("anch", "4 anchors (secondary)"), ("q1ref", "Q1 reference")):
            for kn, blk in e15[cn].items():
                L.append(f"- **{lab}, key {kn}**: rule {_f(blk['rule_acc']['point'])}; n* = **{blk['n_star'] or 'not reached'}** "
                         f"(per donor {blk['n_star_per_donor']}; bootstrap share of n*: "
                         + ", ".join(f"{k} {_f(v, 2)}" for k, v in blk["n_star_bootstrap_distribution"].items() if v > 0) + ")")
                L.append("  " + "; ".join(f"n={n}: {_f(r['median_acc']['point'])} (gap {_f(r['gap_vs_rule']['point'])})" for n, r in blk["curve"].items()))
        L.append("")
        c2 = R["C2"]
        L += [f"### C2 nested readouts (c* = {_f(c2['c_star']['point'])}; margin {c2['margin']}; in-scope guard ≥ -0.005)", "",
              "| key | coverage | OUT decline: closure / mean rule / classifier | closure minus rule [95% CI] | per donor | in-scope acc diff | verdict |",
              "|---|---|---|---|---|---|---|"]
        for kn in K1:
            for cn_, cb in c2["keys"][kn].items():
                d = cb["closure_minus_rule_out_decline"]
                L.append(f"| {kn} | {cn_} | {_f(cb['closure']['out_decline']['point'])} / {_f(cb['top']['out_decline']['point'])} / "
                         f"{_f(cb['classifier']['out_decline']['point'])} | {_f(d['point'])} {_ci(d['ci95'])} | {_pd(d['per_donor'])} | "
                         f"{_f(cb['closure_minus_rule_inscope_acc']['point'])} | {cb['verdict']} |")
        n_ = c2["nested_pair_at_bars"]
        L += ["", "Nested pair at the deployed bars: " + "; ".join(
            f"{arm}: Q2 nested in Q1 {r['Q2_calls_nested_in_Q1_same_class']}, P_f Q1 {_f(r['Q1']['P_f'])}, Q2 {_f(r['Q2']['P_f'])}, "
            f"Q_f Q1 {_f(r['Q1']['Q_f_among_scored'])}, Q2 {_f(r['Q2']['Q_f_among_scored'])}" for arm, r in n_.items()), ""]
        ops = R["operating_points"]
        L += ["### Deployed operating points (bars from val)", "",
              "| method | coverage | selective acc (key) | decision acc | OUT decline | calls on unscored |", "|---|---:|---:|---:|---:|---:|"]
        for nm, op in ops.items():
            if "coverage" not in op or "@" in nm:
                continue
            kn = "q3" if nm.endswith("Q3") else "primary"
            r = op[kn]
            L.append(f"| {nm} | {_f(op['coverage'])} | {_f(r['accuracy_of_calls'])} | {_f(r['decision_accuracy'])} | "
                     f"{_f(r['out_decline_rate'])} | {sum(v for k, v in r['calls_on_unscored'].items() if k != 'n_unscored')} |")
        L.append("")
    L += ["## Implementation choices (declared in registration/addenda/E1.json before site4)", ""]
    L += [f"- {c}" for c in res["implementation_choices"]]
    L += ["", "Full numbers: `E1_results.json`; key validity: `key_validity.json`; log: `progress.log`.", ""]
    path.write_text("\n".join(L))


# ============================================================================ main
def main(argv=None) -> int:
    global _LOG
    a = parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG = a.out_dir / "progress.log"
    import torch

    torch.set_num_threads(int(a.threads))
    if a.anm_root is None and os.environ.get("ANM_ROOT"):
        a.anm_root = Path(os.environ["ANM_ROOT"])
    reg = va.load_registration_amended(a.registration, a.amendment)
    cfg = run_config(a, reg)
    st = registration_state(a)
    say(f"E1 {VERSION} stage={a.stage} smoke={a.smoke} out={a.out_dir}")
    if not a.smoke:
        guard_site4(a, st)
    if a.stage in ("all", "prepare"):
        st["addendum_reproduces"] = check_addendum_reproduces(a, st)
        say(f"addendum core reproduces from train/val: {st['addendum_reproduces']['matches']}")
    if not a.addendum.exists():
        raise SystemExit("E1 addendum missing (bridge_anm/v3_e1_addendum.py)")
    addendum = json.loads(a.addendum.read_text())["computed"]
    fp = file_fingerprint(a, st, cfg)
    meta = load_meta(a)
    splits = eval_splits(a, reg, meta)
    write_json(a.out_dir / "run_state.json", {**st, "fingerprint": fp, "config": cfg,
                                              "splits": {r: {"name": s["name"], "n": int(s["idx"].size),
                                                             "groups": s["groups"]} for r, s in splits.items()}})
    t0 = time.time()
    if a.stage in ("all", "prepare"):
        stage_prepare(a, reg, cfg, fp, splits, meta)
    if a.stage in ("all", "anm"):
        stage_anm(a, reg, cfg, fp, splits)
    if a.stage in ("all", "boot"):
        stage_boot(a, reg, cfg, fp, splits, addendum)
    if a.stage in ("all", "report"):
        res = stage_report(a, reg, cfg, fp, splits, addendum, st)
        res["seconds_this_invocation"] = time.time() - t0
        write_json(a.out_dir / "E1_results.json", res)
        write_report(a.out_dir / "REPORT.md", res)
        say(f"wrote {a.out_dir / 'E1_results.json'} and REPORT.md")
        for k, v in res["summary"].items():
            say(f"  {k}: {v}")
    say(f"done in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
