#!/usr/bin/env python3
"""E6 (v3): external confirmation of E1's conclusions on Hao et al. 2021 PBMC CITE-seq (GEO GSE164378, 3').

Every E6 output is labelled 'RNA possibly seen by TEDDY in pretraining' (registration experiments.E6.contamination;
evidence in registration/addenda/E6.json).

Registration: registration/registration_v3.json (experiments.E6) with amendments A1-A3
(bridge_anm/lib/v3_amend.load_registration_amended()) and the E6 addendum registration/addenda/E6.json (hash in
registration/addenda/HASHES.txt). Everything is frozen: TEDDY and the phase-1 head (L2-normalised official z), the
evidence normaliser (registered training q95, clip [0, 1]), bars, the registered classifier (refit exactly as E1 on
BMMC training cells), the ANM engine (finite_graph_scalar via ANM_ROOT, bridge_anm/lib/v3_e1.AnmBridge), the
panels and the gate logic; the addendum fixes the external map, the missing-protein reduction of the gate, the
per-donor gate estimator, the endpoints, the statistics and the replication rule.

Panels (addendum version 2): PRIMARY = the full frozen registered panels for every arm, identical to E1 (CD5 and
CD94 are not measured externally, but their evidence is the head's prediction from RNA and exists for every cell),
so every primary arm reads the same 12 evidence values and ANM = the mean rule. SENSITIVITY = the literal
reduced-panel reading (T on CD3, CD2; NK on CD122, CD56; arms suffixed _red, plus classifier10), reported beside
the primary and never deciding.

Endpoints (per external donor and pooled; two-stage donor -> cell bootstrap over the 8 donors, B = 2000):
  E1.1a  as E1: ANM engine vs the mean rule (Q1, Q2, Q3) and engine closure vs re-coded closure, registered 0;
         the bridge check (engine vs its exact re-coded closed form) kept beside it
  Q1 selective accuracy at matched coverage (grid and c* = the rule's realised Q1 coverage); E1.4b-style differences
  E1.4a  AURC: ANM top score vs margin (primary); vs classifier / entropy, rule vs margin (secondary)
  E1.C2  key-OUT decline at c*, 0.80, 0.70: ANM closure vs mean rule, in-scope guard
  replication of each v3 (E1 site4) conclusion: sign in a majority of donors and the pooled interval excluding 0
  sensitivity: the same rows on the reduced panels (own c*), with classifier10; never decides
Key validity is computed and written first; if annotation vs per-donor gate kappa < 0.70 the results are reported
on the annotation-only key and labelled 'key not validated'.

Stages (resumable; caches carry a fingerprint of the registration, amendments, addendum, code and inputs):
  prepare  keys and key validity (written first), evidence, scores, classifiers (BMMC training cells only)
  anm      ANM engine per cell (registered panels Q1/Q2 + closure; reduced panels; Q3 anchors)
  boot     pooled two-stage bootstrap replicates in chunks
  report   E6_results.json + REPORT.md
  all      prepare, anm, boot, report

The real run refuses unless the registration, A1-A3 and the E6 addendum are committed with matching hashes, the
external pack and its manifest match the addendum, and the official external embedding is complete (every shard and
the assembled z_rna.npy with a matching z_rna_manifest.json). CPU only.

--smoke runs the same code path on val donor 18303 cells (BMMC; 8 label-free pseudo-donors; the externally absent
proteins blanked; small B). It never opens the external pack or embedding.

Full run (repo root), once the official external embedding is complete (all 8 shards + z_rna.npy):
  ANM_ROOT=<anm_v2_fix> OMP_NUM_THREADS=4 <venv>/bin/python scripts/v3_e6_external.py --stage all --threads 4
Expected about 2-5 min on 4 CPU threads (prepare ~1 min, ANM engine ~0.5 min for 3 x 40,000 instances, bootstrap
~0.5 min for 2,000 two-stage replicates of 40,000 cells at ~8 ms each, report < 1 min); --stage report alone
re-reads the caches (e.g. after E1's site4 results change).
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
sys.path.insert(0, str(ROOT / "bridge_anm"))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e1 as e1  # noqa: E402
from lib import v3_e6 as e6  # noqa: E402
from lib import v3_key as vk  # noqa: E402

VERSION = "v3_e6_external 2.0"
KEYS = ("primary", "annotation_only")
LIN = np.asarray(vk.LINEAGES)
SCHEMA = ROOT / "bridge_anm" / "schemas" / "cite_lineage_finite_field_v0.json"
OUT_BASE = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3")
DEFAULT_OUT = OUT_BASE / "E6"
SMOKE_OUT = OUT_BASE / "E6_smoke"
E1_RESULTS = OUT_BASE / "E1" / "E1_results.json"
E1_SMOKE_RESULTS = OUT_BASE / "E1_smoke" / "E1_results.json"
SMOKE = {"n_pseudo_donors": 8, "seed": 0, "n_boot": 60}
C2_COVS = ("c1", "g0.80", "g0.70")          # primary: c1 = the mean rule's realised Q1 coverage (registered panels)
C2_COVS_RED = ("c1r", "g0.80", "g0.70")     # sensitivity: c1r = the reduced rule's realised Q1 coverage
SENS = "sensitivity: reduced panels (literal missing-protein reading); never decides"

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


# ============================================================================ args, guard
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", default="all", choices=["all", "prepare", "anm", "boot", "report"])
    p.add_argument("--registration", type=Path, default=ROOT / "registration/registration_v3.json")
    p.add_argument("--addendum", type=Path, default=ROOT / "registration/addenda/E6.json")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite", help="BMMC pack (classifier training)")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy", help="BMMC official z")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--external", type=Path, default=None, help="default: the addendum's pack directory")
    p.add_argument("--external-embedding", type=Path, default=None, help="default: the addendum's embedding directory")
    p.add_argument("--e1-results", type=Path, default=None, help="E1 site4 results (v3 direction for replication)")
    p.add_argument("--anm-root", type=Path, default=None, help="default: $ANM_ROOT")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--smoke", action="store_true", help="val donor 18303 cells through the E6 code path, small B")
    p.add_argument("--n-boot", type=int, default=None, help="override B (smoke / tests only; the real run uses 2000)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--chunk", type=int, default=100, help="bootstrap replicates per checkpoint file")
    a = p.parse_args(argv)
    a.out_dir = a.out_dir or (SMOKE_OUT if a.smoke else DEFAULT_OUT)
    a.e1_results = a.e1_results or (E1_SMOKE_RESULTS if a.smoke else E1_RESULTS)
    if a.n_boot is not None and not a.smoke:
        raise SystemExit("--n-boot is for smoke runs only (registered B = 2000)")
    return a


def _git(*args) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    except OSError:
        return None


def committed(path: Path) -> bool:
    r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
    r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
    return r1.returncode == 0 and r2.returncode == 0


CODE = [Path(__file__).resolve(), ROOT / "bridge_anm/lib/v3_e6.py", ROOT / "bridge_anm/lib/v3_e1.py",
        ROOT / "bridge_anm/lib/v3_key.py", ROOT / "bridge_anm/lib/v3_amend.py"]


def registration_state(a) -> dict:
    """Hashes and commit state of everything E6 reads before any external value."""
    rd = a.registration.parent
    st: dict = {"registration_sha256": vk.sha256_file(a.registration),
                "amendments_sha256": {f: vk.sha256_file(rd / f) for f in ("amendment_A1.json", "amendment_A2.json", "amendment_A3.json")}}
    files = [a.registration, a.registration.with_name(a.registration.name + ".sha256")]
    for f in st["amendments_sha256"]:
        files += [rd / f, rd / (f + ".sha256")]
    st["registration_and_amendments_committed"] = all(committed(f) for f in files)
    if a.addendum.exists():
        sha = vk.sha256_file(a.addendum)
        hashes = a.addendum.parent / "HASHES.txt"
        doc = json.loads(a.addendum.read_text())
        st.update({"addendum_sha256": sha,
                   "addendum_hash_recorded": hashes.exists() and f"{sha}  {a.addendum.name}" in hashes.read_text().splitlines(),
                   "addendum_committed": committed(a.addendum) and committed(hashes),
                   "addendum_names_these_files": doc.get("registration_sha256") == st["registration_sha256"]
                   and doc.get("amendments_sha256") == st["amendments_sha256"]})
    else:
        st.update({"addendum_sha256": None, "addendum_hash_recorded": False, "addendum_committed": False,
                   "addendum_names_these_files": False})
    st["code_sha256"] = {str(c.relative_to(ROOT)): vk.sha256_file(c) for c in CODE}
    st["code_committed_clean"] = all(committed(c) for c in CODE)
    st["git_head"] = _git("rev-parse", "HEAD")
    return st


def guard_registration(st: dict, smoke: bool) -> None:
    need = {"addendum E6.json present and recorded in addenda/HASHES.txt": st["addendum_hash_recorded"],
            "addendum names the registration and A1-A3 on disk": st["addendum_names_these_files"]}
    if not smoke:
        need.update({"registration and A1-A3 committed": st["registration_and_amendments_committed"],
                     "addendum E6.json and HASHES.txt committed": st["addendum_committed"],
                     "E6 runner and libraries committed and unmodified": st["code_committed_clean"]})
    bad = [k for k, ok in need.items() if not ok]
    if bad:
        raise SystemExit("E6 refused (commit the E6 addendum before any external outcome): " + "; ".join(bad))


FROZEN_INPUTS = {"ckpt_sha256": "ckpt", "bmmc_cite_arrays_sha256": "processed", "bmmc_z_sha256": "z"}


def guard_frozen_inputs(a, add: dict) -> dict:
    """The phase-1 head checkpoint, the BMMC pack and the BMMC official z must be the files the addendum was built
    on (addendum provenance.inputs); their sha256 are recorded in run_state.json and E6_results.json."""
    want = add.get("provenance", {}).get("inputs", {})
    paths = {"ckpt_sha256": a.ckpt, "bmmc_cite_arrays_sha256": a.processed / "cite_arrays.npz", "bmmc_z_sha256": a.z}
    got = {k: (vk.sha256_file(p) if p.exists() else None) for k, p in paths.items()}
    bad = [f"{k} {got[k]} != addendum {want.get(k)} ({paths[k]})" for k in FROZEN_INPUTS if got[k] != want.get(k)]
    if bad:
        raise SystemExit("E6 refused (frozen inputs differ from the addendum's provenance): " + "; ".join(bad))
    return got


def guard_external(a, add: dict) -> dict:
    """Pack, manifest, map and the official embedding must match the addendum (no value is read here)."""
    ed = add["external_data"]
    pack = a.external / "cite_arrays.npz"
    problems, info = [], {}
    info["pack_sha256"] = vk.sha256_file(pack) if pack.exists() else None
    if info["pack_sha256"] != ed["pack_sha256"]:
        problems.append(f"pack sha256 {info['pack_sha256']} != addendum {ed['pack_sha256']}")
    man_p = a.external / "external_manifest.json"
    if not man_p.exists() or vk.sha256_file(man_p) != ed["manifest_sha256"]:
        problems.append("external_manifest.json missing or its sha256 differs from the addendum")
    else:
        man = json.loads(man_p.read_text())
        m = man["annotation_map_addendum"]["map"]
        if m != add["annotation_map"]["map"] or e6.map_sha256(m) != add["annotation_map"]["sha256_canonical_json"]:
            problems.append("the manifest's annotation map differs from the addendum's")
    emb = a.external_embedding
    req = ed["embedding"]["manifest_must_match"]
    zf, mf = emb / "z_rna.npy", emb / "z_rna_manifest.json"
    shards = sorted((emb / "shards").glob("shard_*.npz")) if (emb / "shards").exists() else []
    n_shards = int(req["n_shards"])
    want = [emb / "shards" / f"shard_{k:04d}.npz" for k in range(n_shards)]
    info["shards_present"] = [s.name for s in shards]
    if [s for s in want if not s.exists()]:
        problems.append(f"official embedding incomplete: {len([s for s in want if s.exists()])} of {n_shards} shards")
    if not (zf.exists() and mf.exists()):
        problems.append("assembled z_rna.npy / z_rna_manifest.json missing (03_embed_rna.py writes them after the last shard)")
    if not problems:
        zm = json.loads(mf.read_text())
        bm = json.loads((a.z.parent / "z_rna_manifest.json").read_text())
        for k in ("preprocessing", "seq_len", "pooling", "normalize_total", "autocast", "length_buckets", "cells", "n_cells", "n_shards"):
            if zm.get(k) != req[k]:
                problems.append(f"embedding manifest {k} = {zm.get(k)!r}, addendum requires {req[k]!r}")
        if zm.get("medians_sha256") != bm.get("medians_sha256"):
            problems.append("embedding gene-median sha256 differs from the BMMC official embedding's")
        for k in ("preprocessing", "seq_len", "pooling", "normalize_total", "autocast", "length_buckets", "batch_size"):
            if zm.get(k) != bm.get(k):
                problems.append(f"embedding recipe {k} differs from the BMMC official embedding ({zm.get(k)!r} vs {bm.get(k)!r})")
        if Path(zm.get("processed_input", "")).resolve() != Path(req["processed_input"]).resolve():
            problems.append(f"embedding processed_input {zm.get('processed_input')} != {req['processed_input']}")
        for s in want:
            with np.load(s) as sh:
                if str(sh["cfg_hash"]) != zm.get("config_hash"):
                    problems.append(f"{s.name} cfg_hash differs from the manifest config_hash")
        if (emb / "z_rna_cells.npy").exists():
            problems.append("z_rna_cells.npy present: the embedding is not of all pack cells in order")
        info["embedding_manifest"] = {k: zm.get(k) for k in ("config_hash", "n_cells", "n_shards", "git_commit", "device",
                                                            "embed_runtime_sec", "medians_sha256")}
        info["z_sha256"] = vk.sha256_file(zf)
    if problems:
        raise SystemExit("E6 refused (external inputs not ready or not as registered): " + "; ".join(problems))
    return info


def fingerprint(a, st: dict, cfg: dict, ext_info: dict) -> str:
    h = hashlib.sha256()
    for k in ("registration_sha256", "addendum_sha256"):
        h.update(str(st.get(k)).encode())
    h.update(json.dumps(st["amendments_sha256"], sort_keys=True).encode())
    for k, v in sorted(st["code_sha256"].items()):
        h.update(f"{k}={v}".encode())
    h.update(json.dumps(st.get("frozen_inputs_sha256"), sort_keys=True).encode())
    for p in (a.processed / "cite_arrays.npz", a.z, a.ckpt):
        s = p.stat()
        h.update(f"{p}:{s.st_size}:{int(s.st_mtime)}".encode())
    h.update(json.dumps(ext_info, sort_keys=True, default=str).encode())
    h.update(json.dumps(cfg, sort_keys=True).encode())
    return h.hexdigest()[:20]


def run_config(a, reg) -> dict:
    B = int(a.n_boot) if a.n_boot is not None else (SMOKE["n_boot"] if a.smoke else 2000)
    return {"smoke": bool(a.smoke), "n_boot": B, "boot_seed": int(reg["seeds"]["bootstrap"]),
            "gate_seed": int(reg["gate"]["threshold_seed"]), **({"smoke_cells": SMOKE} if a.smoke else {})}


# ============================================================================ inputs
def load_inputs(a, reg, add: dict) -> dict:
    """The evaluated cells as one dict: ids (A1.1 tie-break ids), donors, labels + map, measured ADT (NaN = absent),
    the measured-protein list, and the official z. Smoke: BMMC val cells; real: the external pack + embedding."""
    absent_reg = list(add["external_data"]["adt_absent_externally"])
    if a.smoke:
        npz = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
        split, sites, donors = npz["split"].astype(str), npz["sites"].astype(str), npz["donors"].astype(str)
        idx = vk.split_indices(split, sites, donors, reg, "val")
        names = [str(x) for x in npz["adt_names"]]
        absent = absent_reg                       # blank exactly the proteins the external panel lacks
        present = [p for p in names if p not in absent]
        adt = np.asarray(npz["adt"], dtype=np.float32)[idx].copy()
        adt[:, [names.index(p) for p in absent]] = np.nan
        rng = np.random.default_rng(SMOKE["seed"])
        pseudo = np.asarray([f"valP{i + 1}" for i in rng.integers(0, SMOKE["n_pseudo_donors"], idx.size)])
        z = np.asarray(np.load(a.z, mmap_mode="r")[idx], dtype=np.float32)
        return {"name": "smoke_val_18303", "ids": idx.astype(np.int64), "donors": pseudo, "cell_types": npz["cell_types"].astype(str)[idx],
                "amap": dict(reg["annotation_map"]), "adt": adt, "adt_names": names, "present": present, "z": z,
                "absent": absent}
    pack = np.load(a.external / "cite_arrays.npz", allow_pickle=False)
    names = [str(x) for x in pack["adt_names"]]
    pres_mask = np.asarray(pack["adt_present"], dtype=bool)
    present = [n for n, m in zip(names, pres_mask) if m]
    absent = [n for n, m in zip(names, pres_mask) if not m]
    if sorted(absent) != sorted(absent_reg):
        raise SystemExit(f"pack adt_present disagrees with the addendum's absent proteins {absent_reg}")
    if names != list(reg["evidence"]["teddy_head"]["adt_names"]):
        raise SystemExit("external adt_names differ from the registered 134-protein order")
    n = int(pack["donors"].shape[0])
    z = np.load(a.external_embedding / "z_rna.npy", mmap_mode="r")
    if z.shape[0] != n:
        raise SystemExit(f"external z has {z.shape[0]} rows, pack has {n} cells")
    return {"name": "external_hao2021_d5k", "ids": np.arange(n, dtype=np.int64), "donors": pack["donors"].astype(str),
            "cell_types": pack["cell_types"].astype(str), "amap": dict(add["annotation_map"]["map"]),
            "adt": np.asarray(pack["adt"], dtype=np.float32), "adt_names": names, "present": present,
            "z": np.asarray(z, dtype=np.float32), "absent": absent}


# ============================================================================ stage: prepare
def per_type_table(ct: np.ndarray, annot: np.ndarray, gated: np.ndarray, prim: np.ndarray) -> dict:
    out = {}
    for t in sorted(set(ct.tolist())):
        m = ct == t
        out[t] = {"class": str(annot[m][0]), "n": int(m.sum()),
                  "gated": {c: int(np.sum(gated[m] == c)) for c in vk.CLASSES},
                  "in_primary_key": int(np.sum(prim[m] != vk.UNSCORED))}
    return out


def stage_prepare(a, reg, add, cfg, fp, D) -> None:
    cache = a.out_dir / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / "prepare.npz"
    if out.exists() and str(np.load(out)["fingerprint"]) == fp and (cache / "prepare_meta.json").exists():
        say("prepare: cached")
        return
    t0 = time.time()
    gate = reg["gate"]
    spec = e6.reduced_gate_spec(gate, D["present"])
    if not a.smoke and spec["clauses"] != add["key"]["clauses_remaining_in_order"]:
        raise SystemExit("reduced gate recomputed from the pack's adt_present differs from the addendum's")
    # ---- keys and key validity first (verifier side; before any method result)
    annot = e6.annotation_class(D["cell_types"], D["amap"])
    thr = e6.per_donor_thresholds(D["adt"], D["adt_names"], D["donors"], spec["proteins_used"], cfg["gate_seed"])
    nan_used = {p: int(np.isnan(D["adt"][:, D["adt_names"].index(p)]).sum()) for p in spec["proteins_used"]}
    gated = e6.gate_per_donor(D["adt"], D["adt_names"], D["donors"], spec, thr)
    keys = e6.e6_keys(annot, gated)
    kv = e6.key_validation(annot, gated, D["donors"])
    kv.update({"cells": D["name"], "n_cells": int(annot.size), "gate_spec": spec, "per_donor_thresholds": thr,
               "nan_values_in_used_gate_proteins": nan_used,
               "per_type": per_type_table(D["cell_types"], annot, gated, keys["primary"]),
               "key_counts": {k: {c: int(np.sum(keys[k] == c)) for c in list(vk.CLASSES) + [vk.UNSCORED]} for k in KEYS},
               "per_donor_key_counts": {d: {c: int(np.sum(keys["primary"][D["donors"] == d] == c))
                                            for c in list(vk.CLASSES) + [vk.UNSCORED]} for d in sorted(set(D["donors"].tolist()))}})
    write_json(a.out_dir / "key_validity.json", {"registration_sha256": vk.sha256_file(a.registration), "smoke": cfg["smoke"],
                                                 "contamination_label": e6.CONTAMINATION_LABEL, **kv})
    say(f"key validity ({D['name']}): pooled kappa annotation vs per-donor gate {kv['kappa5']:.4f} -> "
        f"{'validated' if kv['validated'] else e6.KEY_NOT_VALIDATED}; decision key {kv['decision_key']}")

    # ---- evidence (same for every method)
    names = list(reg["evidence"]["teddy_head"]["adt_names"])
    if D["adt_names"] != names:
        raise SystemExit("ADT name order differs from the registered head output order")
    v = vk.evidence(vk.head_predict(D["z"], a.ckpt, device="cpu", size_factor=1.0), reg)
    pan = vk.question_panel(reg, "Q1")                    # primary: the full frozen registered panels (as E1)
    pan_red = e6.reduced_panel(pan, D["present"])         # sensitivity: absent panel proteins dropped
    if not a.smoke and (pan != add["panels"]["primary"] or pan_red != add["panels"]["sensitivity_reduced"]):
        raise SystemExit("primary or reduced panel differs from the addendum's")
    p3 = vk.question_panel(reg, "Q3")
    fld, rc = reg["anm"]["field_representation"], reg["anm"]["closure_readout"]
    arr = {"fingerprint": np.asarray(fp), "ids": D["ids"], "donor": D["donors"].astype(str), "ct": D["cell_types"].astype(str),
           "v": v,
           "S": e6.class_scores_panel(v, names, pan), "S_red": e6.class_scores_panel(v, names, pan_red),
           "S3": e6.class_scores_panel(v, names, p3),
           "rec_q1": e6.anm_recoded_scores(v, names, pan, fld), "rec_cl": e6.anm_recoded_scores(v, names, pan, fld, rc),
           "rec_q1_red": e6.anm_recoded_scores(v, names, pan_red, fld),
           "rec_cl_red": e6.anm_recoded_scores(v, names, pan_red, fld, rc),
           "rec_q3": e6.anm_recoded_scores(v, names, p3, fld)}
    for k in ("annotation", "gated_per_donor", "primary", "annotation_only"):
        arr[f"key_{k}"] = keys[k].astype(str)

    # ---- classifiers: BMMC training cells only (as E1); frozen registered 12-feature model + classifier10
    npz = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
    split, sites, donors = npz["split"].astype(str), npz["sites"].astype(str), npz["donors"].astype(str)
    tr = vk.split_indices(split, sites, donors, reg, "train")
    vl = vk.split_indices(split, sites, donors, reg, "val")
    bnames = [str(x) for x in npz["adt_names"]]
    adt_b = np.asarray(npz["adt"], dtype=np.float32)
    ct_b = npz["cell_types"].astype(str)
    zmm = np.load(a.z, mmap_mode="r")
    j = {n: i for i, n in enumerate(names)}

    def bmmc(idx):
        k = vk.build_keys(ct_b[idx], adt_b[idx], bnames, reg)["primary"]
        ev = vk.evidence(vk.head_predict(np.asarray(zmm[idx], dtype=np.float32), a.ckpt, device="cpu", size_factor=1.0), reg)
        return ev, k

    v_tr, y_tr = bmmc(tr)
    v_vl, y_vl = bmmc(vl)
    mt, mv = y_tr != vk.UNSCORED, y_vl != vk.UNSCORED
    cp = reg["classifier"]["primary"]
    f12 = list(cp["feature_proteins"])
    if f12 != [p for k in vk.LINEAGES for p in pan[k]]:
        raise SystemExit("the registered classifier's features are not the primary panel's 12 proteins")
    c10 = add["computed"]["classifier10"]
    f10 = list(c10["features"])
    if f10 != [p for k in vk.LINEAGES for p in pan_red[k]]:
        raise SystemExit("classifier10 features differ from the reduced panel")
    from sklearn.metrics import log_loss

    meta: dict = {"fingerprint": fp, "seconds": None, "panel_primary": pan, "panel_reduced": pan_red, "q3_anchors": p3,
                  "evidence_proteins_primary": f12}
    for nm, feats, C in (("classifier", f12, float(cp["C"])), ("classifier10", f10, float(c10["C"]))):
        mdl = e1.fit_logreg(v_tr[mt][:, [j[p] for p in feats]], y_tr[mt], C)
        P = e1.proba5(mdl, v_vl[mv][:, [j[p] for p in feats]])[:, [vk.CLASSES.index(str(c)) for c in mdl.classes_]]
        ll = round(float(log_loss(y_vl[mv], P, labels=mdl.classes_)), 6)
        want = ({r["C"]: r["val_log_loss"] for r in cp["val_log_loss"]}[cp["C"]] if nm == "classifier"
                else {r["C"]: r["val_log_loss"] for r in c10["val_log_loss"]}[c10["C"]])
        if abs(ll - float(want)) > 1.5e-6:
            raise SystemExit(f"{nm}: val log-loss {ll} does not reproduce the registered / addendum value {want}")
        if nm == "classifier" and int(mt.sum()) != int(cp["n_train_labelled"]):
            raise SystemExit("classifier training set differs from the registration")
        arr[f"P5_{nm}"] = e1.proba5(mdl, v[:, [j[p] for p in feats]])
        meta[nm] = {"features": feats, "C": C, "n_train_labelled": int(mt.sum()), "val_log_loss": ll,
                    "classes": [str(c) for c in mdl.classes_]}
        say(f"{nm}: {len(feats)} features, C {C}, BMMC training cells {int(mt.sum())}, val log-loss {ll} (reproduced)")
    atomic_savez(out, **arr)
    meta["seconds"] = time.time() - t0
    write_json(cache / "prepare_meta.json", meta)
    say(f"prepare done ({annot.size} cells) in {time.time() - t0:.0f}s")


# ============================================================================ stage: anm
def stage_anm(a, reg, fp, chunk_cells: int = 2000) -> None:
    cache = a.out_dir / "cache"
    out = cache / "anm.npz"
    if out.exists() and str(np.load(out)["fingerprint"]) == fp:
        say("anm: cached")
        return
    bridge = e1.AnmBridge(reg, SCHEMA, a.anm_root)
    meta = json.loads((cache / "prepare_meta.json").read_text())
    pan, pan_r, p3 = meta["panel_primary"], meta["panel_reduced"], meta["q3_anchors"]
    names = reg["evidence"]["teddy_head"]["adt_names"]
    j = {n: i for i, n in enumerate(names)}
    P = np.load(cache / "prepare.npz")
    v = P["v"]
    N = v.shape[0]
    parts = []
    for c0 in range(0, N, chunk_cells):
        cf = cache / f"anm_part{c0 // chunk_cells:04d}.npz"
        if cf.exists() and str(np.load(cf)["fingerprint"]) == fp:
            parts.append(cf)
            continue
        t0 = time.time()
        rows = range(c0, min(N, c0 + chunk_cells))
        n = len(rows)
        # tag "": primary (registered panels, 3 events per class); tag "r": reduced-panel sensitivity
        R = {f"{p}_scores": np.zeros((n, 4)) for p in ("q1", "cl", "q1r", "clr", "q3")}
        for p in ("q1", "q2", "cl_q1", "cl_q2", "q1r", "q2r", "clr_q1", "clr_q2", "q3"):
            R[f"{p}_call"] = np.full(n, -1, np.int8)
        R["n_rejected"] = np.zeros(n, np.int16)        # primary runs + Q3 anchors
        R["n_rejected_red"] = np.zeros(n, np.int16)    # reduced-panel runs
        for r, i in enumerate(rows):
            for tag, pn in (("", pan), ("r", pan_r)):
                vals = {k: [v[i, j[p]] for p in pn[k]] for k in vk.LINEAGES}
                o = bridge.run(pn, vals, readouts=("Q1", "Q2"), closure_readouts=("closure_Q1", "closure_Q2"))
                R[f"q1{tag}_scores"][r], R[f"q1{tag}_call"][r] = o["Q1"]
                R[f"q2{tag}_call"][r] = o["Q2"][1]
                R[f"cl{tag}_scores"][r], R[f"cl{tag}_q1_call"][r] = o["closure_Q1"]
                R[f"cl{tag}_q2_call"][r] = o["closure_Q2"][1]
                R["n_rejected" if tag == "" else "n_rejected_red"][r] += o["n_rejected"]
            o3 = bridge.run(p3, {k: [v[i, j[p]] for p in p3[k]] for k in vk.LINEAGES}, readouts=("Q3",))
            R["q3_scores"][r], R["q3_call"][r] = o3["Q3"]
            R["n_rejected"][r] += o3["n_rejected"]
        atomic_savez(cf, fingerprint=np.asarray(fp), **R)
        parts.append(cf)
        say(f"anm: cells {c0}-{c0 + n - 1} of {N} in {time.time() - t0:.1f}s")
    merged = {k: np.concatenate([np.load(p)[k] for p in parts]) for k in np.load(parts[0]).files if k != "fingerprint"}
    atomic_savez(out, fingerprint=np.asarray(fp), anm_commit=np.asarray(str(bridge.anm_commit())),
                 thresholds=np.asarray(json.dumps(bridge.thresholds)), **merged)
    for p in parts:
        p.unlink()
    say(f"anm: done ({N} cells)")


# ============================================================================ evaluator
def build_evaluator(P, A, reg) -> e6.Evaluator:
    amend = reg["amendment_A1"]
    grid = list(reg["experiments"]["common"]["matched_coverage"]["grid"])
    gn = [f"g{c:.2f}" for c in grid]
    ids = P["ids"]
    S, Sr = P["S"], P["S_red"]
    bar = float(reg["questions"]["Q1"]["bar"])
    keys = {k: P[f"key_{k}"].astype(str) for k in KEYS}
    M = e6.Method
    c12, k12 = e1.lineage_conf_call(P["P5_classifier"])
    c10, k10 = e1.lineage_conf_call(P["P5_classifier10"])
    idx = lambda calls: np.asarray([vk.LINEAGES.index(x) for x in calls])  # noqa: E731
    methods = [
        # primary: the full frozen registered panels; every arm reads the same 12 evidence values (as E1)
        M("anm", A["q1_scores"].max(axis=1), A["q1_scores"].argmax(axis=1), ids, amend, gn + ["c1"], KEYS),
        M("rule", S.max(axis=1), S.argmax(axis=1), ids, amend, gn + ["c1"], KEYS, C2_COVS),
        M("margin", vk.margin(S), S.argmax(axis=1), ids, amend, gn, KEYS),
        M("entropy", e1.entropy_confidence(S), S.argmax(axis=1), ids, amend, gn, KEYS),
        M("classifier", c12, idx(k12), ids, amend, gn + ["c1"], KEYS, C2_COVS),
        M("closure", A["cl_scores"].max(axis=1), A["cl_scores"].argmax(axis=1), ids, amend, ["c1"], KEYS, C2_COVS),
        # sensitivity: the literal reduced-panel reading (never decides)
        M("anm_red", A["q1r_scores"].max(axis=1), A["q1r_scores"].argmax(axis=1), ids, amend, gn + ["c1r"], KEYS),
        M("rule_red", Sr.max(axis=1), Sr.argmax(axis=1), ids, amend, gn + ["c1r"], KEYS, C2_COVS_RED),
        M("margin_red", vk.margin(Sr), Sr.argmax(axis=1), ids, amend, gn, KEYS),
        M("entropy_red", e1.entropy_confidence(Sr), Sr.argmax(axis=1), ids, amend, gn, KEYS),
        M("classifier10", c10, idx(k10), ids, amend, gn + ["c1r"], KEYS),
        M("closure_red", A["clr_scores"].max(axis=1), A["clr_scores"].argmax(axis=1), ids, amend, ["c1r"], KEYS, C2_COVS_RED),
    ]
    realised = {"c1": (S.max(axis=1) >= bar).astype(np.float64), "c1r": (Sr.max(axis=1) >= bar).astype(np.float64)}
    return e6.Evaluator(methods, keys, grid, realised, KEYS)


def groups_of(donor: np.ndarray) -> tuple[list[str], list[np.ndarray]]:
    names = sorted(set(donor.tolist()))
    return names, [np.where(donor == d)[0] for d in names]


# ============================================================================ stage: boot
def stage_boot(a, reg, cfg, fp) -> None:
    cache, bdir = a.out_dir / "cache", a.out_dir / "boot"
    bdir.mkdir(parents=True, exist_ok=True)
    final, ffp = bdir / "pooled.npy", bdir / "pooled.fp"
    if final.exists() and ffp.exists() and ffp.read_text() == fp:
        say("boot: cached")
        return
    P, A = np.load(cache / "prepare.npz"), np.load(cache / "anm.npz")
    ev = build_evaluator(P, A, reg)
    names = sorted(ev.stats(np.ones(ev.N)))
    (bdir / "names.json").write_text(json.dumps(names))
    _, groups = groups_of(P["donor"].astype(str))
    B, ch = int(cfg["n_boot"]), int(a.chunk)
    gen = e1.replicate_counts(groups, ev.N, B, cfg["boot_seed"], two_stage=True)
    t0 = time.time()
    rows = []
    for c0 in range(0, B, ch):
        cf = bdir / f"pooled_chunk{c0 // ch:04d}.npy"
        cfp = cf.with_suffix(".fp")
        ws = [next(gen) for _ in range(min(ch, B - c0))]   # always drawn: the sequence is the same on resume
        if cf.exists() and cfp.exists() and cfp.read_text() == fp:
            rows.append(np.load(cf))
            continue
        Mx = np.asarray([[st[k] for k in names] for st in (ev.stats(w) for w in ws)])
        np.save(cf, Mx)
        cfp.write_text(fp)
        rows.append(Mx)
        say(f"boot: {c0 + len(ws)}/{B} replicates, {time.time() - t0:.0f}s elapsed")
    np.save(final, np.concatenate(rows, axis=0))
    ffp.write_text(fp)
    for f in bdir.glob("pooled_chunk*"):
        f.unlink()
    say("boot: done")


# ============================================================================ stage: report
E1_PATHS = {
    "E1.4a anm - margin AURC": ("E1.4", "keys", "primary", "top_vs_margin", "E1.4a_AURC", "point"),
    "E1.4a anm - classifier AURC": ("E1.4", "keys", "primary", "top_vs_classifier", "E1.4a_AURC", "point"),
    "E1.4a anm - entropy AURC": ("E1.4", "keys", "primary", "top_vs_entropy", "E1.4a_AURC", "point"),
    "E1.4a rule - margin AURC": ("E1.4", "keys", "primary", "top_vs_margin", "E1.4a_AURC", "point"),
    "Q1 acc anm - margin @0.90": ("E1.4", "keys", "primary", "top_vs_margin", "E1.4b_acc@0.90", "point"),
    "Q1 acc anm - margin @0.70": ("E1.4", "keys", "primary", "top_vs_margin", "E1.4b_acc@0.70", "point"),
    "Q1 acc anm - classifier @0.90": ("E1.4", "keys", "primary", "top_vs_classifier", "E1.4b_acc@0.90", "point"),
    "Q1 acc anm - classifier @0.70": ("E1.4", "keys", "primary", "top_vs_classifier", "E1.4b_acc@0.70", "point"),
    "E1.C2 closure - rule OUT decline @c*": ("C2", "keys", "primary", "c_star", "closure_minus_rule_out_decline", "point"),
    "E1.C2 closure - rule OUT decline @0.80": ("C2", "keys", "primary", "0.80", "closure_minus_rule_out_decline", "point"),
    "E1.C2 closure - rule OUT decline @0.70": ("C2", "keys", "primary", "0.70", "closure_minus_rule_out_decline", "point"),
}
PRIMARY_ROWS = ("E1.4a anm - margin AURC", "E1.C2 closure - rule OUT decline @c*")
C2_LABELS = ("c*", "0.80", "0.70")


def primary_rows(Rk: dict) -> dict:
    """The replication rows of one key's primary reading (registered panels)."""
    return {**{f"E1.4a {k}": v for k, v in Rk["E1.4a"].items() if "replication" in v}, **Rk["E1.4b"],
            **{f"E1.C2 closure - rule OUT decline @{lab}": Rk["E1.C2"][lab]["closure - rule OUT decline"] for lab in C2_LABELS}}


def sensitivity_rows(Rk: dict) -> dict:
    """The same rows on the reduced-panel sensitivity reading (never decides)."""
    S = Rk["sensitivity_reduced_panels"]
    return {**{f"E1.4a {k}": v for k, v in S["E1.4a"].items() if "replication" in v}, **S["E1.4b"],
            **{f"E1.C2 closure_red - rule_red OUT decline @{lab}": S["E1.C2"][lab]["closure_red - rule_red OUT decline"]
               for lab in C2_LABELS}}


def load_v3_directions(a) -> dict:
    p = a.e1_results
    if not p.exists():
        return {"file": str(p), "available": False, "points": {}, "note": "E1 results missing: replication pending"}
    r = json.loads(p.read_text())
    if bool(r.get("smoke")) and not a.smoke:
        raise SystemExit(f"{p} is a smoke E1 run; the v3 direction must come from E1's site4 run")
    sp = r.get("splits", {}).get("primary", {})
    pts = {k: e6.dig(sp, path) for k, path in E1_PATHS.items()}
    return {"file": str(p), "sha256": vk.sha256_file(p), "available": True, "smoke": bool(r.get("smoke")),
            "split": sp.get("split"), "registration_sha256": r.get("registration_sha256"), "points": pts,
            "E1.1a_passed": e6.dig(sp, ("E1.1a", "passed"))}


class View:
    def __init__(self, point, donors: dict, boot: np.ndarray, names: list):
        self.point, self.donors = point, donors
        self.col = {k: i for i, k in enumerate(names)}
        self.boot = boot

    def g_point(self, k):
        return self.point[k]

    def g_boot(self, k):
        return self.boot[:, self.col[k]]


def block(V: View, fn, margin: float | None = None, v3_point=None, replicate: bool = False) -> dict:
    point = float(fn(V.g_point))
    per = {d: float(fn(lambda k, d=d: V.donors[d][k])) for d in V.donors}
    boot = np.asarray(fn(V.g_boot), dtype=np.float64)
    lo, hi, nf = e1.percentile_ci(boot)
    out = {"point": point, "ci95": [lo, hi], "n_boot_finite": nf, "per_donor": per}
    if margin is not None:
        ver = e1.verdict(point, boot, per, margin)
        out["e1_decision_rule"] = {"verdict": ver["verdict"], "margin": margin,
                                   "note": "registered E1 rule with each of the external donors as 'each primary donor'"}
    if replicate:
        out["replication"] = e6.replication(v3_point, point, [lo, hi], per)
    return out


def call_table(calls: np.ndarray, ct: np.ndarray) -> dict:
    out = {}
    for t in sorted(set(ct.tolist())):
        m = ct == t
        out[t] = {"n": int(m.sum()), **{c: int(np.sum(calls[m] == c)) for c in vk.LINEAGES},
                  "no_call": int(np.sum(calls[m] == vk.NO_CALL))}
    return out


def operating_point(calls: np.ndarray, keys: dict, donor: np.ndarray, ct: np.ndarray) -> dict:
    out = {"coverage": float(np.mean(calls != vk.NO_CALL)), "per_type_calls": call_table(calls, ct)}
    for kn in KEYS:
        r = vk.selective_accuracy(calls, keys[kn])
        r["per_donor"] = {d: vk.selective_accuracy(calls[donor == d], keys[kn][donor == d]) for d in sorted(set(donor.tolist()))}
        out[kn] = r
    return out


def stage_report(a, reg, add, cfg, fp, st, ext_info) -> dict:
    cache, bdir = a.out_dir / "cache", a.out_dir / "boot"
    P, A = np.load(cache / "prepare.npz"), np.load(cache / "anm.npz")
    kv = json.loads((a.out_dir / "key_validity.json").read_text())
    ev = build_evaluator(P, A, reg)
    names = json.loads((bdir / "names.json").read_text())
    donor = P["donor"].astype(str)
    dnames, _ = groups_of(donor)
    point = ev.stats(np.ones(ev.N))
    per = {d: ev.stats((donor == d).astype(np.float64)) for d in dnames}
    V = View(point, per, np.load(bdir / "pooled.npy"), names)
    grid = list(reg["experiments"]["common"]["matched_coverage"]["grid"])
    gn = [f"g{c:.2f}" for c in grid]
    E = reg["experiments"]["E1"]
    m4 = float(E["exp4_which_calls_to_trust"]["margin"])
    m2 = float(E["c2_nested_readouts"]["margin"])
    v3 = load_v3_directions(a)
    v3p = v3["points"]
    dk = kv["decision_key"]
    keys = {k: P[f"key_{k}"].astype(str) for k in KEYS}
    ct = P["ct"].astype(str)

    res: dict = {"experiment": "E6", "builder": VERSION, "smoke": cfg["smoke"], "contamination_label": e6.CONTAMINATION_LABEL,
                 "registration_sha256": st["registration_sha256"], "amendments_sha256": st["amendments_sha256"],
                 "addendum_E6_sha256": st["addendum_sha256"], "state": st, "external_inputs": ext_info, "config": cfg,
                 "cells": {"n": int(ev.N), "per_donor": {d: int(np.sum(donor == d)) for d in dnames}},
                 "key_validity": {k: kv[k] for k in ("kappa5", "threshold", "validated", "decision_key", "label",
                                                     "n_cells_with_annotation_class", "n_annotation_unscored")}
                 | {"pooled": kv["pooled"], "per_donor_kappa": {d: r["kappa5"] for d, r in kv["per_donor"].items()},
                    "key_counts": kv["key_counts"]},
                 "decision_key": dk, "decision_key_label": kv["label"], "v3_directions": v3}
    if cfg["smoke"]:
        res["SMOKE"] = "val donor 18303 cells with 8 label-free pseudo-donors and small B; not a result"

    # ---------------- E1.1a: primary on the registered panels (as E1); reduced-panel sensitivity beside it
    thr = json.loads(str(A["thresholds"]))
    bars = {q: float(reg["questions"][q]["bar"]) for q in ("Q1", "Q2", "Q3")}
    rule = {"Q1": vk.rule_calls(P["S"], bars["Q1"]), "Q2": vk.rule_calls(P["S"], bars["Q2"]),
            "Q3": vk.rule_calls(P["S3"], bars["Q3"])}
    rule_r = {q: vk.rule_calls(P["S_red"], bars[q]) for q in ("Q1", "Q2")}
    eng = {"Q1": e6.engine_calls(A["q1_call"]), "Q2": e6.engine_calls(A["q2_call"]), "Q3": e6.engine_calls(A["q3_call"]),
           "closure_Q1": e6.engine_calls(A["cl_q1_call"]), "closure_Q2": e6.engine_calls(A["cl_q2_call"]),
           "Q1_red": e6.engine_calls(A["q1r_call"]), "Q2_red": e6.engine_calls(A["q2r_call"]),
           "closure_Q1_red": e6.engine_calls(A["clr_q1_call"]), "closure_Q2_red": e6.engine_calls(A["clr_q2_call"])}
    rec = {"Q1": e6.calls_at(P["rec_q1"], thr["Q1"]), "Q2": e6.calls_at(P["rec_q1"], thr["Q2"]),
           "Q3": e6.calls_at(P["rec_q3"], thr["Q3"]),
           "closure_Q1": e6.calls_at(P["rec_cl"], thr["closure_Q1"]), "closure_Q2": e6.calls_at(P["rec_cl"], thr["closure_Q2"]),
           "Q1_red": e6.calls_at(P["rec_q1_red"], thr["Q1"]), "Q2_red": e6.calls_at(P["rec_q1_red"], thr["Q2"]),
           "closure_Q1_red": e6.calls_at(P["rec_cl_red"], thr["closure_Q1"]),
           "closure_Q2_red": e6.calls_at(P["rec_cl_red"], thr["closure_Q2"])}
    # E1's closure check: engine closure calls vs the re-coded closure readout at the registered closure bars
    cl_rule = {q: vk.rule_calls(P["rec_cl"], float(reg["anm"][f"closure_bar_{q}"])) for q in ("Q1", "Q2")}

    def counts(x, y):
        return {"pooled": int(np.sum(x != y)), "per_donor": {d: int(np.sum((x != y)[donor == d])) for d in dnames},
                "share": float(np.mean(x != y))}

    fld = reg["anm"]["field_representation"]
    G1, G3 = vk.field_gain(1, fld), vk.field_gain(3, fld)
    meth = {m.name: m for m in ev.methods}
    e11a = {"registered_value": 0,
            "panels": "registered panels (3 events per class), as E1",
            "engine_vs_mean_rule": {q: counts(eng[q], rule[q]) for q in ("Q1", "Q2", "Q3")},
            "closure_engine_vs_recoded": {f"closure_{q}": counts(eng[f"closure_{q}"], cl_rule[q]) for q in ("Q1", "Q2")},
            "bridge_engine_vs_recoded": {q: counts(eng[q], rec[q]) for q in ("Q1", "Q2", "Q3", "closure_Q1", "closure_Q2")},
            "max_abs_score_diff_engine_vs_recoded": {nm: float(np.max(np.abs(A[f"{ak}_scores"] - P[pk]))) for nm, ak, pk in
                                                     (("Q1", "q1", "rec_q1"), ("closure", "cl", "rec_cl"), ("Q3", "q3", "rec_q3"))},
            "max_rel_diff_engine_vs_Gn_times_S": {
                "Q1": float(np.max(np.abs(A["q1_scores"] - G3 * 3 * P["S"]) / np.maximum(1e-12, np.abs(A["q1_scores"])))),
                "Q3": float(np.max(np.abs(A["q3_scores"] - G1 * P["S3"]) / np.maximum(1e-12, np.abs(A["q3_scores"]))))},
            "engine_events_rejected": int(A["n_rejected"].sum()),
            "ranking_identical_anm_vs_rule": bool(np.array_equal(meth["anm"].order, meth["rule"].order)),
            "anm_commit": str(A["anm_commit"]), "thresholds": thr}
    e11a["passed"] = bool(all(v["pooled"] == 0 for grp in ("engine_vs_mean_rule", "closure_engine_vs_recoded",
                                                          "bridge_engine_vs_recoded") for v in e11a[grp].values())
                          and e11a["engine_events_rejected"] == 0)
    e11a["replication"] = {"v3_E1.1a_passed": v3.get("E1.1a_passed"),
                           "e6_zero_mismatches_every_donor": e11a["passed"],
                           "replicates": (None if v3.get("E1.1a_passed") is None else bool(v3["E1.1a_passed"] and e11a["passed"]))}
    e11a_r = {"status": SENS, "panels": "reduced panels (B 3, T 2, NK 2, myeloid 3 events)",
              "bridge_engine_vs_recoded": {q: counts(eng[q], rec[q]) for q in ("Q1_red", "Q2_red", "closure_Q1_red", "closure_Q2_red")},
              "max_abs_score_diff_engine_vs_recoded": {nm: float(np.max(np.abs(A[f"{ak}_scores"] - P[pk]))) for nm, ak, pk in
                                                       (("Q1_red", "q1r", "rec_q1_red"), ("closure_red", "clr", "rec_cl_red"))},
              "engine_events_rejected": int(A["n_rejected_red"].sum()),
              "engine_vs_mean_rule": {f"{q}_red": counts(eng[f"{q}_red"], rule_r[q]) for q in ("Q1", "Q2")},
              "ranking_identical_anm_red_vs_rule_red": bool(np.array_equal(meth["anm_red"].order, meth["rule_red"].order)),
              "anm_unequal_panels_note": add["sensitivity"]["reduced_panel_reading"]["arms"]["anm_unequal_panels"]}
    e11a_r["bridge_passed"] = bool(all(v["pooled"] == 0 for v in e11a_r["bridge_engine_vs_recoded"].values())
                                   and e11a_r["engine_events_rejected"] == 0)
    e11a["sensitivity_reduced_panels"] = e11a_r
    res["E1.1a"] = e11a

    # ---------------- per key: Q1 selective accuracy, E1.4a/b, C2 (primary), reduced-panel sensitivity
    arms_p = ("anm", "rule", "margin", "entropy", "classifier")
    arms_r = ("anm_red", "rule_red", "margin_red", "entropy_red", "classifier10")
    per_key: dict = {}
    for kn in KEYS:
        R: dict = {"key": kn, "is_decision_key": kn == dk}

        def dA(x, y, kn=kn):
            return lambda g: e6.aurc_from(g, x, kn, gn, grid) - e6.aurc_from(g, y, kn, gn, grid)

        def dacc(x, y, c, kn=kn):
            return lambda g: np.asarray(g(f"{x}|g{c}|{kn}|acc")) - np.asarray(g(f"{y}|g{c}|{kn}|acc"))

        def selacc(arms, cstar, kn=kn):
            return {arm: {c: block(V, lambda g, arm=arm, c=c: g(f"{arm}|{c}|{kn}|acc"))
                          for c in ([cstar] if cstar in meth[arm].acc_covs else []) + gn} for arm in arms}

        def aurcs(arms, kn=kn):
            return {arm: block(V, lambda g, arm=arm: e6.aurc_from(g, arm, kn, gn, grid)) for arm in arms}

        def c2_rows(closure, rule_, covs, with_classifier, kn=kn):
            c2: dict = {}
            for c, lab in zip(covs, ("c*", "0.80", "0.70")):
                arms = (closure, rule_) + (("classifier",) if with_classifier else ())
                cb: dict = {arm: {"out_decline": block(V, lambda g, arm=arm, c=c: g(f"{arm}|{c}|{kn}|outdecl")),
                                  "inscope_acc": block(V, lambda g, arm=arm, c=c: g(f"{arm}|{c}|{kn}|inscope"))}
                            for arm in arms}
                d_out, d_in = f"{closure} - {rule_} OUT decline", f"{closure} - {rule_} in-scope acc"
                cb[d_out] = block(V, lambda g, c=c: np.asarray(g(f"{closure}|{c}|{kn}|outdecl"))
                                  - np.asarray(g(f"{rule_}|{c}|{kn}|outdecl")), m2,
                                  v3p.get(f"E1.C2 closure - rule OUT decline @{lab}"), True)
                cb[d_in] = block(V, lambda g, c=c: np.asarray(g(f"{closure}|{c}|{kn}|inscope"))
                                 - np.asarray(g(f"{rule_}|{c}|{kn}|inscope")))
                if with_classifier:
                    cb["classifier - rule OUT decline (secondary)"] = block(
                        V, lambda g, c=c: np.asarray(g(f"classifier|{c}|{kn}|outdecl")) - np.asarray(g(f"{rule_}|{c}|{kn}|outdecl")))
                guard = cb[d_in]["point"] >= -0.005
                cb["guard_inscope_ok"] = bool(guard)
                vo = cb[d_out]["e1_decision_rule"]["verdict"]
                cb["e1_decision_rule_with_guard"] = ("win" if guard else "not a win (in-scope guard failed)") if vo == "win" else vo
                c2[lab] = cb
            return c2

        # primary: registered panels, the same 12 evidence values for every arm (as E1)
        R["coverage_c_star"] = block(V, lambda g: g("cov|c1"))
        R["Q1_selective_accuracy"] = selacc(arms_p, "c1")
        R["AURC"] = aurcs(arms_p)
        R["E1.4a"] = {
            "anm - margin AURC": block(V, dA("anm", "margin"), m4, v3p.get("E1.4a anm - margin AURC"), True),
            "anm - classifier AURC": block(V, dA("anm", "classifier"), m4, v3p.get("E1.4a anm - classifier AURC"), True),
            "anm - entropy AURC": block(V, dA("anm", "entropy"), m4, v3p.get("E1.4a anm - entropy AURC"), True),
            "rule - margin AURC": block(V, dA("rule", "margin"), m4, v3p.get("E1.4a rule - margin AURC"), True),
            "anm - rule AURC": block(V, dA("anm", "rule"), m4),
        }
        R["E1.4b"] = {f"Q1 acc {x} - {y} @{c}": block(V, dacc(x, y, c), m4, v3p.get(f"Q1 acc {x} - {y} @{c}"), True)
                      for x, y in (("anm", "margin"), ("anm", "classifier")) for c in ("0.90", "0.70")}
        R["E1.C2"] = c2_rows("closure", "rule", C2_COVS, True)

        # sensitivity: the literal reduced-panel reading (never decides); v3 directions of the primary counterparts
        R["sensitivity_reduced_panels"] = {
            "status": SENS,
            "coverage_c_star_reduced": block(V, lambda g: g("cov|c1r")),
            "Q1_selective_accuracy": selacc(arms_r, "c1r"),
            "AURC": aurcs(arms_r),
            "E1.4a": {
                "anm_red - margin_red AURC": block(V, dA("anm_red", "margin_red"), m4, v3p.get("E1.4a anm - margin AURC"), True),
                "anm_red - classifier AURC": block(V, dA("anm_red", "classifier"), m4, v3p.get("E1.4a anm - classifier AURC"), True),
                "anm_red - entropy_red AURC": block(V, dA("anm_red", "entropy_red"), m4, v3p.get("E1.4a anm - entropy AURC"), True),
                "rule_red - margin_red AURC": block(V, dA("rule_red", "margin_red"), m4, v3p.get("E1.4a rule - margin AURC"), True),
                "anm_red - rule_red AURC": block(V, dA("anm_red", "rule_red"), m4),
                "anm_red - classifier10 AURC": block(V, dA("anm_red", "classifier10"), m4,
                                                     v3p.get("E1.4a anm - classifier AURC"), True),
            },
            "E1.4b": {f"Q1 acc {x} - {y} @{c}": block(V, dacc(x, y, c), m4, v3p.get(f"Q1 acc anm - {yv} @{c}"), True)
                      for x, y, yv in (("anm_red", "margin_red", "margin"), ("anm_red", "classifier", "classifier"),
                                       ("anm_red", "classifier10", "classifier")) for c in ("0.90", "0.70")},
            "E1.C2": c2_rows("closure_red", "rule_red", C2_COVS_RED, False),
        }
        per_key[kn] = R
    res["keys"] = per_key

    # ---------------- deployed operating points (descriptive)
    cbars = reg["classifier"]["primary"]["bars"]
    cc, ck = e1.lineage_conf_call(P["P5_classifier"])
    ops = {}
    for q in ("Q1", "Q2"):
        ops[f"rule_{q}"] = operating_point(rule[q], keys, donor, ct)
        ops[f"anm_{q}"] = operating_point(eng[q], keys, donor, ct)
        ops[f"anm_closure_{q}"] = operating_point(eng[f"closure_{q}"], keys, donor, ct)
        ops[f"classifier_{q}"] = operating_point(np.where(cc >= float(cbars[q]), ck, vk.NO_CALL), keys, donor, ct)
    for q in ("Q1", "Q2"):   # sensitivity: reduced panels
        ops[f"rule_red_{q}"] = operating_point(rule_r[q], keys, donor, ct)
        ops[f"anm_red_{q}"] = operating_point(eng[f"{q}_red"], keys, donor, ct)
        ops[f"anm_closure_red_{q}"] = operating_point(eng[f"closure_{q}_red"], keys, donor, ct)
    res["operating_points"] = ops

    def same_acc(nm, cov):
        acc = ops[nm][dk]["accuracy_of_calls"]
        return None if acc is None else bool(abs(point[f"{nm.rsplit('_', 1)[0]}|{cov}|{dk}|acc"] - acc) < 1e-12)

    res["consistency_checks"] = {"rule_Q1_at_c_star_equals_operating_point": same_acc("rule_Q1", "c1"),
                                 "rule_red_Q1_at_c_star_reduced_equals_operating_point": same_acc("rule_red_Q1", "c1r"),
                                 "anm_Q1_calls_equal_rule_Q1_calls": bool(np.array_equal(eng["Q1"], rule["Q1"]))}

    # ---------------- summary (decision key)
    Rk = per_key[dk]
    status_of = {True: "replicates", False: "does not replicate", None: "no v3 direction"}
    summ: dict = {"contamination": e6.CONTAMINATION_LABEL,
                  "key": ("validated" if kv["validated"] else e6.KEY_NOT_VALIDATED)
                         + f" (pooled kappa {kv['kappa5']:.4f}; decision key {dk})",
                  "panels": "primary: the full frozen registered panels for every arm (addendum v2, as E1); "
                            "sensitivity: reduced panels (never decides)",
                  "E1.1a": "pass" if e11a["passed"] else "FAIL: E6 not interpretable (registered: stop)",
                  "E1.1a_engine_vs_mean_rule_Q1": e11a["engine_vs_mean_rule"]["Q1"]["pooled"],
                  "replication": {}}

    def rep_rows(rows: dict, primary: tuple) -> dict:
        out = {}
        for k, b in rows.items():
            r = b["replication"]
            out[k] = {"status": r["status"], "v3_point": r["v3_point"], "e6_point": r["e6_point"],
                      "e6_ci95": r["e6_ci95"], "donors_same_sign": r.get("n_donors_same_sign"),
                      "e1_decision_rule": b["e1_decision_rule"]["verdict"], "primary_row": k in primary}
        return out

    summ["replication"] = rep_rows(primary_rows(Rk), PRIMARY_ROWS)
    summ["replication"]["E1.1a"] = {"status": status_of[e11a["replication"]["replicates"]]}
    summ["sensitivity_reduced_panels"] = {
        "status": SENS,
        "E1.1a_bridge_reduced": "pass" if e11a_r["bridge_passed"] else "discrepancy (reduced reading only)",
        "E1.1a_engine_vs_mean_rule_Q1_reduced": e11a_r["engine_vs_mean_rule"]["Q1_red"]["pooled"],
        "replication": rep_rows(sensitivity_rows(Rk), ())}
    if not v3["available"]:
        summ["replication_note"] = "E1 site4 results missing: replication pending (rerun --stage report later)"
    if not kv["validated"]:
        summ["label"] = e6.KEY_NOT_VALIDATED
    res["summary"] = summ
    res["addendum_disclosures"] = add["disclosures"]
    res["addendum_revision"] = {k: add.get("revision", {}).get(k) for k in ("addendum_version", "change", "reason", "timing")}
    return res


# ============================================================================ REPORT.md
def _f(x, n=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{n}f}"


def _ci(c, n=3):
    return f"[{_f(c[0], n)}, {_f(c[1], n)}]"


def _sign_donors(b):
    pd = b["per_donor"]
    return f"{sum(1 for v in pd.values() if v is not None and np.isfinite(v) and v > 0)}+ / {sum(1 for v in pd.values() if v is not None and np.isfinite(v) and v < 0)}-"


def write_report(path: Path, res: dict) -> None:
    L = []
    kvd = res["key_validity"]
    lab = res["contamination_label"]
    title = f"E6 (v3): external confirmation on Hao et al. 2021 PBMC CITE-seq — {lab}"
    if res["smoke"]:
        title += " — SMOKE TEST, not a result"
    L += [f"# {title}", ""]
    tags = [f"**{lab}.**"]
    if not kvd["validated"]:
        tags.append(f"**{e6.KEY_NOT_VALIDATED}**: every result below is on the annotation-only key.")
    L += [" ".join(tags), "",
          f"Registration `{res['registration_sha256'][:16]}`, E6 addendum `{(res['addendum_E6_sha256'] or 'none')[:16]}`; "
          f"builder {res['builder']}; git `{(res['state'].get('git_head') or '')[:10]}`; "
          f"{res['cells']['n']:,} cells, {len(res['cells']['per_donor'])} donors.", ""]
    if res["smoke"]:
        L += [f"**Smoke run:** {res['SMOKE']}.", ""]
    L += ["## Key validity (computed before any method result)", "",
          f"Pooled 5-class kappa, annotation vs per-donor gate (registered estimator within each donor, reduced gate), on "
          f"{kvd['n_cells_with_annotation_class']:,} cells with an annotation class ({kvd['n_annotation_unscored']:,} Doublet "
          f"cells unscored): **{_f(kvd['kappa5'], 4)}** (threshold {kvd['threshold']}) → "
          f"**{'validated' if kvd['validated'] else e6.KEY_NOT_VALIDATED}**; decision key **{res['decision_key']}**.", "",
          "| donor | kappa |", "|---|---:|"]
    L += [f"| {d} | {_f(k, 4)} |" for d, k in kvd["per_donor_kappa"].items()]
    L += ["", "| class | annotation | gated | agree | recall | precision |", "|---|---:|---:|---:|---:|---:|"]
    for c, r in kvd["pooled"]["per_class"].items():
        L.append(f"| {c} | {r['n_annotation']:,} | {r['n_gated']:,} | {r['n_agree']:,} | {_f(r['recall_vs_annotation'])} | "
                 f"{_f(r['precision_vs_annotation'])} |")
    L += ["", "Key counts: " + "; ".join(f"{k}: " + ", ".join(f"{c} {n:,}" for c, n in v.items()) for k, v in kvd["key_counts"].items()),
          "", "Per-type gate table and per-donor thresholds: `key_validity.json`.", ""]

    s = res["summary"]
    Rk = res["keys"][res["decision_key"]]

    def rep_table(rows: dict, primary: tuple) -> list:
        T = ["| endpoint | v3 point | E6 pooled [95% CI] | donors +/- | replication | E1 decision rule (8 donors) |",
             "|---|---:|---|---|---|---|"]
        for k, b in rows.items():
            r = b["replication"]
            star = " **(primary)**" if k in primary else ""
            T.append(f"| {k}{star} | {_f(r['v3_point'], 4)} | {_f(b['point'], 4)} {_ci(b['ci95'], 4)} | {_sign_donors(b)} | "
                     f"{r['status']} | {b.get('e1_decision_rule', {}).get('verdict', '')} |")
        return T

    a = res["E1.1a"]
    ar = a["sensitivity_reduced_panels"]
    n_mm = sum(v["pooled"] for g in ("engine_vs_mean_rule", "closure_engine_vs_recoded", "bridge_engine_vs_recoded")
               for v in a[g].values())
    L += ["## Replication of the v3 (E1 site4) conclusions", "",
          f"Decision key {res['decision_key']}; panels: the full frozen registered panels for every arm (addendum v2, as E1), "
          "so every arm reads the same 12 evidence values. A conclusion replicates if the E6 difference has the v3 sign in a "
          "majority of the external donors and the pooled two-stage 95% interval excludes 0 on that side. "
          f"v3 direction from `{res['v3_directions']['file']}`" + ("" if res["v3_directions"]["available"] else " (missing: pending)") + ".", ""]
    L += rep_table(primary_rows(Rk), PRIMARY_ROWS)
    L.append(f"| E1.1a (ANM = mean rule, closure, bridge) | E1 passed: {a['replication']['v3_E1.1a_passed']} | "
             f"mismatches {n_mm} | | {s['replication']['E1.1a']['status']} | |")
    L += [""]
    L += ["## E1.1a ANM engine vs the mean rule (registered panels, as E1)", "",
          "Engine vs mean rule (registered 0): " + ", ".join(f"{q} {v['pooled']}" for q, v in a["engine_vs_mean_rule"].items())
          + "; engine closure vs re-coded closure: " + ", ".join(f"{q} {v['pooled']}" for q, v in a["closure_engine_vs_recoded"].items())
          + "; bridge (engine vs exact closed form): " + ", ".join(f"{q} {v['pooled']}" for q, v in a["bridge_engine_vs_recoded"].items())
          + f"; rejected events {a['engine_events_rejected']}; max |engine - closed form| "
          + ", ".join(f"{k} {v:.1e}" for k, v in a["max_abs_score_diff_engine_vs_recoded"].items())
          + f"; ranking identical ANM vs rule: {a['ranking_identical_anm_vs_rule']}. "
          + f"**{'pass' if a['passed'] else 'FAIL: E6 not interpretable'}**.", ""]
    for kn in (res["decision_key"],) + tuple(k for k in KEYS if k != res["decision_key"]):
        R = res["keys"][kn]
        L += [f"## Key `{kn}`" + (" (decision key)" if R["is_decision_key"] else " (other key)"), ""]
        L += [f"c* (rule's realised Q1 coverage, registered panels) {_f(R['coverage_c_star']['point'])} {_ci(R['coverage_c_star']['ci95'])}.", "",
              "### Q1 selective accuracy at matched coverage and AURC", "",
              "| arm | AURC [95% CI] | acc @c* | acc @0.90 | acc @0.70 |", "|---|---|---:|---:|---:|"]
        for arm, b in R["AURC"].items():
            qa = R["Q1_selective_accuracy"][arm]
            L.append(f"| {arm} | {_f(b['point'])} {_ci(b['ci95'])} | {_f(qa.get('c1', {}).get('point'))} | "
                     f"{_f(qa['g0.90']['point'])} | {_f(qa['g0.70']['point'])} |")
        L += ["", "| difference | pooled [95% CI] | per donor | E1 decision rule |", "|---|---|---|---|"]
        for k, b in {**R["E1.4a"], **R["E1.4b"]}.items():
            L.append(f"| {k} | {_f(b['point'], 4)} {_ci(b['ci95'], 4)} | " + ", ".join(f"{d.replace('hao2021_', '')} {_f(v, 3)}" for d, v in b["per_donor"].items())
                     + f" | {b.get('e1_decision_rule', {}).get('verdict', '')} |")
        L += ["", "### E1.C2 nested readouts: key-OUT decline at matched coverage", "",
              "| coverage | OUT decline closure / rule / classifier | closure - rule [95% CI] | in-scope acc diff | E1 rule (with guard) |",
              "|---|---|---|---:|---|"]
        for lab_, cb in R["E1.C2"].items():
            d = cb["closure - rule OUT decline"]
            L.append(f"| {lab_} | {_f(cb['closure']['out_decline']['point'])} / {_f(cb['rule']['out_decline']['point'])} / "
                     f"{_f(cb['classifier']['out_decline']['point'])} | {_f(d['point'], 4)} {_ci(d['ci95'], 4)} | "
                     f"{_f(cb['closure - rule in-scope acc']['point'], 4)} | {cb['e1_decision_rule_with_guard']} |")
        L += [""]

    # ---------------- sensitivity: reduced panels, beside the primary, never deciding
    Sk = Rk["sensitivity_reduced_panels"]
    L += ["## Sensitivity: reduced panels (literal missing-protein reading; never decides)", "",
          "T on CD3, CD2 and NK on CD122, CD56 (CD5 and CD94 dropped); arms `_red` on those panels, the frozen 12-feature "
          "classifier as the version-1 comparator, and classifier10 on the 10 reduced-panel features. Shown beside the "
          "primary; it never changes a primary verdict or replication status.", "",
          f"E1.1a on the reduced panels: bridge (engine vs exact closed form) "
          + ", ".join(f"{q} {v['pooled']}" for q, v in ar["bridge_engine_vs_recoded"].items())
          + f", rejected events {ar['engine_events_rejected']} (**{'pass' if ar['bridge_passed'] else 'discrepancy'}**); "
          f"engine vs mean rule Q1 {ar['engine_vs_mean_rule']['Q1_red']['pooled']:,} cells "
          f"({_f(ar['engine_vs_mean_rule']['Q1_red']['share'])}), Q2 {ar['engine_vs_mean_rule']['Q2_red']['pooled']:,} "
          f"({_f(ar['engine_vs_mean_rule']['Q2_red']['share'])}); ranking identical: {ar['ranking_identical_anm_red_vs_rule_red']}.",
          "", f"_{ar['anm_unequal_panels_note']}_", "",
          f"Decision key {res['decision_key']}; c*_red {_f(Sk['coverage_c_star_reduced']['point'])} "
          f"{_ci(Sk['coverage_c_star_reduced']['ci95'])}.", "",
          "| arm | AURC [95% CI] | acc @c*_red | acc @0.90 | acc @0.70 |", "|---|---|---:|---:|---:|"]
    for arm, b in Sk["AURC"].items():
        qa = Sk["Q1_selective_accuracy"][arm]
        L.append(f"| {arm} | {_f(b['point'])} {_ci(b['ci95'])} | {_f(qa.get('c1r', {}).get('point'))} | "
                 f"{_f(qa['g0.90']['point'])} | {_f(qa['g0.70']['point'])} |")
    L += ["", "Replication rule applied to the reduced rows (sensitivity; v3 directions of the primary counterparts):", ""]
    L += rep_table(sensitivity_rows(Rk), ())
    L += ["", f"Other rows (anm_red - rule_red AURC, in-scope differences, other key): `E6_results.json` "
          "(keys.*.sensitivity_reduced_panels).", ""]
    ops = res["operating_points"]
    dk = res["decision_key"]
    L += ["## Deployed operating points (bars from val; descriptive)", "",
          "| method | coverage | selective acc | decision acc | OUT decline |", "|---|---:|---:|---:|---:|"]
    for nm, op in ops.items():
        r = op[dk]
        L.append(f"| {nm} | {_f(op['coverage'])} | {_f(r['accuracy_of_calls'])} | {_f(r['decision_accuracy'])} | {_f(r['out_decline_rate'])} |")
    pt = ops["rule_Q1"]["per_type_calls"]
    if "Platelet" in pt:
        p = pt["Platelet"]
        L += ["", f"Platelet (keyed OUT; CD62P is in the myeloid panel, A1 review note): rule Q1 calls on {p['n']:,} platelets: "
              + ", ".join(f"{c} {p[c]:,}" for c in vk.LINEAGES) + f", no call {p['no_call']:,}."]
    L += ["", "Per-type call tables for every arm: `E6_results.json` (operating_points.*.per_type_calls).", ""]
    L += ["## Disclosures", ""] + [f"- **{k}**: {v}" for k, v in res["addendum_disclosures"].items()]
    L += ["", "Full numbers: `E6_results.json`; key validity: `key_validity.json`; log: `progress.log`; decisions: "
          "`registration/addenda/E6.json`.", ""]
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
    reg = va.load_registration_amended(a.registration)
    st = registration_state(a)
    guard_registration(st, a.smoke)
    add = json.loads(a.addendum.read_text())
    st["frozen_inputs_sha256"] = guard_frozen_inputs(a, add)
    a.external = a.external or Path(add["external_data"]["pack"]).parent
    a.external_embedding = a.external_embedding or Path(add["external_data"]["embedding"]["dir"])
    cfg = run_config(a, reg)
    say(f"E6 {VERSION} stage={a.stage} smoke={a.smoke} out={a.out_dir} ({e6.CONTAMINATION_LABEL})")
    ext_info = {"smoke": "BMMC val donor 18303; external pack and embedding not opened"} if a.smoke else guard_external(a, add)
    fp = fingerprint(a, st, cfg, ext_info)
    write_json(a.out_dir / "run_state.json", {**st, "fingerprint": fp, "config": cfg, "external_inputs": ext_info})
    t0 = time.time()
    if a.stage in ("all", "prepare"):
        D = load_inputs(a, reg, add)
        stage_prepare(a, reg, add, cfg, fp, D)
        del D
    if a.stage in ("all", "anm"):
        stage_anm(a, reg, fp)
    if a.stage in ("all", "boot"):
        stage_boot(a, reg, cfg, fp)
    if a.stage in ("all", "report"):
        res = stage_report(a, reg, add, cfg, fp, st, ext_info)
        res["seconds_this_invocation"] = time.time() - t0
        write_json(a.out_dir / "E6_results.json", res)
        write_report(a.out_dir / "REPORT.md", res)
        say(f"wrote {a.out_dir / 'E6_results.json'} and REPORT.md")
        for k, v in res["summary"].items():
            say(f"  {k}: {v}")
    say(f"done in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
