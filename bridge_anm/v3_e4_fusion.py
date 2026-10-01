#!/usr/bin/env python3
"""E4 (C6, Experiment 2 redesigned): key-disjoint two-channel fusion on the v3 registration.

Registered design: registration/registration_v3.json -> experiments.E4, loaded with
v3_amend.load_registration_amended() (A1.1 tie-break and per-replicate selection, A1.2 panel-free
channel-2 inputs, A1.9 comparator). The numbers the registration leaves open are fixed on train/val by
the ``prepare`` stage and written to registration/addenda/E4.json (+ E4_models.json, the fitted
coefficients); ``leakage`` proves they do not depend on site4; ``commit-addendum`` commits exactly those
files. No site4 row is read (RNA, protein, labels or embedding) before that commit.

  key        E4 key from the measured primary-panel proteins (Q1 scoring on measured evidence, OUT below
             tau_K; tau_K on val by kappa against the annotation). Secondary: the v3 primary key.
  channel 1  TEDDY + head prediction of the 12 panel proteins from RNA (registered evidence)
  channel 2  ridge from the 122 non-panel proteins, panel-free renormalised (A1.2), to the measured panel
  noise      L0 none; L1 RNA thinned to 0.2 (channel 1 re-embedded); L2 50% of channel-2 inputs set to 0;
             L3 both
  trust      per channel and level on val: mean max(0, Pearson(channel evidence, measured)) over 12 proteins
  methods    F0 channel 1 alone; F0b channel 2 alone; F1 simple average and F2 trust-weighted average
             (declared rules); F3 stacker (val); F4 learned fusion (training cells, L0); F5 TEDDY + ANM fusion
             (support + contradiction events, per-channel source scale = trust / max trust, all at t = 0)
  controls   C1 F5 without propagation; C2 F5 with channel 2 at t = 2 (declared order); C3 F5 without
             contradiction events (= F2); the re-coded closed form of F5 (must equal F5 on every cell)
  endpoint   selective accuracy at matched coverage 0.8 (0.9 secondary), mean over L0-L3, E4 key,
             test_primary: F5 vs the best of F1-F4 on val (A1.9) and F5 vs F2; margin 0.01; two-stage
             donor-then-cell bootstrap (B 2000, seed 1)

Stages
  embed --pool val    MPS/GPU: thin val RNA, re-embed with the official TEDDY-G pipeline (resumable shards)
  prepare             CPU, train/val only: key, channel 2, trust, F3/F4 fits, bars, comparator -> addendum
  leakage             CPU: rebuild the addendum core with site4 rows poisoned (must be identical) and with val
                      rows poisoned (must change)
  commit-addendum     git: commit E4.json, E4_models.json, E4_leakage_check.json and their HASHES.txt lines only
  embed --pool test   MPS/GPU: as embed for every site4 cell (refused before the addendum is committed)
  evaluate            CPU: every method on site4, bootstrap, verdicts -> E4_results.json + REPORT.md (refused
                      before the addendum is committed)
Smoke runs (--smoke NOTE) use val donor 18303 cells only: half of a small val subset plays val, the other
half plays the test split (two pseudo-donors); nothing is written to registration/.

Full run (from the dev clone root; R = the main checkout holding data/ and outputs/):
  PY=<venv>/bin/python; R=/Users/tianchichen/Documents/GitHub/teddy_mm; A=<ANM v2 fix checkout>
  $PY bridge_anm/v3_e4_fusion.py --stage embed --pool val --device mps --threads 4
  ANM_ROOT=$A $PY bridge_anm/v3_e4_fusion.py --stage prepare --threads 4
  ANM_ROOT=$A $PY bridge_anm/v3_e4_fusion.py --stage leakage --threads 4
  ANM_ROOT=$A $PY bridge_anm/v3_e4_fusion.py --stage commit-addendum
  $PY bridge_anm/v3_e4_fusion.py --stage embed --pool test --device mps --threads 4
  ANM_ROOT=$A $PY bridge_anm/v3_e4_fusion.py --stage evaluate --threads 4
Every stage appends to <out-dir>/progress.log; the embed stage writes shards of 512 cells and logs ms/cell and
an ETA after each; rerunning a stage skips finished shards (each carries a configuration hash).
"""
from __future__ import annotations

import os
import sys


def _early_threads(argv) -> int:
    n = 4
    for i, a in enumerate(argv):
        if a == "--threads" and i + 1 < len(argv):
            n = int(argv[i + 1])
        elif a.startswith("--threads="):
            n = int(a.split("=", 1)[1])
    return max(1, n)


_THREADS = _early_threads(sys.argv[1:])
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = str(_THREADS)

import argparse  # noqa: E402
import contextlib  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bridge_anm"))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e4 as e4  # noqa: E402
from lib import v3_key as vk  # noqa: E402

SCRIPT_VERSION = "v3_e4_fusion 1.0"
MAIN = Path(os.environ.get("TEDDY_MM_MAIN", "/Users/tianchichen/Documents/GitHub/teddy_mm"))
SHARD = 512
FLOOR_CHECK_CELLS = 64
THIN_FRACTION = 0.2
DROPOUT_FRACTION = 0.5
NO_CALL_RATE = 0.15
COVERAGES = (0.8, 0.9)  # 0.8 decides the verdict; 0.9 is reported with the same rule
C_GRID = (0.01, 0.1, 1.0, 10.0, 100.0)
N_FOLDS = 5
N_BOOT = 2000
METHODS = ("F0", "F0b", "F1", "F2", "F3", "F4", "F5", "C1", "C2", "C3")
COMPARATOR_CANDIDATES = ("F1", "F2", "F3", "F4")
ANM_VARIANTS = {"F5": {"t2": 0, "contradiction": True, "propagate": True},
                "C1": {"t2": 0, "contradiction": True, "propagate": False},
                "C2": {"t2": 2, "contradiction": True, "propagate": True},
                "C3": {"t2": 0, "contradiction": False, "propagate": True}}
ARM_NAMES = {
    "F0": "channel 1 alone: TEDDY + head evidence, Q1 rule",
    "F0b": "channel 2 alone: protein regressor evidence, Q1 rule",
    "F1": "simple average of the two channels' class scores (declared rule)",
    "F2": "trust-weighted average of the class scores (declared rule)",
    "F3": "stacker: logistic regression on the 8 class scores, trained on val per noise level",
    "F4": "learned fusion: logistic regression on the 24 evidence values, trained on training cells at L0",
    "F5": "TEDDY + ANM fusion (ANM field runs: per-channel source scale, contradiction events, all at t = 0)",
    "C1": "control: F5 without propagation (ANM no_propagation)",
    "C2": "control: F5 with channel 2 at t = 2 (declared order instead of time-blind)",
    "C3": "control: F5 without contradiction events (equals F2 by construction)",
}
ADDENDUM_FILES = ("E4.json", "E4_models.json", "E4_leakage_check.json")


# ============================================================================ helpers
_LOG: Path | None = None


def say(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG is not None:
        _LOG.parent.mkdir(parents=True, exist_ok=True)
        with _LOG.open("a") as f:
            f.write(line + "\n")


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def r6(x) -> float:
    return float(round(float(x), 6))


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(e4.jsonable(obj), indent=1, sort_keys=True))
    os.replace(tmp, path)


def atomic_savez(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def git(*args, env=None, input_bytes: bytes | None = None) -> str:
    r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, env=env, input=input_bytes)
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.decode(errors='replace')}")
    return r.stdout.decode()


def git_committed(path: Path) -> bool | None:
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


# ============================================================================ args
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", required=True,
                   choices=("embed", "prepare", "leakage", "commit-addendum", "evaluate"))
    p.add_argument("--pool", choices=("val", "test"), default=None, help="embed stage: which cells to re-embed")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official",
                   help="official embedding (z_rna.npy + z_rna_manifest.json)")
    p.add_argument("--ckpt", type=Path, default=MAIN.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=MAIN / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--anm-root", type=Path, default=Path(os.environ.get("ANM_ROOT", str(ROOT.parent / "ANM"))))
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--out-dir", type=Path, default=None, help=f"default {MAIN}/outputs/v3/E4 (smoke: required)")
    p.add_argument("--device", default="mps")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=32, help="embedding batch (official run: 32)")
    p.add_argument("--n-boot", type=int, default=N_BOOT)
    p.add_argument("--smoke", default=None, metavar="NOTE", help="smoke run on val cells only")
    p.add_argument("--smoke-n-val", type=int, default=320, help="smoke: val cells re-embedded (half val, half pseudo-test)")
    p.add_argument("--smoke-n-train", type=int, default=0, help="smoke: subsample of training cells (0 = all)")
    a = p.parse_args(argv)
    if a.stage == "embed" and a.pool is None:
        p.error("--stage embed needs --pool val|test")
    if a.smoke:
        if a.out_dir is None:
            p.error("smoke runs need an explicit --out-dir")
        if a.stage == "embed" and a.pool == "test":
            p.error("smoke runs never embed site4 cells")
        if a.stage == "commit-addendum":
            p.error("smoke runs never commit")
    else:
        if a.n_boot != N_BOOT or a.smoke_n_train:
            p.error("--n-boot / --smoke-n-train are for smoke runs only (give --smoke NOTE)")
    if a.out_dir is None:
        a.out_dir = MAIN / "outputs/v3/E4"
    return a


def addenda_dir(a) -> Path:
    return (a.out_dir / "smoke_addenda") if a.smoke else (a.registration_dir / "addenda")


# ============================================================================ data
def load_full(processed: Path, need_adt: bool = True) -> dict:
    """Metadata + protein of every row (the stage decides which rows it reads)."""
    z = np.load(processed / "cite_arrays.npz", allow_pickle=False)
    keys = ["split", "sites", "donors", "cell_types", "adt_names"] + (["adt"] if need_adt else [])
    out = {k: z[k] for k in keys}
    out["split"] = out["split"].astype(str)
    out["donors"] = out["donors"].astype(str)
    out["sites"] = out["sites"].astype(str)
    return out


def smoke_val_ids(meta: dict, reg: dict, n: int) -> np.ndarray:
    """Smoke only: a seeded subset of the val donor's cells (sorted global ids)."""
    vidx = vk.split_indices(meta["split"], meta["sites"], meta["donors"], reg, "val")
    return np.sort(np.random.default_rng(0).permutation(vidx)[: int(n)])


def pool_ids(meta: dict, reg: dict, pool: str) -> np.ndarray:
    if pool == "val":
        return vk.split_indices(meta["split"], meta["sites"], meta["donors"], reg, "val")
    tp = vk.split_indices(meta["split"], meta["sites"], meta["donors"], reg, "test_primary")
    ts = vk.split_indices(meta["split"], meta["sites"], meta["donors"], reg, "test_secondary")
    return np.sort(np.concatenate([tp, ts]))


# ============================================================================ registration gate
def addendum_status(a) -> dict:
    d = addenda_dir(a)
    hf = d / "HASHES.txt"
    lines = hf.read_text().splitlines() if hf.exists() else []
    info = {"addenda_dir": str(d), "files": {}}
    for name in ADDENDUM_FILES:
        f = d / name
        h = sha_bytes(f.read_bytes()) if f.exists() else None
        info["files"][name] = {"exists": f.exists(), "sha256": h, "hash_recorded": bool(h) and f"{h}  {name}" in lines,
                               "committed": (git_committed(f) if f.exists() else False) if not a.smoke else None}
    cf = a.registration_dir / "registration_v3.json"
    info["registration_sha256"] = vk.sha256_file(cf)
    info["registration_hash_file_matches"] = (a.registration_dir / "registration_v3.json.sha256").read_text().split()[0] == info["registration_sha256"]
    info["registration_committed"] = git_committed(cf)
    am = a.registration_dir / "amendment_A1.json"
    info["amendment_A1_sha256"] = vk.sha256_file(am)
    info["amendment_A1_hash_file_matches"] = (a.registration_dir / "amendment_A1.json.sha256").read_text().split()[0] == info["amendment_A1_sha256"]
    info["amendment_A1_committed"] = git_committed(am)
    return info


def require_committed_addendum(a) -> dict:
    """Refuse any site4 work unless the registration, A1 and the E4 addendum are committed and consistent."""
    info = addendum_status(a)
    if a.smoke:
        return info
    need = {
        "registration_v3.json committed, hash file matching": info["registration_committed"] is True and info["registration_hash_file_matches"],
        "amendment_A1.json committed, hash file matching": info["amendment_A1_committed"] is True and info["amendment_A1_hash_file_matches"],
    }
    for name in ADDENDUM_FILES:
        fi = info["files"][name]
        need[f"addenda/{name} exists, sha256 in addenda/HASHES.txt, committed"] = fi["exists"] and fi["hash_recorded"] and fi["committed"] is True
    if all(info["files"][n]["exists"] for n in ADDENDUM_FILES):
        add = json.loads((addenda_dir(a) / "E4.json").read_text())
        lk = json.loads((addenda_dir(a) / "E4_leakage_check.json").read_text())
        need["leakage check passed for this addendum core"] = bool(lk.get("passed")) and lk.get("addendum_core_sha256") == add["core_sha256"]
        need["E4_models.json is the one the addendum pins"] = info["files"]["E4_models.json"]["sha256"] == add["models_sha256"]
    bad = [k for k, ok in need.items() if not ok]
    if bad:
        raise SystemExit("site4 work refused (run prepare, leakage, commit-addendum first): " + "; ".join(bad))
    return info


# ============================================================================ embed (L1 thinning)
def official_manifest(a) -> dict:
    man = json.loads((a.embed_dir / "z_rna_manifest.json").read_text())
    msha = vk.sha256_file(a.medians)
    want = {"preprocessing": "official", "seq_len": 2048, "medians_sha256": msha, "normalize_total": 10000.0,
            "pooling": "gene-mean"}
    bad = {k: (man.get(k), v) for k, v in want.items() if man.get(k) != v}
    if bad:
        raise SystemExit(f"official embedding manifest differs from this re-embedding: {bad}")
    return man


def embed_block(model, csr, ids, token_ids, pad_id, factors, device, amp, batch, fraction, seed):
    """Official TEDDY-G re-embedding of the cells ``ids``; with ``fraction`` the recovered integer counts are
    binomially thinned first (per-cell generator). Returns z [n, d], ntokens, per-cell count stats."""
    import torch

    from teddy_mm.teddy_encoder import official_values, pool_hidden, rank_encode_official

    data, indices, indptr = csr
    n = len(ids)
    n_genes = token_ids.shape[0]
    L = 2048
    tokens = np.full((n, L), pad_id, dtype=np.int64)
    attn = np.zeros((n, L), dtype=np.int64)
    dev_max = np.zeros(n, dtype=np.float64)
    tot_before = np.zeros(n, dtype=np.int64)
    tot_after = np.zeros(n, dtype=np.int64)
    for s in range(0, n, 256):
        sub = ids[s:s + 256]
        block = np.zeros((len(sub), n_genes), dtype=np.float32)
        for r, cid in enumerate(sub):
            lo, hi = int(indptr[cid]), int(indptr[cid + 1])
            vals, cols = data[lo:hi], indices[lo:hi]
            if fraction is None:
                block[r, cols] = vals
                continue
            counts, dv, _ = e4.recover_counts(vals)
            thin = e4.thin_counts(counts, fraction, seed, int(cid))
            dev_max[s + r] = dv
            tot_before[s + r] = int(counts.sum())
            tot_after[s + r] = int(thin.sum())
            block[r, cols] = thin.astype(np.float32)
        v = official_values(block, factors)
        tokens[s:s + len(sub)], attn[s:s + len(sub)] = rank_encode_official(v, token_ids, max_len=L, pad_id=pad_id, pad_to=L)
    ntok = attn.sum(1).astype(np.int32)
    z = np.zeros((n, model.d_model), dtype=np.float32)
    live = np.where(ntok > 0)[0]  # a cell with no counts left keeps the zero embedding (declared)
    order = live[np.argsort(ntok[live], kind="stable")]
    with torch.no_grad():
        for s in range(0, order.size, batch):
            rows = order[s:s + batch]
            w = int(max(1, ntok[rows].max()))
            idt = torch.from_numpy(tokens[rows, :w]).to(device)
            mk = torch.from_numpy(attn[rows, :w]).to(device)
            with (torch.autocast(device_type=device.type, dtype=amp) if amp is not None else contextlib.nullcontext()):
                h = model.hidden_states(idt, mk)
            z[rows] = pool_hidden(h.float(), mk, "gene-mean").float().cpu().numpy()
    return z, ntok, {"dev_max": dev_max, "tot_before": tot_before, "tot_after": tot_after}


def stage_embed(a, reg) -> None:
    import torch

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import load_gene_medians, load_pad_id, load_teddy, load_vocab, median_factors

    torch.set_num_threads(a.threads)
    if a.pool == "test":
        require_committed_addendum(a)  # before any site4 row is read
    npz = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
    meta = {k: npz[k].astype(str) for k in ("split", "sites", "donors")}  # metadata only: no labels, no protein
    ids = pool_ids(meta, reg, a.pool) if not a.smoke else smoke_val_ids(meta, reg, a.smoke_n_val)
    man = official_manifest(a)
    amp = {"fp16": torch.float16, "bf16": torch.bfloat16}.get(man.get("autocast"))
    out = a.out_dir / "embed" / a.pool
    out.mkdir(parents=True, exist_ok=True)
    seed = int(reg["seeds"]["e4_noise"])
    cfg = {"script": SCRIPT_VERSION, "pool": a.pool, "cells_sha256": vk.index_hash(ids), "n_cells": int(ids.size),
           "fraction": THIN_FRACTION, "seed": seed, "rng": "numpy default_rng([seed, 1, cell_id]).binomial(count, fraction)",
           "preprocessing": man["preprocessing"], "seq_len": man["seq_len"], "autocast": man.get("autocast"),
           "medians_sha256": man["medians_sha256"], "batch_size": a.batch_size, "length_buckets": True,
           "pooling": "gene-mean", "ckpt": a.ckpt.name, "shard": SHARD, "smoke": a.smoke}
    cfg_hash = sha_bytes(json.dumps(cfg, sort_keys=True).encode())[:16]
    device = resolve_device(a.device)
    say(f"embed {a.pool}: {ids.size} cells, thinning {THIN_FRACTION} (seed {seed}), device {device}, cfg {cfg_hash}")
    z = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
    csr = (z["rna_data"], z["rna_indices"], z["rna_indptr"])
    token_ids = z["token_ids"]
    factors = median_factors(z["rna_names"], load_gene_medians(a.medians))
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    model = None
    n_sh = int(math.ceil(ids.size / SHARD))
    t_start = time.time()
    done_now = 0
    for k in range(n_sh):
        f = out / f"shard_{k:04d}.npz"
        if f.exists():
            with np.load(f) as old:
                if str(old["cfg_hash"]) != cfg_hash:
                    raise SystemExit(f"{f} was written under another configuration; use a fresh --out-dir")
            continue
        if model is None:
            model = load_teddy(a.ckpt, device)
        sub = ids[k * SHARD:(k + 1) * SHARD]
        t0 = time.time()
        zz, ntok, st = embed_block(model, csr, sub, token_ids, pad_id, factors, device, amp, a.batch_size, THIN_FRACTION, seed)
        atomic_savez(f, idx=sub, z=zz, ntokens=ntok, cfg_hash=np.array(cfg_hash), runtime_sec=np.array(time.time() - t0), **st)
        done_now += 1
        el = time.time() - t_start
        left = sum(1 for j in range(n_sh) if not (out / f"shard_{j:04d}.npz").exists())
        say(f"embed {a.pool} shard {k + 1}/{n_sh}: {sub.size} cells in {time.time() - t0:.0f}s "
            f"({1000 * (time.time() - t0) / max(sub.size, 1):.0f} ms/cell); ETA {el / done_now * left / 60:.1f} min")
    # fp16 floor: the first cells of the pool re-embedded without thinning vs the stored official z
    fc = out / "check_unthinned.npz"
    if not fc.exists():
        if model is None:
            model = load_teddy(a.ckpt, device)
        sub = ids[:min(FLOOR_CHECK_CELLS, ids.size)]
        zz, _, _ = embed_block(model, csr, sub, token_ids, pad_id, factors, device, amp, a.batch_size, None, seed)
        zref = np.asarray(np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")[sub], dtype=np.float32)
        cos = (zz * zref).sum(1) / (np.linalg.norm(zz, axis=1) * np.linalg.norm(zref, axis=1) + 1e-12)
        atomic_savez(fc, idx=sub, cos=cos, cfg_hash=np.array(cfg_hash))
    with np.load(fc) as cf:
        cos = cf["cos"]
    parts = [np.load(out / f"shard_{k:04d}.npz") for k in range(n_sh)]
    idx_all = np.concatenate([p["idx"] for p in parts])
    assert np.array_equal(idx_all, ids), "shard order mismatch"
    zall = np.concatenate([p["z"] for p in parts]).astype(np.float32)
    ntok = np.concatenate([p["ntokens"] for p in parts])
    dev_max = np.concatenate([p["dev_max"] for p in parts])
    tb = np.concatenate([p["tot_before"] for p in parts])
    ta = np.concatenate([p["tot_after"] for p in parts])
    np.save(out / "z_L1.npy", zall)
    np.save(out / "cells.npy", ids.astype(np.int64))
    stats = {**cfg, "config_hash": cfg_hash, "z_L1_sha256": sha_bytes(zall.tobytes()),
             "count_recovery_max_dev": float(dev_max.max()) if dev_max.size else 0.0,
             "n_cells_recovery_dev_gt_1e-3": int((dev_max > 1e-3).sum()),
             "median_counts_before": float(np.median(tb)), "median_counts_after": float(np.median(ta)),
             "median_kept_fraction": float(np.median(ta / np.maximum(tb, 1))),
             "n_cells_zero_tokens": int((ntok == 0).sum()), "median_ntokens": float(np.median(ntok)),
             "unthinned_check": {"n": int(cos.size), "cos_min": float(cos.min()), "cos_median": float(np.median(cos))},
             "runtime_sec": float(sum(float(p["runtime_sec"]) for p in parts))}
    write_json(out / "embed_manifest.json", stats)
    say(f"embed {a.pool} done: {json.dumps({k: stats[k] for k in ('n_cells', 'n_cells_zero_tokens', 'n_cells_recovery_dev_gt_1e-3', 'median_kept_fraction', 'unthinned_check', 'runtime_sec')})}")


def load_l1(a, pool: str) -> tuple[np.ndarray, np.ndarray, dict]:
    out = a.out_dir / "embed" / pool
    if not (out / "z_L1.npy").exists():
        raise SystemExit(f"{out}/z_L1.npy missing: run --stage embed --pool {pool} first")
    return np.load(out / "cells.npy"), np.load(out / "z_L1.npy"), json.loads((out / "embed_manifest.json").read_text())


# ============================================================================ shared method machinery
def softmax_proba(X: np.ndarray, coef, intercept) -> np.ndarray:
    s = np.asarray(X, dtype=np.float64) @ np.asarray(coef, dtype=np.float64).T + np.asarray(intercept, dtype=np.float64)[None, :]
    s -= s.max(axis=1, keepdims=True)
    p = np.exp(s)
    return p / p.sum(axis=1, keepdims=True)


def lineage_conf_calls(P: np.ndarray, classes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Registered classifier decision: the lineage with the highest probability; confidence = that probability."""
    lin = [classes.index(k) for k in vk.LINEAGES]
    Pl = P[:, lin]
    return Pl.max(axis=1), np.asarray(vk.LINEAGES, dtype=object)[Pl.argmax(axis=1)].astype(str)


def fit_lr(X, y, C):
    from sklearn.linear_model import LogisticRegression

    return LogisticRegression(C=float(C), max_iter=3000).fit(X, y)


def log_loss(y, P, classes):
    from sklearn.metrics import log_loss as ll

    return float(ll(y, P, labels=list(classes)))


def method_scores(level: str, e1: np.ndarray, e2: np.ndarray, reg: dict, comp: dict, models: dict,
                  ffr=None, schema_base=None, run_engine: bool = True) -> dict:
    """Confidence and call of every method at one level for one set of cells."""
    S1, S2 = e4.class_scores12(e1, reg), e4.class_scores12(e2, reg)
    tr = comp["trust"][level]
    sc = comp["source_scale"][level]
    out: dict[str, dict[str, np.ndarray]] = {}

    def from_scores(S):
        return {"conf": S.max(axis=1), "call": e4.calls_from_scores(S), "scores": S}

    out["F0"] = from_scores(S1)
    out["F0b"] = from_scores(S2)
    out["F1"] = from_scores(e4.fuse_average(S1, S2))
    out["F2"] = from_scores(e4.fuse_trust_weighted(S1, S2, tr[e4.CHANNELS[0]], tr[e4.CHANNELS[1]]))
    m3 = models["F3"][level]
    c3, k3 = lineage_conf_calls(softmax_proba(np.c_[S1, S2], m3["coef"], m3["intercept"]), m3["classes"])
    out["F3"] = {"conf": c3, "call": k3}
    m4 = models["F4"]
    c4, k4 = lineage_conf_calls(softmax_proba(np.c_[e1, e2], m4["coef"], m4["intercept"]), m4["classes"])
    out["F4"] = {"conf": c4, "call": k4}
    s1, s2 = sc[e4.CHANNELS[0]], sc[e4.CHANNELS[1]]
    checks = {}
    for name, v in ANM_VARIANTS.items():
        closed = e4.anm_fusion_closed(e1, e2, s1, s2, reg, t2=v["t2"], contradiction=v["contradiction"], propagate=v["propagate"])
        if run_engine:
            bar = comp["bars"][name][level] if comp.get("bars") else 0.0
            schema = e4.anm_schema(schema_base, reg, bar)
            A = np.zeros_like(closed)
            rec = np.empty(e1.shape[0], dtype=object)
            n_rej = 0
            for i in range(e1.shape[0]):
                ev = e4.fusion_events(e1[i], e2[i], s1, s2, reg, t2=v["t2"], contradiction=v["contradiction"])
                A[i], rec[i], nr = e4.anm_fusion_engine_cell(ffr, schema, ev, propagate=v["propagate"])
                n_rej += nr
            calls_bar = e4.calls_from_scores(A, bar)
            rec_s = np.asarray(["" if r is None else str(r) for r in rec])
            checks[name] = {"max_abs_engine_minus_closed": float(np.abs(A - closed).max()) if A.size else 0.0,
                            "n_cells_argmax_differs": int(np.sum(e4.top_class(A) != e4.top_class(closed))),
                            "n_rejected_events": int(n_rej),
                            "n_cells_engine_readout_differs_from_bar_rule": int(np.sum(rec_s != calls_bar)) if comp.get("bars") else None}
        else:
            A = closed
        out[name] = {"conf": A.max(axis=1), "call": e4.calls_from_scores(A), "scores": A}
    out["_checks"] = checks
    out["_disagree"] = e4.top_class(S1) != e4.top_class(S2)
    return out


def correct_of(call: np.ndarray, key: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    key = np.asarray(key).astype(str)
    scored = key != vk.UNSCORED
    return (scored & (np.asarray(call).astype(str) == key)), scored


def sa_at(conf, call, key, ids, coverage, amend) -> float | None:
    sel = e4.select_at_coverage(conf, ids, coverage, amend)
    cor, scored = correct_of(call, key)
    ns = int((sel & scored).sum())
    return (int((sel & cor).sum()) / ns) if ns else None


# ============================================================================ prepare (train / val only)
def prepare_core(full: dict, z0: np.ndarray, z1_ids: np.ndarray, z1: np.ndarray, reg: dict, amend: dict,
                 head_ckpt: Path, ffr, schema_base: dict, smoke: dict | None = None) -> tuple[dict, dict]:
    """Everything the addendum fixes, from training and val rows only. ``full`` holds every row as read from
    disk; the first step keeps the non-test rows and nothing below reads another row."""
    split = np.asarray(full["split"]).astype(str)
    nt = np.where(split != "test")[0]
    adt = np.asarray(full["adt"], dtype=np.float64)[nt]
    ct = np.asarray(full["cell_types"])[nt].astype(str)
    zz = np.asarray(z0[nt], dtype=np.float32)
    sp = split[nt]
    gid = nt.astype(np.int64)
    names = [str(x) for x in full["adt_names"]]
    tr = np.where(sp == "train")[0]
    vl = np.where(sp == "val")[0]
    if smoke is not None:
        vl = vl[np.isin(gid[vl], smoke["val_ids"])]
        if smoke.get("n_train"):
            tr = np.sort(np.random.default_rng(0).choice(tr, size=min(int(smoke["n_train"]), tr.size), replace=False))
    flat = e4.panel_flat(reg)
    pj = [names.index(p) for p in flat]
    seed = int(reg["seeds"]["e4_noise"])
    comp: dict = {"n_train": int(tr.size), "n_val": int(vl.size)}
    models: dict = {}

    # ---- E4 key (measured panel evidence)
    q95m = [r6(x) for x in np.percentile(adt[tr][:, pj], 95, axis=0)]
    m_tr = e4.measured_evidence(adt[tr], names, flat, q95m)
    m_vl = e4.measured_evidence(adt[vl], names, flat, q95m)
    annot_vl = vk.annotation_class(ct[vl], reg)
    tau, tau_rows = e4.choose_tau(m_vl, annot_vl, reg)
    key_tr, key_vl = e4.e4_key(m_tr, reg, tau), e4.e4_key(m_vl, reg, tau)
    comp["key"] = {"measured_q95_train": dict(zip(flat, q95m)), "tau_K": tau, "tau_grid_val_kappa": tau_rows,
                   "val_counts": {c: int((key_vl == c).sum()) for c in vk.CLASSES},
                   "train_counts": {c: int((key_tr == c).sum()) for c in vk.CLASSES},
                   "train_kappa_vs_annotation": r6(vk.cohen_kappa(key_tr, vk.annotation_class(ct[tr], reg)))}
    sites_nt = np.asarray(full["sites"]).astype(str)[nt]
    donors_nt = np.asarray(full["donors"]).astype(str)[nt]
    bk = []
    for rows, kk, nm in ((tr, key_tr, "train"), (vl, key_vl, "val")):
        bt = np.asarray([f"{a}|{b}" for a, b in zip(sites_nt[rows], donors_nt[rows])])
        an = vk.annotation_class(ct[rows], reg)
        for b in sorted(set(bt.tolist())):
            mb = bt == b
            bk.append({"split": nm, "batch": b, "n": int(mb.sum()), "kappa_vs_annotation": r6(vk.cohen_kappa(kk[mb], an[mb]))})
    comp["key"]["batch_kappa"] = bk

    # ---- channel 1 (registered evidence) at L0 for train and val; L1 for val from the thinned re-embedding
    e1_tr = vk.evidence(vk.head_predict(zz[tr], head_ckpt), reg)[:, pj]
    e1_vl = vk.evidence(vk.head_predict(zz[vl], head_ckpt), reg)[:, pj]
    pos = {int(c): i for i, c in enumerate(np.asarray(z1_ids, dtype=np.int64))}
    miss = [int(c) for c in gid[vl] if int(c) not in pos]
    if miss:
        raise SystemExit(f"{len(miss)} val cells have no L1 re-embedding (run --stage embed --pool val)")
    z1_vl = np.asarray(z1, dtype=np.float32)[[pos[int(c)] for c in gid[vl]]]
    e1_vl_L1 = vk.evidence(vk.head_predict(z1_vl, head_ckpt), reg)[:, pj]

    # ---- channel 2 (A1.2 inputs) -> ridge (alpha on val)
    w_tr, cols = va.channel2_inputs(adt[tr], names, flat)
    w_vl, _ = va.channel2_inputs(adt[vl], names, flat)
    alpha, alpha_rows, ridge = e4.ridge_fit_select(w_tr, m_tr, w_vl, m_vl)
    coef, icpt = np.asarray(ridge.coef_, dtype=np.float64), np.asarray(ridge.intercept_, dtype=np.float64)
    p_tr = e4.linear_predict(w_tr, coef, icpt)
    q95c2 = [r6(x) for x in np.percentile(p_tr, 95, axis=0)]
    e2_tr = e4.normalised_evidence(p_tr, q95c2)
    e2_vl = e4.normalised_evidence(e4.linear_predict(w_vl, coef, icpt), q95c2)
    w_vl_drop = e4.dropout_inputs(w_vl, gid[vl], DROPOUT_FRACTION, seed)
    e2_vl_L2 = e4.normalised_evidence(e4.linear_predict(w_vl_drop, coef, icpt), q95c2)
    models["channel2_ridge"] = {"inputs": cols, "targets": flat, "alpha": alpha, "coef": coef, "intercept": icpt,
                                "q95_train_pred": dict(zip(flat, q95c2))}
    comp["channel2"] = {"alpha": alpha, "alpha_grid_val_mean_r2": alpha_rows, "q95_train_pred": dict(zip(flat, q95c2)),
                        "n_inputs": len(cols),
                        "val_dropout_mean_inputs_zeroed": r6(np.mean(np.sum((w_vl_drop == 0) & (w_vl != 0), axis=1)))}

    E1 = {"L0": e1_vl, "L1": e1_vl_L1, "L2": e1_vl, "L3": e1_vl_L1}
    E2 = {"L0": e2_vl, "L1": e2_vl, "L2": e2_vl_L2, "L3": e2_vl_L2}

    # ---- trust and source scales per level (val)
    comp["trust"], comp["trust_per_protein"], comp["source_scale"] = {}, {}, {}
    for L in e4.LEVELS:
        t1, r1 = e4.channel_trust(E1[L], m_vl)
        t2, r2 = e4.channel_trust(E2[L], m_vl)
        comp["trust"][L] = {e4.CHANNELS[0]: r6(t1), e4.CHANNELS[1]: r6(t2)}
        comp["trust_per_protein"][L] = {e4.CHANNELS[0]: dict(zip(flat, map(r6, r1))), e4.CHANNELS[1]: dict(zip(flat, map(r6, r2)))}
        comp["source_scale"][L] = e4.source_scales(comp["trust"][L])
    gaps = {L: r6(abs(comp["trust"][L][e4.CHANNELS[0]] - comp["trust"][L][e4.CHANNELS[1]])) for L in e4.LEVELS}
    comp["trust_gap"] = gaps
    comp["trust_untestable"] = bool(all(g < 0.01 for g in gaps.values()))

    # ---- F4: training cells at L0, C by val log-loss at L0
    X4_tr, X4_vl = np.c_[e1_tr, e2_tr], np.c_[e1_vl, e2_vl]
    rows4, fits4 = [], {}
    for C in C_GRID:
        m = fit_lr(X4_tr, key_tr, C)
        fits4[C] = m
        rows4.append({"C": C, "val_log_loss": r6(log_loss(key_vl, m.predict_proba(X4_vl), m.classes_))})
    C4 = sorted(rows4, key=lambda r: (r["val_log_loss"], r["C"]))[0]["C"]
    m4 = fits4[C4]
    models["F4"] = {"classes": [str(c) for c in m4.classes_], "C": C4, "coef": m4.coef_, "intercept": m4.intercept_,
                    "features": [f"{e4.CHANNELS[0]}:{p}" for p in flat] + [f"{e4.CHANNELS[1]}:{p}" for p in flat]}
    comp["F4"] = {"C": C4, "val_log_loss": rows4, "n_train": int(tr.size),
                  "max_abs_proba_check": float(np.abs(softmax_proba(X4_vl, m4.coef_, m4.intercept_) - m4.predict_proba(X4_vl)).max())}

    # ---- F3: per level, val, C by 5-fold cross-fitted log-loss; cross-fitted predictions for the val endpoint
    from sklearn.model_selection import StratifiedKFold

    folds = list(StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed).split(np.zeros(vl.size), key_vl))
    models["F3"], comp["F3"] = {}, {}
    f3_cross = {}
    for L in e4.LEVELS:
        X3 = np.c_[e4.class_scores12(E1[L], reg), e4.class_scores12(E2[L], reg)]
        rows3, cross = [], {}
        for C in C_GRID:
            P = np.zeros((vl.size, 5))
            classes = None
            for trn, tst in folds:
                m = fit_lr(X3[trn], key_vl[trn], C)
                classes = [str(c) for c in m.classes_]
                if len(classes) != 5:
                    raise SystemExit(f"F3 fold has {len(classes)} classes; val too small")
                P[tst] = m.predict_proba(X3[tst])
            cross[C] = (P, classes)
            rows3.append({"C": C, "cv_log_loss": r6(log_loss(key_vl, P, classes))})
        C3 = sorted(rows3, key=lambda r: (r["cv_log_loss"], r["C"]))[0]["C"]
        mfin = fit_lr(X3, key_vl, C3)
        models["F3"][L] = {"classes": [str(c) for c in mfin.classes_], "C": C3, "coef": mfin.coef_, "intercept": mfin.intercept_,
                           "features": [f"{e4.CHANNELS[0]}:S_{k}" for k in vk.LINEAGES] + [f"{e4.CHANNELS[1]}:S_{k}" for k in vk.LINEAGES]}
        comp["F3"][L] = {"C": C3, "cv_log_loss": rows3}
        f3_cross[L] = lineage_conf_calls(*cross[C3])

    # ---- every method on val per level; bars; val endpoint; comparator
    comp_nobars = {"trust": comp["trust"], "source_scale": comp["source_scale"]}
    val_scores = {}
    checks = {}
    for L in e4.LEVELS:
        ms = method_scores(L, E1[L], E2[L], reg, comp_nobars, models, ffr, schema_base, run_engine=True)
        ms["F3"] = {"conf": f3_cross[L][0], "call": f3_cross[L][1]}  # cross-fitted (A1.9)
        checks[L] = ms.pop("_checks")
        ms.pop("_disagree")
        val_scores[L] = ms
    comp["bars"] = {m: {L: e4.val_bar(val_scores[L][m]["conf"], NO_CALL_RATE) for L in e4.LEVELS} for m in METHODS}
    comp["anm_engine_check_val"] = {L: {k: {"max_abs_engine_minus_closed_lt_1e-9": v["max_abs_engine_minus_closed"] < 1e-9,
                                             "n_cells_argmax_differs": v["n_cells_argmax_differs"],
                                             "n_rejected_events": v["n_rejected_events"]} for k, v in checks[L].items()}
                                    for L in e4.LEVELS}
    vend = {m: {} for m in METHODS}
    for m in METHODS:
        for c in COVERAGES:
            per = {L: sa_at(val_scores[L][m]["conf"], val_scores[L][m]["call"], key_vl, gid[vl], c, amend) for L in e4.LEVELS}
            if any(v is None for v in per.values()):
                raise SystemExit(f"val endpoint undefined for {m} at {c}")
            vend[m][str(c)] = {"per_level": {L: r6(v) for L, v in per.items()}, "mean_L0_L3": r6(np.mean(list(per.values())))}
    comp["val_endpoint"] = vend
    cand = sorted(COMPARATOR_CANDIDATES, key=lambda m: (-vend[m]["0.8"]["mean_L0_L3"], COMPARATOR_CANDIDATES.index(m)))
    comp["best_comparator"] = cand[0]
    comp["comparator_ranking_val"] = [{"method": m, "val_mean_sa_0.8": vend[m]["0.8"]["mean_L0_L3"]} for m in cand]
    # sanity: channel 1 at L0 is the registered Q1 rule, so on the full val split its bar is the registered Q1 bar
    comp["F0_L0_bar_vs_registered_Q1_bar"] = {"F0_L0_bar": comp["bars"]["F0"]["L0"],
                                              "registered_Q1_bar": float(reg["questions"]["Q1"]["bar"]),
                                              "applies": smoke is None}
    comp["val_channel_disagreement_share"] = {L: r6(np.mean(e4.top_class(e4.class_scores12(E1[L], reg)) != e4.top_class(e4.class_scores12(E2[L], reg)))) for L in e4.LEVELS}
    return comp, models


ADDENDUM_TEXT = {
    "addendum_to": "registration/registration_v3.json experiments.E4 as amended by A1 (A1.1, A1.2, A1.9); the registered "
                   "key rule, channels, noise levels, trust rule, methods, controls, endpoint, margin and outcomes are not "
                   "changed here",
    "experiment": "E4",
    "addendum_version": 1,
    "fixed_before_site4": "computed by bridge_anm/v3_e4_fusion.py --stage prepare from training (split == train) and val "
                          "(donor 18303) rows only; checked by --stage leakage (site4-poisoned rebuild identical, "
                          "val-poisoned rebuild different) and committed with --stage commit-addendum before any site4 "
                          "row (RNA, protein, labels, embedding) is read by this script",
    "open_choices_fixed": {
        "key": "E4 key as registered: measured evidence = stored ADT (per-cell CLR, as in cite_arrays.npz) / training q95 "
               "of the stored value (6 decimals), clipped to [0, 1]; class score = equal-weight mean over the primary "
               "panel; argmax ties go to the first class in B, T, NK, myeloid order; Cohen kappa over B, T, NK, myeloid, "
               "OUT against the registration annotation map; tau_K = the grid value with the highest val kappa "
               "(6 decimals), ties to the smaller tau; the E4 key has no unscored cells",
        "channel1": "the registered evidence (official z_rna.npy -> phase-1 head, size factor 1, registered q95_train_pred, "
                    "clipped) of the 12 primary-panel proteins",
        "channel2": "A1.2 inputs (v3_amend.channel2_inputs, 122 values); sklearn Ridge(alpha, fit_intercept=True), no "
                    "input scaling, one alpha for the 12 targets; target = the measured panel evidence of the key; "
                    "alpha = the grid value with the highest mean val R2 over the 12 targets (6 decimals), ties to the "
                    "smaller alpha; refit on all training cells (coefficients in E4_models.json); evidence = prediction / "
                    "training q95 of its own predictions (6 decimals), clipped to [0, 1] like channel 1 (ANM event "
                    "values must lie in [0, 1])",
        "L1_thinning": "integer counts = stored RNA row / its smallest non-zero value, rounded (the largest relative "
                       "deviation from an integer is reported); Binomial(count, 0.2) per gene with numpy "
                       "default_rng([seeds.e4_noise, 1, cell_id]); re-embedded with the settings of "
                       "z_rna_manifest.json (official preprocessing, gene medians, torch.topk 2048, fp16 autocast, "
                       "length-bucketed batches of 32, gene-mean); then the frozen head and the registered evidence "
                       "normaliser; a cell with no counts left keeps the zero embedding (counted)",
        "L2_dropout": "per cell, 61 of the 122 channel-2 inputs (after the A1.2 transform) set to 0: numpy "
                      "default_rng([seeds.e4_noise, 2, cell_id]).permutation(122)[:61]; the ridge and its normaliser are "
                      "the clean ones",
        "trust": "Pearson over val cells at the level between the channel evidence (clipped) and the measured panel "
                 "evidence; a constant column gives r = 0; mean of max(0, r) over the 12 proteins, 6 decimals; ANM "
                 "source scale = trust / max of the two channels' trust at the level, 6 decimals",
        "untestable": "the trust mechanism is reported untestable if |trust_1 - trust_2| < 0.01 (6 decimals) at every "
                      "level; the endpoints are still computed and reported",
        "F1_F2": "class scores = equal-weight panel means per channel; F1 = their mean; F2 = (t1 S1 + t2 S2) / (t1 + t2) "
                 "with the level's val trust",
        "F3": "sklearn LogisticRegression(lbfgs, max_iter 3000, multinomial) on the 8 class scores [S1, S2] of the val "
              "cells at the level, labels = E4 key on val (5 classes); C from [0.01, 0.1, 1, 10, 100] by the pooled "
              "log-loss of 5-fold cross-fitted probabilities (StratifiedKFold, shuffle, random_state seeds.e4_noise), "
              "ties to the smaller C; the val endpoint (A1.9) uses those cross-fitted predictions; the test model is "
              "the refit on all val cells at the level",
        "F4": "the same model on the 24 evidence values [e1, e2] of the training cells at L0, labels = E4 key on training "
              "cells; C from the same grid by val log-loss at L0 (the registered classifier procedure); applied "
              "unchanged at every level",
        "F3_F4_decision": "registered classifier decision: the lineage with the highest probability; confidence = that "
                          "probability",
        "F5": "events as documented in bridge_anm/lib/v3_e4.py: support events s_c * e_c[p] for every channel, class and "
              "panel protein; for each channel, its top class (first in B, T, NK, myeloid on ties) contradicts every "
              "other class with value s_c * 0.5 * top score; channel 1 at t = 0, channel 2 at t = 0; ANM default field "
              "(registered: steps 4, retention 0.82, diffusion 0.16, schema source_scale 1.0); action readout; "
              "readout_threshold = the F5 bar of the level. ANM's finite_graph_scalar has one schema source_scale, so "
              "the per-channel source scale multiplies the event value (identical, since the injection is source_scale "
              "* amplitude and amplitude = +/- value)",
        "controls": "C1 = F5 with ANM disable_propagation (every action stays 0, so its selection is the A1.1 "
                    "permutation); C2 = F5 with channel 2 at t = 2; C3 = F5 without contradiction events (proportional "
                    "to F2 cell by cell); closed form = v3_e4.anm_fusion_closed, checked against the engine on every cell",
        "bars": "each method's operating bar = val quantile at no-call 0.15 of its confidence at the level, 6 decimals; "
                "reported with realised site4 coverage only (the endpoint is matched coverage)",
        "confidence": "F0, F0b, F1, F2: top class score; F3, F4: highest lineage probability; F5, C1-C3: top ANM action "
                      "score; rounded to 12 decimals before ranking; ties by the A1.1 permutation; at coverage c a method "
                      "calls its ceil(c N) most confident cells with its argmax lineage",
        "endpoint": "selective accuracy on the E4 key at matched coverage 0.8, mean over L0-L3, test_primary; coverage 0.9 "
                    "is reported with the same rule and never changes the verdict; the differences F5 - comparator and "
                    "F5 - F2 are computed inside each bootstrap replicate on the same resampled cells for every method "
                    "and level",
        "comparator": "A1.9: the F1-F4 method with the highest val mean selective accuracy at 0.8 (F3 cross-fitted); "
                      "ties to the earlier of F1, F2, F3, F4",
        "bootstrap": "two-stage (donors with replacement, then cells with replacement within each drawn donor), B 2000, "
                     "seed 1, one numpy generator, donors in sorted order; selection recomputed in every replicate (A1.1); "
                     "percentile 95% interval; per-donor points (and a within-donor cell bootstrap, descriptive)",
        "decision": "registered: win = point >= 0.01, lower bound > 0, and >= 0.005 in each primary donor; loss = mirror; "
                    "equivalent = interval inside (-0.01, 0.01); else inconclusive",
        "outcomes": "primary (F5 vs comparator): win -> 'field adds value', loss -> 'field hurts', equivalent -> 'F5 "
                    "equivalent to the best non-ANM fusion', inconclusive; field effect (F5 vs F2): equivalent -> 'field "
                    "adds nothing' (acceptable), win -> 'contradiction events and field gain add over the trust-weighted "
                    "rule', loss -> 'they hurt', inconclusive",
        "key_validity_reporting": "reporting rule only (mirrors registration section 4.3 for the v3 key; no E4 key-validity "
                                  "rule is registered): the E4 key's kappa against the annotation is reported per site4 "
                                  "split before any method result; if it is below 0.85 on test_primary, the v3 primary-key "
                                  "results are reported with equal prominence beside the E4-key verdict. The key is not "
                                  "switched and the verdict is not changed. Why: with one training q95 normaliser and one "
                                  "tau_K, the E4 key agrees well with the annotation at site1 (val 0.91) but poorly at site3 "
                                  "(training batches 0.29-0.58; computed.key.batch_kappa)",
        "secondary": "the same endpoints on the v3 primary key (unscored cells count for coverage, not accuracy; partly "
                     "circular for channel 2, disclosed) and on test_secondary (donor 15078); never change the verdict",
    },
}


def stage_prepare(a, reg, amend) -> None:
    ffr = e4.import_anm(a.anm_root)
    schema_base = json.loads((ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    t0 = time.time()
    full = load_full(a.processed)
    z0 = np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")
    z1_ids, z1, emb_man = load_l1(a, "val")
    smoke = None
    if a.smoke:
        ids = smoke_val_ids(full, reg, a.smoke_n_val)
        smoke = {"val_ids": ids[: ids.size // 2], "n_train": a.smoke_n_train}
    say(f"prepare: data loaded ({time.time() - t0:.0f}s); computing the addendum core from train/val rows")
    comp, models = prepare_core(full, z0, z1_ids, z1, reg, amend, a.head_ckpt, ffr, schema_base, smoke)
    models_bytes = e4.json_dump(e4.jsonable(models)).encode()
    core = {"computed": comp, "models_sha256": sha_bytes(models_bytes)}
    core_bytes = e4.json_dump(e4.jsonable(core)).encode()
    st = addendum_status(a)
    add = {**ADDENDUM_TEXT, "computed": comp, "models_file": "E4_models.json", "models_sha256": core["models_sha256"],
           "core_sha256": sha_bytes(core_bytes),
           "registration_sha256": st["registration_sha256"], "amendment_A1_sha256": st["amendment_A1_sha256"],
           "script": "bridge_anm/v3_e4_fusion.py", "script_version": SCRIPT_VERSION,
           "script_sha256": vk.sha256_file(Path(__file__)), "lib_sha256": vk.sha256_file(ROOT / "bridge_anm/lib/v3_e4.py"),
           "val_L1_embedding": {k: emb_man.get(k) for k in ("config_hash", "z_L1_sha256", "n_cells", "n_cells_zero_tokens",
                                                           "n_cells_recovery_dev_gt_1e-3", "median_kept_fraction",
                                                           "unthinned_check")},
           "smoke": a.smoke}
    d = addenda_dir(a)
    d.mkdir(parents=True, exist_ok=True)
    if not a.smoke and (d / "E4.json").exists() and git_committed(d / "E4.json"):
        old = json.loads((d / "E4.json").read_text())
        if old.get("core_sha256") != add["core_sha256"]:
            raise SystemExit(f"a committed E4 addendum (core {old.get('core_sha256', '')[:12]}) exists and this rebuild gives "
                             f"core {add['core_sha256'][:12]}; a frozen addendum is never overwritten (write an amendment)")
        say("prepare: the committed addendum has the same core; nothing rewritten")
        return
    (d / "E4_models.json").write_bytes(models_bytes)
    (d / "E4.json").write_bytes(e4.json_dump(e4.jsonable(add)).encode())
    update_hash_lines(d, {"E4.json": vk.sha256_file(d / "E4.json"), "E4_models.json": vk.sha256_file(d / "E4_models.json")})
    write_json(a.out_dir / "prepare" / "E4_prepare.json", {"addendum": add})
    say(f"prepare done ({time.time() - t0:.0f}s): tau_K {comp['key']['tau_K']}, alpha {comp['channel2']['alpha']}, "
        f"trust {comp['trust']}, comparator {comp['best_comparator']}, untestable {comp['trust_untestable']}; "
        f"wrote {d / 'E4.json'} (core {add['core_sha256'][:12]})")


def update_hash_lines(d: Path, entries: dict[str, str]) -> None:
    """Replace only this experiment's lines in addenda/HASHES.txt (other builders' lines are kept)."""
    hf = d / "HASHES.txt"
    lines = [ln for ln in (hf.read_text().splitlines() if hf.exists() else []) if not any(ln.rstrip().endswith(f"  {n}") for n in entries)]
    hf.write_text("\n".join(lines + [f"{h}  {n}" for n, h in entries.items()]) + "\n")


# ============================================================================ leakage
def stage_leakage(a, reg, amend) -> None:
    """Rebuild the addendum core (a) from the real inputs, (b) with every site4 row of protein, cell type and
    embedding replaced by random values (must be byte-identical to (a) and to the written addendum), (c) with the
    val rows poisoned the same way, including the val L1 embedding (must differ)."""
    ffr = e4.import_anm(a.anm_root)
    schema_base = json.loads((ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    d = addenda_dir(a)
    add = json.loads((d / "E4.json").read_text())
    full = load_full(a.processed)
    z0 = np.asarray(np.load(a.embed_dir / "z_rna.npy", mmap_mode="r"), dtype=np.float32)
    z1_ids, z1, _ = load_l1(a, "val")
    smoke = None
    if a.smoke:
        ids = smoke_val_ids(full, reg, a.smoke_n_val)
        smoke = {"val_ids": ids[: ids.size // 2], "n_train": a.smoke_n_train}
    types = sorted(reg["annotation_map"])

    def core_sha(f, zz0, zz1):
        comp, models = prepare_core(f, zz0, z1_ids, zz1, reg, amend, a.head_ckpt, ffr, schema_base, smoke)
        mb = e4.json_dump(e4.jsonable(models)).encode()
        return sha_bytes(e4.json_dump(e4.jsonable({"computed": comp, "models_sha256": sha_bytes(mb)})).encode())

    def poisoned(rows: np.ndarray, seed: int, z1_rows: np.ndarray | None = None):
        rng = np.random.default_rng(seed)
        f = dict(full)
        f["adt"] = np.array(full["adt"], copy=True)
        f["adt"][rows] = rng.gamma(1.0, 1.0, size=(rows.size, f["adt"].shape[1])).astype(np.float32)
        f["cell_types"] = np.array(full["cell_types"], copy=True)
        f["cell_types"][rows] = rng.choice(types, size=rows.size)
        zz0 = np.array(z0, copy=True)
        zz0[rows] = rng.normal(size=(rows.size, zz0.shape[1])).astype(np.float32)
        zz1 = np.array(z1, copy=True)
        if z1_rows is not None:
            zz1[z1_rows] = rng.normal(size=(z1_rows.size, zz1.shape[1])).astype(np.float32)
        return f, zz0, zz1

    t0 = time.time()
    real = core_sha(full, z0, z1)
    say(f"leakage: real core {real[:12]} ({time.time() - t0:.0f}s)")
    test_rows = np.where(full["split"] == "test")[0]
    s4 = core_sha(*poisoned(test_rows, 101))
    say(f"leakage: site4-poisoned core {s4[:12]} ({time.time() - t0:.0f}s)")
    val_rows = np.where(full["split"] == "val")[0]
    s_val = core_sha(*poisoned(val_rows, 202, np.arange(z1.shape[0])))
    say(f"leakage: val-poisoned core {s_val[:12]} ({time.time() - t0:.0f}s)")
    rep = {"addendum_core_sha256": add["core_sha256"], "real_inputs_core_sha256": real,
           "site4_poisoned_core_sha256": s4, "val_poisoned_core_sha256": s_val,
           "checks": {"real equals the written addendum core": real == add["core_sha256"],
                      "site4-poisoned equals real (no site4 dependence)": s4 == real,
                      "val-poisoned differs (positive control)": s_val != real},
           "poisoned": "protein (all 134), cell type and official embedding rows; the val control also poisons the val "
                       "L1 re-embedding; site4 RNA is never read by prepare (its only RNA-derived input is the val L1 "
                       "embedding)", "script_sha256": vk.sha256_file(Path(__file__)), "smoke": a.smoke}
    rep["passed"] = all(rep["checks"].values())
    (d / "E4_leakage_check.json").write_bytes(e4.json_dump(rep).encode())
    update_hash_lines(d, {"E4_leakage_check.json": vk.sha256_file(d / "E4_leakage_check.json")})
    write_json(a.out_dir / "leakage" / "E4_leakage_check.json", rep)
    say(f"leakage check {'PASSED' if rep['passed'] else 'FAILED'}: {rep['checks']}")
    if not rep["passed"]:
        raise SystemExit(1)


# ============================================================================ commit the addendum
def stage_commit_addendum(a) -> None:
    """Commit exactly the E4 addendum files and their HASHES.txt lines (a temporary index; other builders'
    staged or uncommitted changes are left untouched)."""
    d = addenda_dir(a)
    missing = [n for n in ADDENDUM_FILES if not (d / n).exists()]
    if missing:
        raise SystemExit(f"missing {', '.join(missing)} in {d}: run --stage prepare and --stage leakage first")
    add = json.loads((d / "E4.json").read_text())
    lk = json.loads((d / "E4_leakage_check.json").read_text())
    if not (lk.get("passed") and lk.get("addendum_core_sha256") == add["core_sha256"]):
        raise SystemExit("leakage check missing, failed or for another addendum core; run --stage leakage")
    if vk.sha256_file(d / "E4_models.json") != add["models_sha256"]:
        raise SystemExit("E4_models.json is not the file the addendum pins")
    if add.get("smoke"):
        raise SystemExit("refusing to commit a smoke addendum")
    rel = {n: str((d / n).resolve().relative_to(ROOT.resolve())) for n in ADDENDUM_FILES}
    hrel = str((d / "HASHES.txt").resolve().relative_to(ROOT.resolve()))
    head = git("rev-parse", "HEAD").strip()
    branch = git("symbolic-ref", "--short", "HEAD").strip()
    try:
        old = git("show", f"{head}:{hrel}")
    except RuntimeError:
        old = ""
    keep = [ln for ln in old.splitlines() if not any(ln.rstrip().endswith(f"  {n}") for n in ADDENDUM_FILES)]
    new_hashes = ("\n".join(keep + [f"{vk.sha256_file(d / n)}  {n}" for n in ADDENDUM_FILES]) + "\n").encode()
    files = {rel[n]: (d / n).read_bytes() for n in ADDENDUM_FILES}
    files[hrel] = new_hashes
    msg = (f"E4 addendum: train/val numbers fixed before site4 (core {add['core_sha256'][:12]}, models "
           f"{add['models_sha256'][:12]}), leakage check passed\n\n"
           f"tau_K {add['computed']['key']['tau_K']}, ridge alpha {add['computed']['channel2']['alpha']}, comparator "
           f"{add['computed']['best_comparator']}, trust untestable {add['computed']['trust_untestable']}.\n\n"
           "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n")
    with tempfile.TemporaryDirectory() as td:
        env = {**os.environ, "GIT_INDEX_FILE": str(Path(td) / "index")}
        git("read-tree", head, env=env)
        for path, data in files.items():
            blob = git("hash-object", "-w", "--stdin", input_bytes=data).strip()
            git("update-index", "--add", "--cacheinfo", f"100644,{blob},{path}", env=env)
        tree = git("write-tree", env=env).strip()
    commit = git("commit-tree", tree, "-p", head, "-m", msg).strip()
    git("update-ref", f"refs/heads/{branch}", commit, head)
    git("reset", "-q", "--", *files.keys())
    say(f"committed the E4 addendum as {commit[:10]} on {branch}: {', '.join(files)}")


# ============================================================================ evaluate (site4)
def stage_evaluate(a, reg, amend) -> None:
    gate = require_committed_addendum(a)
    d = addenda_dir(a)
    add = json.loads((d / "E4.json").read_text())
    models = json.loads((d / "E4_models.json").read_text())
    comp = add["computed"]
    ffr = e4.import_anm(a.anm_root)
    schema_base = json.loads((ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    t0 = time.time()
    full = load_full(a.processed)
    names = [str(x) for x in full["adt_names"]]
    flat = e4.panel_flat(reg)
    pj = [names.index(p) for p in flat]
    seed = int(reg["seeds"]["e4_noise"])
    if a.smoke:  # the held-out half of the smoke val subset plays the test split, as two pseudo-donors
        ids_all = smoke_val_ids(full, reg, a.smoke_n_val)
        ev_ids = ids_all[ids_all.size // 2:]
        splits = {"test_primary": ev_ids}
        donor_of = {int(c): ("pseudoA" if i % 2 == 0 else "pseudoB") for i, c in enumerate(ev_ids)}
        primary_donors = ["pseudoA", "pseudoB"]
        l1_pool = "val"
    else:
        splits = {s: vk.split_indices(full["split"], full["sites"], full["donors"], reg, s) for s in ("test_primary", "test_secondary")}
        ev_ids = np.sort(np.concatenate(list(splits.values())))
        donor_of = {int(c): str(full["donors"][c]) for c in ev_ids}
        primary_donors = sorted(reg["splits"]["test_primary"]["donors"])
        l1_pool = "test"
    z1_ids, z1, emb_man = load_l1(a, l1_pool)
    pos = {int(c): i for i, c in enumerate(z1_ids)}
    if any(int(c) not in pos for c in ev_ids):
        raise SystemExit("some evaluated cells have no L1 re-embedding")
    # ---- keys (verifier only: measured protein and annotation of the evaluated cells)
    adt = np.asarray(full["adt"], dtype=np.float64)[ev_ids]
    ct = full["cell_types"][ev_ids].astype(str)
    m12 = e4.measured_evidence(adt, names, flat, [comp["key"]["measured_q95_train"][p] for p in flat])
    key_e4 = e4.e4_key(m12, reg, comp["key"]["tau_K"])
    keys_v3 = vk.build_keys(ct, adt, names, reg)
    annot = keys_v3["annotation"]
    keys = {"E4": key_e4, "v3_primary": keys_v3["primary"]}
    # ---- channels per level
    z0 = np.asarray(np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")[ev_ids], dtype=np.float32)
    e1_L0 = vk.evidence(vk.head_predict(z0, a.head_ckpt), reg)[:, pj]
    e1_L1 = vk.evidence(vk.head_predict(np.asarray(z1, dtype=np.float32)[[pos[int(c)] for c in ev_ids]], a.head_ckpt), reg)[:, pj]
    cm = models["channel2_ridge"]
    w, cols = va.channel2_inputs(adt, names, flat)
    if cols != cm["inputs"]:
        raise SystemExit("channel-2 input order differs from the addendum model")
    q2 = [cm["q95_train_pred"][p] for p in flat]
    e2_L0 = e4.normalised_evidence(e4.linear_predict(w, cm["coef"], cm["intercept"]), q2)
    e2_L2 = e4.normalised_evidence(e4.linear_predict(e4.dropout_inputs(w, ev_ids, DROPOUT_FRACTION, seed), cm["coef"], cm["intercept"]), q2)
    E1 = {"L0": e1_L0, "L1": e1_L1, "L2": e1_L0, "L3": e1_L1}
    E2 = {"L0": e2_L0, "L1": e2_L0, "L2": e2_L2, "L3": e2_L2}
    say(f"evaluate: {ev_ids.size} cells, channels ready ({time.time() - t0:.0f}s); running every method and the ANM engine")
    scores, checks, disagree = {}, {}, {}
    cache = a.out_dir / "e4_cells_scores.npz"
    add_sha = gate["files"]["E4.json"]["sha256"]
    code_sha = sha_bytes(Path(__file__).read_bytes() + (ROOT / "bridge_anm/lib/v3_e4.py").read_bytes())
    for L in e4.LEVELS:  # resumable: a level whose scores were saved under this addendum is not recomputed
        part = a.out_dir / "levels" / f"{L}.npz"
        if part.exists():
            with np.load(part, allow_pickle=False) as P:
                if (str(P["addendum_sha256"]) == add_sha and str(P["code_sha256"]) == code_sha
                        and np.array_equal(P["cell_ids"], ev_ids)):
                    scores[L] = {m: {"conf": P[f"{m}_conf"], "call": P[f"{m}_call"]} for m in METHODS}
                    checks[L] = json.loads(str(P["checks"]))
                    disagree[L] = P["disagree"]
                    say(f"evaluate {L}: scores loaded from {part}")
                    continue
        tl = time.time()
        ms = method_scores(L, E1[L], E2[L], reg, comp, models, ffr, schema_base, run_engine=True)
        checks[L] = ms.pop("_checks")
        disagree[L] = ms.pop("_disagree")
        scores[L] = ms
        atomic_savez(part, cell_ids=ev_ids, addendum_sha256=np.array(add_sha), code_sha256=np.array(code_sha), checks=np.array(json.dumps(checks[L])),
                     disagree=disagree[L], **{f"{m}_conf": ms[m]["conf"] for m in METHODS},
                     **{f"{m}_call": ms[m]["call"] for m in METHODS})
        say(f"evaluate {L}: done in {time.time() - tl:.0f}s; engine vs closed form {json.dumps(checks[L])}")
    atomic_savez(cache, cell_ids=ev_ids,
                 **{f"{L}_{m}_conf": scores[L][m]["conf"] for L in e4.LEVELS for m in METHODS},
                 **{f"{L}_{m}_call": scores[L][m]["call"] for L in e4.LEVELS for m in METHODS},
                 key_E4=key_e4, key_v3_primary=keys["v3_primary"])
    comparator = comp["best_comparator"]
    R = {"experiment": "E4", "script": SCRIPT_VERSION, "registration_sha256": gate["registration_sha256"],
         "amendment_A1_sha256": gate["amendment_A1_sha256"],
         "addendum_sha256": gate["files"]["E4.json"]["sha256"], "addendum_core_sha256": add["core_sha256"],
         "registration_gate": gate, "smoke": a.smoke, "arm_names": ARM_NAMES, "comparator": comparator,
         "margin": float(reg["experiments"]["E4"]["margin"]), "coverages": list(COVERAGES), "n_boot": a.n_boot,
         "trust_val": comp["trust"], "source_scale": comp["source_scale"], "trust_untestable": comp["trust_untestable"],
         "val_endpoint": comp["val_endpoint"], "L1_embedding": {k: emb_man.get(k) for k in (
             "config_hash", "n_cells", "n_cells_zero_tokens", "n_cells_recovery_dev_gt_1e-3", "median_kept_fraction",
             "unthinned_check", "runtime_sec")}}
    # ---- checks: closed form, F5 without contradiction = F2
    R["checks"] = {"anm_engine_vs_closed_form": checks}
    f2c3 = {}
    for L in e4.LEVELS:
        for c in COVERAGES:
            s2 = e4.select_at_coverage(scores[L]["F2"]["conf"], ev_ids, c, amend)
            s3 = e4.select_at_coverage(scores[L]["C3"]["conf"], ev_ids, c, amend)
            f2c3[f"{L}@{c}"] = {"n_selection_differs": int((s2 != s3).sum()),
                                "n_argmax_call_differs": int((scores[L]["F2"]["call"] != scores[L]["C3"]["call"]).sum())}
    R["checks"]["C3_equals_F2"] = f2c3
    R["checks"]["anm_bridge_failure"] = any(v["max_abs_engine_minus_closed"] >= 1e-9 or v["n_cells_argmax_differs"] or v["n_rejected_events"]
                                            for L in e4.LEVELS for v in checks[L].values())
    # ---- key validity and composition (descriptive)
    R["key"] = {}
    for sname, sids in splits.items():
        m = np.isin(ev_ids, sids)
        R["key"][sname] = {"E4_counts": {c: int((key_e4[m] == c).sum()) for c in vk.CLASSES},
                           "E4_kappa_vs_annotation": r6(vk.cohen_kappa(key_e4[m], annot[m])),
                           "E4_vs_v3_primary_agreement_on_scored": r6(np.mean((key_e4[m] == keys["v3_primary"][m])[keys["v3_primary"][m] != vk.UNSCORED])),
                           "v3_primary_counts": {c: int((keys["v3_primary"][m] == c).sum()) for c in (*vk.CLASSES, vk.UNSCORED)}}
    # ---- endpoints
    R["splits"] = {}
    for sname, sids in splits.items():
        m = np.isin(ev_ids, sids)
        R["splits"][sname] = evaluate_split(scores, keys, ev_ids, m, donor_of, sname, primary_donors, comparator, reg, amend,
                                            a.n_boot, comp["bars"], annot, disagree)
        say(f"evaluate {sname}: bootstrap done ({time.time() - t0:.0f}s)")
    prim = R["splits"]["test_primary"]["E4"]["comparisons"]
    R["verdict"] = verdict_text(prim, comparator, comp["trust_untestable"], R["checks"]["anm_bridge_failure"])
    kap = R["key"]["test_primary"]["E4_kappa_vs_annotation"]
    R["verdict"]["E4_key_kappa_test_primary"] = kap
    R["verdict"]["key_validity_flag"] = bool(kap < 0.85)
    if kap < 0.85:
        sec = verdict_text(R["splits"]["test_primary"]["v3_primary"]["comparisons"], comparator, comp["trust_untestable"],
                           R["checks"]["anm_bridge_failure"])
        R["verdict"]["v3_primary_key_equal_prominence"] = sec
        R["verdict"]["summary"] += (f". E4 key kappa vs annotation on test_primary {kap:.3f} < 0.85, so the v3 primary key is "
                                    f"reported with equal prominence: {sec['summary']}")
    write_json(a.out_dir / "E4_results.json", R)
    write_report(a.out_dir / "REPORT.md", R, add)
    say(f"evaluate done ({time.time() - t0:.0f}s): {R['verdict']['summary']}")


def evaluate_split(scores, keys, ev_ids, m, donor_of, sname, primary_donors, comparator, reg, amend, n_boot, bars,
                   annot, disagree) -> dict:
    ids = ev_ids[m]
    donors = np.asarray([str(donor_of[int(c)]) for c in ids])
    dlist = sorted(str(x) for x in set(donors.tolist()))
    margin = float(reg["experiments"]["E4"]["margin"])
    out = {"n_cells": int(ids.size), "donors": {d: int((donors == d).sum()) for d in dlist}}
    W = e4.two_stage_weights(donors, n_boot, int(reg["seeds"]["bootstrap"]), dlist)
    for kname, key_all in keys.items():
        key = key_all[m]
        res = {"per_method": {}, "comparisons": {}}
        boot = {}  # (method, coverage) -> [B] mean over levels
        point = {}
        per_donor = {}
        for meth in METHODS:
            res["per_method"][meth] = {"per_level": {}}
            accB = {c: np.zeros(W.shape[0]) for c in COVERAGES}
            for L in e4.LEVELS:
                conf, call = scores[L][meth]["conf"][m], scores[L][meth]["call"][m]
                cor, scored = correct_of(call, key)
                order = e4.ranking_order(conf, ids, amend)
                row = {"selective_accuracy": {}, "grid": {}, "per_donor": {}}
                for c in COVERAGES:
                    v = sa_at(conf, call, key, ids, c, amend)
                    row["selective_accuracy"][str(c)] = v
                    for ch in range(0, W.shape[0], 250):
                        accB[c][ch:ch + 250] += e4.weighted_selective_accuracy(order, W[ch:ch + 250], cor, scored, c) / len(e4.LEVELS)
                    for dn in dlist:
                        dm = donors == dn
                        row["per_donor"].setdefault(dn, {})[str(c)] = sa_at(conf[dm], call[dm], key[dm], ids[dm], c, amend)
                for c in reg["experiments"]["common"]["matched_coverage"]["grid"]:
                    row["grid"][str(c)] = sa_at(conf, call, key, ids, c, amend)
                res["per_method"][meth]["per_level"][L] = row
            for c in COVERAGES:
                vals = [res["per_method"][meth]["per_level"][L]["selective_accuracy"][str(c)] for L in e4.LEVELS]
                point[(meth, c)] = float(np.mean(vals))
                boot[(meth, c)] = accB[c]
                per_donor[(meth, c)] = {dn: float(np.mean([res["per_method"][meth]["per_level"][L]["per_donor"][dn][str(c)] for L in e4.LEVELS])) for dn in dlist}
                res["per_method"][meth][f"mean_L0_L3@{c}"] = {"point": point[(meth, c)], "ci95": [float(x) for x in np.nanpercentile(boot[(meth, c)], [2.5, 97.5])],
                                                             "per_donor": per_donor[(meth, c)]}
        pairs = {"F5_vs_comparator": ("F5", comparator), "F5_vs_F2": ("F5", "F2"), "F5_vs_F0": ("F5", "F0"),
                 "F2_vs_F0": ("F2", "F0"), "C2_vs_F5": ("C2", "F5"), "C1_vs_F5": ("C1", "F5"), "F5_vs_F1": ("F5", "F1"),
                 "F5_vs_F3": ("F5", "F3"), "F5_vs_F4": ("F5", "F4")}
        for pname, (x, y) in pairs.items():
            for c in COVERAGES:
                dB = boot[(x, c)] - boot[(y, c)]
                pt = point[(x, c)] - point[(y, c)]
                lo, hi = (float(v) for v in np.nanpercentile(dB, [2.5, 97.5]))
                pdn = {dn: per_donor[(x, c)][dn] - per_donor[(y, c)][dn] for dn in dlist}
                prim_d = [pdn[dn] for dn in primary_donors if dn in pdn]
                res["comparisons"][f"{pname}@{c}"] = {"a": x, "b": y, "point": pt, "ci95": [lo, hi], "per_donor": pdn,
                                                      "decision": e4.decide(pt, lo, hi, prim_d, margin) if sname == "test_primary" else
                                                      e4.decide(pt, lo, hi, [], margin)}
        out[kname] = res
    # deployed bars, call tables, channel disagreement (E4 key, descriptive)
    key = keys["E4"][m]
    ops = {}
    for meth in METHODS:
        ops[meth] = {}
        for L in e4.LEVELS:
            conf, call = scores[L][meth]["conf"][m], scores[L][meth]["call"][m]
            bar = bars[meth][L]
            called = conf >= bar
            cor, scored = correct_of(np.where(called, call, vk.NO_CALL), key)
            ops[meth][L] = {"bar": bar, "coverage": float(called.mean()),
                            "selective_accuracy": (float(cor[called].sum() / max(1, (called & scored).sum())) if called.any() else None),
                            "out_decline_rate": float((~called[key == "OUT"]).mean()) if np.any(key == "OUT") else None}
    out["E4"]["operating_points"] = ops
    tables = {}
    for meth in ("F0", "F0b", "F2", comparator, "F5"):
        conf, call = scores["L0"][meth]["conf"][m], scores["L0"][meth]["call"][m]
        sel = e4.select_at_coverage(conf, ids, 0.8, amend)
        t = {}
        for typ in sorted(set(annot[m])):
            tm = annot[m] == typ
            ct_ = {k: int(np.sum(sel & tm & (call == k))) for k in vk.LINEAGES}
            ct_["no_call"] = int(np.sum(~sel & tm))
            t[typ] = ct_
        rec = {k: float(np.sum(sel & (key == k) & (call == k)) / max(1, int(np.sum(key == k)))) for k in vk.LINEAGES}
        tables[meth] = {"call_table_by_annotated_type_L0@0.8": t, "recall_by_E4_key_class_L0@0.8": rec}
    out["E4"]["call_tables"] = tables
    dis = {}
    for L in e4.LEVELS:
        dm = disagree[L][m]
        row = {"share_channels_disagree": float(dm.mean())}
        for meth in ("F2", comparator, "F5"):
            conf, call = scores[L][meth]["conf"][m], scores[L][meth]["call"][m]
            sel = e4.select_at_coverage(conf, ids, 0.8, amend)
            cor, _ = correct_of(call, key)
            row[meth] = {"acc_called_disagree@0.8": float(cor[sel & dm].mean()) if np.any(sel & dm) else None,
                         "acc_called_agree@0.8": float(cor[sel & ~dm].mean()) if np.any(sel & ~dm) else None,
                         "share_disagree_called@0.8": float(np.mean(sel[dm])) if dm.any() else None}
        dis[L] = row
    out["E4"]["channel_disagreement"] = dis
    return out


def verdict_text(prim: dict, comparator: str, untestable: bool, bridge_failure: bool) -> dict:
    p = prim[f"F5_vs_comparator@{COVERAGES[0]}"]
    f = prim[f"F5_vs_F2@{COVERAGES[0]}"]
    prim_map = {"win": "field adds value", "loss": "field hurts", "equivalent": "F5 equivalent to the best non-ANM fusion",
                "inconclusive": "inconclusive"}
    fe_map = {"equivalent": "field adds nothing (acceptable outcome)",
              "win": "contradiction events and field gain add over the trust-weighted rule",
              "loss": "contradiction events and field gain hurt relative to the trust-weighted rule", "inconclusive": "inconclusive"}
    s = (f"primary F5 vs {comparator} at 0.8: {p['decision']} ({prim_map[p['decision']]}), point {p['point']:+.4f} "
         f"[{p['ci95'][0]:+.4f}, {p['ci95'][1]:+.4f}]; field effect F5 vs F2: {f['decision']} ({fe_map[f['decision']]}), "
         f"point {f['point']:+.4f} [{f['ci95'][0]:+.4f}, {f['ci95'][1]:+.4f}]")
    if untestable:
        s += "; trust mechanism untestable (val trust gap < 0.01 at every level)"
    if bridge_failure:
        s += "; ANM BRIDGE FAILURE (engine != closed form): interpretation stops"
    return {"primary": {"decision": p["decision"], "outcome": prim_map[p["decision"]]},
            "field_effect": {"decision": f["decision"], "outcome": fe_map[f["decision"]]},
            "trust_untestable": untestable, "anm_bridge_failure": bridge_failure, "summary": s}


def _f(x, n=4):
    return "n/a" if x is None else f"{x:.{n}f}"


def write_report(path: Path, R: dict, add: dict) -> None:
    comp = add["computed"]
    L = ["# E4: two-channel fusion (C6, Experiment 2 redesigned)", "",
         f"Registration sha256 `{R['registration_sha256']}`; amendment A1 `{R['amendment_A1_sha256']}`; addendum E4 "
         f"`{R['addendum_sha256']}` (core `{R['addendum_core_sha256'][:16]}`).{' SMOKE RUN: ' + R['smoke'] if R['smoke'] else ''}", "",
         "## Verdict", "", R["verdict"]["summary"], "",
         f"Comparator (A1.9, chosen on val before site4): **{R['comparator']}** ({R['arm_names'][R['comparator']]}). "
         f"Margin {R['margin']}. A win needs the point difference >= margin, a bootstrap lower bound > 0 and >= margin/2 in "
         "each primary donor. With 2 primary donors the interval describes these donors, not a population.", "",
         "## Key validity (reported before any method result)", "",
         "| split | E4 key B / T / NK / myeloid / OUT | E4 key kappa vs annotation | agreement with the v3 primary key (scored cells) |",
         "|---|---|---:|---:|"]
    for sname, kk in R["key"].items():
        cnt = kk["E4_counts"]
        L.append(f"| {sname} | {' / '.join(str(cnt[c]) for c in vk.CLASSES)} | {kk['E4_kappa_vs_annotation']:.4f} | "
                 f"{kk['E4_vs_v3_primary_agreement_on_scored']:.4f} |")
    L += ["", f"Val kappa at tau_K: {max(r['val_kappa'] for r in comp['key']['tau_grid_val_kappa']):.4f}; training cells "
          f"{comp['key']['train_kappa_vs_annotation']:.4f} (per batch in the addendum). Key-validity flag (test_primary "
          f"kappa < 0.85): {R['verdict'].get('key_validity_flag')}.", "",
          "## Val trust (fixed in the addendum)", "", "| level | channel 1 (TEDDY + head) | channel 2 (protein regressor) | "
          "source scale 1 | source scale 2 |", "|---|---:|---:|---:|---:|"]
    for lv in e4.LEVELS:
        t, s = R["trust_val"][lv], R["source_scale"][lv]
        L.append(f"| {lv} | {t[e4.CHANNELS[0]]:.4f} | {t[e4.CHANNELS[1]]:.4f} | {s[e4.CHANNELS[0]]:.4f} | {s[e4.CHANNELS[1]]:.4f} |")
    L += ["", f"Trust mechanism untestable: {R['trust_untestable']}. E4 key tau_K {comp['key']['tau_K']}; ridge alpha "
          f"{comp['channel2']['alpha']}; F4 C {comp['F4']['C']}; F3 C per level "
          f"{ {lv: comp['F3'][lv]['C'] for lv in e4.LEVELS} }.", ""]
    for sname, S in R["splits"].items():
        for kname in ("E4", "v3_primary"):
            res = S[kname]
            tag = "primary" if (sname == "test_primary" and kname == "E4") else "secondary"
            L += [f"## {sname}, {kname} key ({tag})", "", f"Cells {S['n_cells']}; donors {S['donors']}.", "",
                  "| method | L0 @0.8 | L1 @0.8 | L2 @0.8 | L3 @0.8 | mean @0.8 [95% CI] | mean @0.9 |", "|---|---:|---:|---:|---:|---:|---:|"]
            for meth in METHODS:
                pm = res["per_method"][meth]
                lv = [pm["per_level"][x]["selective_accuracy"]["0.8"] for x in e4.LEVELS]
                m8, m9 = pm["mean_L0_L3@0.8"], pm["mean_L0_L3@0.9"]
                L.append(f"| {meth} | " + " | ".join(_f(v) for v in lv) +
                         f" | {_f(m8['point'])} [{_f(m8['ci95'][0])}, {_f(m8['ci95'][1])}] | {_f(m9['point'])} |")
            L += ["", "| comparison | coverage | difference [95% CI] | per donor | decision |", "|---|---:|---:|---|---|"]
            for cname, c in res["comparisons"].items():
                nm, cov = cname.split("@")
                L.append(f"| {c['a']} - {c['b']} ({nm}) | {cov} | {c['point']:+.4f} [{c['ci95'][0]:+.4f}, {c['ci95'][1]:+.4f}] | "
                         + ", ".join(f"{k} {v:+.4f}" for k, v in c["per_donor"].items()) + f" | {c['decision']} |")
            L.append("")
        if "operating_points" in S["E4"] and S["E4"]["operating_points"]:
            L += [f"### {sname}: operating points at the val bars (E4 key, descriptive)", "",
                  "| method | " + " | ".join(f"{x} coverage / accuracy" for x in e4.LEVELS) + " |", "|---|" + "---:|" * 4]
            for meth in METHODS:
                op = S["E4"]["operating_points"][meth]
                L.append(f"| {meth} | " + " | ".join(f"{op[x]['coverage']:.3f} / {_f(op[x]['selective_accuracy'], 3)}" for x in e4.LEVELS) + " |")
            L += ["", f"### {sname}: channel disagreement (E4 key, descriptive)", "",
                  "| level | share of cells whose channel top classes differ |", "|---|---:|"]
            for x in e4.LEVELS:
                L.append(f"| {x} | {S['E4']['channel_disagreement'][x]['share_channels_disagree']:.4f} |")
            L.append("")
    ck = R["checks"]
    L += ["## Checks", "", f"- ANM engine vs re-coded closed form, every cell and variant: bridge failure = "
          f"{ck['anm_bridge_failure']} (max |difference| per level and variant in E4_results.json).",
          f"- C3 (F5 without contradiction events) vs F2 at matched coverage: cells whose selection differs, per level and "
          f"coverage: { {k: v['n_selection_differs'] for k, v in ck['C3_equals_F2'].items()} }.",
          f"- L1 re-embedding: {R['L1_embedding']}.", "",
          "## Arms", ""] + [f"- **{k}**: {v}" for k, v in R["arm_names"].items()] + [
          "", "## Notes", "",
          "- F2 is a declared rule, not ANM. F5 is the only arm where ANM's field runs; C3 shows that without "
          "contradiction events F5 ranks cells exactly as F2, so F5 - F2 is the effect of the contradiction events "
          "together with the field's event-count-dependent gain.",
          "- In this setup every event enters at t = 0 (time-blind), so ANM's retention contributes only a constant "
          "gain per star size; C2 (channel 2 at t = 2) shows what a declared order changes.",
          "- The E4 key is fixed by measured panel proteins that neither channel reads (channel 2 reads the 122 other "
          "proteins through the A1.2 panel-free transform). The v3 primary key is secondary and partly circular for "
          "channel 2.",
          "- Coverage 0.9 and test_secondary are reported with the same rules and never change the verdict."]
    path.write_text("\n".join(L) + "\n")


# ============================================================================ main
def main(argv=None) -> int:
    global _LOG
    a = parse_args(argv)
    _LOG = a.out_dir / "progress.log"
    a.out_dir.mkdir(parents=True, exist_ok=True)
    say(f"{SCRIPT_VERSION} stage {a.stage}{' pool ' + a.pool if a.pool else ''}{' SMOKE: ' + a.smoke if a.smoke else ''}")
    reg = va.load_registration_amended(a.registration_dir / "registration_v3.json", a.registration_dir / "amendment_A1.json")
    amend = reg["amendment_A1"]
    if a.stage == "embed":
        stage_embed(a, reg)
    elif a.stage == "prepare":
        stage_prepare(a, reg, amend)
    elif a.stage == "leakage":
        stage_leakage(a, reg, amend)
    elif a.stage == "commit-addendum":
        stage_commit_addendum(a)
    elif a.stage == "evaluate":
        stage_evaluate(a, reg, amend)
    return 0


if __name__ == "__main__":
    sys.exit(main())
