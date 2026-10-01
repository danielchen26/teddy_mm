#!/usr/bin/env python3
"""E5 (C7), Mode A layer dynamics: follow gene-token perturbations through TEDDY-G's 12 layers with
Jacobian-vector products and find the layer at which their NK-T component is attenuated.

Registered design (registration/registration_v3.json, experiments.E5; this script's addendum
registration/addenda/E5.json fixes only what the registration leaves open, adds pre-specified secondary
analyses and states two infeasibilities of the registered procedure):

    h_0 = token embedding + position embedding  (official TEDDY-G input, no CLS)
    h_l = f_l(h_{l-1}),  l = 1..12               (the checkpoint's own encoder layers, post-norm)
    dh_l = J_l dh_{l-1},  J_l = df_l/dh at h_{l-1}   (torch.func.jvp through each layer in turn)

For a post-norm layer f_l(h) = LN2(u + FF(u)), u = LN1(h + SA(h)), so J_l = dLN2 (I + dFF) dLN1 (I + dSA):
the "(I + df/dh)" residual step sits inside the two LayerNorms. PyTorch's fused attention kernels have no
forward-mode derivative, so the JVP runs through explicit attention arithmetic with each layer's own weights
(checked against the module forward for every cell); finite differences run through the module itself.

  cells       E5 cells = unique test_primary cells (donors 13272, 19593) of E3's NK-T pairs under the v3
              primary key, at most 200 per donor (registered); plus a matched random set (secondary)
  d_l         unit(mean training NK - mean training T) of the depth-l gene-mean, v3 primary key
  pushes      self: every E2 gene present, along its own embedding (registered); rand: a random unit vector
              at the same position (registered control); other: 2 random non-E2 tokens (registered control);
              nkt: the 7 declared NK/T genes along d_0 (secondary)
  check       central FD at eps 1e-3 / 1e-2 on the unit tangent, float32 (registered gate), and float64 on
              the CPU for a declared subset (secondary S5)
  endpoint    l* = argmin_l median log g_l, log g_l = log|r_l| - log|r_(l-1)|, r_l = <dz_l, d_l>, over the
              self pushes of the E5 cells; compared with the layer of the largest drop of the probe gap ratio
              (outputs/mode_a_official); specificity vs rand at l* (margin 0.1, two-stage donor bootstrap)
  secondary   S1 nkt-push gain G_l and step change dG_l (localisation), S2 look-alike vs random cells,
              S3 head transmission (head-score lens), S4 linearity on an eps ladder (one-sided slopes),
              S5 float64 FD check, S6 float64 JVP

Stages
  register    write registration/addenda/E5.json + its line in addenda/HASHES.txt
  select      rebuild the E3 pairs on test_primary, choose the cells (no forward pass)
  directions  CPU: d_l per depth, probe directions, head-score gaps (training + val cells)
  respond     GPU/MPS: per-cell JVP + finite differences, one file per cell (resumable)
  report      CPU: two-stage bootstrap, registered endpoint, secondary analyses; E5_results.json + REPORT.md
  all         select, directions, respond, report
A site4 run (any stage but select) is refused unless registration_v3.json and the addendum are committed and
match. Smoke runs use val-donor cells (--cell-pool val) only.

Full run (from the repo root; not run inside the workflow). Paths to the main checkout's data and outputs:
  PY=<venv>/bin/python; R=/Users/tianchichen/Documents/GitHub/teddy_mm
  PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_layer_dynamics.py --stage all \
      --processed $R/data/processed/cite --embed-dir $R/data/processed/cite_official \
      --ckpt /Users/tianchichen/Documents/GitHub/teddy_mwe/ckpt/teddy_g_70M \
      --medians $R/data/reference/teddy_gene_medians.json \
      --head-ckpt $R/outputs/cite_phase1_official/best.pt --probe-dir $R/outputs/mode_a_official \
      --out-dir $R/outputs/v3/E5 --device auto --threads 4
Every stage appends to <out-dir>/progress.log; the respond stage logs seconds per cell and an ETA after every
cell (also in e5_progress.json); rerunning the same command skips finished cells (each cell file carries the
addendum hash and is refused under another).
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
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("mode_a_layer_probes", HERE / "mode_a_layer_probes.py")
lp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lp)

SCRIPT_VERSION = "mode_a_layer_dynamics v3"
PROTEINS, GDT158 = lp.PROTEINS, lp.GDT158  # CD56, CD94, CD335, CD3
rnd = lp.rnd
N_DEPTHS = 13  # 0 = input state, 1..12 = layer outputs
# directions: "self" (registered perturbation: the gene's own embedding), "rand" (registered control: a random unit
# vector at the same position), "other" (registered control: a random other gene of the cell, own embedding),
# "nkt" (secondary: the training NK-T direction of the input gene-mean, declared NK/T genes only)
DIRECTIONS = ("self", "rand", "other", "nkt")
# secondary nkt pushes: symbol -> class (Ensembl ids come from the gene reference like every E2 gene)
NKT_GENES = {"CD3E": "T", "CD3D": "T", "NCAM1": "NK", "KLRD1": "NK", "NCR1": "NK", "FCGR3A": "NK", "KLRF1": "NK"}
# the 13 NK/T genes of E2 must map to the ids of scripts/mode_a_token_probes.py (checked at run time)
E2_NKT_CLASS = {"NCAM1": "NK", "KLRD1": "NK", "NCR1": "NK", "FCGR3A": "NK", "KLRF1": "NK", "KLRC1": "NK", "SH2D1B": "NK",
                "CD3E": "T", "CD3D": "T", "CD3G": "T", "CD5": "T", "CD6": "T", "CD28": "T"}  # as mode_a_token_probes.py
TOKEN_PROBE_IDS = {"NCAM1": "ENSG00000149294", "KLRD1": "ENSG00000134539", "NCR1": "ENSG00000189430",
                   "CD3E": "ENSG00000198851", "CD3D": "ENSG00000167286", "CD3G": "ENSG00000160654",
                   "CD5": "ENSG00000110448", "CD6": "ENSG00000013725", "CD28": "ENSG00000178562",
                   "KLRF1": "ENSG00000150045", "FCGR3A": "ENSG00000203747", "KLRC1": "ENSG00000134545",
                   "SH2D1B": "ENSG00000198574"}
GENE_REFERENCE = "data/raw/external_hao2021_gse164378/ref_10x_grch38_3.0.0/features_grch38_3.0.0.tsv.gz"
E5_CAP_PER_DONOR = 200
N_RANDOM_PER_DONOR = 100
N_OTHER = 2
N_BOOT = 2000
EPS_ABS_CHECK = (1e-3, 1e-2)
REGISTRATION_EPS = (0.02, 0.1, 1.0)  # secondary linearity ladder
EPS_PRIMARY = 0.02
N_FULL_LADDER = 100  # the first cells of the seeded processing order get every eps rung (secondary)
N_FD64_CELLS = 20  # the first cells of the order get the float64 CPU finite-difference check (secondary)
_LOG = {"path": None}


def say(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG["path"] is not None:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


REGISTRATION = {
    "experiment": "E5",
    "addendum_to": "registration/registration_v3.json, experiments.E5 (primary endpoint, cells, directions, "
                   "perturbation, check, margins, falsification are the registered ones and are not changed here)",
    "addendum_version": 3,
    "script": "scripts/mode_a_layer_dynamics.py (" + SCRIPT_VERSION + ")",
    "fixed_before_site4": "fixed on training and validation cells only and committed before any site4 forward pass of "
                          "this script; smoke runs use val donor 18303 cells only",
    "open_choices_fixed": {
        "cells": "E3 pairs rebuilt by registration e3.pairs_rule as amended by A1.8 (k = 10 cosine neighbours of the raw "
                 "L2-normalised z_rna.npy within each test_primary donor, each edge once, one primary-key NK and one "
                 "primary-key T cell; variant 'all'); unique cells; per primary donor in registration order, if more than "
                 f"{E5_CAP_PER_DONOR}, a draw without replacement with numpy default_rng(seeds.e5_cells), one generator "
                 "for both donors",
        "keys": "primary key from bridge_anm/lib/v3_key.build_keys with the frozen registration (training NK / T for "
                "d_l, test_primary NK / T for the pairs)",
        "genes": "E2 genes = panel_coding_genes (flattened) then nk_t_genes, first occurrence kept; Ensembl ids from "
                 f"{GENE_REFERENCE} (10x GRCh38 3.0.0; checked to give the mode_a_token_probes ids for the 13 NK/T genes); "
                 "genes absent from the processed genes or the TEDDY vocabulary are dropped and listed",
        "tangent": "delta h_0 = ||E[gene]|| * v at the gene's position (v unit; the log-gain is scale-free); finite "
                   "differences use c = eps / ||E[gene]||, i.e. eps times the unit vector, as registered",
        "rand_control": "one Gaussian unit vector per (cell, gene) case at the same position, drawn in E2 gene order "
                        "with numpy default_rng([seeds.e5_random_directions, cell index])",
        "other_control": f"{N_OTHER} tokens of the cell that are not E2 genes, drawn with the same generator before the "
                         "rand vectors, each pushed along its own embedding",
        "check_scope": "(cell, gene) cases = the self pushes; pass = relative error <= 0.05 and cosine >= 0.99 at every "
                       "depth 1..12 at eps = 1e-3 (float32, the registered precision); eps = 1e-2 reported",
        "endpoint_scope": "l* over the self pushes in the E5 cells; log g_l = log|r_l| - log|r_(l-1)|, r_l = <dz_l, "
                          "unit(d_l)>, l = 1..12",
        "comparator_layer": "the layer whose step has the largest drop of the mean over CD56 / CD94 / CD335 / CD3 of the "
                            "no_gdT158 gap ratio in outputs/mode_a_official/mode_a_results.json (the 'all' variant reported)",
        "specificity": "D = median log g_{l*}(self) - median log g_{l*}(rand) over the E5 cells; two-stage bootstrap "
                       "(common.statistics as amended by A1: B = 2000, seed 1, l* recomputed inside each replicate); "
                       "'no layer localises the NK-T loss' if the 95% interval "
                       "lies inside (-0.1, 0.1) (common 'equivalent'); 'localised at l*' if D <= -0.1, upper bound < 0 "
                       "and D <= -0.05 in each primary donor (common 'win'); otherwise inconclusive",
        "precision_and_path": "float32; every layer evaluated with its own weights and explicit attention arithmetic "
                              "(see infeasibilities); finite differences through the module's own layers",
    },
    "secondary_analyses": {
        "note": "pre-specified, reported separately, never change the registered E5 verdict",
        "random_cells": f"per primary donor {N_RANDOM_PER_DONOR} test_primary NK / T cells in no E3 pair, class "
                        "proportions and token-count quintiles matched to the donor's E5 cells (seeds.e5_cells generator); a "
                        "quintile with too few cells is filled at random from the donor's other cells of the class (logged)",
        "nkt_push": "for the 7 NK/T genes " + ", ".join(NKT_GENES) + ": a push along unit(mean_NK - mean_T) of the "
                    "training input gene-mean; G_l = (d_l . dz_l / ||d_l||^2) / A0, A0 = ||E[gene]|| / (n ||d_0||), "
                    "dG_l = G_l - G_{l-1}; L_E5 = argmin median dG_l",
        "S1_localised_nkt_attenuation": "L_E5 has a two-stage upper bound < 0, is the argmin in >= 50% of the resamples and "
                                        "has median dG < 0 in both primary donors",
        "S2_lookalike_specificity": "median G_12(E5 cells) - median G_12(random cells), nkt: upper bound < 0 and < 0 in "
                                    "both donors -> look-alike cells damp the push more",
        "S3_head_transmission": "median (H_12 - G_12), H = head-score lens (gradient of mean(CD56, CD94, CD335) - CD3, "
                                "predicted / train p95) over the training NK-T head-score gap, per A0",
        "S4_linearity": "one-sided slopes D+-(e) vs the JVP and vs e/2 (tau 0.1) at e = 0.02 on every self and nkt case, "
                        f"and at 0.1 and 1.0 on the first {N_FULL_LADDER} cells of the seeded order; linear range = largest "
                        "e with pass rate >= 0.90 at every depth",
        "S5_float64_check": "the registered central-FD check in float64 on the CPU through the module's own layers, one "
                            f"case per direction of the first {N_FD64_CELLS} cells of the order, against the device float32 JVP",
        "S6_float64_jvp": "the JVP of the first 2 cases of the first 2 cells recomputed in float64 on the CPU",
    },
    "infeasibilities_needing_an_amendment": [
        "'standard Python path': forward-mode AD is not implemented for PyTorch's fused attention kernels "
        "(_transformer_encoder_layer_fwd, _native_multi_head_attention, _scaled_dot_product_flash_attention_for_cpu, "
        "_scaled_dot_product_attention_math_for_mps; torch 2.14, checked), so the JVP runs through explicit attention "
        "arithmetic with each layer's own weights; its forward equals the module's to about 6e-6 (every cell), its JVP "
        "equals a float64 CPU JVP to about 4e-6 relative and float64 central FD of the module to about 6e-5 at 1e-3 "
        "(val smoke)",
        "'CPU': a CPU run of about 600 cells x ~20 cases is days; the device is MPS in float32 (same precision)",
        "the float32 check at eps = 1e-3 reaches the float32 rounding floor for weak responses at depths 9-12: val smoke "
        "4 of 36 cases with relative error 0.05-0.12 and cosine >= 0.993 there (|eps dz| / |z| about 2e-6), every case "
        "passing at 1e-2 and in float64 (S5). If the registered gate fails for this reason, the registered verdict is "
        "'JVP not validated' unless an amendment adopting S5 is committed before site4",
    ],
    "outputs": "outputs/v3/E5/: E5_results.json, REPORT.md, progress.log, cells/ (one file per cell); every result JSON "
               "carries registration_sha256 and addendum_sha256",
}


def registration_bytes() -> bytes:
    return (json.dumps(REGISTRATION, indent=1, sort_keys=True) + "\n").encode()


def registration_sha() -> str:
    return hashlib.sha256(registration_bytes()).hexdigest()


_vspec = importlib.util.spec_from_file_location("v3_key", ROOT / "bridge_anm" / "lib" / "v3_key.py")
vk = importlib.util.module_from_spec(_vspec)
_vspec.loader.exec_module(vk)


# ============================================================================ args

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "select", "directions", "respond", "report", "all"), default="all")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=ROOT / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--probe-dir", type=Path, default=ROOT / "outputs/mode_a_official",
                   help="mode_a_layer_probes.py out-dir (pairs check, probe curves)")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration",
                   help="holds registration_v3.json(.sha256) and addenda/E5.json + addenda/HASHES.txt")
    p.add_argument("--gene-reference", type=Path, default=ROOT / GENE_REFERENCE,
                   help="GRCh38 features file (Ensembl id, symbol) for the E2 genes")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--cell-pool", choices=("site4", "val"), default="site4",
                   help="val = val donor 18303 NK/T cells (smoke tests; no site4 forward)")
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--jvp-chunk", type=int, default=4, help="tangents per vmap call")
    p.add_argument("--fd-batch", type=int, default=4, help="perturbed sequences per forward call")
    p.add_argument("--cpu64-check", type=int, default=2,
                   help="for the first N processed cells, recompute the JVP of their first 2 cases in float64 on the CPU")
    p.add_argument("--limit", type=int, default=0, help="smoke only: process just the first N cells of the order")
    p.add_argument("--max-genes", type=int, default=0, help="smoke only: at most N genes per cell")
    p.add_argument("--dir-max-train", type=int, default=0,
                   help="smoke only: at most N training (and N val) cells per class for the directions")
    p.add_argument("--n-boot", type=int, default=N_BOOT, help="registered: 2000; smaller only for smoke runs")
    p.add_argument("--n-full-ladder", type=int, default=None, help="smoke only: override the registered n_full_ladder")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    if a.stage != "register" and a.out_dir is None:
        p.error("--out-dir is required")
    if (a.limit or a.max_genes or a.dir_max_train or a.n_boot != N_BOOT or a.n_full_ladder is not None) and not a.smoke:
        p.error("--limit / --max-genes / --dir-max-train / --n-boot / --n-full-ladder are for smoke runs only (give --smoke NOTE)")
    if a.cell_pool == "val" and not a.smoke:
        p.error("--cell-pool val is for smoke runs only")
    return a


# ============================================================================ helpers

def load_reg(a) -> dict:
    """The frozen registration (hash verified against registration_v3.json.sha256)."""
    return vk.load_registration(a.registration_dir / "registration_v3.json", verify_hash=True)


def primary_keys(meta: dict, reg: dict) -> np.ndarray:
    keys = vk.build_keys(meta["cell_types"], meta["adt"], meta["adt_names"], reg)
    return np.asarray(keys["primary"]).astype(str)


def e2_genes(reg: dict, gene_ref: Path, rna_names, vocab: dict) -> tuple[list, list, dict]:
    """E2 genes in registration order (panel coding genes, then NK/T genes; first occurrence kept), mapped to Ensembl
    ids with the gene reference; returns (kept [(symbol, ensembl, vocab id)], dropped [symbol], info)."""
    import gzip
    G = reg["experiments"]["E2"]["perturbations"]["genes"]
    syms = []
    for v in G["panel_coding_genes"].values():
        syms += v if isinstance(v, list) else [v]
    syms += list(G["nk_t_genes"])
    syms = list(dict.fromkeys(syms))
    sym2ens = {}
    with gzip.open(gene_ref, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in syms and parts[1] not in sym2ens:
                sym2ens[parts[1]] = parts[0]
    bad = {s_: (sym2ens.get(s_), e) for s_, e in TOKEN_PROBE_IDS.items() if sym2ens.get(s_) != e}
    if bad:
        raise SystemExit(f"gene reference disagrees with mode_a_token_probes ids: {bad}")
    rn = set(map(str, rna_names))
    kept, dropped = [], []
    for s_ in syms:
        e = sym2ens.get(s_)
        if e is not None and e in rn and e in vocab:
            kept.append((s_, e, int(vocab[e])))
        else:
            dropped.append(s_)
    missing_nkt = [g for g in NKT_GENES if g not in [k[0] for k in kept]]
    if missing_nkt:
        raise SystemExit(f"declared NK/T genes not available: {missing_nkt}")
    info = {"reference": str(gene_ref), "reference_sha256": vk.sha256_file(gene_ref),
            "kept": {k[0]: k[1] for k in kept}, "dropped": dropped}
    return kept, dropped, info


def load_donors(processed: Path) -> np.ndarray:
    with np.load(processed / "cite_arrays.npz", allow_pickle=True) as d:
        return d["donors"].astype(str)


def atomic_savez(path: Path, **arrays):
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def write_json(path: Path, obj):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1))
    os.replace(tmp, path)


# ============================================================================ register

def stage_register(a):
    d = a.registration_dir / "addenda"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "E5.json"
    f.write_bytes(registration_bytes())
    sha = registration_sha()
    hf = d / "HASHES.txt"
    lines = [ln for ln in (hf.read_text().splitlines() if hf.exists() else []) if not ln.rstrip().endswith("  E5.json")]
    hf.write_text("\n".join(lines + [f"{sha}  E5.json"]) + "\n")
    say(f"wrote {f} and its line in {hf} (sha256 {sha})")


def _git_committed(path: Path) -> bool | None:
    import subprocess
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def check_registration(a) -> dict:
    sha = registration_sha()
    f = a.registration_dir / "addenda" / "E5.json"
    hf = a.registration_dir / "addenda" / "HASHES.txt"
    info = {"addendum_sha256": sha, "addendum_file": str(f),
            "addendum_matches_script": f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() == sha,
            "addendum_hash_recorded": hf.exists() and f"{sha}  E5.json" in hf.read_text().splitlines(),
            "addendum_committed": _git_committed(f) if f.exists() else False}
    cf = a.registration_dir / "registration_v3.json"
    ch = a.registration_dir / "registration_v3.json.sha256"
    info["registration_sha256"] = None
    if cf.exists():
        h = hashlib.sha256(cf.read_bytes()).hexdigest()
        info["registration_sha256"] = h
        info["registration_hash_file_matches"] = ch.exists() and ch.read_text().split()[0] == h
        info["registration_committed"] = _git_committed(cf)
    info["amendments"] = {}
    for am in sorted(a.registration_dir.glob("amendment_*.json")):
        if am.name.endswith("_core.json"):
            continue
        h = hashlib.sha256(am.read_bytes()).hexdigest()
        hf_ = am.with_name(am.name + ".sha256")
        info["amendments"][am.name] = {"sha256": h, "committed": _git_committed(am),
                                       "hash_file_matches": hf_.exists() and hf_.read_text().split()[0] == h}
    if a.cell_pool == "site4" and a.stage != "select":  # select draws cells only (no forward pass, no E5 outcome)
        need = {"addendum file equals the script's REGISTRATION": info["addendum_matches_script"],
                "addendum sha256 recorded in addenda/HASHES.txt": info["addendum_hash_recorded"],
                "addendum committed": info["addendum_committed"] is True,
                "registration_v3.json present, hash file matching, committed":
                    bool(info["registration_sha256"]) and info.get("registration_hash_file_matches") is True
                    and info.get("registration_committed") is True,
                "every amendment committed with a matching hash file":
                    all(v["committed"] is True and v["hash_file_matches"] for v in info["amendments"].values())}
        bad = [k for k, ok in need.items() if not ok]
        if bad:
            raise SystemExit("site4 run refused (register and commit before any site4 forward): " + "; ".join(bad))
    return info


# ============================================================================ cells

def stage_select(a, meta, emb, key, reg) -> dict:
    """E5 cells (registered) + the secondary random set; key = registration primary key."""
    out = a.out_dir
    split = meta["split"]
    donors = load_donors(a.processed)
    if emb.partial:
        raise SystemExit("needs the full official embedding")
    ntok = emb.ntokens
    rng = np.random.default_rng(int(reg["seeds"]["e5_cells"]))
    if a.cell_pool == "val":
        sel, role = [], []
        for cls in ("NK", "T"):
            pool = np.where((split == "val") & (key == cls))[0]
            m = min(pool.size, 300)
            sel.append(np.sort(rng.choice(pool, m, replace=False)))
            role += ["val"] * m
        sel = np.concatenate(sel)
        role = np.array(role)
        info = {"pool": "val donor 18303 primary-key NK/T cells (smoke)"}
    else:
        spec = reg["experiments"]["E5"]["cells"]
        if f"at most {E5_CAP_PER_DONOR} per primary donor" not in spec or "k = 10" not in reg["e3"]["pairs_rule"]:
            raise SystemExit("registered E5 cell rule / E3 pairs rule differ from what this script implements")
        prim = list(reg["splits"]["test_primary"]["donors"])
        ev = np.where((split == "test") & np.isin(donors, prim))[0]
        zall = np.asarray(emb.z, dtype=np.float32)
        t0 = time.time()
        P = []  # amendment A1.8: cosine k-NN within each primary donor (knn_pairs L2-normalises), each edge once
        for dn in prim:
            evd = ev[donors[ev] == dn]
            P.append(evd[lp.knn_pairs(zall[evd], 10)])
        P = np.concatenate(P)
        ka, kb = key[P[:, 0]], key[P[:, 1]]
        nkt = ((ka == "NK") & (kb == "T")) | ((ka == "T") & (kb == "NK"))
        Q = P[nkt]
        say(f"E3 pairs rebuilt within each test_primary donor ({time.time() - t0:.0f}s): {len(P)} edges, {len(Q)} NK-T pairs")
        cells_in = np.unique(Q.ravel())
        e5 = []
        for dn in prim:
            c = cells_in[donors[cells_in] == dn]
            if c.size > E5_CAP_PER_DONOR:
                c = np.sort(rng.choice(c, E5_CAP_PER_DONOR, replace=False))
            e5.append(c)
        e5 = np.concatenate(e5)
        in_pair = np.zeros(len(split), bool)
        in_pair[cells_in] = True
        rand, shortfalls = [], {}
        for dn in prim:
            ref_d = e5[donors[e5] == dn]
            for cls in ("NK", "T"):
                ref = ref_d[key[ref_d] == cls]
                want = int(round(N_RANDOM_PER_DONOR * ref.size / max(ref_d.size, 1)))
                if ref.size == 0 or want == 0:
                    continue
                pool = ev[(donors[ev] == dn) & (key[ev] == cls) & ~in_pair[ev]]
                edges = np.quantile(ntok[ref], [0.2, 0.4, 0.6, 0.8])
                sp, sr = np.digitize(ntok[pool], edges), np.digitize(ntok[ref], edges)
                shares = np.bincount(sr, minlength=5) / sr.size
                cnt = np.floor(shares * want).astype(int)
                for q in np.argsort(-(shares * want - cnt))[: want - cnt.sum()]:
                    cnt[q] += 1
                got = []
                for q in range(5):
                    cand = pool[sp == q]
                    m = min(int(cnt[q]), cand.size)
                    if m:
                        got.append(rng.choice(cand, m, replace=False))
                got = np.concatenate(got) if got else np.zeros(0, np.int64)
                short = want - got.size
                if short > 0:  # a quintile ran short: fill at random from the donor's other cells of the class
                    rest = np.setdiff1d(pool, got)
                    if rest.size < short:
                        raise SystemExit(f"random pool too small: donor {dn} class {cls}")
                    got = np.concatenate([got, rng.choice(rest, short, replace=False)])
                    shortfalls[f"{dn}/{cls}"] = int(short)
                rand.append(got)
        rand = np.sort(np.concatenate(rand)) if rand else np.zeros(0, np.int64)
        sel = np.concatenate([e5, rand])
        role = np.array(["pair"] * e5.size + ["random"] * rand.size)
        info = {"pool": "site4 test_primary donors " + ", ".join(prim), "n_e3_edges": int(len(P)), "n_e3_nkt_pairs": int(len(Q)),
                "n_e3_pair_cells_per_donor": {dn: int(np.sum(donors[cells_in] == dn)) for dn in prim},
                "random_quintile_shortfalls_filled_at_random": shortfalls,
                "e3_pairs_sha256": vk.index_hash(sorted(int(x) * 1_000_000 + int(y) for x, y in np.sort(Q, axis=1)))}
    if np.unique(sel).size != sel.size:
        raise SystemExit("duplicate cells")
    order = rng.permutation(sel.size)  # processing order: roles interleaved, so a partial run is a random subset
    info.update({"n_cells": int(sel.size),
                 "by_role_donor_class": {f"{r}/{dn}/{c}": int(np.sum((role == r) & (donors[sel] == dn) & (key[sel] == c)))
                                         for r in np.unique(role) for dn in np.unique(donors[sel]) for c in ("NK", "T")},
                 "ntokens_median_by_role": {str(r): float(np.median(ntok[sel[role == r]])) for r in np.unique(role)}})
    atomic_savez(out / "e5_cells.npz", cells=sel, role=role, order=order, cls=key[sel], donor=donors[sel])
    write_json(out / "e5_cells.json", info)
    say(f"cells: {info['by_role_donor_class']}; token-count medians {info['ntokens_median_by_role']}")
    return info


def load_cells(out: Path):
    with np.load(out / "e5_cells.npz", allow_pickle=False) as d:
        return d["cells"], d["role"].astype(str), d["order"], d["cls"].astype(str)


# ============================================================================ head

def load_head(path: Path, meta: dict, device):
    import torch

    from teddy_mm.models import MLP, AdtDecoder
    blob = torch.load(path, map_location="cpu", weights_only=False)
    w0 = blob["mlp"]["net.0.weight"]
    hidden, z_dim = int(w0.shape[0]), int(w0.shape[1])
    n_adt = meta["adt"].shape[1]
    mlp, dec = MLP(z_dim, z_dim, hidden=hidden), AdtDecoder(z_dim, n_adt, hidden=hidden)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval().to(device)
    dec.eval().to(device)
    sf = float(np.median(meta["adt_size_factor"][meta["split"] == "train"]))
    p95 = lp.p95_replicate(meta)
    col = {n: j for j, n in enumerate(meta["adt_names"])}
    cols = torch.tensor([col[p] for p in PROTEINS], device=device)
    p95v = torch.tensor([p95[p] for p in PROTEINS], dtype=torch.float32, device=device)

    def yhat(z):  # z: [B, d] raw gene-mean -> [B, 4] predicted ADT / train p95 (unclipped)
        zn = z / (torch.linalg.vector_norm(z, dim=-1, keepdim=True) + 1e-6)
        mu = dec(mlp(zn), torch.full((z.shape[0],), sf, device=z.device))[0]
        return mu.index_select(-1, cols) / p95v

    def score(y):  # [B, 4] -> [B]
        return y[:, :3].mean(-1) - y[:, 3]

    info = {"path": str(path), "hidden": hidden, "z_dim": z_dim, "train_median_size_factor": sf,
            "p95": {p: float(p95[p]) for p in PROTEINS}}
    return yhat, score, info


# ============================================================================ directions

def stage_directions(a, meta, emb, key) -> dict:
    import torch
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss
    torch.set_num_threads(_THREADS)
    out = a.out_dir
    split = meta["split"]
    rng = np.random.default_rng(0)
    rows = {}
    for sp in ("train", "val"):
        r = []
        for cls in ("NK", "T"):
            pool = np.where((split == sp) & (key == cls))[0]
            if a.dir_max_train and pool.size > a.dir_max_train:
                pool = np.sort(rng.choice(pool, a.dir_max_train, replace=False))
            r.append(pool)
        rows[sp] = np.sort(np.concatenate(r))
    both = np.concatenate([rows["train"], rows["val"]])
    o = np.argsort(both)
    both_sorted = both[o]
    t0 = time.time()
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    X0, info0 = lp.input_layer_means(meta, both_sorted, a.ckpt, a.medians, seq_len, norm, emb.ntokens[both_sorted])
    say(f"directions: input gene-means of {both.size} train/val NK-T cells in {time.time() - t0:.0f}s")
    pos = {c: i for i, c in enumerate(both_sorted)}
    D = emb.d
    dmu = np.zeros((N_DEPTHS, D), np.float64)
    wraw = np.zeros((N_DEPTHS, D), np.float64)
    probe = []
    ytr = (key[rows["train"]] == "NK").astype(int)
    yva = (key[rows["val"]] == "NK").astype(int)
    for l in range(N_DEPTHS):
        if l == 0:
            Xtr = X0[[pos[c] for c in rows["train"]]]
            Xva = X0[[pos[c] for c in rows["val"]]]
        else:
            Xtr = np.asarray(emb.layer_means[l - 1][rows["train"]], np.float32)
            Xva = np.asarray(emb.layer_means[l - 1][rows["val"]], np.float32)
        dmu[l] = Xtr[ytr == 1].mean(0, dtype=np.float64) - Xtr[ytr == 0].mean(0, dtype=np.float64)
        mu, sd = Xtr.mean(0, dtype=np.float64), Xtr.std(0, dtype=np.float64)
        sd[sd < 1e-8] = 1.0
        best = None
        for C in lp.LR_C_GRID:
            clf = LogisticRegression(C=C, max_iter=3000, class_weight="balanced", random_state=0)
            clf.fit((Xtr - mu) / sd, ytr)
            ll = log_loss(yva, clf.predict_proba((Xva - mu) / sd)[:, 1], labels=[0, 1])
            if best is None or ll < best[0]:
                best = (ll, C, clf)
        w = best[2].coef_[0] / sd
        wraw[l] = w
        probe.append({"depth": l, "C": best[1], "val_logloss": rnd(best[0]), "logit_gap": float(w @ dmu[l]),
                      "dmu_norm": float(np.linalg.norm(dmu[l])),
                      "cos_probe_vs_mudiff": float(w @ dmu[l] / (np.linalg.norm(w) * np.linalg.norm(dmu[l])))})
        say(f"  depth {l:2d}: ||dmu|| {probe[-1]['dmu_norm']:.4f}  probe C {best[1]}  cos(probe, dmu) "
            f"{probe[-1]['cos_probe_vs_mudiff']:.3f}")
    # head-score gap on training NK / T cells (official z, the head's real input)
    yhat, score, hinfo = load_head(a.head_ckpt, meta, torch.device("cpu"))
    z = np.asarray(emb.z[rows["train"]], np.float32)
    with torch.no_grad():
        Y = torch.cat([yhat(torch.from_numpy(z[s:s + 4096])) for s in range(0, z.shape[0], 4096)]).numpy().astype(np.float64)
    S = Y[:, :3].mean(1) - Y[:, 3]
    head_gap_s = float(S[ytr == 1].mean() - S[ytr == 0].mean())
    head_gap_p = (Y[ytr == 1].mean(0) - Y[ytr == 0].mean(0))
    u0 = dmu[0] / np.linalg.norm(dmu[0])
    info = {"n_train": {"NK": int(ytr.sum()), "T": int((1 - ytr).sum())}, "n_val": {"NK": int(yva.sum()), "T": int((1 - yva).sum())},
            "subsampled": bool(a.dir_max_train), "probe": probe, "head": hinfo, "head_gap_s": head_gap_s,
            "head_gap_per_protein": dict(zip(PROTEINS, map(float, head_gap_p))), "input_layer": info0,
            "runtime_sec": round(time.time() - t0, 1)}
    atomic_savez(out / "e5_directions.npz", dmu=dmu, wraw=wraw, u_nkt=u0,
                 head_gap_s=np.array(head_gap_s), head_gap_p=head_gap_p)
    write_json(out / "e5_directions.json", info)
    say(f"directions done ({info['runtime_sec']}s); head-score train NK-T gap {head_gap_s:.4f}")
    return info


def load_directions(out: Path) -> dict:
    with np.load(out / "e5_directions.npz") as d:
        return {k: d[k] for k in d.files}


# ============================================================================ respond

def layer_fn(layer, x):
    """One nn.TransformerEncoderLayer (eval: dropout inactive), same weights, explicit attention
    softmax(QK^T / sqrt(d_h)) V so that forward-mode AD runs on every device. x: [B, L, d], no padding."""
    import torch
    import torch.nn.functional as F
    mha = layer.self_attn
    B, L, d = x.shape
    H = mha.num_heads
    dh = d // H

    def sa(y):
        q, k, v = F.linear(y, mha.in_proj_weight, mha.in_proj_bias).split(d, dim=-1)
        q = q.reshape(B, L, H, dh).transpose(1, 2)
        k = k.reshape(B, L, H, dh).transpose(1, 2)
        v = v.reshape(B, L, H, dh).transpose(1, 2)
        w = torch.softmax(torch.matmul(q, k.transpose(-1, -2)) * (dh ** -0.5), dim=-1)
        o = torch.matmul(w, v).transpose(1, 2).reshape(B, L, d)
        return F.linear(o, mha.out_proj.weight, mha.out_proj.bias)

    def ff(y):
        return layer.linear2(layer.activation(layer.linear1(y)))

    if layer.norm_first:
        x = x + sa(layer.norm1(x))
        return x + ff(layer.norm2(x))
    x = layer.norm1(x + sa(x))
    return layer.norm2(x + ff(x))


def tokenise(meta, sel, a, emb, vocab_pad):
    from teddy_mm.teddy_encoder import load_gene_medians, median_factors, official_values, rank_encode_official
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    factors = median_factors(meta["rna_names"], load_gene_medians(a.medians))
    toks, ntok = [], []
    for s in range(0, len(sel), 256):
        vals = official_values(meta["rna"][sel[s:s + 256]].toarray(), factors, norm)
        t, m = rank_encode_official(vals, meta["token_ids"], max_len=seq_len, pad_id=vocab_pad, pad_to=seq_len)
        for i in range(t.shape[0]):
            n = int(m[i].sum())
            toks.append(t[i, :n].copy())
            ntok.append(n)
    return toks, np.array(ntok)


def stage_respond(a, meta, emb, dirs, reg):
    import copy

    import torch
    from torch.func import jvp, vmap

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab
    torch.set_num_threads(_THREADS)
    out = a.out_dir
    cdir = out / "cells"
    cdir.mkdir(exist_ok=True)
    sel, role, order, cls = load_cells(out)
    n_full = a.n_full_ladder if a.n_full_ladder is not None else N_FULL_LADDER
    full_set = {int(i) for i in order[:n_full]}  # cells with every eps rung (the first of the seeded order)
    fd64_set = {int(i) for i in order[:N_FD64_CELLS]}
    if a.limit:
        order = order[: a.limit]
    device = resolve_device(a.device)
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    genes, dropped, ginfo = e2_genes(reg, a.gene_reference, meta["rna_names"], vocab)
    write_json(out / "e5_genes.json", ginfo)
    tok2gene = {g[2]: j for j, g in enumerate(genes)}
    nkt_idx = {j for j, g in enumerate(genes) if g[0] in NKT_GENES}
    seed_rand = int(reg["seeds"]["e5_random_directions"])
    model = load_teddy(a.ckpt, device)
    for p_ in model.parameters():
        p_.requires_grad_(False)
    layers = list(model.encoder.layers)
    d_model = model.d_model
    E = model.embeddings.weight
    Pm = model.position_embeddings.weight
    yhat, score, _ = load_head(a.head_ckpt, meta, device)
    f32 = lambda x: torch.as_tensor(np.asarray(x, np.float32), device=device)  # noqa: E731
    dmu, wraw, u_nkt = f32(dirs["dmu"]), f32(dirs["wraw"]), f32(dirs["u_nkt"])
    eps_list = [float(e) for e in REGISTRATION_EPS]
    si0 = eps_list.index(EPS_PRIMARY)
    sha = registration_sha()
    todo = []
    for i in order:
        f = cdir / f"cell_{int(sel[i])}.npz"
        if f.exists():
            with np.load(f, allow_pickle=False) as d:
                if json.loads(str(d["meta"])).get("registration_sha256") != sha:
                    raise SystemExit(f"{f} was written under another addendum; use a fresh --out-dir")
        else:
            todo.append(int(i))
    say(f"respond: {len(order)} cells in the order, {len(order) - len(todo)} already done, {len(todo)} to go; "
        f"device {device}, float32, {len(genes)} E2 genes ({len(dropped)} dropped: {dropped}), jvp chunk {a.jvp_chunk}, "
        f"fd batch {a.fd_batch}")
    if not todo:
        return
    toks, ntok = tokenise(meta, sel[todo], a, emb, pad_id)
    ref_ntok = emb.ntokens[sel[todo]]
    if np.any(ntok != ref_ntok):
        raise SystemExit(f"{int(np.sum(ntok != ref_ntok))} cells get another token count than the official run")
    t_start = time.time()
    n_done = n_cpu64 = 0
    layers64 = None

    def pooled_and_own(h, pos_k):  # h: [K, 1, L, d] -> pooled [K, d], own [K, d]
        return h[:, 0].mean(1), h[torch.arange(h.shape[0], device=h.device), 0, pos_k]

    def forward_all(h0):  # explicit arithmetic (the JVP's primal): h0 [B, L, d] -> states per depth
        hs = [h0]
        h = h0
        for layer in layers:
            h = layer_fn(layer, h)
            hs.append(h)
        return hs

    def forward_module(h0):  # the module's own layer call (finite differences): pooled per depth [B, 13, d]
        zs = [h0.mean(1)]
        h = h0
        for layer in layers:
            h = layer(h)  # every position is a real token, so there is no padding mask
            zs.append(h.mean(1))
        return torch.stack(zs, 1)

    for jj, i in enumerate(todo):
        tc = time.time()
        g_id = int(sel[i])
        ids_np = toks[jj]
        L = ids_np.size
        ids = torch.from_numpy(ids_np).to(device)
        posmap = {int(t): q for q, t in enumerate(ids_np)}
        present = [(j, posmap[g[2]]) for j, g in enumerate(genes) if g[2] in posmap]
        if a.max_genes:
            present = present[: a.max_genes]
        rng_c = np.random.default_rng([seed_rand, g_id])
        cand = np.array([q for q in range(L) if int(ids_np[q]) not in tok2gene])
        other_pos = np.sort(rng_c.choice(cand, min(N_OTHER, cand.size), replace=False))
        rv = rng_c.standard_normal((len(present), d_model))
        rv = rv / np.linalg.norm(rv, axis=1, keepdims=True)
        full = int(i) in full_set
        meta_cell = {"cell": g_id, "role": str(role[i]), "class": str(cls[i]), "ntokens": int(L), "registration_sha256": sha,
                     "genes_present": [genes[j][0] for j, _ in present], "full_ladder": full}
        with torch.no_grad():
            h0 = (E[ids] + Pm[:L])[None]  # [1, L, d]
            hs = forward_all(h0)
            zb = torch.stack([h.mean(1)[0] for h in hs])  # [13, d]
            ref_h, ref_lm = model.hidden_states(ids[None], torch.ones(1, L, dtype=torch.long, device=device), return_layer_means=True)
            chk_module = float((ref_lm[:, 0] - zb[1:]).abs().max())
            off = np.stack([np.asarray(emb.layer_means[l_][g_id], np.float32) for l_ in range(emb.n_layers)])
            dif = np.abs(off - zb[1:].cpu().numpy())
            chk_official, chk_official_rel = float(dif.max()), float(dif.max() / (np.abs(off).max() + 1e-12))
        z12 = zb[12].detach()
        y0 = yhat(z12[None])
        g_head = torch.func.grad(lambda z: score(yhat(z[None]))[0])(z12).detach()
        # cases: (gene index or -1, token id, position, direction index, unit vector)
        unit = lambda tid: E[tid] / E[tid].norm()  # noqa: E731
        cases = []
        for k_, (j, t) in enumerate(present):
            cases.append((j, genes[j][2], t, 0, unit(genes[j][2])))
            cases.append((j, genes[j][2], t, 1, f32(rv[k_])))
        for q in other_pos:
            cases.append((-1, int(ids_np[q]), int(q), 2, unit(int(ids_np[q]))))
        for j, t in present:
            if j in nkt_idx:
                cases.append((j, genes[j][2], t, 3, u_nkt))
        K = len(cases)
        sg = torch.stack([E[c[1]].norm() for c in cases])  # [K]
        vdir = torch.stack([c[4] for c in cases])
        pos_k = torch.tensor([c[2] for c in cases], device=device)
        T0 = torch.zeros(K, 1, L, d_model, device=device)
        T0[torch.arange(K, device=device), 0, pos_k] = sg[:, None] * vdir
        J_pool = torch.zeros(K, N_DEPTHS, d_model, device=device)
        J_own = torch.zeros_like(J_pool)
        with torch.no_grad():
            # ---------------- JVP chain
            for c0 in range(0, K, a.jvp_chunk):
                c1 = min(c0 + a.jvp_chunk, K)
                dh = T0[c0:c1]
                pk = pos_k[c0:c1]
                J_pool[c0:c1, 0], J_own[c0:c1, 0] = pooled_and_own(dh, pk)
                for li, layer in enumerate(layers):
                    f = lambda h, _layer=layer: layer_fn(_layer, h)  # noqa: E731
                    dh = vmap(lambda t, _h=hs[li], _f=f: jvp(_f, (_h,), (t,))[1])(dh)
                    J_pool[c0:c1, li + 1], J_own[c0:c1, li + 1] = pooled_and_own(dh, pk)
                del dh
            dY = vmap(lambda t: jvp(yhat, (z12[None],), (t[None],))[1][0])(J_pool[:, 12])  # [K, 4]
            dS = dY[:, :3].mean(-1) - dY[:, 3]
            # ---------------- finite differences through the module's own layers, one batch shape for every sequence.
            # self: registered check (+-1e-3 on the unit tangent; +-1e-2 on full-ladder cells) and the eps_primary rung;
            # nkt: the eps_primary rung; full-ladder cells: every rung for self and nkt; rand / other: JVP only.
            jobs = [(0, ("base",), 0.0)]
            for k, c in enumerate(cases):
                di = c[3]
                if di not in (0, 3):
                    continue
                if di == 0:
                    for ei, e in enumerate(EPS_ABS_CHECK):
                        if ei == 0 or full:
                            jobs += [(k, ("abs", ei, q), sgn * e / float(sg[k])) for q, sgn in enumerate((1, -1))]
                for si, e in enumerate(eps_list):
                    if si == si0 or full:
                        jobs += [(k, ("lad", si, q), cc) for q, cc in enumerate((e, -e, e / 2, -e / 2))]
            nb = a.fd_batch
            while len(jobs) % nb:
                jobs.append(jobs[-1])
            Zd, Yd = {}, {}
            for b0 in range(0, len(jobs), nb):
                bj = jobs[b0:b0 + nb]
                kk = torch.tensor([k for k, _, _ in bj], device=device)
                cc = torch.tensor([c for _, _, c in bj], device=device, dtype=torch.float32)
                zz = forward_module(h0.expand(nb, L, d_model) + cc[:, None, None] * T0[kk, 0])  # [nb, 13, d]
                yy = yhat(zz[:, 12])
                for r_, (k, key_, _) in enumerate(bj):
                    Zd[(k, key_)] = zz[r_]
                    Yd[(k, key_)] = yy[r_]
            z0, yb = Zd[(0, ("base",))], Yd[(0, ("base",))]
            chk_base = float((z0 - zb).abs().max())  # module path vs explicit path, unperturbed
        # statistics (NaN where not run). lin[:, s, j, l]: 0 err_plus ||D+(e) - J|| / ||J||, D+ = (z(+e) - z0) / e;
        # 1 err_minus (D- = (z0 - z(-e)) / e); 2 / 3 slope_plus / minus ||D(e) - D(e/2)|| / ||D(e/2)||; 4 err_central
        # (Dc = (D+ + D-) / 2); 5 slope_central; 6 kappa ||D+ - D-|| / ||D+ + D-||
        S = len(eps_list)
        lin = np.full((K, S, 7, N_DEPTHS), np.nan, np.float32)
        fd_proj = np.full((K, S, 4, N_DEPTHS, 3), np.nan)  # slopes D+(e), D-(e), D+(e/2), D-(e/2) on (d_l, probe, head lens)
        fd_head = np.full((K, S, 4), np.nan)
        jchk = np.full((K, len(EPS_ABS_CHECK), 2, N_DEPTHS), np.nan, np.float32)  # registered check: rel. error, cosine
        nrm = lambda x: x.norm(dim=-1).clamp(min=1e-30)  # noqa: E731
        sc = lambda y: (y[..., :3].mean(-1) - y[..., 3])  # noqa: E731
        s0 = sc(yb)
        for k in range(K):
            J = J_pool[k]
            Jn = nrm(J)
            for si, e in enumerate(eps_list):
                if (k, ("lad", si, 0)) not in Zd:
                    continue
                zp, zm, zph, zmh = (Zd[(k, ("lad", si, q))] for q in range(4))
                Dp, Dm, Dph, Dmh = (zp - z0) / e, (z0 - zm) / e, (zph - z0) / (e / 2), (z0 - zmh) / (e / 2)
                Dc, Dch = (Dp + Dm) / 2, (Dph + Dmh) / 2
                st = [nrm(Dp - J) / Jn, nrm(Dm - J) / Jn, nrm(Dp - Dph) / nrm(Dph), nrm(Dm - Dmh) / nrm(Dmh),
                      nrm(Dc - J) / Jn, nrm(Dc - Dch) / nrm(Dch), nrm(Dp - Dm) / nrm(Dp + Dm)]
                lin[k, si] = torch.stack(st).cpu().numpy()
                for hi, D in enumerate((Dp, Dm, Dph, Dmh)):
                    fd_proj[k, si, hi, :, 0] = (D * dmu).sum(-1).cpu().double().numpy()
                    fd_proj[k, si, hi, :, 1] = (D * wraw).sum(-1).cpu().double().numpy()
                    fd_proj[k, si, hi, :, 2] = (D @ g_head).cpu().double().numpy()
                for hi, (q, sgn, ee) in enumerate(((0, 1, e), (1, -1, e), (2, 1, e / 2), (3, -1, e / 2))):
                    fd_head[k, si, hi] = float(sgn * (sc(Yd[(k, ("lad", si, q))]) - s0) / ee)
            for ei, e in enumerate(EPS_ABS_CHECK):
                if (k, ("abs", ei, 0)) not in Zd:
                    continue
                c = e / float(sg[k])
                Dc = (Zd[(k, ("abs", ei, 0))] - Zd[(k, ("abs", ei, 1))]) / (2 * c)
                jchk[k, ei, 0] = (nrm(Dc - J) / Jn).cpu().numpy()
                jchk[k, ei, 1] = ((Dc * J).sum(-1) / (nrm(Dc) * Jn)).cpu().numpy()
        cpu64 = None
        if n_cpu64 < a.cpu64_check:
            cpu64 = cpu64_jvp_check(layers, h0.detach().cpu().double(), T0[:2].detach().cpu().double(),
                                    J_pool[:2].detach().cpu().double())
            n_cpu64 += 1
        fd64 = None
        if int(i) in fd64_set:
            if layers64 is None:
                layers64 = [copy.deepcopy(l_).cpu().double() for l_ in layers]
            pick = [k for k in (next((k for k, c in enumerate(cases) if c[3] == di), None) for di in range(len(DIRECTIONS)))
                    if k is not None]
            fd64 = fd64_check(layers64, h0.detach().cpu().double(), T0[pick].detach().cpu().double(),
                              sg[pick].detach().cpu().double(), J_pool[pick].detach().cpu().double())
            fd64["cases"] = pick
        meta_cell.update({"check_explicit_vs_module_layer_means_max_abs": chk_module,
                          "check_vs_official_fp16_layer_means_max_abs": chk_official,
                          "check_vs_official_fp16_layer_means_max_rel": chk_official_rel,
                          "check_fd_base_vs_jvp_base_max_abs": chk_base,
                          "cpu64": cpu64, "fd64": fd64, "sec": round(time.time() - tc, 2)})
        atomic_savez(cdir / f"cell_{g_id}.npz",
                     meta=np.array(json.dumps(meta_cell)), n_cases=np.array(K),
                     gene=np.array([c[0] for c in cases], np.int16), token_id=np.array([c[1] for c in cases], np.int32),
                     direction=np.array([c[3] for c in cases], np.int8), pos=np.array([c[2] for c in cases], np.int32),
                     s_gene=sg.cpu().numpy().astype(np.float64), jvp_check=jchk,
                     J_pool=J_pool.cpu().numpy(), J_own=J_own.cpu().numpy(),
                     dY=dY.cpu().double().numpy(), dS=dS.cpu().double().numpy(),
                     g_head=g_head.cpu().numpy(), y0=y0[0].detach().cpu().numpy(), z_base=zb.cpu().numpy(),
                     lin=lin, fd_proj=fd_proj, fd_head=fd_head)
        n_done += 1
        rate = (time.time() - t_start) / n_done
        say(f"  cell {jj + 1}/{len(todo)} ({role[i]}, {cls[i]}, {L} tokens, {K} cases, full {full}) {time.time() - tc:.1f}s; "
            f"module check {chk_module:.1e}"
            + (f", cpu64 rel {cpu64['max_rel_diff_pooled']:.1e}" if cpu64 else "")
            + f"; mean {rate:.1f}s/cell, ETA {rate * (len(todo) - n_done) / 60:.1f} min")
        write_json(out / "e5_progress.json", {"done_this_call": n_done, "to_do_this_call": len(todo),
                                               "sec_per_cell": round(rate, 2), "eta_min": round(rate * (len(todo) - n_done) / 60, 1),
                                               "updated": time.strftime("%Y-%m-%d %H:%M:%S")})
        if device.type == "mps":
            torch.mps.empty_cache()
    say(f"respond done: {n_done} cells in {(time.time() - t_start) / 60:.1f} min")


def fd64_check(layers64, h0, T, sg, J32):
    """experiments_v3 JVP check in float64 on the CPU: central FD through the module's own layers on the unit tangent
    (c = eps_abs / ||E[token]||, c * T = eps_abs * v) against the device float32 JVP; rel. error and cosine per depth."""
    import torch
    t0 = time.time()

    def pooled(x):
        zs = [x.mean(1)]
        h = x
        for l_ in layers64:
            h = l_(h)
            zs.append(h.mean(1))
        return torch.stack(zs, 1)  # [B, 13, d]
    out = np.zeros((T.shape[0], len(EPS_ABS_CHECK), 2, N_DEPTHS))
    with torch.no_grad():
        for ei, e in enumerate(EPS_ABS_CHECK):
            c = (e / sg)[:, None, None]
            x = torch.cat([h0 + c * T[:, 0], h0 - c * T[:, 0]])
            Z = pooled(x)
            n = T.shape[0]
            Dc = (Z[:n] - Z[n:]) / (2 * c)
            Jn = J32.norm(dim=-1).clamp(min=1e-300)
            out[:, ei, 0] = ((Dc - J32).norm(dim=-1) / Jn).numpy()
            out[:, ei, 1] = ((Dc * J32).sum(-1) / (Dc.norm(dim=-1).clamp(min=1e-300) * Jn)).numpy()
    return {"rel_err_cos": out.tolist(), "sec": round(time.time() - t0, 1)}


def cpu64_jvp_check(layers, h0, T, J32):
    """float64 CPU JVP of the same cases; relative difference to the device float32 JVP (pooled, per depth)."""
    import copy

    import torch
    from torch.func import jvp, vmap
    L64 = [copy.deepcopy(l_).cpu().double() for l_ in layers]
    t0 = time.time()
    with torch.no_grad():
        dh = T
        hs = [h0]
        h = h0
        for l_ in L64:
            h = layer_fn(l_, h)
            hs.append(h)
        pooled = [dh[:, 0].mean(1)]
        for li, l_ in enumerate(L64):
            dh = vmap(lambda t, _h=hs[li], _l=l_: jvp(lambda x: layer_fn(_l, x), (_h,), (t,))[1])(dh)
            pooled.append(dh[:, 0].mean(1))
        P = torch.stack(pooled, 1)  # [2, 13, d]
    rel = ((P - J32).norm(dim=-1) / P.norm(dim=-1).clamp(min=1e-300)).numpy()
    return {"max_rel_diff_pooled": float(rel.max()), "rel_diff_by_depth": [float(x) for x in rel.max(0)],
            "sec": round(time.time() - t0, 1)}


# ============================================================================ report

def load_cases(out: Path, sel, role, cls, donors, dirs, gsyms):
    files = sorted((out / "cells").glob("cell_*.npz"))
    want = {int(c) for c in sel}
    keys = ("cell", "role", "cls", "donor", "gene", "direction", "ntok", "pos", "G", "P", "H", "Gown", "Nrel", "cosmu",
            "dS", "dY", "lin", "G_fd", "jchk", "own_norm_gain", "logr", "gsign")
    rows = {k: [] for k in keys}
    metas = []
    dmu, wraw = dirs["dmu"], dirs["wraw"]
    dmu2 = (dmu ** 2).sum(1)
    dmun = np.sqrt(dmu2)
    logit_gap = (wraw * dmu).sum(1)
    head_gap = float(dirs["head_gap_s"])
    head_gap_p = dirs["head_gap_p"]
    rl = {int(c): (r, k) for c, r, k in zip(sel, role, cls)}
    sha = registration_sha()
    for f in files:
        with np.load(f, allow_pickle=False) as d:
            m = json.loads(str(d["meta"]))
            if m["cell"] not in want:
                continue
            if m.get("registration_sha256") != sha:
                raise SystemExit(f"{f} was written under another registration")
            metas.append(m)
            K = int(d["n_cases"])
            if K == 0:
                continue
            n = m["ntokens"]
            A0 = d["s_gene"] / (n * dmun[0])  # [K]
            Jp, Jo = d["J_pool"].astype(np.float64), d["J_own"].astype(np.float64)
            proj = np.einsum("kld,ld->kl", Jp, dmu)
            Jn = np.linalg.norm(Jp, axis=-1)
            own_n = np.linalg.norm(Jo, axis=-1)
            gene = d["gene"].astype(int)
            r, k = rl[m["cell"]]
            add = {"G": proj / dmu2[None] / A0[:, None],
                   "Gown": np.einsum("kld,ld->kl", Jo, dmu) / dmu2[None] / n / A0[:, None],
                   "P": np.einsum("kld,ld->kl", Jp, wraw) / logit_gap[None] / A0[:, None],
                   "H": (Jp @ d["g_head"].astype(np.float64)) / head_gap / A0[:, None],
                   "Nrel": Jn / dmun[None] / A0[:, None],
                   "cosmu": proj / (Jn * dmun[None] + 1e-300),
                   "dS": d["dS"] / head_gap / A0,
                   "dY": d["dY"] / head_gap_p[None] / A0[:, None],
                   "lin": d["lin"], "jchk": d["jvp_check"],
                   "G_fd": d["fd_proj"][..., 0] / dmu2[None, None, None] / A0[:, None, None, None],  # [K, S, 4, 13]
                   "own_norm_gain": own_n / own_n[:, :1],
                   "logr": np.log(np.abs(proj / dmun[None]) + 1e-30),
                   "gsign": np.array([{"NK": 1.0, "T": -1.0}.get(E2_NKT_CLASS.get(gsyms[g]), 0.0) if g >= 0 else 0.0 for g in gene]),
                   "gene": gene, "direction": d["direction"].astype(int), "pos": d["pos"]}
            for kk, v in add.items():
                rows[kk].append(v)
            rows["cell"] += [m["cell"]] * K
            rows["role"] += [r] * K
            rows["cls"] += [k] * K
            rows["donor"] += [donors[m["cell"]]] * K
            rows["ntok"] += [n] * K
    C = {}
    for kk, v in rows.items():
        C[kk] = (np.concatenate(v) if isinstance(v[0], np.ndarray) else np.array(v)) if v else np.zeros(0)
    return C, metas


class Boot:
    """Resampling of cells (all of a cell's cases kept): 'cell' = cells with replacement within each stratum;
    'donor' = two-stage, donors with replacement, then cells with replacement within each drawn donor and stratum."""

    def __init__(self, cell, donor, stratum, n_boot: int, seed: int):
        o = np.argsort(cell, kind="stable")
        cs, start = np.unique(cell[o], return_index=True)
        self.groups = np.split(o, start[1:])
        first = o[start]
        cd, cst = donor[first], stratum[first]
        self.donors = sorted(set(cd.tolist()))
        rng = np.random.default_rng(seed)
        strata = [np.where(cst == s)[0] for s in np.unique(cst)]
        bydon = {d_: [np.where((cd == d_) & (cst == s))[0] for s in np.unique(cst[cd == d_])] for d_ in self.donors}
        self.samples = {"cell": [], "donor": []}
        for _ in range(n_boot):
            self.samples["cell"].append(np.concatenate([g[rng.integers(0, g.size, g.size)] for g in strata]))
        for _ in range(n_boot):
            ds = rng.integers(0, len(self.donors), len(self.donors))
            self.samples["donor"].append(np.concatenate([g[rng.integers(0, g.size, g.size)]
                                                         for x in ds for g in bydon[self.donors[x]]]))
        self.n = sum(g.size for g in self.groups)

    def cases(self, cell_sample):
        return np.concatenate([self.groups[g] for g in cell_sample])


def pct(v):
    v = np.asarray(v, np.float64)
    return np.percentile(v, [2.5, 97.5], axis=0)


def subset_summary(X: dict, donor: np.ndarray, boot: Boot) -> dict:
    """X: name -> [cases, D] per-case quantities of one subset. Median over cases with two-stage ('donor') and
    cell-cluster ('cell') 95% intervals, per-donor medians; for 'dG' also the bootstrap share of each argmin."""
    names = list(X)
    widths = [X[k].shape[1] for k in names]
    M = np.concatenate([X[k] for k in names], 1)
    full = np.median(M, 0)
    res = {}
    boots = {}
    for lev in ("donor", "cell"):
        B = np.array([np.median(M[boot.cases(s)], 0) for s in boot.samples[lev]])
        boots[lev] = B
    per_d = {d_: np.median(M[donor == d_], 0) for d_ in sorted(set(donor.tolist()))}
    o = 0
    for k, w in zip(names, widths):
        sl = slice(o, o + w)
        res[k] = {"median": full[sl].tolist(),
                  "ci_two_stage": pct(boots["donor"][:, sl]).T.tolist(),
                  "ci_cell": pct(boots["cell"][:, sl]).T.tolist(),
                  "per_donor": {d_: v[sl].tolist() for d_, v in per_d.items()}}
        if k in ("dG", "logg"):
            for lev in ("donor", "cell"):
                A = np.argmin(boots[lev][:, sl], 1) + 1
                res[k][f"argmin_share_{'two_stage' if lev == 'donor' else 'cell'}"] = {int(v): float(np.mean(A == v)) for v in np.unique(A)}
        o += w
    return res


def joint_diff(xa, xb, cell, donor, stratum, boot_seed, n_boot):
    """median(xa over stratum 'a') - median(xb over stratum 'b') with a joint two-stage bootstrap (one donor draw,
    cells resampled within donor and stratum) and per-donor values. xa, xb: per-case values over the union."""
    b = Boot(cell, donor, stratum, n_boot, boot_seed)
    A, Bm = stratum == "a", stratum == "b"
    full = float(np.median(xa[A]) - np.median(xb[Bm]))

    def st(ix):
        ia, ib = ix[A[ix]], ix[Bm[ix]]
        return np.median(xa[ia]) - np.median(xb[ib]) if ia.size and ib.size else np.nan
    out = {"diff": full}
    for lev, nm in (("donor", "ci_two_stage"), ("cell", "ci_cell")):
        v = np.array([st(b.cases(s)) for s in b.samples[lev]])
        v = v[np.isfinite(v)]
        out[nm] = pct(v).tolist() if v.size else None
    out["per_donor"] = {d_: float(np.median(xa[A & (donor == d_)]) - np.median(xb[Bm & (donor == d_)]))
                        for d_ in sorted(set(donor.tolist())) if (A & (donor == d_)).any() and (Bm & (donor == d_)).any()}
    return out


def specificity_boot(LG, cell, donor, stratum, lstar, seed, n_boot):
    """Registered E5 specificity D = median log g_{l*}(self) - median log g_{l*}(control). LG: [cases, 12] log-gains
    over the union of self ('a') and control ('b') cases. Two-stage bootstrap with l* recomputed inside each replicate
    (amendment A1, common.statistics.bootstrap); per-donor D at the full-sample l*."""
    b = Boot(cell, donor, stratum, n_boot, seed)
    A, Bm = stratum == "a", stratum == "b"

    def st(ix):
        ia, ib = ix[A[ix]], ix[Bm[ix]]
        if not ia.size or not ib.size:
            return np.nan, -1
        ls = int(np.argmin(np.median(LG[ia], 0)))
        return float(np.median(LG[ia, ls]) - np.median(LG[ib, ls])), ls + 1
    full = float(np.median(LG[A, lstar - 1]) - np.median(LG[Bm, lstar - 1]))
    out = {"diff": full, "l_star": lstar}
    for lev, nm in (("donor", "ci_two_stage"), ("cell", "ci_cell")):
        v = [st(b.cases(s_)) for s_ in b.samples[lev]]
        d = np.array([x[0] for x in v])
        ok = np.isfinite(d)
        out[nm] = pct(d[ok]).tolist() if ok.any() else None
        ls = np.array([x[1] for x in v])[ok]
        out[f"l_star_share_{'two_stage' if lev == 'donor' else 'cell'}"] = {int(x): float(np.mean(ls == x)) for x in np.unique(ls)}
    out["per_donor"] = {d_: float(np.median(LG[A & (donor == d_), lstar - 1]) - np.median(LG[Bm & (donor == d_), lstar - 1]))
                        for d_ in sorted(set(donor.tolist())) if (A & (donor == d_)).any() and (Bm & (donor == d_)).any()}
    return out


def probe_curves(probe_dir: Path) -> dict:
    R = json.loads((probe_dir / "mode_a_results.json").read_text())
    by = {l_["name"]: l_ for l_ in R["layers"]}
    names = ["input"] + [f"layer_{i}" for i in range(1, 13)]
    cur = {"depth_names": names, "order_rate_no_gdT158": [], "auc": [], "gap_ratio_mean_no_gdT158": [],
           "gap_ratio_mean_all": [], "gap_ratio_no_gdT158": {p: [] for p in PROTEINS}, "centred_distance_ratio_no_gdT158": []}
    for nm in names:
        l_ = by[nm]
        cur["order_rate_no_gdT158"].append(l_["nk_t_probe"]["pairs"]["no_gdT158"]["order_rate"])
        cur["auc"].append(l_["nk_t_probe"]["site4_auc"])
        g = l_["protein_probe"]["gap"]["no_gdT158"]
        rs = [g[p]["ratio"] for p in PROTEINS]
        for p, r_ in zip(PROTEINS, rs):
            cur["gap_ratio_no_gdT158"][p].append(r_)
        cur["gap_ratio_mean_no_gdT158"].append(float(np.mean(rs)))
        cur["gap_ratio_mean_all"].append(float(np.mean([l_["protein_probe"]["gap"]["all"][p]["ratio"] for p in PROTEINS])))
        cur["centred_distance_ratio_no_gdT158"].append(l_["cosine"]["centred"]["no_gdT158"]["distance_ratio_neighbour_over_random"])
    gm = np.array(cur["gap_ratio_mean_no_gdT158"])
    cur["L_probe_argmin_gap_ratio"] = int(np.argmin(gm[1:]) + 1)
    cur["L_probe_largest_drop_gap_ratio"] = int(np.argmin(np.diff(gm)) + 1)
    cur["L_probe_largest_drop_gap_ratio_all"] = int(np.argmin(np.diff(np.array(cur["gap_ratio_mean_all"]))) + 1)
    cur["L_probe_largest_drop_distance_ratio"] = int(np.argmin(np.diff(np.array(cur["centred_distance_ratio_no_gdT158"], float))) + 1)
    cur["source"] = str(probe_dir / "mode_a_results.json")
    tp = probe_dir / "token" / "mode_a_token_results.json"
    if tp.exists():
        T = json.loads(tp.read_text())
        cur["token_probes_ridge_gap_ratio_no_gdT158"] = {dn: {fs: {p: T["feature_sets"][dn][fs]["ridge"]["gap_no_gdT158"][p]["ratio"]
                                                                  for p in PROTEINS} for fs in ("genemean", "tokens")}
                                                         for dn in T["feature_sets"]}
    return cur


def spearman(a, b):
    from scipy.stats import spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[ok], b[ok]).statistic) if ok.sum() >= 3 else None


def _all_donors(per_donor: dict, fn) -> bool:
    return bool(per_donor) and all(fn(v) for v in per_donor.values())


def stage_report(a, meta, reg_info, reg):
    out = a.out_dir
    sel, role, _, cls = load_cells(out)
    donors = load_donors(a.processed)
    dirs = load_directions(out)
    dinfo = json.loads((out / "e5_directions.json").read_text())
    ginfo = json.loads((out / "e5_genes.json").read_text())
    gsyms = list(ginfo["kept"])
    C, metas = load_cases(out, sel, role, cls, donors, dirs, gsyms)
    tau = 0.1
    eps = list(REGISTRATION_EPS)
    si0 = eps.index(EPS_PRIMARY)
    seed = int(reg["seeds"]["bootstrap"])
    margins = reg["experiments"]["E5"]["margins"]
    R = {"experiment": "E5", "smoke_test": bool(a.smoke), "smoke_note": a.smoke, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
         "script_version": SCRIPT_VERSION, "registration_sha256": reg_info.get("registration_sha256"),
         "addendum_sha256": reg_info["addendum_sha256"], "registration_check": reg_info, "cell_pool": a.cell_pool,
         "n_boot": a.n_boot, "n_cells_selected": int(sel.size), "n_cells_done": len(metas),
         "complete": len(metas) == sel.size, "n_cases": int(C["G"].shape[0]) if C["G"].size else 0,
         "genes": ginfo, "directions_info": dinfo}
    if not C["G"].size:
        write_json(out / "E5_results.json", R)
        say("no cases yet")
        return
    roles = [r for r in ("pair", "random", "val") if np.any(C["role"] == r)]
    main = "pair" if "pair" in roles else roles[0]
    R["cases_by_role_direction"] = {f"{r}/{d}": int(np.sum((C["role"] == r) & (C["direction"] == di)))
                                    for r in roles for di, d in enumerate(DIRECTIONS)}
    R["checks"] = {
        "explicit_vs_module_layer_means_max_abs": max(m["check_explicit_vs_module_layer_means_max_abs"] for m in metas),
        "vs_official_fp16_layer_means_max_rel": max(m["check_vs_official_fp16_layer_means_max_rel"] for m in metas),
        "fd_base_module_vs_explicit_max_abs": max(m["check_fd_base_vs_jvp_base_max_abs"] for m in metas),
        "head_lens_12_vs_head_jvp_max_abs": float(np.max(np.abs(C["H"][:, 12] - C["dS"]))),
        "cells_without_e2_genes": int(sum(1 for m in metas if not m["genes_present"])),
        "sec_per_cell_median": float(np.median([m["sec"] for m in metas])),
    }
    # ---------------- registered JVP check: float32 central FD on the self cases
    jc, selfm = C["jchk"], C["direction"] == 0
    chk = {"scope": "self cases (cell, gene); every depth 1..12; float32", "eps_abs": list(EPS_ABS_CHECK)}
    for ei, e in enumerate(EPS_ABS_CHECK):
        h = selfm & np.isfinite(jc[:, ei, 0, 1])
        ok = (jc[h, ei, 0, 1:].max(1) <= 0.05) & (jc[h, ei, 1, 1:].min(1) >= 0.99)
        chk[str(e)] = {"n_cases": int(h.sum()), "pass_share": float(ok.mean()) if h.any() else None,
                       "median_rel_err_by_depth": np.median(jc[h, ei, 0], 0).tolist() if h.any() else None,
                       "p05_cos_by_depth": np.percentile(jc[h, ei, 1], 5, axis=0).tolist() if h.any() else None,
                       "pass_share_by_depth": ((jc[h, ei, 0] <= 0.05) & (jc[h, ei, 1] >= 0.99)).mean(0).tolist() if h.any() else None}
    chk["validated"] = bool(chk[str(EPS_ABS_CHECK[0])]["pass_share"] is not None and chk[str(EPS_ABS_CHECK[0])]["pass_share"] >= 0.95)
    R["registered_jvp_check"] = chk
    # ---------------- S5 / S6 float64 checks
    f64 = [np.array(m["fd64"]["rel_err_cos"]) for m in metas if m.get("fd64")]
    if f64:
        F = np.concatenate(f64)
        fok = [(F[:, ei, 0, 1:].max(1) <= 0.05) & (F[:, ei, 1, 1:].min(1) >= 0.99) for ei in range(len(EPS_ABS_CHECK))]
        R["S5_float64_check"] = {"n_cases": int(F.shape[0]), "pass_share": {str(e): float(fok[ei].mean()) for ei, e in enumerate(EPS_ABS_CHECK)},
                                 "max_rel_err_by_depth": {str(e): F[:, ei, 0].max(0).tolist() for ei, e in enumerate(EPS_ABS_CHECK)},
                                 "min_cos_by_depth": {str(e): F[:, ei, 1].min(0).tolist() for ei, e in enumerate(EPS_ABS_CHECK)},
                                 "would_validate": bool(fok[0].mean() >= 0.95)}
    R["S6_float64_jvp_max_rel_diff"] = max((m["cpu64"]["max_rel_diff_pooled"] for m in metas if m.get("cpu64")), default=None)
    # ---------------- S4 linearity (self and nkt cases)
    LS = C["lin"]
    has = np.isfinite(LS[:, :, 0, 1])
    os_pass = np.all(LS[:, :, :4] <= tau, axis=2)
    c_pass = np.all(LS[:, :, 4:6] <= tau, axis=2)
    names = ("err_plus", "err_minus", "slope_plus", "slope_minus", "err_central", "slope_central", "kappa")
    lin = {}
    for si, e in enumerate(eps):
        h = has[:, si]
        lin[str(e)] = None if not h.any() else {
            "n_cases": int(h.sum()), "one_sided_pass_rate_by_depth": os_pass[h, si].mean(0).tolist(),
            "central_pass_rate_by_depth": c_pass[h, si].mean(0).tolist(),
            "median_by_depth": {nm: np.median(LS[h, si, j], 0).tolist() for j, nm in enumerate(names)},
            "one_sided_pass_rate_by_direction_depth": {d: os_pass[h & (C["direction"] == di), si].mean(0).tolist()
                                                       for di, d in enumerate(DIRECTIONS) if np.any(h & (C["direction"] == di))}}
    h1 = np.array(lin[str(eps[si0])]["one_sided_pass_rate_by_depth"][1:]) if lin[str(eps[si0])] else np.zeros(12)
    ok_eps = [e for e in eps if lin[str(e)] and np.all(np.array(lin[str(e)]["one_sided_pass_rate_by_depth"][1:]) >= 0.90)]
    R["S4_linearity"] = {"by_eps": lin, "holds_at_eps_primary": bool(np.all(h1 >= 0.90)), "eps_primary": EPS_PRIMARY,
                         "min_one_sided_pass_rate_eps_primary": float(h1.min()),
                         "linear_range_largest_eps": max(ok_eps) if ok_eps else None}
    Gp, Gm = C["G_fd"][:, :, 0], C["G_fd"][:, :, 1]
    # ---------------- per subset summaries
    curves, idx = {}, {}
    for r in roles:
        for di, d in enumerate(DIRECTIONS):
            m_ = (C["role"] == r) & (C["direction"] == di)
            if d == "self":
                m_ &= C["gsign"] != 0  # sign-aligned G only for the 13 NK/T genes; the log-gain uses every E2 gene
            ix = np.where(m_)[0]
            ixl = np.where((C["role"] == r) & (C["direction"] == di))[0]
            if ixl.size >= 2:
                idx[(r, d)] = ixl
            if ix.size < 2:
                continue
            sg = C["gsign"][ix][:, None] if d == "self" else np.ones((ix.size, 1))
            G = C["G"][ix] * sg
            X = {"G": G, "dG": np.diff(G, axis=1), "P": C["P"][ix] * sg, "H": C["H"][ix] * sg,
                 "G_own": C["Gown"][ix] * sg, "norm_rel": C["Nrel"][ix], "cos_with_d": C["cosmu"][ix] * sg,
                 "head_minus_z_12": (C["H"][ix, 12:13] - C["G"][ix, 12:13]) * sg}
            b = Boot(C["cell"][ix], C["donor"][ix], np.zeros(ix.size, int), a.n_boot, seed)
            blk = {"n_cases": int(ix.size), "n_cells": int(np.unique(C["cell"][ix]).size), **subset_summary(X, C["donor"][ix], b)}
            hx = has[ix]
            blk["G_fd_plus_by_eps"] = {str(e): np.median((Gp[ix, si] * sg)[hx[:, si]], 0).tolist() if hx[:, si].any() else None
                                       for si, e in enumerate(eps)}
            blk["G_fd_minus_by_eps"] = {str(e): np.median((Gm[ix, si] * sg)[hx[:, si]], 0).tolist() if hx[:, si].any() else None
                                        for si, e in enumerate(eps)}
            if d == "self":
                blk["sign_consistency_by_depth"] = (G > 0).mean(0).tolist()
            blk["by_gene"] = {gsyms[g]: {"n_cases": int((C["gene"][ix] == g).sum()), "G_12_median": float(np.median(G[C["gene"][ix] == g, 12])),
                                         "median_rank_position": float(np.median(C["pos"][ix][C["gene"][ix] == g]))}
                              for g in np.unique(C["gene"][ix]) if g >= 0}
            curves[f"{r}/{d}"] = blk
    R["curves"] = curves
    # log-gain curves (registered quantity) for every subset, all E2 genes
    lg_curves = {}
    for (r, d), ixl in idx.items():
        X = {"logg": np.diff(C["logr"][ixl], axis=1)}
        b = Boot(C["cell"][ixl], C["donor"][ixl], np.zeros(ixl.size, int), a.n_boot, seed)
        lg_curves[f"{r}/{d}"] = {"n_cases": int(ixl.size), **subset_summary(X, C["donor"][ixl], b)}
    R["log_gain_curves"] = lg_curves
    PC = probe_curves(a.probe_dir)
    R["probe_curves"] = PC
    # ---------------- registered E5 endpoint
    E5 = {"cells": main, "jvp_check_validated": chk["validated"]}
    ls = lg_curves.get(f"{main}/self")
    if ls:
        med = np.array(ls["logg"]["median"])
        lstar = int(np.argmin(med) + 1)
        Lc = PC["L_probe_largest_drop_gap_ratio"]
        E5.update({"l_star": lstar, "median_log_gain_at_l_star": float(med[lstar - 1]),
                   "l_star_bootstrap_share_two_stage": ls["logg"]["argmin_share_two_stage"].get(lstar, 0.0),
                   "comparator_layer_no_gdT158": Lc, "comparator_layer_all": PC["L_probe_largest_drop_gap_ratio_all"],
                   "agreement_within_margin": abs(lstar - Lc) <= int(margins["localisation_agreement_layers"])})
        mg = float(margins["direction_specificity_log_gain"])
        for ctl in ("rand", "other"):
            if (main, ctl) in idx:
                ia, ib = idx[(main, "self")], idx[(main, ctl)]
                u = np.concatenate([ia, ib])
                st = np.array(["a"] * ia.size + ["b"] * ib.size)
                dd = specificity_boot(np.diff(C["logr"][u], axis=1), C["cell"][u], C["donor"][u], st, lstar, seed, a.n_boot)
                c = dd["ci_two_stage"]
                if c and -mg < c[0] and c[1] < mg:
                    v_ = "equivalent: no layer localises the NK-T loss"
                elif c and dd["diff"] <= -mg and c[1] < 0 and _all_donors(dd["per_donor"], lambda v: v <= -mg / 2):
                    v_ = f"localised: the NK-T component loses more than the control at layer {lstar}"
                else:
                    v_ = "inconclusive"
                E5[f"specificity_vs_{ctl}"] = {**dd, "verdict": v_}
    if not chk["validated"]:
        ver = "JVP not validated (registered check failed): E5 stops; registered reading: a linear-response description is not valid at this scale"
        if R.get("S5_float64_check", {}).get("would_validate"):
            ver += ". Note: the float64 version of the same check passes (S5), so the failure is float32 rounding; see the addendum's infeasibilities"
    elif "specificity_vs_rand" in E5:
        sv = E5["specificity_vs_rand"]["verdict"]
        ver = ("no layer localises the NK-T loss (registered falsification)" if sv.startswith("equivalent")
               else f"NK-T loss localised at layer {E5['l_star']} (probe-curve layer {E5['comparator_layer_no_gdT158']}, "
                    f"agreement within 1 layer: {E5['agreement_within_margin']})" if sv.startswith("localised") else "inconclusive")
    else:
        ver = "no self / rand cases"
    E5["verdict"] = ver if a.cell_pool == "site4" else f"(smoke on {a.cell_pool} cells; not a registered result) {ver}"
    R["registered_endpoint"] = E5
    # ---------------- secondary S1-S3 (nkt pushes)
    sec = {}
    cm = curves.get(f"{main}/nkt")
    if cm:
        dG = cm["dG"]
        L_e5 = int(np.argmin(dG["median"]) + 1)
        ub = dG["ci_two_stage"][L_e5 - 1][1]
        share = dG["argmin_share_two_stage"].get(L_e5, 0.0)
        pdn = {d_: v[L_e5 - 1] for d_, v in dG["per_donor"].items()}
        sec["S1_localised_nkt_attenuation"] = {"holds": bool(ub < 0 and share >= 0.5 and _all_donors(pdn, lambda v: v < 0)),
                                               "L_E5": L_e5, "median_dG": dG["median"][L_e5 - 1], "ci_two_stage": dG["ci_two_stage"][L_e5 - 1],
                                               "argmin_share_two_stage": share, "per_donor_dG": pdn,
                                               "L_probe_argmin_gap_ratio": PC["L_probe_argmin_gap_ratio"],
                                               "spearman_G_vs_probe_gap_ratio": spearman(cm["G"]["median"], PC["gap_ratio_mean_no_gdT158"])}
        if ("pair", "nkt") in idx and ("random", "nkt") in idx:
            ia, ib = idx[("pair", "nkt")], idx[("random", "nkt")]
            u = np.concatenate([ia, ib])
            st = np.array(["a"] * ia.size + ["b"] * ib.size)
            dd = joint_diff(C["G"][u, 12], C["G"][u, 12], C["cell"][u], C["donor"][u], st, seed, a.n_boot)
            c = dd["ci_two_stage"]
            sec["S2_lookalike_specificity"] = {**dd, "holds": bool(c and c[1] < 0 and _all_donors(dd["per_donor"], lambda v: v < 0)),
                                               "opposite": bool(c and c[0] > 0 and _all_donors(dd["per_donor"], lambda v: v > 0))}
        hm = cm["head_minus_z_12"]
        pdh = {d_: v[0] for d_, v in hm["per_donor"].items()}
        sec["S3_head_transmission"] = {"median_H12_minus_G12": hm["median"][0], "ci_two_stage": hm["ci_two_stage"][0],
                                       "per_donor": pdh, "head_passes_less": bool(hm["ci_two_stage"][0][1] < 0 and _all_donors(pdh, lambda v: v < 0))}
    R["secondary"] = sec
    write_json(out / "E5_results.json", R)
    write_report(out, R)
    say(f"wrote {out / 'E5_results.json'} and REPORT.md ({R['n_cases']} cases, {len(metas)} cells)")


def _f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def _ci(c, n=2):
    return "" if c is None else f" [{_f(c[0], n)}, {_f(c[1], n)}]"


def write_report(out: Path, R: dict):
    L = []
    eps = list(REGISTRATION_EPS)
    PC = R["probe_curves"]
    E5 = R["registered_endpoint"]
    chk = R["registered_jvp_check"]
    flag = f"SMOKE TEST: {R['smoke_note']}\n\n" if R["smoke_test"] else ""
    ck = R["checks"]
    s5 = R.get("S5_float64_check")
    L += ["# E5 Mode A layer dynamics (C7)", "", flag +
          f"{SCRIPT_VERSION}; cells {R['n_cells_done']} of {R['n_cells_selected']} ({R['cell_pool']}); cases {R['n_cases']}; "
          f"registration_sha256 {R['registration_sha256']}; addendum_sha256 {R['addendum_sha256']}.", "",
          f"**Registered verdict:** {E5['verdict']}", "",
          "## Registered JVP check (float32 central FD on the unit tangent; self cases; every depth)", "",
          "| eps | cases | pass share (rel. err <= 0.05 and cos >= 0.99) |", "|---|---|---|"]
    for e in EPS_ABS_CHECK:
        c = chk[str(e)]
        L.append(f"| {e} | {c['n_cases']} | {_f(c['pass_share'])} |")
    L += ["", f"Validated (>= 0.95 at 1e-3): {chk['validated']}. S5 float64 version: "
          + (", ".join(f"{e}: {_f(v)}" for e, v in s5["pass_share"].items()) + f" of {s5['n_cases']} cases" if s5 else "not run")
          + f". S6 float64 JVP vs device JVP max rel. {_f(R['S6_float64_jvp_max_rel_diff'], 7) if R['S6_float64_jvp_max_rel_diff'] is not None else 'n/a'}. "
          f"Explicit layers vs module forward max |diff| {ck['explicit_vs_module_layer_means_max_abs']:.1e}; vs the official fp16 "
          f"layer means max rel. {ck['vs_official_fp16_layer_means_max_rel']:.1e}; median {_f(ck['sec_per_cell_median'], 1)} s per cell.", ""]
    lg = R["log_gain_curves"]
    main = E5["cells"]
    if f"{main}/self" in lg:
        L += [f"## Registered endpoint: median log-gain of the NK-T component per layer ({main} cells)", "",
              "log g_l = log|r_l| - log|r_(l-1)|, r_l = <dz_l, unit(d_l)>. Intervals: two-stage bootstrap.", "",
              "| layer | self (E2 genes) | rand (same positions) | other genes | nkt (secondary) | probe gap-ratio step (no_gdT158) |",
              "|---|---|---|---|---|---|"]
        gm = np.array(PC["gap_ratio_mean_no_gdT158"])
        for l_ in range(1, N_DEPTHS):
            row = [f"{_f(lg[f'{main}/{d}']['logg']['median'][l_ - 1])}" + (_ci(lg[f'{main}/{d}']['logg']['ci_two_stage'][l_ - 1]) if d == "self" else "")
                   if f"{main}/{d}" in lg else "" for d in ("self", "rand", "other", "nkt")]
            L.append(f"| {l_} | " + " | ".join(row) + f" | {_f(float(gm[l_] - gm[l_ - 1]))} |")
        L += ["", "Endpoint: `" + json.dumps({k: v for k, v in E5.items() if not k.startswith("specificity")}, default=float) + "`", ""]
        for k in ("specificity_vs_rand", "specificity_vs_other"):
            if k in E5:
                L.append(f"- {k}: D = {_f(E5[k]['diff'])}{_ci(E5[k]['ci_two_stage'])}, per donor "
                         + ", ".join(f"{d_} {_f(v)}" for d_, v in E5[k]["per_donor"].items()) + f" -> {E5[k]['verdict']}")
    cm, cr = R["curves"].get(f"{main}/nkt"), R["curves"].get("random/nkt")
    S4 = R["S4_linearity"]
    if cm:
        L += ["", f"## Secondary: NK-T push (nkt) through the layers ({main} cells) next to the static probes", "",
              "G = NK-T component (fraction of the class-mean distance per unit of an aligned input push; G_0 = 1); dG = median "
              "step change; H = head lens (depth 12 = the head's own response); pass = one-sided linearity pass rate (self and "
              f"nkt cases). Linear range: eps <= {S4['linear_range_largest_eps']}.", "",
              "| depth | G | dG | G (random cells) | H | " + " | ".join(f"pass {e}" for e in eps)
              + " | probe gap ratio | probe order rate | centred dist. ratio |",
              "|---|---|---|---|---|" + "---|" * len(eps) + "---|---|---|"]
        g, dG = cm["G"], cm["dG"]
        for l_ in range(N_DEPTHS):
            L.append(f"| {PC['depth_names'][l_]} | {_f(g['median'][l_])}{_ci(g['ci_two_stage'][l_])} | "
                     + (f"{_f(dG['median'][l_ - 1])}{_ci(dG['ci_two_stage'][l_ - 1])}" if l_ else "")
                     + f" | {_f(cr['G']['median'][l_]) if cr else ''} | {_f(cm['H']['median'][l_])} | "
                     + " | ".join(_f(S4["by_eps"][str(e)]["one_sided_pass_rate_by_depth"][l_], 2) if S4["by_eps"][str(e)] else "n/a" for e in eps)
                     + f" | {_f(PC['gap_ratio_mean_no_gdT158'][l_])} | {_f(PC['order_rate_no_gdT158'][l_])} | "
                     f"{_f(PC['centred_distance_ratio_no_gdT158'][l_])} |")
        L += ["", f"Gain vs push size ({main}/nkt, median G at depths 3, 6, 9, 12):", "", "| push | G_3 | G_6 | G_9 | G_12 |", "|---|---|---|---|---|",
              "| JVP | " + " | ".join(_f(g["median"][l_]) for l_ in (3, 6, 9, 12)) + " |"]
        for e in eps:
            for nm, k_ in (("+", "G_fd_plus_by_eps"), ("-", "G_fd_minus_by_eps")):
                v = cm[k_][str(e)]
                L.append(f"| {nm}{e} | " + (" | ".join(_f(v[l_]) for l_ in (3, 6, 9, 12)) if v else "n/a | | |") + " |")
    L += ["", "## Secondary hypotheses (addendum; never change the registered verdict)", ""]
    for k, v in R["secondary"].items():
        L.append(f"- **{k}**: `" + json.dumps(v, default=float)[:1200] + "`")
    L.append(f"- **S4_linearity**: holds at eps {S4['eps_primary']}: {S4['holds_at_eps_primary']} (min pass {_f(S4['min_one_sided_pass_rate_eps_primary'])})")
    L += ["", "## All subsets (median G at depths 0, 6, 9, 12; self = the 13 NK/T genes, sign-aligned)", "",
          "| subset | cases | cells | G_0 | G_6 | G_9 | G_12 | H_12 |", "|---|---|---|---|---|---|---|---|"]
    for key_, blk in R["curves"].items():
        gm_ = blk["G"]["median"]
        L.append(f"| {key_} | {blk['n_cases']} | {blk['n_cells']} | {_f(gm_[0])} | {_f(gm_[6])} | {_f(gm_[9])} | {_f(gm_[12])} | "
                 f"{_f(blk['H']['median'][12])} |")
    L += ["", f"Genes: {len(R['genes']['kept'])} E2 genes kept; dropped: {R['genes']['dropped']}."]
    (out / "REPORT.md").write_text("\n".join(L) + "\n")


# ============================================================================ main

def main(argv=None) -> int:
    a = parse_args(argv)
    if a.stage == "register":
        stage_register(a)
        return 0
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG["path"] = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage} out {a.out_dir}" + (f" SMOKE: {a.smoke}" if a.smoke else ""))
    reg_info = check_registration(a)
    write_json(a.out_dir / "e5_registration_used.json", {**reg_info, "registration": REGISTRATION})
    need_rna = a.stage in ("directions", "respond", "all")
    meta = lp.load_meta(a.processed, need_rna=need_rna)
    reg = load_reg(a)
    key = primary_keys(meta, reg)
    if a.stage in ("select", "directions", "respond", "all"):
        emb = lp.Embedding(a.embed_dir, None, len(meta["split"]))
    if a.stage in ("select", "all"):
        if (a.out_dir / "e5_cells.npz").exists() and a.stage == "all":
            say("reusing e5_cells.npz (the selection is deterministic; delete the file to rebuild it)")
        else:
            stage_select(a, meta, emb, key, reg)
    if a.stage in ("directions", "all"):
        dj = a.out_dir / "e5_directions.json"
        if (a.out_dir / "e5_directions.npz").exists() and dj.exists() and a.stage == "all":
            if json.loads(dj.read_text()).get("subsampled") and not a.smoke:
                raise SystemExit("e5_directions.npz in this out-dir was built from a smoke subsample; delete it first")
            say("reusing e5_directions.npz")
        else:
            stage_directions(a, meta, emb, key)
    if a.stage in ("respond", "all"):
        stage_respond(a, meta, emb, load_directions(a.out_dir), reg)
    if a.stage in ("report", "all"):
        stage_report(a, meta, reg_info, reg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
