#!/usr/bin/env python3
"""E2 (C3 response decomposition, C5 linearity): perturb the RNA of the registered E2 cells, re-embed with
the official TEDDY-G preprocessing, and follow each perturbation through the frozen pipeline
(z -> phase-1 head -> registered evidence -> Q1 rule and ANM) to say where a lost call was lost.

Registered design: registration/registration_v3.json experiments.E2, as amended by A1 (A1.7 falsification,
A1.1 bootstrap), loaded with bridge_anm/lib/v3_amend.load_registration_amended(). This script's addendum
registration/addenda/E2.json fixes only what the registration leaves open (the z-probe C, the gene ids, the
count recovery and random-number scheme, unit and statistic definitions; see SPEC below) and is committed
before any site4 forward pass.

  cells        registration e2_subset: 1,000 test_primary cells (500 per donor 13272, 19593) + 250 of donor
               15078 (secondary); smoke runs use val donor 18303 cells only (--cell-pool val)
  perturb      floor (one unperturbed re-embedding); RNA binomial thinning to 0.95/0.9/0.8/0.6/0.2/0.05 of the
               UMIs, seeds 11 and 12; each E2 gene among the cell's tokens scaled by 1 - eps, eps 1.0/0.5/0.25
  measure      dz, dv (12 panel evidences), dS and dmargin, dANM (ANM engine), dcall
  linearity    (eps, eps/2) slope ratio and cosine per stage (z, evidence, score), floor-guarded
  decompose    lost call -> representation (z-probe wrong on z') / head (rule argmax wrong) / decision
  comparator   the plain perturbation -> accuracy curve (Q1 selective and decision accuracy)
  falsified    unless (i) or (ii) of A1's E2 text holds (two-stage donor bootstrap, B = 2000, seed 1)

Stages
  register     train/val only: gene ids, z-probe C on val, count-recovery check, leakage self-check (site4
               rows poisoned -> identical addendum core; val rows poisoned -> changed); writes
               registration/addenda/E2.json, its line in addenda/HASHES.txt and <out-dir>/z_probe.npz
  embed        GPU (MPS) fp16: every variant of every cell, chunk files <out-dir>/chunks/chunk_NNNN.npz
               (resumable: finished chunks are skipped; progress and ETA in progress.log, e2_progress.json)
  report       CPU: head, evidence, rule, ANM engine, z-probe, linearity, decomposition, curve, falsification,
               bootstrap; writes E2_results.json, REPORT.md, e2_cases.npz
  all          embed then report
A site4 run (embed / report / all on --cell-pool e2_subset) is refused unless registration_v3.json,
amendment_A1.json and addenda/E2.json are committed, hash-checked and equal to this script's SPEC.

Full run (not run inside the workflow):
  R=/Users/tianchichen/Documents/GitHub/teddy_mm; PY=<venv>/bin/python
  ANM_ROOT=<anm_v2_fix> $PY scripts/e2_response_decomposition.py --stage all --cell-pool e2_subset \\
      --processed $R/data/processed/cite --embed-dir $R/data/processed/cite_official \\
      --ckpt /Users/tianchichen/Documents/GitHub/teddy_mwe/ckpt/teddy_g_70M \\
      --medians $R/data/reference/teddy_gene_medians.json \\
      --head-ckpt $R/outputs/cite_phase1_official/best.pt --out-dir $R/outputs/v3/E2 --device auto --threads 4
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
import ast  # noqa: E402
import contextlib  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bridge_anm.lib import e2_response as er  # noqa: E402
from bridge_anm.lib import v3_amend as va  # noqa: E402
from bridge_anm.lib import v3_key as vk  # noqa: E402

SCRIPT_VERSION = "e2_response_decomposition v1"
GENE_REFERENCE = "data/raw/external_hao2021_gse164378/ref_10x_grch38_3.0.0/features_grch38_3.0.0.tsv.gz"
TOKEN_PROBES = ROOT / "scripts" / "mode_a_token_probes.py"
C_GRID = (0.01, 0.1, 1.0, 10.0, 100.0)
N_COUNT_CHECK = 2000
COUNT_TOL = 1.0  # deviation ratio (e2_response.recover_counts): valid when <= 1
N_BOOT = 2000
EPS_SMALL, EPS_TARGET = 0.05, 0.8  # criterion (ii): thinning eps (1 - kept fraction)
MIN_LOST = 30  # A1.7
LOSS_TOL, SHARE_MARGIN, AUROC_MARGIN = 0.02, 0.15, 0.02
KIND_FLOOR, KIND_THIN, KIND_GENE = 0, 1, 2
_LOG = {"path": None}

SPEC = {
    "experiment": "E2",
    "addendum_version": 1,
    "addendum_to": "registration/registration_v3.json experiments.E2 as amended by A1 (A1.7 falsification text, A1.1 "
                   "bootstrap); the registered cells, perturbations, measurements, linearity test, decomposition, "
                   "comparator, margins and falsification are not changed here",
    "script": "scripts/e2_response_decomposition.py (" + SCRIPT_VERSION + ")",
    "fixed_before_site4": "fixed on training and validation cells only and committed before any site4 forward pass of "
                          "this script; smoke runs use val donor 18303 cells only",
    "open_choices_fixed": {
        "genes": "E2 genes = panel_coding_genes flattened, then nk_t_genes, first occurrence kept; Ensembl ids from the "
                 f"10x GRCh38 3.0.0 features file {GENE_REFERENCE} (the gene reference E5 uses; a 10x reference "
                 "table, not the external dataset), checked to give the ids of scripts/mode_a_token_probes.py GENES "
                 "for the 13 NK/T genes; genes absent from the processed gene list or the TEDDY vocabulary are "
                 "dropped and listed (computed.genes)",
        "z_probe": "sklearn LogisticRegression(solver lbfgs, max_iter 3000, multinomial, no feature scaling) on "
                   "L2-normalised z (z / (||z|| + 1e-6), the head's input) of training primary-key cells (5 classes, "
                   "unscored excluded); C = argmin of the val primary-key log-loss over [0.01, 0.1, 1, 10, 100] (ties: "
                   "the smaller C); the selected fit is saved to <out-dir>/z_probe.npz and the sha256 of its "
                   "coefficients is recorded in computed.z_probe",
        "baseline": "each cell's stored official embedding (data/processed/cite_official/z_rna.npy) -> head -> "
                    "evidence -> Q1 rule (the baseline E1 uses); every response is delta(eps) = x(perturbed) - "
                    "x(baseline)",
        "floor": "one unperturbed re-embedding per cell through this script's pipeline; floor_x = ||x(re-embedding) - "
                 "x(baseline)|| per cell and stage; the floor's own call changes are reported",
        "counts": "the stored RNA values are UMI counts times one per-cell factor; counts = rint(q), q = value / "
                  "smallest non-zero value of the row; a row is valid when every |q - rint(q)| <= 1e-3 + 4e-7 * "
                  "rint(q) (float32 rounding of the value and of the unit); invalid rows are flagged and reported "
                  f"(computed.count_check: {N_COUNT_CHECK} train/val cells drawn with numpy "
                  "default_rng(seeds.e2_thinning[0]))",
        "thinning": "x' ~ Binomial(counts, kept fraction) per gene with numpy default_rng([seed, global cell id, level "
                    "index]), level index in the registered order [0.95, 0.9, 0.8, 0.6, 0.2, 0.05], seeds 11 and 12 "
                    "(seeds.e2_thinning): independent draws per level and seed; thinned counts then go through the "
                    "official preprocessing (counts/total x 1e4, / gene medians, top 2,048, zeros dropped)",
        "gene_scaling": "for each E2 gene among the cell's baseline tokens, its official value (after /total x 1e4 and "
                        "/ median) is multiplied by (1 - eps), eps in [1.0, 0.5, 0.25], before torch.topk; eps = 1 "
                        "removes the token",
        "embedding": "the official pipeline of scripts/03_embed_rna.py: model.hidden_states with the bool key-padding "
                     "mask under torch.autocast fp16 on the device (MPS), gene-mean pooling; the variants of a chunk "
                     "of cells are length-sorted and run in batches of 32 padded to the batch's longest sequence",
        "empty_variant": "a variant left with no tokens gets no embedding: no call, z-probe class undefined (a lost case "
                         "is then 'representation'), excluded from the linearity rows (counted)",
        "measurements": "dz = ||l2(z') - l2(z)||; dv = ||v' - v|| over the 12 primary-panel evidences (clipped, as "
                        "registered; the unclipped norm is reported beside it); dS = S'[c0] - S[c0] for the baseline "
                        "argmax class c0 and dmargin = m(S') - m(S), m = top1 - top2; dANM = A'[c0] - A[c0], A = the "
                        "ANM engine's action score; dcall = Q1 call changed",
        "anm": "ANM finite_graph_scalar from ANM_ROOT on the baseline and every case: 4 actions, 3 support events per "
               "class (value = evidence), all at t = 0, registered field; readout threshold "
               "v3_key.anm_readout_threshold(Q1); ANM calls are compared with the rule's on every case (a mismatch is "
               "reported as a bridge failure) and dANM with G(3) * 3 * dS",
        "linearity": "stages z (512-d l2(z)), evidence (12 panel values), score (4 rule scores); thinning uses the "
                     "seed-averaged response mean_s delta_s(eps) per cell; gene scaling uses each (cell, gene) case; "
                     "per pair and stage: share linear, share below floor, median rho and cosine over the test_primary "
                     "e2_subset cells (thinning) or cases (gene scaling); also per gene",
        "units": "E2.dec and criterion (i) units = perturbation family x level: thinning per kept fraction (6 units; "
                 "cases = cell x 2 seeds) and gene scaling per eps pooled over genes (3 units; cases = (cell, gene), "
                 "each case weighted equally); secondary rows apply the same procedure to the thinning units plus the "
                 "finer gene x eps units and never decide the verdict",
        "accuracy_loss_for_i": "Q1 decision accuracy on the primary key (correct call, or no call on a key-OUT case) "
                               "over the unit's scored cases, baseline minus perturbed; the selective-accuracy version "
                               "is reported beside it",
        "lost_cells_for_i": "distinct test_primary e2_subset cells with at least one lost case in the unit (>= 30 in "
                            "both units, A1.7); share of a label = lost cases with the label / lost cases",
        "interval_for_i": "pairs qualify on the point estimates; the share difference's interval is the registered "
                          "two-stage bootstrap (B = 2000, seed 1, donors in sorted order; cells carry all their "
                          "cases); replicates in which a unit has no lost case are dropped from that interval and "
                          "counted; every (pair, label) test is listed; a Bonferroni-adjusted interval and the "
                          "per-donor clause of the common win rule are reported beside the verdict (secondary)",
        "criterion_ii": "cases = (test_primary e2_subset cell, seed) at kept fraction 0.2 (eps 0.8); y = the Q1 call "
                        "differs from the baseline call; comparator score = -m (baseline margin top1 - top2); "
                        "predictor score = -m_hat, m_hat = m + (0.8 / 0.05) * (mean_s msigned_s - m), msigned_s = "
                        "S'[c0] - max_{k != c0} S'[k] at kept fraction 0.95, seed s (the comparator's baseline "
                        "information plus the small-eps response); dAUROC = AUROC(-m_hat) - AUROC(-m), weighted "
                        "Mann-Whitney (ties 1/2); (ii) holds if dAUROC >= 0.02 and the bootstrap lower bound > 0. "
                        "Secondary rows: the bar-aware distance d = min(signed margin, distance of S[c0] to the bar "
                        "on the side of the baseline call), extrapolated the same way, against baseline d and against m",
        "verdict": "'adds information beyond the curve' if (i) or (ii) holds (A1's E2 text); otherwise 'adds nothing "
                   "beyond the curve (falsified)'; whether the per-donor clause of the common win rule also holds is "
                   "reported with each criterion",
        "secondary_split": "donor 15078 (250 cells) is reported with the same procedure and never changes the verdict",
        "smoke": "smoke runs use val donor 18303 cells only (--cell-pool val); no site4 RNA, embedding, protein or "
                 "label is read before this addendum is committed",
    },
    "outputs": "outputs/v3/E2/: E2_results.json, REPORT.md, progress.log, e2_progress.json, z_probe.npz, chunks/, "
               "e2_cases.npz; every result JSON carries registration_sha256, amendment_sha256 and addendum_sha256",
}


def spec_sha() -> str:
    return hashlib.sha256((json.dumps(SPEC, indent=1, sort_keys=True) + "\n").encode()).hexdigest()


def say(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG["path"] is not None:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


def write_json(path: Path, obj) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default))
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def clean(x):
    """JSON-safe copy: NaN/inf -> None, numpy -> python, floats rounded to 6 decimals."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return round(float(x), 6) if np.isfinite(x) else None
    return x


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


# ============================================================================ args

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "embed", "report", "all"), default="all")
    p.add_argument("--cell-pool", choices=("e2_subset", "val"), default="e2_subset",
                   help="val = val donor 18303 cells (smoke only; no site4 data)")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=ROOT / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--gene-reference", type=Path, default=ROOT / GENE_REFERENCE)
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--z-probe", type=Path, default=None,
                   help="z-probe file written by the register stage (default <out-dir>/z_probe.npz, else the file "
                        "named in the addendum)")
    p.add_argument("--anm-root", type=Path, default=Path(os.environ.get("ANM_ROOT", str(ROOT.parent / "ANM"))))
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--chunk-cells", type=int, default=10, help="cells per resumable chunk file")
    p.add_argument("--limit", type=int, default=0, help="smoke only: the first N cells of the pool")
    p.add_argument("--n-boot", type=int, default=N_BOOT, help="registered 2000; smaller only for smoke runs")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    if (a.limit or a.n_boot != N_BOOT) and not a.smoke:
        p.error("--limit / --n-boot are for smoke runs only (give --smoke NOTE)")
    if a.cell_pool == "val" and not a.smoke and a.stage != "register":
        p.error("--cell-pool val is for smoke runs only")
    if a.cell_pool == "e2_subset" and a.smoke and a.stage != "register":
        p.error("smoke runs use --cell-pool val (no site4 forward inside smoke tests)")
    return a


# ============================================================================ registration

def _git_committed(path: Path) -> bool | None:
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def addendum_path(a) -> Path:
    return a.registration_dir / "addenda" / "E2.json"


def check_registration(a) -> dict:
    """Hash checks; refuses a site4 stage unless every registration file is committed and matches."""
    regp = a.registration_dir / "registration_v3.json"
    amp = a.registration_dir / "amendment_A1.json"
    reg = va.load_registration_amended(regp, amp)  # verifies both sha256 files
    info = {"registration_sha256": vk.sha256_file(regp), "amendment_sha256": va.amendment_sha256(amp),
            "registration_committed": _git_committed(regp), "amendment_committed": _git_committed(amp),
            "spec_sha256": spec_sha()}
    f = addendum_path(a)
    hf = a.registration_dir / "addenda" / "HASHES.txt"
    info["addendum_file"] = str(f)
    info["addendum_exists"] = f.exists()
    if f.exists():
        sha = vk.sha256_file(f)
        add = json.loads(f.read_text())
        info["addendum_sha256"] = sha
        info["addendum_hash_recorded"] = hf.exists() and f"{sha}  E2.json" in hf.read_text().splitlines()
        info["addendum_committed"] = _git_committed(f)
        info["addendum_spec_matches_script"] = add.get("spec_sha256") == spec_sha() and add.get("spec") == SPEC
        info["addendum_registration_matches"] = (add.get("registration_sha256") == info["registration_sha256"]
                                                 and add.get("amendment_sha256") == info["amendment_sha256"])
    site4 = a.cell_pool == "e2_subset" and a.stage in ("embed", "report", "all")
    if site4:
        need = {"registration_v3.json committed": info["registration_committed"] is True,
                "amendment_A1.json committed": info["amendment_committed"] is True,
                "addenda/E2.json exists": info["addendum_exists"],
                "addenda/E2.json sha256 in addenda/HASHES.txt": info.get("addendum_hash_recorded") is True,
                "addenda/E2.json committed": info.get("addendum_committed") is True,
                "addenda/E2.json spec equals this script's SPEC": info.get("addendum_spec_matches_script") is True,
                "addenda/E2.json names this registration and A1": info.get("addendum_registration_matches") is True}
        bad = [k for k, ok in need.items() if not ok]
        if bad:
            raise SystemExit("site4 run refused (register and commit the E2 addendum first): " + "; ".join(bad))
    return {"info": info, "reg": reg}


def load_addendum(a) -> dict:
    f = addendum_path(a)
    if not f.exists():
        raise SystemExit(f"missing {f}: run --stage register first")
    return json.loads(f.read_text())


# ============================================================================ data

def load_pack(processed: Path, need_rna: bool) -> dict:
    try:
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=False)
        _ = npz["adt_names"]
    except ValueError:
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=True)
    out = {k: npz[k].astype(str) for k in ("split", "sites", "donors", "cell_types")}
    out["adt_names"] = [str(x) for x in npz["adt_names"]]
    out["rna_names"] = npz["rna_names"].astype(str)
    out["token_ids"] = npz["token_ids"].astype(np.int64)
    out["_npz"] = npz
    if need_rna:
        from scipy import sparse
        out["rna"] = sparse.csr_matrix((npz["rna_data"], npz["rna_indices"], npz["rna_indptr"]),
                                       shape=tuple(npz["rna_shape"]))
    return out


def adt_rows(pack: dict, rows: np.ndarray) -> np.ndarray:
    """Measured ADT of the given rows only (verifier side: keys)."""
    return np.asarray(pack["_npz"]["adt"], dtype=np.float32)[np.asarray(rows, dtype=np.int64)]


def token_probe_ids() -> dict[str, str]:
    """GENES (symbol -> Ensembl) of scripts/mode_a_token_probes.py, read from its source with ast."""
    tree = ast.parse(TOKEN_PROBES.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "GENES" for t in node.targets):
            return {str(k): str(v) for k, v in ast.literal_eval(node.value).items()}
    raise SystemExit(f"GENES not found in {TOKEN_PROBES}")


def map_genes(reg: dict, gene_ref: Path, rna_names: np.ndarray, vocab: dict, token_ids: np.ndarray) -> dict:
    G = reg["experiments"]["E2"]["perturbations"]["genes"]
    syms = er.flatten_genes(G["panel_coding_genes"], G["nk_t_genes"])
    sym2ens: dict[str, str] = {}
    with gzip.open(gene_ref, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in syms and parts[1] not in sym2ens:
                sym2ens[parts[1]] = parts[0]
    probe_ids = token_probe_ids()
    bad = {s: (sym2ens.get(s), e) for s, e in probe_ids.items() if sym2ens.get(s) != e}
    if bad:
        raise SystemExit(f"gene reference disagrees with mode_a_token_probes ids: {bad}")
    col = {str(g): i for i, g in enumerate(map(str, rna_names))}
    kept, dropped = [], []
    for s in syms:
        e = sym2ens.get(s)
        if e is not None and e in col and e in vocab:
            if int(token_ids[col[e]]) != int(vocab[e]):
                raise SystemExit(f"{s}: processed column {col[e]} has token id {token_ids[col[e]]}, vocab says {vocab[e]}")
            kept.append({"symbol": s, "ensembl": e, "column": int(col[e]), "vocab_id": int(vocab[e])})
        else:
            dropped.append({"symbol": s, "ensembl": e, "in_processed": bool(e in col) if e else False,
                            "in_vocab": bool(e in vocab) if e else False})
    nkt = set(map(str, G["nk_t_genes"]))
    panel = {s for v in G["panel_coding_genes"].values() for s in (v if isinstance(v, list) else [v])}
    return {"reference": GENE_REFERENCE, "reference_sha256": vk.sha256_file(gene_ref), "symbols": syms,
            "kept": kept, "dropped": dropped, "panel_genes": sorted(panel), "nk_t_genes": sorted(nkt)}


def cell_pool(a, reg: dict, pack: dict) -> tuple[np.ndarray, np.ndarray]:
    """(global cell ids, role) in pool order; role 'primary' / 'secondary' (e2_subset) or 'val_smoke'."""
    if a.cell_pool == "e2_subset":
        e = reg["e2_subset"]
        prim = np.asarray(e["primary_indices"], dtype=np.int64)
        sec = np.asarray(e["secondary_indices"], dtype=np.int64)
        if vk.index_hash(prim) != e["sha256_primary"] or vk.index_hash(sec) != e["sha256_secondary"]:
            raise SystemExit("e2_subset indices differ from their registered hashes")
        cells = np.concatenate([prim, sec])
        role = np.array(["primary"] * prim.size + ["secondary"] * sec.size)
    else:
        va_idx = vk.split_indices(pack["split"], pack["sites"], pack["donors"], reg, "val")
        rng = np.random.default_rng(int(reg["seeds"]["e2_subset"]))
        n = a.limit if a.limit else 40
        cells = np.sort(rng.choice(va_idx, size=min(n, va_idx.size), replace=False)).astype(np.int64)
        role = np.array(["val_smoke"] * cells.size)
    if a.limit:
        cells, role = cells[: a.limit], role[: a.limit]
    return cells, role


# ============================================================================ register (train/val only)

def l2_rows(z: np.ndarray) -> np.ndarray:
    return er.l2n(np.asarray(z, dtype=np.float64))


def coef_sha(coef: np.ndarray, intercept: np.ndarray, classes: np.ndarray) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(coef, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(intercept, dtype=np.float64).tobytes())
    h.update(",".join(map(str, classes)).encode())
    return h.hexdigest()


def register_core(z: np.ndarray, adt: np.ndarray, cell_types: np.ndarray, rna_row, split: np.ndarray,
                  sites: np.ndarray, donors: np.ndarray, adt_names: list, reg: dict, *, return_model: bool = False):
    """The computed part of the addendum, from training and val rows only.

    z, adt, cell_types are full-length arrays; only the train and val rows are ever indexed. ``rna_row(i)``
    returns a dense RNA row and is only called for train/val rows (the count check)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss

    tr = vk.split_indices(split, sites, donors, reg, "train")
    vl = vk.split_indices(split, sites, donors, reg, "val")
    ktr = vk.build_keys(cell_types[tr], adt[tr], adt_names, reg)["primary"]
    kvl = vk.build_keys(cell_types[vl], adt[vl], adt_names, reg)["primary"]
    mtr, mvl = ktr != vk.UNSCORED, kvl != vk.UNSCORED
    Xtr, ytr = l2_rows(z[tr[mtr]]), ktr[mtr]
    Xvl, yvl = l2_rows(z[vl[mvl]]), kvl[mvl]
    grid, models = [], {}
    for C in C_GRID:
        m = LogisticRegression(C=C, solver="lbfgs", max_iter=3000)
        m.fit(Xtr, ytr)
        p = m.predict_proba(Xvl)
        ll = float(log_loss(yvl, p, labels=m.classes_))
        grid.append({"C": C, "val_log_loss": round(ll, 6), "val_accuracy": round(float(np.mean(m.predict(Xvl) == yvl)), 6),
                     "n_iter": int(np.max(m.n_iter_))})
        models[C] = m
    best = min(grid, key=lambda r: (r["val_log_loss"], r["C"]))
    m = models[best["C"]]
    # count recovery check on seeded train/val rows
    rng = np.random.default_rng(int(reg["seeds"]["e2_thinning"][0]))
    pool = np.concatenate([tr, vl])
    sample = np.sort(rng.choice(pool, size=min(N_COUNT_CHECK, pool.size), replace=False))
    devs = np.array([er.recover_counts(rna_row(int(i)))[2] for i in sample])
    core = {
        "z_probe": {"C_grid": grid, "C": best["C"], "classes": [str(c) for c in m.classes_],
                    "n_train": int(Xtr.shape[0]), "n_val": int(Xvl.shape[0]),
                    "coef_sha256": coef_sha(m.coef_, m.intercept_, m.classes_)},
        "count_check": {"n_cells": int(sample.size), "share_passing": round(float(np.mean(devs <= COUNT_TOL)), 6),
                        "max_deviation_ratio": round(float(np.max(devs)), 6),
                        "rule": "valid when max |q - rint(q)| / (1e-3 + 4e-7 rint(q)) <= 1"},
    }
    return (core, m) if return_model else core


def stage_register(a) -> None:
    reg = va.load_registration_amended(a.registration_dir / "registration_v3.json",
                                       a.registration_dir / "amendment_A1.json")
    say("register: train/val only")
    pack = load_pack(a.processed, need_rna=True)
    from teddy_mm.teddy_encoder import load_vocab
    vocab = load_vocab(a.ckpt)
    genes = map_genes(reg, a.gene_reference, pack["rna_names"], vocab, pack["token_ids"])
    say(f"genes kept {[g['symbol'] for g in genes['kept']]}, dropped {[g['symbol'] for g in genes['dropped']]}")
    z = np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")
    adt = np.asarray(pack["_npz"]["adt"], dtype=np.float32)
    ct = pack["cell_types"]
    rna = pack["rna"]
    args = (pack["split"], pack["sites"], pack["donors"], pack["adt_names"], reg)
    t0 = time.time()
    core, model = register_core(z, adt, ct, lambda i: rna[i].toarray().ravel(), *args, return_model=True)
    say(f"z-probe C = {core['z_probe']['C']} ({time.time() - t0:.0f}s); count check {core['count_check']}")
    # leakage self-check: site4 rows poisoned -> identical core; val rows poisoned -> changed core
    test = np.where(pack["split"] == "test")[0]
    valr = np.where(pack["split"] == "val")[0]
    rng = np.random.default_rng(int(reg["seeds"]["global"]))

    def poisoned(rows):
        zp = np.array(z, dtype=np.float32, copy=True)
        zp[rows] = rng.normal(size=(rows.size, zp.shape[1])).astype(np.float32)
        ap = adt.copy()
        ap[rows] = rng.normal(size=(rows.size, ap.shape[1])).astype(np.float32)
        cp = ct.copy()
        cp[rows] = rng.permutation(np.asarray(sorted(set(ct.tolist()))))[rng.integers(0, len(set(ct.tolist())), rows.size)]
        rs = set(rows.tolist())
        return zp, ap, cp, (lambda i: rng.poisson(3.0, size=rna.shape[1]).astype(np.float32) * 0.7 if i in rs
                            else rna[i].toarray().ravel())

    zp, ap, cp, rr = poisoned(test)
    core_t = register_core(zp, ap, cp, rr, *args)
    zp, ap, cp, rr = poisoned(valr)
    core_v = register_core(zp, ap, cp, rr, *args)
    leak = {"site4_poisoned_core_identical": core_t == core, "val_poisoned_core_changed": core_v != core}
    say(f"leakage self-check {leak}")
    if not (leak["site4_poisoned_core_identical"] and leak["val_poisoned_core_changed"]):
        raise SystemExit(f"leakage self-check failed: {leak}")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    zp_file = a.out_dir / "z_probe.npz"
    atomic_savez(zp_file, coef=model.coef_.astype(np.float64), intercept=model.intercept_.astype(np.float64),
                 classes=np.asarray(model.classes_).astype(str), C=np.array(core["z_probe"]["C"]))
    core["z_probe"]["file"] = str(zp_file)
    add = {"experiment": "E2", "addendum_version": SPEC["addendum_version"],
           "registration_sha256": vk.sha256_file(a.registration_dir / "registration_v3.json"),
           "amendment_sha256": va.amendment_sha256(a.registration_dir / "amendment_A1.json"),
           "spec_sha256": spec_sha(), "spec": SPEC,
           "computed": {**core, "genes": genes, "leakage_self_check": leak},
           "provenance": {"script_sha256": vk.sha256_file(Path(__file__)),
                          "e2_response_sha256": vk.sha256_file(ROOT / "bridge_anm/lib/e2_response.py"),
                          "z_sha256": vk.sha256_file(a.embed_dir / "z_rna.npy"),
                          "cite_arrays_sha256": vk.sha256_file(a.processed / "cite_arrays.npz"),
                          "rows_read": "train and val rows of z, adt, cell types and RNA; split/site/donor metadata of "
                                       "all rows; site4 rows only inside the poisoned rebuild",
                          "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    d = a.registration_dir / "addenda"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "E2.json"
    f.write_text(json.dumps(add, indent=1, sort_keys=True) + "\n")
    sha = vk.sha256_file(f)
    hf = d / "HASHES.txt"
    lines = [ln for ln in (hf.read_text().splitlines() if hf.exists() else []) if not ln.rstrip().endswith("  E2.json")]
    hf.write_text("\n".join(lines + [f"{sha}  E2.json"]) + "\n")
    say(f"wrote {f} (sha256 {sha}) and its line in {hf}; commit both before any site4 run")


# ============================================================================ embed

def resolve_device(name: str):
    from teddy_mm.device import resolve_device as rd
    return rd(name)


def variants_for_cell(cell: int, row: np.ndarray, factors: np.ndarray, token_ids: np.ndarray, pad_id: int,
                      genes: list[dict], reg: dict, seq_len: int, norm: float) -> dict:
    """All registered variants of one cell: token id arrays plus their descriptors."""
    from teddy_mm.teddy_encoder import official_values, rank_encode_official

    P = reg["experiments"]["E2"]["perturbations"]
    kept = [float(f) for f in P["rna_thinning"]["fractions_kept"]]
    seeds = [int(s) for s in reg["seeds"]["e2_thinning"]]
    eps_g = [float(e) for e in P["gene_scaling"]["epsilon"]]
    base_vals = official_values(row[None, :], factors, norm)
    counts, unit, dev = er.recover_counts(row)
    rows, desc = [base_vals[0]], [(KIND_FLOOR, 0.0, -1, -1)]
    for li, f in enumerate(kept):
        for s in seeds:
            c2 = er.thin_counts(counts, f, er.thinning_rng(s, cell, li))
            rows.append(official_values(c2[None, :].astype(np.float32), factors, norm)[0])
            desc.append((KIND_THIN, f, s, -1))
    btok, battn = rank_encode_official(base_vals, token_ids, max_len=seq_len, pad_id=pad_id, pad_to=seq_len)
    base_ids = btok[0, : int(battn[0].sum())]
    present = set(base_ids.tolist())
    for gi, g in enumerate(genes):
        if g["vocab_id"] not in present:
            continue
        for e in eps_g:
            rows.append(er.scale_gene(base_vals[0], g["column"], 1.0 - e))
            desc.append((KIND_GENE, e, -1, gi))
    V = np.stack(rows).astype(np.float32)
    tok, attn = rank_encode_official(V, token_ids, max_len=seq_len, pad_id=pad_id, pad_to=seq_len)
    ids = [tok[i, : int(attn[i].sum())].copy() for i in range(V.shape[0])]
    ident = np.array([x.size == base_ids.size and np.array_equal(x, base_ids) for x in ids])
    return {"ids": ids, "desc": desc, "identical": ident, "count_unit": unit, "count_dev": dev,
            "ntok_base": int(base_ids.size)}


def embed_sequences(model, seqs: list[np.ndarray], pad_id: int, device, bs: int) -> np.ndarray:
    """Gene-mean of the last layer for each token sequence (official path; fp16 autocast on the device)."""
    import torch

    from teddy_mm.teddy_encoder import pool_hidden

    out = np.full((len(seqs), model.d_model), np.nan, dtype=np.float32)
    lens = np.array([s.size for s in seqs])
    order = [i for i in np.argsort(lens, kind="stable") if lens[i] > 0]
    amp = torch.float16
    with torch.no_grad():
        for b in range(0, len(order), bs):
            rows = order[b: b + bs]
            width = int(max(lens[r] for r in rows))
            ids = np.full((len(rows), width), pad_id, dtype=np.int64)
            msk = np.zeros((len(rows), width), dtype=np.int64)
            for j, r in enumerate(rows):
                ids[j, : lens[r]] = seqs[r]
                msk[j, : lens[r]] = 1
            it, mt = torch.from_numpy(ids).to(device), torch.from_numpy(msk).to(device)
            ctx = torch.autocast(device_type=device.type, dtype=amp) if device.type in ("mps", "cuda", "cpu") \
                else contextlib.nullcontext()
            with ctx:
                h = model.hidden_states(it, mt)
            out[rows] = pool_hidden(h.float(), mt, "gene-mean").float().cpu().numpy()
    return out


def stage_embed(a, reg: dict, reg_info: dict) -> None:
    import torch

    from teddy_mm.teddy_encoder import load_gene_medians, load_pad_id, load_teddy, load_vocab, median_factors

    torch.set_num_threads(max(1, a.threads))
    add = load_addendum(a)
    pack = load_pack(a.processed, need_rna=True)
    cells, role = cell_pool(a, reg, pack)
    man = json.loads((a.embed_dir / "z_rna_manifest.json").read_text())
    seq_len, norm = int(man.get("seq_len", 2048)), float(man.get("normalize_total", 1e4))
    if man.get("autocast") != "fp16" or man.get("pooling") != "gene-mean":
        raise SystemExit(f"embedding manifest is not the official fp16 gene-mean run: {man.get('autocast')}, {man.get('pooling')}")
    if vk.sha256_file(a.medians) != man.get("medians_sha256"):
        raise SystemExit(f"{a.medians} differs from the gene medians of the official embedding run")
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    factors = median_factors(pack["rna_names"], load_gene_medians(a.medians))
    genes = add["computed"]["genes"]["kept"]
    ntok_stored = np.load(a.embed_dir / "z_rna_ntokens.npy")
    device = resolve_device(a.device)
    chunk_dir = a.out_dir / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    n_chunks = er.ceil_div(cells.size, a.chunk_cells)
    tag = {"spec_sha256": spec_sha(), "addendum_sha256": reg_info.get("addendum_sha256"), "pool": a.cell_pool}
    model = None
    t_start, done_now, n_seq_now = time.time(), 0, 0
    say(f"embed: {cells.size} cells in {n_chunks} chunks of {a.chunk_cells}; device {device}; seq_len {seq_len}")
    for k in range(n_chunks):
        path = chunk_dir / f"chunk_{k:04d}.npz"
        sel = cells[k * a.chunk_cells: (k + 1) * a.chunk_cells]
        if path.exists():
            with np.load(path) as old:
                if str(old["tag"]) != json.dumps(tag, sort_keys=True) or not np.array_equal(old["cells"], sel):
                    raise SystemExit(f"{path} was written under another addendum / pool; use a new --out-dir")
            continue
        if model is None:
            model = load_teddy(a.ckpt, device)
        t0 = time.time()
        seqs, vcell, vkind, vlevel, vseed, vgene, vident = [], [], [], [], [], [], []
        cunit, cdev, cntok = [], [], []
        for c in sel:
            row = pack["rna"][int(c)].toarray().ravel().astype(np.float32)
            V = variants_for_cell(int(c), row, factors, pack["token_ids"], pad_id, genes, reg, seq_len, norm)
            for ids, (kind, lev, sd, gi), ident in zip(V["ids"], V["desc"], V["identical"]):
                seqs.append(ids)
                vcell.append(int(c))
                vkind.append(kind)
                vlevel.append(lev)
                vseed.append(sd)
                vgene.append(gi)
                vident.append(bool(ident))
            cunit.append(V["count_unit"])
            cdev.append(V["count_dev"])
            cntok.append(V["ntok_base"])
        Z = embed_sequences(model, seqs, pad_id, device, a.batch_size)
        atomic_savez(path, tag=np.array(json.dumps(tag, sort_keys=True)), cells=sel,
                     var_cell=np.array(vcell, dtype=np.int64), var_kind=np.array(vkind, dtype=np.int8),
                     var_level=np.array(vlevel, dtype=np.float64), var_seed=np.array(vseed, dtype=np.int64),
                     var_gene=np.array(vgene, dtype=np.int64), var_ntok=np.array([s.size for s in seqs], dtype=np.int32),
                     var_identical=np.array(vident), z=Z, cell_count_unit=np.array(cunit),
                     cell_count_dev=np.array(cdev), cell_ntok_base=np.array(cntok, dtype=np.int32),
                     cell_ntok_stored=ntok_stored[sel].astype(np.int32), runtime_sec=np.array(time.time() - t0))
        done_now += 1
        n_seq_now += len(seqs)
        el = time.time() - t_start
        left = sum(1 for j in range(k + 1, n_chunks) if not (chunk_dir / f"chunk_{j:04d}.npz").exists())
        eta = el / done_now * left
        say(f"chunk {k + 1}/{n_chunks}: {sel.size} cells, {len(seqs)} sequences in {time.time() - t0:.0f}s "
            f"({1000 * (time.time() - t0) / max(len(seqs), 1):.0f} ms/seq); ETA {eta / 60:.1f} min")
        write_json(a.out_dir / "e2_progress.json", {"chunks_total": n_chunks, "chunk_last": k, "elapsed_sec": el,
                                                    "eta_sec": eta, "sequences_this_run": n_seq_now,
                                                    "sec_per_sequence": el / max(n_seq_now, 1)})
    say("embed: all chunks present")


# ============================================================================ report

def load_chunks(a, cells: np.ndarray) -> dict:
    files = sorted((a.out_dir / "chunks").glob("chunk_[0-9][0-9][0-9][0-9].npz"))
    n_chunks = er.ceil_div(cells.size, a.chunk_cells)
    if len(files) < n_chunks:
        raise SystemExit(f"{len(files)} of {n_chunks} chunks present; run --stage embed first")
    keys = ("var_cell", "var_kind", "var_level", "var_seed", "var_gene", "var_ntok", "var_identical", "z")
    ckeys = ("cells", "cell_count_unit", "cell_count_dev", "cell_ntok_base", "cell_ntok_stored")
    out = {k: [] for k in keys + ckeys}
    rt = 0.0
    for f in files[:n_chunks]:
        with np.load(f) as d:
            for k in keys + ckeys:
                out[k].append(d[k])
            rt += float(d["runtime_sec"])
    res = {k: np.concatenate(v) for k, v in out.items()}
    if not np.array_equal(res["cells"], cells):
        raise SystemExit("chunk cells differ from the pool order")
    res["embed_runtime_sec"] = rt
    return res


class Pipeline:
    """z -> head -> evidence -> Q1 rule scores / calls, z-probe class, ANM engine scores / calls."""

    def __init__(self, a, reg: dict, probe: dict, pack: dict):
        self.a, self.reg = a, reg
        self.names = pack["adt_names"]
        self.bar = float(reg["questions"]["Q1"]["bar"])
        pan = vk.question_panel(reg, "Q1")
        j = {n: i for i, n in enumerate(self.names)}
        self.panel_cols = [j[p] for k in vk.LINEAGES for p in pan[k]]  # 12, class-major
        self.coef, self.icpt, self.classes = probe["coef"], probe["intercept"], probe["classes"]
        self.anm_thr = vk.anm_readout_threshold(reg, "Q1")
        self.G3 = vk.field_gain(3, reg["anm"]["field_representation"])
        self._ffr = None

    def head(self, z: np.ndarray) -> dict:
        ok = np.all(np.isfinite(z), axis=1)
        n = z.shape[0]
        pred = np.full((n, len(self.names)), np.nan)
        if ok.any():
            pred[ok] = vk.head_predict(z[ok], self.a.head_ckpt, device="cpu")
        v_raw = vk.evidence(np.where(ok[:, None], pred, 0.0), self.reg, clip=False)
        v = np.clip(v_raw, 0.0, 1.0)
        S = vk.class_scores(v, self.names, self.reg, "Q1")
        calls = vk.rule_calls(S, self.bar).astype(object)
        calls[~ok] = vk.NO_CALL
        argmax = np.asarray([vk.LINEAGES[i] for i in S.argmax(axis=1)], dtype=object)
        argmax[~ok] = ""
        zl = er.l2n(np.where(ok[:, None], z, 0.0))
        logits = zl @ self.coef.T + self.icpt
        probe = np.asarray(self.classes)[logits.argmax(axis=1)].astype(object)
        probe[~ok] = ""
        return {"ok": ok, "zl": zl, "v12": v[:, self.panel_cols], "v12_raw": v_raw[:, self.panel_cols], "S": S,
                "call": calls.astype(str), "argmax": argmax.astype(str), "probe": probe.astype(str)}

    def anm(self, v12: np.ndarray, ok: np.ndarray) -> dict:
        if self._ffr is None:
            sys.path.insert(0, str(self.a.anm_root))
            from active_neural_matter.field import finite_field_runner as ffr
            base_schema = json.loads((ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
            self._ffr = ffr
            self._schema = er.anm_schema(base_schema, self.reg["anm"]["field_representation"], self.anm_thr)
        A, calls = er.anm_engine(v12, ok, self._ffr, self._schema, vk.LINEAGES)
        return {"A": A, "call": calls}


def boot_ci(rep: np.ndarray) -> list:
    lo, hi, nb = er.percentile_ci(rep)
    return [lo, hi] if nb == 0 else [lo, hi, f"{nb} undefined replicates dropped"]


def unit_block(T, W, per_donor_T: dict) -> dict:
    pt = er.unit_stats(T)
    out = {"point": pt}
    if W is not None:
        bs = er.unit_stats(T, W)
        out["ci95"] = {k: boot_ci(v) for k, v in bs.items() if isinstance(v, np.ndarray)}
        out["_boot"] = bs
    out["per_donor"] = {d: er.unit_stats(t) for d, t in per_donor_T.items()}
    return out


def analyse(H: dict, base: dict, case: dict, key: np.ndarray, donor: np.ndarray, gene_syms: list[str],
            reg: dict, n_boot: int, seed: int, bar: float) -> dict:
    """Every E2 endpoint on one population of cells (positions 0..n-1).

    H: perturbed pipeline outputs per case; base: baseline outputs per cell; case: descriptors
    (cell position, kind, level, seed, gene); key: primary key per cell; donor per cell."""
    n = key.size
    cc = case["cell"]
    kind, lev, sd, gi = case["kind"], case["level"], case["seed"], case["gene"]
    P = reg["experiments"]["E2"]["perturbations"]
    kept_levels = [float(f) for f in P["rna_thinning"]["fractions_kept"]]
    eps_g = [float(e) for e in P["gene_scaling"]["epsilon"]]
    donors = sorted(set(donor.tolist()))
    W = er.two_stage_weights(donor, n_boot, seed) if n_boot > 0 else None
    k_case = key[cc]
    cb_case = base["call"][cc]
    lost = er.lost_cases(cb_case, H["call"], k_case)
    label = er.label_lost(H["probe"], H["argmax"], k_case)

    def tables(mask, cells_mask=None):
        m = mask if cells_mask is None else (mask & cells_mask[cc])
        return er.unit_cell_tables(cc[m], n, key=k_case[m], call_base=cb_case[m], call_pert=H["call"][m],
                                   label=label[m])

    def block(mask):
        T = tables(mask)
        pd = {d: tables(mask, donor == d) for d in donors}
        return unit_block(T, W, pd)

    units, units_fine = {}, {}
    for f in kept_levels:
        units[f"thin_keep{f:g}"] = block((kind == KIND_THIN) & np.isclose(lev, f))
    for e in eps_g:
        units[f"gene_eps{e:g}"] = block((kind == KIND_GENE) & np.isclose(lev, e))
        for g, sym in enumerate(gene_syms):
            m = (kind == KIND_GENE) & np.isclose(lev, e) & (gi == g)
            if m.any():
                units_fine[f"gene_{sym}_eps{e:g}"] = block(m)
    floor_blk = block(kind == KIND_FLOOR)

    # ---- measurements per unit (medians over cases)
    c0 = base["S"].argmax(axis=1)
    dS = H["S"][np.arange(cc.size), c0[cc]] - base["S"][cc, c0[cc]]
    dmargin = er.top_margin(H["S"]) - er.top_margin(base["S"])[cc]
    dA = H["A"][np.arange(cc.size), c0[cc]] - base["A"][cc, c0[cc]]
    dz = np.linalg.norm(H["zl"] - base["zl"][cc], axis=1)
    dv = np.linalg.norm(H["v12"] - base["v12"][cc], axis=1)
    dv_raw = np.linalg.norm(H["v12_raw"] - base["v12_raw"][cc], axis=1)
    okc = H["ok"]
    dz[~okc] = dv[~okc] = dv_raw[~okc] = dS[~okc] = dmargin[~okc] = dA[~okc] = np.nan

    def med(x, m):
        x = x[m & np.isfinite(x)]
        return {"median": float(np.median(x)) if x.size else None, "mean": float(np.mean(x)) if x.size else None,
                "n": int(x.size)}

    def meas(m):
        return {"dz": med(dz, m), "dv_clipped": med(dv, m), "dv_unclipped": med(dv_raw, m), "dS_called": med(dS, m),
                "dmargin": med(dmargin, m), "dANM_called": med(dA, m),
                "dcall_rate": float(np.mean(H["call"][m] != cb_case[m])) if m.any() else None,
                "n_cases": int(m.sum()), "n_empty": int(np.sum(m & ~okc)),
                "n_tokens_identical_to_baseline": int(np.sum(m & case["identical"]))}

    measurements = {"floor": meas(kind == KIND_FLOOR)}
    for f in kept_levels:
        measurements[f"thin_keep{f:g}"] = meas((kind == KIND_THIN) & np.isclose(lev, f))
    for e in eps_g:
        measurements[f"gene_eps{e:g}"] = meas((kind == KIND_GENE) & np.isclose(lev, e))
    anm_identity = float(np.nanmax(np.abs(dA - 3.0 * vk.field_gain(3, reg["anm"]["field_representation"]) * dS))) \
        if np.isfinite(dA).any() else None

    # ---- linearity (C5)
    floor_idx = np.full(n, -1)
    floor_idx[cc[kind == KIND_FLOOR]] = np.where(kind == KIND_FLOOR)[0]
    stages = {"z": ("zl", 512), "evidence": ("v12", 12), "score": ("S", 4)}

    def stage_floor(st):
        fx = np.full(n, np.nan)
        has = floor_idx >= 0
        fx[has] = np.linalg.norm(H[st][floor_idx[has]] - base[st][has], axis=1)
        return fx

    floors = {s: stage_floor(f) for s, (f, _) in stages.items()}

    def thin_mean_delta(f, st):
        """Seed-averaged response per cell at kept fraction f; NaN unless every seed is present and non-empty."""
        seeds = sorted(set(sd[(kind == KIND_THIN)].tolist()))
        acc = np.zeros((n, base[st].shape[1]))
        cnt = np.zeros(n)
        for s in seeds:
            m = (kind == KIND_THIN) & np.isclose(lev, f) & (sd == s) & okc
            acc[cc[m]] += H[st][m] - base[st][cc[m]]
            cnt[cc[m]] += 1
        out = np.full_like(acc, np.nan)
        full = cnt == len(seeds)
        out[full] = acc[full] / len(seeds)
        return out

    def lin_summary(Lr, floor_ok_rows):
        m = floor_ok_rows
        res = {"n": int(m.sum()), "share_linear": float(np.mean(Lr["linear"][m])) if m.any() else None,
               "share_below_floor": float(np.mean(Lr["below_floor"][m])) if m.any() else None}
        for k in ("rho", "cos"):
            x = Lr[k][m]
            x = x[np.isfinite(x)]
            res[f"median_{k}"] = float(np.median(x)) if x.size else None
        return res

    linear = {"thinning": {}, "gene_scaling": {}, "gene_scaling_per_gene": {}}
    thin_pairs = [(0.1, 0.05), (0.2, 0.1), (0.4, 0.2), (0.8, 0.4)]
    for e_full, e_half in thin_pairs:
        row = {}
        for s, (f, _) in stages.items():
            a1 = thin_mean_delta(round(1 - e_full, 6), f)
            a2 = thin_mean_delta(round(1 - e_half, 6), f)
            Lr = er.linearity(a1, a2, floors[s])
            row[s] = lin_summary(Lr, np.all(np.isfinite(a1), axis=1) & np.all(np.isfinite(a2), axis=1))
            row[s]["n_excluded_empty_or_missing"] = int(n - row[s]["n"])
        linear["thinning"][f"eps{e_full:g}_vs_eps{e_half:g}"] = row
    gene_pairs = [(1.0, 0.5), (0.5, 0.25)]
    # (cell, gene) -> case index per eps
    gidx = {}
    for e in eps_g:
        m = np.where((kind == KIND_GENE) & np.isclose(lev, e))[0]
        gidx[e] = {(int(cc[i]), int(gi[i])): i for i in m}
    for e_full, e_half in gene_pairs:
        common = sorted(set(gidx[e_full]) & set(gidx[e_half]))
        i1 = np.array([gidx[e_full][p] for p in common], dtype=np.int64)
        i2 = np.array([gidx[e_half][p] for p in common], dtype=np.int64)
        cells_p = np.array([p[0] for p in common], dtype=np.int64)
        genes_p = np.array([p[1] for p in common], dtype=np.int64)
        row, rowg = {}, {}
        for s, (f, _) in stages.items():
            if i1.size == 0:
                row[s] = {"n": 0}
                continue
            a1 = H[f][i1] - base[f][cells_p]
            a2 = H[f][i2] - base[f][cells_p]
            okp = okc[i1] & okc[i2]
            Lr = er.linearity(np.where(okp[:, None], a1, np.nan), np.where(okp[:, None], a2, np.nan), floors[s][cells_p])
            row[s] = lin_summary(Lr, okp)
            row[s]["n_excluded_empty"] = int((~okp).sum())
            for g, sym in enumerate(gene_syms):
                mg = okp & (genes_p == g)
                if mg.any():
                    rowg.setdefault(sym, {})[s] = lin_summary(Lr, mg)
        linear["gene_scaling"][f"eps{e_full:g}_vs_eps{e_half:g}"] = row
        linear["gene_scaling_per_gene"][f"eps{e_full:g}_vs_eps{e_half:g}"] = rowg

    # ---- decomposition (C3) diagnostics: probe already wrong at baseline among lost cases
    probe_base_wrong = base["probe"][cc] != k_case
    diag = {}
    for name, m in [("thinning", kind == KIND_THIN), ("gene_scaling", kind == KIND_GENE)]:
        lm = lost & m
        diag[name] = {"n_lost_cases": int(lm.sum()),
                      "share_lost_with_probe_wrong_at_baseline": float(np.mean(probe_base_wrong[lm])) if lm.any() else None}

    # ---- criterion (i) on the family x level units, and the fine gene units (secondary)
    def crit_i(U, loss_key):
        pts = {u: b["point"] for u, b in U.items()}
        bts = {u: b["_boot"] for u, b in U.items()} if W is not None else {u: {k: np.full(1, np.nan) for k in b["point"]}
                                                                           for u, b in U.items()}
        pdn = {d: {u: b["per_donor"][d] for u, b in U.items()} for d in donors}
        return er.criterion_i(pts, bts, loss_key=loss_key, loss_tol=LOSS_TOL, min_lost=MIN_LOST,
                              share_margin=SHARE_MARGIN, per_donor=pdn)

    crit1 = crit_i(units, "decision_loss")
    crit1_sel = crit_i(units, "selective_loss")
    crit1_fine = crit_i({**{u: b for u, b in units.items() if u.startswith("thin_")}, **units_fine}, "decision_loss")

    # ---- criterion (ii)
    def crit_ii():
        keep_small, keep_target = round(1 - EPS_SMALL, 6), round(1 - EPS_TARGET, 6)
        m_base = er.top_margin(base["S"])
        seeds = sorted(set(sd[kind == KIND_THIN].tolist()))
        ms = np.zeros(n)
        sc0 = np.zeros(n)
        cnt = np.zeros(n)
        for s in seeds:
            m = (kind == KIND_THIN) & np.isclose(lev, keep_small) & (sd == s) & okc
            ms[cc[m]] += er.signed_margin(H["S"][m], c0[cc[m]])
            sc0[cc[m]] += H["S"][m, c0[cc[m]]]
            cnt[cc[m]] += 1
        full = cnt == len(seeds)
        ms = np.where(full, ms / np.maximum(cnt, 1), np.nan)
        sc0 = np.where(full, sc0 / np.maximum(cnt, 1), np.nan)
        m_hat = er.extrapolate(m_base, ms, EPS_SMALL, EPS_TARGET)
        called0 = base["call"] != vk.NO_CALL
        s_c0 = base["S"][np.arange(n), c0]
        side = np.where(called0, 1.0, -1.0)
        b_base = side * (s_c0 - bar)
        b_hat = side * (er.extrapolate(s_c0, sc0, EPS_SMALL, EPS_TARGET) - bar)
        d_base = np.minimum(m_base, b_base)
        d_hat = np.minimum(m_hat, b_hat)
        mt = (kind == KIND_THIN) & np.isclose(lev, keep_target)
        ci_ = cc[mt]
        y = H["call"][mt] != base["call"][ci_]
        keep = np.isfinite(m_hat[ci_])
        ci_, y = ci_[keep], y[keep]
        scores = {"margin": -m_base[ci_], "extrapolated_margin": -m_hat[ci_], "distance": -d_base[ci_],
                  "extrapolated_distance": -d_hat[ci_]}
        groups = {k: er.auroc_groups(v) for k, v in scores.items()}
        pt = {k: er.auroc(v, y, groups=groups[k]) for k, v in scores.items()}
        comps = {"primary: extrapolated_margin - margin": ("extrapolated_margin", "margin"),
                 "secondary: extrapolated_distance - distance": ("extrapolated_distance", "distance"),
                 "secondary: extrapolated_distance - margin": ("extrapolated_distance", "margin")}
        res = {"n_cases": int(y.size), "n_changed": int(y.sum()), "n_cells_excluded_missing_small_eps": int((~full).sum()),
               "auroc": pt, "differences": {}}
        per_d = {}
        for d in donors:
            md = donor[ci_] == d
            per_d[d] = {k: er.auroc(v[md], y[md]) for k, v in scores.items()}
        res["auroc_per_donor"] = per_d
        for name, (p_, q_) in comps.items():
            point = pt[p_] - pt[q_]
            rep = np.full(W.shape[0] if W is not None else 0, np.nan)
            if W is not None:
                for b in range(W.shape[0]):
                    w = W[b, ci_].astype(np.float64)
                    rep[b] = er.auroc(scores[p_], y, w, groups[p_]) - er.auroc(scores[q_], y, w, groups[q_])
            lo, hi, nb = er.percentile_ci(rep) if W is not None else (None, None, 0)
            pdd = {d: (per_d[d][p_] - per_d[d][q_]) for d in donors}
            res["differences"][name] = {
                "point": point, "ci95": [lo, hi], "n_undefined_replicates": nb, "per_donor": pdd,
                "holds": bool(np.isfinite(point) and point >= AUROC_MARGIN and lo is not None and lo > 0),
                "per_donor_clause": bool(all(v is not None and np.isfinite(v) and v >= AUROC_MARGIN / 2 for v in pdd.values())),
                "common_rule": er.win_rule(point, [lo, hi], pdd, AUROC_MARGIN)}
        res["holds"] = res["differences"]["primary: extrapolated_margin - margin"]["holds"]
        return res

    crit2 = crit_ii()
    verdict = ("adds information beyond the curve" if (crit1["holds"] or crit2["holds"])
               else "adds nothing beyond the curve (falsified)")

    def strip(U):
        return {u: {k: v for k, v in b.items() if k != "_boot"} for u, b in U.items()}

    return {
        "n_cells": int(n), "donors": {d: int(np.sum(donor == d)) for d in donors},
        "n_cases": int(cc.size), "n_boot": int(n_boot),
        "baseline": {"q1_coverage": float(np.mean(base["call"] != vk.NO_CALL)),
                     **{k: v for k, v in vk.selective_accuracy(base["call"], key).items()}},
        "curve_comparator": strip(units), "curve_fine_gene_units": strip(units_fine), "floor_unit": strip({"floor": floor_blk}),
        "measurements": measurements, "anm_identity_max_abs_dANM_minus_3G3_dS": anm_identity,
        "linearity": linear, "decomposition_diagnostics": diag,
        "falsification": {"criterion_i": crit1, "criterion_i_selective_accuracy_secondary": crit1_sel,
                          "criterion_i_with_gene_units_secondary": {k: v for k, v in crit1_fine.items() if k != "rows"}
                          | {"rows_passing": [r for r in crit1_fine["rows"] if r["passes"]]},
                          "criterion_ii": crit2, "verdict": verdict,
                          "verdict_rule": reg["experiments"]["E2"]["falsification"]},
    }


def stage_report(a, reg: dict, info: dict) -> None:
    add = load_addendum(a)
    pack = load_pack(a.processed, need_rna=False)
    cells, role = cell_pool(a, reg, pack)
    ch = load_chunks(a, cells)
    zp_path = a.z_probe or (a.out_dir / "z_probe.npz")
    if not zp_path.exists():
        zp_path = Path(add["computed"]["z_probe"]["file"])
    with np.load(zp_path) as d:
        probe = {"coef": d["coef"], "intercept": d["intercept"], "classes": d["classes"].astype(str)}
    if coef_sha(probe["coef"], probe["intercept"], probe["classes"]) != add["computed"]["z_probe"]["coef_sha256"]:
        raise SystemExit(f"{zp_path}: z-probe coefficients differ from the addendum's sha256")
    inputs = {"head_ckpt_sha256": vk.sha256_file(a.head_ckpt), "z_rna_sha256": vk.sha256_file(a.embed_dir / "z_rna.npy"),
              "teddy_weights_sha256": vk.sha256_file(a.ckpt / "model.safetensors")}
    regin = reg["provenance"]["inputs"]
    if inputs["head_ckpt_sha256"] != regin["ckpt_sha256"] or inputs["z_rna_sha256"] != regin["z_sha256"]:
        raise SystemExit(f"head checkpoint or stored embedding differs from the registered inputs: {inputs}")
    say(f"report: {cells.size} cells, {ch['z'].shape[0]} cases; z-probe C {add['computed']['z_probe']['C']}")
    gene_syms = [g["symbol"] for g in add["computed"]["genes"]["kept"]]
    # verifier side: keys of the pool cells
    adt = adt_rows(pack, cells)
    keys = vk.build_keys(pack["cell_types"][cells], adt, pack["adt_names"], reg)
    key = np.asarray(keys["primary"]).astype(str)
    donor = pack["donors"][cells]
    pipe = Pipeline(a, reg, probe, pack)
    z_store = np.asarray(np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")[cells], dtype=np.float32)
    t0 = time.time()
    base = pipe.head(z_store)
    base_anm = pipe.anm(base["v12"], base["ok"])
    base["A"], base["anm_call"] = base_anm["A"], base_anm["call"]
    H = pipe.head(ch["z"].astype(np.float32))
    an = pipe.anm(H["v12"], H["ok"])
    H["A"], H["anm_call"] = an["A"], an["call"]
    say(f"head + ANM engine on {cells.size + ch['z'].shape[0]} rows in {time.time() - t0:.0f}s")
    pos = {int(c): i for i, c in enumerate(cells)}
    case = {"cell": np.array([pos[int(c)] for c in ch["var_cell"]], dtype=np.int64), "kind": ch["var_kind"].astype(int),
            "level": ch["var_level"], "seed": ch["var_seed"], "gene": ch["var_gene"], "identical": ch["var_identical"]}
    # ANM vs rule, cell by cell (baseline and every case)
    mism_base = int(np.sum(base_anm["call"] != base["call"]))
    mism_case = int(np.sum(H["anm_call"] != H["call"]))
    G3 = vk.field_gain(3, reg["anm"]["field_representation"])
    cf_err = float(np.nanmax(np.abs(H["A"] - 3 * G3 * H["S"])))
    checks = {
        "count_recovery_cells_flagged": int(np.sum(ch["cell_count_dev"] > COUNT_TOL)),
        "count_recovery_max_deviation_ratio": float(np.max(ch["cell_count_dev"])),
        "baseline_ntokens_mismatch_vs_official_run": int(np.sum(ch["cell_ntok_base"] != ch["cell_ntok_stored"])),
        "anm_engine_call_mismatches_baseline": mism_base, "anm_engine_call_mismatches_cases": mism_case,
        "anm_engine_vs_closed_form_max_abs": cf_err,
        "empty_variants": int(np.sum(~H["ok"])), "embed_runtime_sec": ch["embed_runtime_sec"],
        "floor_dz": _quant(np.linalg.norm(H["zl"][case["kind"] == KIND_FLOOR] - base["zl"][case["cell"][case["kind"] == KIND_FLOOR]], axis=1)),
        "floor_call_changes": int(np.sum(H["call"][case["kind"] == KIND_FLOOR] != base["call"][case["cell"][case["kind"] == KIND_FLOOR]])),
    }
    if mism_base or mism_case:
        say(f"WARN bridge failure: ANM engine calls differ from the rule on {mism_base} baseline cells and {mism_case} cases")
    say(f"checks {checks}")
    seed = int(reg["seeds"]["bootstrap"])
    pops = {}
    if a.cell_pool == "e2_subset":
        pops["test_primary"] = role == "primary"
        pops["test_secondary_15078"] = role == "secondary"
    else:
        pops["val_smoke"] = role == "val_smoke"
    results = {}
    for name, pm in pops.items():
        idx = np.where(pm)[0]
        remap = np.full(cells.size, -1)
        remap[idx] = np.arange(idx.size)
        cm = pm[case["cell"]]
        sub_case = {k: (v[cm] if k != "cell" else remap[v[cm]]) for k, v in case.items()}
        Hs = {k: v[cm] for k, v in H.items()}
        bs = {k: v[idx] for k, v in base.items()}
        t1 = time.time()
        results[name] = analyse(Hs, bs, sub_case, key[idx], donor[idx], gene_syms, reg, a.n_boot, seed, pipe.bar)
        say(f"{name}: analysed in {time.time() - t1:.0f}s; verdict: {results[name]['falsification']['verdict']}")
    out = {"experiment": "E2", "script": SCRIPT_VERSION, "smoke": a.smoke, "cell_pool": a.cell_pool,
           "registration_sha256": info["registration_sha256"],
           "amendment_sha256": info["amendment_sha256"],
           "addendum_sha256": info.get("addendum_sha256"), "spec_sha256": spec_sha(),
           "z_probe": {k: add["computed"]["z_probe"][k] for k in ("C", "C_grid", "coef_sha256")},
           "genes": {"kept": gene_syms, "dropped": [g["symbol"] for g in add["computed"]["genes"]["dropped"]]},
           "inputs": inputs, "checks": checks, "results": results,
           "primary_population": "test_primary" if a.cell_pool == "e2_subset" else "val_smoke",
           "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out = clean(out)
    write_json(a.out_dir / "E2_results.json", out)
    np.savez_compressed(a.out_dir / "e2_cases.npz", cells=cells, role=role, key=key, donor=donor,
                        case_cell=case["cell"], case_kind=case["kind"], case_level=case["level"], case_seed=case["seed"],
                        case_gene=case["gene"], case_identical=case["identical"], call=H["call"], probe=H["probe"],
                        argmax=H["argmax"], S=H["S"], v12=H["v12"], A=H["A"], base_call=base["call"],
                        base_S=base["S"], base_probe=base["probe"], base_v12=base["v12"])
    write_report(a.out_dir, out)
    say(f"wrote {a.out_dir / 'E2_results.json'} and REPORT.md")


def _quant(x) -> dict:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if not x.size:
        return {"n": 0}
    return {"n": int(x.size), "median": float(np.median(x)), "p95": float(np.quantile(x, 0.95)), "max": float(x.max())}


# ============================================================================ report markdown

def _f(x, n=3):
    return "n/a" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))


def _ci(c, n=3):
    if not c or c[0] is None:
        return "n/a"
    return f"[{c[0]:.{n}f}, {c[1]:.{n}f}]" + (f" ({c[2]})" if len(c) > 2 else "")


def write_report(out_dir: Path, R: dict) -> None:
    L = []
    L.append("# E2: response decomposition (C3) and linearity (C5)\n")
    if R.get("smoke"):
        L.append(f"**SMOKE RUN** ({R['smoke']}): val donor cells only, reduced bootstrap. Not a result.\n")
    L.append(f"Registration `{R['registration_sha256'][:16]}`, amendment A1 `{R['amendment_sha256'][:16]}`, "
             f"E2 addendum `{(R.get('addendum_sha256') or 'none')[:16]}`. z-probe C = {R['z_probe']['C']}. "
             f"Genes kept: {', '.join(R['genes']['kept'])}; dropped: {', '.join(R['genes']['dropped']) or 'none'}.\n")
    c = R["checks"]
    L.append("## Checks\n")
    cf = c["anm_engine_vs_closed_form_max_abs"]
    L.append(f"- ANM engine (finite_graph_scalar, events at t = 0) vs the Q1 rule: {c['anm_engine_call_mismatches_baseline']} "
             f"baseline and {c['anm_engine_call_mismatches_cases']} case mismatches (registered value 0); engine vs closed "
             f"form max |A - 3 G(3) S| = {'n/a' if cf is None else f'{cf:.1e}'}.")
    L.append(f"- Floor (unperturbed re-embedding vs stored z): median dz {_f(c['floor_dz'].get('median'), 5)}, p95 "
             f"{_f(c['floor_dz'].get('p95'), 5)}; call changes from the floor alone: {c['floor_call_changes']}.")
    L.append(f"- Count recovery flagged cells: {c['count_recovery_cells_flagged']}; baseline token-count mismatches vs the "
             f"official run: {c['baseline_ntokens_mismatch_vs_official_run']}; empty variants: {c['empty_variants']}.\n")
    for pop, res in R["results"].items():
        L.append(f"## {pop} ({res['n_cells']} cells, {res['n_cases']} cases; donors {res['donors']})\n")
        b = res["baseline"]
        L.append(f"Baseline Q1: coverage {_f(b['q1_coverage'])}, selective accuracy {_f(b['accuracy_of_calls'])}, "
                 f"decision accuracy {_f(b['decision_accuracy'])} (primary key).\n")
        L.append("### Comparator: perturbation -> accuracy curve, and E2.dec localisation\n")
        L.append("| unit | cases | decision acc (base -> pert) | loss [95% CI] | selective acc (pert) | coverage (pert) | "
                 "lost cells | representation | head | decision |")
        L.append("|---|---:|---|---|---:|---:|---:|---|---|---|")
        for u, blk in res["curve_comparator"].items():
            p, ci = blk["point"], blk.get("ci95", {})
            L.append(f"| {u} | {p['n_cases']} | {_f(p['decision_acc_base'])} -> {_f(p['decision_acc_pert'])} | "
                     f"{_f(p['decision_loss'])} {_ci(ci.get('decision_loss'))} | {_f(p['selective_acc_pert'])} | "
                     f"{_f(p['coverage_pert'])} | {p['n_lost_cells']} | "
                     + " | ".join(f"{_f(p['share_' + k])} {_ci(ci.get('share_' + k), 2)}" for k in er.LABELS) + " |")
        L.append("")
        L.append("### Measurements (median over cases)\n")
        L.append("| unit | dz | dv (clipped) | dv (unclipped) | dS called | dmargin | dANM called | call change rate | "
                 "identical tokens |")
        L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for u, m in res["measurements"].items():
            L.append(f"| {u} | {_f(m['dz']['median'], 4)} | {_f(m['dv_clipped']['median'], 4)} | "
                     f"{_f(m['dv_unclipped']['median'], 4)} | {_f(m['dS_called']['median'], 4)} | "
                     f"{_f(m['dmargin']['median'], 4)} | {_f(m['dANM_called']['median'], 4)} | {_f(m['dcall_rate'])} | "
                     f"{m['n_tokens_identical_to_baseline']} |")
        L.append("")
        L.append("### Linearity (C5): share of rows in the linear regime (|rho - 1| <= 0.1, cos >= 0.9, above 3 x floor)\n")
        L.append("| family | pair | stage | n | share linear | share below floor | median rho | median cos |")
        L.append("|---|---|---|---:|---:|---:|---:|---:|")
        for fam in ("thinning", "gene_scaling"):
            for pair, row in res["linearity"][fam].items():
                for st, s in row.items():
                    L.append(f"| {fam} | {pair} | {st} | {s.get('n')} | {_f(s.get('share_linear'))} | "
                             f"{_f(s.get('share_below_floor'))} | {_f(s.get('median_rho'))} | {_f(s.get('median_cos'))} |")
        L.append("")
        F = res["falsification"]
        c1, c2 = F["criterion_i"], F["criterion_ii"]
        L.append("### Falsification\n")
        L.append(f"**Verdict: {F['verdict']}.**\n")
        L.append(f"- (i) {c1['n_pairs_compared']} unit pairs qualify (decision-accuracy losses within {LOSS_TOL}, >= "
                 f"{MIN_LOST} lost cells each), {c1['n_tests']} (pair, label) tests; holds: {c1['holds']} (with the "
                 f"per-donor clause: {c1['holds_with_per_donor_clause']}; Bonferroni: {c1['holds_bonferroni']}).")
        for r in c1["rows"]:
            if r["passes"] or abs(r["diff"] or 0) >= SHARE_MARGIN:
                L.append(f"  - {r['unit_a']} vs {r['unit_b']}, {r['label']}: {_f(r['share_a'])} vs {_f(r['share_b'])}, "
                         f"diff {_f(r['diff'])} {_ci(r['ci95'])}; passes {r['passes']}")
        d = c2["differences"]["primary: extrapolated_margin - margin"]
        L.append(f"- (ii) {c2['n_cases']} cases, {c2['n_changed']} change call at eps 0.8. AUROC: margin "
                 f"{_f(c2['auroc']['margin'])}, extrapolated margin {_f(c2['auroc']['extrapolated_margin'])}; difference "
                 f"{_f(d['point'])} {_ci(d['ci95'])}; holds: {d['holds']} (per-donor clause {d['per_donor_clause']}; "
                 f"common rule: {d['common_rule']}).")
        for name, dd in c2["differences"].items():
            if name.startswith("secondary"):
                L.append(f"  - {name}: {_f(dd['point'])} {_ci(dd['ci95'])} ({dd['common_rule']})")
        L.append(f"- Rule: {F['verdict_rule']}\n")
    L.append("ANM here is the response-theory verifier and readout form: with every event at t = 0 its action score is "
             "G(3) x 3 x the rule score, so its calls and responses equal the declared rule's (checked above).\n")
    (out_dir / "REPORT.md").write_text("\n".join(L) + "\n")


# ============================================================================ main

def main(argv=None) -> int:
    a = parse_args(argv)
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG["path"] = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage} pool {a.cell_pool} out {a.out_dir}" + (f" SMOKE: {a.smoke}" if a.smoke else ""))
    if a.stage == "register":
        stage_register(a)
        return 0
    rc = check_registration(a)
    write_json(a.out_dir / "e2_registration_used.json", {**rc["info"], "spec": SPEC})
    if a.stage in ("embed", "all"):
        stage_embed(a, rc["reg"], rc["info"])
    if a.stage in ("report", "all"):
        stage_report(a, rc["reg"], rc["info"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
