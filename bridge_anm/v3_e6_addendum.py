#!/usr/bin/env python3
"""E6 addendum (registration/addenda/E6.json): everything E6 needs that the registration leaves to the addendum,
fixed before any external outcome (key, kappa, accuracy) is computed.

Contents (registration experiments.E6 + the decisions this builder had to take, each by a registered principle):

* the external annotation map (celltype.l2 -> class) with its sha256, and the independent re-derivation of the map
  from the 31 label names (``INDEPENDENT_MAP``, written before the committed map was opened) with every disagreement
  (there are none) and the four items prep flagged for confirmation;
* the contamination statement and its sources;
* the external pack (d5k: 5,000 label-free cells per donor) and its reason; the official embedding it needs;
* the gate with the missing-protein rule (exact clauses that remain), the per-donor estimator, the keys and the
  key-validation failure rule; the panels that remain;
* arms, endpoints, statistics, the replication rule and the sensitivity rows;
* one number fixed from train/val by a registered procedure: the C of the secondary same-evidence classifier
  (registered C grid, lowest val log-loss on val primary-key cells, ties -> smaller C);
* disclosures (prep-order deviation, what the orchestrator and this builder have seen).

Rows read: BMMC split train and val only (protein, cell types, embedding) for the secondary classifier's C. Of the
external data only metadata is read: external_manifest.json (label names, ADT name mapping, the committed map),
meta.json and prepare_stats.json donor/site counts (label-free), and the pack file's sha256 (bytes hashed, not parsed).
No external ADT value, embedding value, label count by prediction or site4 value is read.

Usage (repo root):
  python bridge_anm/v3_e6_addendum.py            # writes registration/addenda/E6.json + its HASHES.txt line
  python bridge_anm/v3_e6_addendum.py --print    # the document on stdout, nothing written
"""
from __future__ import annotations

import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

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
from lib import v3_e6 as e6  # noqa: E402
from lib import v3_key as vk  # noqa: E402

BUILDER_VERSION = "v3_e6_addendum 1.0"
ADDENDUM = ROOT / "registration" / "addenda" / "E6.json"
HASHES = ROOT / "registration" / "addenda" / "HASHES.txt"
EXT = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/data/processed/external_hao2021_d5k")
EXT_EMB = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/data/processed/external_hao2021_d5k_official")
PREP_REPORT = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/e6_external_prep/PREP_external_hao2021_d5k.md")
E1_SITE4_RESULTS = Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E1/E1_results.json")
COMMITTED_MAP_SHA = "5492a4176e1c661abefa3d287e32fa8dfcf2b12b6edfaeef93546e3c9f5e0660"

# Written 2026-10-01 by the E6 builder from the 31 celltype.l2 names (keys of the pack manifest's class map, read
# without their values) and registration_v3.json annotation_map / annotation_map_notes, before the committed map was
# opened; then compared with it (31 of 31 agree; same sha256).
INDEPENDENT_MAP = {
    "ASDC": "OUT", "B intermediate": "B", "B memory": "B", "B naive": "B", "CD14 Mono": "myeloid",
    "CD16 Mono": "myeloid", "CD4 CTL": "T", "CD4 Naive": "T", "CD4 Proliferating": "T", "CD4 TCM": "T",
    "CD4 TEM": "T", "CD8 Naive": "T", "CD8 Proliferating": "T", "CD8 TCM": "T", "CD8 TEM": "T",
    "Doublet": "unscored", "Eryth": "OUT", "HSPC": "OUT", "ILC": "OUT", "MAIT": "T", "NK": "NK",
    "NK Proliferating": "NK", "NK_CD56bright": "NK", "Plasmablast": "OUT", "Platelet": "OUT", "Treg": "T",
    "cDC1": "myeloid", "cDC2": "myeloid", "dnT": "T", "gdT": "T", "pDC": "OUT",
}
INDEPENDENT_REASONS = {
    "B intermediate, B memory, B naive": "mature B; registered principle: naive / transitional / B1 B -> B",
    "CD14 Mono, CD16 Mono, cDC1, cDC2": "monocytes and conventional DCs; registered: CD14+ / CD16+ Mono, cDC1, cDC2 -> myeloid",
    "CD4 CTL, CD4 Naive, CD4 TCM, CD4 TEM, CD8 Naive, CD8 TCM, CD8 TEM, Treg, MAIT, gdT, dnT":
        "T lineage; registered: CD4+ / CD8+ T, T reg, MAIT, gdT TCRVD2+ / gdT CD158b+, dnT -> T",
    "CD4 Proliferating, CD8 Proliferating": "lineage-committed cycling T; registered: T prog cycling -> T",
    "NK, NK_CD56bright": "conventional NK (CD56bright CD16- is an NK subset, not an ILC); registered: NK, NK CD158e1+ -> NK",
    "NK Proliferating": "lineage-committed cycling NK; by the registered cycling principle (T prog cycling -> its lineage) -> NK",
    "ILC": "registered: ILC / ILC1 -> OUT (no TCR, not conventional NK)",
    "Plasmablast": "registered: plasma cells / plasmablasts -> OUT (terminal B state outside the mature-B question)",
    "HSPC": "stem / progenitor; registered: HSC, Lymph prog, G/M prog, MK/E prog -> OUT",
    "Eryth": "erythroid; registered: Proerythroblast, Erythroblast, Normoblast, Reticulocyte -> OUT",
    "Platelet": "megakaryocyte lineage (anucleate fragments); registered: MK/E prog -> OUT; no B / T / NK / myeloid call is correct",
    "pDC": "registered: pDC -> OUT (not a monocyte or conventional DC)",
    "ASDC": "AXL+ SIGLEC6+ DC: neither a monocyte nor a conventional DC (the registered myeloid class is monocytes + "
            "cDC1/cDC2) and CD123+ pDC-like; by the registered pDC principle -> OUT",
    "Doublet": "not a cell type: none of the five classes applies (OUT would assert a single cell outside the four "
               "lineages) -> unscored in every key, kept in coverage denominators like any unscored cell",
}
FLAGGED_FOR_OWNER = {
    "Platelet": {"committed": "OUT", "independent": "OUT", "decision": "OUT",
                 "principle": "registered MK/E prog -> OUT: a megakaryocyte-lineage population outside the four lineages"},
    "ASDC": {"committed": "OUT", "independent": "OUT", "decision": "OUT",
             "principle": "registered pDC -> OUT: myeloid = monocytes and conventional DCs only"},
    "NK Proliferating": {"committed": "NK", "independent": "NK", "decision": "NK",
                         "principle": "registered T prog cycling -> T: a lineage-committed cycling population keeps its lineage"},
    "Doublet": {"committed": "unscored", "independent": "unscored", "decision": "unscored",
                "principle": "not a cell type; the five-class key has no class for it (registered: cells without a key "
                             "class are unscored and stay in coverage denominators)"},
}

CONTAMINATION = {
    "label": e6.CONTAMINATION_LABEL,
    "answer": "likely yes",
    "evidence": [
        "TEDDY paper (arXiv 2503.03485, section on data): the training corpus is derived from CELLxGENE, 1,399 "
        "single-cell RNA-seq datasets downloaded on June 11, 2024 (160M cells, 70M from primary datasets); QC removed "
        "cells with fewer than 225 gene counts or more than 10% mitochondrial reads and excluded 10x Chromium v1 "
        "studies; 116M cells remained. Held out: five named disease datasets (chronic kidney disease, Alzheimer's, "
        "gastric cancer, rheumatoid arthritis, pediatric Crohn's) and 82 unnamed donors for a held-out-donor task.",
        "CELLxGENE curation API (collections/b0cf0afa-ec40-4d65-b570-ed4ceacc6813, read 2026-10-01): collection "
        "'Integrated analysis of multimodal single-cell data' (doi 10.1016/j.cell.2021.04.048), published "
        "2022-07-15; dataset 'nygc multimodal pbmc' (ed5d841d-6346-47d4-ab2f-7119ad7e3a35), 161,764 cells, all "
        "is_primary_data, assay 10x 3' v3 - the cells of GSE164378's 3' set.",
        "Merck/TEDDY model card (huggingface.co/Merck/TEDDY): 116 million cells; no dataset list or exclusion of this set.",
    ],
    "reading": "the dataset was public in CELLxGENE two years before TEDDY's download, is 10x v3 (not excluded) and is "
               "not among the named held-out datasets; unless its donors happen to be among the 82 unnamed held-out "
               "donors, TEDDY-G's pretraining saw these cells' RNA. Our heads, thresholds, panels and rules were never "
               "fit on it. Every compared arm reads the same TEDDY embedding, so exposure does not favour one arm, but "
               "absolute accuracies may be optimistic. Every E6 output carries the label.",
    "sources": ["https://arxiv.org/html/2503.03485", "https://api.cellxgene.cziscience.com/curation/v1/collections/b0cf0afa-ec40-4d65-b570-ed4ceacc6813",
                "https://huggingface.co/Merck/TEDDY"],
    "checked_utc": "2026-10-01",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def content_bytes(d: dict) -> bytes:
    return (json.dumps(d, indent=1, ensure_ascii=False) + "\n").encode()


def _git(*args) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    except OSError:
        return None


# ============================================================================ external metadata (no values)
def external_metadata(ext: Path) -> dict:
    man = json.loads((ext / "external_manifest.json").read_text())
    meta = json.loads((ext / "meta.json").read_text())
    stats = json.loads((ext / "prepare_stats.json").read_text())
    names = [r["ours"] for r in man["adt_mapping"]]
    present = [r["ours"] for r in man["adt_mapping"] if r["source"]]
    absent = [r["ours"] for r in man["adt_mapping"] if not r["source"]]
    if len(present) != int(meta["n_adt_present"]):
        raise SystemExit(f"manifest maps {len(present)} proteins, meta says {meta['n_adt_present']} present")
    amap = man["annotation_map_addendum"]["map"]
    if e6.map_sha256(amap) != man["annotation_map_addendum"]["sha256_canonical_json"]:
        raise SystemExit("committed map does not hash to its recorded sha256")
    return {"manifest": man, "meta": meta, "donors": stats["donors"], "sites": stats["sites"],
            "adt_names": names, "present": present, "absent": absent, "map": amap,
            "absent_rules": {r["ours"]: r["rule"] for r in man["adt_mapping"] if not r["source"]},
            "two_clone": [r["ours"] for r in man["adt_mapping"] if r["rule"].startswith("mean of two clones")]}


# ============================================================================ train/val number (secondary classifier C)
def classifier10_core(args, reg, features: list[str]) -> dict:
    """Registered classifier procedure (section 8) on the remaining panel proteins: training-cell primary key,
    C = lowest val log-loss on val primary-key cells over the registered grid (ties -> smaller C). Train/val rows only."""
    from sklearn.metrics import log_loss

    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    split = npz["split"].astype(str)
    names = [str(x) for x in npz["adt_names"]]
    nt = np.where(split != "test")[0]                      # the only rows whose protein / labels / z are read
    adt = np.asarray(npz["adt"], dtype=np.float32)[nt]
    ct = npz["cell_types"].astype(str)[nt]
    zmm = np.load(args.z, mmap_mode="r")
    z = np.asarray(zmm[nt], dtype=np.float32)
    sp = split[nt]
    tr, vl = sp == "train", sp == "val"
    j = {n: i for i, n in enumerate(names)}
    prim = vk.build_keys(ct, adt, names, reg)["primary"]
    v = vk.evidence(vk.head_predict(z, args.ckpt, device="cpu", size_factor=1.0), reg)
    X = v[:, [j[p] for p in features]]
    mt, mv = tr & (prim != vk.UNSCORED), vl & (prim != vk.UNSCORED)
    grid = list(reg["classifier"]["primary"]["C_grid"])
    res = []
    for C in grid:
        mdl = e1.fit_logreg(X[mt], prim[mt], C)
        P = e1.proba5(mdl, X[mv])[:, [vk.CLASSES.index(str(c)) for c in mdl.classes_]]
        res.append({"C": C, "val_log_loss": round(float(log_loss(prim[mv], P, labels=mdl.classes_)), 6)})
    best = min(res, key=lambda r: (r["val_log_loss"], r["C"]))["C"]
    # the frozen 12-feature classifier must reproduce its registered val log-loss (same pipeline)
    f12 = list(reg["classifier"]["primary"]["feature_proteins"])
    m12 = e1.fit_logreg(v[mt][:, [j[p] for p in f12]], prim[mt], reg["classifier"]["primary"]["C"])
    P12 = e1.proba5(m12, v[mv][:, [j[p] for p in f12]])[:, [vk.CLASSES.index(str(c)) for c in m12.classes_]]
    ll12 = round(float(log_loss(prim[mv], P12, labels=m12.classes_)), 6)
    reg_ll12 = {r["C"]: r["val_log_loss"] for r in reg["classifier"]["primary"]["val_log_loss"]}[reg["classifier"]["primary"]["C"]]
    if abs(ll12 - reg_ll12) > 1.5e-6:
        raise SystemExit(f"registered classifier val log-loss not reproduced: {ll12} vs {reg_ll12}")
    log(f"classifier10 C grid {res} -> C {best}; registered 12-feature val log-loss reproduced ({ll12})")
    return {"features": features, "labels": "training-cell primary key (5 classes incl. OUT; BMMC global gate)",
            "C_grid": grid, "val_log_loss": res, "C": best,
            "C_selection": "lowest val log-loss on val primary-key cells, ties -> smaller C (registration classifier.primary.C_selection)",
            "n_train_labelled": int(mt.sum()), "n_val_labelled": int(mv.sum()),
            "registered_12_feature_val_log_loss_reproduced": {"registered": reg_ll12, "recomputed": ll12},
            "rows_read": "BMMC split train and val only"}


# ============================================================================ document
def build(args) -> dict:
    reg = va.load_registration_amended(args.registration)
    E6 = reg["experiments"]["E6"]
    ext = external_metadata(args.external)
    if sorted(ext["map"]) != sorted(INDEPENDENT_MAP):
        raise SystemExit("label sets differ between the committed map and the independent re-derivation")
    disagreements = [{"label": k, "committed": ext["map"][k], "independent": INDEPENDENT_MAP[k]}
                     for k in sorted(INDEPENDENT_MAP) if ext["map"][k] != INDEPENDENT_MAP[k]]
    final_map = dict(sorted(ext["map"].items()))
    for d in disagreements:   # none expected; any would be decided by principle here before the document is written
        raise SystemExit(f"map disagreement must be decided by principle before the addendum is written: {d}")
    map_sha = e6.map_sha256(final_map)
    if map_sha != COMMITTED_MAP_SHA or e6.map_sha256(INDEPENDENT_MAP) != COMMITTED_MAP_SHA:
        raise SystemExit("map sha256 differs from the committed 5492a417...0660")

    gate = reg["gate"]
    spec = e6.reduced_gate_spec(gate, ext["present"])
    pan_reg = vk.question_panel(reg, "Q1")
    pan = e6.reduced_panel(pan_reg, ext["present"])
    panel_absent = sorted({p for k in vk.LINEAGES for p in pan_reg[k] if p not in ext["present"]})
    anchors = vk.question_panel(reg, "Q3")
    if e6.reduced_panel(anchors, ext["present"]) != anchors:
        raise SystemExit("a Q3 anchor is not measured externally")
    f10 = [p for k in vk.LINEAGES for p in pan[k]]
    clf10 = classifier10_core(args, reg, f10)
    fld = reg["anm"]["field_representation"]
    thr = {"Q1": vk.anm_readout_threshold(reg, "Q1"), "Q2": vk.anm_readout_threshold(reg, "Q2"),
           "Q3": vk.anm_readout_threshold(reg, "Q3")}
    gains = {str(n): round(vk.field_gain(n, fld), 6) for n in (1, 2, 3)}
    task_said = {"CD335": "said absent in the task text; measured externally (exact name)",
                 "CD71": "said absent in the task text; measured externally (exact name)",
                 "CD33": "not named in the task text; NOT measured externally (not in Hao's 3' panel)"}

    reg_sha = vk.sha256_file(args.registration)
    am = {f: vk.sha256_file(args.registration.parent / f) for f in ("amendment_A1.json", "amendment_A2.json", "amendment_A3.json")}
    pack = args.external / "cite_arrays.npz"
    doc = {
        "addendum_to": "registration/registration_v3.json experiments.E6, as amended by A1, A2 and A3",
        "experiment": "E6",
        "addendum_version": 1,
        "title": "E6: external confirmation on Hao et al. 2021 PBMC CITE-seq (GEO GSE164378, 3' data) - " + e6.CONTAMINATION_LABEL,
        "registration_sha256": reg_sha,
        "amendments_sha256": am,
        "registered_E6": E6,
        "fixed_before_external_outcomes": "every choice here was written from the registration, the external label "
            "names and ADT names, and BMMC train/val cells, and committed before any external key, kappa, accuracy or "
            "other outcome was computed; smoke runs of the E6 runner use val donor 18303 cells only",
        "contamination": CONTAMINATION,
        "external_data": {
            "dataset": ext["meta"]["dataset"],
            "pack": str(pack), "pack_sha256": vk.sha256_file(pack),
            "manifest_sha256": vk.sha256_file(args.external / "external_manifest.json"),
            "meta_sha256": vk.sha256_file(args.external / "meta.json"),
            "prep_commit": ext["manifest"]["git_commit"], "prep_script": ext["manifest"]["script"],
            "n_cells": int(ext["meta"]["n_cells"]), "donors": ext["donors"], "sites": ext["sites"],
            "cell_selection": "d5k: 5,000 cells per donor, uniformly at random within donor (numpy default_rng(0)), "
                              "chosen by donor id only (label-free); 40,000 cells, 8 donors (hao2021_P1..P8)",
            "why_d5k": "compute: the official TEDDY-G embedding ran at 117.6 ms per cell on MPS for our 90,261 cells, and "
                       "this set tokenises longer (median about 1,950 tokens per cell, 46% at the 2,048 cap, vs 1,272 and "
                       "15% on our pack), so the full 161,764 cells would need more than 5.3 h of a shared GPU against "
                       "about 1.3 h or more for 40,000; 5,000 per donor keeps every donor equally weighted for the "
                       "donor-level replication rule. No label or outcome entered the choice.",
            "gene_universe": "train (the 12,052 training genes, as the head was fit on embeddings over that universe)",
            "time_points": "each donor contributes 3 time points (metadata 0 / 2 / 7); they are pooled within donor and "
                           "the bootstrap resamples donors (8)",
            "adt_absent_externally": ext["absent"],
            "adt_absent_reasons": ext["absent_rules"],
            "n_adt_measured": len(ext["present"]),
            "adt_two_clone_means": ext["two_clone"],
            "adt": "per-cell CLR log1p(c / g_i), g_i over all 228 measured antibodies; NaN = not measured; CD3, CD4, "
                   "CD56, CD45, CD11b, CD38, CD44, CD26 are the mean of two clones' raw counts (prep rule, data-free)",
            "embedding": {
                "dir": str(args.external_embedding),
                "required": "all shards (8 = ceil(40,000 / 5,024)) under shards/ and the assembled z_rna.npy + "
                            "z_rna_manifest.json, written by scripts/03_embed_rna.py; the runner refuses otherwise",
                "manifest_must_match": {"preprocessing": "official", "seq_len": 2048, "pooling": "gene-mean",
                                        "normalize_total": 10000.0, "autocast": "fp16", "length_buckets": True,
                                        "cells": "all", "n_cells": int(ext["meta"]["n_cells"]), "n_shards": 8,
                                        "medians_sha256": "same as the BMMC official embedding (data/processed/cite_official/z_rna_manifest.json)",
                                        "processed_input": str(args.external), "each shard cfg_hash": "= manifest config_hash"},
            },
        },
        "annotation_map": {
            "map": final_map, "sha256_canonical_json": map_sha,
            "canonical_json": "json.dumps(map, sort_keys=True, separators=(',', ':'))",
            "n_labels": len(final_map),
            "independent_rederivation": {
                "how": "the E6 builder listed the 31 celltype.l2 names (keys of the manifest's class map, values not "
                       "printed), wrote a map from registration_v3.json annotation_map and annotation_map_notes only, "
                       "saved it, and then opened the committed map",
                "map": INDEPENDENT_MAP, "sha256_canonical_json": e6.map_sha256(INDEPENDENT_MAP),
                "reasons": INDEPENDENT_REASONS,
                "n_agree": sum(INDEPENDENT_MAP[k] == ext["map"][k] for k in INDEPENDENT_MAP),
                "disagreements": disagreements,
                "flagged_by_prep_for_owner_confirmation": FLAGGED_FOR_OWNER,
            },
            "doublet": "Doublet -> unscored in every key (annotation, primary, annotation_only); Doublet cells stay in "
                       "every coverage denominator (registration: coverage over all cells, unscored included)",
        },
        "key": {
            "gate_proteins_registered": list(gate["proteins"]),
            "gate_proteins_absent_externally": spec["absent_gate_proteins"],
            "gate_proteins_used": spec["proteins_used"],
            "missing_protein_rule": "registered E6.missing_proteins, read as: a condition on an absent protein is removed "
                                    "from its clause (the registered example: CD335 absent -> NK count over the available "
                                    "NK markers, nk_min 2 kept); a clause whose defining condition is gone is dropped "
                                    "(the registered example: CD71 absent -> erythroid clause dropped)",
            "task_text_check": task_said,
            "clauses_remaining_in_order": spec["clauses"],
            "clauses_dropped": spec["dropped_clauses"],
            "conditions_removed": spec["removed_conditions"],
            "clauses_plain": [
                "1. OUT (erythroid) if CD71 high (kept: CD71 measured)",
                "2. T if CD3 high and CD19 not high and CD14 not high (unchanged)",
                "3. B if CD3 not high and CD19 high and CD14 not high (unchanged; b_markers = [CD19] as registered)",
                "4. NK if CD3, CD19, CD14 not high and at least 2 of CD56, CD335, CD16 high (CD94 removed from the NK "
                "count, nk_min 2 kept; the 'CD33 not high' condition removed)",
                "5. myeloid if CD3 and CD19 not high and any of CD14, CD11c high (CD33 removed from the any-of set)",
                "6. OUT otherwise",
            ],
            "consequence_note": "without CD33, a cell high only for CD33 can no longer gate myeloid (it gates OUT), and the "
                                "NK gate no longer excludes CD33-high cells; without CD94, CD56bright CD16- NK cells need "
                                "CD335 to reach 2 NK markers",
            "thresholds": "registered estimator gmm_nonzero_post50 (v3_key.gmm_nonzero_threshold: 2-component GMM on "
                          "the donor's non-zero CLR values, first grid point between the means where the upper "
                          "component's posterior reaches 0.5; random_state = gate.threshold_seed 0; all of the donor's "
                          "cells, no subsampling), estimated per measured gate protein within each external donor, "
                          "label-free; high = CLR > threshold (strict); NaN never high. BMMC thresholds are not used.",
            "keys": {"primary": "annotation class and per-donor gated class agree -> that class; else unscored",
                     "annotation_only": "annotation class (Doublet unscored)"},
            "validation": {"statistic": "Cohen's kappa (5 classes) of annotation class vs per-donor gated class, pooled "
                                        "over the 8 donors, on cells with an annotation class (Doublet excluded); per-donor "
                                        "kappas are descriptive",
                           "rule": "kappa < 0.70 -> every E6 result is reported on the annotation-only key and labelled "
                                   "'key not validated' (decision key = annotation_only); otherwise decision key = primary",
                           "both_keys_always_reported": True},
            "not_used": "the registered sensitivity key primary_no_gdT158 has no external analogue (Hao has one 'gdT' "
                        "label) and the Q3 key is not needed (no Q3 accuracy endpoint)",
        },
        "evidence": {
            "embedding": "official TEDDY-G gene-mean z of the external cells (z_rna.npy above), L2-normalised",
            "head": "phase-1 MLP + NB decoder mean, outputs/cite_phase1_official/best.pt (frozen; sha256 recorded by the runner)",
            "normaliser": "the registered training q95 per protein (registration evidence.teddy_head.q95_train_pred), clip [0, 1]",
            "measured_protein": "no measured external protein enters the evidence or any arm (key only)",
        },
        "panels": {
            "registered_primary": pan_reg,
            "absent_externally": panel_absent,
            "remaining_E6_primary": pan,
            "rule": "registered E6.missing_proteins: a class whose panel protein is absent is scored on its remaining "
                    "panel proteins: T on CD3, CD2 (CD5 absent); NK on CD122, CD56 (CD94 absent); B and myeloid unchanged",
            "bars": {"Q1": reg["questions"]["Q1"]["bar"], "Q2": reg["questions"]["Q2"]["bar"],
                     "Q3": reg["questions"]["Q3"]["bar"], "note": "registered bars, unchanged (the score stays a mean in [0, 1])"},
            "q3_anchors": anchors, "q3_note": "all four anchors are measured; Q3 is used only for E1.1a",
            "platelet_note": "A1 review: CD62P (myeloid panel) is a platelet protein and Platelet is keyed OUT, so myeloid "
                             "calls on platelets are a predictable error mode; the per-type call table reports them; the "
                             "panel is not changed",
        },
        "arms": {
            "rule": "mean of the evidence over each class's remaining panel; argmax (B, T, NK, myeloid order on ties); "
                    "call if the top score reaches the registered bar; confidence = top score",
            "anm": "ANM finite_graph_scalar engine via finite_field_runner (ANM_ROOT), bridge_anm/lib/v3_e1.AnmBridge "
                   "unchanged: one support event per remaining panel protein at t = 0 (B 3, T 2, NK 2, myeloid 3 events), "
                   "registered field and closure weights; readout thresholds as registered (Q1 %.6f, Q2 %.6f = G(3) * 3 * bar; "
                   "Q3 %.6f); confidence = top action score; ANM's own argmax" % (thr["Q1"], thr["Q2"], thr["Q3"]),
            "anm_unequal_panels": "with unequal panel sizes ANM's action score is G(n_k) * n_k * S_k (G(1) %s, G(2) %s, "
                                  "G(3) %s), so ANM is no longer the mean rule: a 2-event class scores about 2/3 of a 3-event "
                                  "class at the same mean. The registered equivalence (registration anm.equivalence) assumes "
                                  "equal panel sizes, which the registered missing-protein rule breaks; ANM is run as frozen, "
                                  "not re-weighted, and the full-panel sensitivity shows the equal-size case"
                                  % (gains["1"], gains["2"], gains["3"]),
            "closure": "ANM closure readout (readout_coordinates over each action's event sites, registered weights); "
                       "closure bars Q1 %s, Q2 %s; confidence = top closure score; closure argmax"
                       % (reg["anm"]["closure_bar_Q1"], reg["anm"]["closure_bar_Q2"]),
            "margin": "top-1 minus top-2 mean-rule score; ranks the rule's argmax calls",
            "entropy": "1 - H(q)/log 4, q = S / sum(S) of the mean-rule scores (0 when sum S = 0); ranks the rule's argmax calls",
            "classifier": "frozen registered classifier: logistic regression on the 12 registered primary-panel features "
                          "(C = 100, lbfgs, max_iter 3000, seed 0), refit exactly as E1 on the BMMC training-cell primary "
                          "key (58,657 cells, checked); it reads CD5 and CD94 evidence (predicted from RNA, exists for all "
                          "cells) that the reduced rule panels no longer read - an asymmetry in the classifier's favour; "
                          "its own lineage argmax, confidence = that probability; bars Q1 %s, Q2 %s for operating points"
                          % (reg["classifier"]["primary"]["bars"]["Q1"], reg["classifier"]["primary"]["bars"]["Q2"]),
            "classifier10": "secondary, same evidence as the rule: the registered classifier procedure on the 10 remaining "
                            "panel proteins (C fixed below from BMMC train/val); matched-coverage rows only (no bar)",
        },
        "endpoints": {
            "E1.1a": {
                "bridge": "registered E1.1a for E6: engine calls vs the exact re-coded ANM on the same panel (each class a "
                          "star with its own n_k events, v3_key.field_star; same thresholds): Q1, Q2, Q3 and closure Q1, Q2; "
                          "registered value 0 per donor and pooled; any mismatch or rejected event -> E6 not interpretable "
                          "(registered failure: stop)",
                "vs_mean_rule": "engine calls vs the mean rule's calls (Q1, Q2): count and share per donor and pooled; "
                                "expected nonzero under the reduced panels (above); reported, not a bridge failure. In the "
                                "full-panel sensitivity the two checks coincide and 0 is required",
            },
            "Q1_selective_accuracy_at_matched_coverage": "every arm at each registered grid coverage (0.95 ... 0.50) and at "
                "c* = the mean rule's realised Q1 coverage at its bar on the evaluated cells (recomputed per replicate), "
                "decision key and the other key, per donor and pooled; difference rows (E1.4b): anm - margin and "
                "anm - classifier at 0.90 and 0.70; deployed operating points (bars from val) descriptive",
            "E1.4a": "AURC (trapezoid over the registered grid / 0.45). Primary: anm - margin (registered margin 0.005). "
                     "Secondary: anm - classifier, anm - entropy, rule - margin (E1's 'top' form), anm - rule, "
                     "anm - classifier10",
            "E1.C2": "OUT decline (share of key-OUT cells not among the selected cells) at matched coverage c in {c*, 0.80, "
                     "0.70}: ANM closure - mean rule (registered margin 0.05) with the in-scope selective accuracy "
                     "difference beside it (guard >= -0.005, as E1); classifier - rule secondary",
            "per_donor": "every endpoint is reported per external donor (point) and pooled (point + two-stage interval)",
            "not_run": "E1.1b-d, E1.3, E1.5, Q2 / Q3 accuracy endpoints, the MLP, E1's sensitivity panels",
        },
        "statistics": {
            "bootstrap": "two-stage: draw 8 donors with replacement from the 8 (sorted), then cells with replacement within "
                         "each drawn donor; B = 2000; numpy default_rng(seeds.bootstrap = %d); every endpoint recomputed per "
                         "replicate, including the matched-coverage selection and c* (A1.1); percentile 95%% interval; "
                         "undefined replicates excluded and counted (v3_e1.replicate_counts / percentile_ci)"
                         % int(reg["seeds"]["bootstrap"]),
            "tie_break": "A1.1: rank of the external row index (0..39,999) in numpy default_rng(29).permutation(90261) "
                         "(v3_amend.tie_break_rank), the same for every arm; label-free",
            "replication_rule": "registered: a v3 conclusion replicates if its sign holds in a majority of the external donors "
                                "(at least 5 of 8 donor points with the v3 sign; an undefined donor value counts against) "
                                "and the pooled two-stage 95% interval excludes 0 on the same side",
            "v3_direction": "the sign of the E1 site4 point difference on E1's primary split (test_primary) and primary key, "
                            "read by the runner at report time from outputs/v3/E1/E1_results.json (its sha256 is recorded; a "
                            "smoke E1 file is refused in the real run). The same direction is used for both E6 keys. A v3 "
                            "point of exactly 0 gives 'no v3 direction' (the E6 interval is reported)",
            "replication_rows": {
                "E1.1a": "E1 splits.primary.E1.1a.passed; replicates if the E6 bridge check has 0 mismatches in every donor",
                "E1.4a anm - margin AURC": "E1 E1.4.keys.primary.top_vs_margin.E1.4a_AURC (E1's 'top' = rule = ANM)",
                "E1.4a anm - classifier AURC": "E1 top_vs_classifier.E1.4a_AURC (secondary)",
                "E1.4a anm - entropy AURC": "E1 top_vs_entropy.E1.4a_AURC (secondary)",
                "E1.4a rule - margin AURC": "same E1 number as anm - margin (secondary: E1's form on the E6 mean rule)",
                "Q1 acc anm - margin @0.90 / @0.70": "E1 top_vs_margin.E1.4b_acc@0.90 / @0.70",
                "Q1 acc anm - classifier @0.90 / @0.70": "E1 top_vs_classifier.E1.4b_acc@0.90 / @0.70 (secondary)",
                "E1.C2 closure - rule OUT decline @c*, 0.80, 0.70": "E1 C2.keys.primary.{c_star, 0.80, 0.70}.closure_minus_rule_out_decline",
            },
            "e1_decision_rule_beside": "the registered E1 decision rule (win / loss / equivalent / inconclusive with the "
                                       "endpoint's margin, 'each primary donor' read as each of the 8 external donors) is "
                                       "reported beside the replication result and never replaces it",
        },
        "sensitivity": {
            "full_registered_panels": "rule, ANM (engine, 3 events per class), closure, margin on the registered 12-protein "
                                      "panel including CD5 and CD94 evidence (predicted; no measured protein enters evidence): "
                                      "E1.1a, AURC anm - margin, C2 closure - rule at its own c*; never changes a verdict",
            "annotation_only_key": "always reported beside the decision key",
        },
        "computed": {"classifier10": clf10,
                     "anm_readout_thresholds": thr, "field_gain": gains},
        "disclosures": {
            "prep_order_deviation": "the external map (scripts/prepare_external_cite.py L2_CLASS, prep commit "
                                    f"{ext['manifest']['git_commit'][:7]}) was written from label names, but its first version "
                                    "followed the answer-key prototype (plasma -> B, pDC -> myeloid, ILC -> NK) and was revised "
                                    "to the registered principles after the prep agent's smoke checks had read external values "
                                    "(per-class medians of 8 measured proteins on 2,000 cells, ADT totals, token counts), and it "
                                    "was committed with the code that reads external values (A1 review note). The registered "
                                    "order 'committed before any external ADT or RNA value is read' was therefore not met in "
                                    "time. No value changed a mapping: the independent re-derivation above, from names only, "
                                    "reproduces the committed map exactly (31 of 31, same sha256).",
            "orchestrator_has_seen": "the orchestrating session has already seen site4 results of E1, E3 and E4. E6 takes no "
                                     "number from them; E1's site4 point differences enter only as the v3 direction of the "
                                     "replication rule, read by the runner at report time.",
            "this_builder_read": "label names (31) before writing the map; then the committed map, the prep report "
                                 "PREP_external_hao2021_d5k.md (per-donor annotation-class counts, family counts inside OUT, "
                                 "ADT coverage by name, gene, UMI and token statistics, integrity summaries with aggregate "
                                 "ADT-total ratios; no per-protein ADT value and no prediction), the manifest's ADT name mapping, prepare_stats.json donor and site counts, "
                                 "and the timing lines of E1's progress log (engine and bootstrap seconds); its smoke runs (val donor 18303) read "
                                 "E1's smoke results (val) to exercise the replication code, and it checked that the E1 "
                                 "site4 result paths the runner reads exist (types only, no value printed). It read no "
                                 "external ADT value, embedding value, label-vs-prediction comparison or site4 outcome.",
            "task_text_vs_panel": "the task text named CD335 and CD71 as absent (the registration's examples); both are "
                                  "measured. The absent gate proteins are CD94 and CD33, the absent panel proteins CD5 and "
                                  "CD94; the registered rule is applied to the actual absent set.",
            "contamination": e6.CONTAMINATION_LABEL,
        },
        "outputs": "/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E6/: E6_results.json, REPORT.md, key_validity.json, "
                   "progress.log (smoke: outputs/v3/E6_smoke)",
        "runner": "scripts/v3_e6_external.py (library bridge_anm/lib/v3_e6.py)",
        "provenance": {
            "builder": BUILDER_VERSION, "builder_sha256": vk.sha256_file(Path(__file__)),
            "inputs": {"bmmc_cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz"),
                       "bmmc_z_sha256": vk.sha256_file(args.z), "ckpt_sha256": vk.sha256_file(args.ckpt),
                       "prep_report_sha256": vk.sha256_file(PREP_REPORT) if PREP_REPORT.exists() else None},
            "git_head_at_build": _git("rev-parse", "HEAD"),
            "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "rows_read": "BMMC train + val (classifier10 C); external: metadata only",
        },
    }
    return doc


def write(doc: dict, path: Path = ADDENDUM, hashes: Path = HASHES) -> str:
    b = content_bytes(doc)
    path.write_bytes(b)
    sha = hashlib.sha256(b).hexdigest()
    lines = [ln for ln in (hashes.read_text().splitlines() if hashes.exists() else []) if not ln.rstrip().endswith("  E6.json")]
    hashes.write_text("\n".join(lines + [f"{sha}  E6.json"]) + "\n")
    return sha


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--registration", type=Path, default=ROOT / "registration/registration_v3.json")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--external", type=Path, default=EXT)
    p.add_argument("--external-embedding", type=Path, default=EXT_EMB)
    p.add_argument("--print", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    import torch

    torch.set_num_threads(4)
    a = parse_args(argv)
    doc = build(a)
    if a.print:
        print(json.dumps(doc, indent=1))
        return 0
    sha = write(doc)
    log(f"wrote {ADDENDUM} sha256 {sha} (HASHES.txt updated)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
