#!/usr/bin/env python3
"""Build amendment A3 to the v3 registration (registration/amendment_A3.json) from training rows only.

A3 follows A1 and A2 and is written before any site4 evaluation of E3. registration_v3.json, A1, A2 and
the committed addenda are not changed; A3 overrides one experiment field in ``declared()``:

* A3.1 (E3.H3a_pooling): the pooling claim (D_R2 - D_R1 >= 0.05, "the loss is in mean pooling") gets the
  scale guard A2.1 gave E3.H3a: a registered pooling win also needs (b') GR_R2(flagged) - GR_R1(flagged) > 0
  and (c') Dlog_R2 - Dlog_R1 > 0, each with its two-stage 95% lower bound > 0 and > 0 in each primary donor,
  otherwise the verdict is 'not supported (scale artefact)'. Same algebra as A2.1 (R1 in the head's place);
  stricter only. The two statistics are the ones commit b4ecad4 added to the E3 evaluate code as a secondary
  that never decides; under A3 they decide.

The rule is fixed by the algebra of the statistic; the only numbers computed here are a point-estimate
illustration on NK-T look-alike pairs of the training donors (pairs as E3 builds them, the registered
phase-1 head's evidence): the head's flagged / unflagged gap-ratio difference, which sets how much uniform
shrinkage reaches the registered pooling margin, and D_R2 - D_R1, (b') and (c') for a pure rescaling and for
a flagged-pair repair of a head-like R1. Nothing computed here decides anything or sets a number of A3.
Protein, cell type and embedding values are read on split=train rows only. ``--leakage-check`` rebuilds the
core from copies whose site4 rows (protein, cell type, embedding) are random: the core must be
byte-identical; positive control: poisoning the training rows must change it.

Usage (from the repo root):
  python bridge_anm/v3_build_amendment_A3.py                 # core + leakage check + amendment_A3.json(.sha256) + .md
  python bridge_anm/v3_build_amendment_A3.py --core-only --out-dir D
  python bridge_anm/v3_build_amendment_A3.py --record-hashes # A3 lines in registration/HASHES.txt
"""
from __future__ import annotations

import os
import sys

_THREADS = "4"
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, _THREADS)

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va_  # noqa: E402
from lib import v3_e3 as e3  # noqa: E402
from lib import v3_key as vk  # noqa: E402

BUILDER_VERSION = "v3_build_amendment_A3 1.0"
POOL_MARGIN = 0.05             # registered E3.H3a_pooling margin (unchanged)
SHRINK_S = 0.5                 # illustration only: R2 keeps half of R1's NK-T evidence differences (declared, not fitted)
POISON_SEEDS = {"test": 101, "train": 202, "val": 303}   # as A2
STATISTICS_COMMIT = "b4ecad4"  # E3 evaluate code that already computes (b') and (c') (decides: False before A3)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def r6(x):
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return float(round(x, 6)) if np.isfinite(x) else None
    if isinstance(x, dict):
        return {k: r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r6(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def content_bytes(d: dict) -> bytes:
    return (json.dumps(d, indent=1, sort_keys=False, ensure_ascii=False) + "\n").encode()


# ============================================================================ inputs
class Inputs:
    """The arrays the core reads. ``split`` and ``donors`` are read on every row; adt, cell types and z on
    split=train rows only."""

    def __init__(self, split, donors, adt, adt_names, cell_types, z):
        self.split, self.donors, self.adt, self.adt_names = split, donors, adt, list(adt_names)
        self.cell_types, self.z = cell_types, z

    @classmethod
    def load(cls, processed: Path, z_path: Path) -> "Inputs":
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=False)
        z = np.load(z_path, mmap_mode="r")
        return cls(npz["split"].astype(str), npz["donors"].astype(str), np.asarray(npz["adt"], dtype=np.float32),
                   [str(x) for x in npz["adt_names"]], npz["cell_types"].astype(str), np.asarray(z, dtype=np.float32))

    def poisoned(self, which: str, seed: int) -> "Inputs":
        """Copy whose rows of one split carry random protein, cell types and embedding (as A2)."""
        rng = np.random.default_rng(seed)
        m = self.split == which
        adt, ct, z = self.adt.copy(), self.cell_types.copy(), self.z.copy()
        adt[m] = rng.gamma(1.0, 1.0, size=(int(m.sum()), adt.shape[1])).astype(adt.dtype)
        ct[m] = rng.choice(np.unique(self.cell_types), size=int(m.sum()))
        z[m] = rng.normal(size=(int(m.sum()), z.shape[1])).astype(np.float32)
        return Inputs(self.split, self.donors, adt, self.adt_names, ct, z)


# ============================================================================ statistics (as the E3 evaluate code)
def gr_from(dE_r: np.ndarray, dM: np.ndarray, sel: np.ndarray) -> float:
    """scripts/v3_e3_nkt_repair.gr_from: mean over the gap proteins of median |evidence diff| / median |measured diff|."""
    if sel.size == 0:
        return float("nan")
    num = np.median(dE_r[sel], axis=0)
    den = np.median(dM[sel], axis=0)
    r = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return float(np.mean(r)) if np.all(np.isfinite(r)) else float("nan")


def pooling_stats(dR1: np.ndarray, dR2: np.ndarray, dM: np.ndarray, flag: np.ndarray, ix: np.ndarray) -> dict:
    """The E3 evaluate quantities of E3.H3a_pooling on pair indices ix: D_R2 - D_R1, (b') and (c'), with the head
    and the null set to R1 (they cancel in D_R2 - D_R1 and in Dlog_R2 - Dlog_R1)."""
    gr = {}
    for r, d in (("R1", dR1), ("R2", dR2)):
        gr[r] = {nm: gr_from(d, dM, ix[fm]) for nm, fm in (("flagged", flag[ix]), ("unflagged", ~flag[ix]))}
    gr["head"] = gr["null"] = gr["R1"]
    return {"GR_R1_flagged": gr["R1"]["flagged"], "GR_R1_unflagged": gr["R1"]["unflagged"],
            "GR_R2_flagged": gr["R2"]["flagged"], "GR_R2_unflagged": gr["R2"]["unflagged"],
            "D_R2_minus_D_R1": e3.d_statistic(gr, "R2") - e3.d_statistic(gr, "R1"),
            "term_flagged_R2_minus_R1": gr["R2"]["flagged"] - gr["R1"]["flagged"],
            "Dlog_R2_minus_Dlog_R1": e3.d_log_statistic(gr, "R2") - e3.d_log_statistic(gr, "R1")}


def scenario_points(dR1, dR2, dM, flag, donor, order) -> dict:
    """Point values (pooled and per training donor) of D_R2 - D_R1, (b') and (c') for one constructed readout pair."""
    allix = np.arange(flag.size)
    keys = ("D_R2_minus_D_R1", "term_flagged_R2_minus_R1", "Dlog_R2_minus_Dlog_R1")
    pt = pooling_stats(dR1, dR2, dM, flag, allix)
    pdn = {d: pooling_stats(dR1, dR2, dM, flag, allix[donor == d]) for d in order}
    return {k: {"point": pt[k], "per_donor": {d: pdn[d][k] for d in order}} for k in keys}


# ============================================================================ computed core
def build_core(X: Inputs, reg: dict, e3add: dict, head_ckpt: Path) -> dict:
    t0 = time.time()
    split = X.split
    n = int(split.size)
    tr = np.where(split == "train")[0].astype(np.int64)   # the only rows whose protein / cell type / embedding are read
    spec = reg["splits"]["train"]
    if tr.size != int(spec["n"]) or vk.index_hash(tr) != spec["sha256_indices"]:
        raise SystemExit("split=train rows differ from the registration")
    flag_cos = float(reg["e3"]["flag_cosine"])
    targets = list(e3add["declared"]["targets"])
    j = {nm: i for i, nm in enumerate(X.adt_names)}
    tc = [j[p] for p in targets]
    gidx = [targets.index(p) for p in e3.GAP_PROTEINS]
    # registered primary key of the training cells (verifier), as E3 keys_for
    k = vk.build_keys(X.cell_types[tr], X.adt[tr], X.adt_names, reg)
    key = np.full(n, vk.UNSCORED, dtype=object)
    key[tr] = k["primary"]
    key = key.astype(str)
    P = e3.nkt_pairs_within_donors(X.z, tr, X.donors, {"all": key}, 10, flag_cos)
    pv = P["variants"]["all"]
    order = sorted(set(X.donors[tr].tolist()))
    log(f"core: {pv['nk'].size} training NK-T pairs ({int(pv['flag'].sum())} flagged) ({time.time() - t0:.0f}s)")
    cells = np.unique(np.concatenate([pv["nk"], pv["t"]])) if pv["nk"].size else np.zeros(0, np.int64)
    dH = dM = np.zeros((0, len(gidx)))
    if cells.size:
        qh = np.asarray(reg["evidence"]["teddy_head"]["q95_train_pred"], np.float64)
        ph = vk.head_predict(np.asarray(X.z[cells], dtype=np.float32), head_ckpt, device="cpu", size_factor=1.0)
        ev = e3.normalise(ph[:, tc], qh[tc])                               # head evidence (registered normaliser)
        q95m = np.asarray([e3add["computed_train_val"]["q95_measured_train"][p] for p in targets])
        meas = e3.normalise(X.adt[cells][:, tc], q95m)                      # measured evidence (addendum q95)
        a_, b_ = np.searchsorted(cells, pv["nk"]), np.searchsorted(cells, pv["t"])
        dH = np.abs(ev[a_][:, gidx] - ev[b_][:, gidx])
        dM = np.abs(meas[a_][:, gidx] - meas[b_][:, gidx])
    flag, donor = pv["flag"], pv["donor"]
    allix = np.arange(flag.size)
    gh = {nm: gr_from(dH, dM, allix[fm]) for nm, fm in (("flagged", flag), ("unflagged", ~flag))}
    gh_d = {d: {nm: gr_from(dH, dM, allix[(donor == d) & fm]) for nm, fm in (("flagged", flag), ("unflagged", ~flag))}
            for d in order}
    delta = gh["unflagged"] - gh["flagged"]
    scen = {"uniform_shrinkage": scenario_points(dH, SHRINK_S * dH, dM, flag, donor, order),
            "flagged_repair": scenario_points(dH, np.where(flag[:, None], dM, dH), dM, flag, donor, order)}
    us = scen["uniform_shrinkage"]
    ident = (abs(us["D_R2_minus_D_R1"]["point"] - (1 - SHRINK_S) * delta) if np.isfinite(delta) else None)
    log(f"core: illustration done ({time.time() - t0:.0f}s)")
    core = {
        "amendment_id": va_.A3_ID,
        "builder": BUILDER_VERSION,
        "train_pairs": {
            "rule": "as E3 evaluate (addendum declared.pairs): k = 10 cosine neighbours of the raw official z within each "
                    "training donor (v3_e3.nkt_pairs_within_donors), NK-T pair = one key-NK and one key-T cell (registration "
                    f"primary key) of the same donor, flagged if cosine >= {flag_cos}",
            "donors": order, "n_pairs": int(flag.size), "n_flagged": int(flag.sum()),
            "per_donor": {d: {"n_pairs": int(np.sum(donor == d)), "n_flagged": int(np.sum(flag & (donor == d)))} for d in order},
            "pairs_sha256": e3.pairs_sha256(pv["nk"], pv["t"]),
        },
        "head_gap_ratio": {
            "what": "the registered phase-1 head's gap ratio GR (E3 addendum declared.gap_ratio; evidence = prediction / "
                    "registered q95_train_pred, clipped; measured = ADT / addendum q95_measured_train, clipped) on the "
                    "training pairs",
            "GR_head_flagged": gh["flagged"], "GR_head_unflagged": gh["unflagged"],
            "unflagged_minus_flagged": delta,
            "per_donor": {d: {"GR_head_flagged": v["flagged"], "GR_head_unflagged": v["unflagged"],
                              "unflagged_minus_flagged": v["unflagged"] - v["flagged"]} for d, v in gh_d.items()},
        },
        "scale_artefact_size": {
            "statement": "if R2 keeps a common fraction s of R1's NK-T evidence differences, D_R2 - D_R1 = (1 - s) "
                         "(GR_R1(unflagged) - GR_R1(flagged)); with a readout like the head as R1 the registered pooling "
                         "margin 0.05 is reached by a pure rescaling with 1 - s = 0.05 / (GR_head(unflagged) - "
                         "GR_head(flagged)) when that difference is positive",
            "one_minus_s_reaching_margin": (POOL_MARGIN / delta) if (np.isfinite(delta) and delta > 0) else None,
            "registered_margin": POOL_MARGIN,
        },
        "illustration": {
            "what": "point values (pooled over the training pairs and per training donor; no bootstrap, no verdict) of "
                    "D_R2 - D_R1, (b') term_flagged_R2_minus_R1 and (c') Dlog_R2_minus_Dlog_R1, computed as E3 evaluate "
                    "does, for two constructed readout pairs with R1 = the head: (uniform_shrinkage) R2 = the head's NK-T "
                    "evidence differences times s, which repairs nothing; (flagged_repair) R2 = the measured differences "
                    "on flagged pairs and the head's on unflagged pairs; head and null set to R1 (they cancel)",
            "s": SHRINK_S, "donor_order": order,
            "scenarios": scen,
            "uniform_shrinkage_identity_abs_error": ident,
        },
        "rows_read": {"split_and_donor_columns": n, "protein_cell_type_embedding_rows": int(tr.size),
                      "val_rows_values_read": 0, "site4_rows_values_read": 0,
                      "n_val_rows": int(np.sum(split == "val")), "n_site4_rows": int(np.sum(split == "test"))},
    }
    return r6(core)


# ============================================================================ leakage check
def leakage_check(X: Inputs, reg: dict, e3add: dict, head_ckpt: Path, real_bytes: bytes) -> dict:
    out = {"pattern": "as A2 (v3_build_amendment_A2.py) and A1 (v3_leakage_check.py): rebuild the core from copies of the "
                      "inputs whose rows of one split carry random protein, cell types and embedding, and compare bytes",
           "real_core_sha256": sha256_bytes(real_bytes)}
    for which, role in (("test", "site4 rows poisoned: core must be byte-identical"),
                        ("val", "val rows poisoned: A3 reads no val value, so identical as well (information)"),
                        ("train", "positive control, training rows poisoned: core must change")):
        b = content_bytes(build_core(X.poisoned(which, POISON_SEEDS[which]), reg, e3add, head_ckpt))
        out[f"poison_{which}"] = {"seed": POISON_SEEDS[which], "role": role, "core_sha256": sha256_bytes(b),
                                  "identical_to_real": b == real_bytes}
        log(f"leakage: {which} rows poisoned -> identical {b == real_bytes}")
    out["passed"] = bool(out["poison_test"]["identical_to_real"] and not out["poison_train"]["identical_to_real"])
    return out


# ============================================================================ declared content
REG_POOLING = "D_R2 - D_R1 >= 0.05: the loss is in mean pooling"
COND_B = ("(b') recovery on flagged pairs beyond R1: term_flagged_R2_minus_R1 = GR_R2(flagged) - GR_R1(flagged) > 0, its "
          "two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor")
COND_C = ("(c') the scale-free contrast Dlog_R2 - Dlog_R1 = [ln GR_R2(fl) - ln GR_R1(fl)] - [ln GR_R2(unfl) - ln GR_R1(unfl)] "
          "(Dlog_R as the E3 builder computes it, addendum secondary_H3a S1), > 0 with its two-stage 95% interval "
          "excluding 0 (lower bound > 0) and > 0 in each primary donor")


def declared(core: dict, reg: dict, a2: dict, e3add_sha: str, a2_sha: str) -> dict:
    sa = core["scale_artefact_size"]
    hg = core["head_gap_ratio"]
    il = core["illustration"]["scenarios"]
    tp = core["train_pairs"]
    findings = [
        {"id": "A3.1", "area": "E3.H3a_pooling (D_R2 - D_R1, 'the loss is in mean pooling'): scale artefact of the registered contrast",
         "found_by": "review of E3 before any site4 evaluation of E3: A2 left E3.H3a_pooling on the registered rule (A2.1 "
                     "'H3a_pooling (D_R2 - D_R1) is unchanged'), and commit b4ecad4 then added the A2.1-style conditions "
                     "(b') and (c') beside the registered pooling verdict as a secondary that never decides "
                     "(scripts/v3_e3_nkt_repair.py, run_h3a / pooling_scale_check, decides: False)",
         "finding": "the head and the null cancel in the registered pooling contrast: D_R2 - D_R1 = [GR_R2(fl) - GR_R1(fl)] "
                    "- [GR_R2(unfl) - GR_R1(unfl)], which is not scale-free: if R2 keeps a common fraction s < 1 of R1's "
                    "NK-T evidence differences (GR_R2 = s GR_R1 on every pair set), D_R2 - D_R1 = (1 - s) * (GR_R1(unfl) - "
                    "GR_R1(fl)) > 0 whenever R1 keeps less of the gap on flagged pairs. Such an R2 recovers less of the "
                    "measured gap than R1 on flagged pairs, yet can pass the registered win rule and support 'the loss is "
                    "in mean pooling'. This is the algebra of A2.1 with R1 in the head's place",
         "fix": "the registered D_R2 - D_R1 computation, margin 0.05 and win / loss / equivalent rules are kept. The pooling "
                "claim wins only if the registered rule gives a win AND, under the same registered two-stage donor "
                "bootstrap (B = 2000, seed 1, the same replicates as D_R2 - D_R1) and per-donor rule: "
                f"{COND_B}; {COND_C}. An undefined value fails. If the registered pooling verdict is a win but (b') or "
                "(c') fails, the E3.H3a_pooling verdict is 'not supported (scale artefact)'; a loss, equivalent or "
                "inconclusive registered verdict is unchanged. (b') and (c') are the statistics term_flagged_R2_minus_R1 "
                f"and Dlog_R2_minus_Dlog_R1 that the E3 evaluate code computes since commit {STATISTICS_COMMIT} with "
                "decides = False; under A3 they decide. E3.H3a (as amended by A2.1), E3.H3b and the E3 falsification "
                "(A2.1) are unchanged",
         "effect": "stricter: every A3 pooling win is a registered pooling win (v3_amend.a3_pooling_verdict). Training-pair "
                   f"illustration (computed; {tp['n_pairs']} NK-T pairs of the {len(tp['donors'])} training donors, "
                   f"{tp['n_flagged']} flagged): the head's GR is {hg['GR_head_flagged']} on flagged and "
                   f"{hg['GR_head_unflagged']} on unflagged pairs, so with a head-like R1 a pure rescaling reaches the "
                   f"margin at 1 - s = {sa['one_minus_s_reaching_margin']} (point estimate). R2 = {core['illustration']['s']} "
                   f"x R1 has D_R2 - D_R1 = {il['uniform_shrinkage']['D_R2_minus_D_R1']['point']} with (b') "
                   f"{il['uniform_shrinkage']['term_flagged_R2_minus_R1']['point']} and (c') "
                   f"{il['uniform_shrinkage']['Dlog_R2_minus_Dlog_R1']['point']}: positive D, nothing repaired; a "
                   f"flagged-pair repair has (b') {il['flagged_repair']['term_flagged_R2_minus_R1']['point']} and (c') "
                   f"{il['flagged_repair']['Dlog_R2_minus_Dlog_R1']['point']}. The head's flagged / unflagged difference "
                   "varies in sign across training donors (computed.head_gap_ratio.per_donor), so the size of the "
                   "artefact depends on R1; the rule does not"},
    ]
    overrides = {
        "E3/endpoints/E3.H3a_pooling":
            REG_POOLING + f" (registered; margin 0.05 with the registration win / loss / equivalent rules). A3.1: the pooling "
            f"claim wins only if the registered rule wins and, under the same two-stage donor bootstrap and per-donor rule, "
            f"{COND_B} and {COND_C}; a registered pooling win without (b') and (c') is 'not supported (scale artefact)'",
    }
    e3 = {"h3a_pooling_win": {
        "registered_D": "unchanged: D_R2 - D_R1 as registered, margin 0.05, registration win / loss / equivalent rules "
                        "(E3 addendum declared.H3a: 'pooling = D_R2 - D_R1')",
        "conditions": {"b_prime": {"statistic": "term_flagged_R2_minus_R1 = GR_R2(flagged) - GR_R1(flagged)",
                                   "rule": "point > 0, two-stage 95% lower bound > 0, > 0 in each primary donor"},
                       "c_prime": {"statistic": "Dlog_R2_minus_Dlog_R1 = Dlog_R2 - Dlog_R1 (E3 addendum secondary_H3a S1)",
                                   "rule": "point > 0, two-stage 95% lower bound > 0, > 0 in each primary donor"}},
        "bootstrap": "the registered two-stage bootstrap of E3.H3a (A1.1): donors then pairs within donors, B = 2000, seed 1, "
                     "the same replicates as D_R2 - D_R1; percentile interval over the defined replicates",
        "verdict": "win if the registered pooling verdict is win and (b') and (c') hold; 'not supported (scale artefact)' if "
                   "it is a win and (b') or (c') fails; otherwise the registered pooling verdict",
        "decides": "test_primary, variant all; other splits, variants and sensitivities are reported with the same rule",
        "statistics_in_code": f"scripts/v3_e3_nkt_repair.py h3a_stats: term_flagged_R2_minus_R1, Dlog_R2_minus_Dlog_R1 "
                              f"(added by commit {STATISTICS_COMMIT} as a secondary with decides = False)",
        "functions": "v3_amend.a2_positive_condition, v3_amend.a3_pooling_verdict",
    }}
    disclosure = {
        "E1_site4_results_seen": "before A3 was written, the orchestrating session had already seen the E1 site4 (test) "
                                 "results (E1 is the Mode B rerun). A3 concerns only E3 (the E3.H3a_pooling verdict) and "
                                 "follows from the algebra of the registered statistic (A2.1's algebra with R1 in the "
                                 "head's place); it reads no E1 output and nothing in it depends on an E1 result",
        "E3_site4_state": "no E3 site4 evaluation had been run: E3 evaluate, the only E3 stage that reads site4 protein or "
                          "cell types, needs R2's site4 predictions, which E3 r2_predict (a TEDDY + R2 forward of site4 RNA; "
                          "no protein or label) was still writing while A3 was written; A3 and its builder read none of them. "
                          "From A3 on, E3 evaluate on site4 is refused unless A3 is committed and matches",
        "E3_val_smoke_seen": "the val-only smoke runs of E3 evaluate (val donor 18303, small subsets, labelled smoke) had "
                             "been seen; A3's rule does not depend on them",
        "numbers_computed_for_A3": "the only numbers computed for A3 use split=train rows (protein, cell types, embedding) and the "
                         "split and donor columns (leakage_check)",
    }
    supersedes_add = [
        {"file": "registration/addenda/E3.json", "sha256": e3add_sha, "unchanged_file": True,
         "fields": ["declared.H3a: 'pooling = D_R2 - D_R1; margin 0.05 with the registration win / loss / equivalent rules' "
                    "-> the A3.1 pooling win (registered win plus (b') and (c'))",
                    "declared.secondary_H3a: S1 (Dlog_R) enters the A3.1 pooling condition (c') as Dlog_R2 - Dlog_R1"]},
    ]
    supersedes_am = [
        {"file": "registration/amendment_A2.json", "sha256": a2_sha, "unchanged_file": True,
         "fields": ["findings A2.1 fix: 'H3a_pooling (D_R2 - D_R1) is unchanged' -> A3.1"]},
    ]
    change_log = [
        {"file": "registration/amendment_A3.json (+ .sha256, amendment_A3_core.json, AMENDMENT_A3.md)",
         "change": "new amendment A3", "reason": "A3.1 above; written before any site4 evaluation of E3"},
        {"file": "bridge_anm/lib/v3_amend.py",
         "change": "load_registration_amended() applies A1, A2, then A3, verifying the registration, A1, A2 and A3 hashes "
                   "and that A3 names the registration, A1 and A2 on disk; amendment_A3_status() for the E3 site4 evaluate "
                   "guard; a3_pooling_verdict()",
         "reason": "every v3 builder reads the amended registration through this loader"},
        {"file": "scripts/v3_e3_nkt_repair.py",
         "change": "evaluate on site4 refuses unless A3 is committed and matches; the E3.H3a_pooling verdict is the A3.1 "
                   f"verdict, computed on unrounded values and the same replicates (the {STATISTICS_COMMIT} secondary block "
                   "becomes the deciding A3 block; the registered pooling verdict is reported beside it); DECLARED (the "
                   "addendum's declared part) is unchanged",
         "reason": "A3.1; earlier stages (null, r1, r2_states, r2_fit, r2_predict) are unchanged"},
        {"file": "tests/test_v3_amend.py, tests/test_e2_response.py",
         "change": "A3 loading, tamper / chain / missing-file detection, status, the E3 evaluate guard, A3.1 end to end "
                   "through run_h3a / verdicts; the E2 guard test copies A3 (the loader needs it)",
         "reason": "A3"},
        {"file": "registration/HASHES.txt", "change": "A3 lines; v3_amend.py line moved to the A3 version (the A2 version's "
                                                       "sha256 is in amendment_A2.json provenance)",
         "reason": "hash record of the registration files"},
    ]
    return {"title": "Amendment A3 to the v3 pre-registration: E3 H3a_pooling scale artefact",
            "written_before": "any site4 (test) evaluation of E3; the only numbers computed here use split=train rows "
                              "(protein, cell types, embedding) and the split and donor columns (leakage_check)",
            "scope": "registration_v3.json, amendment_A1.json, amendment_A2.json and the addenda are unchanged; builders "
                     "apply A1, A2, then A3 (v3_amend.load_registration_amended()). Where A3 differs from the registration "
                     "as amended by A1 and A2, from A2 (supersedes_in_amendments) or from the E3 addendum "
                     "(supersedes_in_addenda), A3 wins. A3 only makes one verdict stricter (A3.1) and concerns E3 only",
            "disclosure": disclosure,
            "findings": findings, "experiments_overrides": overrides, "e3": e3,
            "supersedes_in_addenda": supersedes_add, "supersedes_in_amendments": supersedes_am, "change_log": change_log}


# ============================================================================ render
def render_md(full: dict) -> str:
    cap = lambda t: t[:1].upper() + t[1:]  # noqa: E731
    am = full["amends"]
    L = [f"# {full['title']}", "",
         f"**Status:** written before {full['written_before']}. File `registration/amendment_A3.json`, sha256 in "
         "`amendment_A3.json.sha256`. It follows amendment A2 (`registration/amendment_A2.json` sha256 "
         f"`{am['amendment_A2_sha256']}`) and A1 (`registration/amendment_A1.json` sha256 `{am['amendment_A1_sha256']}`) "
         f"and amends `registration/registration_v3.json` sha256 `{am['registration_sha256']}`; all three are unchanged. "
         "This page is rendered from the JSON by `bridge_anm/v3_build_amendment_A3.py`.", "",
         f"**Scope.** {full['scope']}.", "", "## Disclosure", ""]
    for k, v in full["disclosure"].items():
        L.append(f"- **{k.replace('_', ' ')}.** {cap(v)}.")
    L += ["", "## Findings and fixes", ""]
    for f in full["findings"]:
        L += [f"### {f['id']} {f['area']}", "", f"- **Found by.** {cap(f['found_by'])}.",
              f"- **Finding.** {cap(f['finding'])}.", f"- **Fix.** {cap(f['fix'])}.", f"- **Effect.** {cap(f['effect'])}.", ""]
    c = full["computed"]
    tp, hg, sa, il = c["train_pairs"], c["head_gap_ratio"], c["scale_artefact_size"], c["illustration"]
    L += ["## Illustration on training pairs (training rows only; point values; decides nothing)", "",
          f"Pairs: {cap(tp['rule'])}. {tp['n_pairs']} NK-T pairs ({tp['n_flagged']} flagged) in the training donors "
          + ", ".join(f"{d} {x['n_pairs']} ({x['n_flagged']})" for d, x in tp["per_donor"].items())
          + f"; pairs sha256 `{tp['pairs_sha256'][:16]}...`.", "",
          f"{cap(hg['what'])}: flagged {hg['GR_head_flagged']}, unflagged {hg['GR_head_unflagged']}, difference "
          f"{hg['unflagged_minus_flagged']}. {cap(sa['statement'])}: 1 - s = {sa['one_minus_s_reaching_margin']}.", "",
          "| training donor | GR_head flagged | GR_head unflagged | unflagged - flagged |", "|---|---:|---:|---:|"]
    for d, v in hg["per_donor"].items():
        L.append(f"| {d} | {v['GR_head_flagged']} | {v['GR_head_unflagged']} | {v['unflagged_minus_flagged']} |")
    L += ["", f"{cap(il['what'])}. s = {il['s']}.", "",
          "| scenario | statistic | pooled | " + " | ".join(il["donor_order"]) + " |",
          "|---|---|---:|" + "---:|" * len(il["donor_order"])]
    for nm, sc in il["scenarios"].items():
        for k, row in sc.items():
            L.append(f"| {nm} | {k} | {row['point']} | " + " | ".join(str(row["per_donor"][d]) for d in il["donor_order"]) + " |")
    L += ["", f"Uniform shrinkage: |D_R2 - D_R1 - (1 - s)(GR_head(unfl) - GR_head(fl))| = "
          f"{il['uniform_shrinkage_identity_abs_error']} (the identity above).", ""]
    lk = full["leakage_check"]
    L += ["## Leakage check", "", f"{cap(lk['pattern'])}.", "", "| poisoned rows | role | core identical to real |", "|---|---|---|"]
    for k in ("poison_test", "poison_val", "poison_train"):
        L.append(f"| {k.split('_')[1]} (seed {lk[k]['seed']}) | {lk[k]['role']} | {lk[k]['identical_to_real']} |")
    L += ["", f"Result: **{'PASS' if lk['passed'] else 'FAIL'}**; real core sha256 `{lk['real_core_sha256']}`.", "",
          "## Experiment fields replaced (path in registration_v3.json → experiments, after A1 and A2)", ""]
    for p_, v in full["experiments_overrides"].items():
        L.append(f"- `{p_}`: {v}")
    L += ["", "## Addendum and amendment fields superseded (files unchanged)", ""]
    for sp in full["supersedes_in_addenda"] + full["supersedes_in_amendments"]:
        L.append(f"- `{sp['file']}` (sha256 `{sp['sha256'][:16]}...`): " + "; ".join(sp["fields"]))
    L += ["", "## Change log", "", "| file | change | reason |", "|---|---|---|"]
    for ch in full["change_log"]:
        L.append(f"| `{ch['file']}` | {ch['change']} | {ch['reason']} |")
    pv = full["provenance"]
    L += ["", "## Provenance", "",
          f"Builder `{pv['builder']}` sha256 `{pv['builder_sha256']}`; `v3_amend.py` sha256 `{pv['v3_amend_sha256']}`; "
          f"`v3_e3.py` sha256 `{pv['v3_e3_sha256']}`; computed core sha256 `{pv['core_sha256']}`; built {pv['built_utc']}.", ""]
    return "\n".join(L)


# ============================================================================ hashes
A3_HASH_FILES = ("amendment_A3.json", "amendment_A3_core.json", "AMENDMENT_A3.md",
                 "../bridge_anm/lib/v3_amend.py", "../bridge_anm/v3_build_amendment_A3.py")


def record_hashes(reg_dir: Path) -> None:
    hf = reg_dir / "HASHES.txt"
    lines = hf.read_text().splitlines()
    names = set(A3_HASH_FILES)
    keep = [ln for ln in lines if not (len(ln.split()) == 2 and ln.split()[1] in names) and not ln.startswith("# A3:")]
    if keep and keep[0].startswith("#"):
        keep[0] = ("# sha256 of the v3 registration files (recorded at commit; registration_v3.json is the frozen "
                   "registration, amendment_A1.json amends it before site4, amendment_A2.json follows A1 before any site4 "
                   "evaluation of E2 and E3, amendment_A3.json follows A2 before any site4 evaluation of E3)")
    keep.append("# A3: v3_amend.py extended by A3 (its A2-era sha256 59089ef5... is recorded in amendment_A2.json provenance)")
    keep += [f"{vk.sha256_file(reg_dir / n)}  {n}" for n in A3_HASH_FILES]
    hf.write_text("\n".join(keep) + "\n")
    log(f"recorded {len(A3_HASH_FILES)} A3 lines in {hf}")


# ============================================================================ main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--out-dir", type=Path, default=None, help="default: --registration-dir")
    p.add_argument("--core-only", action="store_true", help="write only amendment_A3_core.json")
    p.add_argument("--skip-leakage-check", action="store_true", help="development only; the amendment needs the check")
    p.add_argument("--record-hashes", action="store_true", help="only rewrite the A3 lines of registration/HASHES.txt")
    args = p.parse_args(argv)
    rd = args.registration_dir
    if args.record_hashes:
        record_hashes(rd)
        return
    out_dir = args.out_dir or rd
    out_dir.mkdir(parents=True, exist_ok=True)
    import torch
    torch.set_num_threads(int(_THREADS))
    reg_f, a1_f, a2_f = rd / "registration_v3.json", rd / "amendment_A1.json", rd / "amendment_A2.json"
    reg = vk.load_registration(reg_f)
    va_.load_amendment(a1_f, registration_path=reg_f)
    a2 = va_.load_amendment_A2(a2_f, registration_path=reg_f, amendment_A1_path=a1_f)
    e3f = rd / "addenda" / "E3.json"
    if f"{vk.sha256_file(e3f)}  {e3f.name}" not in (rd / "addenda" / "HASHES.txt").read_text().splitlines():
        raise SystemExit(f"{e3f} does not match addenda/HASHES.txt")
    e3add = json.loads(e3f.read_text())
    if vk.sha256_file(args.ckpt) != reg["provenance"]["inputs"]["ckpt_sha256"]:
        raise SystemExit("head checkpoint differs from the registered input")
    X = Inputs.load(args.processed, args.z)
    core = build_core(X, reg, e3add, args.ckpt)
    core_b = content_bytes(core)
    (out_dir / "amendment_A3_core.json").write_bytes(core_b)
    log(f"core sha256 {sha256_bytes(core_b)}")
    if args.core_only:
        return
    leak = (leakage_check(X, reg, e3add, args.ckpt, core_b) if not args.skip_leakage_check
            else {"passed": False, "skipped": True})
    if not leak["passed"] and not args.skip_leakage_check:
        raise SystemExit(f"leakage check failed: {leak}")
    full = {"amendment_id": va_.A3_ID,
            "amends": {"file": "registration/registration_v3.json", "registration_sha256": vk.sha256_file(reg_f),
                       "after": "registration/amendment_A2.json", "amendment_A1_sha256": vk.sha256_file(a1_f),
                       "amendment_A2_sha256": vk.sha256_file(a2_f)},
            **declared(core, reg, a2, vk.sha256_file(e3f), vk.sha256_file(a2_f)),
            "computed": core, "leakage_check": leak,
            "provenance": {"builder": BUILDER_VERSION, "builder_sha256": vk.sha256_file(Path(__file__)),
                           "v3_amend_sha256": vk.sha256_file(Path(va_.__file__)),
                           "v3_e3_sha256": vk.sha256_file(Path(e3.__file__)),
                           "v3_key_sha256": vk.sha256_file(Path(vk.__file__)),
                           "core_sha256": sha256_bytes(core_b),
                           "inputs": {"cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz"),
                                      "z_sha256": vk.sha256_file(args.z), "ckpt_sha256": vk.sha256_file(args.ckpt),
                                      "addendum_E3_sha256": vk.sha256_file(e3f)},
                           "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                           "test_rows_read": "the split and donor columns only (leakage_check)"}}
    out = out_dir / "amendment_A3.json"
    out.write_bytes(content_bytes(full))
    h = vk.sha256_file(out)
    (out_dir / "amendment_A3.json.sha256").write_text(f"{h}  amendment_A3.json\n")
    (out_dir / "AMENDMENT_A3.md").write_text(render_md(full))
    log(f"wrote {out} sha256 {h}")


if __name__ == "__main__":
    main()
