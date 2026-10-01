#!/usr/bin/env python3
"""E1 addendum (registration/addenda/E1.json): the numbers E1 needs that the registration leaves open,
fixed by the procedures the registration writes, from training and validation cells only.

What it fixes (registration section 9 'Addenda'; A1.10: an addendum may only fix open numbers):

* sensitivity panels (canonical9_plus_nk, canonical9_variant_plus_nk, val_k2, val_k5), which the
  registration lists for 'rule and classifier' rows without bars or classifier settings:
  - rule bars for Q1 and Q2 = the val quantile of the top panel score at the declared val no-call
    rates 0.15 / 0.30 (section 7, the procedure that set the primary-panel bars);
  - the classifier's C = lowest val log-loss on the registered grid (section 8), and its Q1 / Q2 bars
    = the val quantile of the highest lineage probability at the same no-call rates (section 8 bar_rule);
* Q3 operating-point bars (val no-call 0.15, section 8 bar_rule) for the four A1.6 classifiers and the
  registered Q3 classifier, used only to report each method's deployed operating point (every E1
  comparison is at matched coverage and needs no bar).

It also recomputes the registered numbers E1 reads (rule bars, closure bars, classifier C, val
log-loss and bars, the A1.6 C values) and stops if any differs from the registration: E1 then runs
exactly the registered methods.

The implementation choices of the E1 builder that the registration does not fix are declared here
as text (``IMPLEMENTATION_CHOICES``), before any site4 evaluation.

As in v3_build_registration.py, the protein matrix, the cell types and the embedding are read only
on non-test rows (``nt``). ``--leakcheck`` rebuilds the computed core from inputs whose site4 rows
are replaced by random values (must be byte-identical) and from val-poisoned inputs (must change),
with the poisoning helpers of bridge_anm/v3_leakage_check.py.

Usage (repo root):
  ANM_ROOT=... python bridge_anm/v3_e1_addendum.py               # writes registration/addenda/E1.json + HASHES.txt line
  python bridge_anm/v3_e1_addendum.py --core-only --out-dir D      # computed core only
  python bridge_anm/v3_e1_addendum.py --leakcheck                  # static + poisoning checks
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
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va_  # noqa: E402
from lib import v3_e1 as e1  # noqa: E402
from lib import v3_key as vk  # noqa: E402

BUILDER_VERSION = "v3_e1_addendum 1.0"
CORE_NAME = "E1_addendum_core.json"
ADDENDUM = ROOT / "registration" / "addenda" / "E1.json"
HASHES = ROOT / "registration" / "addenda" / "HASHES.txt"
NO_CALL_TARGETS = {"Q1": 0.15, "Q2": 0.30, "Q3": 0.15}

IMPLEMENTATION_CHOICES = [
    "Classifier decision: the most probable of the four lineages, confidence = its probability (registration "
    "classifier.primary.decision); OUT is never called. A label-cost draw with one class predicts that class.",
    "E1.4: the top score (rule = ANM), the margin and the entropy score rank the rule's argmax calls; a classifier "
    "ranks its own argmax calls. The classifier probability ranking the rule's calls is a labelled diagnostic.",
    "Entropy confidence 1 - H(q)/log 4 with q = S / sum(S); when sum(S) = 0 it is 0.",
    "Argmax ties between classes go to the first class in the order B, T, NK, myeloid (v3_key.rule_calls; ANM's "
    "readout takes the first maximum in the same action order); the share of calls with tied top classes is reported.",
    "ANM runs through the engine path of finite_field_runner.run_instances (preparation apparatus, "
    "validate_source_events, build_graph, evolve_field, readout) with the registered field; event ids "
    "'<lineage>:<protein>'; Q1 and Q2 share one field per cell (different readout thresholds); Q3 has its own "
    "instance (one anchor event per class); the closure readout uses readout_coordinates over each action's event "
    "sites with the registered closure weights. E1.1a counts call mismatches engine vs re-coded rule for Q1, Q2, Q3 "
    "and, for C2, engine closure vs re-coded closure (Q1, Q2 bars).",
    "E1.3 engine leave-one-out: each event of the called class is deleted in turn and the field re-evolved; the "
    "top-1 set is every event whose deletion lowers the called action's score by at least the maximum decrease "
    "minus 1e-12 * max(1, |maximum decrease|) (float tolerance); the engine flip is the call after deleting a top-1 "
    "event, read with the question's threshold. E1.3c uses the exact closed-form flip (A1.3); E1.3a checks that the "
    "engine agrees with it.",
    "E1.3b: shares are reported among all calls of a donor and within the called class; the permutation null uses "
    "1000 permutations seeded with seeds.global (no E1.3-specific seed is registered); two-sided permutation p-value "
    "around the null mean.",
    "E1.3c: m = number of flip-sensitive Q1 calls (scored or not); the comparison set is the m Q1 calls with the "
    "lowest margin (ties by the A1.1 permutation); error rate = wrong among scored cells of each set (a call on a "
    "key-OUT cell is wrong).",
    "E1.5 draws: for draw seed s, the labelled training pool (training cells the key scores) is sorted by cell id "
    "and permuted with numpy default_rng(s); a draw of size n is its first n cells (draws are nested across n). "
    "n = 0 and n = all are deterministic and repeated for every draw. The Q1 reference curve (primary classifier "
    "features and C, primary-key labels) has n >= 25 only (there is no label-free Q1 classifier).",
    "Q3 annotation-only sensitivity key: annotated CD14+ Mono -> myeloid, the other annotated myeloid types "
    "(CD16+ Mono, cDC2, cDC1) -> OUT, every other cell its annotation class.",
    "Key-validity flag (registration 4.3): decided by the test_primary 5-class kappa (annotation vs gate); the "
    "secondary split's flag by its own kappa. Annotation-only rows are reported next to every primary-key row.",
    "Bootstrap: count vectors; per replicate, donors drawn with numpy integers over the sorted donor list, then cells "
    "with replacement within each drawn donor; replicates generated in sequence from default_rng(seeds.bootstrap) "
    "for the pooled primary split; per-donor rows: the point within the donor plus a one-stage (cells within the "
    "donor) bootstrap from default_rng(seeds.bootstrap); the secondary split (one donor) likewise. Percentile "
    "intervals by numpy.percentile; replicates where an endpoint is undefined are excluded and counted.",
    "Secondary classifier (registration classifier.secondary, MLP, random_state = seeds.classifier 0-4): E1.4 rows "
    "only, mean over the five seeds; it has no registered bar and no operating point.",
    "Sensitivity panels: rule and classifier rows for E1.4 and the deployed operating points (bars and C fixed in "
    "this addendum); ANM on the primary panel only (registration E1.panels). Sensitivity keys annotation_only, "
    "primary_per_batch, primary_no_gdT158 for E1.4, E1.3c and C2; q3_annotation_only for E1.1c/d and E1.5.",
    "C2: the classifier arm (primary logistic classifier, its own argmax lineage) is reported with its "
    "differences from the mean rule as secondary rows; the verdict compares ANM closure with the mean rule only.",
    "Signs: E1.1c/d difference = rule minus classifier (win = the rule is better); E1.3c = error of flip-sensitive "
    "calls minus error of as many lowest-margin calls; E1.4 = top score minus comparator; C2 = ANM closure minus mean "
    "rule. E1.1b counts a change to or from no call as a change.",
    "The A1.1 tie-break applies to every matched-coverage selection (including E1.3c's lowest-margin set and the "
    "label-cost classifiers); an internal check confirms that the selection at the rule's own realised coverage "
    "reproduces its deployed calls.",
    "Smoke runs (inside the workflow) evaluate val donor 18303 cells only, never site4.",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def rnd(x) -> float:
    return float(round(float(x), 6))


def r6(x):
    if isinstance(x, (float, np.floating)):
        return float(round(float(x), 6))
    if isinstance(x, dict):
        return {k: r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r6(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def content_bytes(d: dict) -> bytes:
    return (json.dumps(d, indent=1, sort_keys=False, ensure_ascii=False) + "\n").encode()


# ============================================================================ computed core (train + val only)
def build_core(args) -> dict:
    t0 = time.time()
    reg = va_.load_registration_amended(args.registration, args.amendment)
    a1 = reg["amendment_A1"]["computed"]["classifiers"]
    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    split = npz["split"].astype(str)
    adt_names = [str(x) for x in npz["adt_names"]]
    n_all = int(split.size)

    # ---- the only rows whose protein / annotation / embedding are read: split in {train, val}
    nt = np.where(split != "test")[0]
    adt = np.asarray(npz["adt"], dtype=np.float32)[nt]
    ct = npz["cell_types"].astype(str)[nt]
    zmm = np.load(args.z, mmap_mode="r")
    if zmm.shape[0] != n_all:
        raise SystemExit(f"z rows {zmm.shape[0]} != cells {n_all}")
    z = np.asarray(zmm[nt], dtype=np.float32)
    del zmm
    sp = split[nt]
    tr, va = sp == "train", sp == "val"
    j = {n: i for i, n in enumerate(adt_names)}
    log(f"loaded {nt.size} train+val rows of {n_all}")

    keys = vk.build_keys(ct, adt, adt_names, reg)
    prim, q3k = keys["primary"], keys["q3"]
    v = vk.evidence(vk.head_predict(z, args.ckpt, device="cpu", size_factor=1.0), reg)
    from sklearn.metrics import log_loss

    def val_ll(model, X, y):
        m = va & (y != vk.UNSCORED)
        return float(log_loss(y[m], e1.proba5(model, X[m])[:, [vk.CLASSES.index(str(c)) for c in model.classes_]],
                              labels=model.classes_))

    def fit(feats, y, C):
        X = v[:, [j[p] for p in feats]]
        m = tr & (y != vk.UNSCORED)
        return e1.fit_logreg(X[m], y[m], C), X

    def lineage_bar(model, X, q):
        conf, _ = e1.lineage_conf_call(e1.proba5(model, X[va]))
        return rnd(np.quantile(conf, NO_CALL_TARGETS[q]))

    core: dict = {"experiment": "E1", "builder": BUILDER_VERSION}

    # ---------------------------------------------------------------- reproduce the registered numbers E1 reads
    rep: dict = {}
    for q in ("Q1", "Q2", "Q3"):
        top = vk.class_scores(v[va], adt_names, reg, q).max(axis=1)
        rep[f"rule_bar_{q}"] = {"registered": reg["questions"][q]["bar"], "recomputed": rnd(np.quantile(top, NO_CALL_TARGETS[q]))}
    for q in ("Q1", "Q2"):
        cl = vk.anm_action_scores(v[va], adt_names, reg, q, readout="closure").max(axis=1)
        rep[f"closure_bar_{q}"] = {"registered": reg["anm"][f"closure_bar_{q}"], "recomputed": rnd(np.quantile(cl, NO_CALL_TARGETS[q]))}
    cp = reg["classifier"]["primary"]
    m_p, X_p = fit(cp["feature_proteins"], prim, cp["C"])
    reg_ll = {r["C"]: r["val_log_loss"] for r in cp["val_log_loss"]}
    rep["classifier_primary_val_log_loss"] = {"registered": reg_ll[cp["C"]], "recomputed": rnd(val_ll(m_p, X_p, prim))}
    for q in ("Q1", "Q2"):
        rep[f"classifier_primary_bar_{q}"] = {"registered": cp["bars"][q], "recomputed": lineage_bar(m_p, X_p, q)}
    cq = reg["classifier"]["q3"]
    m_q3, X_q3 = fit(cq["feature_proteins"], q3k, cq["C"])
    reg_ll3 = {r["C"]: r["val_log_loss"] for r in cq["val_log_loss"]}
    rep["classifier_q3_val_log_loss"] = {"registered": reg_ll3[cq["C"]], "recomputed": rnd(val_ll(m_q3, X_q3, q3k))}
    a1_models = {}
    for name, labels in (("q1labels_q3features", prim), ("q1labels_anchors", prim), ("q3labels_anchors", q3k)):
        spec = a1[name]
        mdl, Xm = fit(spec["features"], labels, spec["C"])
        a1_models[name] = (mdl, Xm)
        ll = {r["C"]: r["val_log_loss"] for r in spec["val_log_loss"]}
        rep[f"A1_{name}_val_log_loss"] = {"registered": ll[spec["C"]], "recomputed": rnd(val_ll(mdl, Xm, labels))}
    for r in rep.values():
        r["match"] = bool(abs(float(r["registered"]) - float(r["recomputed"])) <= 1.5e-6)
    rep_ok = all(r["match"] for r in rep.values())
    core["reproduced_registered_numbers"] = {"all_match": rep_ok, "items": rep}
    if not rep_ok:
        bad = {k: r for k, r in rep.items() if not r["match"]}
        raise SystemExit(f"registered numbers not reproduced (E1 would not run the registered methods): {bad}")
    log("registered bars, closure bars, classifier log-loss and bars reproduce")

    # ---------------------------------------------------------------- sensitivity panels (section 7 / 8 procedures)
    grid = list(reg["classifier"]["primary"]["C_grid"])
    sens = {}
    for name in reg["panels"]["sensitivity"]:
        S = vk.class_scores(v[va], adt_names, reg, "Q1", panel=name)
        top = S.max(axis=1)
        bars = {q: rnd(np.quantile(top, NO_CALL_TARGETS[q])) for q in ("Q1", "Q2")}
        feats = [p for k in vk.LINEAGES for p in reg["panels"][name]["proteins"][k]]
        res = []
        for C in grid:
            mdl, Xs = fit(feats, prim, C)
            res.append({"C": C, "val_log_loss": rnd(val_ll(mdl, Xs, prim))})
        bestC = min(res, key=lambda r: (r["val_log_loss"], r["C"]))["C"]
        mdl, Xs = fit(feats, prim, bestC)
        sens[name] = {
            "proteins": reg["panels"][name]["proteins"],
            "rule_bars": bars,
            "rule_val_no_call_realised": {q: rnd(np.mean(top < bars[q])) for q in bars},
            "classifier": {"features": feats, "labels": "training-cell primary key (5 classes incl. OUT)",
                           "C_grid": grid, "val_log_loss": res, "C": bestC,
                           "C_selection": "lowest val log-loss on val primary-key cells, ties -> smaller C",
                           "n_train_labelled": int((tr & (prim != vk.UNSCORED)).sum()),
                           "bars": {q: lineage_bar(mdl, Xs, q) for q in ("Q1", "Q2")}},
        }
        log(f"sensitivity panel {name}: bars {bars}, C {bestC}")
    core["sensitivity_panels"] = sens

    # ---------------------------------------------------------------- Q3 operating-point bars (reporting only)
    q3b = {"registered_q3_classifier": {"C": cq["C"], "features": cq["feature_proteins"], "labels": "q3 key",
                                        "bar_Q3": lineage_bar(m_q3, X_q3, "Q3")}}
    for name, (mdl, Xm) in a1_models.items():
        q3b[name] = {"C": a1[name]["C"], "features": a1[name]["features"], "bar_Q3": lineage_bar(mdl, Xm, "Q3")}
    core["q3_operating_point_bars"] = {"no_call_target": NO_CALL_TARGETS["Q3"],
                                       "use": "deployed operating point only; every E1 comparison is at matched coverage",
                                       "classifiers": q3b}
    log(f"Q3 operating-point bars: { {k: v_['bar_Q3'] for k, v_ in q3b.items()} }")
    log(f"core built in {time.time() - t0:.1f}s")
    return r6(core)


# ============================================================================ declared content
def _git_head() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
    except OSError:
        return None


def assemble(core: dict, args) -> dict:
    reg_sha = vk.sha256_file(args.registration)
    am_sha = vk.sha256_file(args.amendment)
    return {
        "addendum_to": "registration/registration_v3.json experiments.E1, as amended by registration/amendment_A1.json",
        "experiment": "E1",
        "addendum_version": 1,
        "registration_sha256": reg_sha,
        "amendment_A1_sha256": am_sha,
        "fixed_before_site4": "computed from training and validation cells only and committed before any site4 "
                              "evaluation of E1; smoke runs of the E1 builder use val donor 18303 cells only",
        "scope": "numbers the registration leaves open, fixed by the procedures it writes (bars: section 7 val "
                 "quantile at the declared no-call rate; classifier C: section 8 lowest val log-loss on the "
                 "registered grid; classifier bars: section 8 bar_rule), plus the implementation choices of the E1 "
                 "builder that the registration does not fix. No registered endpoint, key, cell set, margin or "
                 "comparator is changed (A1.10).",
        "computed": core,
        "implementation_choices": IMPLEMENTATION_CHOICES,
        "provenance": {
            "builder": BUILDER_VERSION,
            "builder_sha256": vk.sha256_file(Path(__file__)),
            "v3_key_sha256": vk.sha256_file(ROOT / "bridge_anm/lib/v3_key.py"),
            "v3_amend_sha256": vk.sha256_file(ROOT / "bridge_anm/lib/v3_amend.py"),
            "v3_e1_sha256": vk.sha256_file(ROOT / "bridge_anm/lib/v3_e1.py"),
            "core_sha256": sha256_bytes(content_bytes(core)),
            "inputs": {"cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz"),
                       "z_sha256": vk.sha256_file(args.z), "ckpt_sha256": vk.sha256_file(args.ckpt)},
            "git_head_at_build": _git_head(),
            "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "test_rows_read": "split metadata only (protein, cell types and embedding read on non-test rows)",
        },
    }


def write_addendum(doc: dict, path: Path = ADDENDUM, hashes: Path = HASHES) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    b = content_bytes(doc)
    path.write_bytes(b)
    sha = sha256_bytes(b)
    lines = [ln for ln in (hashes.read_text().splitlines() if hashes.exists() else []) if not ln.rstrip().endswith("  E1.json")]
    hashes.write_text("\n".join(lines + [f"{sha}  E1.json"]) + "\n")
    return sha


# ============================================================================ leakage check
def leakcheck(args) -> int:
    import tempfile

    import v3_leakage_check as lc

    lc.BUILDERS["E1_addendum"] = (Path(__file__).resolve(), CORE_NAME)
    py = sys.executable
    rep: dict = {"script": "bridge_anm/v3_e1_addendum.py --leakcheck",
                 "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    rep["static"] = lc.static_check(Path(__file__).resolve(),
                                    (lc.KEYLIB, lc.AMENDLIB, ROOT / "bridge_anm/lib/v3_e1.py"))
    log(f"static: {'PASS' if rep['static']['passed'] else rep['static']['issues']}")
    tmp = Path(tempfile.mkdtemp(prefix="v3_e1_leak_", dir=args.work_dir))
    try:
        real = lc.build_core(args.processed, args.z, args.ckpt, tmp / "real", py, "E1_addendum")
        rep["real_core_sha256"] = sha256_bytes(real)
        if ADDENDUM.exists():
            rec = json.loads(ADDENDUM.read_text())["provenance"]["core_sha256"]
            rep["addendum_core_sha256"] = rec
            rep["real_matches_addendum"] = rec == sha256_bytes(real)
        src_npz, src_z = args.processed / "cite_arrays.npz", args.z
        npz_p, z_p, info = lc.write_poisoned(src_npz, src_z, tmp / "poison_test", "test", 101)
        b = lc.build_core(npz_p.parent, z_p, args.ckpt, tmp / "core_poison_test", py, "E1_addendum")
        rep["poison_test"] = {**info, "core_sha256": sha256_bytes(b), "identical_to_real": b == real}
        log(f"site4 poisoned: identical = {b == real}")
        npz_m, z_m, info_m = lc.write_poisoned(src_npz, src_z, tmp / "poison_val_mild", "val", 303, fraction=0.05,
                                               protein_only=True)
        try:
            bm = lc.build_core(npz_m.parent, z_m, args.ckpt, tmp / "core_poison_val_mild", py, "E1_addendum")
            rep["poison_val_mild_control"] = {**info_m, "core_sha256": sha256_bytes(bm), "differs_from_real": bm != real,
                                              "builder_outcome": "completed"}
        except RuntimeError as exc:
            rep["poison_val_mild_control"] = {**info_m, "core_sha256": None, "differs_from_real": True,
                                              "builder_outcome": "stopped: " + str(exc).strip().splitlines()[-1][:300]}
        log(f"val (5% protein) poisoned: differs = {rep['poison_val_mild_control']['differs_from_real']}")
    finally:
        import shutil

        shutil.rmtree(tmp, ignore_errors=True)
    rep["passed"] = bool(rep["static"]["passed"] and rep["poison_test"]["identical_to_real"]
                         and rep["poison_val_mild_control"]["differs_from_real"]
                         and rep.get("real_matches_addendum", True))
    rep["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    js = json.dumps(rep, indent=1) + "\n"
    args.report_dir.mkdir(parents=True, exist_ok=True)
    (args.report_dir / "E1_addendum_leakage_check.json").write_text(js)
    (ROOT / "registration" / "addenda" / "E1_leakage_check.json").write_text(js)
    log("PASS" if rep["passed"] else "FAIL")
    return 0 if rep["passed"] else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration", type=Path, default=ROOT / "registration/registration_v3.json")
    p.add_argument("--amendment", type=Path, default=ROOT / "registration/amendment_A1.json")
    p.add_argument("--core-only", action="store_true")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--leakcheck", action="store_true")
    p.add_argument("--work-dir", type=Path, default=None, help="parent of the temporary poisoned copies")
    p.add_argument("--report-dir", type=Path, default=Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/E1"))
    args = p.parse_args(argv)
    if args.leakcheck:
        return leakcheck(args)
    core = build_core(args)
    if args.core_only:
        out = args.out_dir or Path(".")
        out.mkdir(parents=True, exist_ok=True)
        (out / CORE_NAME).write_bytes(content_bytes(core))
        log(f"wrote {out / CORE_NAME} (sha256 {sha256_bytes(content_bytes(core))})")
        return 0
    doc = assemble(core, args)
    sha = write_addendum(doc)
    log(f"wrote {ADDENDUM} (sha256 {sha}) and its line in {HASHES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
