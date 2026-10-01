#!/usr/bin/env python3
"""Prove that no site4 (test) label or measured site4 protein entered any v3 registration choice.

Three checks; all must pass.

1. Static: in the registration builder, every read of the protein matrix, the cell
   types and the embedding is restricted to the non-test rows (``[nt]`` with
   ``nt = split != "test"``) and nothing reads an earlier site4 result file.
2. Poisoning (invariance): rebuild the registration core from inputs whose site4
   rows of protein, cell type and embedding are replaced by random values. The core
   must be byte-identical to the core built from the real inputs (and to the core
   whose sha256 the frozen registration records).
3. Sensitivity (positive control): poison the *val* rows the same way. The core
   must change -- this shows check 2 would have caught a dependence.

Writes leakage_check.json and LEAKAGE_CHECK.md to --report-dir (default
outputs/v3/registration) and a copy of the JSON next to the registration.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "bridge_anm" / "v3_build_registration.py"
KEYLIB = ROOT / "bridge_anm" / "lib" / "v3_key.py"
NEEDED_KEYS = ("split", "sites", "donors", "adt_names", "adt", "cell_types")
FORBIDDEN_READS = ("test_per_protein", "adt_true", "hard_proof", "anm_cite_bridge", "scope_refine",
                   "mode_a_official", "keys_site4", "pred_test", "idx_test")


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def static_check() -> dict:
    src = BUILDER.read_text()
    res = {"file": str(BUILDER.relative_to(ROOT)), "issues": []}
    m = re.search(r'nt = np\.where\(split != "test"\)\[0\]', src)
    res["nontest_selector_present"] = bool(m)
    if not m:
        res["issues"].append("non-test selector nt not found")
    for pat in (r'npz\["adt"\]', r'npz\["cell_types"\]', r"zmm\["):
        for mm in re.finditer(pat, src):
            tail = src[mm.end(): mm.end() + 60]
            if not re.match(r'(?:\.astype\(str\)|, dtype=np\.float32\))?\[nt\]|nt\]', tail) and "[nt]" not in tail[:40]:
                res["issues"].append(f"{pat} read without [nt] at offset {mm.start()}: {tail[:40]!r}")
    for word in FORBIDDEN_READS:
        for f in (BUILDER, KEYLIB):
            if word in f.read_text():
                res["issues"].append(f"{f.name} mentions {word!r}")
    res["passed"] = not res["issues"]
    return res


def write_poisoned(src_npz: Path, src_z: Path, out_dir: Path, rows: str, seed: int, fraction: float = 1.0,
                   protein_only: bool = False) -> tuple[Path, Path, dict]:
    """Copy the needed arrays; replace the chosen rows' protein (and, unless protein_only, cell type and
    embedding) by random values. ``fraction`` < 1 poisons a seeded subset of those rows."""
    out_dir.mkdir(parents=True, exist_ok=True)
    z = np.load(src_npz, allow_pickle=False)
    arrs = {k: z[k] for k in NEEDED_KEYS}
    split = arrs["split"].astype(str)
    sel = np.where(split == rows)[0]
    rng = np.random.default_rng(seed)
    if fraction < 1.0:
        sel = np.sort(rng.choice(sel, size=max(1, int(round(fraction * sel.size))), replace=False))
    adt = np.array(arrs["adt"], dtype=np.float32, copy=True)
    adt[sel] = rng.uniform(0.0, 8.0, size=(sel.size, adt.shape[1])).astype(np.float32)
    ct = arrs["cell_types"].copy()
    if not protein_only:
        types = np.array(sorted(set(arrs["cell_types"].astype(str))), dtype=ct.dtype)
        ct[sel] = types[rng.integers(0, types.size, size=sel.size)]
    arrs["adt"], arrs["cell_types"] = adt, ct
    np.savez(out_dir / "cite_arrays.npz", **arrs)
    zz = np.load(src_z, mmap_mode="r")
    zp = np.lib.format.open_memmap(out_dir / "z_rna.npy", mode="w+", dtype=zz.dtype, shape=zz.shape)
    zp[:] = zz[:]
    if not protein_only:
        zp[sel] = rng.normal(0.0, 1.0, size=(sel.size, zz.shape[1])).astype(zz.dtype)
    zp.flush()
    del zp
    info = {"rows_poisoned": rows, "n_rows": int(sel.size), "fraction": fraction, "seed": seed,
            "what": "adt rows -> U(0, 8)" + ("" if protein_only else "; cell_types -> random registered type; z rows -> N(0, 1)")}
    return out_dir / "cite_arrays.npz", out_dir / "z_rna.npy", info


def build_core(processed_dir: Path, z: Path, ckpt: Path, out_dir: Path, python: str) -> bytes:
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "4")
    cmd = [python, str(BUILDER), "--core-only", "--processed", str(processed_dir), "--z", str(z),
           "--ckpt", str(ckpt), "--out-dir", str(out_dir)]
    t = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"builder failed ({r.returncode}):\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    print(f"  built core in {time.time() - t:.0f}s -> {out_dir}", flush=True)
    return (out_dir / "registration_v3_core.json").read_bytes()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration", type=Path, default=ROOT / "registration/registration_v3.json")
    p.add_argument("--work-dir", type=Path, default=None, help="scratch dir for poisoned copies (deleted after)")
    p.add_argument("--report-dir", type=Path,
                   default=Path("/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/registration"))
    p.add_argument("--skip-sensitivity", action="store_true")
    args = p.parse_args(argv)
    py = sys.executable
    report = {"script": str(Path(__file__).relative_to(ROOT)), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    print("1. static check", flush=True)
    report["static"] = static_check()
    print("  ", "PASS" if report["static"]["passed"] else f"FAIL {report['static']['issues']}", flush=True)

    reg = json.loads(args.registration.read_text()) if args.registration.exists() else None
    recorded_core = (reg or {}).get("provenance", {}).get("core_sha256")
    tmp_parent = args.work_dir or Path(tempfile.mkdtemp(prefix="v3_leak_"))
    tmp_parent.mkdir(parents=True, exist_ok=True)
    try:
        print("2. real inputs", flush=True)
        real = build_core(args.processed, args.z, args.ckpt, tmp_parent / "core_real", py)
        report["real_core_sha256"] = sha(real)
        report["registration_core_sha256"] = recorded_core
        report["real_matches_registration"] = (recorded_core == sha(real)) if recorded_core else None

        print("3. site4 rows poisoned", flush=True)
        npz_p, z_p, info = write_poisoned(args.processed / "cite_arrays.npz", args.z, tmp_parent / "poison_test", "test", 101)
        poisoned = build_core(npz_p.parent, z_p, args.ckpt, tmp_parent / "core_poison_test", py)
        report["poison_test"] = {**info, "core_sha256": sha(poisoned), "identical_to_real": poisoned == real}
        print("   identical:", poisoned == real, flush=True)

        if not args.skip_sensitivity:
            print("4. val rows poisoned (positive control)", flush=True)
            npz_v, z_v, info_v = write_poisoned(args.processed / "cite_arrays.npz", args.z, tmp_parent / "poison_val", "val", 202)
            try:
                pv = build_core(npz_v.parent, z_v, args.ckpt, tmp_parent / "core_poison_val", py)
                report["poison_val_control"] = {**info_v, "core_sha256": sha(pv), "differs_from_real": pv != real,
                                                "builder_outcome": "completed"}
            except RuntimeError as exc:
                # random val rows can leave a class with no protein above the AUROC floor; the build then
                # stops. That is also a change of output, i.e. the check detects the dependence.
                report["poison_val_control"] = {**info_v, "core_sha256": None, "differs_from_real": True,
                                                "builder_outcome": "stopped: " + str(exc).strip().splitlines()[-1][:300]}
            print("   differs:", report["poison_val_control"]["differs_from_real"],
                  report["poison_val_control"]["builder_outcome"][:120], flush=True)
            print("5. 5% of val rows, protein only (mild positive control)", flush=True)
            npz_m, z_m, info_m = write_poisoned(args.processed / "cite_arrays.npz", args.z, tmp_parent / "poison_val_mild",
                                                "val", 303, fraction=0.05, protein_only=True)
            pm = build_core(npz_m.parent, z_m, args.ckpt, tmp_parent / "core_poison_val_mild", py)
            report["poison_val_mild_control"] = {**info_m, "core_sha256": sha(pm), "differs_from_real": pm != real}
            print("   differs:", pm != real, flush=True)
    finally:
        if args.work_dir is None:
            import shutil

            shutil.rmtree(tmp_parent, ignore_errors=True)

    report["passed"] = bool(report["static"]["passed"] and report["poison_test"]["identical_to_real"]
                            and (report.get("poison_val_control", {}).get("differs_from_real", True))
                            and (report.get("poison_val_mild_control", {}).get("differs_from_real", True))
                            and report["real_matches_registration"] in (True, None))
    report["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.report_dir.mkdir(parents=True, exist_ok=True)
    js = json.dumps(report, indent=1) + "\n"
    (args.report_dir / "leakage_check.json").write_text(js)
    (ROOT / "registration" / "leakage_check_report.json").write_text(js)
    md = [
        "# v3 leakage check", "",
        f"Result: **{'PASS' if report['passed'] else 'FAIL'}** ({report['finished_utc']})", "",
        "| check | result |", "|---|---|",
        f"| static: builder reads protein / cell types / embedding only on non-test rows; no earlier site4 result file | {'pass' if report['static']['passed'] else 'fail: ' + '; '.join(report['static']['issues'])} |",
        f"| core built from the real inputs equals the core recorded in the registration | {report['real_matches_registration']} |",
        f"| site4 protein, cell types and embedding replaced by random values: core byte-identical | {report['poison_test']['identical_to_real']} |",
    ]
    if "poison_val_control" in report:
        md.append(f"| positive control, all val rows replaced: output changes | {report['poison_val_control']['differs_from_real']} ({report['poison_val_control']['builder_outcome'][:80]}) |")
    if "poison_val_mild_control" in report:
        md.append(f"| mild positive control, protein of 5% of val rows replaced: core changes | {report['poison_val_mild_control']['differs_from_real']} |")
    md += ["", f"Real core sha256: `{report['real_core_sha256']}`", "",
           "The poisoned build reads the same split / site / donor metadata, so the registered splits and the "
           "label-free E2 subset are unchanged by construction; everything else in the core (gate thresholds, "
           "key quality, evidence normaliser, panels, bars, classifier settings, E3 flag threshold) is shown not "
           "to depend on any site4 protein, label or embedding value.", ""]
    (args.report_dir / "LEAKAGE_CHECK.md").write_text("\n".join(md))
    print("PASS" if report["passed"] else "FAIL", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
