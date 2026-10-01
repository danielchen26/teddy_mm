#!/usr/bin/env python3
"""Render registration/REGISTRATION_v3.md from REGISTRATION_v3.template.md and registration_v3.json.

Every number in the rendered document comes from the registration JSON, so the prose
cannot drift from the frozen values. Placeholders:

  {{a.b.c}}            value at that path (floats: 4 significant decimals)
  {{a.b.c|.6f}}        with a format spec
  {{TABLE:name}}       a generated table (see TABLES)
  {{SHA256}}           sha256 of registration_v3.json (from its .sha256 file)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REG_DIR = ROOT / "registration"


def get(reg, path):
    cur = reg
    for part in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        else:
            cur = cur[part]
    return cur


def fmt(v, spec=None):
    if spec:
        return format(v, spec)
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.4f}"
    if isinstance(v, list):
        return ", ".join(str(x) for x in v)
    return str(v)


def t_annotation(reg):
    counts = reg["data"]["cell_types_per_split"]
    notes = reg["annotation_map_notes"]
    order = {c: i for i, c in enumerate(reg["classes"])}
    rows = sorted(reg["annotation_map"].items(), key=lambda kv: (order[kv[1]], kv[0]))
    out = ["| annotated cell type | class | train | val | note |", "|---|---|---:|---:|---|"]
    for t, c in rows:
        out.append(f"| {t} | {c} | {counts[t]['train']:,} | {counts[t]['val']:,} | {notes.get(t, '')} |")
    return "\n".join(out)


def t_gate_thresholds(reg):
    g = reg["gate"]
    used = set(g["b_markers"]) | set(g["nk_markers"]) | {"CD3", "CD14", "CD33", "CD11c", "CD71"}
    out = ["| protein | threshold (CLR) | role in the chosen gate |", "|---|---:|---|"]
    role = {"CD3": "T", "CD19": "B marker", "CD20": "B marker", "CD56": "NK marker", "CD94": "NK marker",
            "CD335": "NK marker", "CD16": "NK marker", "CD14": "myeloid; excludes T, B, NK; Q3 key",
            "CD33": "myeloid; excludes NK", "CD11c": "myeloid", "CD71": "erythroid -> OUT"}
    for p in g["proteins"]:
        r = role[p] if p in used else "not used (B-marker variant not chosen)"
        out.append(f"| {p} | {g['thresholds'][p]:.6f} | {r} |")
    return "\n".join(out)


def t_gate_grid(reg):
    out = ["| threshold method | B markers | NK markers | mean train-batch kappa | min train-batch kappa | val kappa |",
           "|---|---|---|---:|---:|---:|"]
    for r in sorted(reg["gate"]["selection"]["grid_results"], key=lambda r: -r["mean_train_batch_kappa"]):
        out.append(f"| {r['threshold_method']} | {', '.join(r['b_markers'])} | {', '.join(r['nk_markers'])} | "
                   f"{r['mean_train_batch_kappa']:.4f} | {r['min_train_batch_kappa']:.4f} | {r['val_kappa5']:.4f} |")
    return "\n".join(out)


def t_key_quality(reg):
    kq = reg["key_quality"]
    out = ["| split | cells | agreement | kappa (5 classes) | B recall / precision | T | NK | myeloid | OUT |",
           "|---|---:|---:|---:|---|---|---|---|---|"]
    for sp in ("train", "val"):
        r = kq[sp]
        cells = [f"{r['per_class'][c]['recall_vs_annotation']:.3f} / {r['per_class'][c]['precision_vs_annotation']:.3f}"
                 for c in ("B", "T", "NK", "myeloid", "OUT")]
        out.append(f"| {sp} | {r['n']:,} | {r['agreement']:.4f} | {r['kappa5']:.4f} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def t_key_batches(reg):
    out = ["| batch (site, donor) | cells | agreement | kappa |", "|---|---:|---:|---:|"]
    for b, r in reg["key_quality"]["per_batch"].items():
        out.append(f"| {b.replace('|', ' / ')} | {r['n']:,} | {r['agreement']:.4f} | {r['kappa5']:.4f} |")
    return "\n".join(out)


def t_key_counts(reg):
    kq = reg["key_quality"]
    cols = ["B", "T", "NK", "myeloid", "OUT", "unscored"]
    out = ["| split | " + " | ".join(cols) + " |", "|---|" + "---:|" * len(cols)]
    for sp in ("train", "val"):
        out.append(f"| {sp} | " + " | ".join(f"{kq[sp]['primary_key_counts'][c]:,}" for c in cols) + " |")
    return "\n".join(out)


def t_key_per_type(reg):
    pt = reg["key_quality"]["per_cell_type_train_val"]
    out = ["| annotated type | class | train+val cells | gated B / T / NK / myeloid / OUT | in primary key |",
           "|---|---|---:|---|---:|"]
    order = {c: i for i, c in enumerate(reg["classes"])}
    for t, r in sorted(pt.items(), key=lambda kv: (order[kv[1]["annotation_class"]], kv[0])):
        g = r["gated"]
        out.append(f"| {t} | {r['annotation_class']} | {r['n_train_val']:,} | "
                   f"{g['B']} / {g['T']} / {g['NK']} / {g['myeloid']} / {g['OUT']} | {r['in_primary_key']:,} |")
    return "\n".join(out)


def t_nk(reg):
    nk = reg["panels"]["selection"]["nk_must_consider"]
    out = ["| protein | worst-pair AUROC (NK) | NK vs rest AUROC | NK vs T AUROC | assigned class | rank for NK |",
           "|---|---:|---:|---:|---|---:|"]
    for p, r in nk.items():
        out.append(f"| {p} | {r['worst_pair_auroc_nk']:.4f} | {r['auroc_nk_vs_rest']:.4f} | {r['auroc_nk_vs_t']:.4f} | "
                   f"{r['assigned_class']} | {r['rank_for_nk']} |")
    return "\n".join(out)


def t_top10(reg):
    out = []
    for k, rows in reg["panels"]["selection"]["top10_per_class"].items():
        out.append(f"**{k}**: " + "; ".join(f"{r['protein']} {r['worst_pair_auroc']:.4f}"
                                            + ("" if r["assigned"] == k else f" (assigned {r['assigned']})")
                                            for r in rows))
        out.append("")
    return "\n".join(out).rstrip()


def t_panels(reg):
    names = ["primary"] + reg["panels"]["sensitivity"]
    out = ["| panel | B | T | NK | myeloid | role |", "|---|---|---|---|---|---|"]
    for n in names:
        p = reg["panels"][n]["proteins"]
        role = "primary" if n == "primary" else reg["panels"][n].get("note", "sensitivity")
        out.append(f"| {n} | {', '.join(p['B'])} | {', '.join(p['T'])} | {', '.join(p['NK'])} | {', '.join(p['myeloid'])} | {role} |")
    return "\n".join(out)


def t_questions(reg):
    out = ["| question | evidence per class | bar | val no-call target | val no-call realised | val accuracy of calls | val decision accuracy | val OUT decline |",
           "|---|---|---:|---:|---:|---:|---:|---:|"]
    for q in ("Q1", "Q2", "Q3"):
        r = reg["questions"][q]
        o = r["val_outcome_rule"]
        ev = "primary panel, equal weights" if r["panel"] == "primary" else "anchor: " + ", ".join(f"{k} {p}" for k, p in r["anchors"].items())
        out.append(f"| {q} {r['name']} | {ev} | {r['bar']:.6f} | {r['val_no_call_target']:.2f} | {r['val_no_call_realised']:.4f} | "
                   f"{o['accuracy_of_calls']:.4f} | {o['decision_accuracy']:.4f} | {o['out_decline_rate']:.4f} |")
    return "\n".join(out)


def t_clf(reg):
    out = ["| C | val log-loss (primary classifier) | val log-loss (Q3 classifier) |", "|---:|---:|---:|"]
    for a, b in zip(reg["classifier"]["primary"]["val_log_loss"], reg["classifier"]["q3"]["val_log_loss"]):
        out.append(f"| {a['C']} | {a['val_log_loss']:.6f} | {b['val_log_loss']:.6f} |")
    return "\n".join(out)


def t_splits(reg):
    s = reg["splits"]
    out = ["| split | rule | cells | index sha256 (first 16) |", "|---|---|---:|---|"]
    for n in ("train", "val", "test_primary", "test_secondary"):
        out.append(f"| {n} | {s[n]['rule']} | {s[n]['n']:,} | `{s[n]['sha256_indices'][:16]}` |")
    return "\n".join(out)


TABLES = {"annotation_map": t_annotation, "gate_thresholds": t_gate_thresholds, "gate_grid": t_gate_grid,
          "key_quality": t_key_quality, "key_batches": t_key_batches, "key_counts": t_key_counts,
          "key_per_type": t_key_per_type, "nk_candidates": t_nk, "top10": t_top10, "panels": t_panels,
          "questions": t_questions, "classifier": t_clf, "splits": t_splits}


def render(template: str, reg: dict, sha: str) -> str:
    def sub(m):
        token = m.group(1).strip()
        if token == "SHA256":
            return sha
        if token.startswith("TABLE:"):
            return TABLES[token[6:]](reg)
        path, _, spec = token.partition("|")
        return fmt(get(reg, path.strip()), spec.strip() or None)

    out = re.sub(r"\{\{(.+?)\}\}", sub, template)
    if "{{" in out:
        raise SystemExit("unresolved placeholder")
    return out


def main():
    reg_path = REG_DIR / "registration_v3.json"
    reg = json.loads(reg_path.read_text())
    sha = (REG_DIR / "registration_v3.json.sha256").read_text().split()[0]
    from lib.v3_key import sha256_file  # noqa: E402

    if sha256_file(reg_path) != sha:
        raise SystemExit("registration_v3.json does not match its .sha256 file")
    tpl = (REG_DIR / "REGISTRATION_v3.template.md").read_text()
    (REG_DIR / "REGISTRATION_v3.md").write_text(render(tpl, reg, sha))
    print("wrote", REG_DIR / "REGISTRATION_v3.md")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
