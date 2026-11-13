"""Declaration loading, defaults and hashing.

Schema (version 1). Keys marked * are required.

anm_declaration*: 1
name*: short id
status: exploratory | sealed            (sealed requires data_hashes and code_hash)
kind: sufficiency | intervention        (intervention studies edit the system; not handled by run())
system*:
  id*: e.g. TEDDY-G-70M
  frozen*: true
  inputs*: [measurement ids the system reads]
measurements*:                          (what was measured in each cell: the provenance roots)
  <id>: {modality: rna|protein|..., kind: counts|values}
source*:
  controlled*: true                     (a controlled perturbation; observational contrasts are refused)
  kind: crispri | crispr_ko | drug | ...
  control*: label of control cells
  treated*: regex or list of labels forming the treated class
  cluster: target                       (unit for folds and bootstrap; controls are their own units)
  stratum: a value, or a list of values, of the adapter's stratum column (e.g. a culture condition)
  assign_from_guides: false             (true: label cells from their detected guides; with control_guides,
                                         the gene names of the non-targeting guides. A cell is "control" when
                                         every guide is a control guide, gene G when every guide targets G or is
                                         a control guide, and is left out with no guide or guides for two genes)
  max_control_cells, control_seed       (a seeded subsample of the controls, for compute)
readouts*:
  - id*: name
    measurement*: a measurement id (never a node: a readout is an observation, not a model output)
    statistic*: e.g. mean_log1p_cp10k | adt_clr_mean
    features*: adapter-specific selector (gene regex, protein names, ...)
nodes*:                                 (candidate retained states and probes)
  - id*: name
    from*: a measurement id or another node id
    op*: operation name (adapter registry)
    operator*: system | analysis | approximation
    role*: chain | probe | positive_control | technical
chain*: [node ids in computation order]  (first = reference view of the input; all on the real path)
measurement_sharing: none | count_split | declared_shared   (needed when a readout uses a system input)
readout_classes: {ridge: {...}, hgb: {...}}
inference: {folds: 5, bootstrap: 1000, forward_bootstrap: 300, seed: 0}
decision: {margin_u: 0.10, alpha: 0.05}
"""
import copy
import hashlib
import json
from pathlib import Path

DEFAULTS = {
    "status": "exploratory",
    "kind": "sufficiency",
    "measurement_sharing": "none",
    "readout_classes": {
        "ridge": {"regressor": "RidgeCV", "alphas": [10.0 ** k for k in (-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)],
                  "classifier": "LogisticRegression", "C": 0.05},
        "hgb": {"regressor": "HistGradientBoostingRegressor", "classifier": "HistGradientBoostingClassifier",
                "max_iter": 150, "learning_rate": 0.05, "max_leaf_nodes": 15},
    },
    "inference": {"folds": 5, "bootstrap": 1000, "forward_bootstrap": 300, "seed": 0, "standardize_on": "control"},
    "decision": {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05},
}


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_declaration(path_or_dict):
    """Load a declaration from a .yaml/.yml/.json path or a dict; fill documented defaults."""
    if isinstance(path_or_dict, dict):
        raw = path_or_dict
    else:
        p = Path(path_or_dict)
        text = p.read_text()
        if p.suffix in (".yaml", ".yml"):
            import yaml
            raw = yaml.safe_load(text)
        else:
            raw = json.loads(text)
    return _merge(DEFAULTS, raw)


def declaration_hash(decl):
    """sha256 of the canonical JSON of the filled declaration (what the run actually used)."""
    blob = json.dumps(decl, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(blob).hexdigest()


def node_map(decl):
    return {n["id"]: n for n in decl.get("nodes", [])}


def ancestors(decl, node_id):
    """Measurement ids and node ids upstream of a node (following `from`)."""
    nodes = node_map(decl)
    seen, out_meas, out_nodes = set(), set(), []
    cur = node_id
    while cur is not None and cur not in seen:
        seen.add(cur)
        n = nodes.get(cur)
        if n is None:
            if cur in decl.get("measurements", {}):
                out_meas.add(cur)
            break
        src = n.get("from")
        if isinstance(src, list):          # multi-input node
            for s in src:
                m, ns = ancestors(decl, s)
                out_meas |= m
                out_nodes += ns + ([s] if s in nodes else [])
            break
        if src in nodes:
            out_nodes.append(src)
        cur = src
    return out_meas, out_nodes
