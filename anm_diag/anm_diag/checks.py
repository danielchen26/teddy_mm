"""Legality checks: refuse designs that cannot answer the ANM question.

Each rule cites what it protects. Errors stop run(); warnings are written into the report.

  S1  a controlled source with a control group                  (paper §2: the source is an intervention)
  R1  a readout is a measurement, never a node or model output  (Fig. 1: state and readout measured separately)
  R2  a readout's measurement must not feed any candidate,
      unless it is a declared system input with a sharing policy (else the readout leaks into the state)
  C1  chain nodes follow the real computation, one after another (only on-path stages are compared)
  C2  chain nodes after the reference are the system's own operators (no approximations on the path)
  C3  no analyst compression (PCA, k-means, ...) on the chain or as a candidate (§S1.4: candidates are
      declared maps; compression is a transform the analyst adds and must not hide inside a comparison)
  C4  the system stays fixed: no input masking, weight edits or state patches in a sufficiency study
      (§2: editing the component or its input would test a system intervention)
  C5  probes and controls hang off declared nodes; positive control and reference present
  E1  inverse statistics are continuous conditional tests, not partitions of the state
  K1  margin and alpha declared; coverage warning when there are few treated clusters
  K2  a sealed declaration carries data and code hashes
"""
import re

from .spec import ancestors, node_map

COMPRESSION_OPS = {"pca", "svd", "kmeans", "umap", "tsne", "nmf", "random_projection", "hvg_subset"}
INTERVENTION_OPS = {"mask_input", "drop_genes_from_input", "edit_weights", "clamp_state", "patch_state",
                    "mean_preserving_patch", "finetune"}
PARTITION_STATS = {"kmeans_partition", "cluster_tv", "discrete_partition"}
ALLOWED_INVERSE = {"crossfit_partial_linear"}


class CheckError(Exception):
    def __init__(self, problems):
        self.problems = problems
        super().__init__("\n".join(f"[{p['code']}] {p['message']}" for p in problems if p["severity"] == "error"))


def _p(code, severity, message):
    return {"code": code, "severity": severity, "message": message}


def check(decl, n_treated_clusters=None, raise_on_error=False):
    out = []
    req = ["name", "system", "measurements", "source", "readouts", "nodes", "chain"]
    for k in req:
        if k not in decl or decl[k] in (None, [], {}):
            out.append(_p("D1", "error", f"missing required key '{k}'"))
    if any(p["code"] == "D1" for p in out):
        if raise_on_error:
            raise CheckError(out)
        return out

    meas = decl["measurements"]
    nodes = node_map(decl)
    sysinp = set(decl["system"].get("inputs", []))
    kind = decl.get("kind", "sufficiency")

    # S1 source
    src = decl["source"]
    if not src.get("controlled", False):
        out.append(_p("S1", "error", "no controlled source perturbation: an ANM response test needs a source pushed by "
                                     "the experimenter and a control group (observational contrasts or decision layers "
                                     "without a source cannot test retained information)"))
    if not src.get("control") or not src.get("treated"):
        out.append(_p("S1", "error", "source must name the control label and the treated class"))

    # system
    if not decl["system"].get("frozen", False):
        out.append(_p("C4", "error", "system.frozen must be true for a sufficiency study"))
    for m in sysinp:
        if m not in meas:
            out.append(_p("D2", "error", f"system input '{m}' is not a declared measurement"))

    # R1/R2 readouts
    for r in decl["readouts"]:
        rid, rm = r.get("id"), r.get("measurement")
        if rm in nodes:
            out.append(_p("R1", "error", f"readout '{rid}' is taken from node '{rm}': a readout must be an independent "
                                         f"measurement, not a model output or candidate state"))
            continue
        if rm not in meas:
            out.append(_p("R1", "error", f"readout '{rid}' names unknown measurement '{rm}'"))
            continue
        feeders = [nid for nid in nodes if rm in ancestors(decl, nid)[0] and nodes[nid].get("role") != "positive_control"]
        if feeders:
            if rm in sysinp:
                pol = decl.get("measurement_sharing", "none")
                if pol not in ("count_split", "declared_shared"):
                    out.append(_p("R2", "error", f"readout '{rid}' uses the system input '{rm}' that also feeds "
                                                 f"{sorted(feeders)}; declare measurement_sharing: count_split (preferred) "
                                                 f"or declared_shared"))
                elif pol == "declared_shared":
                    out.append(_p("R2", "warning", f"readout '{rid}' shares measurement '{rm}' with the candidates "
                                                   f"(declared_shared): counting noise is shared; count splitting would remove it"))
            else:
                out.append(_p("R2", "error", f"readout '{rid}' measurement '{rm}' feeds candidate(s) {sorted(feeders)}: "
                                             f"the observation leaks into the state (e.g. a model that takes measured protein "
                                             f"as input, tested against measured protein)"))

    # C1-C4 chain
    chain = decl["chain"]
    for i, cid in enumerate(chain):
        n = nodes.get(cid)
        if n is None:
            out.append(_p("C1", "error", f"chain node '{cid}' is not declared"))
            continue
        if n.get("role") != "chain":
            out.append(_p("C1", "error", f"node '{cid}' is in the chain but has role '{n.get('role')}'"))
        src_ = n.get("from")
        if i == 0:
            if src_ not in sysinp:
                out.append(_p("C1", "error", f"first chain node '{cid}' must be read from a system input, not '{src_}'"))
        else:
            if src_ != chain[i - 1]:
                out.append(_p("C1", "error", f"chain node '{cid}' is computed from '{src_}', not from the previous stage "
                                             f"'{chain[i - 1]}': only stages on one computation path can be compared"))
            if n.get("operator") != "system":
                out.append(_p("C2", "error", f"chain node '{cid}' uses operator '{n.get('operator')}': stages after the "
                                             f"reference must be the system's own operators (an approximation is a probe, "
                                             f"not a stage)"))
    for nid, n in nodes.items():
        op = str(n.get("op", "")).split(":")[0].lower()
        if op in COMPRESSION_OPS and n.get("operator") != "system":
            out.append(_p("C3", "error", f"node '{nid}' applies analyst compression '{op}': candidates enter the test "
                                         f"with full features; a compression step must not hide inside a comparison"))
        if op in INTERVENTION_OPS and kind == "sufficiency":
            out.append(_p("C4", "error", f"node '{nid}' uses '{op}', which edits the system or its input: that tests a "
                                         f"system intervention, not whether a retained state is sufficient"))
        if n.get("role") not in ("chain", "probe", "positive_control", "technical"):
            out.append(_p("C5", "error", f"node '{nid}' has unknown role '{n.get('role')}'"))
        f = n.get("from")
        for s in (f if isinstance(f, list) else [f]):
            if s not in nodes and s not in meas:
                out.append(_p("C5", "error", f"node '{nid}' reads unknown '{s}'"))
    if not any(n.get("role") == "positive_control" for n in nodes.values()):
        out.append(_p("C5", "warning", "no positive control node: support verdicts will be marked uninformative"))

    # E1 statistics
    stat = decl.get("inference", {}).get("inverse", "crossfit_partial_linear")
    if stat in PARTITION_STATS or stat not in ALLOWED_INVERSE:
        out.append(_p("E1", "error", f"inverse statistic '{stat}' is not allowed: use the continuous conditional test "
                                     f"(crossfit_partial_linear); coarse partitions merge readout differences by construction"))

    # K1/K2 decision and sealing
    d = decl.get("decision", {})
    if "margin_u" not in d or "alpha" not in d:
        out.append(_p("K1", "error", "decision.margin_u and decision.alpha must be declared"))
    if n_treated_clusters is not None and n_treated_clusters < 10:
        out.append(_p("K1", "warning", f"only {n_treated_clusters} treated clusters: cluster-bootstrap intervals are "
                                       f"unreliable; read directions only"))
    if decl.get("status") == "sealed":
        if not decl.get("data_hashes") or not decl.get("code_hash"):
            out.append(_p("K2", "error", "a sealed declaration must carry data_hashes and code_hash"))
    if kind != "sufficiency":
        out.append(_p("K3", "error", f"kind '{kind}' is not handled by run(): it changes the system"))

    if raise_on_error and any(p["severity"] == "error" for p in out):
        raise CheckError(out)
    return out


def errors(problems):
    return [p for p in problems if p["severity"] == "error"]
