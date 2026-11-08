"""Three-way verdicts (paper §2.3, Table 2) and the ordered root-cause reading.

Per candidate x readout, over the union of the declared readout classes:
  support      some class has upper bound of u < margin   (bounded support on the tested comparisons)
  reject       every class has lower bound of u > margin  (rejects the declared readout classes on this state, Eq. 5)
  unresolved   otherwise
Root cause, first match wins:
  no_effect             the readout response Delta is not resolved (interval includes 0): data
  control_failed        the positive control is not supported: the test cannot find information
  relative_loss_at      the reference is not supported, but one chain step adds a resolved loss
                        (lower bound of the change in u > step_margin): a relative, weaker claim
  reference_unresolved  the reference stage (first chain node) is not supported and no step adds a resolved
                        loss: no system-specific claim
  loss_at               a chain stage is rejected after a supported stage: the step between them merged
                        the distinction (it names a step of the map, not a faulty component)
  unresolved            stages after the reference are unresolved without a rejection
  no_deficit            every chain stage is supported
"""


def verdict(interval, margin):
    lo, hi = interval
    if hi < margin:
        return "support"
    if lo > margin:
        return "reject"
    return "unresolved"


def verdict_union(by_class, margin):
    """Verdict for the union of the declared readout classes.

    Supported if some declared class reproduces the response within the margin (its upper bound < margin);
    rejected only if every declared class fails (every lower bound > margin); otherwise unresolved.
    Returns (verdict, chosen class): the supporting class with the lowest upper bound, else the class with
    the lowest point estimate.
    """
    sup = [c for c, x in by_class.items() if x["ci"][1] < margin]
    if sup:
        return "support", min(sup, key=lambda c: by_class[c]["ci"][1])
    if all(x["ci"][0] > margin for x in by_class.values()):
        return "reject", min(by_class, key=lambda c: by_class[c]["u"])
    return "unresolved", min(by_class, key=lambda c: by_class[c]["u"])


def root_cause(delta_ci, pc_verdict, chain, chain_verdicts, step_cis=None, step_margin=None):
    """step_cis: {"a->b": [lo, hi]} for consecutive chain stages (change in u, best class)."""
    lo, hi = delta_ci
    if lo <= 0 <= hi:
        return {"cause": "no_effect", "detail": "readout response not resolved: nothing to explain (data)"}
    if pc_verdict is not None and pc_verdict != "support":
        return {"cause": "control_failed", "detail": "positive control not supported: the test cannot recover "
                                                     "information known to be present; no verdict on the chain"}
    ref = chain[0]
    if chain_verdicts[ref] != "support":
        # absolute sufficiency is unresolved; a resolved additional loss along the chain is still reportable
        if step_cis and step_margin is not None:
            for a, b in zip(chain[:-1], chain[1:]):
                s = step_cis.get(f"{a}->{b}")
                if s is not None and s[0] > step_margin:
                    return {"cause": "relative_loss_at", "from": a, "to": b,
                            "detail": f"reference '{ref}' is {chain_verdicts[ref]} (the readout is not fully explained "
                                      f"even by the input view), but '{a}' -> '{b}' adds a resolved loss "
                                      f"(change in u > {step_margin})"}
        if chain_verdicts[ref] == "reject":
            why = (f"reference stage '{ref}' is rejected: the input view itself does not reproduce this response "
                   f"for any declared readout class, so by these classes the response is not in the input "
                   f"(a limit of the data or the readout classes, not a loss inside the system)")
        else:
            why = (f"reference stage '{ref}' is unresolved: the data cannot tell whether the input view "
                   f"reproduces this response")
        return {"cause": "reference_unresolved", "detail": why + "; no claim about later stages"}
    last_support = ref
    pending = []
    for c in chain[1:]:
        v = chain_verdicts[c]
        if v == "support":
            last_support, pending = c, []
        elif v == "reject":
            span = [last_support] + pending + [c]
            return {"cause": "loss_at", "from": last_support, "to": c, "span": span,
                    "detail": f"distinction merged between '{last_support}' and '{c}'"
                              + (f" (via unresolved {pending})" if pending else "")}
        else:
            pending.append(c)
    if pending:
        return {"cause": "unresolved", "detail": f"stages {pending} unresolved after '{last_support}'"}
    return {"cause": "no_deficit", "detail": "every chain stage supported for this readout"}
