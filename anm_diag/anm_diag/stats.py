"""Estimators: forward AUC and the inverse unexplained fraction u, with cluster bootstrap.

Inverse (paper §2.3, stochastic case, mean statistic): with A the source indicator, O the readout and
z a candidate state, cross-fitted residuals give
    theta(z) = sum_i (O_i - m(z_i)) (A_i - e(z_i)) / sum_i (A_i - e(z_i))^2,   u(z) = theta(z) / Delta,
where Delta is the same statistic with no state (the treated-minus-control mean difference). u = 0 means
the readout response is fully accounted for by a readout depending only on z (within the declared
readout class); u = 1 means z accounts for none of it. Folds are grouped by cluster (treated target);
intervals come from a cluster bootstrap (treated clusters resampled whole, controls resampled by cell).
"""
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import roc_auc_score


def _hgb_extra(spec):
    """Optional boosting settings a declaration may fix (e.g. max_features for wide input views)."""
    return {k: spec[k] for k in ("max_features", "l2_regularization", "min_samples_leaf") if k in spec}


def make_learner(spec, kind):
    if kind == "m":
        if spec["regressor"] == "RidgeCV":
            return RidgeCV(alphas=np.asarray(spec["alphas"], float))
        if spec["regressor"] == "HistGradientBoostingRegressor":
            return HistGradientBoostingRegressor(max_iter=spec["max_iter"], learning_rate=spec["learning_rate"],
                                                 max_leaf_nodes=spec["max_leaf_nodes"], random_state=0, **_hgb_extra(spec))
    else:
        if spec["classifier"] == "LogisticRegression":
            return LogisticRegression(C=spec["C"], max_iter=2000)
        if spec["classifier"] == "HistGradientBoostingClassifier":
            return HistGradientBoostingClassifier(max_iter=spec["max_iter"], learning_rate=spec["learning_rate"],
                                                  max_leaf_nodes=spec["max_leaf_nodes"], random_state=0, **_hgb_extra(spec))
    raise ValueError(f"unknown learner spec {spec}")


def folds_by_cluster(clusters, k, seed):
    ug = np.unique(clusters)
    rng = np.random.default_rng(seed)
    rng.shuffle(ug)
    fo = {g: i % k for i, g in enumerate(ug)}
    return np.array([fo[g] for g in clusters])


def fold_is_usable(fold, A, k):
    """Every training part (all folds but one) must hold both treated and control cells."""
    A = np.asarray(A)
    return all(len(np.unique(A[fold != f])) == 2 for f in range(k))


def usable_folds(clusters, A, k, seed, max_tries=100):
    """Folds by cluster with the declared seed; if that split leaves a training part without treated cells
    (or without controls), the next seeds are tried in order. Returns (fold, seed used).

    folds_by_cluster deals whole clusters to folds at random. With few treated clusters they can all land in
    one fold, and the model for every other fold is then trained on controls only. The fallback reads cluster
    labels, the treatment indicator and the seed, never an outcome, and it leaves every split that was usable
    exactly as it was."""
    for s in range(int(seed), int(seed) + max_tries):
        fold = folds_by_cluster(clusters, k, s)
        if fold_is_usable(fold, A, k):
            return fold, s
    raise ValueError(f"no usable fold split within {max_tries} seeds: too few treated clusters for {k} folds "
                     f"(cluster by a finer unit, e.g. source.cluster: guide)")


def cluster_bootstrap(A, clusters, B, seed):
    gx = np.unique(clusters[A == 1])
    gi = {g: np.where(clusters == g)[0] for g in gx}
    ic = np.where(A == 0)[0]
    rng = np.random.default_rng(seed + 1)
    return [np.concatenate([gi[g] for g in rng.choice(gx, len(gx))] + [rng.choice(ic, len(ic))]) for _ in range(B)]


def theta(ry, ra, s=slice(None)):
    return float((ry[s] * ra[s]).sum() / (ra[s] * ra[s]).sum())


class Crossfit:
    """Caches cross-fitted residuals per (state, class)."""

    def __init__(self, A, fold, classes):
        self.A, self.fold, self.classes = A, fold, classes
        self.k = int(fold.max()) + 1
        self._ra = {}

    def ra(self, sid, S, cls):
        key = (sid, cls)
        if key not in self._ra:
            out = np.zeros(len(self.A))
            for f in range(self.k):
                tr, te = self.fold != f, self.fold == f
                mdl = make_learner(self.classes[cls], "e").fit(S[tr], self.A[tr])
                out[te] = self.A[te] - mdl.predict_proba(S[te])[:, 1]
            self._ra[key] = out
        return self._ra[key]

    def ry(self, S, Y, cls):
        out = np.zeros(len(Y))
        for f in range(self.k):
            tr, te = self.fold != f, self.fold == f
            mdl = make_learner(self.classes[cls], "m").fit(S[tr], Y[tr])
            out[te] = Y[te] - mdl.predict(S[te])
        return out


def forward_auc(A, ra, boots):
    score = A - ra
    point = float(roc_auc_score(A, score))
    bs = np.array([roc_auc_score(A[s], score[s]) for s in boots])
    return point, bs


def ci(x, alpha):
    return [float(np.percentile(x, 100 * alpha / 2)), float(np.percentile(x, 100 * (1 - alpha / 2)))]
