"""Fold assignment: the declared split is kept whenever it can be used, and replaced only when it cannot."""
import numpy as np
import pytest

from anm_diag.stats import fold_is_usable, folds_by_cluster, usable_folds


def _case(n_controls, treated=("T1", "T2"), per_cluster=30):
    cl = np.array([t for t in treated for _ in range(per_cluster)] + [f"c{i}" for i in range(n_controls)])
    A = np.array([1] * (len(treated) * per_cluster) + [0] * n_controls)
    return cl, A


def test_a_usable_split_is_left_exactly_as_declared():
    cl = np.array([f"g{i % 40}" for i in range(400)] + [f"c{i}" for i in range(300)])
    A = np.array([1] * 400 + [0] * 300)
    fold, used = usable_folds(cl, A, 5, 0)
    assert used == 0 and (fold == folds_by_cluster(cl, 5, 0)).all() and fold_is_usable(fold, A, 5)


def test_a_split_with_all_treated_clusters_in_one_fold_is_replaced():
    bad = None
    for n in range(10, 400):                      # with two treated clusters about one split in five is unusable
        cl, A = _case(n)
        if not fold_is_usable(folds_by_cluster(cl, 5, 0), A, 5):
            bad = (cl, A)
            break
    assert bad is not None
    cl, A = bad
    fold, used = usable_folds(cl, A, 5, 0)
    assert used > 0 and fold_is_usable(fold, A, 5) and (fold == folds_by_cluster(cl, 5, used)).all()
    assert all(not fold_is_usable(folds_by_cluster(cl, 5, s), A, 5) for s in range(0, used))


def test_one_treated_cluster_cannot_be_split():
    cl, A = _case(50, treated=("T1",))
    with pytest.raises(ValueError, match="too few treated clusters"):
        usable_folds(cl, A, 5, 0)
