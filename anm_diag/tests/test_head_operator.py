"""The head stages must apply the teddy_mm head exactly as teddy_mm does: unit-length embeddings into the MLP
(teddy_mm scripts/04_train.py line 101, 05_eval.py line 55), then the decoder."""
import numpy as np

from anm_diag.adapters.perturb_cite import PerturbCiteAdapter


def _fake(seed=0, d=8, n_adt=3):
    rng = np.random.default_rng(seed)
    W = {}
    for pre, out in (("", d), ("mu.", n_adt)):
        W[pre + "net.0.weight"] = rng.normal(size=(d, d)); W[pre + "net.0.bias"] = rng.normal(size=d)
        W[pre + "net.3.weight"] = rng.normal(size=(d, d)); W[pre + "net.3.bias"] = rng.normal(size=d)
        W[pre + "net.6.weight"] = rng.normal(size=(out, d)); W[pre + "net.6.bias"] = rng.normal(size=out)
    a = object.__new__(PerturbCiteAdapter)
    z = rng.normal(size=(5, d)) * 26.0
    a._cache = {"head": {"mlp": {k: v for k, v in W.items() if not k.startswith("mu.")},
                         "dec": {k: v for k, v in W.items() if k.startswith("mu.")}},
                "teddy": {"z_last": z}, "teddy_idx": np.arange(5)}
    return a, z, W


def test_head_reads_unit_length_embeddings():
    a, z, W = _fake()
    lat = a._head_latent(None, None, None)
    unit = z / np.linalg.norm(z, axis=1, keepdims=True)
    silu = lambda x: x / (1 + np.exp(-x))
    h = silu(unit @ W["net.0.weight"].T + W["net.0.bias"])
    h = silu(h @ W["net.3.weight"].T + W["net.3.bias"])
    ref = h @ W["net.6.weight"].T + W["net.6.bias"]
    assert np.allclose(lat, ref, atol=1e-4)


def test_head_is_invariant_to_embedding_scale():
    a, z, _ = _fake()
    lat1 = a._head_latent(None, None, None)
    a._cache["teddy"]["z_last"] = z * 0.01
    a._cache.pop("teddy_lat", None)
    lat2 = a._head_latent(None, None, None)
    assert np.allclose(lat1, lat2, atol=1e-4)
    out = a._head_output(None, None, None)
    assert out.shape == (5, 3) and np.isfinite(out).all()
