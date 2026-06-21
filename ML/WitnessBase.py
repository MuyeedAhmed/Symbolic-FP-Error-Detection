import numpy as np
from FPReduction import FPReduction


def sq_dist_matrix(X, fold_fn):
    X = X.astype(np.float32)
    n = len(X)
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            sq = ((X[i] - X[j]) ** 2).astype(np.float32)
            D[i, j] = D[j, i] = float(fold_fn(sq))
    return D


def two_cluster_dataset(p, q, eps_sq, n_anchor=5, seed=0):
    p, q = p.astype(np.float32), q.astype(np.float32)
    eps    = float(np.sqrt(float(eps_sq)))
    pq_dir = (q - p) / np.linalg.norm(q - p)
    rng    = np.random.default_rng(seed)

    def _anchors(center, away_dir, n):
        d     = len(p)
        noise = rng.standard_normal((n, d)).astype(np.float32)
        noise -= (noise @ away_dir[:, None]) * away_dir
        noise += away_dir
        noise /= np.linalg.norm(noise, axis=1, keepdims=True)
        return (center + noise * np.float32(eps / 1000)).astype(np.float32)

    X = np.vstack([
        p, q,
        _anchors(p, -pq_dir, n_anchor),
        _anchors(q, +pq_dir, n_anchor),
    ]).astype(np.float32)
    return X, float(eps_sq), n_anchor + 1


def sqrt_candidates(target, n_ulp=2):
    if target <= 0:
        return [np.float32(0.0)] if target == 0 else []
    base = FPReduction.f32_to_bits(np.float32(np.sqrt(float(target))))
    return [
        a for delta in range(-n_ulp, n_ulp + 1)
        if np.isfinite(a := FPReduction.bits_to_f32(base + delta))
        and np.float32(a * a) == target
    ]
