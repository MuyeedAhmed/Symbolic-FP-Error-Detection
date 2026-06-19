import numpy as np
import z3

from FPReduction import FPReduction
from KMeansW_Iter0_Z3 import (
    prove_arithmetic_divergence,
    find_distance_flip,
    blocking_clause,
    margin_constraint,
)


class KMeansWitnessFinder:
    def __init__(self,
                 d: int = 4,
                 hw_a: str = "left",
                 hw_b: str = "tree",
                 n_anchor: int = 10,
                 n_boundary: int = 5,
                 n_iter: int = 50):
        self.d          = d
        self.hw_a       = hw_a
        self.hw_b       = hw_b
        self.fold_a     = FPReduction.get(hw_a)
        self.fold_b     = FPReduction.get(hw_b)
        self.n_anchor   = n_anchor
        self.n_boundary = n_boundary
        self.n_iter     = n_iter
        self._FP32      = z3.Float32()


    def _run_kmeans(self, X: np.ndarray, fold, k: int = 2,
                    init_centroids: np.ndarray | None = None) -> np.ndarray:
        na = self.n_anchor
        if init_centroids is not None:
            centroids = init_centroids.astype(np.float32)
        else:
            centroids = np.vstack([
                X[:na].mean(axis=0),
                X[na:2*na].mean(axis=0),
            ]).astype(np.float32)

        for _ in range(self.n_iter):
            labels = np.zeros(len(X), dtype=int)
            for i, xi in enumerate(X):
                sq = (xi.astype(np.float32) - centroids).astype(np.float32) ** 2
                dists = [fold(sq[k_]) for k_ in range(k)]
                labels[i] = int(np.argmin(dists))
            new_c = np.zeros_like(centroids)
            for ck in range(k):
                m = X[labels == ck]
                new_c[ck] = m.mean(axis=0).astype(np.float32) if len(m) else centroids[ck]
            if np.all(new_c == centroids):
                break
            centroids = new_c

        return labels
