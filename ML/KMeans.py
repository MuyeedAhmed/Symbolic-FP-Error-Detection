import numpy as np
import z3

from FPReduction import FPReduction
from Z3Solver import (
    prove_arithmetic_divergence,
    find_distance_flip,
    blocking_clause,
    margin_constraint,
)


class KMeansWitnessFinder:
    def __init__(self, d=4, hw_a="left", hw_b="tree",
                 n_anchor=10, n_boundary=5, n_iter=50):
        self.d          = d
        self.hw_a       = hw_a
        self.hw_b       = hw_b
        self.fold_a     = FPReduction.get(hw_a)
        self.fold_b     = FPReduction.get(hw_b)
        self.n_anchor   = n_anchor
        self.n_boundary = n_boundary
        self.n_iter     = n_iter
        self._FP32      = z3.Float32()

    def _run_kmeans(self, X, fold, k=2, init_centroids=None):
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
                sq    = (xi.astype(np.float32) - centroids).astype(np.float32) ** 2
                dists = [fold(sq[k_]) for k_ in range(k)]
                labels[i] = int(np.argmin(dists))
            new_c = np.zeros_like(centroids)
            for ck in range(k):
                m = X[labels == ck]
                new_c[ck] = (m.mean(axis=0).astype(np.float32)
                             if len(m) else centroids[ck])
            if np.all(new_c == centroids):
                break
            centroids = new_c
        return labels

    @staticmethod
    def _canonical(labels, X, k=2):
        means   = [X[labels == c, 0].mean() if (labels == c).any() else np.inf for c in range(k)]
        order   = np.argsort(means)
        mapping = {old: new for new, old in enumerate(order)}
        return np.array([mapping[l] for l in labels])

    def test_propagation(self, witness, anchor_sigma=0.05):
        c0  = witness["c0"]
        c1  = witness["c1"]
        x_w = witness["x"].astype(np.float32)
        rng = np.random.default_rng(seed=0)

        anchors0 = (c0 + rng.standard_normal((self.n_anchor, self.d)).astype(np.float32) * np.float32(anchor_sigma))
        anchors1 = (c1 + rng.standard_normal((self.n_anchor, self.d)).astype(np.float32) * np.float32(anchor_sigma))

        normal = np.float32(c1 - c0)
        nn     = np.float32(np.linalg.norm(normal))
        if nn == 0:
            boundary = np.tile(x_w, (self.n_boundary, 1)).astype(np.float32)
        else:
            scales   = np.linspace(-1e-4, 1e-4, self.n_boundary, dtype=np.float32)
            boundary = np.array([x_w + s * normal / nn for s in scales],
                                dtype=np.float32)

        X      = np.vstack([anchors0, anchors1, boundary]).astype(np.float32)
        init_c = np.vstack([c0, c1]).astype(np.float32)
        lab_A  = self._run_kmeans(X, self.fold_a, init_centroids=init_c)
        lab_B  = self._run_kmeans(X, self.fold_b, init_centroids=init_c)
        cA     = self._canonical(lab_A, X)
        cB     = self._canonical(lab_B, X)
        diff   = np.where(cA != cB)[0]

        return {"X": X, "labels_A": cA, "labels_B": cB,
                "diverges": len(diff) > 0, "n_diverging": len(diff),
                "diverging_idx": diff}

    def find_witness(self, max_attempts=20, timeout_ms=30_000, verbose=True):
        extra      = []
        min_margin = 0.0
        last_w     = None

        d_sq_base = None
        c1 = prove_arithmetic_divergence(n=self.d, hw_a=self.hw_a, hw_b=self.hw_b, timeout_ms=timeout_ms)
        if c1 is not None:
            d_sq_base = np.abs(c1).astype(np.float32)

        for attempt in range(max_attempts):
            w = find_distance_flip(
                d=self.d, hw_a=self.hw_a, hw_b=self.hw_b,
                timeout_ms=timeout_ms, extra_constraints=extra,
                d_sq_base=d_sq_base,
            )
            if w is None:
                break

            if verbose:
                print(f"  x  = {w['x']}")
                print(f"  c0 = {w['c0']}")
                print(f"  c1 = {w['c1']}")
                print(f"  {self.hw_a}: dist(x,c0)={w['d0_A']:.10g}  dist(x,c1)={w['d1_A']:.10g}  x->cluster {w['label_A']}")
                print(f"  {self.hw_b}: dist(x,c0)={w['d0_B']:.10g}  dist(x,c1)={w['d1_B']:.10g}  x->cluster {w['label_B']}")

            last_w = w
            assert w["label_A"] != w["label_B"]

            prop = self.test_propagation(w, anchor_sigma=0.0)
            last_w["propagation"] = prop

            if prop["diverges"]:
                return last_w


            extra.append(blocking_clause(w, self._FP32))
            min_margin = max(min_margin * 2, float(max(w["margin_A"], w["margin_B"])) * 2)
            if min_margin > 0:
                extra.append(margin_constraint(min_margin, w, self._FP32))

        return last_w


if __name__ == "__main__":
    print("KMeans — float32 fold-order sensitivity\n")

    finder = KMeansWitnessFinder(d=4, hw_a="left", hw_b="tree", n_anchor=10, n_boundary=5, n_iter=50)
    result = finder.find_witness(max_attempts=20, timeout_ms=60_000)

    if result is None:
        raise SystemExit("No witness found.")

    print("\nWitness")
    print(f"  x  = {result['x']}")
    print(f"  c0 = {result['c0']}")
    print(f"  c1 = {result['c1']}")
    print(f"  left_fold: dist(x,c0)={result['d0_A']:.10g}  dist(x,c1)={result['d1_A']:.10g}  x -> cluster {result['label_A']}")
    print(f"  tree_fold: dist(x,c0)={result['d0_B']:.10g}  dist(x,c1)={result['d1_B']:.10g}  x -> cluster {result['label_B']}")

    prop = result.get("propagation", {})
    if prop:
        n = len(prop["X"])
        print(f"\nResults (iterated KMeans, {n} points)")
        print(f"  left_fold: {prop['labels_A']}")
        print(f"  tree_fold: {prop['labels_B']}")
        print(f"  Diverging: {prop['n_diverging']}/{n}")

