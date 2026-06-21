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

    @staticmethod
    def _canonical(labels: np.ndarray, X: np.ndarray, k: int = 2) -> np.ndarray:
        means  = [X[labels == c, 0].mean() if (labels == c).any() else np.inf
                  for c in range(k)]
        order  = np.argsort(means)
        mapping = {old: new for new, old in enumerate(order)}
        return np.array([mapping[l] for l in labels])

    def test_propagation(self,
                         witness: dict,
                         anchor_sigma: float = 0.05
                         ) -> dict:
        c0 = witness["c0"]
        c1 = witness["c1"]

        rng = np.random.default_rng(seed=0)   # deterministic, fixed seed

        anchors0 = (c0 + rng.standard_normal((self.n_anchor, self.d)).astype(np.float32)
                    * np.float32(anchor_sigma))
        anchors1 = (c1 + rng.standard_normal((self.n_anchor, self.d)).astype(np.float32)
                    * np.float32(anchor_sigma))

        
        x_w    = witness["x"].astype(np.float32)
        normal = np.float32(c1 - c0)
        nn     = np.float32(np.linalg.norm(normal))
        if nn == 0:
            boundary = np.tile(x_w, (self.n_boundary, 1)).astype(np.float32)
        else:
            scales   = np.linspace(-1e-4, 1e-4, self.n_boundary, dtype=np.float32)
            boundary = np.array([x_w + s * normal / nn for s in scales],
                                dtype=np.float32)

        X = np.vstack([anchors0, anchors1, boundary]).astype(np.float32)

        init_c = np.vstack([c0, c1]).astype(np.float32)
        lab_A = self._run_kmeans(X, self.fold_a, init_centroids=init_c)
        lab_B = self._run_kmeans(X, self.fold_b, init_centroids=init_c)

        cA = self._canonical(lab_A, X)
        cB = self._canonical(lab_B, X)

        diverging_idx = np.where(cA != cB)[0]
        diverges      = len(diverging_idx) > 0

        return {
            "X":            X,
            "labels_A":     cA,
            "labels_B":     cB,
            "diverges":     diverges,
            "n_diverging":  len(diverging_idx),
            "diverging_idx": diverging_idx,
        }


    def find_witness(self,
                     max_attempts: int = 20,
                     timeout_ms: int = 30_000,
                     verbose: bool = True
                     ) -> dict | None:
        extra   = []
        min_margin = 0.0
        last_witness = None

        d_sq_base = None
        c1 = prove_arithmetic_divergence(n=self.d, hw_a=self.hw_a, hw_b=self.hw_b, timeout_ms=timeout_ms)
        if c1 is not None:
            d_sq_base = np.abs(c1).astype(np.float32)
        if verbose and d_sq_base is not None:
            print(f"  Phase-1 anchor: d_sq = {d_sq_base}")

        for attempt in range(max_attempts):
            if verbose:
                print(f"\n[attempt {attempt}] Z3 search "
                      f"(min_margin={min_margin:.2e}, "
                      f"extra_constraints={len(extra)}) …")

            w = find_distance_flip(
                d=self.d,
                hw_a=self.hw_a,
                hw_b=self.hw_b,
                timeout_ms=timeout_ms,
                extra_constraints=extra,
                d_sq_base=d_sq_base,
            )

            if w is None:
                if verbose:
                    print("  Z3 returned no result (timeout or UNSAT) — stopping.")
                break

            if verbose:
                self._print_flip(w, attempt)

            last_witness = w

            assert w["label_A"] != w["label_B"]

            prop = self.test_propagation(w, anchor_sigma=0.0)
            last_witness["propagation"] = prop

            if prop["diverges"]:
                if verbose:
                    print(f"  Propagates! {prop['n_diverging']} point(s) differ "
                          f"in final KMeans labels.")
                return last_witness

            if verbose:
                print(f"  No cascade ({prop['n_diverging']} diverging points). "
                      f"Tightening constraints …")

            extra.append(blocking_clause(w, self._FP32))
            min_margin = max(min_margin * 2, float(max(w["margin_A"], w["margin_B"])) * 2)
            if min_margin > 0:
                extra.append(margin_constraint(min_margin, w, self._FP32))

        if verbose and last_witness is not None:
            print(f"\n[find_witness] Budget exhausted.  "
                  f"Returning last non-propagating witness.")
        return last_witness

    def _print_flip(self, w: dict, attempt: int) -> None:
        print(f"  Z3 witness found (attempt {attempt}):")
        print(f"    a  = {w['a']}")
        print(f"    b  = {w['b']}")
        print(f"    x  = {w['x']}")
        print(f"    c0 = {w['c0']}")
        print(f"    c1 = {w['c1']}")
        print(f"    dist(x,c0)  {self.hw_a}={w['d0_A']:.10f}  "
              f"{self.hw_b}={w['d0_B']:.10f}")
        print(f"    dist(x,c1)  {self.hw_a}={w['d1_A']:.10f}  "
              f"{self.hw_b}={w['d1_B']:.10f}")
        print(f"    HW_A ({self.hw_a}) → cluster {w['label_A']}")
        print(f"    HW_B ({self.hw_b}) → cluster {w['label_B']}  ← FLIP")
        print(f"    margin_A={w['margin_A']:.3e}  margin_B={w['margin_B']:.3e}")


if __name__ == "__main__":
    print("=" * 65)
    print("KMeans FP tiling witness  —  systematic Z3 + propagation")
    print("=" * 65)

    finder = KMeansWitnessFinder(
        d=4,
        hw_a="left",
        hw_b="tree",
        n_anchor=10,
        n_boundary=5,
        n_iter=50,
    )

    result = finder.find_witness(max_attempts=20, timeout_ms=60_000)

    if result is None:
        print("\nNo witness found.")
    else:
        prop = result.get("propagation", {})
        print("\n" + "=" * 65)
        print("FINAL WITNESS SUMMARY")
        print("=" * 65)
        print(f"  d={finder.d}  HW_A={finder.hw_a}  HW_B={finder.hw_b}")
        print(f"  x  = {result['x']}")
        print(f"  c0 = {result['c0']}")
        print(f"  c1 = {result['c1']}")
        print(f"  Single-point assignment reversal: True")
        if prop:
            print(f"  Propagation to full KMeans:       "
                  f"{'True' if prop['diverges'] else 'False'}")
            if prop["diverges"]:
                print(f"  Points with different final labels: "
                      f"{prop['n_diverging']} / {len(prop['X'])}")
