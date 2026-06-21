import itertools
import numpy as np
from FPReduction import FPReduction
from Z3Solver import prove_arithmetic_divergence
from WitnessBase import sqrt_candidates


def _derive_witness(d_sq):
    val_l = FPReduction.left_fold(d_sq)
    val_t = FPReduction.tree_fold(d_sq)
    if val_l == val_t:
        return None
    threshold = float(min(val_l, val_t))
    hw_in  = "tree" if val_t <= val_l else "left"
    hw_out = "left" if hw_in == "tree" else "tree"
    cands  = [sqrt_candidates(d_sq[j]) for j in range(len(d_sq))]
    if any(not c for c in cands):
        return None
    for combo in itertools.product(*cands):
        a = np.array(combo, dtype=np.float32)
        if np.all((a * a).astype(np.float32) == d_sq):
            return dict(
                d_sq=d_sq, threshold=threshold,
                val_left=float(val_l), val_tree=float(val_t),
                ref=np.zeros(len(d_sq), np.float32), x=a,
                hw_inside=hw_in, hw_outside=hw_out,
            )
    return None


def _random_subtree(X, rng, depth, max_depth):
    n = len(X)
    if n <= 1 or depth >= max_depth:
        return np.full(n, float(depth))
    f    = rng.integers(X.shape[1])
    vals = X[:, f]
    lo, hi = float(vals.min()), float(vals.max())
    if lo == hi:
        return np.full(n, float(depth))
    t  = rng.uniform(lo, hi)
    lm = vals <= t
    rm = ~lm
    result = np.zeros(n)
    if lm.sum() > 0:
        result[lm] = _random_subtree(X[lm], rng, depth + 1, max_depth)
    if rm.sum() > 0:
        result[rm] = _random_subtree(X[rm], rng, depth + 1, max_depth)
    return result


def run_iforest(X, ref, fold_fn, threshold, n_trees=100, seed=0):
    n         = len(X)
    max_depth = int(np.ceil(np.log2(max(n, 2)))) + 1
    rng       = np.random.default_rng(seed)

    scores    = np.array([fold_fn(((pt.astype(np.float32) - ref) ** 2).astype(np.float32))
                          for pt in X])
    left_mask  = scores <= threshold
    right_mask = ~left_mask
    left_idx   = np.where(left_mask)[0]

    total = np.zeros(n)
    for _ in range(n_trees):
        depths = np.zeros(n)
        depths[right_mask] = 1.0
        if len(left_idx) > 1:
            depths[left_idx] = _random_subtree(X[left_idx], rng, 1, max_depth)
        elif len(left_idx) == 1:
            depths[left_idx] = 1.0
        total += depths

    return total / n_trees


if __name__ == "__main__":
    print("Isolation Forest — float32 fold-order sensitivity\n")

    raw  = prove_arithmetic_divergence(n=4, hw_a="left", hw_b="tree")
    d_sq = np.abs(raw).astype(np.float32)

    w = _derive_witness(d_sq)
    if w is None:
        raise SystemExit("witness derivation failed")

    ref, thr = w["ref"], w["threshold"]

    print("\nWitness")
    print(f"  d_sq      = {w['d_sq']}")
    print(f"  left_fold = {w['val_left']:.10g}")
    print(f"  tree_fold = {w['val_tree']:.10g}")
    print(f"  threshold = {w['threshold']:.10g}  (= min of the two)")
    print(f"  ref       = {ref}")
    print(f"  x         = {w['x']}")
    print(f"  {w['hw_inside']:9s}: score(x) = threshold  -> x joins normal subtree")
    print(f"  {w['hw_outside']:9s}: score(x) = threshold + 1 ULP  -> x isolated at root (path=1)")

    N_NORMAL  = 20
    rng_data  = np.random.default_rng(42)
    normals   = (ref + rng_data.uniform(-0.005, 0.005, (N_NORMAL, len(ref)))).astype(np.float32)
    X         = np.vstack([w["x"], normals]).astype(np.float32)
    pt_labels = ["x"] + [f"normal-{i}" for i in range(N_NORMAL)]

    print(f"\nDataset points ({len(X)} total)")
    for i, pt in enumerate(X):
        sl = float(FPReduction.left_fold(((pt.astype(np.float32) - ref) ** 2).astype(np.float32)))
        st = float(FPReduction.tree_fold(((pt.astype(np.float32) - ref) ** 2).astype(np.float32)))
        print(f"  idx {i:2d}  {pt_labels[i]:10s}  {pt}  score_left={sl:.6g}  score_tree={st:.6g}")

    ANOMALY_THRESHOLD = 1.5
    paths_l = run_iforest(X, ref, FPReduction.left_fold,  thr)
    paths_t = run_iforest(X, ref, FPReduction.tree_fold, thr)
    anom_l  = paths_l < ANOMALY_THRESHOLD
    anom_t  = paths_t < ANOMALY_THRESHOLD

    print(f"\nResults (IsolationForest, {len(X)} points, anomaly_threshold={ANOMALY_THRESHOLD})")
    for i in range(len(X)):
        cl = "ANOMALY" if anom_l[i] else "normal"
        ct = "ANOMALY" if anom_t[i] else "normal"
        print(f"  idx {i:2d}  {pt_labels[i]:10s}  left_fold path={paths_l[i]:.2f} {cl:7s}  "
              f"tree_fold path={paths_t[i]:.2f} {ct}")

    print(f"\n  x: left_fold={'ANOMALY' if anom_l[0] else 'normal'}  "
          f"tree_fold={'ANOMALY' if anom_t[0] else 'normal'}")
    print(f"  Diverge: {bool(anom_l[0] != anom_t[0])}")
