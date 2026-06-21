import itertools
import numpy as np
from sklearn.cluster import DBSCAN

from FPReduction import FPReduction
from Z3Solver import prove_arithmetic_divergence
from WitnessBase import sq_dist_matrix, two_cluster_dataset, sqrt_candidates


def _dbscan(D, eps_sq, min_samples):
    return DBSCAN(eps=float(eps_sq), min_samples=min_samples, metric="precomputed").fit_predict(D)

def _n_clusters(labels):
    return len(set(labels[labels >= 0]))

def _derive_witness(d_sq):
    val_l = FPReduction.left_fold(d_sq)
    val_t = FPReduction.tree_fold(d_sq)
    if val_l == val_t:
        return None
    eps_sq = float(min(val_l, val_t))
    hw_in  = "left" if val_l <= val_t else "tree"
    hw_out = "tree" if hw_in == "left" else "left"

    cands = [sqrt_candidates(d_sq[j]) for j in range(len(d_sq))]
    if any(not c for c in cands):
        return None
    for combo in itertools.product(*cands):
        a = np.array(combo, dtype=np.float32)
        if np.all((a * a).astype(np.float32) == d_sq):
            return dict(
                d_sq=d_sq, eps_sq=eps_sq,
                val_left=float(val_l), val_tree=float(val_t),
                p=np.zeros(len(d_sq), np.float32), q=a,
                hw_inside=hw_in, hw_outside=hw_out,
            )
    return None


if __name__ == "__main__":
    print("DBSCAN — float32 fold-order sensitivity\n")
    N_ANCHOR = 5

    raw  = prove_arithmetic_divergence(n=4, hw_a="left", hw_b="tree")
    d_sq = np.abs(raw).astype(np.float32)

    w = _derive_witness(d_sq)
    if w is None:
        raise SystemExit("witness derivation failed")

    print("\nWitness")
    print(f"  d_sq      = {w['d_sq']}")
    print(f"  left_fold = {w['val_left']:.10g}")
    print(f"  tree_fold = {w['val_tree']:.10g}")
    print(f"  eps^2     = {w['eps_sq']:.10g}")
    print(f"  p         = {w['p']}")
    print(f"  q         = {w['q']}")
    print(f"  {w['hw_inside']}_fold places (p,q) on boundary; {w['hw_outside']}_fold places it outside")

    X, eps_sq, min_s = two_cluster_dataset(w["p"], w["q"], w["eps_sq"], N_ANCHOR)

    labels = ["p", "q"] + [f"p-anchor-{i}" for i in range(N_ANCHOR)] + [f"q-anchor-{i}" for i in range(N_ANCHOR)]
    print("\nDataset points")
    for i, pt in enumerate(X):
        d_p = float(np.sum((pt - w["p"]) ** 2))
        d_q = float(np.sum((pt - w["q"]) ** 2))
        print(f"  idx {i:2d}  {labels[i]:12s}  {pt}  dist^2_to_p={d_p:.6g}  dist^2_to_q={d_q:.6g}")

    la = _dbscan(sq_dist_matrix(X, FPReduction.left_fold), eps_sq, min_s)
    lb = _dbscan(sq_dist_matrix(X, FPReduction.tree_fold), eps_sq, min_s)

    print("\nResults")
    print(f"  left_fold (CPU):  {la}  ({_n_clusters(la)} cluster(s))")
    print(f"  tree_fold (GPU):  {lb}  ({_n_clusters(lb)} cluster(s))")
    print(f"  Diverge: {not np.array_equal(la, lb)}")
