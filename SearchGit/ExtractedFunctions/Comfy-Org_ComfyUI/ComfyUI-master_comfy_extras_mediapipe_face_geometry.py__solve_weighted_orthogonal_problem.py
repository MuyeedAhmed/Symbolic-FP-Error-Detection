def _solve_weighted_orthogonal_problem(src: np.ndarray, tgt: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted orthogonal Procrustes (similarity). Returns 4x4 M with
    `target ≈ M @ homogeneous(source)` in the weighted LS sense. fp64 for
    SVD stability. Port of procrustes_solver.cc."""
    sqrt_w = np.sqrt(weights.astype(np.float64))
    w_total = float((sqrt_w ** 2).sum())
    ws = src.astype(np.float64) * sqrt_w
    wt = tgt.astype(np.float64) * sqrt_w

    c_w = (ws @ sqrt_w) / w_total
    centered = ws - np.outer(c_w, sqrt_w)
    U, _S, Vt = np.linalg.svd(wt @ centered.T, full_matrices=True)
    # Disallow reflection: flip the least-significant axis when det(U)·det(V)<0.
    post, pre = U.copy(), Vt.T.copy()
    if np.linalg.det(post) * np.linalg.det(pre) < 0:
        post[:, 2] *= -1.0
    R = post @ pre.T

    denom = float((centered * ws).sum())
    if denom < 1e-12:
        raise ValueError("Procrustes denominator collapsed (degenerate source).")
    scale = float((R @ centered * wt).sum()) / denom
    translation = ((wt - scale * (R @ ws)) @ sqrt_w) / w_total

    M = np.eye(4, dtype=np.float64)
    M[:3, :3] = scale * R
    M[:3, 3] = translation
    return M