def _taubin_smooth(verts, faces, iters, lam=0.5, mu=-0.53):
    # Taubin lambda|mu smoothing: low-pass the mesh surface without the shrinkage of a Laplacian blur
    # (the mu inflation pass cancels the lambda pass's volume loss). Uniform (umbrella) weights.
    if iters <= 0 or len(verts) == 0 or len(faces) == 0:
        return verts
    nv = len(verts)
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [0, 2]]], 0)
    e = np.concatenate([e, e[:, ::-1]], 0)  # symmetric adjacency
    adj = coo_matrix((np.ones(len(e), np.float32), (e[:, 0], e[:, 1])), shape=(nv, nv)).tocsr()
    adj.data[:] = 1.0
    deg = np.clip(np.asarray(adj.sum(1)).ravel(), 1.0, None).astype(np.float32)[:, None]
    v = verts.astype(np.float32)  # fp32 matvec: ~2x faster, sub-micron drift on unit-scale verts
    for _ in range(int(iters)):
        for fac in (lam, mu):
            v = v + np.float32(fac) * ((adj @ v) / deg - v)  # fac * (mean(neighbours) - v)
    return np.ascontiguousarray(v)