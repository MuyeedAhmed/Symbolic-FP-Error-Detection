def _umeyama(source: np.ndarray, destination: np.ndarray, estimate_scale: bool) -> np.ndarray:
    """Estimate N-D similarity transformation with or without scaling.

    Imported, and slightly adapted, directly from:
    https://github.com/scikit-image/scikit-image/blob/master/skimage/transform/_geometric.py


    Parameters
    ----------
    source
        (M, N) array source coordinates.
    destination
        (M, N) array destination coordinates.
    estimate_scale
        Whether to estimate scaling factor.

    Returns
    -------
    (N + 1, N + 1) The homogeneous similarity transformation matrix. The matrix contains NaN values
    only if the problem is not well-conditioned.

    References
    ----------
    .. [1] "Least-squares estimation of transformation parameters between two
            point patterns", Shinji Umeyama, PAMI 1991, :DOI:`10.1109/34.88573`
    """
    # pylint:disable=invalid-name,too-many-locals
    num = source.shape[0]
    dim = source.shape[1]

    # Compute mean of source and destination.
    src_mean = source.mean(axis=0)
    dst_mean = destination.mean(axis=0)

    # Subtract mean from source and destination.
    src_demean = source - src_mean
    dst_demean = destination - dst_mean

    # Eq. (38).
    A = dst_demean.T @ src_demean / num

    # Eq. (39).
    d = np.ones((dim,), dtype=np.double)
    if np.linalg.det(A) < 0:
        d[dim - 1] = -1

    retval = np.eye(dim + 1, dtype=np.double)

    U, S, V = np.linalg.svd(A)

    # Eq. (40) and (43).
    rank = np.linalg.matrix_rank(A)
    if rank == 0:
        return np.nan * retval
    if rank == dim - 1:
        if np.linalg.det(U) * np.linalg.det(V) > 0:
            retval[:dim, :dim] = U @ V
        else:
            s = d[dim - 1]
            d[dim - 1] = -1
            retval[:dim, :dim] = U @ np.diag(d) @ V
            d[dim - 1] = s
    else:
        retval[:dim, :dim] = U @ np.diag(d) @ V

    if estimate_scale:
        # Eq. (41) and (42).
        scale = 1.0 / src_demean.var(axis=0).sum() * (S @ d)
    else:
        scale = 1.0

    retval[:dim, dim] = dst_mean - scale * (retval[:dim, :dim] @ src_mean.T)
    retval[:dim, :dim] *= scale

    return retval