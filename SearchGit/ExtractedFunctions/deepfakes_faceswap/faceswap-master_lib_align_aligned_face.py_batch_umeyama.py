def batch_umeyama(source: np.ndarray, destination: np.ndarray, estimate_scale: bool) -> np.ndarray:
    """A batch implementation to estimate N-D similarity transformation with or without scaling.

    Parameters
    ----------
    source
        (B, M, N) array source coordinates.
    destination
        (M, N) array destination coordinates.
    estimate_scale: bool
        Whether to estimate scaling factor.

    Returns
    -------
    (B, N + 1, N + 1) The homogeneous similarity transformation matrix. The matrix contains NaN
    values only if the problem is not well-conditioned.

    References
    ----------
    .. [1] "Least-squares estimation of transformation parameters between two
            point patterns", Shinji Umeyama, PAMI 1991, :DOI:`10.1109/34.88573`
    """
    # pylint:disable=too-many-locals
    batch_size, num, dim = source.shape  # (B, M, N)

    # Compute mean of source and destination.
    src_mean = source.mean(axis=1)       # (B, N)
    dst_mean = destination.mean(axis=0)  # (N, )

    # Subtract mean from source and destination.
    src_demean = source - src_mean[:, None]  # (B, M, N)
    dst_demean = destination - dst_mean      # (M, N)

    # Eq. (38).
    a = dst_demean.T @ src_demean / num  # (B, N, N)

    # SVD
    u, s, vt = np.linalg.svd(a)

    rot = u @ vt
    det_rot = np.linalg.det(rot)
    # Fix improper rotations
    vt[det_rot < 0, -1, :] *= -1
    rot = u @ vt

    if estimate_scale:
        # Eq. (41) and (42).
        var_src = src_demean.var(axis=1).sum(axis=1)  # (B,)
        scale = s.sum(axis=1) / var_src
    else:
        scale = np.ones(batch_size)

    trans = dst_mean - scale[:, None] * ((rot @ src_mean[..., None])[..., 0])
    retval = np.zeros((batch_size, dim + 1, dim + 1), dtype=source.dtype)
    retval[:, -1, -1] = 1.0

    retval[:, :dim, :dim] = scale[:, None, None] * rot
    retval[:, :dim, dim] = trans
    return retval