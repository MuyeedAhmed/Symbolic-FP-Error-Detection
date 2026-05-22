def _sym_decorrelation(W):
    """Symmetric decorrelation
    i.e. W <- (W * W.T) ^{-1/2} * W
    """
    s, u = linalg.eigh(np.dot(W, W.T))
    # Avoid sqrt of negative values because of rounding errors. Note that
    # np.sqrt(tiny) is larger than tiny and therefore this clipping also
    # prevents division by zero in the next step.
    s = np.clip(s, a_min=np.finfo(W.dtype).tiny, a_max=None)

    # u (resp. s) contains the eigenvectors (resp. square roots of
    # the eigenvalues) of W * W.T
    return np.linalg.multi_dot([u * (1.0 / np.sqrt(s)), u.T, W])