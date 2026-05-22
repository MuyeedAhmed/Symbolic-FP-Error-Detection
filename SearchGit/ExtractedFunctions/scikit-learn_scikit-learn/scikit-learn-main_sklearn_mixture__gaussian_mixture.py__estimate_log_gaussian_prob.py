def _estimate_log_gaussian_prob(X, means, precisions_chol, covariance_type, xp=None):
    """Estimate the log Gaussian probability.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)

    means : array-like of shape (n_components, n_features)

    precisions_chol : array-like
        Cholesky decompositions of the precision matrices.
        'full' : shape of (n_components, n_features, n_features)
        'tied' : shape of (n_features, n_features)
        'diag' : shape of (n_components, n_features)
        'spherical' : shape of (n_components,)

    covariance_type : {'full', 'tied', 'diag', 'spherical'}

    Returns
    -------
    log_prob : array, shape (n_samples, n_components)
    """
    xp, _, device_ = get_namespace_and_device(X, means, precisions_chol, xp=xp)
    n_samples, n_features = X.shape
    n_components, _ = means.shape
    # The determinant of the precision matrix from the Cholesky decomposition
    # corresponds to the negative half of the determinant of the full precision
    # matrix.
    # In short: det(precision_chol) = - det(precision) / 2
    log_det = _compute_log_det_cholesky(precisions_chol, covariance_type, n_features)

    if covariance_type == "full":
        log_prob = xp.empty((n_samples, n_components), dtype=X.dtype, device=device_)
        for k in range(means.shape[0]):
            mu = means[k, :]
            prec_chol = precisions_chol[k, :, :]
            y = (X @ prec_chol) - (mu @ prec_chol)
            log_prob[:, k] = xp.sum(xp.square(y), axis=1)

    elif covariance_type == "tied":
        log_prob = xp.empty((n_samples, n_components), dtype=X.dtype, device=device_)
        for k in range(means.shape[0]):
            mu = means[k, :]
            y = (X @ precisions_chol) - (mu @ precisions_chol)
            log_prob[:, k] = xp.sum(xp.square(y), axis=1)

    elif covariance_type == "diag":
        precisions = precisions_chol**2
        log_prob = (
            xp.sum((means**2 * precisions), axis=1)
            - 2.0 * (X @ (means * precisions).T)
            + (X**2 @ precisions.T)
        )

    elif covariance_type == "spherical":
        precisions = precisions_chol**2
        log_prob = (
            xp.sum(means**2, axis=1) * precisions
            - 2 * (X @ means.T * precisions)
            + xp.linalg.outer(row_norms(X, squared=True), precisions)
        )
    # Since we are using the precision of the Cholesky decomposition,
    # `- 0.5 * log_det_precision` becomes `+ log_det_precision_chol`
    return -0.5 * (n_features * math.log(2 * math.pi) + log_prob) + log_det