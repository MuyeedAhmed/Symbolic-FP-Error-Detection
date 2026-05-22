def _estimate_gaussian_covariances_tied(resp, X, nk, means, reg_covar, xp=None):
    """Estimate the tied covariance matrix.

    Parameters
    ----------
    resp : array-like of shape (n_samples, n_components)

    X : array-like of shape (n_samples, n_features)

    nk : array-like of shape (n_components,)

    means : array-like of shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariance : array, shape (n_features, n_features)
        The tied covariance matrix of the components.
    """
    xp, _ = get_namespace(X, means, xp=xp)
    avg_X2 = X.T @ X
    avg_means2 = nk * means.T @ means
    covariance = avg_X2 - avg_means2
    covariance /= xp.sum(nk)
    _add_to_diagonal(covariance, reg_covar, xp)
    return covariance