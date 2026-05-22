def _estimate_gaussian_covariances_diag(resp, X, nk, means, reg_covar, xp=None):
    """Estimate the diagonal covariance vectors.

    Parameters
    ----------
    responsibilities : array-like of shape (n_samples, n_components)

    X : array-like of shape (n_samples, n_features)

    nk : array-like of shape (n_components,)

    means : array-like of shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components, n_features)
        The covariance vector of the current components.
    """
    xp, _ = get_namespace(X, xp=xp)
    avg_X2 = (resp.T @ (X * X)) / nk[:, xp.newaxis]
    avg_means2 = means**2
    return avg_X2 - avg_means2 + reg_covar