def _estimate_gaussian_covariances_full(resp, X, nk, means, reg_covar, xp=None):
    """Estimate the full covariance matrices.

    Parameters
    ----------
    resp : array-like of shape (n_samples, n_components)

    X : array-like of shape (n_samples, n_features)

    nk : array-like of shape (n_components,)

    means : array-like of shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components, n_features, n_features)
        The covariance matrix of the current components.
    """
    xp, _, device_ = get_namespace_and_device(X, xp=xp)
    n_components, n_features = means.shape
    covariances = xp.empty(
        (n_components, n_features, n_features), device=device_, dtype=X.dtype
    )
    for k in range(n_components):
        diff = X - means[k, :]
        covariances[k, :, :] = ((resp[:, k] * diff.T) @ diff) / nk[k]
        _add_to_diagonal(covariances[k, :, :], reg_covar, xp)
    return covariances