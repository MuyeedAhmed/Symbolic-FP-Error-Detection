def _naive_ledoit_wolf_shrinkage(X):
    # A simple implementation of the formulas from Ledoit & Wolf

    # The computation below achieves the following computations of the
    # "O. Ledoit and M. Wolf, A Well-Conditioned Estimator for
    # Large-Dimensional Covariance Matrices"
    # beta and delta are given in the beginning of section 3.2
    n_samples, n_features = X.shape
    emp_cov = empirical_covariance(X, assume_centered=False)
    mu = np.trace(emp_cov) / n_features
    delta_ = emp_cov.copy()
    delta_.flat[:: n_features + 1] -= mu
    delta = (delta_**2).sum() / n_features
    X2 = X**2
    beta_ = (
        1.0
        / (n_features * n_samples)
        * np.sum(np.dot(X2.T, X2) / n_samples - emp_cov**2)
    )

    beta = min(beta_, delta)
    shrinkage = beta / delta
    return shrinkage