def test_compute_covariance(shape, uniform_weights, csr_container):
    rng = np.random.RandomState(0)
    X = rng.randn(*shape)
    if uniform_weights:
        sw = np.ones(X.shape[0])
    else:
        sw = rng.chisquare(1, shape[0])
    sqrt_sw = np.sqrt(sw)
    X_mean = np.average(X, axis=0, weights=sw)
    X_centered = (X - X_mean) * sqrt_sw[:, None]
    true_covariance = X_centered.T.dot(X_centered)
    X_sparse = csr_container(X * sqrt_sw[:, None])
    gcv = _RidgeGCV(fit_intercept=True)
    computed_cov = gcv._compute_covariance(X_sparse, X_mean, sqrt_sw)
    assert_allclose(true_covariance, computed_cov)