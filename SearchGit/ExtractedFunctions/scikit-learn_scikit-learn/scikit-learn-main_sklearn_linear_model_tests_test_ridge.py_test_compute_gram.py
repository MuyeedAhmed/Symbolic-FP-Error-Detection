def test_compute_gram(shape, uniform_weights, csr_container):
    rng = np.random.RandomState(0)
    X = rng.randn(*shape)
    if uniform_weights:
        sw = np.ones(X.shape[0])
    else:
        sw = rng.chisquare(1, shape[0])
    sqrt_sw = np.sqrt(sw)
    X_mean = np.average(X, axis=0, weights=sw)
    X_centered = (X - X_mean) * sqrt_sw[:, None]
    true_gram = X_centered.dot(X_centered.T)
    X_sparse = csr_container(X * sqrt_sw[:, None])
    gcv = _RidgeGCV(fit_intercept=True)
    computed_gram = gcv._compute_gram(X_sparse, X_mean, sqrt_sw)
    assert_allclose(true_gram, computed_gram)