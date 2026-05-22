def test_elasticnet_precompute_gram():
    # Check the dtype-aware check for a precomputed Gram matrix
    # (see https://github.com/scikit-learn/scikit-learn/pull/22059
    # and https://github.com/scikit-learn/scikit-learn/issues/21997).
    # Here: (X_c.T, X_c)[2, 3] is not equal to np.dot(X_c[:, 2], X_c[:, 3])
    # but within tolerance for np.float32

    rng = np.random.RandomState(58)
    X = rng.binomial(1, 0.25, (1000, 4)).astype(np.float32)
    y = rng.rand(1000).astype(np.float32)

    X_c = X - np.average(X, axis=0)
    gram = np.dot(X_c.T, X_c)

    clf1 = ElasticNet(alpha=0.01, precompute=gram)
    clf1.fit(X_c, y)

    clf2 = ElasticNet(alpha=0.01, precompute=False)
    clf2.fit(X, y)

    assert_allclose(clf1.coef_, clf2.coef_)