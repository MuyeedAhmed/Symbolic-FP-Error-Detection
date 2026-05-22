def test_elasticnet_precompute_gram_weighted_samples():
    # check the equivalence between passing a precomputed Gram matrix and
    # internal computation using sample weights.
    X, y, _, _ = build_dataset()

    rng = np.random.RandomState(0)
    sample_weight = rng.lognormal(size=y.shape)

    w_norm = sample_weight * (y.shape / np.sum(sample_weight))
    X_c = X - np.average(X, axis=0, weights=w_norm)
    X_r = X_c * np.sqrt(w_norm)[:, np.newaxis]
    gram = np.dot(X_r.T, X_r)

    clf1 = ElasticNet(alpha=0.01, precompute=gram)
    clf1.fit(X_c, y, sample_weight=sample_weight)

    clf2 = ElasticNet(alpha=0.01, precompute=False)
    clf2.fit(X, y, sample_weight=sample_weight)

    assert_allclose(clf1.coef_, clf2.coef_)