def test_elasticnet_precompute_incorrect_gram():
    # check that passing an invalid precomputed Gram matrix will raise an
    # error.
    X, y, _, _ = build_dataset()

    rng = np.random.RandomState(0)

    X_centered = X - np.average(X, axis=0)
    garbage = rng.standard_normal(X.shape)
    precompute = np.dot(garbage.T, garbage)

    clf = ElasticNet(alpha=0.01, precompute=precompute)
    msg = "Gram matrix.*did not pass validation.*"
    with pytest.raises(ValueError, match=msg):
        clf.fit(X_centered, y)