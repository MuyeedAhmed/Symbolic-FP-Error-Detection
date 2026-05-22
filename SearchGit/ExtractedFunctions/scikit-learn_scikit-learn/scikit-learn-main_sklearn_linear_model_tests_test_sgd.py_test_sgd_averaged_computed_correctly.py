def test_sgd_averaged_computed_correctly(klass):
    # Tests the average regressor matches the naive implementation

    eta = 0.001
    alpha = 0.01
    n_samples = 20
    n_features = 10
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)

    # simple linear function without noise
    y = np.dot(X, w)

    clf = klass(
        loss="squared_error",
        learning_rate="constant",
        eta0=eta,
        alpha=alpha,
        fit_intercept=True,
        max_iter=1,
        average=True,
        shuffle=False,
    )

    clf.fit(X, y)
    average_weights, average_intercept = asgd(klass, X, y, eta, alpha)

    assert_array_almost_equal(clf.coef_, average_weights, decimal=16)
    assert_almost_equal(clf.intercept_, average_intercept, decimal=16)