def test_sgd_averaged_partial_fit(klass):
    # Tests whether the partial fit yields the same average as the fit
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

    clf.partial_fit(X[: int(n_samples / 2)][:], y[: int(n_samples / 2)])
    clf.partial_fit(X[int(n_samples / 2) :][:], y[int(n_samples / 2) :])
    average_weights, average_intercept = asgd(klass, X, y, eta, alpha)

    assert_array_almost_equal(clf.coef_, average_weights, decimal=16)
    assert_almost_equal(clf.intercept_[0], average_intercept, decimal=16)