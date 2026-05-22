def test_average_binary_computed_correctly(klass):
    # Checks the SGDClassifier correctly computes the average weights
    eta = 0.1
    alpha = 2.0
    n_samples = 20
    n_features = 10
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)

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

    # simple linear function without noise
    y = np.dot(X, w)
    y = np.sign(y)

    clf.fit(X, y)

    average_weights, average_intercept = asgd(klass, X, y, eta, alpha)
    average_weights = average_weights.reshape(1, -1)
    assert_array_almost_equal(clf.coef_, average_weights, decimal=14)
    assert_almost_equal(clf.intercept_, average_intercept, decimal=14)