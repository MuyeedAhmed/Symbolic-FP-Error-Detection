def test_enet_alpha_max(
    estimatorCV, estimator, X_is_sparse, fit_intercept, positive, sample_weight
):
    X = np.array([[3.0, -1.0], [2.0, -5.0], [5.0, -3.0], [1.0, -4.0]])
    beta = np.array([1, -2])
    y = X @ beta
    params = dict(fit_intercept=fit_intercept, positive=positive)
    if issubclass(estimator, MultiTaskElasticNet):
        n_tasks = 3
        y = np.tile(y[:, None], reps=(1, n_tasks))
        params.pop("positive")
        if positive:
            return

    if X_is_sparse:
        X = X_is_sparse(X)
    # Test alpha_max makes coefs zero.
    reg = estimatorCV(alphas=1, cv=2, eps=1, **params)
    reg.fit(X, y, sample_weight=sample_weight)
    assert_allclose(reg.coef_, 0, atol=1e-5)
    alpha_max = reg.alpha_
    # Test smaller alpha makes coefs nonzero.
    reg = estimator(alpha=0.99 * alpha_max, tol=1e-8, **params)
    reg.fit(X, y, sample_weight=sample_weight)
    assert_array_less(1e-3, np.max(np.abs(reg.coef_)))

    if positive:
        # Make sure that the positive constraint changes alpha_max,
        # i.e. test the meaningfulness of the test data.
        not_positive_alpha_max = (
            estimatorCV(alphas=1, cv=2, eps=1, **{**params, "positive": not positive})
            .fit(X, y, sample_weight=sample_weight)
            .alpha_
        )
        assert not np.isclose(alpha_max, not_positive_alpha_max), (
            "Test data cannot distinguish alpha_max between positive=True and False."
        )