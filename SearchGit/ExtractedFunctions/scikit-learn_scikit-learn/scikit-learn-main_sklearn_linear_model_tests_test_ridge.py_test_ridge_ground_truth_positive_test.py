def test_ridge_ground_truth_positive_test(fit_intercept, alpha):
    """Test that Ridge w/wo positive converges to the same solution.

    Ridge with positive=True and positive=False must give the same
    when the ground truth coefs are all positive.
    """
    rng = np.random.RandomState(42)
    X = rng.randn(300, 100)
    coef = rng.uniform(0.1, 1.0, size=X.shape[1])
    if fit_intercept:
        intercept = 1
        y = X @ coef + intercept
    else:
        y = X @ coef
    y += rng.normal(size=X.shape[0]) * 0.01

    results = []
    for positive in [True, False]:
        model = Ridge(
            alpha=alpha, positive=positive, fit_intercept=fit_intercept, tol=1e-10
        )
        results.append(model.fit(X, y).coef_)
    assert_allclose(*results, atol=1e-6, rtol=0)