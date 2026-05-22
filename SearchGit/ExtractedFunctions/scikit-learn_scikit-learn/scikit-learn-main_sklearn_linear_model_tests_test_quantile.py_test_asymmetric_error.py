def test_asymmetric_error(quantile):
    """Test quantile regression for asymmetric distributed targets."""
    n_samples = 1000
    rng = np.random.RandomState(42)
    X = np.concatenate(
        (
            np.abs(rng.randn(n_samples)[:, None]),
            -rng.randint(2, size=(n_samples, 1)),
        ),
        axis=1,
    )
    intercept = 1.23
    coef = np.array([0.5, -2])
    #  Take care that X @ coef + intercept > 0
    assert np.min(X @ coef + intercept) > 0
    # For an exponential distribution with rate lambda, e.g. exp(-lambda * x),
    # the quantile at level q is:
    #   quantile(q) = - log(1 - q) / lambda
    #   scale = 1/lambda = -quantile(q) / log(1 - q)
    y = rng.exponential(
        scale=-(X @ coef + intercept) / np.log(1 - quantile), size=n_samples
    )
    model = QuantileRegressor(
        quantile=quantile,
        alpha=0,
    ).fit(X, y)
    # This test can be made to pass with any solver but in the interest
    # of sparing continuous integration resources, the test is performed
    # with the fastest solver only.

    assert model.intercept_ == approx(intercept, rel=0.2)
    assert_allclose(model.coef_, coef, rtol=0.6)
    assert_allclose(np.mean(model.predict(X) > y), quantile, atol=1e-2)

    # Now compare to Nelder-Mead optimization with L1 penalty
    alpha = 0.01
    model.set_params(alpha=alpha).fit(X, y)
    model_coef = np.r_[model.intercept_, model.coef_]

    def func(coef):
        loss = mean_pinball_loss(y, X @ coef[1:] + coef[0], alpha=quantile)
        L1 = np.sum(np.abs(coef[1:]))
        return loss + alpha * L1

    res = minimize(
        fun=func,
        x0=[1, 0, -1],
        method="Nelder-Mead",
        tol=1e-12,
        options={"maxiter": 2000},
    )

    assert func(model_coef) == approx(func(res.x))
    assert_allclose(model.intercept_, res.x[0])
    assert_allclose(model.coef_, res.x[1:])
    assert_allclose(np.mean(model.predict(X) > y), quantile, atol=1e-2)