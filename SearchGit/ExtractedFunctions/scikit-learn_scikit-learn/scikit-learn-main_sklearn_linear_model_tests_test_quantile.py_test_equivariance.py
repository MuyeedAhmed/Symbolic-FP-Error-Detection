def test_equivariance(quantile):
    """Test equivariace of quantile regression.

    See Koenker (2005) Quantile Regression, Chapter 2.2.3.
    """
    rng = np.random.RandomState(42)
    n_samples, n_features = 100, 5
    X, y = make_regression(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_features,
        noise=0,
        random_state=rng,
        shuffle=False,
    )
    # make y asymmetric
    y += rng.exponential(scale=100, size=y.shape)
    params = dict(alpha=0)
    model1 = QuantileRegressor(quantile=quantile, **params).fit(X, y)

    # coef(q; a*y, X) = a * coef(q; y, X)
    a = 2.5
    model2 = QuantileRegressor(quantile=quantile, **params).fit(X, a * y)
    assert model2.intercept_ == approx(a * model1.intercept_, rel=1e-5)
    assert_allclose(model2.coef_, a * model1.coef_, rtol=1e-5)

    # coef(1-q; -a*y, X) = -a * coef(q; y, X)
    model2 = QuantileRegressor(quantile=1 - quantile, **params).fit(X, -a * y)
    assert model2.intercept_ == approx(-a * model1.intercept_, rel=1e-5)
    assert_allclose(model2.coef_, -a * model1.coef_, rtol=1e-5)

    # coef(q; y + X @ g, X) = coef(q; y, X) + g
    g_intercept, g_coef = rng.randn(), rng.randn(n_features)
    model2 = QuantileRegressor(quantile=quantile, **params)
    model2.fit(X, y + X @ g_coef + g_intercept)
    assert model2.intercept_ == approx(model1.intercept_ + g_intercept)
    assert_allclose(model2.coef_, model1.coef_ + g_coef, rtol=1e-6)

    # coef(q; y, X @ A) = A^-1 @ coef(q; y, X)
    A = rng.randn(n_features, n_features)
    model2 = QuantileRegressor(quantile=quantile, **params)
    model2.fit(X @ A, y)
    assert model2.intercept_ == approx(model1.intercept_, rel=1e-5)
    assert_allclose(model2.coef_, np.linalg.solve(A, model1.coef_), rtol=1e-5)