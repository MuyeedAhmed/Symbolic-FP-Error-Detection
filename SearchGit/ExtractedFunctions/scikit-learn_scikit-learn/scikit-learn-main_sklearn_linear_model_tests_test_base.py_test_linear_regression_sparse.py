def test_linear_regression_sparse(global_random_seed):
    # Test that linear regression also works with sparse data
    rng = np.random.RandomState(global_random_seed)
    n = 100
    X = _sparse_eye_array(n, n)
    beta = rng.rand(n)
    y = X @ beta

    ols = LinearRegression()
    ols.fit(X, y.ravel())
    assert_array_almost_equal(beta, ols.coef_ + ols.intercept_)

    assert_array_almost_equal(ols.predict(X) - y.ravel(), 0)