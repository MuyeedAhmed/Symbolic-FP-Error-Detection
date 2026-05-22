def test_linear_regression_vs_lstsq(dtype):
    """
    Check that LinearRegression is as good as `scipy.linalg.lstsq`.
    Non regression test for issue #33032.
    """
    rng = np.random.RandomState(1137)
    n_samples = 500_000

    x1 = rng.rand(n_samples)
    x2 = 0.3 * x1 + 0.1 * rng.rand(n_samples)
    X = np.column_stack([x1, x2])
    y = X @ [0.5, 2.0] + 0.1 * rng.rand(n_samples)

    X = X.astype(dtype)
    y = y.astype(dtype)

    coef_scipy = linalg.lstsq(X, y)[0]
    coef_sklearn = LinearRegression(fit_intercept=False).fit(X, y).coef_

    rmse_scipy = np.linalg.norm(y - X @ coef_scipy)
    rmse_sklearn = np.linalg.norm(y - X @ coef_sklearn)
    assert rmse_sklearn == pytest.approx(rmse_scipy, rel=1e-6)