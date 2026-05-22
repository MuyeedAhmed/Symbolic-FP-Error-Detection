def test_ridge_positive_error_test(solver):
    """Test input validation for positive argument in Ridge."""
    alpha = 0.1
    X = np.array([[1, 2], [3, 4]])
    coef = np.array([1, -1])
    y = X @ coef

    model = Ridge(alpha=alpha, positive=True, solver=solver, fit_intercept=False)
    with pytest.raises(ValueError, match="does not support positive"):
        model.fit(X, y)

    with pytest.raises(ValueError, match="only 'lbfgs' solver can be used"):
        _, _ = ridge_regression(
            X, y, alpha, positive=True, solver=solver, return_intercept=False
        )