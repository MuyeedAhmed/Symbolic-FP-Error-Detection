def test_ridge_positive_regression_test(solver, fit_intercept, alpha):
    """Test that positive Ridge finds true positive coefficients."""
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    coef = np.array([1, -10])
    if fit_intercept:
        intercept = 20
        y = X.dot(coef) + intercept
    else:
        y = X.dot(coef)

    model = Ridge(
        alpha=alpha, positive=True, solver=solver, fit_intercept=fit_intercept
    )
    model.fit(X, y)
    assert np.all(model.coef_ >= 0)