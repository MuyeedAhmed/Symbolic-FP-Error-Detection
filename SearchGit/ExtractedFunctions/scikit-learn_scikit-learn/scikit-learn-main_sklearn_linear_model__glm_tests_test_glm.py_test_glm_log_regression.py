def test_glm_log_regression(solver, fit_intercept, estimator):
    """Test GLM regression with log link on a simple dataset."""
    coef = [0.2, -0.1]
    X = np.array([[0, 1, 2, 3, 4], [1, 1, 1, 1, 1]]).T
    y = np.exp(np.dot(X, coef))
    glm = clone(estimator).set_params(
        alpha=0,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-8,
    )
    if fit_intercept:
        res = glm.fit(X[:, :-1], y)
        assert_allclose(res.coef_, coef[:-1], rtol=1e-6)
        assert_allclose(res.intercept_, coef[-1], rtol=1e-6)
    else:
        res = glm.fit(X, y)
        assert_allclose(res.coef_, coef, rtol=2e-6)