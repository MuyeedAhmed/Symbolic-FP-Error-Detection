def test_graphical_lasso_when_alpha_equals_0(global_random_seed):
    """Test graphical_lasso's early return condition when alpha=0."""
    X = np.random.RandomState(global_random_seed).randn(100, 10)
    emp_cov = empirical_covariance(X, assume_centered=True)

    model = GraphicalLasso(alpha=0, covariance="precomputed").fit(emp_cov)
    assert_allclose(model.precision_, np.linalg.inv(emp_cov))

    _, precision = graphical_lasso(emp_cov, alpha=0)
    assert_allclose(precision, np.linalg.inv(emp_cov))