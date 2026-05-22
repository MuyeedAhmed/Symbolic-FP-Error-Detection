def test_bayesian_covariance_matrix(n_samples, n_features, global_random_seed):
    """Check the posterior covariance matrix sigma_

    Non-regression test for https://github.com/scikit-learn/scikit-learn/issues/31093
    """
    X, y = datasets.make_regression(
        n_samples, n_features, random_state=global_random_seed
    )
    reg = BayesianRidge(fit_intercept=False).fit(X, y)
    covariance_matrix = np.linalg.inv(
        reg.lambda_ * np.identity(n_features) + reg.alpha_ * np.dot(X.T, X)
    )
    assert_allclose(reg.sigma_, covariance_matrix, rtol=1e-6)