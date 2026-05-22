def test_lasso_dual_gap():
    """
    Check that Lasso.dual_gap_ matches its objective formulation, with the
    datafit normalized by n_samples
    """
    X, y, _, _ = build_dataset(n_samples=10, n_features=30)
    n_samples = len(y)
    alpha = 0.01 * np.max(np.abs(X.T @ y)) / n_samples
    clf = Lasso(alpha=alpha, fit_intercept=False).fit(X, y)
    w = clf.coef_
    R = y - X @ w
    primal = 0.5 * np.mean(R**2) + clf.alpha * np.sum(np.abs(w))
    # dual pt: R / n_samples, dual constraint: norm(X.T @ theta, inf) <= alpha
    R /= np.max(np.abs(X.T @ R) / (n_samples * alpha))
    dual = 0.5 * (np.mean(y**2) - np.mean((y - R) ** 2))
    assert_allclose(clf.dual_gap_, primal - dual)