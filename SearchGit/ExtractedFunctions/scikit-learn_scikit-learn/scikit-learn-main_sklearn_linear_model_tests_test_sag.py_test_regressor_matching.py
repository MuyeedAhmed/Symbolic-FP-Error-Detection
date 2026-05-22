def test_regressor_matching():
    n_samples = 10
    n_features = 5

    rng = np.random.RandomState(10)
    X = rng.normal(size=(n_samples, n_features))
    true_w = rng.normal(size=n_features)
    y = X.dot(true_w)

    alpha = 1.0
    n_iter = 100
    fit_intercept = True

    step_size = get_step_size(X, alpha, fit_intercept, classification=False)
    clf = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00000000001,
        solver="sag",
        alpha=alpha * n_samples,
        max_iter=n_iter,
    )
    clf.fit(X, y)

    weights1, intercept1 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
    )
    weights2, intercept2 = sag(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
    )

    assert_allclose(weights1, clf.coef_)
    assert_allclose(intercept1, clf.intercept_)
    assert_allclose(weights2, clf.coef_)
    assert_allclose(intercept2, clf.intercept_)