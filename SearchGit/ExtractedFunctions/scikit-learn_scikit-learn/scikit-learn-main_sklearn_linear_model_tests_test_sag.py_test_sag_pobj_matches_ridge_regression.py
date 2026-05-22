def test_sag_pobj_matches_ridge_regression(csr_container):
    """tests if the sag pobj matches ridge reg"""
    n_samples = 100
    n_features = 10
    alpha = 1.0
    n_iter = 100
    fit_intercept = False
    rng = np.random.RandomState(10)
    X = rng.normal(size=(n_samples, n_features))
    true_w = rng.normal(size=n_features)
    y = X.dot(true_w)

    clf1 = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00000000001,
        solver="sag",
        alpha=alpha,
        max_iter=n_iter,
        random_state=42,
    )
    clf2 = clone(clf1)
    clf3 = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00001,
        solver="lsqr",
        alpha=alpha,
        max_iter=n_iter,
        random_state=42,
    )

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, squared_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, squared_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, squared_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj1, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj2, decimal=4)