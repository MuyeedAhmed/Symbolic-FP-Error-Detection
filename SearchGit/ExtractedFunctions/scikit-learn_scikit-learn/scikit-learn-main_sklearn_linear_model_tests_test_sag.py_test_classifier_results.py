def test_classifier_results(csr_container):
    """tests if classifier results match target"""
    alpha = 0.1
    n_features = 20
    n_samples = 10
    tol = 0.01
    max_iter = 200
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)
    y = np.sign(y)
    clf1 = LogisticRegression(
        solver="sag",
        C=1.0 / alpha / n_samples,
        max_iter=max_iter,
        tol=tol,
        random_state=77,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    pred1 = clf1.predict(X)
    pred2 = clf2.predict(X)
    assert_almost_equal(pred1, y, decimal=12)
    assert_almost_equal(pred2, y, decimal=12)