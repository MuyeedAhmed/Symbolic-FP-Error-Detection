def test_enet_toy():
    # Test ElasticNet for various parameters of alpha and l1_ratio.
    # Actually, the parameters alpha = 0 should not be allowed. However,
    # we test it as a border case.
    # ElasticNet is tested with and without precomputed Gram matrix

    X = np.array([[-1.0], [0.0], [1.0]])
    Y = [-1, 0, 1]  # just a straight line
    T = [[2.0], [3.0], [4.0]]  # test sample

    # this should be the same as lasso
    clf = ElasticNet(alpha=1e-8, l1_ratio=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, l1_ratio=0.3, max_iter=100, precompute=False)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=True)
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=np.dot(X.T, X))
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, l1_ratio=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090, 1.3636, 1.8181], 3)
    assert_almost_equal(clf.dual_gap_, 0)