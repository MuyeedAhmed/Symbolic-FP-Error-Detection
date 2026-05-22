def test_error():
    # Test that it gives proper exception on deficient input.
    for name, TreeEstimator in CLF_TREES.items():
        # predict before fit
        est = TreeEstimator()
        with pytest.raises(NotFittedError):
            est.predict_proba(X)

        est.fit(X, y)
        X2 = [[-2, -1, 1]]  # wrong feature shape for sample
        with pytest.raises(ValueError):
            est.predict_proba(X2)

        # Wrong dimensions
        est = TreeEstimator()
        y2 = y[:-1]
        with pytest.raises(ValueError):
            est.fit(X, y2)

        # Test with arrays that are non-contiguous.
        Xf = np.asfortranarray(X)
        est = TreeEstimator()
        est.fit(Xf, y)
        assert_almost_equal(est.predict(T), true_result)

        # predict before fitting
        est = TreeEstimator()
        with pytest.raises(NotFittedError):
            est.predict(T)

        # predict on vector with different dims
        est.fit(X, y)
        t = np.asarray(T)
        with pytest.raises(ValueError):
            est.predict(t[:, 1:])

        # wrong sample shape
        Xt = np.array(X).T

        est = TreeEstimator()
        est.fit(np.dot(X, Xt), y)
        with pytest.raises(ValueError):
            est.predict(X)
        with pytest.raises(ValueError):
            est.apply(X)

        clf = TreeEstimator()
        clf.fit(X, y)
        with pytest.raises(ValueError):
            clf.predict(Xt)
        with pytest.raises(ValueError):
            clf.apply(Xt)

        # apply before fitting
        est = TreeEstimator()
        with pytest.raises(NotFittedError):
            est.apply(T)

    # non positive target for Poisson splitting Criterion
    est = DecisionTreeRegressor(criterion="poisson")
    with pytest.raises(ValueError, match="y is not positive.*Poisson"):
        est.fit([[0, 1, 2]], [0, 0, 0])
    with pytest.raises(ValueError, match="Some.*y are negative.*Poisson"):
        est.fit([[0, 1, 2]], [5, -0.1, 2])