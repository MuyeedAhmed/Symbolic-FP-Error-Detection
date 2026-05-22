def test_svr_coef_sign(global_random_seed):
    # Test that SVR(kernel="linear") has coef_ with the right sign.
    # Non-regression test for #2933.
    rng = np.random.RandomState(global_random_seed)
    X = rng.randn(10, 3)
    y = rng.randn(10)

    for svr in [
        svm.SVR(kernel="linear"),
        svm.NuSVR(kernel="linear"),
        svm.LinearSVR(),
    ]:
        svr.fit(X, y)
        assert_array_almost_equal(
            svr.predict(X), np.dot(X, svr.coef_.ravel()) + svr.intercept_
        )