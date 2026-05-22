def test_enet_float_precision():
    # Generate dataset
    X, y, X_test, y_test = build_dataset(n_samples=20, n_features=10)
    # Here we have a small number of iterations, and thus the
    # ElasticNet might not converge. This is to speed up tests

    for fit_intercept in [True, False]:
        coef = {}
        intercept = {}
        for dtype in [np.float64, np.float32]:
            clf = ElasticNet(
                alpha=0.5,
                max_iter=100,
                precompute=False,
                fit_intercept=fit_intercept,
            )

            X = dtype(X)
            y = dtype(y)
            ignore_warnings(clf.fit)(X, y)

            coef[("simple", dtype)] = clf.coef_
            intercept[("simple", dtype)] = clf.intercept_

            assert clf.coef_.dtype == dtype

            # test precompute Gram array
            Gram = X.T.dot(X)
            clf_precompute = ElasticNet(
                alpha=0.5,
                max_iter=100,
                precompute=Gram,
                fit_intercept=fit_intercept,
            )
            ignore_warnings(clf_precompute.fit)(X, y)
            assert_array_almost_equal(clf.coef_, clf_precompute.coef_)
            assert_array_almost_equal(clf.intercept_, clf_precompute.intercept_)

            # test multi task enet
            multi_y = np.hstack((y[:, np.newaxis], y[:, np.newaxis]))
            clf_multioutput = MultiTaskElasticNet(
                alpha=0.5,
                max_iter=100,
                fit_intercept=fit_intercept,
            )
            clf_multioutput.fit(X, multi_y)
            coef[("multi", dtype)] = clf_multioutput.coef_
            intercept[("multi", dtype)] = clf_multioutput.intercept_
            assert clf.coef_.dtype == dtype

        for v in ["simple", "multi"]:
            assert_array_almost_equal(
                coef[(v, np.float32)], coef[(v, np.float64)], decimal=4
            )
            assert_array_almost_equal(
                intercept[(v, np.float32)], intercept[(v, np.float64)], decimal=4
            )