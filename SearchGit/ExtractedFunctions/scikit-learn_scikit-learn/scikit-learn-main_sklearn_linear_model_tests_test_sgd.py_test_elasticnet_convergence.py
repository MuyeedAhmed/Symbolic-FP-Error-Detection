def test_elasticnet_convergence(klass):
    # Check that the SGD output is consistent with coordinate descent

    n_samples, n_features = 1000, 5
    rng = np.random.RandomState(0)
    X = rng.randn(n_samples, n_features)
    # ground_truth linear model that generate y from X and to which the
    # models should converge if the regularizer would be set to 0.0
    ground_truth_coef = rng.randn(n_features)
    y = np.dot(X, ground_truth_coef)

    # XXX: alpha = 0.1 seems to cause convergence problems
    for alpha in [0.01, 0.001]:
        for l1_ratio in [0.5, 0.8, 1.0]:
            cd = linear_model.ElasticNet(
                alpha=alpha, l1_ratio=l1_ratio, fit_intercept=False
            )
            cd.fit(X, y)
            sgd = klass(
                penalty="elasticnet",
                max_iter=50,
                alpha=alpha,
                l1_ratio=l1_ratio,
                fit_intercept=False,
            )
            sgd.fit(X, y)
            err_msg = (
                "cd and sgd did not converge to comparable "
                "results for alpha=%f and l1_ratio=%f" % (alpha, l1_ratio)
            )
            assert_almost_equal(cd.coef_, sgd.coef_, decimal=2, err_msg=err_msg)