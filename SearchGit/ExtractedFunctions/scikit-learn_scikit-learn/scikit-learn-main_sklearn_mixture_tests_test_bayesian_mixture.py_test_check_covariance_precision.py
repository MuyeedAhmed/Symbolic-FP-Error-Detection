def test_check_covariance_precision():
    # We check that the dot product of the covariance and the precision
    # matrices is identity.
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components, n_features = 2 * rand_data.n_components, 2

    # Computation of the full_covariance
    bgmm = BayesianGaussianMixture(
        n_components=n_components, max_iter=100, random_state=rng, tol=1e-3, reg_covar=0
    )
    for covar_type in COVARIANCE_TYPE:
        bgmm.covariance_type = covar_type
        bgmm.fit(rand_data.X[covar_type])

        if covar_type == "full":
            for covar, precision in zip(bgmm.covariances_, bgmm.precisions_):
                assert_almost_equal(np.dot(covar, precision), np.eye(n_features))
        elif covar_type == "tied":
            assert_almost_equal(
                np.dot(bgmm.covariances_, bgmm.precisions_), np.eye(n_features)
            )

        elif covar_type == "diag":
            assert_almost_equal(
                bgmm.covariances_ * bgmm.precisions_,
                np.ones((n_components, n_features)),
            )

        else:
            assert_almost_equal(
                bgmm.covariances_ * bgmm.precisions_, np.ones(n_components)
            )