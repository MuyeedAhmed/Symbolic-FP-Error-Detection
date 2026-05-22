def test_fitted_precision_covariance_concistency(covar_type, global_dtype):
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7, dtype=global_dtype)
    n_components = rand_data.n_components

    X = rand_data.X[covar_type]
    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type=covar_type,
        random_state=rng,
        n_init=5,
    )
    gmm.fit(X)
    assert gmm.precisions_.dtype == global_dtype
    assert gmm.covariances_.dtype == global_dtype
    if covar_type == "full":
        for prec, covar in zip(gmm.precisions_, gmm.covariances_):
            assert_array_almost_equal(linalg.inv(prec), covar)
    elif covar_type == "tied":
        assert_array_almost_equal(linalg.inv(gmm.precisions_), gmm.covariances_)
    else:
        assert_array_almost_equal(gmm.precisions_, 1.0 / gmm.covariances_)