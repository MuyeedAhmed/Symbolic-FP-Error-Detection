def test_gaussian_suffstat_sk_spherical(global_dtype):
    # computing spherical covariance equals to the variance of one-dimension
    # data after flattening, n_components=1
    rng = np.random.RandomState(0)
    n_samples, n_features = 500, 2

    X = rng.rand(n_samples, n_features).astype(global_dtype)
    X = X - X.mean()
    resp = np.ones((n_samples, 1), dtype=global_dtype)
    nk = np.array([n_samples], dtype=global_dtype)
    xk = X.mean()
    covars_pred_spherical = _estimate_gaussian_covariances_spherical(resp, X, nk, xk, 0)
    covars_pred_spherical2 = np.dot(X.flatten().T, X.flatten()) / (
        n_features * n_samples
    )
    assert_almost_equal(covars_pred_spherical, covars_pred_spherical2)
    assert covars_pred_spherical.dtype == global_dtype

    # check the precision computation
    precs_chol_pred = _compute_precision_cholesky(covars_pred_spherical, "spherical")
    assert_almost_equal(covars_pred_spherical, 1.0 / precs_chol_pred**2)
    assert precs_chol_pred.dtype == global_dtype