def test_suffstat_sk_diag():
    # test against 'full' case
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 2, 2

    resp = rng.rand(n_samples, n_components)
    resp = resp / resp.sum(axis=1)[:, np.newaxis]
    X = rng.rand(n_samples, n_features)
    nk = resp.sum(axis=0)
    xk = np.dot(resp.T, X) / nk[:, np.newaxis]
    covars_pred_full = _estimate_gaussian_covariances_full(resp, X, nk, xk, 0)
    covars_pred_diag = _estimate_gaussian_covariances_diag(resp, X, nk, xk, 0)

    ecov = EmpiricalCovariance()
    for cov_full, cov_diag in zip(covars_pred_full, covars_pred_diag):
        ecov.covariance_ = np.diag(np.diag(cov_full))
        cov_diag = np.diag(cov_diag)
        assert_almost_equal(ecov.error_norm(cov_diag, norm="frobenius"), 0)
        assert_almost_equal(ecov.error_norm(cov_diag, norm="spectral"), 0)

    # check the precision computation
    precs_chol_pred = _compute_precision_cholesky(covars_pred_diag, "diag")
    assert_almost_equal(covars_pred_diag, 1.0 / precs_chol_pred**2)