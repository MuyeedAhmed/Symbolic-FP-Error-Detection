def test_suffstat_sk_tied():
    # use equation Nk * Sk / N = S_tied
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 2, 2

    resp = rng.rand(n_samples, n_components)
    resp = resp / resp.sum(axis=1)[:, np.newaxis]
    X = rng.rand(n_samples, n_features)
    nk = resp.sum(axis=0)
    xk = np.dot(resp.T, X) / nk[:, np.newaxis]

    covars_pred_full = _estimate_gaussian_covariances_full(resp, X, nk, xk, 0)
    covars_pred_full = (
        np.sum(nk[:, np.newaxis, np.newaxis] * covars_pred_full, 0) / n_samples
    )

    covars_pred_tied = _estimate_gaussian_covariances_tied(resp, X, nk, xk, 0)

    ecov = EmpiricalCovariance()
    ecov.covariance_ = covars_pred_full
    assert_almost_equal(ecov.error_norm(covars_pred_tied, norm="frobenius"), 0)
    assert_almost_equal(ecov.error_norm(covars_pred_tied, norm="spectral"), 0)

    # check the precision computation
    precs_chol_pred = _compute_precision_cholesky(covars_pred_tied, "tied")
    precs_pred = np.dot(precs_chol_pred, precs_chol_pred.T)
    precs_est = linalg.inv(covars_pred_tied)
    assert_array_almost_equal(precs_est, precs_pred)