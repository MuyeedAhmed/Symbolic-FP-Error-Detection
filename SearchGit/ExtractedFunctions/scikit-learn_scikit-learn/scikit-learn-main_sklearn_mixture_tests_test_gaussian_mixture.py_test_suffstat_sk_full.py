def test_suffstat_sk_full():
    # compare the precision matrix compute from the
    # EmpiricalCovariance.covariance fitted on X*sqrt(resp)
    # with _sufficient_sk_full, n_components=1
    rng = np.random.RandomState(0)
    n_samples, n_features = 500, 2

    # special case 1, assuming data is "centered"
    X = rng.rand(n_samples, n_features)
    resp = rng.rand(n_samples, 1)
    X_resp = np.sqrt(resp) * X
    nk = np.array([n_samples])
    xk = np.zeros((1, n_features))
    covars_pred = _estimate_gaussian_covariances_full(resp, X, nk, xk, 0)
    ecov = EmpiricalCovariance(assume_centered=True)
    ecov.fit(X_resp)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm="frobenius"), 0)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm="spectral"), 0)

    # check the precision computation
    precs_chol_pred = _compute_precision_cholesky(covars_pred, "full")
    precs_pred = np.array([np.dot(prec, prec.T) for prec in precs_chol_pred])
    precs_est = np.array([linalg.inv(cov) for cov in covars_pred])
    assert_array_almost_equal(precs_est, precs_pred)

    # special case 2, assuming resp are all ones
    resp = np.ones((n_samples, 1))
    nk = np.array([n_samples])
    xk = X.mean(axis=0).reshape((1, -1))
    covars_pred = _estimate_gaussian_covariances_full(resp, X, nk, xk, 0)
    ecov = EmpiricalCovariance(assume_centered=False)
    ecov.fit(X)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm="frobenius"), 0)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm="spectral"), 0)

    # check the precision computation
    precs_chol_pred = _compute_precision_cholesky(covars_pred, "full")
    precs_pred = np.array([np.dot(prec, prec.T) for prec in precs_chol_pred])
    precs_est = np.array([linalg.inv(cov) for cov in covars_pred])
    assert_array_almost_equal(precs_est, precs_pred)