def test_center_kernel():
    # Test that KernelCenterer is equivalent to StandardScaler
    # in feature space
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    scaler = StandardScaler(with_std=False)
    scaler.fit(X_fit)
    X_fit_centered = scaler.transform(X_fit)
    K_fit = np.dot(X_fit, X_fit.T)

    # center fit time matrix
    centerer = KernelCenterer()
    K_fit_centered = np.dot(X_fit_centered, X_fit_centered.T)
    K_fit_centered2 = centerer.fit_transform(K_fit)
    assert_array_almost_equal(K_fit_centered, K_fit_centered2)

    # center predict time matrix
    X_pred = rng.random_sample((2, 4))
    K_pred = np.dot(X_pred, X_fit.T)
    X_pred_centered = scaler.transform(X_pred)
    K_pred_centered = np.dot(X_pred_centered, X_fit_centered.T)
    K_pred_centered2 = centerer.transform(K_pred)
    assert_array_almost_equal(K_pred_centered, K_pred_centered2)

    # check the results coherence with the method proposed in:
    # B. Schölkopf, A. Smola, and K.R. Müller,
    # "Nonlinear component analysis as a kernel eigenvalue problem"
    # equation (B.3)

    # K_centered3 = (I - 1_M) K (I - 1_M)
    #             =  K - 1_M K - K 1_M + 1_M K 1_M
    ones_M = np.ones_like(K_fit) / K_fit.shape[0]
    K_fit_centered3 = K_fit - ones_M @ K_fit - K_fit @ ones_M + ones_M @ K_fit @ ones_M
    assert_allclose(K_fit_centered, K_fit_centered3)

    # K_test_centered3 = (K_test - 1'_M K)(I - 1_M)
    #                  = K_test - 1'_M K - K_test 1_M + 1'_M K 1_M
    ones_prime_M = np.ones_like(K_pred) / K_fit.shape[0]
    K_pred_centered3 = (
        K_pred - ones_prime_M @ K_fit - K_pred @ ones_M + ones_prime_M @ K_fit @ ones_M
    )
    assert_allclose(K_pred_centered, K_pred_centered3)