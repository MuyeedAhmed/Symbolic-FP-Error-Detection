def test_kernelcenterer_non_linear_kernel():
    """Check kernel centering for non-linear kernel."""
    rng = np.random.RandomState(0)
    X, X_test = rng.randn(100, 50), rng.randn(20, 50)

    def phi(X):
        """Our mapping function phi."""
        return np.vstack(
            [
                np.clip(X, a_min=0, a_max=None),
                -np.clip(X, a_min=None, a_max=0),
            ]
        )

    phi_X = phi(X)
    phi_X_test = phi(X_test)

    # centered the projection
    scaler = StandardScaler(with_std=False)
    phi_X_center = scaler.fit_transform(phi_X)
    phi_X_test_center = scaler.transform(phi_X_test)

    # create the different kernel
    K = phi_X @ phi_X.T
    K_test = phi_X_test @ phi_X.T
    K_center = phi_X_center @ phi_X_center.T
    K_test_center = phi_X_test_center @ phi_X_center.T

    kernel_centerer = KernelCenterer()
    kernel_centerer.fit(K)

    assert_allclose(kernel_centerer.transform(K), K_center)
    assert_allclose(kernel_centerer.transform(K_test), K_test_center)

    # check the results coherence with the method proposed in:
    # B. Schölkopf, A. Smola, and K.R. Müller,
    # "Nonlinear component analysis as a kernel eigenvalue problem"
    # equation (B.3)

    # K_centered = (I - 1_M) K (I - 1_M)
    #            =  K - 1_M K - K 1_M + 1_M K 1_M
    ones_M = np.ones_like(K) / K.shape[0]
    K_centered = K - ones_M @ K - K @ ones_M + ones_M @ K @ ones_M
    assert_allclose(kernel_centerer.transform(K), K_centered)

    # K_test_centered = (K_test - 1'_M K)(I - 1_M)
    #                 = K_test - 1'_M K - K_test 1_M + 1'_M K 1_M
    ones_prime_M = np.ones_like(K_test) / K.shape[0]
    K_test_centered = (
        K_test - ones_prime_M @ K - K_test @ ones_M + ones_prime_M @ K @ ones_M
    )
    assert_allclose(kernel_centerer.transform(K_test), K_test_centered)