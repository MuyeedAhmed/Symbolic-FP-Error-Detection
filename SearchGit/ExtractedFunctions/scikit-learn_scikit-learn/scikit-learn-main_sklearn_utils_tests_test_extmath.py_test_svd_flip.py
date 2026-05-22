def test_svd_flip():
    # Check that svd_flip works in both situations, and reconstructs input.
    rs = np.random.RandomState(1999)
    n_samples = 20
    n_features = 10
    X = rs.randn(n_samples, n_features)

    # Check matrix reconstruction
    U, S, Vt = linalg.svd(X, full_matrices=False)
    U1, V1 = svd_flip(U, Vt, u_based_decision=False)
    assert_almost_equal(np.dot(U1 * S, V1), X, decimal=6)

    # Check transposed matrix reconstruction
    XT = X.T
    U, S, Vt = linalg.svd(XT, full_matrices=False)
    U2, V2 = svd_flip(U, Vt, u_based_decision=True)
    assert_almost_equal(np.dot(U2 * S, V2), XT, decimal=6)

    # Check that different flip methods are equivalent under reconstruction
    U_flip1, V_flip1 = svd_flip(U, Vt, u_based_decision=True)
    assert_almost_equal(np.dot(U_flip1 * S, V_flip1), XT, decimal=6)
    U_flip2, V_flip2 = svd_flip(U, Vt, u_based_decision=False)
    assert_almost_equal(np.dot(U_flip2 * S, V_flip2), XT, decimal=6)