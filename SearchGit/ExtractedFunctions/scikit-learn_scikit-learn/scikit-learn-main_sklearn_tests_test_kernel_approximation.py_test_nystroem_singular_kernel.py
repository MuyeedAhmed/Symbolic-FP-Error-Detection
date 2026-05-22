def test_nystroem_singular_kernel():
    # test that nystroem works with singular kernel matrix
    rng = np.random.RandomState(0)
    X = rng.rand(10, 20)
    X = np.vstack([X] * 2)  # duplicate samples

    gamma = 100
    N = Nystroem(gamma=gamma, n_components=X.shape[0]).fit(X)
    X_transformed = N.transform(X)

    K = rbf_kernel(X, gamma=gamma)

    assert_array_almost_equal(K, np.dot(X_transformed, X_transformed.T))
    assert np.all(np.isfinite(Y))