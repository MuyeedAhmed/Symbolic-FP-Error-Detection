def test_implicit_center_rmatmat(global_random_seed, centered_matrices):
    X_sparse_centered, X_dense_centered = centered_matrices
    rng = np.random.default_rng(global_random_seed)
    Y = rng.standard_normal((X_dense_centered.shape[0], 50))
    assert_allclose(X_dense_centered.T @ Y, X_sparse_centered.rmatmat(Y))
    assert_allclose(X_dense_centered.T @ Y, X_sparse_centered.T @ Y)