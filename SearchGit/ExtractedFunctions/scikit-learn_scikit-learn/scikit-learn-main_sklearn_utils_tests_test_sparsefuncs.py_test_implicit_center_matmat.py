def test_implicit_center_matmat(global_random_seed, centered_matrices):
    X_sparse_centered, X_dense_centered = centered_matrices
    rng = np.random.default_rng(global_random_seed)
    Y = rng.standard_normal((X_dense_centered.shape[1], 50))
    assert_allclose(X_dense_centered @ Y, X_sparse_centered.matmat(Y))
    assert_allclose(X_dense_centered @ Y, X_sparse_centered @ Y)