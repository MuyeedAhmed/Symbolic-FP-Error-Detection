def test_implicit_center_matvec(global_random_seed, centered_matrices):
    X_sparse_centered, X_dense_centered = centered_matrices
    rng = np.random.default_rng(global_random_seed)
    y = rng.standard_normal(X_dense_centered.shape[1])
    assert_allclose(X_dense_centered @ y, X_sparse_centered.matvec(y))
    assert_allclose(X_dense_centered @ y, X_sparse_centered @ y)