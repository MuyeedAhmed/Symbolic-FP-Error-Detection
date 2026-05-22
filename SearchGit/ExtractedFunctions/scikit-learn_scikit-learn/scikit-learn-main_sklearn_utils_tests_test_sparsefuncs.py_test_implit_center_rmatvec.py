def test_implit_center_rmatvec(global_random_seed, centered_matrices):
    X_sparse_centered, X_dense_centered = centered_matrices
    rng = np.random.default_rng(global_random_seed)
    y = rng.standard_normal(X_dense_centered.shape[0])
    assert_allclose(X_dense_centered.T @ y, X_sparse_centered.rmatvec(y))
    assert_allclose(X_dense_centered.T @ y, X_sparse_centered.T @ y)