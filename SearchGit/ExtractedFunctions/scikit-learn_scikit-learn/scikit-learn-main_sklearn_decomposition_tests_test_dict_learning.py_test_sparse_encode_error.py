def test_sparse_encode_error():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V**2, axis=1)[:, np.newaxis]
    code = sparse_encode(X, V, alpha=0.001)
    assert not np.all(code == 0)
    assert np.sqrt(np.sum((np.dot(code, V) - X) ** 2)) < 0.1