def test_sparse_coder_estimator():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V**2, axis=1)[:, np.newaxis]
    coder = SparseCoder(
        dictionary=V, transform_algorithm="lasso_lars", transform_alpha=0.001
    )
    code = coder.fit_transform(X)
    Xr = coder.inverse_transform(code)
    assert not np.all(code == 0)
    assert np.sqrt(np.sum((np.dot(code, V) - X) ** 2)) < 0.1
    np.testing.assert_allclose(Xr, np.dot(code, V))