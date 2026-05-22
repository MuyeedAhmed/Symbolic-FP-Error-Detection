def test_nystroem_approximation():
    # some basic tests
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))

    # With n_components = n_samples this is exact
    X_transformed = Nystroem(n_components=X.shape[0]).fit_transform(X)
    K = rbf_kernel(X)
    assert_array_almost_equal(np.dot(X_transformed, X_transformed.T), K)

    trans = Nystroem(n_components=2, random_state=rnd)
    X_transformed = trans.fit(X).transform(X)
    assert X_transformed.shape == (X.shape[0], 2)

    # test callable kernel
    trans = Nystroem(n_components=2, kernel=_linear_kernel, random_state=rnd)
    X_transformed = trans.fit(X).transform(X)
    assert X_transformed.shape == (X.shape[0], 2)

    # test that available kernels fit and transform
    kernels_available = kernel_metrics()
    for kern in kernels_available:
        trans = Nystroem(n_components=2, kernel=kern, random_state=rnd)
        X_transformed = trans.fit(X).transform(X)
        assert X_transformed.shape == (X.shape[0], 2)