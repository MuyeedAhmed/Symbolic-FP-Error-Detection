def test_nystroem_default_parameters():
    rnd = np.random.RandomState(42)
    X = rnd.uniform(size=(10, 4))

    # rbf kernel should behave as gamma=None by default
    # aka gamma = 1 / n_features
    nystroem = Nystroem(n_components=10)
    X_transformed = nystroem.fit_transform(X)
    K = rbf_kernel(X, gamma=None)
    K2 = np.dot(X_transformed, X_transformed.T)
    assert_array_almost_equal(K, K2)

    # chi2 kernel should behave as gamma=1 by default
    nystroem = Nystroem(kernel="chi2", n_components=10)
    X_transformed = nystroem.fit_transform(X)
    K = chi2_kernel(X, gamma=1)
    K2 = np.dot(X_transformed, X_transformed.T)
    assert_array_almost_equal(K, K2)