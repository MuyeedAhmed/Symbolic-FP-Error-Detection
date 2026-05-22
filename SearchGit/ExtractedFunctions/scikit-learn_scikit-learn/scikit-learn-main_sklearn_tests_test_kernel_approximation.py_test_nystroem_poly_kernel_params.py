def test_nystroem_poly_kernel_params():
    # Non-regression: Nystroem should pass other parameters beside gamma.
    rnd = np.random.RandomState(37)
    X = rnd.uniform(size=(10, 4))

    K = polynomial_kernel(X, degree=3.1, coef0=0.1)
    nystroem = Nystroem(
        kernel="polynomial", n_components=X.shape[0], degree=3.1, coef0=0.1
    )
    X_transformed = nystroem.fit_transform(X)
    assert_array_almost_equal(np.dot(X_transformed, X_transformed.T), K)