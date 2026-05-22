def test_pca_singular_values(svd_solver):
    rng = np.random.RandomState(0)
    n_samples, n_features = 100, 80
    X = rng.randn(n_samples, n_features)

    pca = PCA(n_components=2, svd_solver=svd_solver, random_state=rng)
    X_trans = pca.fit_transform(X)

    # compare to the Frobenius norm
    assert_allclose(
        np.sum(pca.singular_values_**2), np.linalg.norm(X_trans, "fro") ** 2
    )
    # Compare to the 2-norms of the score vectors
    assert_allclose(pca.singular_values_, np.sqrt(np.sum(X_trans**2, axis=0)))

    # set the singular values and see what er get back
    n_samples, n_features = 100, 110
    X = rng.randn(n_samples, n_features)

    pca = PCA(n_components=3, svd_solver=svd_solver, random_state=rng)
    X_trans = pca.fit_transform(X)
    X_trans /= np.sqrt(np.sum(X_trans**2, axis=0))
    X_trans[:, 0] *= 3.142
    X_trans[:, 1] *= 2.718
    X_hat = np.dot(X_trans, pca.components_)
    pca.fit(X_hat)
    assert_allclose(pca.singular_values_, [3.142, 2.718, 1.0])