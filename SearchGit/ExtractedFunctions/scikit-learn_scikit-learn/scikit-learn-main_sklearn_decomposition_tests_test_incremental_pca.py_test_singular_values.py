def test_singular_values(global_random_seed):
    # Check that the IncrementalPCA output has the correct singular values

    rng = np.random.RandomState(global_random_seed)
    n_samples = 1000
    n_features = 100

    X = datasets.make_low_rank_matrix(
        n_samples, n_features, tail_strength=0.0, effective_rank=10, random_state=rng
    )

    pca = PCA(n_components=10, svd_solver="full", random_state=rng).fit(X)
    ipca = IncrementalPCA(n_components=10, batch_size=150).fit(X)
    assert_array_almost_equal(pca.singular_values_, ipca.singular_values_, 2)

    # Compare to the Frobenius norm
    X_pca = pca.transform(X)
    X_ipca = ipca.transform(X)
    assert_array_almost_equal(
        np.sum(pca.singular_values_**2.0), np.linalg.norm(X_pca, "fro") ** 2.0, 12
    )
    assert_array_almost_equal(
        np.sum(ipca.singular_values_**2.0), np.linalg.norm(X_ipca, "fro") ** 2.0, 2
    )

    # Compare to the 2-norms of the score vectors
    assert_array_almost_equal(
        pca.singular_values_, np.sqrt(np.sum(X_pca**2.0, axis=0)), 12
    )
    assert_array_almost_equal(
        ipca.singular_values_, np.sqrt(np.sum(X_ipca**2.0, axis=0)), 2
    )

    # Set the singular values and see what we get back
    rng = np.random.RandomState(global_random_seed)
    n_samples = 100
    n_features = 110

    X = datasets.make_low_rank_matrix(
        n_samples, n_features, tail_strength=0.0, effective_rank=3, random_state=rng
    )

    pca = PCA(n_components=3, svd_solver="full", random_state=rng)
    ipca = IncrementalPCA(n_components=3, batch_size=100)

    X_pca = pca.fit_transform(X)
    X_pca /= np.sqrt(np.sum(X_pca**2.0, axis=0))
    X_pca[:, 0] *= 3.142
    X_pca[:, 1] *= 2.718

    X_hat = np.dot(X_pca, pca.components_)
    pca.fit(X_hat)
    ipca.fit(X_hat)
    assert_array_almost_equal(pca.singular_values_, [3.142, 2.718, 1.0], 14)
    assert_array_almost_equal(ipca.singular_values_, [3.142, 2.718, 1.0], 14)