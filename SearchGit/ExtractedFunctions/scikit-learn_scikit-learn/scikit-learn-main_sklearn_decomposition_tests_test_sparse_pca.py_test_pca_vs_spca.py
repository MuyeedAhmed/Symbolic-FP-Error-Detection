def test_pca_vs_spca(global_random_seed):
    rng = np.random.RandomState(global_random_seed)
    Y, _, _ = generate_toy_data(3, 1000, (8, 8), random_state=rng)
    Z, _, _ = generate_toy_data(3, 10, (8, 8), random_state=rng)
    spca = SparsePCA(alpha=0, ridge_alpha=0, n_components=2, random_state=rng)
    pca = PCA(n_components=2, random_state=rng)
    pca.fit(Y)
    spca.fit(Y)
    results_test_pca = pca.transform(Z)
    results_test_spca = spca.transform(Z)
    assert_allclose(
        np.abs(spca.components_.dot(pca.components_.T)), np.eye(2), atol=1e-4
    )
    results_test_pca *= np.sign(results_test_pca[0, :])
    results_test_spca *= np.sign(results_test_spca[0, :])
    assert_allclose(results_test_pca, results_test_spca, atol=1e-4)