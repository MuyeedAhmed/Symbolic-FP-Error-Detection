def test_whitening(solver, copy):
    # Check that PCA output has unit-variance
    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 80
    n_components = 30
    rank = 50

    # some low rank data with correlated features
    X = np.dot(
        rng.randn(n_samples, rank),
        np.dot(np.diag(np.linspace(10.0, 1.0, rank)), rng.randn(rank, n_features)),
    )
    # the component-wise variance of the first 50 features is 3 times the
    # mean component-wise variance of the remaining 30 features
    X[:, :50] *= 3

    assert X.shape == (n_samples, n_features)

    # the component-wise variance is thus highly varying:
    assert X.std(axis=0).std() > 43.8

    # whiten the data while projecting to the lower dim subspace
    X_ = X.copy()  # make sure we keep an original across iterations.
    pca = PCA(
        n_components=n_components,
        whiten=True,
        copy=copy,
        svd_solver=solver,
        random_state=0,
        iterated_power=7,
    )
    # test fit_transform
    X_whitened = pca.fit_transform(X_.copy())
    assert X_whitened.shape == (n_samples, n_components)
    X_whitened2 = pca.transform(X_)
    assert_allclose(X_whitened, X_whitened2, rtol=5e-4)

    assert_allclose(X_whitened.std(ddof=1, axis=0), np.ones(n_components))
    assert_allclose(X_whitened.mean(axis=0), np.zeros(n_components), atol=1e-12)

    X_ = X.copy()
    pca = PCA(
        n_components=n_components, whiten=False, copy=copy, svd_solver=solver
    ).fit(X_.copy())
    X_unwhitened = pca.transform(X_)
    assert X_unwhitened.shape == (n_samples, n_components)

    # in that case the output components still have varying variances
    assert X_unwhitened.std(axis=0).std() == pytest.approx(74.1, rel=1e-1)