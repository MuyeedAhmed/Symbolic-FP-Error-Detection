def test_lle_simple_grid(global_dtype):
    # note: ARPACK is numerically unstable, so this test will fail for
    #       some random seeds.  We choose 42 because the tests pass.
    #       for arm64 platforms 2 makes the test fail.
    # TODO: rewrite this test to make less sensitive to the random seed,
    # irrespective of the platform.
    rng = np.random.RandomState(42)

    # grid of equidistant points in 2D, n_components = n_dim
    X = np.array(list(product(range(5), repeat=2)))
    X = X + 1e-10 * rng.uniform(size=X.shape)
    X = X.astype(global_dtype, copy=False)

    n_components = 2
    clf = manifold.LocallyLinearEmbedding(
        n_neighbors=5, n_components=n_components, random_state=rng
    )
    tol = 0.1

    N = barycenter_kneighbors_graph(X, clf.n_neighbors).toarray()
    reconstruction_error = linalg.norm(np.dot(N, X) - X, "fro")
    assert reconstruction_error < tol

    for solver in eigen_solvers:
        clf.set_params(eigen_solver=solver)
        clf.fit(X)
        assert clf.embedding_.shape[1] == n_components
        reconstruction_error = (
            linalg.norm(np.dot(N, clf.embedding_) - clf.embedding_, "fro") ** 2
        )

        assert reconstruction_error < tol
        assert_allclose(clf.reconstruction_error_, reconstruction_error, atol=1e-1)

    # re-embed a noisy version of X using the transform method
    noise = rng.randn(*X.shape).astype(global_dtype, copy=False) / 100
    X_reembedded = clf.transform(X + noise)
    assert linalg.norm(X_reembedded - clf.embedding_) < tol