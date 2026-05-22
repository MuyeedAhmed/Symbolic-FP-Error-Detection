def test_ridge_regression_dtype_stability(solver, seed):
    random_state = np.random.RandomState(seed)
    n_samples, n_features = 6, 5
    X = random_state.randn(n_samples, n_features)
    coef = random_state.randn(n_features)
    y = np.dot(X, coef) + 0.01 * random_state.randn(n_samples)
    alpha = 1.0
    positive = solver == "lbfgs"
    results = dict()
    # XXX: Sparse CG seems to be far less numerically stable than the
    # others, maybe we should not enable float32 for this one.
    atol = 1e-3 if solver == "sparse_cg" else 1e-5
    for current_dtype in (np.float32, np.float64):
        results[current_dtype] = ridge_regression(
            X.astype(current_dtype),
            y.astype(current_dtype),
            alpha=alpha,
            solver=solver,
            random_state=random_state,
            sample_weight=None,
            positive=positive,
            max_iter=500,
            tol=1e-10,
            return_n_iter=False,
            return_intercept=False,
        )

    assert results[np.float32].dtype == np.float32
    assert results[np.float64].dtype == np.float64
    assert_allclose(results[np.float32], results[np.float64], atol=atol)