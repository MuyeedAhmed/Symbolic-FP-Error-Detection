def test_iterative_imputer_transform_recovery(rank):
    rng = np.random.RandomState(0)
    n = 70
    d = 70
    A = rng.rand(n, rank)
    B = rng.rand(rank, d)
    X_filled = np.dot(A, B)
    nan_mask = rng.rand(n, d) < 0.5
    X_missing = X_filled.copy()
    X_missing[nan_mask] = np.nan

    # split up data in half
    n = n // 2
    X_train = X_missing[:n]
    X_test_filled = X_filled[n:]
    X_test = X_missing[n:]

    imputer = IterativeImputer(
        max_iter=5, imputation_order="descending", verbose=1, random_state=rng
    ).fit(X_train)
    X_test_est = imputer.transform(X_test)
    assert_allclose(X_test_filled, X_test_est, atol=0.1)