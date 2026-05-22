def test_iterative_imputer_early_stopping():
    rng = np.random.RandomState(0)
    n = 50
    d = 5
    A = rng.rand(n, 1)
    B = rng.rand(1, d)
    X = np.dot(A, B)
    nan_mask = rng.rand(n, d) < 0.5
    X_missing = X.copy()
    X_missing[nan_mask] = np.nan

    imputer = IterativeImputer(
        max_iter=100, tol=1e-2, sample_posterior=False, verbose=1, random_state=rng
    )
    X_filled_100 = imputer.fit_transform(X_missing)
    assert len(imputer.imputation_sequence_) == d * imputer.n_iter_

    imputer = IterativeImputer(
        max_iter=imputer.n_iter_, sample_posterior=False, verbose=1, random_state=rng
    )
    X_filled_early = imputer.fit_transform(X_missing)
    assert_allclose(X_filled_100, X_filled_early, atol=1e-7)

    imputer = IterativeImputer(
        max_iter=100, tol=0, sample_posterior=False, verbose=1, random_state=rng
    )
    imputer.fit(X_missing)
    assert imputer.n_iter_ == imputer.max_iter