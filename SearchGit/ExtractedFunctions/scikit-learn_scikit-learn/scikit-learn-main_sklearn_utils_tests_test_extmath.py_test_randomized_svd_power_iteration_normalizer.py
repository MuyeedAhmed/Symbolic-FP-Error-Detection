def test_randomized_svd_power_iteration_normalizer():
    # randomized_svd with power_iteration_normalized='none' diverges for
    # large number of power iterations on this dataset
    rng = np.random.RandomState(42)
    X = make_low_rank_matrix(100, 500, effective_rank=50, random_state=rng)
    X += 3 * rng.randint(0, 2, size=X.shape)
    n_components = 50

    # Check that it diverges with many (non-normalized) power iterations
    U, s, Vt = randomized_svd(
        X, n_components, n_iter=2, power_iteration_normalizer="none", random_state=0
    )
    A = X - U.dot(np.diag(s).dot(Vt))
    error_2 = linalg.norm(A, ord="fro")
    U, s, Vt = randomized_svd(
        X, n_components, n_iter=20, power_iteration_normalizer="none", random_state=0
    )
    A = X - U.dot(np.diag(s).dot(Vt))
    error_20 = linalg.norm(A, ord="fro")
    assert np.abs(error_2 - error_20) > 100

    for normalizer in ["LU", "QR", "auto"]:
        U, s, Vt = randomized_svd(
            X,
            n_components,
            n_iter=2,
            power_iteration_normalizer=normalizer,
            random_state=0,
        )
        A = X - U.dot(np.diag(s).dot(Vt))
        error_2 = linalg.norm(A, ord="fro")

        for i in [5, 10, 50]:
            U, s, Vt = randomized_svd(
                X,
                n_components,
                n_iter=i,
                power_iteration_normalizer=normalizer,
                random_state=0,
            )
            A = X - U.dot(np.diag(s).dot(Vt))
            error = linalg.norm(A, ord="fro")
            assert 15 > np.abs(error_2 - error)