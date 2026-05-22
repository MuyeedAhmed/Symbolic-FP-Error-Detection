def test_make_regression(global_random_seed):
    X, y, c = make_regression(
        n_samples=200,
        n_features=10,
        n_informative=3,
        effective_rank=5,
        coef=True,
        bias=0.0,
        noise=1.0,
        random_state=global_random_seed,
    )

    assert X.shape == (200, 10), "X shape mismatch"
    assert y.shape == (200,), "y shape mismatch"
    assert c.shape == (10,), "coef shape mismatch"
    assert sum(c != 0.0) == 3, "Unexpected number of informative features"

    # Test that y ~= np.dot(X, c) + bias + N(0, 1.0).
    assert_almost_equal(np.std(y - np.dot(X, c)), 1.0, decimal=1)

    # Test with small number of features.
    X, y = make_regression(n_samples=100, n_features=1)  # n_informative=3
    assert X.shape == (100, 1)