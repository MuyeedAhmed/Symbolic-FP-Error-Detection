def test_make_regression_multitarget(global_random_seed):
    X, y, c = make_regression(
        n_samples=100,
        n_features=10,
        n_informative=3,
        n_targets=3,
        coef=True,
        noise=1.0,
        random_state=global_random_seed,
    )

    assert X.shape == (100, 10), "X shape mismatch"
    assert y.shape == (100, 3), "y shape mismatch"
    assert c.shape == (10, 3), "coef shape mismatch"
    assert_array_equal(sum(c != 0.0), 3, "Unexpected number of informative features")

    # Test that y ~= np.dot(X, c) + bias + N(0, 1.0)
    assert_almost_equal(np.std(y - np.dot(X, c)), 1.0, decimal=1)