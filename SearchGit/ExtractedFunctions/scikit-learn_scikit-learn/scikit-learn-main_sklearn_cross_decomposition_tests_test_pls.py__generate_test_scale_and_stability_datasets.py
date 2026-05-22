def _generate_test_scale_and_stability_datasets():
    """Generate dataset for test_scale_and_stability"""
    # dataset for non-regression 7818
    rng = np.random.RandomState(0)
    n_samples = 1000
    n_targets = 5
    n_features = 10
    Q = rng.randn(n_targets, n_features)
    y = rng.randn(n_samples, n_targets)
    X = np.dot(y, Q) + 2 * rng.randn(n_samples, n_features) + 1
    X *= 1000
    yield X, y

    # Data set where one of the features is constraint
    X, y = load_linnerud(return_X_y=True)
    # causes X[:, -1].std() to be zero
    X[:, -1] = 1.0
    yield X, y

    X = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [2.0, 2.0, 2.0], [3.0, 5.0, 4.0]])
    y = np.array([[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]])
    yield X, y

    # Seeds that provide a non-regression test for #18746, where CCA fails
    seeds = [530, 741]
    for seed in seeds:
        rng = np.random.RandomState(seed)
        X = rng.randn(4, 3)
        y = rng.randn(4, 2)
        yield X, y