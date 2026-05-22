def make_sparse_data(
    sparse_container,
    n_samples=100,
    n_features=100,
    n_informative=10,
    seed=42,
    positive=False,
    n_targets=1,
):
    random_state = np.random.RandomState(seed)

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples

    # generate a ground truth model
    w = random_state.randn(n_features, n_targets)
    w[n_informative:] = 0.0  # only the top features are impacting the model
    if positive:
        w = np.abs(w)

    X = random_state.randn(n_samples, n_features)
    rnd = random_state.uniform(size=(n_samples, n_features))
    X[rnd > 0.5] = 0.0  # 50% of zeros in input signal

    # generate training ground truth labels
    y = np.dot(X, w)
    X = sparse_container(X)
    if n_targets == 1:
        y = np.ravel(y)
    return X, y