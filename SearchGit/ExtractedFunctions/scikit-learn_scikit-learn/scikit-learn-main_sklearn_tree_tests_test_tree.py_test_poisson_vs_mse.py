def test_poisson_vs_mse():
    # For a Poisson distributed target, Poisson loss should give better results
    # than squared error measured in Poisson deviance as metric.
    # We have a similar test, test_poisson(), in
    # sklearn/ensemble/_hist_gradient_boosting/tests/test_gradient_boosting.py
    rng = np.random.RandomState(42)
    n_train, n_test, n_features = 500, 500, 10
    X = datasets.make_low_rank_matrix(
        n_samples=n_train + n_test, n_features=n_features, random_state=rng
    )
    # We create a log-linear Poisson model and downscale coef as it will get
    # exponentiated.
    coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)
    y = rng.poisson(lam=np.exp(X @ coef))
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=n_test, random_state=rng
    )
    # We prevent some overfitting by setting min_samples_split=10.
    tree_poi = DecisionTreeRegressor(
        criterion="poisson", min_samples_split=10, random_state=rng
    )
    tree_mse = DecisionTreeRegressor(
        criterion="squared_error", min_samples_split=10, random_state=rng
    )

    tree_poi.fit(X_train, y_train)
    tree_mse.fit(X_train, y_train)
    dummy = DummyRegressor(strategy="mean").fit(X_train, y_train)

    for X, y, val in [(X_train, y_train, "train"), (X_test, y_test, "test")]:
        metric_poi = mean_poisson_deviance(y, tree_poi.predict(X))
        # squared_error might produce non-positive predictions => clip
        metric_mse = mean_poisson_deviance(y, np.clip(tree_mse.predict(X), 1e-15, None))
        metric_dummy = mean_poisson_deviance(y, dummy.predict(X))
        # As squared_error might correctly predict 0 in train set, its train
        # score can be better than Poisson. This is no longer the case for the
        # test set.
        if val == "test":
            assert metric_poi < 0.5 * metric_mse
        assert metric_poi < 0.75 * metric_dummy