def test_poisson_vs_mse():
    """Test that random forest with poisson criterion performs better than
    mse for a poisson target.

    There is a similar test for DecisionTreeRegressor.
    """
    rng = np.random.RandomState(42)
    n_train, n_test, n_features = 500, 500, 10
    X = datasets.make_low_rank_matrix(
        n_samples=n_train + n_test, n_features=n_features, random_state=rng
    )
    #  We create a log-linear Poisson model and downscale coef as it will get
    # exponentiated.
    coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)
    y = rng.poisson(lam=np.exp(X @ coef))
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=n_test, random_state=rng
    )
    # We prevent some overfitting by setting min_samples_split=10.
    forest_poi = RandomForestRegressor(
        criterion="poisson", min_samples_leaf=10, max_features="sqrt", random_state=rng
    )
    forest_mse = RandomForestRegressor(
        criterion="squared_error",
        min_samples_leaf=10,
        max_features="sqrt",
        random_state=rng,
    )

    forest_poi.fit(X_train, y_train)
    forest_mse.fit(X_train, y_train)
    dummy = DummyRegressor(strategy="mean").fit(X_train, y_train)

    for X, y, data_name in [(X_train, y_train, "train"), (X_test, y_test, "test")]:
        metric_poi = mean_poisson_deviance(y, forest_poi.predict(X))
        # squared_error forest might produce non-positive predictions => clip
        # If y = 0 for those, the poisson deviance gets too good.
        # If we drew more samples, we would eventually get y > 0 and the
        # poisson deviance would explode, i.e. be undefined. Therefore, we do
        # not clip to a tiny value like 1e-15, but to 1e-6. This acts like a
        # small penalty to the non-positive predictions.
        metric_mse = mean_poisson_deviance(
            y, np.clip(forest_mse.predict(X), 1e-6, None)
        )
        metric_dummy = mean_poisson_deviance(y, dummy.predict(X))
        # As squared_error might correctly predict 0 in train set, its train
        # score can be better than Poisson. This is no longer the case for the
        # test set. But keep the above comment for clipping in mind.
        if data_name == "test":
            assert metric_poi < metric_mse
        assert metric_poi < 0.8 * metric_dummy