def test_poisson():
    # For Poisson distributed target, Poisson loss should give better results
    # than least squares measured in Poisson deviance as metric.
    rng = np.random.RandomState(42)
    n_train, n_test, n_features = 500, 100, 100
    X = make_low_rank_matrix(
        n_samples=n_train + n_test, n_features=n_features, random_state=rng
    )
    # We create a log-linear Poisson model and downscale coef as it will get
    # exponentiated.
    coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)
    y = rng.poisson(lam=np.exp(X @ coef))
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=n_test, random_state=rng
    )
    gbdt_pois = HistGradientBoostingRegressor(loss="poisson", random_state=rng)
    gbdt_ls = HistGradientBoostingRegressor(loss="squared_error", random_state=rng)
    gbdt_pois.fit(X_train, y_train)
    gbdt_ls.fit(X_train, y_train)
    dummy = DummyRegressor(strategy="mean").fit(X_train, y_train)

    for X, y in [(X_train, y_train), (X_test, y_test)]:
        metric_pois = mean_poisson_deviance(y, gbdt_pois.predict(X))
        # squared_error might produce non-positive predictions => clip
        metric_ls = mean_poisson_deviance(y, np.clip(gbdt_ls.predict(X), 1e-15, None))
        metric_dummy = mean_poisson_deviance(y, dummy.predict(X))
        assert metric_pois < metric_ls
        assert metric_pois < metric_dummy