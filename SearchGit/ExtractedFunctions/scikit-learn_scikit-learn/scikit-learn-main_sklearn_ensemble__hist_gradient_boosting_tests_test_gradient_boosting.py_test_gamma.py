def test_gamma():
    # For a Gamma distributed target, we expect an HGBT trained with the Gamma deviance
    # (loss) to give better results than an HGBT with any other loss function, measured
    # in out-of-sample Gamma deviance as metric/score.
    # Note that squared error could potentially predict negative values which is
    # invalid (np.inf) for the Gamma deviance. A Poisson HGBT (having a log link)
    # does not have that defect.
    # Important note: It seems that a Poisson HGBT almost always has better
    # out-of-sample performance than the Gamma HGBT, measured in Gamma deviance.
    # LightGBM shows the same behaviour. Hence, we only compare to a squared error
    # HGBT, but not to a Poisson deviance HGBT.
    rng = np.random.RandomState(42)
    n_train, n_test, n_features = 500, 100, 20
    X = make_low_rank_matrix(
        n_samples=n_train + n_test,
        n_features=n_features,
        random_state=rng,
    )
    # We create a log-linear Gamma model. This gives y.min ~ 1e-2, y.max ~ 1e2
    coef = rng.uniform(low=-10, high=20, size=n_features)
    # Numpy parametrizes gamma(shape=k, scale=theta) with mean = k * theta and
    # variance = k * theta^2. We parametrize it instead with mean = exp(X @ coef)
    # and variance = dispersion * mean^2 by setting k = 1 / dispersion,
    # theta =  dispersion * mean.
    dispersion = 0.5
    y = rng.gamma(shape=1 / dispersion, scale=dispersion * np.exp(X @ coef))
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=n_test, random_state=rng
    )
    gbdt_gamma = HistGradientBoostingRegressor(loss="gamma", random_state=123)
    gbdt_mse = HistGradientBoostingRegressor(loss="squared_error", random_state=123)
    dummy = DummyRegressor(strategy="mean")
    for model in (gbdt_gamma, gbdt_mse, dummy):
        model.fit(X_train, y_train)

    for X, y in [(X_train, y_train), (X_test, y_test)]:
        loss_gbdt_gamma = mean_gamma_deviance(y, gbdt_gamma.predict(X))
        # We restrict the squared error HGBT to predict at least the minimum seen y at
        # train time to make it strictly positive.
        loss_gbdt_mse = mean_gamma_deviance(
            y, np.maximum(np.min(y_train), gbdt_mse.predict(X))
        )
        loss_dummy = mean_gamma_deviance(y, dummy.predict(X))
        assert loss_gbdt_gamma < loss_dummy
        assert loss_gbdt_gamma < loss_gbdt_mse