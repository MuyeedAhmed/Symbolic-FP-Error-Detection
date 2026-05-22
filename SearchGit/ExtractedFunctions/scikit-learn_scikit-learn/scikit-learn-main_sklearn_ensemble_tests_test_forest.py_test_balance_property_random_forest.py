def test_balance_property_random_forest(criterion):
    """ "Test that sum(y_pred)==sum(y_true) on the training set."""
    rng = np.random.RandomState(42)
    n_train, n_test, n_features = 500, 500, 10
    X = datasets.make_low_rank_matrix(
        n_samples=n_train + n_test, n_features=n_features, random_state=rng
    )

    coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)
    y = rng.poisson(lam=np.exp(X @ coef))

    reg = RandomForestRegressor(
        criterion=criterion, n_estimators=10, bootstrap=False, random_state=rng
    )
    reg.fit(X, y)

    assert np.sum(reg.predict(X)) == pytest.approx(np.sum(y))