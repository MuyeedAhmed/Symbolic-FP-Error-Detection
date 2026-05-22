def test_huber_better_r2_score():
    # Test that huber returns a better r2 score than non-outliers"""
    X, y = make_regression_with_outliers()
    huber = HuberRegressor(alpha=0.01)
    huber.fit(X, y)
    linear_loss = np.dot(X, huber.coef_) + huber.intercept_ - y
    mask = np.abs(linear_loss) < huber.epsilon * huber.scale_
    huber_score = huber.score(X[mask], y[mask])
    huber_outlier_score = huber.score(X[~mask], y[~mask])

    # The Ridge regressor should be influenced by the outliers and hence
    # give a worse score on the non-outliers as compared to the huber
    # regressor.
    ridge = Ridge(alpha=0.01)
    ridge.fit(X, y)
    ridge_score = ridge.score(X[mask], y[mask])
    ridge_outlier_score = ridge.score(X[~mask], y[~mask])
    assert huber_score > ridge_score

    # The huber model should also fit poorly on the outliers.
    assert ridge_outlier_score > huber_outlier_score