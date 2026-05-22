def test_pls_prediction(PLSEstimator, scale):
    """Check the behaviour of the prediction function."""
    d = load_linnerud()
    X = d.data
    y = d.target

    pls = PLSEstimator(copy=True, scale=scale).fit(X, y)
    y_pred = pls.predict(X, copy=True)

    y_mean = y.mean(axis=0)
    X_trans = X - X.mean(axis=0)

    assert_allclose(pls.intercept_, y_mean)
    assert_allclose(y_pred, X_trans @ pls.coef_.T + pls.intercept_)