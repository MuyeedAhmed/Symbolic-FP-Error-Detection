def test_quantile_asymmetric_error(quantile):
    """Test quantile regression for asymmetric distributed targets."""
    n_samples = 10_000
    rng = np.random.RandomState(42)
    # take care that X @ coef + intercept > 0
    X = np.concatenate(
        (
            np.abs(rng.randn(n_samples)[:, None]),
            -rng.randint(2, size=(n_samples, 1)),
        ),
        axis=1,
    )
    intercept = 1.23
    coef = np.array([0.5, -2])
    # For an exponential distribution with rate lambda, e.g. exp(-lambda * x),
    # the quantile at level q is:
    #   quantile(q) = - log(1 - q) / lambda
    #   scale = 1/lambda = -quantile(q) / log(1-q)
    y = rng.exponential(
        scale=-(X @ coef + intercept) / np.log(1 - quantile), size=n_samples
    )
    model = HistGradientBoostingRegressor(
        loss="quantile",
        quantile=quantile,
        max_iter=25,
        random_state=0,
        max_leaf_nodes=10,
    ).fit(X, y)
    assert_allclose(np.mean(model.predict(X) > y), quantile, rtol=1e-2)

    pinball_loss = PinballLoss(quantile=quantile)
    loss_true_quantile = pinball_loss(y, X @ coef + intercept)
    loss_pred_quantile = pinball_loss(y, model.predict(X))
    # we are overfitting
    assert loss_pred_quantile <= loss_true_quantile