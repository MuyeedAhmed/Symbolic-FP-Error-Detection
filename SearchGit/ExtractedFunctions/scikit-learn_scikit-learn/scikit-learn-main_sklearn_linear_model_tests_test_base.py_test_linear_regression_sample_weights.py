def test_linear_regression_sample_weights(
    sparse_container, fit_intercept, global_random_seed
):
    rng = np.random.RandomState(global_random_seed)

    # It would not work with under-determined systems
    n_samples, n_features = 6, 5

    X = rng.normal(size=(n_samples, n_features))
    if sparse_container is not None:
        X = sparse_container(X)
    y = rng.normal(size=n_samples)

    sample_weight = 1.0 + rng.uniform(size=n_samples)

    # LinearRegression with explicit sample_weight
    reg = LinearRegression(fit_intercept=fit_intercept, tol=1e-16)
    reg.fit(X, y, sample_weight=sample_weight)
    coefs1 = reg.coef_
    inter1 = reg.intercept_

    assert reg.coef_.shape == (X.shape[1],)  # sanity checks

    # Closed form of the weighted least square
    # theta = (X^T W X)^(-1) @ X^T W y
    W = np.diag(sample_weight)
    X_aug = X if not fit_intercept else add_dummy_feature(X)

    Xw = X_aug.T @ W @ X_aug
    yw = X_aug.T @ W @ y
    coefs2 = linalg.solve(Xw, yw)

    if not fit_intercept:
        assert_allclose(coefs1, coefs2)
    else:
        assert_allclose(coefs1, coefs2[1:])
        assert_allclose(inter1, coefs2[0])