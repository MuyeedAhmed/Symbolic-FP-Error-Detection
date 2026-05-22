def test_pls_regression_scaling_coef():
    """Check that when using `scale=True`, the coefficients are using the std. dev. from
    both `X` and `y`.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/27964
    """
    # handcrafted data where we can predict y from X with an additional scaling factor
    rng = np.random.RandomState(0)
    coef = rng.uniform(size=(3, 5))
    X = rng.normal(scale=10, size=(30, 5))  # add a std of 10
    y = X @ coef.T

    # we need to make sure that the dimension of the latent space is large enough to
    # perfectly predict `y` from `X` (no information loss)
    pls = PLSRegression(n_components=5, scale=True).fit(X, y)
    assert_allclose(pls.coef_, coef)

    # we therefore should be able to predict `y` from `X`
    assert_allclose(pls.predict(X), y)