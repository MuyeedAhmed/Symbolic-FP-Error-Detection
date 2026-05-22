def test_collinearity():
    # Check that lars_path is robust to collinearity in input
    X = np.array([[3.0, 3.0, 1.0], [2.0, 2.0, 0.0], [1.0, 1.0, 0]])
    y = np.array([1.0, 0.0, 0])
    rng = np.random.RandomState(0)

    f = ignore_warnings
    _, _, coef_path_ = f(linear_model.lars_path)(X, y, alpha_min=0.01)
    assert not np.isnan(coef_path_).any()
    residual = np.dot(X, coef_path_[:, -1]) - y
    assert (residual**2).sum() < 1.0  # just make sure it's bounded

    n_samples = 10
    X = rng.rand(n_samples, 5)
    y = np.zeros(n_samples)
    _, _, coef_path_ = linear_model.lars_path(
        X,
        y,
        Gram="auto",
        copy_X=False,
        copy_Gram=False,
        alpha_min=0.0,
        method="lasso",
        verbose=0,
        max_iter=500,
    )
    assert_array_almost_equal(coef_path_, np.zeros_like(coef_path_))