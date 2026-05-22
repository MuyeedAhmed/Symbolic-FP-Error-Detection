def test_no_path_all_precomputed():
    # Test that the ``return_path=False`` option with Gram and Xy remains
    # correct
    X, y = 3 * diabetes.data, diabetes.target
    G = np.dot(X.T, X)
    Xy = np.dot(X.T, y)
    alphas_, _, coef_path_ = linear_model.lars_path(
        X, y, method="lasso", Xy=Xy, Gram=G, alpha_min=0.9
    )
    alpha_, _, coef = linear_model.lars_path(
        X, y, method="lasso", Gram=G, Xy=Xy, alpha_min=0.9, return_path=False
    )

    assert_array_almost_equal(coef, coef_path_[:, -1])
    assert alpha_ == alphas_[-1]