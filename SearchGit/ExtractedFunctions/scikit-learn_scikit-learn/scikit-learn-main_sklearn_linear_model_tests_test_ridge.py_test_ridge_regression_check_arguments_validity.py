def test_ridge_regression_check_arguments_validity(
    return_intercept, sample_weight, container, solver
):
    """check if all combinations of arguments give valid estimations"""

    # test excludes 'svd' solver because it raises exception for sparse inputs

    rng = check_random_state(42)
    X = rng.rand(1000, 3)
    true_coefs = [1, 2, 0.1]
    y = np.dot(X, true_coefs)
    true_intercept = 0.0
    if return_intercept:
        true_intercept = 10000.0
    y += true_intercept
    X_testing = container(X)

    alpha, tol = 1e-3, 1e-6
    atol = 1e-3 if _IS_32BIT else 1e-4

    positive = solver == "lbfgs"

    if solver not in ["sag", "auto"] and return_intercept:
        with pytest.raises(ValueError, match="In Ridge, only 'sag' solver"):
            ridge_regression(
                X_testing,
                y,
                alpha=alpha,
                solver=solver,
                sample_weight=sample_weight,
                return_intercept=return_intercept,
                positive=positive,
                tol=tol,
                random_state=rng,
            )
        return

    out = ridge_regression(
        X_testing,
        y,
        alpha=alpha,
        solver=solver,
        sample_weight=sample_weight,
        positive=positive,
        return_intercept=return_intercept,
        tol=tol,
        random_state=rng,
    )

    if return_intercept:
        coef, intercept = out
        assert_allclose(coef, true_coefs, rtol=0, atol=atol)
        assert_allclose(intercept, true_intercept, rtol=0, atol=atol)
    else:
        assert_allclose(out, true_coefs, rtol=0, atol=atol)