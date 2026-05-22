def test_linalg_warning_with_newton_solver(global_random_seed):
    """
    Test that the Newton solver raises a warning and falls back to LBFGS when
    encountering a singular or ill-conditioned Hessian matrix.

    This test assess the behavior of `PoissonRegressor` with the "newton-cholesky"
    solver.
    It verifies the following:-
    - The model significantly improves upon the constant baseline deviance.
    - LBFGS remains robust on collinear data.
    - The Newton solver raises a `LinAlgWarning` on collinear data and falls
      back to LBFGS.
    """
    newton_solver = "newton-cholesky"
    rng = np.random.RandomState(global_random_seed)
    # Use at least 20 samples to reduce the likelihood of getting a degenerate
    # dataset for any global_random_seed.
    X_orig = rng.normal(size=(20, 3))
    y = rng.poisson(
        np.exp(X_orig @ np.ones(X_orig.shape[1])), size=X_orig.shape[0]
    ).astype(np.float64)

    # Collinear variation of the same input features.
    X_collinear = np.hstack([X_orig] * 10)

    # Let's consider the deviance of a constant baseline on this problem.
    baseline_pred = np.full_like(y, y.mean())
    constant_model_deviance = mean_poisson_deviance(y, baseline_pred)
    assert constant_model_deviance > 1.0

    # No warning raised on well-conditioned design, even without regularization.
    tol = 1e-10
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reg = PoissonRegressor(solver=newton_solver, alpha=0.0, tol=tol).fit(X_orig, y)
    original_newton_deviance = mean_poisson_deviance(y, reg.predict(X_orig))

    # On this dataset, we should have enough data points to not make it
    # possible to get a near zero deviance (for the any of the admissible
    # random seeds). This will make it easier to interpret meaning of rtol in
    # the subsequent assertions:
    assert original_newton_deviance > 0.2

    # We check that the model could successfully fit information in X_orig to
    # improve upon the constant baseline by a large margin (when evaluated on
    # the training set).
    assert constant_model_deviance - original_newton_deviance > 0.1

    # LBFGS is robust to a collinear design because its approximation of the
    # Hessian is Symmeric Positive Definite by construction. Let's record its
    # solution
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reg = PoissonRegressor(solver="lbfgs", alpha=0.0, tol=tol).fit(X_collinear, y)
    collinear_lbfgs_deviance = mean_poisson_deviance(y, reg.predict(X_collinear))

    # The LBFGS solution on the collinear is expected to reach a comparable
    # solution to the Newton solution on the original data.
    rtol = 1e-6
    assert collinear_lbfgs_deviance == pytest.approx(original_newton_deviance, rel=rtol)

    # Fitting a Newton solver on the collinear version of the training data
    # without regularization should raise an informative warning and fallback
    # to the LBFGS solver.
    msg = (
        "The inner solver of .*Newton.*Solver stumbled upon a singular or very "
        "ill-conditioned Hessian matrix"
    )
    with pytest.warns(scipy.linalg.LinAlgWarning, match=msg):
        reg = PoissonRegressor(solver=newton_solver, alpha=0.0, tol=tol).fit(
            X_collinear, y
        )
    # As a result we should still automatically converge to a good solution.
    collinear_newton_deviance = mean_poisson_deviance(y, reg.predict(X_collinear))
    assert collinear_newton_deviance == pytest.approx(
        original_newton_deviance, rel=rtol
    )

    # Increasing the regularization slightly should make the problem go away:
    with warnings.catch_warnings():
        warnings.simplefilter("error", scipy.linalg.LinAlgWarning)
        reg = PoissonRegressor(solver=newton_solver, alpha=1e-10).fit(X_collinear, y)

    # The slightly penalized model on the collinear data should be close enough
    # to the unpenalized model on the original data.
    penalized_collinear_newton_deviance = mean_poisson_deviance(
        y, reg.predict(X_collinear)
    )
    assert penalized_collinear_newton_deviance == pytest.approx(
        original_newton_deviance, rel=rtol
    )