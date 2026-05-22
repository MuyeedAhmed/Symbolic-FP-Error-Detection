def test_poisson_regressor_array_api_compliance(
    use_sample_weight,
    array_namespace,
    device_name,
    dtype_name,
    global_random_seed,
):
    xp, device = _array_api_for_tests(array_namespace, device_name, dtype_name)
    rng = np.random.default_rng(global_random_seed)
    n_samples = 1000
    n_features = 3
    X_np = rng.normal(size=(n_samples, n_features))
    beta = np.array([0.5, -0.3, 0.8])  # true coefficients
    intercept = 1.0
    mu = np.exp(X_np @ beta + intercept)  # Poisson mean with log-link.
    y_np = rng.poisson(mu)
    # Ensure that we have non-zero targets for meaningful testing:
    assert (y_np > 0).mean() > 0.1

    X_np = X_np.astype(dtype_name, copy=False)
    y_np = y_np.astype(dtype_name, copy=False)
    X_xp = xp.asarray(X_np, device=device)
    y_xp = xp.asarray(y_np, device=device)

    if use_sample_weight:
        sample_weight = (
            rng.uniform(-1, 5, size=n_samples)
            .clip(0, None)  # over-represent null weights to cover edge-cases.
            .astype(dtype_name)
        )
    else:
        sample_weight = None

    params = dict(alpha=1, solver="lbfgs", max_iter=500)
    params["tol"] = 3e-6 if dtype_name == "float32" else 1e-13
    glm_np = PoissonRegressor(**params).fit(X_np, y_np, sample_weight=sample_weight)
    assert glm_np.n_iter_ < glm_np.max_iter

    # Test that alpha was not too large for meaningful testing.
    assert np.abs(glm_np.coef_).max() > 0.1

    predict_np = glm_np.predict(X_np)
    atol = _atol_for_type(dtype_name)
    rtol = 2e-3 if dtype_name == "float32" else 3e-7

    with config_context(array_api_dispatch=True):
        glm_xp = PoissonRegressor(**params).fit(X_xp, y_xp, sample_weight=sample_weight)
        if dtype_name == "float64":
            assert abs(glm_xp.n_iter_ - glm_np.n_iter_) <= 1

        for attr_name in ("coef_", "intercept_"):
            attr_xp = getattr(glm_xp, attr_name)
            attr_np = getattr(glm_np, attr_name)
            assert_allclose(
                move_to(attr_xp, xp=np, device="cpu"), attr_np, rtol=rtol, atol=atol
            )
            assert attr_xp.dtype == X_xp.dtype
            assert array_api_device(attr_xp) == array_api_device(X_xp)

        predict_xp = glm_xp.predict(X_xp)
        assert_allclose(
            move_to(predict_xp, xp=np, device="cpu"),
            predict_np,
            rtol=rtol,
            atol=atol,
        )
        assert predict_xp.dtype == X_xp.dtype
        assert array_api_device(predict_xp) == array_api_device(X_xp)