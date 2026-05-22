def test_loss_grad_hess_are_the_same(
    base_loss,
    fit_intercept,
    sample_weight,
    l2_reg_strength,
    csr_container,
    global_random_seed,
):
    """Test that loss and gradient are the same across different functions."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=10, n_features=5, seed=global_random_seed
    )
    X_old, y_old, coef_old = X.copy(), y.copy(), coef.copy()

    if sample_weight == "range":
        sample_weight = np.linspace(1, y.shape[0], num=y.shape[0])

    l1 = loss.loss(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g1 = loss.gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l2, g2 = loss.loss_gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g3, h3 = loss.gradient_hessian_product(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g4, h4, _ = loss.gradient_hessian(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    assert_allclose(l1, l2)
    assert_allclose(g1, g2)
    assert_allclose(g1, g3)
    assert_allclose(g1, g4)
    # The ravelling only takes effect for multiclass.
    assert_allclose(h4 @ g4.ravel(order="F"), h3(g3).ravel(order="F"))
    # Test that gradient_out and hessian_out are considered properly.
    g_out = np.empty_like(coef)
    h_out = np.empty_like(coef, shape=(coef.size, coef.size))
    g5, h5, _ = loss.gradient_hessian(
        coef,
        X,
        y,
        sample_weight=sample_weight,
        l2_reg_strength=l2_reg_strength,
        gradient_out=g_out,
        hessian_out=h_out,
    )
    assert np.shares_memory(g5, g_out)
    assert np.shares_memory(h5, h_out)
    assert_allclose(g5, g_out)
    assert_allclose(h5, h_out)
    assert_allclose(g1, g5)
    assert_allclose(h5, h4)

    # same for sparse X
    Xs = csr_container(X)
    l1_sp = loss.loss(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g1_sp = loss.gradient(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l2_sp, g2_sp = loss.loss_gradient(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g3_sp, h3_sp = loss.gradient_hessian_product(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g4_sp, h4_sp, _ = loss.gradient_hessian(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    assert_allclose(l1, l1_sp)
    assert_allclose(l1, l2_sp)
    assert_allclose(g1, g1_sp)
    assert_allclose(g1, g2_sp)
    assert_allclose(g1, g3_sp)
    assert_allclose(h3(g1), h3_sp(g1_sp))
    assert_allclose(g1, g4_sp)
    assert_allclose(h4, h4_sp)

    # X, y and coef should not have changed
    assert_allclose(X, X_old)
    assert_allclose(Xs.toarray(), X_old)
    assert_allclose(y, y_old)
    assert_allclose(coef, coef_old)