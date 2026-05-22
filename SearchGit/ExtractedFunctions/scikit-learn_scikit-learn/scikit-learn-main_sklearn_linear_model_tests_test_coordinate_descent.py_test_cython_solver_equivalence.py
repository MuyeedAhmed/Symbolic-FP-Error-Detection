def test_cython_solver_equivalence(sparse_csc_type):
    """Test that all 3 Cython solvers for 1-d targets give same results."""
    X, y = make_regression()
    X_mean = X.mean(axis=0)
    X_centered = np.asfortranarray(X - X_mean)
    y -= y.mean()
    alpha_max = np.linalg.norm(X.T @ y, ord=np.inf)
    alpha = alpha_max / 10
    params = {
        "beta": 0,
        "max_iter": 100,
        "tol": 1e-10,
        "rng": np.random.RandomState(0),  # not used, but needed as argument
        "random": False,
        "positive": False,
    }

    def zc():
        """Create a new zero coefficient array (zc)."""
        return np.zeros(X.shape[1])

    # For alpha_max, coefficients must all be zero.
    coef_1 = zc()
    for do_screening in [True, False]:
        cd_fast.enet_coordinate_descent(
            w=coef_1,
            alpha=alpha_max,
            X=X_centered,
            y=y,
            **params,
            do_screening=do_screening,
        )
        assert_allclose(coef_1, 0)

    # Without gap safe screening rules
    coef_1 = zc()
    cd_fast.enet_coordinate_descent(
        w=coef_1, alpha=alpha, X=X_centered, y=y, **params, do_screening=False
    )
    # At least 2 coefficients are non-zero
    assert 2 <= np.sum(np.abs(coef_1) > 1e-8) < X.shape[1]

    # With gap safe screening rules
    coef_2 = zc()
    cd_fast.enet_coordinate_descent(
        w=coef_2, alpha=alpha, X=X_centered, y=y, **params, do_screening=True
    )
    assert_allclose(coef_2, coef_1)

    # Sparse
    Xs = sparse_csc_type(X)
    for do_screening in [True, False]:
        coef_3 = zc()
        cd_fast.sparse_enet_coordinate_descent(
            w=coef_3,
            alpha=alpha,
            X_data=Xs.data,
            X_indices=Xs.indices,
            X_indptr=Xs.indptr,
            y=y,
            sample_weight=None,
            X_mean=X_mean,
            **params,
            do_screening=do_screening,
        )
        assert_allclose(coef_3, coef_1)

    # Gram
    for do_screening in [True, False]:
        coef_4 = zc()
        cd_fast.enet_coordinate_descent_gram(
            w=coef_4,
            alpha=alpha,
            Q=X_centered.T @ X_centered,
            q=X_centered.T @ y,
            y=y,
            **params,
            do_screening=do_screening,
        )
        assert_allclose(coef_4, coef_1)