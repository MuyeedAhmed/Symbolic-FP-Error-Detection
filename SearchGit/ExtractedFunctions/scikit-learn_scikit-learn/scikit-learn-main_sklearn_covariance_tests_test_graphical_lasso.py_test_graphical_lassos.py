def test_graphical_lassos(global_random_seed):
    """Test the graphical lasso solvers."""
    # Sample data from a sparse multivariate normal
    dim = 10
    n_samples = 100
    random_state = check_random_state(global_random_seed)
    prec = make_sparse_spd_matrix(dim, alpha=0.95, random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    emp_cov = empirical_covariance(X)

    for alpha in (0.0, 0.1, 0.25):
        covs = dict()
        icovs = dict()
        for method in ("cd", "lars"):
            cov_, icov_, costs = graphical_lasso(
                emp_cov,
                return_costs=True,
                alpha=alpha,
                mode=method,
                tol=1e-7,
                enet_tol=1e-11,
                max_iter=100,
            )
            covs[method] = cov_
            icovs[method] = icov_
            costs, dual_gap = np.array(costs).T
            # Check that the costs always decrease (doesn't hold if alpha == 0)
            if not alpha == 0:
                # use 1e-10 since the cost can be exactly 0
                assert_array_less(np.diff(costs), 1e-10)
        # Check that the 2 approaches give similar results
        assert_allclose(covs["cd"], covs["lars"], atol=2e-3)
        assert_allclose(icovs["cd"], icovs["lars"], atol=2e-3)

    # Smoke test the estimator
    model = GraphicalLasso(alpha=0.25, tol=1e-7, enet_tol=1e-11, max_iter=100).fit(X)
    model.score(X)
    assert_allclose(model.covariance_, covs["cd"], rtol=1e-6)

    # For a centered matrix, assume_centered could be chosen True or False
    # Check that this returns indeed the same result for centered data
    Z = X - X.mean(0)
    precs = list()
    for assume_centered in (False, True):
        prec_ = GraphicalLasso(assume_centered=assume_centered).fit(Z).precision_
        precs.append(prec_)
    assert_array_almost_equal(precs[0], precs[1])