def test_graphical_lasso_cv(global_random_seed):
    # Sample data from a sparse multivariate normal
    dim = 5
    n_samples = 6
    random_state = np.random.RandomState(global_random_seed)
    prec = make_sparse_spd_matrix(dim, alpha=0.96, random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    # Capture stdout, to smoke test the verbose mode
    orig_stdout = sys.stdout
    try:
        sys.stdout = StringIO()
        # We need verbose very high so that Parallel prints on stdout
        GraphicalLassoCV(verbose=100, alphas=5, tol=1e-1).fit(X)
    finally:
        sys.stdout = orig_stdout