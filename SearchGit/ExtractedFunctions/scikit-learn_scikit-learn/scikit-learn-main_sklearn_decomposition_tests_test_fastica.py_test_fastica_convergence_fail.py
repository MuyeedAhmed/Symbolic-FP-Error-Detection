def test_fastica_convergence_fail(global_random_seed):
    # Test the FastICA algorithm on very simple data
    # (see test_non_square_fastica).
    # Ensure a ConvergenceWarning raised if the tolerance is sufficiently low.
    rng = np.random.RandomState(global_random_seed)

    n_samples = 1000
    # Generate two sources:
    t = np.linspace(0, 100, n_samples)
    s1 = np.sin(t)
    s2 = np.ceil(np.sin(np.pi * t))
    s = np.c_[s1, s2].T
    center_and_norm(s)

    # Mixing matrix
    mixing = rng.randn(6, 2)
    m = np.dot(mixing, s)

    # Do fastICA with tolerance 0. to ensure failing convergence
    warn_msg = (
        "FastICA did not converge. Consider increasing tolerance "
        "or the maximum number of iterations."
    )
    with pytest.warns(ConvergenceWarning, match=warn_msg):
        ica = FastICA(
            algorithm="parallel", n_components=2, random_state=rng, max_iter=2, tol=0.0
        )
        ica.fit(m.T)