def test_fastica_eigh_low_rank_warning(global_random_seed):
    """Test FastICA eigh solver raises warning for low-rank data."""
    rng = np.random.RandomState(global_random_seed)
    A = rng.randn(10, 2)
    X = A @ A.T
    ica = FastICA(random_state=0, whiten="unit-variance", whiten_solver="eigh")
    msg = "There are some small singular values"

    with pytest.warns(UserWarning, match=msg):
        with ignore_warnings(category=ConvergenceWarning):
            # The FastICA solver may not converge for some data with specific
            # random seeds but this happens after the whiten step so this is
            # not want we want to test here.
            ica.fit(X)