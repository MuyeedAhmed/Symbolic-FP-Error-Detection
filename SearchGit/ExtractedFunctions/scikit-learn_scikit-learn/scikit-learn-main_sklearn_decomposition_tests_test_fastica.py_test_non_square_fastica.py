def test_non_square_fastica(global_random_seed, add_noise):
    # Test the FastICA algorithm on very simple data.
    rng = np.random.RandomState(global_random_seed)

    n_samples = 1000
    # Generate two sources:
    t = np.linspace(0, 100, n_samples)
    s1 = np.sin(t)
    s2 = np.ceil(np.sin(np.pi * t))
    s = np.c_[s1, s2].T
    center_and_norm(s)
    s1, s2 = s

    # Mixing matrix
    mixing = rng.randn(6, 2)
    m = np.dot(mixing, s)

    if add_noise:
        m += 0.1 * rng.randn(6, n_samples)

    center_and_norm(m)

    k_, mixing_, s_ = fastica(
        m.T, n_components=2, whiten="unit-variance", random_state=rng
    )
    s_ = s_.T

    # Check that the mixing model described in the docstring holds:
    assert_allclose(s_, np.dot(np.dot(mixing_, k_), m))

    center_and_norm(s_)
    s1_, s2_ = s_
    # Check to see if the sources have been estimated
    # in the wrong order
    if abs(np.dot(s1_, s2)) > abs(np.dot(s1_, s1)):
        s2_, s1_ = s_
    s1_ *= np.sign(np.dot(s1_, s1))
    s2_ *= np.sign(np.dot(s2_, s2))

    # Check that we have estimated the original sources
    if not add_noise:
        assert_allclose(np.dot(s1_, s1) / n_samples, 1, atol=1e-3)
        assert_allclose(np.dot(s2_, s2) / n_samples, 1, atol=1e-3)