def test_gs(global_random_seed):
    # Test gram schmidt orthonormalization
    # generate a random orthogonal  matrix
    rng = np.random.RandomState(global_random_seed)
    W, _, _ = np.linalg.svd(rng.randn(10, 10))
    w = rng.randn(10)
    _gs_decorrelation(w, W, 10)
    assert (w**2).sum() < 1.0e-10
    w = rng.randn(10)
    u = _gs_decorrelation(w, W, 5)
    tmp = np.dot(u, W.T)
    assert (tmp[:5] ** 2).sum() < 1.0e-10