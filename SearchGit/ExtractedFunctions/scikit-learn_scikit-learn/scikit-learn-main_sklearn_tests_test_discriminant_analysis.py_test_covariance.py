def test_covariance():
    x, y = make_blobs(n_samples=100, n_features=5, centers=1, random_state=42)

    # make features correlated
    x = np.dot(x, np.arange(x.shape[1] ** 2).reshape(x.shape[1], x.shape[1]))

    c_e = _cov(x, "empirical")
    assert_almost_equal(c_e, c_e.T)

    c_s = _cov(x, "auto")
    assert_almost_equal(c_s, c_s.T)