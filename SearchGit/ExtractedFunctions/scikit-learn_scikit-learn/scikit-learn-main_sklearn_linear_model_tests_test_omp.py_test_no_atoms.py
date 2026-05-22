def test_no_atoms():
    y_empty = np.zeros_like(y)
    Xy_empty = np.dot(X.T, y_empty)
    gamma_empty = ignore_warnings(orthogonal_mp)(X, y_empty, n_nonzero_coefs=1)
    gamma_empty_gram = ignore_warnings(orthogonal_mp)(G, Xy_empty, n_nonzero_coefs=1)
    assert np.all(gamma_empty == 0)
    assert np.all(gamma_empty_gram == 0)