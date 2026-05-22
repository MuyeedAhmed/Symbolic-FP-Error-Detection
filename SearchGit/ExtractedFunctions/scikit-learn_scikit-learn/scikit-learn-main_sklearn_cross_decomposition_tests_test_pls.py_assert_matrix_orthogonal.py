def assert_matrix_orthogonal(M):
    K = np.dot(M.T, M)
    assert_array_almost_equal(K, np.diag(np.diag(K)))