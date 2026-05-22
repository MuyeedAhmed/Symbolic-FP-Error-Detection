def normalize_adjacent_matrix(A):
    assert A.ndim == 2
    assert A.shape[0] == A.shape[1]

    A = A + np.eye(A.shape[0])
    d = np.sum(A, axis=0)
    d = np.clip(d, 0, None)
    d_inv = np.power(d, -0.5).flatten()
    d_inv[np.isinf(d_inv)] = 0.0
    d_inv = np.diag(d_inv)
    G = A.dot(d_inv).transpose().dot(d_inv)
    return G