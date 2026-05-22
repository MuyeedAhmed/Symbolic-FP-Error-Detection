def euclidean_distance_matrix(A, B):
    """Calculate the Euclidean distance matrix.

    Args:
        A (ndarray): The point sequence.
        B (ndarray): The point sequence with the same dimensions as A.

    returns:
        D (ndarray): The Euclidean distance matrix.
    """
    assert A.ndim == 2
    assert B.ndim == 2
    assert A.shape[1] == B.shape[1]

    m = A.shape[0]
    n = B.shape[0]

    A_dots = (A * A).sum(axis=1).reshape((m, 1)) * np.ones(shape=(1, n))
    B_dots = (B * B).sum(axis=1) * np.ones(shape=(m, 1))
    D_squared = A_dots + B_dots - 2 * A.dot(B.T)

    zero_mask = np.less(D_squared, 0.0)
    D_squared[zero_mask] = 0.0
    D = np.sqrt(D_squared)
    return D