def compute_low_dim_affinities(embedding_matrix: ndarray) -> tuple[ndarray, ndarray]:
    """
    Compute low-dimensional affinities (Q matrix) using a Student-t distribution.

    Args:
        embedding_matrix: Low-dimensional embedding of shape (n_samples, n_components).

    Returns:
        tuple[ndarray, ndarray]: (Q probability matrix, numerator matrix).

    >>> y = np.array([[0.0, 0.0], [1.0, 0.0]])
    >>> q_matrix, numerators = compute_low_dim_affinities(y)
    >>> q_matrix.shape
    (2, 2)
    """
    squared_sum = np.sum(np.square(embedding_matrix), axis=1)
    numerator_matrix = 1 / (
        1
        + np.add(
            np.add(-2 * np.dot(embedding_matrix, embedding_matrix.T), squared_sum).T,
            squared_sum,
        )
    )
    np.fill_diagonal(numerator_matrix, 0)

    q_matrix = numerator_matrix / np.sum(numerator_matrix)
    return q_matrix, numerator_matrix