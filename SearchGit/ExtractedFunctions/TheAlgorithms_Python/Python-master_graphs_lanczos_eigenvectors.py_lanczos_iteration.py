def lanczos_iteration(
    graph: list[list[int | None]], num_eigenvectors: int
) -> tuple[np.ndarray, np.ndarray]:
    """Constructs the tridiagonal matrix and orthonormal basis vectors using the
    Lanczos method.

    Args:
        graph: The graph represented as a list of adjacency lists.
        num_eigenvectors: The number of largest eigenvalues and eigenvectors
                          to approximate.

    Returns:
        A tuple containing:
            - tridiagonal_matrix: A (num_eigenvectors x num_eigenvectors) symmetric
                                  matrix.
            - orthonormal_basis: A (num_nodes x num_eigenvectors) matrix of orthonormal
                                 basis vectors.

    Raises:
        ValueError: If num_eigenvectors is less than 1 or greater than the number of
                    nodes.

    >>> graph = [[1, 2], [0, 2], [0, 1]]
    >>> T, Q = lanczos_iteration(graph, 2)
    >>> T.shape == (2, 2) and Q.shape == (3, 2)
    True
    """
    num_nodes: int = len(graph)
    if not (1 <= num_eigenvectors <= num_nodes):
        raise ValueError(
            "Number of eigenvectors must be between 1 and the number of "
            "nodes in the graph."
        )

    orthonormal_basis: np.ndarray = np.zeros((num_nodes, num_eigenvectors))
    tridiagonal_matrix: np.ndarray = np.zeros((num_eigenvectors, num_eigenvectors))

    rng = np.random.default_rng()
    initial_vector: np.ndarray = rng.random(num_nodes)
    initial_vector /= np.sqrt(np.dot(initial_vector, initial_vector))
    orthonormal_basis[:, 0] = initial_vector

    prev_beta: float = 0.0
    for iter_index in range(num_eigenvectors):
        result_vector: np.ndarray = multiply_matrix_vector(
            graph, orthonormal_basis[:, iter_index]
        )
        if iter_index > 0:
            result_vector -= prev_beta * orthonormal_basis[:, iter_index - 1]
        alpha_value: float = np.dot(orthonormal_basis[:, iter_index], result_vector)
        result_vector -= alpha_value * orthonormal_basis[:, iter_index]

        prev_beta = np.sqrt(np.dot(result_vector, result_vector))
        if iter_index < num_eigenvectors - 1 and prev_beta > 1e-10:
            orthonormal_basis[:, iter_index + 1] = result_vector / prev_beta
        tridiagonal_matrix[iter_index, iter_index] = alpha_value
        if iter_index < num_eigenvectors - 1:
            tridiagonal_matrix[iter_index, iter_index + 1] = prev_beta
            tridiagonal_matrix[iter_index + 1, iter_index] = prev_beta
    return tridiagonal_matrix, orthonormal_basis