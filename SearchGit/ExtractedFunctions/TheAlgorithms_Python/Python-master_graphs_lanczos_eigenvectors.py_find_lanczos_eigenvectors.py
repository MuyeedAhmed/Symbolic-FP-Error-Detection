def find_lanczos_eigenvectors(
    graph: list[list[int | None]], num_eigenvectors: int
) -> tuple[np.ndarray, np.ndarray]:
    """Computes the largest eigenvalues and their corresponding eigenvectors using the
    Lanczos method.

    Args:
        graph: The graph as a list of adjacency lists.
        num_eigenvectors: Number of largest eigenvalues and eigenvectors to compute.

    Returns:
        A tuple containing:
            - eigenvalues: 1D array of the largest eigenvalues in descending order.
            - eigenvectors: 2D array where each column is an eigenvector corresponding
                            to an eigenvalue.

    Raises:
        ValueError: If the graph format is invalid or num_eigenvectors is out of bounds.

    >>> eigenvalues, eigenvectors = find_lanczos_eigenvectors(
    ...     [[1, 2], [0, 2], [0, 1]], 2
    ... )
    >>> len(eigenvalues) == 2 and eigenvectors.shape[1] == 2
    True
    """
    validate_adjacency_list(graph)
    tridiagonal_matrix, orthonormal_basis = lanczos_iteration(graph, num_eigenvectors)
    eigenvalues, eigenvectors = np.linalg.eigh(tridiagonal_matrix)
    return eigenvalues[::-1], np.dot(orthonormal_basis, eigenvectors[:, ::-1])