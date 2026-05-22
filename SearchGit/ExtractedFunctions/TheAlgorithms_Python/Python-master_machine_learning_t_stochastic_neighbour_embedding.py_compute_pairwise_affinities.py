def compute_pairwise_affinities(data_matrix: ndarray, sigma: float = 1.0) -> ndarray:
    """
    Compute high-dimensional affinities (P matrix) using a Gaussian kernel.

    Args:
        data_matrix: Input data of shape (n_samples, n_features).
        sigma: Gaussian kernel bandwidth.

    Returns:
        ndarray: Symmetrized probability matrix.

    >>> x = np.array([[0.0, 0.0], [1.0, 0.0]])
    >>> probabilities = compute_pairwise_affinities(x)
    >>> float(round(probabilities[0, 1], 3))
    0.25
    """
    n_samples = data_matrix.shape[0]
    squared_sum = np.sum(np.square(data_matrix), axis=1)
    squared_distance = np.add(
        np.add(-2 * np.dot(data_matrix, data_matrix.T), squared_sum).T, squared_sum
    )

    affinity_matrix = np.exp(-squared_distance / (2 * sigma**2))
    np.fill_diagonal(affinity_matrix, 0)

    affinity_matrix /= np.sum(affinity_matrix)
    return (affinity_matrix + affinity_matrix.T) / (2 * n_samples)