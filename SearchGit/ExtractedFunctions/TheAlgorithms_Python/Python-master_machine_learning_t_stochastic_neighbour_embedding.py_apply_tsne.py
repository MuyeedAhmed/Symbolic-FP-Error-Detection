def apply_tsne(
    data_matrix: ndarray,
    n_components: int = 2,
    learning_rate: float = 200.0,
    n_iter: int = 500,
) -> ndarray:
    """
    Apply t-SNE for dimensionality reduction.

    Args:
        data_matrix: Original dataset (features).
        n_components: Target dimension (2D or 3D).
        learning_rate: Step size for gradient descent.
        n_iter: Number of iterations.

    Returns:
        ndarray: Low-dimensional embedding of the data.

    >>> features, _ = collect_dataset()
    >>> embedding = apply_tsne(features, n_components=2, n_iter=50)
    >>> embedding.shape
    (150, 2)
    """
    if n_components < 1 or n_iter < 1:
        raise ValueError("n_components and n_iter must be >= 1")

    n_samples = data_matrix.shape[0]
    rng = np.random.default_rng()
    embedding = rng.standard_normal((n_samples, n_components)) * 1e-4

    high_dim_affinities = compute_pairwise_affinities(data_matrix)
    high_dim_affinities = np.maximum(high_dim_affinities, 1e-12)

    embedding_increment = np.zeros_like(embedding)
    momentum = 0.5

    for iteration in range(n_iter):
        low_dim_affinities, numerator_matrix = compute_low_dim_affinities(embedding)
        low_dim_affinities = np.maximum(low_dim_affinities, 1e-12)

        affinity_diff = high_dim_affinities - low_dim_affinities

        gradient = 4 * (
            np.dot((affinity_diff * numerator_matrix), embedding)
            - np.multiply(
                np.sum(affinity_diff * numerator_matrix, axis=1)[:, np.newaxis],
                embedding,
            )
        )

        embedding_increment = momentum * embedding_increment - learning_rate * gradient
        embedding += embedding_increment

        if iteration == int(n_iter / 4):
            momentum = 0.8

    return embedding