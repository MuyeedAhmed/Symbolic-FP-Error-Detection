def local_weight(
    point: np.ndarray, x_train: np.ndarray, y_train: np.ndarray, tau: float
) -> np.ndarray:
    """
    Calculate the local weights at a given prediction point using the weight
    matrix for that point

    Args:
        point: x-value at which the prediction is being made
        x_train: ndarray of x-values for training
        y_train: ndarray of y-values for training
        tau: bandwidth value, controls how quickly the weight of training values
            decreases as the distance from the prediction point increases
    Returns:
        ndarray of local weights
    >>> local_weight(
    ...     np.array([1., 1.]),
    ...     np.array([[16.99, 10.34], [21.01,23.68], [24.59,25.69]]),
    ...     np.array([[1.01, 1.66, 3.5]]),
    ...     0.6
    ... )
    array([[0.00873174],
           [0.08272556]])
    """
    weight_mat = weight_matrix(point, x_train, tau)
    weight = np.linalg.inv(x_train.T @ weight_mat @ x_train) @ (
        x_train.T @ weight_mat @ y_train.T
    )

    return weight