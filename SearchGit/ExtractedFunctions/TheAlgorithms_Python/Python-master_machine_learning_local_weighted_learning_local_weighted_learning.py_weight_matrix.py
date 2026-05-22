def weight_matrix(point: np.ndarray, x_train: np.ndarray, tau: float) -> np.ndarray:
    """
    Calculate the weight of every point in the training data around a given
    prediction point

    Args:
        point: x-value at which the prediction is being made
        x_train: ndarray of x-values for training
        tau: bandwidth value, controls how quickly the weight of training values
            decreases as the distance from the prediction point increases

    Returns:
        m x m weight matrix around the prediction point, where m is the size of
        the training set
    >>> weight_matrix(
    ...     np.array([1., 1.]),
    ...     np.array([[16.99, 10.34], [21.01,23.68], [24.59,25.69]]),
    ...     0.6
    ... )
    array([[1.43807972e-207, 0.00000000e+000, 0.00000000e+000],
           [0.00000000e+000, 0.00000000e+000, 0.00000000e+000],
           [0.00000000e+000, 0.00000000e+000, 0.00000000e+000]])
    """
    m = len(x_train)  # Number of training samples
    weights = np.eye(m)  # Initialize weights as identity matrix
    for j in range(m):
        diff = point - x_train[j]
        weights[j, j] = np.exp(diff @ diff.T / (-2.0 * tau**2))

    return weights