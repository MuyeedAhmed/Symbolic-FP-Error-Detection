def local_weight_regression(
    x_train: np.ndarray, y_train: np.ndarray, tau: float
) -> np.ndarray:
    """
    Calculate predictions for each point in the training data

    Args:
        x_train: ndarray of x-values for training
        y_train: ndarray of y-values for training
        tau: bandwidth value, controls how quickly the weight of training values
            decreases as the distance from the prediction point increases

    Returns:
        ndarray of predictions
    >>> local_weight_regression(
    ...     np.array([[16.99, 10.34], [21.01, 23.68], [24.59, 25.69]]),
    ...     np.array([[1.01, 1.66, 3.5]]),
    ...     0.6
    ... )
    array([1.07173261, 1.65970737, 3.50160179])
    """
    y_pred = np.zeros(len(x_train))  # Initialize array of predictions
    for i, item in enumerate(x_train):
        y_pred[i] = np.dot(item, local_weight(item, x_train, y_train, tau)).item()

    return y_pred