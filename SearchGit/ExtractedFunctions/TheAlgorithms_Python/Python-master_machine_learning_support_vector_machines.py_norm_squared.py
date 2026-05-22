def norm_squared(vector: ndarray) -> float:
    """
    Return the squared second norm of vector
    norm_squared(v) = sum(x * x for x in v)

    Args:
        vector (ndarray): input vector

    Returns:
        float: squared second norm of vector

    >>> int(norm_squared([1, 2]))
    5
    >>> int(norm_squared(np.asarray([1, 2])))
    5
    >>> int(norm_squared([0, 0]))
    0
    """
    return np.dot(vector, vector)