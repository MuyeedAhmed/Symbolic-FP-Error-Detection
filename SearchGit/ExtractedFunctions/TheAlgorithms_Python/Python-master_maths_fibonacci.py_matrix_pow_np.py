def matrix_pow_np(m: ndarray, power: int) -> ndarray:
    """
    Raises a matrix to the power of 'power' using binary exponentiation.

    Args:
        m: Matrix as a numpy array.
        power: The power to which the matrix is to be raised.

    Returns:
        The matrix raised to the power.

    Raises:
        ValueError: If power is negative.

    >>> m = np.array([[1, 1], [1, 0]], dtype=int)
    >>> matrix_pow_np(m, 0)  # Identity matrix when raised to the power of 0
    array([[1, 0],
           [0, 1]])

    >>> matrix_pow_np(m, 1)  # Same matrix when raised to the power of 1
    array([[1, 1],
           [1, 0]])

    >>> matrix_pow_np(m, 5)
    array([[8, 5],
           [5, 3]])

    >>> matrix_pow_np(m, -1)
    Traceback (most recent call last):
        ...
    ValueError: power is negative
    """
    result = np.array([[1, 0], [0, 1]], dtype=int)  # Identity Matrix
    base = m
    if power < 0:  # Negative power is not allowed
        raise ValueError("power is negative")
    while power:
        if power % 2 == 1:
            result = np.dot(result, base)
        base = np.dot(base, base)
        power //= 2
    return result