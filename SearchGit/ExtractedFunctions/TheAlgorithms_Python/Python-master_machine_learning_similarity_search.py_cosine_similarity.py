def cosine_similarity(input_a: np.ndarray, input_b: np.ndarray) -> float:
    """
    Calculates cosine similarity between two data.
    :param input_a: ndarray of first vector.
    :param input_b: ndarray of second vector.
    :return: Cosine similarity of input_a and input_b. By using math.sqrt(),
             result will be float.

    >>> cosine_similarity(np.array([1]), np.array([1]))
    1.0
    >>> cosine_similarity(np.array([1, 2]), np.array([6, 32]))
    0.9615239476408232
    """
    return float(np.dot(input_a, input_b) / (norm(input_a) * norm(input_b)))