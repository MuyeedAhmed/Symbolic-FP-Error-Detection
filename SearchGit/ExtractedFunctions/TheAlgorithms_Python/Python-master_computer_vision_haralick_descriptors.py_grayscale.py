def grayscale(image: np.ndarray) -> np.ndarray:
    """
    Uses luminance weights to transform RGB channel to greyscale, by
    taking the dot product between the channel and the weights.

    Example:
        >>> grayscale(np.array([[[108, 201, 72], [255, 11,  127]],
        ...                     [[56,  56,  56], [128, 255, 107]]]))
        array([[158,  97],
               [ 56, 200]], dtype=uint8)
    """
    return np.dot(image[:, :, 0:3], [0.299, 0.587, 0.114]).astype(np.uint8)