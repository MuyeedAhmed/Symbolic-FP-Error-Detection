def cosine_similarity(X: Matrix, Y: Matrix) -> np.ndarray:  # noqa: N803
    """Row-wise cosine similarity between two equal-width matrices."""
    if len(X) == 0 or len(Y) == 0:
        return np.array([])

    x: np.ndarray = np.array(X)
    y: np.ndarray = np.array(Y)
    if x.shape[1] != y.shape[1]:
        msg = (
            f"Number of columns in X and Y must be the same. X has shape {x.shape} "
            f"and Y has shape {y.shape}."
        )
        raise ValueError(msg)
    try:
        import simsimd as simd  # noqa: PLC0415

        x = np.array(x, dtype=np.float32)
        y = np.array(y, dtype=np.float32)
        return 1 - np.array(simd.cdist(x, y, metric="cosine"))
    except ImportError:
        x_norm = np.linalg.norm(x, axis=1)
        y_norm = np.linalg.norm(y, axis=1)
        # Ignore divide by zero errors run time warnings as those are handled below.
        with np.errstate(divide="ignore", invalid="ignore"):
            similarity = np.dot(x, y.T) / np.outer(x_norm, y_norm)
        similarity[np.isnan(similarity) | np.isinf(similarity)] = 0.0
        return similarity