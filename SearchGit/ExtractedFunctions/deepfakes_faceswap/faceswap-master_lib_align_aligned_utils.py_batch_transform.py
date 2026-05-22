def batch_transform(matrices: npt.NDArray[np.float32],
                    points: npt.NDArray[np.float32],
                    in_place: bool = False) -> npt.NDArray[np.float32]:
    """Batch transform an array of (N, M, 2) points by the given (N, 3, 3) affine matrices

    Parameters
    ----------
    matrices
        The matrices to use to transform the points
    points
        The points to be transformed
    in_place
        ``True`` to directly transform the given points in place. ``False`` to return a new array

    Returns
    -------
    The transformed points
    """
    retval = points if in_place else np.empty_like(points)
    linear = matrices[:, :2, :2]
    translation = matrices[:, :2, 2]
    retval[:] = points @ linear.transpose(0, 2, 1) + translation[:, None, :]
    return retval