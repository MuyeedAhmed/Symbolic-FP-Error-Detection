def _check_precomputed_gram_matrix(
    X, precompute, X_offset, X_scale, rtol=None, atol=1e-5
):
    """Computes a single element of the gram matrix and compares it to
    the corresponding element of the user supplied gram matrix.

    If the values do not match a ValueError will be thrown.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Data array.

    precompute : array-like of shape (n_features, n_features)
        User-supplied gram matrix.

    X_offset : ndarray of shape (n_features,)
        Array of feature means used to center design matrix.

    X_scale : ndarray of shape (n_features,)
        Array of feature scale factors used to normalize design matrix.

    rtol : float, default=None
        Relative tolerance; see numpy.allclose
        If None, it is set to 1e-4 for arrays of dtype numpy.float32 and 1e-7
        otherwise.

    atol : float, default=1e-5
        absolute tolerance; see :func`numpy.allclose`. Note that the default
        here is more tolerant than the default for
        :func:`numpy.testing.assert_allclose`, where `atol=0`.

    Raises
    ------
    ValueError
        Raised when the provided Gram matrix is not consistent.
    """

    n_features = X.shape[1]
    f1 = n_features // 2
    f2 = min(f1 + 1, n_features - 1)

    v1 = (X[:, f1] - X_offset[f1]) * X_scale[f1]
    v2 = (X[:, f2] - X_offset[f2]) * X_scale[f2]

    expected = np.dot(v1, v2)
    actual = precompute[f1, f2]

    dtypes = [precompute.dtype, expected.dtype]
    if rtol is None:
        rtols = [1e-4 if dtype == np.float32 else 1e-7 for dtype in dtypes]
        rtol = max(rtols)

    if not np.isclose(expected, actual, rtol=rtol, atol=atol):
        raise ValueError(
            "Gram matrix passed in via 'precompute' parameter "
            "did not pass validation when a single element was "
            "checked - please check that it was computed "
            f"properly. For element ({f1},{f2}) we computed "
            f"{expected} but the user-supplied value was "
            f"{actual}."
        )