def _randomized_range_finder(
    A, *, size, n_iter, power_iteration_normalizer="auto", random_state=None
):
    """Body of randomized_range_finder without input validation."""
    xp, is_array_api_compliant = get_namespace(A)
    random_state = check_random_state(random_state)

    # Generating normal random vectors with shape: (A.shape[1], size)
    # XXX: generate random number directly from xp if it's possible
    # one day.
    Q = random_state.normal(size=(A.shape[1], size))
    if A.dtype == xp.float32 or (
        is_array_api_compliant
        and _max_precision_float_dtype(xp, device=device(A)) == xp.float32
    ):
        # Use float32 computation and components if A has a float32 dtype
        # or if A has integer dtype and device doesn't not support float64.

        # Downcast while Q is still a NumPy array to avoid allocating float64
        # on devices that do not support it. The Array API does not require
        # xp.asarray(..., dtype=..., device=device) to accept such a downcast.
        Q = Q.astype(np.float32, copy=False)

    if is_array_api_compliant:
        Q = xp.asarray(Q, device=device(A))
    else:
        Q = xp.asarray(Q)

    # Deal with "auto" mode
    if power_iteration_normalizer == "auto":
        if n_iter <= 2:
            power_iteration_normalizer = "none"
        elif is_array_api_compliant:
            # XXX: https://github.com/data-apis/array-api/issues/627
            warnings.warn(
                "Array API does not support LU factorization, falling back to QR"
                " instead. Set `power_iteration_normalizer='QR'` explicitly to silence"
                " this warning."
            )
            power_iteration_normalizer = "QR"
        else:
            power_iteration_normalizer = "LU"
    elif power_iteration_normalizer == "LU" and is_array_api_compliant:
        raise ValueError(
            "Array API does not support LU factorization. Set "
            "`power_iteration_normalizer='QR'` instead."
        )

    if is_array_api_compliant:
        qr_normalizer = partial(xp.linalg.qr, mode="reduced")
    else:
        # Use scipy.linalg instead of numpy.linalg when not explicitly
        # using the Array API.
        qr_normalizer = partial(linalg.qr, mode="economic", check_finite=False)

    if power_iteration_normalizer == "QR":
        normalizer = qr_normalizer
    elif power_iteration_normalizer == "LU":
        normalizer = partial(linalg.lu, permute_l=True, check_finite=False)
    else:
        normalizer = lambda x: (x, None)

    # Perform power iterations with Q to further 'imprint' the top
    # singular vectors of A in Q
    for _ in range(n_iter):
        Q, _ = normalizer(A @ Q)
        Q, _ = normalizer(A.T @ Q)

    # Sample the range of A using by linear projection of Q
    # Extract an orthonormal basis
    Q, _ = qr_normalizer(A @ Q)

    return Q