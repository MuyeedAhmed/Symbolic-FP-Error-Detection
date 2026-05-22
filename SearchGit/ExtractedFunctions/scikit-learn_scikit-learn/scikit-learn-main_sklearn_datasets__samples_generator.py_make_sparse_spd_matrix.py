def make_sparse_spd_matrix(
    n_dim=1,
    *,
    alpha=0.95,
    norm_diag=False,
    smallest_coef=0.1,
    largest_coef=0.9,
    sparse_format=None,
    random_state=None,
):
    """Generate a sparse symmetric definite positive matrix.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_dim : int, default=1
        The size of the random matrix to generate.

        .. versionchanged:: 1.4
            Renamed from ``dim`` to ``n_dim``.

    alpha : float, default=0.95
        The probability that a coefficient is zero (see notes). Larger values
        enforce more sparsity. The value should be in the range 0 and 1.

    norm_diag : bool, default=False
        Whether to normalize the output matrix to make the leading diagonal
        elements all 1.

    smallest_coef : float, default=0.1
        The value of the smallest coefficient between 0 and 1.

    largest_coef : float, default=0.9
        The value of the largest coefficient between 0 and 1.

    sparse_format : str, default=None
        String representing the output sparse format, such as 'csc', 'csr', etc.
        If ``None``, return a dense numpy ndarray.

        .. versionadded:: 1.4

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for dataset creation. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    prec : ndarray or sparse matrix of shape (dim, dim)
        The generated matrix. If ``sparse_format=None``, this would be an ndarray.
        Otherwise, this will be a sparse matrix of the specified format.

    See Also
    --------
    make_spd_matrix : Generate a random symmetric, positive-definite matrix.

    Notes
    -----
    The sparsity is actually imposed on the cholesky factor of the matrix.
    Thus alpha does not translate directly into the filling fraction of
    the matrix itself.

    Examples
    --------
    >>> from sklearn.datasets import make_sparse_spd_matrix
    >>> make_sparse_spd_matrix(n_dim=4, norm_diag=False, random_state=42)
    array([[1., 0., 0., 0.],
           [0., 1., 0., 0.],
           [0., 0., 1., 0.],
           [0., 0., 0., 1.]])
    """
    random_state = check_random_state(random_state)

    chol = -_sparse_eye_array(n_dim)
    aux = _sparse_random_array(
        shape=(n_dim, n_dim),
        density=1 - alpha,
        data_sampler=lambda size: random_state.uniform(
            low=smallest_coef, high=largest_coef, size=size
        ),
        random_state=random_state,
    )
    # We need to avoid "coo" format because it does not support slicing
    aux = sp.tril(aux, k=-1, format="csc")

    # Permute the lines: we don't want to have asymmetries in the final
    # SPD matrix
    permutation = random_state.permutation(n_dim)
    aux = aux[permutation].T[permutation]
    chol += aux
    prec = chol.T @ chol

    if norm_diag:
        # Form the diagonal vector into a row matrix
        d = _sparse_diags_array(1.0 / np.sqrt(prec.diagonal()))
        prec = d @ prec @ d

    if sparse_format is None:
        return prec.toarray()
    else:
        return _align_api_if_sparse(prec.asformat(sparse_format))