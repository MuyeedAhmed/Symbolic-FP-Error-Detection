def make_spd_matrix(n_dim, *, random_state=None):
    """Generate a random symmetric, positive-definite matrix.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_dim : int
        The matrix dimension.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for dataset creation. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    X : ndarray of shape (n_dim, n_dim)
        The random symmetric, positive-definite matrix.

    See Also
    --------
    make_sparse_spd_matrix: Generate a sparse symmetric definite positive matrix.

    Examples
    --------
    >>> from sklearn.datasets import make_spd_matrix
    >>> make_spd_matrix(n_dim=2, random_state=42)
    array([[2.093, 0.346],
           [0.346, 0.218]])
    """
    generator = check_random_state(random_state)

    A = generator.uniform(size=(n_dim, n_dim))
    U, _, Vt = linalg.svd(np.dot(A.T, A), check_finite=False)
    X = np.dot(np.dot(U, 1.0 + np.diag(generator.uniform(size=n_dim))), Vt)

    return X