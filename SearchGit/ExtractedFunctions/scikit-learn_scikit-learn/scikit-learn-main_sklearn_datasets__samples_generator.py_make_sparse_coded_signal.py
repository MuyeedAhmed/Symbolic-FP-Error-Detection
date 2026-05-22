def make_sparse_coded_signal(
    n_samples,
    *,
    n_components,
    n_features,
    n_nonzero_coefs,
    random_state=None,
):
    """Generate a signal as a sparse combination of dictionary elements.

    Returns matrices `Y`, `D` and `X` such that `Y = XD` where `X` is of shape
    `(n_samples, n_components)`, `D` is of shape `(n_components, n_features)`, and
    each row of `X` has exactly `n_nonzero_coefs` non-zero elements.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int
        Number of samples to generate.

    n_components : int
        Number of components in the dictionary.

    n_features : int
        Number of features of the dataset to generate.

    n_nonzero_coefs : int
        Number of active (non-zero) coefficients in each sample.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for dataset creation. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    data : ndarray of shape (n_samples, n_features)
        The encoded signal (Y).

    dictionary : ndarray of shape (n_components, n_features)
        The dictionary with normalized components (D).

    code : ndarray of shape (n_samples, n_components)
        The sparse code such that each column of this matrix has exactly
        n_nonzero_coefs non-zero items (X).

    Examples
    --------
    >>> from sklearn.datasets import make_sparse_coded_signal
    >>> data, dictionary, code = make_sparse_coded_signal(
    ...     n_samples=50,
    ...     n_components=100,
    ...     n_features=10,
    ...     n_nonzero_coefs=4,
    ...     random_state=0
    ... )
    >>> data.shape
    (50, 10)
    >>> dictionary.shape
    (100, 10)
    >>> code.shape
    (50, 100)
    """
    generator = check_random_state(random_state)

    # generate dictionary
    D = generator.standard_normal(size=(n_features, n_components))
    D /= np.sqrt(np.sum((D**2), axis=0))

    # generate code
    X = np.zeros((n_components, n_samples))
    for i in range(n_samples):
        idx = np.arange(n_components)
        generator.shuffle(idx)
        idx = idx[:n_nonzero_coefs]
        X[idx, i] = generator.standard_normal(size=n_nonzero_coefs)

    # encode signal
    Y = np.dot(D, X)

    # Transpose to have shapes consistent with the rest of the API
    Y, D, X = Y.T, D.T, X.T

    return map(np.squeeze, (Y, D, X))