def make_low_rank_matrix(
    n_samples=100,
    n_features=100,
    *,
    effective_rank=10,
    tail_strength=0.5,
    random_state=None,
):
    """Generate a mostly low rank matrix with bell-shaped singular values.

    Most of the variance can be explained by a bell-shaped curve of width
    effective_rank: the low rank part of the singular values profile is::

        (1 - tail_strength) * exp(-1.0 * (i / effective_rank) ** 2)

    The remaining singular values' tail is fat, decreasing as::

        tail_strength * exp(-0.1 * i / effective_rank).

    The low rank part of the profile can be considered the structured
    signal part of the data while the tail can be considered the noisy
    part of the data that cannot be summarized by a low number of linear
    components (singular vectors).

    This kind of singular profiles is often seen in practice, for instance:
     - gray level pictures of faces
     - TF-IDF vectors of text documents crawled from the web

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, default=100
        The number of samples.

    n_features : int, default=100
        The number of features.

    effective_rank : int, default=10
        The approximate number of singular vectors required to explain most of
        the data by linear combinations.

    tail_strength : float, default=0.5
        The relative importance of the fat noisy tail of the singular values
        profile. The value should be between 0 and 1.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for dataset creation. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    X : ndarray of shape (n_samples, n_features)
        The matrix.

    Examples
    --------
    >>> from numpy.linalg import svd
    >>> from sklearn.datasets import make_low_rank_matrix
    >>> X = make_low_rank_matrix(
    ...     n_samples=50,
    ...     n_features=25,
    ...     effective_rank=5,
    ...     tail_strength=0.01,
    ...     random_state=0,
    ... )
    >>> X.shape
    (50, 25)
    """
    generator = check_random_state(random_state)
    n = min(n_samples, n_features)

    # Random (ortho normal) vectors
    u, _ = linalg.qr(
        generator.standard_normal(size=(n_samples, n)),
        mode="economic",
        check_finite=False,
    )
    v, _ = linalg.qr(
        generator.standard_normal(size=(n_features, n)),
        mode="economic",
        check_finite=False,
    )

    # Index of the singular values
    singular_ind = np.arange(n, dtype=np.float64)

    # Build the singular profile by assembling signal and noise components
    low_rank = (1 - tail_strength) * np.exp(-1.0 * (singular_ind / effective_rank) ** 2)
    tail = tail_strength * np.exp(-0.1 * singular_ind / effective_rank)
    s = np.identity(n) * (low_rank + tail)

    return np.dot(np.dot(u, s), v.T)