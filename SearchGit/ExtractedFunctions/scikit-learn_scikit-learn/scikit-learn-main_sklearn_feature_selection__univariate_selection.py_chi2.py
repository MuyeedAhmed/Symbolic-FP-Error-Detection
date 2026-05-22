def chi2(X, y):
    """Compute chi-squared stats between each non-negative feature and class.

    This score can be used to select the `n_features` features with the
    highest values for the test chi-squared statistic from X, which must
    contain only **non-negative integer feature values** such as booleans or frequencies
    (e.g., term counts in document classification), relative to the classes.

    If some of your features are continuous, you need to bin them, for
    example by using :class:`~sklearn.preprocessing.KBinsDiscretizer`.

    Recall that the chi-square test measures dependence between stochastic
    variables, so using this function "weeds out" the features that are the
    most likely to be independent of class and therefore irrelevant for
    classification.

    Read more in the :ref:`User Guide <univariate_feature_selection>`.

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Sample vectors.

    y : array-like of shape (n_samples,)
        Target vector (class labels).

    Returns
    -------
    chi2 : ndarray of shape (n_features,)
        Chi2 statistics for each feature.

    p_values : ndarray of shape (n_features,)
        P-values for each feature.

    See Also
    --------
    f_classif : ANOVA F-value between label/feature for classification tasks.
    f_regression : F-value between label/feature for regression tasks.

    Notes
    -----
    Complexity of this algorithm is O(n_classes * n_features).

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.feature_selection import chi2
    >>> X = np.array([[1, 1, 3],
    ...               [0, 1, 5],
    ...               [5, 4, 1],
    ...               [6, 6, 2],
    ...               [1, 4, 0],
    ...               [0, 0, 0]])
    >>> y = np.array([1, 1, 0, 0, 2, 2])
    >>> chi2_stats, p_values = chi2(X, y)
    >>> chi2_stats
    array([15.3,  6.5       ,  8.9])
    >>> p_values
    array([0.000456, 0.0387, 0.0116 ])
    """

    # XXX: we might want to do some of the following in logspace instead for
    # numerical stability.
    # Converting X to float allows getting better performance for the
    # safe_sparse_dot call made below.
    X = check_array(X, accept_sparse="csr", dtype=(np.float64, np.float32))
    if np.any((X.data if issparse(X) else X) < 0):
        raise ValueError("Input X must be non-negative.")

    # Use a sparse representation for Y by default to reduce memory usage when
    # y has many unique classes.
    Y = LabelBinarizer(sparse_output=True).fit_transform(y)
    if Y.shape[1] == 1:
        Y = Y.toarray()
        Y = np.append(1 - Y, Y, axis=1)

    observed = safe_sparse_dot(Y.T, X)  # n_classes * n_features

    if issparse(observed):
        # convert back to a dense array before calling _chisquare
        # XXX: could _chisquare be reimplement to accept sparse matrices for
        # cases where both n_classes and n_features are large (and X is
        # sparse)?
        observed = observed.toarray()

    feature_count = X.sum(axis=0).reshape(1, -1)
    class_prob = Y.mean(axis=0).reshape(1, -1)
    expected = np.dot(class_prob.T, feature_count)

    return _chisquare(observed, expected)