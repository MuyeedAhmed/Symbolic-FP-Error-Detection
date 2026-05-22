def fowlkes_mallows_score(labels_true, labels_pred):
    """Measure the similarity of two clusterings of a set of points.

    .. versionadded:: 0.18

    The Fowlkes-Mallows index (FMI) is defined as the geometric mean of
    the precision and recall::

        FMI = TP / sqrt((TP + FP) * (TP + FN))

    Where ``TP`` is the number of **True Positive** (i.e. the number of pairs of
    points that belong to the same cluster in both ``labels_true`` and
    ``labels_pred``), ``FP`` is the number of **False Positive** (i.e. the
    number of pairs of points that belong to the same cluster in
    ``labels_pred`` but not in ``labels_true``) and ``FN`` is the number of
    **False Negative** (i.e. the number of pairs of points that belong to the
    same cluster in ``labels_true`` but not in ``labels_pred``).

    The score ranges from 0 to 1. A high value indicates a good similarity
    between two clusters.

    Read more in the :ref:`User Guide <fowlkes_mallows_scores>`.

    Parameters
    ----------
    labels_true : array-like of shape (n_samples,)
        A clustering of the data into disjoint subsets.

    labels_pred : array-like of shape (n_samples,)
        A clustering of the data into disjoint subsets.

    Returns
    -------
    score : float
       The resulting Fowlkes-Mallows score.

    References
    ----------
    .. [1] `E. B. Fowkles and C. L. Mallows, 1983. "A method for comparing two
       hierarchical clusterings". Journal of the American Statistical
       Association
       <https://www.tandfonline.com/doi/abs/10.1080/01621459.1983.10478008>`_

    .. [2] `Wikipedia entry for the Fowlkes-Mallows Index
           <https://en.wikipedia.org/wiki/Fowlkes-Mallows_index>`_

    Examples
    --------

    Perfect labelings are both homogeneous and complete, hence have
    score 1.0::

      >>> from sklearn.metrics.cluster import fowlkes_mallows_score
      >>> fowlkes_mallows_score([0, 0, 1, 1], [0, 0, 1, 1])
      1.0
      >>> fowlkes_mallows_score([0, 0, 1, 1], [1, 1, 0, 0])
      1.0

    If classes members are completely split across different clusters,
    the assignment is totally random, hence the FMI is null::

      >>> fowlkes_mallows_score([0, 0, 0, 0], [0, 1, 2, 3])
      0.0
    """

    labels_true, labels_pred = check_clusterings(labels_true, labels_pred)
    (n_samples,) = labels_true.shape

    c = contingency_matrix(labels_true, labels_pred, sparse=True)
    c = c.astype(np.int64, copy=False)
    tk = np.dot(c.data, c.data) - n_samples
    pk = np.sum(np.asarray(c.sum(axis=0)).ravel() ** 2) - n_samples
    qk = np.sum(np.asarray(c.sum(axis=1)).ravel() ** 2) - n_samples
    return float(np.sqrt(tk / pk) * np.sqrt(tk / qk)) if tk != 0.0 else 0.0