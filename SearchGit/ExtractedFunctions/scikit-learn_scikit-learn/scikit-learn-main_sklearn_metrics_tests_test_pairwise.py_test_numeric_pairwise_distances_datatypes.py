def test_numeric_pairwise_distances_datatypes(metric, global_dtype, y_is_x):
    # Check that pairwise distances gives the same result as pdist and cdist
    # regardless of input datatype when using any scipy metric for comparing
    # numeric vectors
    #
    # This test is necessary because pairwise_distances used to throw an
    # error when using metric='seuclidean' and the input data was not
    # of type np.float64 (#15730)

    rng = np.random.RandomState(0)

    X = rng.random_sample((5, 4)).astype(global_dtype, copy=False)

    params = {}
    if y_is_x:
        Y = X
        expected_dist = squareform(pdist(X, metric=metric))
    else:
        Y = rng.random_sample((5, 4)).astype(global_dtype, copy=False)
        expected_dist = cdist(X, Y, metric=metric)
        # precompute parameters for seuclidean & mahalanobis when x is not y
        if metric == "seuclidean":
            params = {"V": np.var(np.vstack([X, Y]), axis=0, ddof=1, dtype=np.float64)}
        elif metric == "mahalanobis":
            params = {"VI": np.linalg.inv(np.cov(np.vstack([X, Y]).T)).T}

    dist = pairwise_distances(X, Y, metric=metric, **params)

    assert_allclose(dist, expected_dist)