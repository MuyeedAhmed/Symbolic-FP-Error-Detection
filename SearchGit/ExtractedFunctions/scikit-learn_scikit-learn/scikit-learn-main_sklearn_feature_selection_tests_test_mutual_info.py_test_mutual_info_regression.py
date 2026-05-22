def test_mutual_info_regression(global_dtype):
    # We generate sample from multivariate normal distribution, using
    # transformation from initially uncorrelated variables. The zero
    # variables after transformation is selected as the target vector,
    # it has the strongest correlation with the variable 2, and
    # the weakest correlation with the variable 1.
    T = np.array([[1, 0.5, 2, 1], [0, 1, 0.1, 0.0], [0, 0.1, 1, 0.1], [0, 0.1, 0.1, 1]])
    cov = T.dot(T.T)
    mean = np.zeros(4)

    rng = check_random_state(0)
    Z = rng.multivariate_normal(mean, cov, size=1000).astype(global_dtype, copy=False)
    X = Z[:, 1:]
    y = Z[:, 0]

    mi = mutual_info_regression(X, y, random_state=0)
    assert_array_equal(np.argsort(-mi), np.array([1, 2, 0]))
    # XXX: should mutual_info_regression be fixed to avoid
    # up-casting float32 inputs to float64?
    assert mi.dtype == np.float64