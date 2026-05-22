def test_cross_val_score_precomputed():
    # test for svm with precomputed kernel
    svm = SVC(kernel="precomputed")
    iris = load_iris()
    X, y = iris.data, iris.target
    linear_kernel = np.dot(X, X.T)
    score_precomputed = cross_val_score(svm, linear_kernel, y)
    svm = SVC(kernel="linear")
    score_linear = cross_val_score(svm, X, y)
    assert_array_almost_equal(score_precomputed, score_linear)

    # test with callable
    svm = SVC(kernel=lambda x, y: np.dot(x, y.T))
    score_callable = cross_val_score(svm, X, y)
    assert_array_almost_equal(score_precomputed, score_callable)

    # Error raised for non-square X
    svm = SVC(kernel="precomputed")
    with pytest.raises(ValueError):
        cross_val_score(svm, X, y)

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    with pytest.raises(ValueError):
        cross_val_score(svm, linear_kernel.tolist(), y)