def test_all_precomputed():
    # Test that lars_path with precomputed Gram and Xy gives the right answer
    G = np.dot(X.T, X)
    Xy = np.dot(X.T, y)
    for method in "lar", "lasso":
        output = linear_model.lars_path(X, y, method=method)
        output_pre = linear_model.lars_path(X, y, Gram=G, Xy=Xy, method=method)
        for expected, got in zip(output, output_pre):
            assert_array_almost_equal(expected, got)