def test_lars_precompute(classifier):
    # Check for different values of precompute
    G = np.dot(X.T, X)

    clf = classifier(precompute=G)
    output_1 = ignore_warnings(clf.fit)(X, y).coef_
    for precompute in [True, False, "auto", None]:
        clf = classifier(precompute=precompute)
        output_2 = clf.fit(X, y).coef_
        assert_array_almost_equal(output_1, output_2, decimal=8)