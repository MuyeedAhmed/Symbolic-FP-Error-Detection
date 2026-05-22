def test_decision_function(global_random_seed):
    iris = get_iris_dataset(global_random_seed)
    # Test decision_function
    # Sanity check, test that decision_function implemented in python
    # returns the same as the one in libsvm
    # multi class:
    clf = svm.SVC(kernel="linear", C=0.1, decision_function_shape="ovo").fit(
        iris.data, iris.target
    )

    dec = np.dot(iris.data, clf.coef_.T) + clf.intercept_

    assert_array_almost_equal(dec, clf.decision_function(iris.data))

    # binary:
    clf.fit(X, Y)
    dec = np.dot(X, clf.coef_.T) + clf.intercept_
    prediction = clf.predict(X)
    assert_array_almost_equal(dec.ravel(), clf.decision_function(X))
    assert_array_almost_equal(
        prediction, clf.classes_[(clf.decision_function(X) > 0).astype(int)]
    )
    expected = np.array([-1.0, -0.66, -1.0, 0.66, 1.0, 1.0])
    assert_array_almost_equal(clf.decision_function(X), expected, 2)

    # kernel binary:
    clf = svm.SVC(kernel="rbf", gamma=1, decision_function_shape="ovo")
    clf.fit(X, Y)

    rbfs = rbf_kernel(X, clf.support_vectors_, gamma=clf.gamma)
    dec = np.dot(rbfs, clf.dual_coef_.T) + clf.intercept_
    assert_array_almost_equal(dec.ravel(), clf.decision_function(X))