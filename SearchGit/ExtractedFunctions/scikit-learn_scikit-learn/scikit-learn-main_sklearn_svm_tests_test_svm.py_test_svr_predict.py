def test_svr_predict(global_random_seed):
    # Test SVR's decision_function
    # Sanity check, test that predict implemented in python
    # returns the same as the one in libsvm
    iris = get_iris_dataset(global_random_seed)

    X = iris.data
    y = iris.target

    # linear kernel
    reg = svm.SVR(kernel="linear", C=0.1).fit(X, y)

    dec = np.dot(X, reg.coef_.T) + reg.intercept_
    assert_array_almost_equal(dec.ravel(), reg.predict(X).ravel())

    # rbf kernel
    reg = svm.SVR(kernel="rbf", gamma=1).fit(X, y)

    rbfs = rbf_kernel(X, reg.support_vectors_, gamma=reg.gamma)
    dec = np.dot(rbfs, reg.dual_coef_.T) + reg.intercept_
    assert_array_almost_equal(dec.ravel(), reg.predict(X).ravel())