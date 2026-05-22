def test_svc_clone_with_callable_kernel():
    iris = get_iris_dataset(42)

    # create SVM with callable linear kernel, check that results are the same
    # as with built-in linear kernel
    svm_callable = svm.SVC(
        kernel=lambda x, y: np.dot(x, y.T),
        random_state=0,
        decision_function_shape="ovr",
    )
    # clone for checking clonability with lambda functions..
    svm_cloned = base.clone(svm_callable)
    svm_cloned.fit(iris.data, iris.target)

    svm_builtin = svm.SVC(
        kernel="linear",
        random_state=0,
        decision_function_shape="ovr",
    )

    svm_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(svm_cloned.dual_coef_, svm_builtin.dual_coef_)
    assert_array_almost_equal(svm_cloned.intercept_, svm_builtin.intercept_)
    assert_array_equal(svm_cloned.predict(iris.data), svm_builtin.predict(iris.data))

    assert_array_almost_equal(
        svm_cloned.decision_function(iris.data),
        svm_builtin.decision_function(iris.data),
    )