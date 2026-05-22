def test_linearsvc_crammer_singer(global_random_seed):
    # Test LinearSVC with crammer_singer multi-class svm
    iris = get_iris_dataset(global_random_seed)

    ovr_clf = svm.LinearSVC(random_state=global_random_seed).fit(iris.data, iris.target)
    cs_clf = svm.LinearSVC(
        multi_class="crammer_singer", random_state=global_random_seed
    )
    cs_clf.fit(iris.data, iris.target)

    # similar prediction for ovr and crammer-singer:
    assert (ovr_clf.predict(iris.data) == cs_clf.predict(iris.data)).mean() > 0.9

    # classifiers shouldn't be the same
    assert (ovr_clf.coef_ != cs_clf.coef_).all()

    # test decision function
    assert_array_equal(
        cs_clf.predict(iris.data),
        np.argmax(cs_clf.decision_function(iris.data), axis=1),
    )
    dec_func = np.dot(iris.data, cs_clf.coef_.T) + cs_clf.intercept_
    assert_array_almost_equal(dec_func, cs_clf.decision_function(iris.data))