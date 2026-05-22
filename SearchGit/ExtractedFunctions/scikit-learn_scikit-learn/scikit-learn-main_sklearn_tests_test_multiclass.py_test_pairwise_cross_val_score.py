def test_pairwise_cross_val_score(MultiClassClassifier):
    clf_precomputed = svm.SVC(kernel="precomputed")
    clf_notprecomputed = svm.SVC(kernel="linear")

    X, y = iris.data, iris.target

    multiclass_clf_notprecomputed = MultiClassClassifier(clf_notprecomputed)
    multiclass_clf_precomputed = MultiClassClassifier(clf_precomputed)

    linear_kernel = np.dot(X, X.T)
    score_not_precomputed = cross_val_score(
        multiclass_clf_notprecomputed, X, y, error_score="raise"
    )
    score_precomputed = cross_val_score(
        multiclass_clf_precomputed, linear_kernel, y, error_score="raise"
    )
    assert_array_equal(score_precomputed, score_not_precomputed)