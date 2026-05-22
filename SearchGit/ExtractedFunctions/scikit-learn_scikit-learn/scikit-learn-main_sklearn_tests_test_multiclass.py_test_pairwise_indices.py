def test_pairwise_indices():
    clf_precomputed = svm.SVC(kernel="precomputed")
    X, y = iris.data, iris.target

    ovr_false = OneVsOneClassifier(clf_precomputed)
    linear_kernel = np.dot(X, X.T)
    ovr_false.fit(linear_kernel, y)

    n_estimators = len(ovr_false.estimators_)
    precomputed_indices = ovr_false.pairwise_indices_

    for idx in precomputed_indices:
        assert (
            idx.shape[0] * n_estimators / (n_estimators - 1) == linear_kernel.shape[0]
        )