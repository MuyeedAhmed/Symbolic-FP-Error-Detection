def test_bad_input(lil_container, global_random_seed):
    # Test dimensions for labels
    Y2 = Y[:-1]  # wrong dimensions for labels
    with pytest.raises(ValueError):
        svm.SVC().fit(X, Y2)

    # Test with arrays that are non-contiguous.
    for clf in (svm.SVC(), svm.LinearSVC(random_state=global_random_seed)):
        Xf = np.asfortranarray(X)
        assert not Xf.flags["C_CONTIGUOUS"]
        yf = np.ascontiguousarray(np.tile(Y, (2, 1)).T)
        yf = yf[:, -1]
        assert not yf.flags["F_CONTIGUOUS"]
        assert not yf.flags["C_CONTIGUOUS"]
        clf.fit(Xf, yf)
        assert_array_equal(clf.predict(T), true_result)

    # error for precomputed kernelsx
    clf = svm.SVC(kernel="precomputed")
    with pytest.raises(ValueError):
        clf.fit(X, Y)

    # predict with sparse input when trained with dense
    clf = svm.SVC().fit(X, Y)
    with pytest.raises(ValueError):
        clf.predict(lil_container(X))

    Xt = np.array(X).T
    clf.fit(np.dot(X, Xt), Y)
    with pytest.raises(ValueError):
        clf.predict(X)

    clf = svm.SVC()
    clf.fit(X, Y)
    with pytest.raises(ValueError):
        clf.predict(Xt)