def test_label_spreading_closed_form(global_dtype, Estimator, parameters, alpha):
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200, random_state=0)
    X = X.astype(global_dtype, copy=False)
    y[::3] = -1

    gamma = 0.1
    clf = label_propagation.LabelSpreading(gamma=gamma).fit(X, y)
    # adopting notation from Zhou et al (2004):
    S = clf._build_graph()
    Y = np.zeros((len(y), n_classes + 1), dtype=X.dtype)
    Y[np.arange(len(y)), y] = 1
    Y = Y[:, :-1]

    expected = np.dot(np.linalg.inv(np.eye(len(S), dtype=S.dtype) - alpha * S), Y)
    expected /= expected.sum(axis=1)[:, np.newaxis]

    clf = label_propagation.LabelSpreading(
        max_iter=100, alpha=alpha, tol=1e-10, gamma=gamma
    )
    clf.fit(X, y)

    assert_allclose(expected, clf.label_distributions_)