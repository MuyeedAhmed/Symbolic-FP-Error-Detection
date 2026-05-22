def test_predict_proba_multilabel():
    # Test that predict_proba works as expected for multilabel.
    # Multilabel should not use softmax which makes probabilities sum to 1
    X, Y = make_multilabel_classification(
        n_samples=50, random_state=0, return_indicator=True
    )
    n_samples, n_classes = Y.shape

    clf = MLPClassifier(solver="lbfgs", hidden_layer_sizes=30, random_state=0)
    clf.fit(X, Y)
    y_proba = clf.predict_proba(X)

    assert y_proba.shape == (n_samples, n_classes)
    assert_array_equal(y_proba > 0.5, Y)

    y_log_proba = clf.predict_log_proba(X)
    proba_max = y_proba.argmax(axis=1)
    proba_log_max = y_log_proba.argmax(axis=1)

    assert (y_proba.sum(1) - 1).dot(y_proba.sum(1) - 1) > 1e-10
    assert_array_equal(proba_max, proba_log_max)
    assert_allclose(y_log_proba, np.log(y_proba))