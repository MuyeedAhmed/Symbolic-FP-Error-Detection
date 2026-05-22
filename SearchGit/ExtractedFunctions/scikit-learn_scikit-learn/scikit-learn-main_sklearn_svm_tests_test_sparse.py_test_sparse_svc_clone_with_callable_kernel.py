def test_sparse_svc_clone_with_callable_kernel(lil_container):
    # Test that the "dense_fit" is called even though we use sparse input
    # meaning that everything works fine.
    a = svm.SVC(C=1, kernel=lambda x, y: x @ y.T, probability=True, random_state=0)
    b = base.clone(a)

    X_sp = lil_container(X)
    b.fit(X_sp, Y)
    pred = b.predict(X_sp)
    b.predict_proba(X_sp)

    dense_svm = svm.SVC(
        C=1, kernel=lambda x, y: np.dot(x, y.T), probability=True, random_state=0
    )
    pred_dense = dense_svm.fit(X, Y).predict(X)
    assert_array_equal(pred_dense, pred)