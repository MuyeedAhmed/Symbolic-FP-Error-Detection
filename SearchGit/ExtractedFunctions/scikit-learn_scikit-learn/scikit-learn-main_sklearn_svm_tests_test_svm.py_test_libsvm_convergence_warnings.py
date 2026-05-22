def test_libsvm_convergence_warnings(global_random_seed):
    a = svm.SVC(
        kernel=lambda x, y: np.dot(x, y.T),
        random_state=global_random_seed,
        max_iter=2,
    )
    warning_msg = (
        r"Solver terminated early \(max_iter=2\).  Consider pre-processing "
        r"your data with StandardScaler or MinMaxScaler."
    )
    with pytest.warns(ConvergenceWarning, match=warning_msg):
        a.fit(np.array(X), Y)
    assert np.all(a.n_iter_ == 2)