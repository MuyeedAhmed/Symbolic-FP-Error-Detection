def test_timeout(lil_container):
    sp = svm.SVC(C=1, kernel=lambda x, y: x @ y.T, random_state=0, max_iter=1)
    warning_msg = (
        r"Solver terminated early \(max_iter=1\).  Consider pre-processing "
        r"your data with StandardScaler or MinMaxScaler."
    )
    with pytest.warns(ConvergenceWarning, match=warning_msg):
        sp.fit(lil_container(X), Y)