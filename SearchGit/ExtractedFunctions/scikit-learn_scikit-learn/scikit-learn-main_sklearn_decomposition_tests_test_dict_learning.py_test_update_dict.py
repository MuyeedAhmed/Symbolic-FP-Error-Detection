def test_update_dict():
    # Check the dict update in batch mode vs online mode
    # Non-regression test for #4866
    rng = np.random.RandomState(0)

    code = np.array([[0.5, -0.5], [0.1, 0.9]])
    dictionary = np.array([[1.0, 0.0], [0.6, 0.8]])

    X = np.dot(code, dictionary) + rng.randn(2, 2)

    # full batch update
    newd_batch = dictionary.copy()
    _update_dict(newd_batch, X, code)

    # online update
    A = np.dot(code.T, code)
    B = np.dot(X.T, code)
    newd_online = dictionary.copy()
    _update_dict(newd_online, X, code, A, B)

    assert_allclose(newd_batch, newd_online)