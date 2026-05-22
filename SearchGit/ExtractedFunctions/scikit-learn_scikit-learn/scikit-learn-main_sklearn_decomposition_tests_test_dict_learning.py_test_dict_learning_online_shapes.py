def test_dict_learning_online_shapes():
    rng = np.random.RandomState(0)
    n_components = 8

    code, dictionary = dict_learning_online(
        X,
        n_components=n_components,
        batch_size=4,
        max_iter=10,
        method="cd",
        random_state=rng,
        return_code=True,
    )
    assert code.shape == (n_samples, n_components)
    assert dictionary.shape == (n_components, n_features)
    assert np.dot(code, dictionary).shape == X.shape

    dictionary = dict_learning_online(
        X,
        n_components=n_components,
        batch_size=4,
        max_iter=10,
        method="cd",
        random_state=rng,
        return_code=False,
    )
    assert dictionary.shape == (n_components, n_features)