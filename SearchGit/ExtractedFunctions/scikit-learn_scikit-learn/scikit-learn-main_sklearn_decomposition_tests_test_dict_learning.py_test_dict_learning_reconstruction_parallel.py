def test_dict_learning_reconstruction_parallel():
    # regression test that parallel reconstruction works with n_jobs>1
    n_components = 12
    dico = DictionaryLearning(
        n_components,
        transform_algorithm="omp",
        transform_alpha=0.001,
        random_state=0,
        n_jobs=4,
    )
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.set_params(transform_algorithm="lasso_lars")
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)