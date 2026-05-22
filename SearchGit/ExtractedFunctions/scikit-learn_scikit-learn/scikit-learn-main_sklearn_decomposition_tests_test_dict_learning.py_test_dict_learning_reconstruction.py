def test_dict_learning_reconstruction():
    n_components = 12
    dico = DictionaryLearning(
        n_components, transform_algorithm="omp", transform_alpha=0.001, random_state=0
    )
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)
    assert_array_almost_equal(dico.inverse_transform(code), X)

    dico.set_params(transform_algorithm="lasso_lars")
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)
    assert_array_almost_equal(dico.inverse_transform(code), X, decimal=2)

    # test error raised for wrong code size
    with pytest.raises(ValueError, match="Expected 12, got 11."):
        dico.inverse_transform(code[:, :-1])