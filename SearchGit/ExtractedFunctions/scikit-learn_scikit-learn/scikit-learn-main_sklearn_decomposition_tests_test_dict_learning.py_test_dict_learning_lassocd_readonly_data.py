def test_dict_learning_lassocd_readonly_data():
    n_components = 12
    with TempMemmap(X) as X_read_only:
        dico = DictionaryLearning(
            n_components,
            transform_algorithm="lasso_cd",
            transform_alpha=0.001,
            random_state=0,
            n_jobs=4,
        )
        with ignore_warnings(category=ConvergenceWarning):
            code = dico.fit(X_read_only).transform(X_read_only)
        assert_array_almost_equal(
            np.dot(code, dico.components_), X_read_only, decimal=2
        )