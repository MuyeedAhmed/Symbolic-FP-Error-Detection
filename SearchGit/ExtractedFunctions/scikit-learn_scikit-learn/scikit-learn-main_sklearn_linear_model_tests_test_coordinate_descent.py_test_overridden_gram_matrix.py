def test_overridden_gram_matrix():
    X, y, _, _ = build_dataset(n_samples=20, n_features=10)
    Gram = X.T.dot(X)
    clf = ElasticNet(selection="cyclic", tol=1e-8, precompute=Gram)
    warning_message = (
        "Gram matrix was provided but X was centered"
        " to fit intercept: recomputing Gram matrix."
    )
    with pytest.warns(UserWarning, match=warning_message):
        clf.fit(X, y)