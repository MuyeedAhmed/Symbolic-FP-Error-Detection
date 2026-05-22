def test_pls_canonical_basics():
    # Basic checks for PLSCanonical
    d = load_linnerud()
    X = d.data
    y = d.target

    pls = PLSCanonical(n_components=X.shape[1])
    pls.fit(X, y)

    assert_matrix_orthogonal(pls.x_weights_)
    assert_matrix_orthogonal(pls.y_weights_)
    assert_matrix_orthogonal(pls._x_scores)
    assert_matrix_orthogonal(pls._y_scores)

    # Check X = TP' and y = UQ'
    T = pls._x_scores
    P = pls.x_loadings_
    U = pls._y_scores
    Q = pls.y_loadings_
    # Need to scale first
    Xc, yc, x_mean, y_mean, x_std, y_std = _center_scale_xy(
        X.copy(), y.copy(), scale=True
    )
    assert_array_almost_equal(Xc, np.dot(T, P.T))
    assert_array_almost_equal(yc, np.dot(U, Q.T))

    # Check that rotations on training data lead to scores
    Xt = pls.transform(X)
    assert_array_almost_equal(Xt, pls._x_scores)
    Xt, yt = pls.transform(X, y)
    assert_array_almost_equal(Xt, pls._x_scores)
    assert_array_almost_equal(yt, pls._y_scores)

    # Check that inverse_transform works
    X_back = pls.inverse_transform(Xt)
    assert_array_almost_equal(X_back, X)
    _, y_back = pls.inverse_transform(Xt, yt)
    assert_array_almost_equal(y_back, y)