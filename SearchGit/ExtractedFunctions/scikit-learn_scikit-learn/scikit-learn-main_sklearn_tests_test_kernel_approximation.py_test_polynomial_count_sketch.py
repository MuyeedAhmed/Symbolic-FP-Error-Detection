def test_polynomial_count_sketch(gamma, degree, coef0, n_components):
    # test that PolynomialCountSketch approximates polynomial
    # kernel on random data

    # compute exact kernel
    kernel = polynomial_kernel(X, Y, gamma=gamma, degree=degree, coef0=coef0)

    # approximate kernel mapping
    ps_transform = PolynomialCountSketch(
        n_components=n_components,
        gamma=gamma,
        coef0=coef0,
        degree=degree,
        random_state=42,
    )
    X_trans = ps_transform.fit_transform(X)
    Y_trans = ps_transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)

    error = kernel - kernel_approx
    assert np.abs(np.mean(error)) <= 0.05  # close to unbiased
    np.abs(error, out=error)
    assert np.max(error) <= 0.1  # nothing too far off
    assert np.mean(error) <= 0.05  # mean is fairly close