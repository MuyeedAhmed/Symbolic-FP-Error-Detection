def test_additive_chi2_sampler(csr_container):
    # test that AdditiveChi2Sampler approximates kernel on random data

    # compute exact kernel
    # abbreviations for easier formula
    X_ = X[:, np.newaxis, :].copy()
    Y_ = Y[np.newaxis, :, :].copy()

    large_kernel = 2 * X_ * Y_ / (X_ + Y_)

    # reduce to n_samples_x x n_samples_y by summing over features
    kernel = large_kernel.sum(axis=2)

    # approximate kernel mapping
    transform = AdditiveChi2Sampler(sample_steps=3)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)

    kernel_approx = np.dot(X_trans, Y_trans.T)

    assert_array_almost_equal(kernel, kernel_approx, 1)

    X_sp_trans = transform.fit_transform(csr_container(X))
    Y_sp_trans = transform.transform(csr_container(Y))

    assert_array_equal(X_trans, X_sp_trans.toarray())
    assert_array_equal(Y_trans, Y_sp_trans.toarray())

    # test error is raised on negative input
    Y_neg = Y.copy()
    Y_neg[0, 0] = -1
    msg = "Negative values in data passed to"
    with pytest.raises(ValueError, match=msg):
        transform.fit(Y_neg)