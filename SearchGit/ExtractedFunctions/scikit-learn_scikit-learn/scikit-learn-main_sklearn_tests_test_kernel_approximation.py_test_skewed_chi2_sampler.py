def test_skewed_chi2_sampler():
    # test that RBFSampler approximates kernel on random data

    # compute exact kernel
    c = 0.03
    # set on negative component but greater than c to ensure that the kernel
    # approximation is valid on the group (-c; +\infty) endowed with the skewed
    # multiplication.
    Y_ = Y.copy()
    Y_[0, 0] = -c / 2.0

    # abbreviations for easier formula
    X_c = (X + c)[:, np.newaxis, :]
    Y_c = (Y_ + c)[np.newaxis, :, :]

    # we do it in log-space in the hope that it's more stable
    # this array is n_samples_x x n_samples_y big x n_features
    log_kernel = (
        (np.log(X_c) / 2.0) + (np.log(Y_c) / 2.0) + np.log(2.0) - np.log(X_c + Y_c)
    )
    # reduce to n_samples_x x n_samples_y by summing over features in log-space
    kernel = np.exp(log_kernel.sum(axis=2))

    # approximate kernel mapping
    transform = SkewedChi2Sampler(skewedness=c, n_components=1000, random_state=42)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y_)

    kernel_approx = np.dot(X_trans, Y_trans.T)
    assert_array_almost_equal(kernel, kernel_approx, 1)
    assert np.isfinite(kernel).all(), "NaNs found in the Gram matrix"
    assert np.isfinite(kernel_approx).all(), "NaNs found in the approximate Gram matrix"

    # test error is raised on when inputs contains values smaller than -c
    Y_neg = Y_.copy()
    Y_neg[0, 0] = -c * 2.0
    msg = "X may not contain entries smaller than -skewedness"
    with pytest.raises(ValueError, match=msg):
        transform.transform(Y_neg)