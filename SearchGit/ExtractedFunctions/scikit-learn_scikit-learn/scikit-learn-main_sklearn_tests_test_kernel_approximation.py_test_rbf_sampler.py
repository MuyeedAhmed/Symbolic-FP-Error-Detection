def test_rbf_sampler():
    # test that RBFSampler approximates kernel on random data
    # compute exact kernel
    gamma = 10.0
    kernel = rbf_kernel(X, Y, gamma=gamma)

    # approximate kernel mapping
    rbf_transform = RBFSampler(gamma=gamma, n_components=1000, random_state=42)
    X_trans = rbf_transform.fit_transform(X)
    Y_trans = rbf_transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)

    error = kernel - kernel_approx
    assert np.abs(np.mean(error)) <= 0.01  # close to unbiased
    np.abs(error, out=error)
    assert np.max(error) <= 0.1  # nothing too far off
    assert np.mean(error) <= 0.05  # mean is fairly close