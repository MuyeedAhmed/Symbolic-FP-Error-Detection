def test_gradient():
    # Test gradient of Kullback-Leibler divergence.
    random_state = check_random_state(0)

    n_samples = 50
    n_features = 2
    n_components = 2
    alpha = 1.0

    distances = random_state.randn(n_samples, n_features).astype(np.float32)
    distances = np.abs(distances.dot(distances.T))
    np.fill_diagonal(distances, 0.0)
    X_embedded = random_state.randn(n_samples, n_components).astype(np.float32)

    P = _joint_probabilities(distances, desired_perplexity=25.0, verbose=0)

    def fun(params):
        return _kl_divergence(params, P, alpha, n_samples, n_components)[0]

    def grad(params):
        return _kl_divergence(params, P, alpha, n_samples, n_components)[1]

    assert_almost_equal(check_grad(fun, grad, X_embedded.ravel()), 0.0, decimal=5)