def test_randomized_eigsh_reconst_low_rank(n, rank):
    """Check that randomized_eigsh is able to reconstruct a low rank psd matrix

    Tests that the decomposition provided by `_randomized_eigsh` leads to
    orthonormal eigenvectors, and that a low rank PSD matrix can be effectively
    reconstructed with good accuracy using it.
    """
    assert rank < n

    # create a low rank PSD
    rng = np.random.RandomState(69)
    X = rng.randn(n, rank)
    A = X @ X.T

    # approximate A with the "right" number of components
    S, V = _randomized_eigsh(A, n_components=rank, random_state=rng)
    # orthonormality checks
    assert_array_almost_equal(np.linalg.norm(V, axis=0), np.ones(S.shape))
    assert_array_almost_equal(V.T @ V, np.diag(np.ones(S.shape)))
    # reconstruction
    A_reconstruct = V @ np.diag(S) @ V.T

    # test that the approximation is good
    assert_array_almost_equal(A_reconstruct, A, decimal=6)