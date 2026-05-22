def test_randomized_eigsh(dtype):
    """Test that `_randomized_eigsh` returns the appropriate components"""

    rng = np.random.RandomState(42)
    X = np.diag(np.array([1.0, -2.0, 0.0, 3.0], dtype=dtype))
    # random rotation that preserves the eigenvalues of X
    rand_rot = np.linalg.qr(rng.normal(size=X.shape))[0]
    X = rand_rot @ X @ rand_rot.T

    # with 'module' selection method, the negative eigenvalue shows up
    eigvals, eigvecs = _randomized_eigsh(X, n_components=2, selection="module")
    # eigenvalues
    assert eigvals.shape == (2,)
    assert_array_almost_equal(eigvals, [3.0, -2.0])  # negative eigenvalue here
    # eigenvectors
    assert eigvecs.shape == (4, 2)

    # with 'value' selection method, the negative eigenvalue does not show up
    with pytest.raises(NotImplementedError):
        _randomized_eigsh(X, n_components=2, selection="value")