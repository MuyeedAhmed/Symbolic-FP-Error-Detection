def test_safe_sparse_dot_2d_1d(container):
    rng = np.random.RandomState(0)
    B = rng.random_sample((10))

    # 2D @ 1D
    A = rng.random_sample((30, 10))
    expected = np.dot(A, B)
    actual = safe_sparse_dot(container(A), B)
    assert_allclose(actual, expected)

    # 1D @ 2D
    A = rng.random_sample((10, 30))
    expected = np.dot(B, A)
    actual = safe_sparse_dot(B, container(A))
    assert_allclose(actual, expected)