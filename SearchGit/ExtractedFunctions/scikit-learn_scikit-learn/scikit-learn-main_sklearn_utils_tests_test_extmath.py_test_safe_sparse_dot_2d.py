def test_safe_sparse_dot_2d(A_container, B_container):
    rng = np.random.RandomState(0)

    A = rng.random_sample((30, 10))
    B = rng.random_sample((10, 20))
    expected = np.dot(A, B)

    A = A_container(A)
    B = B_container(B)
    actual = safe_sparse_dot(A, B, dense_output=True)

    assert_allclose(actual, expected)