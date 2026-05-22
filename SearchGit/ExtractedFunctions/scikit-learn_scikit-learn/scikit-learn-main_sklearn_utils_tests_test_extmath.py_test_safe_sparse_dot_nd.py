def test_safe_sparse_dot_nd(csr_container):
    rng = np.random.RandomState(0)

    # dense ND / sparse
    A = rng.random_sample((2, 3, 4, 5, 6))
    B = rng.random_sample((6, 7))
    expected = np.dot(A, B)
    B = csr_container(B)
    actual = safe_sparse_dot(A, B)
    assert_allclose(actual, expected)

    # sparse / dense ND
    A = rng.random_sample((2, 3))
    B = rng.random_sample((4, 5, 3, 6))
    expected = np.dot(A, B)
    A = csr_container(A)
    actual = safe_sparse_dot(A, B)
    assert_allclose(actual, expected)