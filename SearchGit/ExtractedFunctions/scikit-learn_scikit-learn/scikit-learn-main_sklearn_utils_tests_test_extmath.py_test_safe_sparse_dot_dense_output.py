def test_safe_sparse_dot_dense_output(dense_output):
    rng = np.random.RandomState(0)

    A = _sparse_random_array((30, 10), density=0.1, rng=rng)
    B = _sparse_random_array((10, 20), density=0.1, rng=rng)

    expected = A.dot(B)
    actual = safe_sparse_dot(A, B, dense_output=dense_output)

    assert sparse.issparse(actual) == (not dense_output)

    if dense_output:
        expected = expected.toarray()
    assert_allclose_dense_sparse(actual, expected)