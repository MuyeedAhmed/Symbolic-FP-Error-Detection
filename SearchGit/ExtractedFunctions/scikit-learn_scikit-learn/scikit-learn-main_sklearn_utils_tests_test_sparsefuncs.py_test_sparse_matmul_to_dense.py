def test_sparse_matmul_to_dense(
    global_random_seed, out_is_None, a_container, b_container, dtype
):
    """Test that sparse_matmul_to_dense computes correctly."""
    rng = np.random.default_rng(global_random_seed)
    n1, n2, n3 = 10, 19, 13
    a_dense = rng.standard_normal((n1, n2)).astype(dtype)
    b_dense = rng.standard_normal((n2, n3)).astype(dtype)
    a_dense.flat[rng.choice([False, True], size=n1 * n2, p=[0.5, 0.5])] = 0
    b_dense.flat[rng.choice([False, True], size=n2 * n3, p=[0.5, 0.5])] = 0
    a = a_container(a_dense)
    b = b_container(b_dense)
    if out_is_None:
        out = None
    else:
        out = np.empty((n1, n3), dtype=dtype)

    result = sparse_matmul_to_dense(a, b, out=out)
    # Use atol to account for the wide range of values in the computed matrix.
    assert_allclose(result, a_dense @ b_dense, atol=1e-7)
    if not out_is_None:
        assert_allclose(out, result, atol=1e-7)