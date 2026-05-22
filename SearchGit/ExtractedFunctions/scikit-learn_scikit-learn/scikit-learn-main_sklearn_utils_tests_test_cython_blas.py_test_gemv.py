def test_gemv(dtype, opA, transA, order):
    gemv = _gemv_memview[_numpy_to_cython(dtype)]

    rng = np.random.RandomState(0)
    A = np.asarray(
        opA(rng.random_sample((20, 10)).astype(dtype, copy=False)), order=ORDER[order]
    )
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(20).astype(dtype, copy=False)
    alpha, beta = 2.5, -0.5

    expected = alpha * opA(A).dot(x) + beta * y
    gemv(transA, alpha, A, x, beta, y)

    assert_allclose(y, expected, rtol=RTOL[dtype])