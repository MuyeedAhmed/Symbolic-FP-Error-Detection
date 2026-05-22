def test_gemm(dtype, opA, transA, opB, transB, order):
    gemm = _gemm_memview[_numpy_to_cython(dtype)]

    rng = np.random.RandomState(0)
    A = np.asarray(
        opA(rng.random_sample((30, 10)).astype(dtype, copy=False)), order=ORDER[order]
    )
    B = np.asarray(
        opB(rng.random_sample((10, 20)).astype(dtype, copy=False)), order=ORDER[order]
    )
    C = np.asarray(
        rng.random_sample((30, 20)).astype(dtype, copy=False), order=ORDER[order]
    )
    alpha, beta = 2.5, -0.5

    expected = alpha * opA(A).dot(opB(B)) + beta * C
    gemm(transA, transB, alpha, A, B, beta, C)

    assert_allclose(C, expected, rtol=RTOL[dtype])