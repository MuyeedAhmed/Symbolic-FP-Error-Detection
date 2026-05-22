def test_dot(dtype):
    dot = _dot_memview[_numpy_to_cython(dtype)]

    rng = np.random.RandomState(0)
    x = rng.random_sample(10).astype(dtype, copy=False)
    y = rng.random_sample(10).astype(dtype, copy=False)

    expected = x.dot(y)
    actual = dot(x, y)

    assert_allclose(actual, expected, rtol=RTOL[dtype])