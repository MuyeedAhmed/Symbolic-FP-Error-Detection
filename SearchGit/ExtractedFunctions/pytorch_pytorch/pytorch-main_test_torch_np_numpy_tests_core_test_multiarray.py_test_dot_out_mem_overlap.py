    def test_dot_out_mem_overlap(self):
        np.random.seed(1)

        # Test BLAS and non-BLAS code paths, including all dtypes
        # that dot() supports
        dtypes = [np.dtype(code) for code in np.typecodes["All"] if code not in "USVM"]
        for dtype in dtypes:
            a = np.random.rand(3, 3).astype(dtype)

            # Valid dot() output arrays must be aligned
            b = _aligned_zeros((3, 3), dtype=dtype)
            b[...] = np.random.rand(3, 3)

            y = np.dot(a, b)
            x = np.dot(a, b, out=b)
            assert_equal(x, y, err_msg=repr(dtype))

            # Check invalid output array
            assert_raises(ValueError, np.dot, a, b, out=b[::2])
            assert_raises(ValueError, np.dot, a, b, out=b.T)