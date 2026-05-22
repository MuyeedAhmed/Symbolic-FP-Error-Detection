    def test_out_arg(self):
        a = np.ones((5, 2), dtype=float)
        b = np.array([[1, 3], [5, 7]], dtype=float)
        tgt = np.dot(a, b)

        # test as positional argument
        msg = "out positional argument"
        out = np.zeros((5, 2), dtype=float)
        self.matmul(a, b, out)
        assert_array_equal(out, tgt, err_msg=msg)

        # test as keyword argument
        msg = "out keyword argument"
        out = np.zeros((5, 2), dtype=float)
        self.matmul(a, b, out=out)
        assert_array_equal(out, tgt, err_msg=msg)

        # test out with not allowed type cast (safe casting)
        msg = "Cannot cast"
        out = np.zeros((5, 2), dtype=np.int32)
        assert_raises_regex(TypeError, msg, self.matmul, a, b, out=out)

        # test out with type upcast to complex
        out = np.zeros((5, 2), dtype=np.complex128)
        c = self.matmul(a, b, out=out)
        assert_(c is out)
        c = c.astype(tgt.dtype)
        assert_array_equal(c, tgt)