    def test_blas64_dot(self):
        n = 2**32
        a = np.zeros([1, n], dtype=np.float32)
        b = np.ones([1, 1], dtype=np.float32)
        a[0, -1] = 1
        c = np.dot(b, a)
        assert_equal(c[0, -1], 1)