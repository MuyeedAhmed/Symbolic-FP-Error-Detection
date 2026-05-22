    def test_dot_out(self):
        # if HAVE_CBLAS, will use WRITEBACKIFCOPY
        a = np.arange(9, dtype=float).reshape(3, 3)
        b = np.dot(a, a, out=a)
        assert_equal(b, np.array([[15, 18, 21], [42, 54, 66], [69, 90, 111]]))