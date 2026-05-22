    def test_dotvecscalar(self):
        np.random.seed(100)
        # Numpy guarantees the random stream, and we don't. So inline the
        # values from numpy 1.24.1
        # b1 = np.random.rand(1, 1)
        b1 = np.array([[0.54340494]])

        # b2 = np.random.rand(1, 4)
        b2 = np.array([[0.27836939, 0.42451759, 0.84477613, 0.00471886]])

        res = np.dot(b1, b2)
        tgt = np.array([[0.15126730, 0.23068496, 0.45905553, 0.00256425]])
        assert_almost_equal(res, tgt, decimal=self.N)