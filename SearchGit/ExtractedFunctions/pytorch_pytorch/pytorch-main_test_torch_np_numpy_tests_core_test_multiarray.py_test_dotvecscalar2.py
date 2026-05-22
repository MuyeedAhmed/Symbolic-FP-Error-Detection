    def test_dotvecscalar2(self):
        np.random.seed(100)
        # b1 = np.random.rand(4, 1)
        b1 = np.array([[0.54340494], [0.27836939], [0.42451759], [0.84477613]])

        # b2 = np.random.rand(1, 1)
        b2 = np.array([[0.00471886]])

        res = np.dot(b1, b2)
        tgt = np.array([[0.00256425], [0.00131359], [0.00200324], [0.00398638]])
        assert_almost_equal(res, tgt, decimal=self.N)