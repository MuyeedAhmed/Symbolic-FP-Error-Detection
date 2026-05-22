    def test_dotvecvecinner(self):
        b1, b3 = self.b1, self.b3
        res = np.dot(b3, b1)
        tgt = np.array([[0.23129668]])
        assert_almost_equal(res, tgt, decimal=self.N)