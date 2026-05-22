    def test_dotmatmat(self):
        A = self.A
        res = np.dot(A.transpose(), A)
        tgt = np.array([[1.45046013, 0.86323640], [0.86323640, 0.84934569]])
        assert_almost_equal(res, tgt, decimal=self.N)