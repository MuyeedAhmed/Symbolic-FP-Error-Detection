    def test_dotmatvec(self):
        A, b1 = self.A, self.b1
        res = np.dot(A, b1)
        tgt = np.array([[0.32114320], [0.04889721], [0.15696029], [0.33612621]])
        assert_almost_equal(res, tgt, decimal=self.N)