    def test_dynamic_programming_optimization_and_out(self):
        # multi_dot with four or more arguments uses the dynamic programming
        # optimization and therefore deserve a separate test
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))
        D = np.random.random((2, 1))
        out = np.zeros((6, 1))
        ret = multi_dot([A, B, C, D], out=out)
        if out is not ret:
            raise AssertionError("Expected out is ret")
        assert_almost_equal(out, A.dot(B).dot(C).dot(D))