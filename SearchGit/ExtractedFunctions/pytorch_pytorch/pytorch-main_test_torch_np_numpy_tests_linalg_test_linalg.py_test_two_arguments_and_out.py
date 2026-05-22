    def test_two_arguments_and_out(self):
        # separate code path with two arguments
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        out = np.zeros((6, 6))
        ret = multi_dot([A, B], out=out)
        if out is not ret:
            raise AssertionError("Expected out is ret")
        assert_almost_equal(out, A.dot(B))
        assert_almost_equal(out, np.dot(A, B))