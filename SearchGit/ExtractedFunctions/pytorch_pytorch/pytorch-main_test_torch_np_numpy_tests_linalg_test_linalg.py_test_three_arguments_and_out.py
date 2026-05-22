    def test_three_arguments_and_out(self):
        # multi_dot with three arguments uses a fast hand coded algorithm to
        # determine the optimal order. Therefore test it separately.
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))

        out = np.zeros((6, 2))
        ret = multi_dot([A, B, C], out=out)
        if out is not ret:
            raise AssertionError("Expected out is ret")
        assert_almost_equal(out, A.dot(B).dot(C))
        assert_almost_equal(out, np.dot(A, np.dot(B, C)))