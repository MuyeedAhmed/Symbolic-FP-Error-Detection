    def test_basic_function_with_three_arguments(self):
        # multi_dot with three arguments uses a fast hand coded algorithm to
        # determine the optimal order. Therefore test it separately.
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))

        assert_almost_equal(multi_dot([A, B, C]), A.dot(B).dot(C))
        assert_almost_equal(multi_dot([A, B, C]), np.dot(A, np.dot(B, C)))