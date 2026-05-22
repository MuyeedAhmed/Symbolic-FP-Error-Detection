    def test_basic_function_with_two_arguments(self):
        # separate code path with two arguments
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))

        assert_almost_equal(multi_dot([A, B]), A.dot(B))
        assert_almost_equal(multi_dot([A, B]), np.dot(A, B))