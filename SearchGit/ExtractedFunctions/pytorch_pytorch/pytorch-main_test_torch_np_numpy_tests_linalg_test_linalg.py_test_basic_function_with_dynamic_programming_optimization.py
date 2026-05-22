    def test_basic_function_with_dynamic_programming_optimization(self):
        # multi_dot with four or more arguments uses the dynamic programming
        # optimization and therefore deserve a separate
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))
        D = np.random.random((2, 1))
        assert_almost_equal(multi_dot([A, B, C, D]), A.dot(B).dot(C).dot(D))