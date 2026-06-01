    def test_vector_as_first_argument(self):
        # The first argument can be 1-D
        A1d = np.random.random(2)  # 1-D
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))
        D = np.random.random((2, 2))

        # the result should be 1-D
        assert_equal(multi_dot([A1d, B, C, D]).shape, (2,))