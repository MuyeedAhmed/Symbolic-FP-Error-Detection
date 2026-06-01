    def test_vector_as_last_argument(self):
        # The last argument can be 1-D
        A = np.random.random((6, 2))
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))
        D1d = np.random.random(2)  # 1-D

        # the result should be 1-D
        assert_equal(multi_dot([A, B, C, D1d]).shape, (6,))