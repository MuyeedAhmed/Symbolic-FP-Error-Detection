    def test_vector_as_first_and_last_argument(self):
        # The first and last arguments can be 1-D
        A1d = np.random.random(2)  # 1-D
        B = np.random.random((2, 6))
        C = np.random.random((6, 2))
        D1d = np.random.random(2)  # 1-D

        # the result should be a scalar
        assert_equal(multi_dot([A1d, B, C, D1d]).shape, ())