    def test_types(self, dtype):
        x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
        assert_equal(linalg.inv(x).dtype, dtype)