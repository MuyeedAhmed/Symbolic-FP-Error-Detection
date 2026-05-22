    def test_zero_dimension(self):
        # Test resolution to issue #5663
        a = np.zeros((3, 0))
        b = np.zeros((0, 4))
        td = np.tensordot(a, b, (1, 0))
        assert_array_equal(td, np.dot(a, b))