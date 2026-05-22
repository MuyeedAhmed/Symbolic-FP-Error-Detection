    def test_dot_array_order(self):
        a = np.array([[1, 2], [3, 4]], order="C")
        b = np.array([[1, 2], [3, 4]], order="F")
        res = np.dot(a, a)

        # integer arrays are exact
        assert_equal(np.dot(a, b), res)
        assert_equal(np.dot(b, a), res)
        assert_equal(np.dot(b, b), res)