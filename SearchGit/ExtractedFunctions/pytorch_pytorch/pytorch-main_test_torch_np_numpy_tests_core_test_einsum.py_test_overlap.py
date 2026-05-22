    def test_overlap(self):
        a = np.arange(9, dtype=int).reshape(3, 3)
        b = np.arange(9, dtype=int).reshape(3, 3)
        d = np.dot(a, b)
        # sanity check
        c = np.einsum("ij,jk->ik", a, b)
        assert_equal(c, d)
        # gh-10080, out overlaps one of the operands
        c = np.einsum("ij,jk->ik", a, b, out=b)
        assert_equal(c, d)