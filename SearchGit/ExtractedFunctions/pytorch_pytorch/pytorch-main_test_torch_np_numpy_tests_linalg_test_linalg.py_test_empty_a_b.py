    def test_empty_a_b(self, m, n, n_rhs):
        a = np.arange(m * n).reshape(m, n)
        b = np.ones((m, n_rhs))
        x, residuals, rank, s = linalg.lstsq(a, b, rcond=None)
        if m == 0:
            assert_((x == 0).all())
        assert_equal(x.shape, (n, n_rhs))
        assert_equal(residuals.shape, ((n_rhs,) if m > n else (0,)))
        if m > n and n_rhs > 0:
            # residuals are exactly the squared norms of b's columns
            r = b - np.dot(a, x)
            assert_almost_equal(residuals, (r * r).sum(axis=-2))
        assert_equal(rank, min(m, n))
        assert_equal(s.shape, (min(m, n),))