    def do(self, a, b, tags):
        arr = np.asarray(a)
        m, n = arr.shape
        u, s, vt = linalg.svd(a, False)
        x, residuals, rank, sv = linalg.lstsq(a, b, rcond=-1)
        if m == 0:
            assert_((x == 0).all())
        if m <= n:
            assert_almost_equal(b, dot(a, x), single_decimal=5)
            assert_equal(rank, m)
        else:
            assert_equal(rank, n)
        #     assert_almost_equal(sv, sv.__array_wrap__(s))
        if rank == n and m > n:
            expect_resids = (np.asarray(abs(np.dot(a, x) - b)) ** 2).sum(axis=0)
            expect_resids = np.asarray(expect_resids)
            if np.asarray(b).ndim == 1:
                expect_resids = expect_resids.reshape(
                    1,
                )
                assert_equal(residuals.shape, expect_resids.shape)
        else:
            expect_resids = np.array([])  # .view(type(x))
        assert_almost_equal(residuals, expect_resids, single_decimal=5)
        assert_(np.issubdtype(residuals.dtype, np.floating))
        assert_(consistent_subclass(x, b))
        assert_(consistent_subclass(residuals, b))