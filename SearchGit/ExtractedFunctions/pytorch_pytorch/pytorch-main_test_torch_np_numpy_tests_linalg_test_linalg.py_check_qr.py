    def check_qr(self, a):
        # This test expects the argument `a` to be an ndarray or
        # a subclass of an ndarray of inexact type.
        a_type = type(a)
        a_dtype = a.dtype
        m, n = a.shape
        k = min(m, n)

        # mode == 'complete'
        q, r = linalg.qr(a, mode="complete")
        assert_(q.dtype == a_dtype)
        assert_(r.dtype == a_dtype)
        assert_(isinstance(q, a_type))
        assert_(isinstance(r, a_type))
        assert_(q.shape == (m, m))
        assert_(r.shape == (m, n))
        assert_almost_equal(dot(q, r), a, single_decimal=5)
        assert_almost_equal(dot(q.T.conj(), q), np.eye(m))
        assert_almost_equal(np.triu(r), r)

        # mode == 'reduced'
        q1, r1 = linalg.qr(a, mode="reduced")
        assert_(q1.dtype == a_dtype)
        assert_(r1.dtype == a_dtype)
        assert_(isinstance(q1, a_type))
        assert_(isinstance(r1, a_type))
        assert_(q1.shape == (m, k))
        assert_(r1.shape == (k, n))
        assert_almost_equal(dot(q1, r1), a, single_decimal=5)
        assert_almost_equal(dot(q1.T.conj(), q1), np.eye(k))
        assert_almost_equal(np.triu(r1), r1)

        # mode == 'r'
        r2 = linalg.qr(a, mode="r")
        assert_(r2.dtype == a_dtype)
        assert_(isinstance(r2, a_type))
        assert_almost_equal(r2, r1)