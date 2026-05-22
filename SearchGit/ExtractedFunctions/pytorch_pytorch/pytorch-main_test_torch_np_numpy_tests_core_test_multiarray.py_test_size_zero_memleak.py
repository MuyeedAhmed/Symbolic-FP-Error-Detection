    def test_size_zero_memleak(self):
        # Regression test for issue 9615
        # Exercises a special-case code path for dot products of length
        # zero in cblasfuncs (making it is specific to floating dtypes).
        a = np.array([], dtype=np.float64)
        x = np.array(2.0)
        for _ in range(100):
            np.dot(a, a, out=x)
        if HAS_REFCOUNT:
            assert_(sys.getrefcount(x) < 50)