    def test_matmul_bool(self):
        # gh-14439
        a = np.array([[1, 0], [1, 1]], dtype=bool)
        if np.max(a.view(np.uint8)) != 1:
            raise AssertionError(f"max mismatch: {np.max(a.view(np.uint8))} != 1")
        b = np.matmul(a, a)
        # matmul with boolean output should always be 0, 1
        if np.max(b.view(np.uint8)) != 1:
            raise AssertionError(f"max mismatch: {np.max(b.view(np.uint8))} != 1")

        # rg = np.random.default_rng(np.random.PCG64(43))
        # d = rg.integers(2, size=4*5, dtype=np.int8)
        # d = d.reshape(4, 5) > 0
        np.random.seed(1234)
        d = np.random.randint(2, size=(4, 5)) > 0

        out1 = np.matmul(d, d.reshape(5, 4))
        out2 = np.dot(d, d.reshape(5, 4))
        assert_equal(out1, out2)

        c = np.matmul(np.zeros((2, 0), dtype=bool), np.zeros(0, dtype=bool))
        if np.any(c):
            raise AssertionError("c should be all False")