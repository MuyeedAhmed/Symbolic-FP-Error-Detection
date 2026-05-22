    def test_huge_vectordot(self, dtype):
        # Large vector multiplications are chunked with 32bit BLAS
        # Test that the chunking does the right thing, see also gh-22262
        data = np.ones(2**30 + 100, dtype=dtype)
        res = np.dot(data, data)
        if res != 2**30 + 100:
            raise AssertionError(f"result mismatch: {res} != {2**30 + 100}")