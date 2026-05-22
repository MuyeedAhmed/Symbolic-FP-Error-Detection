        def assert_dot_close(A, X, desired):
            assert_allclose(np.dot(A, X), desired, rtol=1e-5, atol=1e-7)