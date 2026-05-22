    def test_matmul(self, dtypes):
        import jax.numpy as jnp

        dtype1, dtype2 = dtypes
        x1 = backend.Variable("ones", shape=(1,), dtype=dtype1, trainable=False)
        x2 = backend.Variable("ones", shape=(1,), dtype=dtype2, trainable=False)
        x1_jax = jnp.ones((1,), dtype=dtype1)
        x2_jax = jnp.ones((1,), dtype=dtype2)
        expected_dtype = standardize_dtype(jnp.matmul(x1_jax, x2_jax).dtype)

        self.assertDType(x1 @ x2, expected_dtype)
        self.assertDType(x1.__rmatmul__(x2), expected_dtype)