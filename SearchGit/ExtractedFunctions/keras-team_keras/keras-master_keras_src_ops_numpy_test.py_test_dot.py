    def test_dot(self, dtypes):
        import jax.numpy as jnp

        dtype1, dtype2 = dtypes
        x1 = knp.ones((2, 3, 4), dtype=dtype1)
        x2 = knp.ones((4, 3), dtype=dtype2)
        x1_jax = jnp.ones((2, 3, 4), dtype=dtype1)
        x2_jax = jnp.ones((4, 3), dtype=dtype2)
        expected_dtype = standardize_dtype(jnp.dot(x1_jax, x2_jax).dtype)

        self.assertEqual(
            standardize_dtype(knp.dot(x1, x2).dtype), expected_dtype
        )
        self.assertEqual(knp.Dot().symbolic_call(x1, x2).dtype, expected_dtype)