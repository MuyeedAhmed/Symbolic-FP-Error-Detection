    def decode(
        self,
        inputs: jax.Array,
    ) -> jax.Array:
        return jnp.dot(inputs, self.embeddings.T.astype(inputs.dtype))