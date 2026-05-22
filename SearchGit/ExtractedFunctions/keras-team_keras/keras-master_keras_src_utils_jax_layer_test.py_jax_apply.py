        def jax_apply(params, inputs):
            return jnp.dot(inputs, params[0])