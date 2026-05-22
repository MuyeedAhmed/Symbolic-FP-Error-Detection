def jax_model_no_state_apply(params, inputs):
    activations = inputs.reshape((inputs.shape[0], -1))  # flatten
    for w, b in params[:-1]:
        outputs = jnp.dot(activations, w) + b
        activations = jnp.tanh(outputs)

    final_w, final_b = params[-1]
    logits = jnp.dot(activations, final_w) + final_b
    return jax.nn.softmax(logits, axis=-1)