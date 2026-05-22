def _compute_posterior_kernel(x_tr, x_ts, kernel_func, ridge_penalty):
  """Computes the posterior covariance matrix of a Gaussian process."""
  num_sample = x_tr.shape[0]

  k_tt_inv = tf.linalg.inv(
      kernel_func(x_tr, x_tr) + ridge_penalty * np.eye(num_sample))
  k_ts = kernel_func(x_tr, x_ts)
  k_ss = kernel_func(x_ts, x_ts)

  return k_ss - tf.matmul(k_ts, tf.matmul(k_tt_inv, k_ts), transpose_a=True)