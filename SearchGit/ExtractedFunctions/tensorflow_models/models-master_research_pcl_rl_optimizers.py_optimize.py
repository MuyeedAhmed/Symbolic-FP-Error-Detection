  def optimize(self, sess, feed_dict):
    reg_input, reg_weight, old_values, targets = sess.run(
        [self.inputs, self.regression_weight, self.values, self.targets],
        feed_dict=feed_dict)

    intended_values = targets * self.mix_frac + old_values * (1 - self.mix_frac)

    # taken from rllab
    reg_coeff = 1e-5
    for _ in range(5):
      best_fit_weight = np.linalg.lstsq(
          reg_input.T.dot(reg_input) +
          reg_coeff * np.identity(reg_input.shape[1]),
          reg_input.T.dot(intended_values))[0]
      if not np.any(np.isnan(best_fit_weight)):
        break
      reg_coeff *= 10

    if len(best_fit_weight.shape) == 1:
      best_fit_weight = np.expand_dims(best_fit_weight, -1)

    sess.run(self.update_regression_weight,
             feed_dict={self.new_regression_weight: best_fit_weight})