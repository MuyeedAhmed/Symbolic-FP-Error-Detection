  def optimize(self, sess, feed_dict):
    old_theta = sess.run(self.flat_vars)
    loss_flat_grad = sess.run(self.loss_flat_gradient,
                              feed_dict=feed_dict)

    def calc_fisher_vector_product(tangent):
      feed_dict[self.flat_tangent] = tangent
      fvp = sess.run(self.fisher_vector_product,
                     feed_dict=feed_dict)
      fvp += self.cg_damping * tangent
      return fvp

    step_dir = conjugate_gradient(calc_fisher_vector_product, -loss_flat_grad)

    shs = 0.5 * step_dir.dot(calc_fisher_vector_product(step_dir))
    lm = np.sqrt(shs / self.max_divergence)
    fullstep = step_dir / lm
    neggdotstepdir = -loss_flat_grad.dot(step_dir)

    def calc_loss(theta):
      sess.run(self.set_vars, feed_dict={self.flat_theta: theta})
      if self.divergence is None:
        return sess.run(self.raw_loss, feed_dict=feed_dict), True
      else:
        raw_loss, divergence = sess.run(
            [self.raw_loss, self.divergence], feed_dict=feed_dict)
        return raw_loss, divergence < self.max_divergence

    # find optimal theta
    theta = linesearch(calc_loss, old_theta, fullstep, neggdotstepdir / lm)
    if self.divergence is not None:
      final_divergence = sess.run(self.divergence, feed_dict=feed_dict)
    else:
      final_divergence = None

    # set vars accordingly
    if final_divergence is None or final_divergence < self.max_divergence:
      sess.run(self.set_vars, feed_dict={self.flat_theta: theta})
    else:
      sess.run(self.set_vars, feed_dict={self.flat_theta: old_theta})