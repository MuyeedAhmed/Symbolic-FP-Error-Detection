  def compute_predictive_covariance(self, gp_feature):
    """Computes posterior predictive variance.

    Approximates the Gaussian process posterior using random features.
    Given training random feature Phi_tr (num_train, num_hidden) and testing
    random feature Phi_ts (batch_size, num_hidden). The predictive covariance
    matrix is computed as (assuming Gaussian likelihood):

    s * Phi_ts @ inv(t(Phi_tr) * Phi_tr + s * I) @ t(Phi_ts),

    where s is the ridge factor to be used for stablizing the inverse, and I is
    the identity matrix with shape (num_hidden, num_hidden).

    Args:
      gp_feature: (tf.Tensor) The random feature of testing data to be used for
        computing the covariance matrix. Shape (batch_size, gp_hidden_size).

    Returns:
      (tf.Tensor) Predictive covariance matrix, shape (batch_size, batch_size).
    """
    # Computes the covariance matrix of the feature coefficient.
    feature_cov_matrix = tf.linalg.inv(self.precision_matrix)

    # Computes the covariance matrix of the gp prediction.
    cov_feature_product = tf.matmul(
        feature_cov_matrix, gp_feature, transpose_b=True) * self.ridge_penalty
    gp_cov_matrix = tf.matmul(gp_feature, cov_feature_product)
    return gp_cov_matrix