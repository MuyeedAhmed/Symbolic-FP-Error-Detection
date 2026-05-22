  def test_laplace_covariance_minibatch(self, generate_orthogonal_data):
    """Tests if model correctly learns population-lvel precision matrix."""
    batch_size = 50
    epochs = 1000
    x_data = _generate_rbf_data(self.x_ts, generate_orthogonal_data)
    data_iterator = _make_minibatch_iterator(x_data, batch_size, epochs)

    # Estimates precision matrix using minibatch.
    cov_estimator = gaussian_process.LaplaceRandomFeatureCovariance(
        momentum=0.999, ridge_penalty=0)

    for minibatch_data in data_iterator:
      _ = cov_estimator(minibatch_data, training=True)

    # Evaluation
    prec_mat_expected = x_data.T.dot(x_data)
    prec_mat_computed = (
        cov_estimator.precision_matrix.numpy() * self.num_test_sample)

    np.testing.assert_allclose(prec_mat_computed, prec_mat_expected,
                               **self.prec_tolerance)