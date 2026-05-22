  def test_random_feature_prior_approximation(self):
    """Tests random feature GP's ability in approximating exact GP prior."""
    num_inducing = 10240
    rfgp_model = gaussian_process.RandomFeatureGaussianProcess(
        units=1,
        num_inducing=num_inducing,
        normalize_input=False,
        gp_kernel_type='gaussian',
        return_random_features=True)

    # Extract random features.
    _, _, gp_feature = rfgp_model(self.x_tr, training=True)
    gp_feature_np = gp_feature.numpy()

    prior_kernel_computed = gp_feature_np.dot(gp_feature_np.T)
    prior_kernel_expected = self.rbf_kern_func(self.x_tr, self.x_tr)
    np.testing.assert_allclose(prior_kernel_computed, prior_kernel_expected,
                               **self.cov_tolerance)