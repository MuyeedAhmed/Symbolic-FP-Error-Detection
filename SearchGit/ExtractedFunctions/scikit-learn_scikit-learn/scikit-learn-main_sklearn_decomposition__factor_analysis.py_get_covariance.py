    def get_covariance(self):
        """Compute data covariance with the FactorAnalysis model.

        ``cov = components_.T * components_ + diag(noise_variance)``

        Returns
        -------
        cov : ndarray of shape (n_features, n_features)
            Estimated covariance of data.
        """
        check_is_fitted(self)

        cov = np.dot(self.components_.T, self.components_)
        cov.flat[:: len(cov) + 1] += self.noise_variance_  # modify diag inplace
        return cov