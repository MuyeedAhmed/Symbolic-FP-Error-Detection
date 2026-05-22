    def get_precision(self):
        """Compute data precision matrix with the FactorAnalysis model.

        Returns
        -------
        precision : ndarray of shape (n_features, n_features)
            Estimated precision of data.
        """
        check_is_fitted(self)

        n_features = self.components_.shape[1]

        # handle corner cases first
        if self.n_components == 0:
            return np.diag(1.0 / self.noise_variance_)
        if self.n_components == n_features:
            return linalg.inv(self.get_covariance())

        # Get precision using matrix inversion lemma
        components_ = self.components_
        precision = np.dot(components_ / self.noise_variance_, components_.T)
        precision.flat[:: len(precision) + 1] += 1.0
        precision = np.dot(components_.T, np.dot(linalg.inv(precision), components_))
        precision /= self.noise_variance_[:, np.newaxis]
        precision /= -self.noise_variance_[np.newaxis, :]
        precision.flat[:: len(precision) + 1] += 1.0 / self.noise_variance_
        return precision