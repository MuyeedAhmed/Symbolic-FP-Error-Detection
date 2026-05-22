    def _estimate_wishart_tied(self, nk, xk, sk):
        """Estimate the tied Wishart distribution parameters.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        nk : array-like of shape (n_components,)

        xk : array-like of shape (n_components, n_features)

        sk : array-like of shape (n_features, n_features)
        """
        _, n_features = xk.shape

        # Warning : in some Bishop book, there is a typo on the formula 10.63
        # `degrees_of_freedom_k = degrees_of_freedom_0 + Nk`
        # is the correct formula
        self.degrees_of_freedom_ = (
            self.degrees_of_freedom_prior_ + nk.sum() / self.n_components
        )

        diff = xk - self.mean_prior_
        self.covariances_ = (
            self.covariance_prior_
            + sk * nk.sum() / self.n_components
            + self.mean_precision_prior_
            / self.n_components
            * np.dot((nk / self.mean_precision_) * diff.T, diff)
        )

        # Contrary to the original bishop book, we normalize the covariances
        self.covariances_ /= self.degrees_of_freedom_