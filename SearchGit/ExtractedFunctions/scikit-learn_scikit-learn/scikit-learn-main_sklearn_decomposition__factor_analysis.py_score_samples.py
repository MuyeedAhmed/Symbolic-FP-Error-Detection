    def score_samples(self, X):
        """Compute the log-likelihood of each sample.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            The data.

        Returns
        -------
        ll : ndarray of shape (n_samples,)
            Log-likelihood of each sample under the current model.
        """
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)
        Xr = X - self.mean_
        precision = self.get_precision()
        n_features = X.shape[1]
        log_like = -0.5 * (Xr * (np.dot(Xr, precision))).sum(axis=1)
        log_like -= 0.5 * (n_features * log(2.0 * np.pi) - fast_logdet(precision))
        return log_like