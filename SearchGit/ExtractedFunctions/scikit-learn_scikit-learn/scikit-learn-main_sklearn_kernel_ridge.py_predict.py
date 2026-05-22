    def predict(self, X):
        """Predict using the kernel ridge model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples. If kernel == "precomputed" this is instead a
            precomputed kernel matrix, shape = [n_samples,
            n_samples_fitted], where n_samples_fitted is the number of
            samples used in the fitting for this estimator.

        Returns
        -------
        C : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Returns predicted values.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse=("csr", "csc"), reset=False)
        K = self._get_kernel(X, self.X_fit_)
        return np.dot(K, self.dual_coef_)